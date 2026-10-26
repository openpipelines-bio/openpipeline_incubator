import json
import sys
from pathlib import Path

import pytest
from insitupy import InSituData, InSituExperiment
from insitupy.io import read_xenium

## VIASH START
meta = {
    "name": "insitupy_calc_distance",
    "resources_dir": "resources_test",
}
## VIASH END
xenium_bundle = Path(meta["resources_dir"]) / "xenium_multicellseg_tiny"

crop_xlim = (250, 600)
crop_ylim = (300, 550)


# QuPath-style GeoJSON with a single rectangular polygon, in um.
def write_geojson(path: Path, name: str, xlim: tuple, ylim: tuple):
    (x0, x1), (y0, y1) = xlim, ylim
    feature = {
        "type": "Feature",
        "id": f"{name}-1",
        "geometry": {
            "type": "Polygon",
            "coordinates": [[[x0, y0], [x1, y0], [x1, y1], [x0, y1], [x0, y0]]],
        },
        "properties": {
            "objectType": "annotation",
            "classification": {"name": name},
        },
    }
    path.write_text(json.dumps({"type": "FeatureCollection", "features": [feature]}))
    return path


# There are no insitupy resources on s3 yet, so the fixture below converts the
# Xenium bundle itself and adds a small "tumor" annotation (inside the crop, so
# both experiment samples have it) and a "tissue" region covering everything.
@pytest.fixture(scope="module")
def annotated_xd(tmp_path_factory):
    shapes = tmp_path_factory.mktemp("shapes")
    annot = write_geojson(shapes / "annot.geojson", "tumor", (350, 450), (350, 450))
    region = write_geojson(shapes / "region.geojson", "tissue", (0, 2000), (0, 2000))

    xd = read_xenium(xenium_bundle, verbose=False)
    xd.import_annotations(files=annot, keys="test_annot", scale_factor=1)
    xd.import_regions(files=region, keys="test_region", scale_factor=1)
    return xd


@pytest.fixture(scope="module")
def insitupy_insitudata(tmp_path_factory, annotated_xd):
    project = tmp_path_factory.mktemp("insitupy_project") / "sample_1"
    annotated_xd.saveas(project, verbose=False)
    return project


@pytest.fixture(scope="module")
def insitupy_experiment(tmp_path_factory, annotated_xd):
    experiment = tmp_path_factory.mktemp("insitupy_experiment") / "experiment"
    exp = InSituExperiment()
    exp.add(data=annotated_xd, metadata={"sample_id": "sample_1"})
    exp.add(
        data=annotated_xd.crop(xlim=crop_xlim, ylim=crop_ylim),
        metadata={"sample_id": "sample_2"},
    )
    exp.saveas(experiment, verbose=False)
    return experiment


def assert_distance(xd: InSituData, key: str):
    obsm = xd.cells.table.obsm
    assert "distance_from" in obsm, "No distances found in .obsm."
    assert key in obsm["distance_from"].columns, (
        f"Key '{key}' not found in .obsm['distance_from']."
    )
    assert obsm["distance_from"][key].notna().any(), f"Distances '{key}' are all NaN."


def test_insitudata(run_component, insitupy_insitudata, tmp_path):
    output = tmp_path / "distance_insitudata"
    run_component(
        [
            "--input",
            str(insitupy_insitudata),
            "--annotation_tuple",
            "test_annot;tumor",
            "--output",
            str(output),
        ]
    )
    assert output.is_dir(), "Output project was not created."
    # Without --key_to_save, the annotation name is used as key.
    assert_distance(InSituData.read(output), "tumor")


def test_insitudata_region(run_component, insitupy_insitudata, tmp_path):
    output = tmp_path / "distance_insitudata"
    run_component(
        [
            "--input",
            str(insitupy_insitudata),
            "--annotation_tuple",
            "test_annot;tumor",
            "--region_tuple",
            "test_region;tissue",
            "--key_to_save",
            "distance_to_tumor",
            "--output",
            str(output),
        ]
    )
    assert output.is_dir(), "Output project was not created."
    assert_distance(InSituData.read(output), "distance_to_tumor")


def test_experiment(run_component, insitupy_experiment, tmp_path):
    output = tmp_path / "distance_experiment"
    run_component(
        [
            "--input",
            str(insitupy_experiment),
            "--annotation_tuple",
            "test_annot;tumor",
            "--output",
            str(output),
        ]
    )
    assert output.is_dir(), "Output experiment was not created."
    exp = InSituExperiment.read(output)
    exp.load_all()
    for _, xd in exp.iterdata():
        assert_distance(xd, "tumor")


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-W", "ignore::pytest.PytestAssertRewriteWarning"]))
