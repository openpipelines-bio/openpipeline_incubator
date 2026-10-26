import json
import sys
from pathlib import Path

import pytest
from insitupy import InSituData, InSituExperiment
from insitupy.io import read_xenium

## VIASH START
meta = {
    "name": "insitupy_annotate_cells",
    "resources_dir": "resources_test",
}
## VIASH END
xenium_bundle = Path(meta["resources_dir"]) / "xenium_multicellseg_tiny"

crop_xlim = (250, 600)
crop_ylim = (300, 550)


# There are no insitupy resources on s3 yet, so the fixtures below convert the
# Xenium bundle themselves.
@pytest.fixture(scope="module")
def insitupy_insitudata(tmp_path_factory):
    project = tmp_path_factory.mktemp("insitupy_project") / "sample_1"
    xd = read_xenium(xenium_bundle, verbose=False)
    xd.saveas(project, verbose=False)
    return project


@pytest.fixture(scope="module")
def insitupy_experiment(tmp_path_factory):
    experiment = tmp_path_factory.mktemp("insitupy_experiment") / "experiment"
    xd = read_xenium(xenium_bundle, verbose=False)

    exp = InSituExperiment()
    exp.add(data=xd, metadata={"sample_id": "sample_1"})
    exp.add(
        data=xd.crop(xlim=crop_xlim, ylim=crop_ylim),
        metadata={"sample_id": "sample_2"},
    )
    exp.saveas(experiment, verbose=False)
    return experiment


# QuPath-style GeoJSON with one polygon covering the whole tiny dataset. The
# coordinates are in um, so the component is run with a scale factor of 1.
@pytest.fixture(scope="module")
def shapes_geojson(tmp_path_factory):
    shapes = tmp_path_factory.mktemp("shapes") / "shapes.geojson"
    feature = {
        "type": "Feature",
        "id": "tissue-1",
        "geometry": {
            "type": "Polygon",
            "coordinates": [[[0, 0], [2000, 0], [2000, 2000], [0, 2000], [0, 0]]],
        },
        "properties": {
            "objectType": "annotation",
            "classification": {"name": "tissue"},
        },
    }
    shapes.write_text(json.dumps({"type": "FeatureCollection", "features": [feature]}))
    return shapes


def annotation_args(shapes: Path):
    return [
        "--input_annotations",
        str(shapes),
        "--annotation_keys",
        "test_annot",
        "--annotation_scale_factor",
        "1",
    ]


def region_args(shapes: Path):
    return [
        "--input_regions",
        str(shapes),
        "--region_keys",
        "test_region",
        "--region_scale_factor",
        "1",
    ]


def assert_annotations(xd: InSituData):
    assert "test_annot" in xd.annotations.keys()
    assert "annotations-test_annot" in xd.cells.table.obs.columns


def assert_regions(xd: InSituData):
    assert "test_region" in xd.regions.keys()
    assert "regions-test_region" in xd.cells.table.obs.columns


def run_insitudata(run_component, input_project: Path, output: Path, args: list):
    run_component(["--input", str(input_project), "--output", str(output), *args])
    assert output.is_dir(), "Output project was not created."
    return InSituData.read(output)


def run_experiment(run_component, input_experiment: Path, output: Path, args: list):
    run_component(["--input", str(input_experiment), "--output", str(output), *args])
    assert output.is_dir(), "Output experiment was not created."
    exp = InSituExperiment.read(output)
    exp.load_all()
    return [xd for _, xd in exp.iterdata()]


def test_insitudata_annotations(
    run_component, insitupy_insitudata, shapes_geojson, tmp_path
):
    xd = run_insitudata(
        run_component,
        insitupy_insitudata,
        tmp_path / "annotated_insitudata",
        annotation_args(shapes_geojson),
    )
    assert_annotations(xd)


def test_insitudata_regions(
    run_component, insitupy_insitudata, shapes_geojson, tmp_path
):
    xd = run_insitudata(
        run_component,
        insitupy_insitudata,
        tmp_path / "annotated_insitudata",
        region_args(shapes_geojson),
    )
    assert_regions(xd)


def test_experiment_annotations(
    run_component, insitupy_experiment, shapes_geojson, tmp_path
):
    samples = run_experiment(
        run_component,
        insitupy_experiment,
        tmp_path / "annotated_experiment",
        annotation_args(shapes_geojson),
    )
    for xd in samples:
        assert_annotations(xd)


def test_experiment_regions(
    run_component, insitupy_experiment, shapes_geojson, tmp_path
):
    samples = run_experiment(
        run_component,
        insitupy_experiment,
        tmp_path / "annotated_experiment",
        region_args(shapes_geojson),
    )
    for xd in samples:
        assert_regions(xd)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-W", "ignore::pytest.PytestAssertRewriteWarning"]))
