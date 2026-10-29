import json
import sys
from pathlib import Path

import numpy as np
import pytest
import scanpy as sc
from insitupy import InSituExperiment
from insitupy.io import read_xenium

## VIASH START
meta = {
    "name": "insitupy_dge",
    "resources_dir": "resources_test",
}
## VIASH END
xenium_bundle = Path(meta["resources_dir"]) / "xenium_multicellseg_tiny"

crop_xlim = (250, 600)
crop_ylim = (300, 550)


# QuPath-style GeoJSON with one rectangular polygon per class, in um.
def write_geojson(path: Path, polygons: dict):
    features = []
    for name, ((x0, x1), (y0, y1)) in polygons.items():
        features.append(
            {
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
        )
    path.write_text(json.dumps({"type": "FeatureCollection", "features": features}))
    return path


# There are no insitupy resources on s3 yet, so the fixture below converts the
# Xenium bundle itself and adds a "tumor" (~50 cells) and "healthy" (~200 cells)
# annotation (both inside the crop, so both experiment samples have them), a
# "tissue" region covering everything and a dummy "cell_type" column (A/B).
# The matrix is log1p-normalized, as required by dge.
@pytest.fixture(scope="module")
def annotated_xd(tmp_path_factory):
    shapes = tmp_path_factory.mktemp("shapes")
    annot = write_geojson(
        shapes / "annot.geojson",
        {"tumor": ((350, 450), (350, 450)), "healthy": ((250, 350), (300, 550))},
    )
    region = write_geojson(
        shapes / "region.geojson", {"tissue": ((0, 2000), (0, 2000))}
    )

    xd = read_xenium(xenium_bundle, verbose=False)
    xd.import_annotations(files=annot, keys="test_annot", scale_factor=1)
    xd.import_regions(files=region, keys="test_region", scale_factor=1)

    adata = xd.cells.table
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    adata.obs["cell_type"] = np.where(np.arange(adata.n_obs) % 2 == 0, "A", "B")
    return xd


@pytest.fixture(scope="module")
def insitupy_insitudata(tmp_path_factory, annotated_xd):
    project = tmp_path_factory.mktemp("insitupy_project") / "sample_1"
    annotated_xd.saveas(project, verbose=False)
    return project


@pytest.fixture(scope="module")
def insitupy_insitudata_ref(tmp_path_factory, annotated_xd):
    project = tmp_path_factory.mktemp("insitupy_project_ref") / "sample_2"
    annotated_xd.crop(xlim=crop_xlim, ylim=crop_ylim).saveas(project, verbose=False)
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


def assert_outputs(output: Path, files=("main.csv", "config.toml")):
    assert output.is_dir(), "Output folder was not created."
    for file in files:
        assert (output / file).is_file(), f"'{file}' was not created."


def test_insitudata(run_component, insitupy_insitudata, tmp_path):
    output = tmp_path / "dge_insitudata"
    run_component(
        [
            "--input",
            str(insitupy_insitudata),
            "--target_annotation_tuple",
            "test_annot;tumor",
            "--ref_annotation_tuple",
            "rest",
            "--output",
            str(output),
        ]
    )
    assert_outputs(output)


def test_insitudata_ref_annotation(run_component, insitupy_insitudata, tmp_path):
    output = tmp_path / "dge_ref_annotation"
    run_component(
        [
            "--input",
            str(insitupy_insitudata),
            "--target_annotation_tuple",
            "test_annot;tumor",
            "--target_region_tuple",
            "test_region;tissue",
            "--target_name",
            "tumor",
            "--target_metadata",
            "condition=treated",
            "--ref_annotation_tuple",
            "test_annot;healthy",
            "--ref_name",
            "healthy",
            "--ref_metadata",
            "condition=control",
            "--output",
            str(output),
        ]
    )
    assert_outputs(output)


def test_insitudata_cell_type(run_component, insitupy_insitudata, tmp_path):
    output = tmp_path / "dge_cell_type"
    run_component(
        [
            "--input",
            str(insitupy_insitudata),
            "--target_annotation_tuple",
            "test_annot;tumor",
            "--target_cell_type_tuple",
            "cell_type;A",
            "--ref_cell_type_tuple",
            "rest",
            "--output",
            str(output),
        ]
    )
    assert_outputs(output)


def test_insitudata_ref_project(
    run_component, insitupy_insitudata, insitupy_insitudata_ref, tmp_path
):
    output = tmp_path / "dge_ref_project"
    run_component(
        [
            "--input",
            str(insitupy_insitudata),
            "--target_annotation_tuple",
            "test_annot;tumor",
            "--ref",
            str(insitupy_insitudata_ref),
            "--output",
            str(output),
        ]
    )
    assert_outputs(output)


def test_insitudata_neighbors(run_component, insitupy_insitudata, tmp_path):
    output = tmp_path / "dge_neighbors"
    run_component(
        [
            "--input",
            str(insitupy_insitudata),
            "--target_annotation_tuple",
            "test_annot;tumor",
            "--ref_annotation_tuple",
            "test_annot;healthy",
            "--consider_neighbors",
            "--output",
            str(output),
        ]
    )
    assert_outputs(
        output,
        ("main.csv", "config.toml", "target_neighborhood.csv", "ref_neighborhood.csv"),
    )


def test_experiment(run_component, insitupy_experiment, tmp_path):
    output = tmp_path / "dge_experiment"
    run_component(
        [
            "--input",
            str(insitupy_experiment),
            "--target_id",
            "0",
            "--ref_id",
            "1",
            "--target_annotation_tuple",
            "test_annot;tumor",
            "--output",
            str(output),
        ]
    )
    assert_outputs(output)


def test_experiment_rest(run_component, insitupy_experiment, tmp_path):
    output = tmp_path / "dge_experiment_rest"
    run_component(
        [
            "--input",
            str(insitupy_experiment),
            "--target_id",
            "0",
            "--ref_id",
            "rest",
            "--target_annotation_tuple",
            "test_annot;tumor",
            "--output",
            str(output),
        ]
    )
    assert_outputs(output)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-W", "ignore::pytest.PytestAssertRewriteWarning"]))
