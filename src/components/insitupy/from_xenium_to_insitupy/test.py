import sys
from pathlib import Path

import pytest
from insitupy import InSituData
from insitupy.io.data import CELLSEG_NAMES


## VIASH START
meta = {
    "name": "from_xenium_to_insitupy",
    "resources_dir": "resources_test/xenium_multichannel",
}
## VIASH END

input_path = Path(meta["resources_dir"]) / "xenium_multicellseg_tiny"


def assert_insitupy_output(path: Path):
    assert path.is_dir(), "Output folder was not created."
    try:
        xd = InSituData.read(path)
    except ValueError as e:
        pytest.fail(f"{path.name} could not be read as an InSituData project: {e}")

    try:
        xd.load_all()
    except RuntimeError as e:
        pytest.fail(f"Not all modalities of {xd} could be loaded: {e}")

    return xd


def run_and_read(run_component, output: Path, extra_args: list):
    run_component(["--input", str(input_path), "--output", str(output), *extra_args])
    return assert_insitupy_output(output)


# The images insitupy treats as cell segmentation stains
CELLSEG_IMAGES = set(CELLSEG_NAMES)


def test_basic_execution(run_component, tmp_path):
    output = tmp_path / "insitupy_project"
    run_component(
        [
            "--input",
            str(input_path),
            "--output",
            str(output),
        ]
    )
    assert_insitupy_output(output)


def test_spatial_data_backend_basic(run_component, tmp_path):
    output = tmp_path / "insitupy_project"
    run_component(
        [
            "--input",
            str(input_path),
            "--backend",
            "spatialdata",
            "--output",
            str(output),
        ]
    )
    assert_insitupy_output(output)


def test_spatial_data_backend_extended(run_component, tmp_path):
    output = tmp_path / "insitupy_project"
    run_component(
        [
            "--input",
            str(input_path),
            "--backend",
            "spatialdata",
            "--slide_id",
            "slide",
            "--sample_id",
            "sample",
            "--output",
            str(output),
        ]
    )
    assert_insitupy_output(output)

    xd = InSituData.read(output)
    assert xd.slide_id == "slide", (
        "xd should contain a slide id when --slide_id is enabled"
    )
    assert xd.sample_id == "sample", (
        "xd should contain a sample id when --sample_id is enabled"
    )


def test_load_cell_segmentation_images(run_component, tmp_path):
    xd_default = run_and_read(run_component, tmp_path / "default", [])
    xd_cellseg = run_and_read(
        run_component, tmp_path / "cellseg", ["--load_cell_segmentation_images"]
    )

    assert not CELLSEG_IMAGES & set(xd_default.images.names), (
        "Cell segmentation images should not be loaded by default."
    )
    assert CELLSEG_IMAGES <= set(xd_cellseg.images.names), (
        "--load_cell_segmentation_images should add the cell segmentation images."
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
