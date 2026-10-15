import sys
from pathlib import Path
import pytest
from insitupy import InSituData


## VIASH START
meta = {
    "name": "from_visium_to_insitupy",
    "resources_dir": "resources_test/visium",
}
## VIASH END

input_path = (
    Path(meta["resources_dir"]) / "Visium_FFPE_Human_Ovarian_Cancer_tiny_spaceranger"
)
pixel_size = 4.035


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


def test_basic_execution(run_component, tmp_path):
    output = tmp_path / "insitupy_project"
    run_component(
        [
            "--input",
            str(input_path),
            "--fullres_pixel_size",
            str(pixel_size),
            "--output",
            str(output),
        ]
    )
    assert_insitupy_output(output)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
