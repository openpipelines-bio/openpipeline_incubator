import sys
from pathlib import Path

import pytest
from insitupy import InSituData


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


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))