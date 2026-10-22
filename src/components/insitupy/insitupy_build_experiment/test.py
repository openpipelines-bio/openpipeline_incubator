import sys
from pathlib import Path
from insitupy import InSituExperiment
import pandas as pd

import pytest
from insitupy.io import read_xenium

## VIASH START
meta = {
    "name": "insitupy_build_experiment",
    "resources_dir": "resources_test",
}
## VIASH END

xenium_bundle = Path(meta["resources_dir"]) / "xenium_multicellseg_tiny"

crop_xlim = (250, 600)
crop_ylim = (300, 550)


def assert_insitupy_output(path: Path):
    assert path.is_dir(), "Output folder was not created."
    exp = InSituExperiment.read(path)
    exp.load_all()


metadata = pd.DataFrame(
    {
        "sample_id": ["sample_1", "sample_2"],
        "therapy": ["drug", "control"],
        "age": [23, 26],
    }
)


# There are no insitupy resources on s3 yet, so this fixture converts the Xenium
# bundle itself into two samples.
# The second sample is a crop of the first.
@pytest.fixture(scope="module")
def insitupy_projects(tmp_path_factory):
    tmp_dir = tmp_path_factory.mktemp("insitupy_projects")

    xd = read_xenium(xenium_bundle, verbose=False)

    sample_1 = tmp_dir / "sample_1"
    xd.saveas(sample_1, verbose=False)

    sample_2 = tmp_dir / "sample_2"
    xd.crop(xlim=crop_xlim, ylim=crop_ylim).saveas(sample_2, verbose=False)

    return [sample_1, sample_2]


def test_basic_execution(run_component, insitupy_projects, tmp_path):
    output = tmp_path / "insitupy_experiment"
    run_component(
        [
            "--input",
            str(insitupy_projects[0]),
            "--input",
            str(insitupy_projects[1]),
            "--sample_id",
            "sample_1",
            "--sample_id",
            "sample_2",
            "--output",
            str(output),
        ]
    )
    assert_insitupy_output(output)


def test_metadata(run_component, insitupy_projects, tmp_path):
    output = tmp_path / "insitupy_experiment"
    metadata.to_csv(tmp_path / "metadata.csv", index=False)
    run_component(
        [
            "--input",
            str(insitupy_projects[0]),
            "--input",
            str(insitupy_projects[1]),
            "--sample_id",
            "sample_1",
            "--sample_id",
            "sample_2",
            "--metadata",
            str(tmp_path / "metadata.csv"),
            "--output",
            str(output),
        ]
    )
    assert_insitupy_output(output)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
