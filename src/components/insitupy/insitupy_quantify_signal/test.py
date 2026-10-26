import sys
from pathlib import Path

import pytest
from insitupy import InSituData, InSituExperiment
from insitupy.io import read_xenium

## VIASH START
meta = {
    "name": "insitupy_quantify_signal",
    "resources_dir": "resources_test",
}
## VIASH END
xenium_bundle = Path(meta["resources_dir"]) / "xenium_multicellseg_tiny"

# The DAPI image is the only image in the tiny Xenium bundle.
image_name = "nuclei"

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


def assert_quantified(xd: InSituData, column: str):
    obs = xd.cells.table.obs
    assert column in obs.columns, f"Column '{column}' not found in .obs."
    assert obs[column].notna().any(), f"Column '{column}' contains only NaN."


def test_insitudata(run_component, insitupy_insitudata, tmp_path):
    output = tmp_path / "quantified_insitudata"
    run_component(
        [
            "--input",
            str(insitupy_insitudata),
            "--image_name",
            image_name,
            "--output",
            str(output),
        ]
    )
    assert output.is_dir(), "Output project was not created."
    assert_quantified(InSituData.read(output), "nuclei_signal_cells_median")


def test_insitudata_parameters(run_component, insitupy_insitudata, tmp_path):
    output = tmp_path / "quantified_insitudata"
    run_component(
        [
            "--input",
            str(insitupy_insitudata),
            "--image_name",
            image_name,
            "--cells_compartment",
            "nuclei",
            "--method",
            "mean",
            "--downsample_factor",
            "2",
            "--output",
            str(output),
        ]
    )
    assert output.is_dir(), "Output project was not created."
    assert_quantified(InSituData.read(output), "nuclei_signal_nuclei_mean")


def test_experiment(run_component, insitupy_experiment, tmp_path):
    output = tmp_path / "quantified_experiment"
    run_component(
        [
            "--input",
            str(insitupy_experiment),
            "--image_name",
            image_name,
            "--output",
            str(output),
        ]
    )
    assert output.is_dir(), "Output experiment was not created."
    exp = InSituExperiment.read(output)
    exp.load_all()
    for _, xd in exp.iterdata():
        assert_quantified(xd, "nuclei_signal_cells_median")


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-W", "ignore::pytest.PytestAssertRewriteWarning"]))
