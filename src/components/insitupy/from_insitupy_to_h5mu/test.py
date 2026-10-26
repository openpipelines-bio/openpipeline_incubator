import sys
from pathlib import Path

import h5py
import mudata as mu
import pytest
from insitupy import InSituExperiment
from insitupy.io import read_xenium

## VIASH START
meta = {
    "name": "from_insitupy_to_h5mu",
    "resources_dir": "resources_test",
}
## VIASH END

mu.set_options(pull_on_update=False)

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


def assert_h5mu_output(output: Path, mod: str = "rna"):
    assert output.is_file(), "Output h5mu was not created."
    mdata = mu.read_h5mu(output)
    assert mod in mdata.mod


def test_insitudata(run_component, insitupy_insitudata, tmp_path):
    output = tmp_path / "insitudata.h5mu"
    run_component(
        [
            "--input",
            str(insitupy_insitudata),
            "--output",
            str(output),
        ]
    )
    assert_h5mu_output(output)


def test_insituexperiment(run_component, insitupy_experiment, tmp_path):
    output = tmp_path / "insitudata.h5mu"
    run_component(
        [
            "--input",
            str(insitupy_experiment),
            "--output",
            str(output),
        ]
    )
    assert_h5mu_output(output)


def test_output_compression(run_component, insitupy_insitudata, tmp_path):
    output = tmp_path / "compressed.h5mu"
    run_component(
        [
            "--input",
            str(insitupy_insitudata),
            "--output",
            str(output),
            "--output_compression",
            "gzip",
        ]
    )
    assert_h5mu_output(output)
    with h5py.File(output, "r") as f:
        x = f["mod/rna/X"]
        x = x["data"] if isinstance(x, h5py.Group) else x
        assert x.compression == "gzip", "Output X is not gzip compressed."


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-W", "ignore::pytest.PytestAssertRewriteWarning"]))
