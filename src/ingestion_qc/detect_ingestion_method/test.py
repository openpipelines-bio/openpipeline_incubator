import pytest
import h5py
import os
import anndata as ad
import sys

## VIASH START
meta = {
    "resources_dir": "resources_test"
}
## VIASH END


def test_cellranger(run_component, tmp_path):
    output = tmp_path / "output_cellranger.h5mu"

    run_component(
        [
            "--input", meta["resources_dir"] + "/sample_one.qc.h5mu",
            "--output", output
        ]
    )

    assert os.path.exists(output), "Output file was not created"
    with h5py.File (output, "r") as out_file:
        uns = ad.experimental.read_elem(out_file["uns"])
        assert uns["ingestion_method"] == "cellranger_multi", "cellranger_multi not detected"


def test_xenium(run_component, tmp_path):
    output = tmp_path / "output_xenium.h5mu"

    run_component(
        [
            "--input", meta["resources_dir"] + "/xenium_tiny.qc.h5mu",
            "--output", output
        ]
    )

    assert os.path.exists(output), "Output file was not created"
    with h5py.File (output, "r") as out_file:
        uns = ad.experimental.read_elem(out_file["uns"])
        assert uns["ingestion_method"] == "xenium", "xenium not detected"


def test_cosmx(run_component, tmp_path):
    output = tmp_path / "output_cosmx.h5mu"

    run_component(
        [
            "--input", meta["resources_dir"] + "/Lung5_Rep2_tiny.qc.h5mu",
            "--output", output
        ]
    )

    assert os.path.exists(output), "Output file was not created"
    with h5py.File (output, "r") as out_file:
        uns = ad.experimental.read_elem(out_file["uns"])
        assert uns["ingestion_method"] == "cosmx", "cosmx not detected"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))