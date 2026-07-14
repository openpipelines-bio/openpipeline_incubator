import os
import subprocess
import sys

import anndata as ad
import mudata as mu
import pytest


def _prepare_sc_input(tmp_path):
    mdata = mu.read_h5mu(meta["resources_dir"] + "/TS_Blood_filtered.h5mu")
    adata = mdata.mod["rna"].copy()
    # Align gene identifiers with the unversioned Ensembl IDs used by the Visium panel.
    adata.var_names = adata.var["ensemblid"].str.split(".").str[0]
    adata.var.index.name = None

    path = tmp_path / "input_sc.h5ad"
    adata.write_h5ad(path)
    return path


def _prepare_spatial_input(tmp_path):
    mdata = mu.read_h5mu(meta["resources_dir"] + "/Visium_FFPE_Human_Ovarian_Cancer_tiny.h5mu")
    adata = mdata.mod["rna"].copy()

    path = tmp_path / "input_spatial.h5ad"
    adata.write_h5ad(path)
    return path


def test_simple_execution(run_component, tmp_path):
    input_sc = _prepare_sc_input(tmp_path)
    input_spatial = _prepare_spatial_input(tmp_path)
    output = tmp_path / "output.h5ad"

    run_component(
        [
            "--input_sc",
            input_sc,
            "--input_spatial",
            input_spatial,
            "--celltype_key",
            "cell_type",
            "--obsm_spatial_coordinates",
            "spatial",
            "--output",
            output,
            "--spot_num",
            "50",
            "--cell_num",
            "5",
            "--k",
            "5",
            "--num_epochs",
            "2",
            "--batch_size",
            "8",
            "--predicted_size",
            "8",
        ]
    )
    assert os.path.exists(output), "output h5ad was not created"

    sp_adata = ad.read_h5ad(output)
    assert "spatial" in sp_adata.obsm, "Predicted coordinates missing from .obsm['spatial']"
    assert sp_adata.obsm["spatial"].shape[1] == 2, "Predicted coordinates should have 2 dimensions"
    assert sp_adata.n_obs > 0, "Output AnnData should contain predicted cells"


def test_missing_celltype_key(run_component, tmp_path):
    input_sc = _prepare_sc_input(tmp_path)
    input_spatial = _prepare_spatial_input(tmp_path)
    output = tmp_path / "output.h5ad"

    with pytest.raises(subprocess.CalledProcessError):
        run_component(
            [
                "--input_sc",
                input_sc,
                "--input_spatial",
                input_spatial,
                "--celltype_key",
                "does_not_exist",
                "--obsm_spatial_coordinates",
                "spatial",
                "--output",
                output,
                "--spot_num",
                "50",
                "--cell_num",
                "5",
                "--k",
                "5",
                "--num_epochs",
                "2",
                "--batch_size",
                "8",
                "--predicted_size",
                "8",
            ]
        )

    assert not os.path.exists(output), "output h5ad should not be created on failure"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))