import os
import subprocess
import sys

import anndata as ad
import pandas as pd
import pytest


def _resource(name):
    return meta["resources_dir"] + "/" + name


def test_simple_execution(run_component, tmp_path):
    output = tmp_path / "output.h5ad"
    output_cell_fractions = tmp_path / "output_cell_fractions.csv"
    output_correlation = tmp_path / "output_correlation.csv"

    run_component(
        [
            "--bulk_data",
            _resource("input_bulk.tsv"),
            "--input_sc",
            _resource("input_sc.h5ad"),
            "--gene_id_mapping",
            _resource("pair_GRCm39.tsv"),
            "--celltype_key",
            "clusters",
            "--output",
            output,
            "--output_cell_fractions",
            output_cell_fractions,
            "--output_correlation",
            output_correlation,
            "--top_marker_num",
            "20",
            "--batch_size",
            "8",
            "--hidden_size",
            "8",
            "--epoch_num",
            "2",
            "--patience",
            "1",
            "--leiden_size",
            "2",
            "--n_comps",
            "10",
        ]
    )

    assert os.path.exists(output), "output h5ad was not created"
    assert os.path.exists(output_cell_fractions), "cell fraction csv was not created"
    assert os.path.exists(output_correlation), "correlation csv was not created"

    generated = ad.read_h5ad(output)
    assert generated.n_obs > 0, "Output AnnData should contain generated cells"
    assert "clusters" in generated.obs, "Cell type labels missing from output"

    cell_fractions = pd.read_csv(output_cell_fractions, index_col=0)
    assert cell_fractions.shape[0] > 0, "Cell fraction table should not be empty"

    correlation = pd.read_csv(output_correlation, index_col=0)
    assert correlation.shape[0] > 0 and correlation.shape[1] > 0, (
        "Correlation table should not be empty"
    )


def test_missing_celltype_key(run_component, tmp_path):
    output = tmp_path / "output.h5ad"

    with pytest.raises(subprocess.CalledProcessError):
        run_component(
            [
                "--bulk_data",
                _resource("input_bulk.tsv"),
                "--input_sc",
                _resource("input_sc.h5ad"),
                "--gene_id_mapping",
                _resource("pair_GRCm39.tsv"),
                "--celltype_key",
                "does_not_exist",
                "--output",
                output,
                "--top_marker_num",
                "20",
                "--batch_size",
                "8",
                "--hidden_size",
                "8",
                "--epoch_num",
                "2",
                "--patience",
                "1",
                "--leiden_size",
                "2",
                "--n_comps",
                "10",
            ]
        )

    assert not os.path.exists(output), "output h5ad should not be created on failure"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))