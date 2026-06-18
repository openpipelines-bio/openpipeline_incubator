import sys
import pytest
from pathlib import Path
import mudata as md
import subprocess
import re

## VIASH START
meta = {
    "executable": "./target/executable/qc/calculate_qc_metrics_scanpy/calculate_qc_metrics_scanpy",
    "resources_dir": "./resources_test/pbmc_1k_protein_v3/",
    "config": "./src/qc/calculate_qc_metrics_scanpy/config.vsh.yaml",
}
## VIASH END


@pytest.fixture
def input_path():
    return Path(
        f"{meta['resources_dir']}/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"
    )


@pytest.fixture
def input_with_mito_ribo(tmp_path, input_path):
    """Input with boolean .var columns flagging mitochondrial and ribosomal genes."""
    output_path = tmp_path / "input_with_mito_ribo.h5mu"
    mdata = md.read_h5mu(input_path)
    rna = mdata.mod["rna"]
    rna.var["mitochondrial"] = rna.var_names.str.match("^[mM][tT]-")
    rna.var["ribosomal"] = rna.var_names.str.match("^RP[SL]")
    mdata.write_h5mu(output_path)
    return output_path


def test_default_metrics(run_component, tmp_path, input_path):
    output_path = tmp_path / "output.h5mu"

    run_component(
        [
            "--input",
            input_path,
            "--modality",
            "rna",
            "--output",
            output_path,
        ]
    )

    assert output_path.exists(), "Output file does not exist"
    adata_out = md.read_h5ad(output_path, mod="rna")

    for col in ("total_counts", "n_genes_by_counts"):
        assert col in adata_out.obs.columns, f"Expected .obs column '{col}'"
    # log1p enabled by default
    assert "log1p_total_counts" in adata_out.obs.columns
    # var metrics
    assert "n_cells_by_counts" in adata_out.var.columns
    assert "total_counts" in adata_out.var.columns


def test_no_log1p(run_component, tmp_path, input_path):
    output_path = tmp_path / "output.h5mu"

    run_component(
        [
            "--input",
            input_path,
            "--modality",
            "rna",
            "--log1p",
            "false",
            "--output",
            output_path,
        ]
    )

    adata_out = md.read_h5ad(output_path, mod="rna")
    assert "total_counts" in adata_out.obs.columns
    assert "log1p_total_counts" not in adata_out.obs.columns


def test_top_n_vars(run_component, tmp_path, input_path):
    output_path = tmp_path / "output.h5mu"

    run_component(
        [
            "--input",
            input_path,
            "--modality",
            "rna",
            "--top_n_vars",
            "50;100",
            "--output",
            output_path,
        ]
    )

    adata_out = md.read_h5ad(output_path, mod="rna")
    assert "pct_counts_in_top_50_genes" in adata_out.obs.columns
    assert "pct_counts_in_top_100_genes" in adata_out.obs.columns


def test_qc_vars(run_component, tmp_path, input_with_mito_ribo):
    output_path = tmp_path / "output.h5mu"

    run_component(
        [
            "--input",
            input_with_mito_ribo,
            "--modality",
            "rna",
            "--qc_vars",
            "mitochondrial",
            "--output",
            output_path,
        ]
    )

    adata_out = md.read_h5ad(output_path, mod="rna")
    assert "pct_counts_mitochondrial" in adata_out.obs.columns
    assert "total_counts_mitochondrial" in adata_out.obs.columns
    pct = adata_out.obs["pct_counts_mitochondrial"]
    assert ((pct >= 0) & (pct <= 100)).all(), "Percentages must be in [0, 100]"


def test_multiple_qc_vars(run_component, tmp_path, input_with_mito_ribo):
    output_path = tmp_path / "output.h5mu"

    run_component(
        [
            "--input",
            input_with_mito_ribo,
            "--modality",
            "rna",
            "--qc_vars",
            "mitochondrial;ribosomal",
            "--output",
            output_path,
        ]
    )

    adata_out = md.read_h5ad(output_path, mod="rna")
    for qc_var in ("mitochondrial", "ribosomal"):
        for prefix in ("pct_counts", "total_counts"):
            col = f"{prefix}_{qc_var}"
            assert col in adata_out.obs.columns, f"Expected .obs column '{col}'"
        pct = adata_out.obs[f"pct_counts_{qc_var}"]
        assert ((pct >= 0) & (pct <= 100)).all(), "Percentages must be in [0, 100]"


def test_missing_qc_var_raises(run_component, tmp_path, input_path):
    output_path = tmp_path / "output.h5mu"

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                input_path,
                "--modality",
                "rna",
                "--qc_vars",
                "does_not_exist",
                "--output",
                output_path,
            ]
        )

    assert re.search(
        r"does_not_exist.*not found in \.var", err.value.stdout.decode("utf-8")
    ), f"Expected error message not found: {err.value.stdout.decode('utf-8')}"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
