import sys
import pytest
import numpy as np
import pandas as pd
from scipy.sparse import issparse
from mudata import read_h5mu

## VIASH START
par = {
    "input": "output.h5mu",
    "og_input": "input.h5mu",
    "cell_count_report": "cell_counts.tsv",
    "rna_modality": "rna",
    "prot_modality": "prot",
    "sample_id_column": "sample_id",
    "expected_sample_id": None,
    "expect_scrublet": True,
}
## VIASH END

# The cell-count report column names match the defaults of the cell_count_report component, which
# process_single_sample does not override.
RNA_COLUMN = "cell_count_after_rna_filter"
SCRUBLET_COLUMN = "cell_count_after_scrublet_filter"
PROT_COLUMN = "cell_count_after_protein_filter"
ALL_COLUMN = "cell_count_after_all_filter"

# Post-filter QC aggregations the report adds, matching the cell_count_report defaults.
QC_AGGREGATIONS = ("median", "mean", "max")
RNA_QC_METRICS = ("total_counts", "pct_counts_mt", "n_genes_by_counts")
PROT_QC_METRICS = ("total_counts", "n_genes_by_counts")


def _sum(matrix, axis):
    """Sum a (possibly sparse) matrix along an axis, returned as a flat float array."""
    return np.asarray(matrix.sum(axis=axis)).ravel().astype(float)


def _nonzero_count(matrix, axis):
    """Count non-zero entries of a (possibly sparse) matrix along an axis, as a flat int array."""
    nonzero = (matrix > 0) if issparse(matrix) else (np.asarray(matrix) > 0)
    return np.asarray(nonzero.sum(axis=axis)).ravel().astype(int)


def test_run():
    og_mudata = read_h5mu(par["og_input"])
    output_mudata = read_h5mu(par["input"])

    rna_modality = par["rna_modality"]
    prot_modality = par["prot_modality"]

    # Every modality present in the input must still be present in the output.
    assert set(og_mudata.mod.keys()) == set(output_mudata.mod.keys()), (
        f"Modalities should be unchanged. Input: {list(og_mudata.mod.keys())}, "
        f"output: {list(output_mudata.mod.keys())}"
    )
    for modality in (rna_modality, prot_modality):
        assert (
            modality in output_mudata.mod
        ), f"Modality '{modality}' should be present in the output."

    og_rna = og_mudata.mod[rna_modality]
    output_rna = output_mudata.mod[rna_modality]
    output_prot = output_mudata.mod[prot_modality]

    # add_id stamped the workflow event id onto a .obs column before the modalities were split,
    # so every modality in the output must carry that column with the expected sample id.
    if par["expected_sample_id"] is not None:
        for modality, mod_data in (
            (rna_modality, output_rna),
            (prot_modality, output_prot),
        ):
            assert par["sample_id_column"] in mod_data.obs.columns, (
                f"add_id should have written the sample id column "
                f"'{par['sample_id_column']}' to .obs of modality '{modality}'. "
                f"Found: {mod_data.obs.columns.to_list()}"
            )
            assert (
                mod_data.obs[par["sample_id_column"]] == par["expected_sample_id"]
            ).all(), (
                f"All cells in modality '{modality}' should carry sample id "
                f"'{par['expected_sample_id']}' in .obs['{par['sample_id_column']}']."
            )

    # Unlike the flag-only QC sub-workflows, process_single_sample subsets: the output must not
    # contain more cells or genes than the input.
    assert output_rna.n_obs <= og_rna.n_obs, (
        f"Output should not contain more cells than the input. Input: {og_rna.n_obs} cells, "
        f"output: {output_rna.n_obs} cells."
    )
    assert output_rna.n_vars <= og_rna.n_vars, (
        f"Rare RNA genes should be removed, so the output should not contain more genes than "
        f"the input. Input: {og_rna.n_vars} genes, output: {output_rna.n_vars} genes."
    )

    # The modalities were intersected, so RNA and protein must share the exact same set of cells.
    pd.testing.assert_index_equal(
        output_rna.obs_names, output_prot.obs_names, check_order=False
    )

    # The cell-count report: one row per sample, the expected stage columns, and counts that are
    # consistent with one another.
    report = pd.read_csv(par["cell_count_report"], sep="\t")

    assert par["sample_id_column"] in report.columns, (
        f"Cell-count report should contain the sample id column "
        f"'{par['sample_id_column']}'. Found: {report.columns.to_list()}"
    )

    # add_id wrote the workflow event id into the sample-id .obs column, so the report groups by
    # that column and the event id must appear as a sample id.
    if par["expected_sample_id"] is not None:
        assert par["expected_sample_id"] in set(report[par["sample_id_column"]]), (
            f"Report should be grouped by the workflow event id "
            f"'{par['expected_sample_id']}'. Found sample ids: "
            f"{report[par['sample_id_column']].to_list()}"
        )

    expected_columns = [RNA_COLUMN, PROT_COLUMN, ALL_COLUMN]
    if par["expect_scrublet"]:
        expected_columns.append(SCRUBLET_COLUMN)
    else:
        assert SCRUBLET_COLUMN not in report.columns, (
            f"Scrublet was skipped, so the report should not contain '{SCRUBLET_COLUMN}'. "
            f"Found: {report.columns.to_list()}"
        )
    for column in expected_columns:
        assert column in report.columns, (
            f"Cell-count report should contain the stage column '{column}'. "
            f"Found: {report.columns.to_list()}"
        )

    assert len(report) >= 1, "Cell-count report should contain at least one sample."

    # The combined "all" count is the intersection of every stage, so it can never exceed any
    # individual stage count.
    per_stage_columns = [RNA_COLUMN, PROT_COLUMN]
    if par["expect_scrublet"]:
        per_stage_columns.append(SCRUBLET_COLUMN)
    for column in per_stage_columns:
        assert (
            report[ALL_COLUMN] <= report[column]
        ).all(), f"'{ALL_COLUMN}' should never exceed '{column}' for any sample."

    # The report is computed before subsetting; the final output is exactly the cells passing
    # every filter in every modality. The two must agree.
    assert report[ALL_COLUMN].sum() == output_rna.n_obs, (
        f"Total cells surviving all filters in the report ({report[ALL_COLUMN].sum()}) should "
        f"equal the number of cells in the subset output ({output_rna.n_obs})."
    )

    # The report carries post-filter QC aggregations for every (modality, aggregation, metric).
    expected_qc_columns = [
        f"{par['rna_modality']}_{agg}_{metric}_after_all_filter"
        for metric in RNA_QC_METRICS
        for agg in QC_AGGREGATIONS
    ] + [
        f"{par['prot_modality']}_{agg}_{metric}_after_all_filter"
        for metric in PROT_QC_METRICS
        for agg in QC_AGGREGATIONS
    ]
    for column in expected_qc_columns:
        assert column in report.columns, (
            f"Cell-count report should contain the QC column '{column}'. "
            f"Found: {report.columns.to_list()}"
        )

    # The QC aggregations are computed over the surviving cells, so for every sample with at
    # least one surviving cell the mean total counts must be positive and the max must not be
    # smaller than the median.
    survived = report[report[ALL_COLUMN] > 0]
    for modality in (par["rna_modality"], par["prot_modality"]):
        median_col = f"{modality}_median_total_counts_after_all_filter"
        max_col = f"{modality}_max_total_counts_after_all_filter"
        assert (
            survived[median_col] > 0
        ).all(), f"Median total counts should be positive for surviving samples in '{median_col}'."
        assert (
            survived[max_col] >= survived[median_col]
        ).all(), f"Max total counts should be >= median in '{max_col}'."


def test_post_filter_qc_metrics_recomputed():
    """The post-filter calculate_qc_metrics steps (rna_/prot_post_filter_qc_metrics in main.nf)
    overwrite the QC slots so the output carries metrics reflecting the *filtered* data, not the
    pre-filter values the sub-workflows wrote. This recomputes the metrics independently from the
    output's own count matrix and asserts the stored .obs/.var values match. Because RNA gene
    filtering dropped genes and cell filtering dropped cells, stale pre-filter values would not
    match these recomputations — so a pass proves the second calculate_qc_metrics ran on the
    subset data.
    """
    output_mudata = read_h5mu(par["input"])
    output_rna = output_mudata.mod[par["rna_modality"]]
    output_prot = output_mudata.mod[par["prot_modality"]]

    # The workflow leaves rna_layer/prot_layer unset in this test, so metrics were computed from
    # .X. Compare against .X here.
    for modality_name, mod_data, metrics in (
        (par["rna_modality"], output_rna, RNA_QC_METRICS),
        (par["prot_modality"], output_prot, PROT_QC_METRICS),
    ):
        for metric in metrics:
            assert metric in mod_data.obs.columns, (
                f"Post-filter QC metric '{metric}' should be present in the '{modality_name}' "
                f".obs. Found: {mod_data.obs.columns.to_list()}"
            )

        matrix = mod_data.X
        stored_total = mod_data.obs["total_counts"].to_numpy(dtype=float)
        expected_total = _sum(matrix, axis=1)
        np.testing.assert_allclose(
            stored_total,
            expected_total,
            rtol=1e-4,
            err_msg=(
                f"'{modality_name}' .obs['total_counts'] should equal the per-cell sum of the "
                f"filtered matrix; a mismatch means the metrics were not recomputed after "
                f"filtering."
            ),
        )

        stored_n_genes = mod_data.obs["n_genes_by_counts"].to_numpy(dtype=int)
        expected_n_genes = _nonzero_count(matrix, axis=1)
        np.testing.assert_array_equal(
            stored_n_genes,
            expected_n_genes,
            err_msg=(
                f"'{modality_name}' .obs['n_genes_by_counts'] should equal the per-cell count of "
                f"detected features in the filtered matrix."
            ),
        )

    # RNA mitochondrial percentage: recompute from the preserved `mt` .var flag. This exercises the
    # `qc_vars: [mt]` path of the post-filter recompute and confirms the flag survived filtering.
    assert "mt" in output_rna.var.columns, (
        f"The 'mt' .var flag should be preserved so pct_counts_mt can be recomputed. "
        f"Found: {output_rna.var.columns.to_list()}"
    )
    mt_indices = np.where(output_rna.var["mt"].to_numpy(dtype=bool))[0]
    mt_counts = _sum(output_rna.X[:, mt_indices], axis=1)
    total_counts = _sum(output_rna.X, axis=1)
    expected_pct_mt = 100.0 * mt_counts / total_counts
    np.testing.assert_allclose(
        output_rna.obs["pct_counts_mt"].to_numpy(dtype=float),
        expected_pct_mt,
        rtol=1e-4,
        atol=1e-6,
        err_msg=(
            "RNA .obs['pct_counts_mt'] should equal the mitochondrial fraction of the filtered "
            "matrix; a mismatch means it was not recomputed after gene/cell filtering."
        ),
    )

    # scanpy also refreshes the per-gene .var metrics. n_cells_by_counts must reflect the surviving
    # cell set, which changed after cell filtering and intersection.
    if "n_cells_by_counts" in output_rna.var.columns:
        np.testing.assert_array_equal(
            output_rna.var["n_cells_by_counts"].to_numpy(dtype=int),
            _nonzero_count(output_rna.X, axis=0),
            err_msg=(
                "RNA .var['n_cells_by_counts'] should equal the per-gene count of cells with "
                "non-zero expression in the filtered matrix."
            ),
        )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "--import-mode=importlib"]))
