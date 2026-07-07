import sys
import pytest
import pandas as pd
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
        assert modality in output_mudata.mod, (
            f"Modality '{modality}' should be present in the output."
        )

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
        assert (report[ALL_COLUMN] <= report[column]).all(), (
            f"'{ALL_COLUMN}' should never exceed '{column}' for any sample."
        )

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
        assert (survived[median_col] > 0).all(), (
            f"Median total counts should be positive for surviving samples in '{median_col}'."
        )
        assert (survived[max_col] >= survived[median_col]).all(), (
            f"Max total counts should be >= median in '{max_col}'."
        )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "--import-mode=importlib"]))
