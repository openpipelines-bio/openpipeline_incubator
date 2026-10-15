import pytest
import uuid
import numpy as np
import pandas as pd
from anndata import AnnData
from mudata import MuData

## VIASH START
meta = {
    "name": "cell_count_report",
    "resources_dir": "resources_test/",
    "executable": "target/executable/metadata/cell_count_report/cell_count_report",
    "config": "src/metadata/cell_count_report/config.vsh.yaml",
}
## VIASH END


@pytest.fixture
def write_temp_h5mu(tmp_path):
    def wrapper(test_h5mu):
        test_h5mu_path = tmp_path / f"{str(uuid.uuid4())}.h5mu"
        test_h5mu.write_h5mu(test_h5mu_path)
        return test_h5mu_path

    return wrapper


@pytest.fixture
def h5mu():
    # Two samples, 4 cells each. Boolean keep-columns (True = keep) chosen so the
    # per-stage and combined counts are easy to verify by hand. Per-cell QC metric columns
    # carry distinctive values so the post-filter aggregations can also be checked by hand.
    obs_names = [f"cell_{i}" for i in range(8)]
    sample_id = ["s1"] * 4 + ["s2"] * 4

    # RNA filter masks (AND of these two columns gives the RNA keep set).
    filter_counts_rna = [True, True, True, False, True, True, False, False]
    filter_quantile_rna = [True, True, False, True, True, False, True, True]
    # rna_keep:                T,    T,    F,     F,    T,    F,     F,     F
    filter_scrublet = [True, False, True, True, True, True, True, False]
    # prot masks
    filter_counts_prot = [True, True, True, True, True, True, True, True]
    filter_quantile_prot = [True, True, True, True, False, True, True, True]
    # all_keep = rna & scrublet & prot:
    #   s1: [T, F, F, F] -> only cell_0 survives
    #   s2: [F, F, F, F] -> no cell survives

    # Per-cell QC metrics (already present in .obs from the upstream QC sub-workflows).
    rna_total_counts = [1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000]
    rna_pct_counts_mt = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
    rna_n_genes_by_counts = [10, 20, 30, 40, 50, 60, 70, 80]
    prot_total_counts = [100, 200, 300, 400, 500, 600, 700, 800]
    prot_n_genes_by_counts = [1, 2, 3, 4, 5, 6, 7, 8]

    rna_obs = pd.DataFrame(
        {
            "sample_id": sample_id,
            "filter_counts_rna": filter_counts_rna,
            "filter_quantile_rna": filter_quantile_rna,
            "filter_scrublet": filter_scrublet,
            "total_counts": rna_total_counts,
            "pct_counts_mt": rna_pct_counts_mt,
            "n_genes_by_counts": rna_n_genes_by_counts,
        },
        index=obs_names,
    )
    prot_obs = pd.DataFrame(
        {
            "sample_id": sample_id,
            "filter_counts_prot": filter_counts_prot,
            "filter_quantile_prot": filter_quantile_prot,
            "total_counts": prot_total_counts,
            "n_genes_by_counts": prot_n_genes_by_counts,
        },
        index=obs_names,
    )
    rna = AnnData(np.ones((8, 3)), obs=rna_obs)
    prot = AnnData(np.ones((8, 2)), obs=prot_obs)
    return MuData({"rna": rna, "prot": prot})


def test_cell_count_report(run_component, h5mu, write_temp_h5mu, tmp_path):
    output = tmp_path / "cell_counts.tsv"
    run_component(
        [
            "--input",
            write_temp_h5mu(h5mu),
            "--modality",
            "rna",
            "--prot_modality",
            "prot",
            "--sample_id_column",
            "sample_id",
            "--rna_filter_columns",
            "filter_counts_rna;filter_quantile_rna",
            "--scrublet_filter_column",
            "filter_scrublet",
            "--prot_filter_columns",
            "filter_counts_prot;filter_quantile_prot",
            "--output",
            output,
        ]
    )
    assert output.is_file(), "Output TSV must have been created."
    report = pd.read_csv(output, sep="\t").set_index("sample_id")

    expected_columns = {
        "cell_count_after_rna_filter",
        "cell_count_after_scrublet_filter",
        "cell_count_after_protein_filter",
        "cell_count_after_all_filter",
    }
    assert expected_columns.issubset(set(report.columns))
    assert set(report.index) == {"s1", "s2"}

    # s1: rna_keep = [T,T,F,F] -> 2 ; scrublet = [T,F,T,T] -> 3 ; prot = [T,T,T,T] -> 4
    #     all = rna & scrublet & prot = [T,F,F,F] -> 1
    assert report.loc["s1", "cell_count_after_rna_filter"] == 2
    assert report.loc["s1", "cell_count_after_scrublet_filter"] == 3
    assert report.loc["s1", "cell_count_after_protein_filter"] == 4
    assert report.loc["s1", "cell_count_after_all_filter"] == 1

    # s2: rna_keep = [T,F,F,F] -> 1 ; scrublet = [T,T,T,F] -> 3 ; prot = [F,T,T,T] -> 3
    #     all = [F,F,F,F]... cell_4: rna T, scrublet T, prot F -> F ; so all -> 0
    assert report.loc["s2", "cell_count_after_rna_filter"] == 1
    assert report.loc["s2", "cell_count_after_scrublet_filter"] == 3
    assert report.loc["s2", "cell_count_after_protein_filter"] == 3
    assert report.loc["s2", "cell_count_after_all_filter"] == 0

    # The combined count never exceeds any single-stage count.
    for sample in report.index:
        for col in (
            "cell_count_after_rna_filter",
            "cell_count_after_scrublet_filter",
            "cell_count_after_protein_filter",
        ):
            assert (
                report.loc[sample, "cell_count_after_all_filter"]
                <= report.loc[sample, col]
            )


def test_cell_count_report_qc_metrics(run_component, h5mu, write_temp_h5mu, tmp_path):
    # The post-filter QC aggregations are computed over the cells surviving every filter
    # (the all_keep set), using the per-cell QC columns already present in .obs.
    output = tmp_path / "cell_counts.tsv"
    run_component(
        [
            "--input",
            write_temp_h5mu(h5mu),
            "--modality",
            "rna",
            "--prot_modality",
            "prot",
            "--sample_id_column",
            "sample_id",
            "--rna_filter_columns",
            "filter_counts_rna;filter_quantile_rna",
            "--scrublet_filter_column",
            "filter_scrublet",
            "--prot_filter_columns",
            "filter_counts_prot;filter_quantile_prot",
            "--output",
            output,
        ]
    )
    report = pd.read_csv(output, sep="\t").set_index("sample_id")

    # One column per (modality, aggregation, metric) combination, suffixed with the filter stage.
    expected_qc_columns = set()
    for fun in ("median", "mean", "max"):
        for metric in ("total_counts", "pct_counts_mt", "n_genes_by_counts"):
            expected_qc_columns.add(f"rna_{fun}_{metric}_after_all_filter")
        for metric in ("total_counts", "n_genes_by_counts"):
            expected_qc_columns.add(f"prot_{fun}_{metric}_after_all_filter")
    assert expected_qc_columns.issubset(set(report.columns)), (
        f"Missing QC columns: {expected_qc_columns - set(report.columns)}"
    )

    # s1: only cell_0 survives, so median == mean == max == cell_0's value for every metric.
    for fun in ("median", "mean", "max"):
        assert report.loc["s1", f"rna_{fun}_total_counts_after_all_filter"] == 1000.0
        assert report.loc["s1", f"rna_{fun}_pct_counts_mt_after_all_filter"] == 1.0
        assert report.loc["s1", f"rna_{fun}_n_genes_by_counts_after_all_filter"] == 10.0
        assert report.loc["s1", f"prot_{fun}_total_counts_after_all_filter"] == 100.0
        assert report.loc["s1", f"prot_{fun}_n_genes_by_counts_after_all_filter"] == 1.0

    # s2: no cell survives, so the QC aggregations are undefined (NaN) but the sample is
    # still reported (its per-stage counts are meaningful).
    for fun in ("median", "mean", "max"):
        assert pd.isna(report.loc["s2", f"rna_{fun}_total_counts_after_all_filter"])
        assert pd.isna(report.loc["s2", f"prot_{fun}_total_counts_after_all_filter"])


def test_cell_count_report_multiple_inputs(
    run_component, h5mu, write_temp_h5mu, tmp_path
):
    # Splitting the two samples into separate input files must yield the same report as
    # passing them in a single file: the report combines the per-cell masks across inputs.
    s1 = h5mu[h5mu["rna"].obs["sample_id"] == "s1", :].copy()
    s2 = h5mu[h5mu["rna"].obs["sample_id"] == "s2", :].copy()
    output = tmp_path / "cell_counts.tsv"
    run_component(
        [
            "--input",
            f"{write_temp_h5mu(s1)};{write_temp_h5mu(s2)}",
            "--modality",
            "rna",
            "--prot_modality",
            "prot",
            "--sample_id_column",
            "sample_id",
            "--rna_filter_columns",
            "filter_counts_rna;filter_quantile_rna",
            "--scrublet_filter_column",
            "filter_scrublet",
            "--prot_filter_columns",
            "filter_counts_prot;filter_quantile_prot",
            "--output",
            output,
        ]
    )
    assert output.is_file(), "Output TSV must have been created."
    report = pd.read_csv(output, sep="\t").set_index("sample_id")

    assert set(report.index) == {"s1", "s2"}
    # Same expected counts as the single-file case.
    assert report.loc["s1", "cell_count_after_rna_filter"] == 2
    assert report.loc["s1", "cell_count_after_scrublet_filter"] == 3
    assert report.loc["s1", "cell_count_after_protein_filter"] == 4
    assert report.loc["s1", "cell_count_after_all_filter"] == 1
    assert report.loc["s2", "cell_count_after_rna_filter"] == 1
    assert report.loc["s2", "cell_count_after_scrublet_filter"] == 3
    assert report.loc["s2", "cell_count_after_protein_filter"] == 3
    assert report.loc["s2", "cell_count_after_all_filter"] == 0
    # The QC aggregations are also identical when the samples are split across files.
    assert report.loc["s1", "rna_max_total_counts_after_all_filter"] == 1000.0
    assert report.loc["s1", "prot_max_n_genes_by_counts_after_all_filter"] == 1.0


def test_cell_count_report_sample_id_fallback(
    run_component, h5mu, write_temp_h5mu, tmp_path
):
    # When the sample id column is absent from .obs, the provided --sample_id is used as the
    # sample label (the process_single_sample workflow passes the event id here).
    no_sample_col = h5mu.copy()
    no_sample_col["rna"].obs = no_sample_col["rna"].obs.drop(columns=["sample_id"])
    no_sample_col["prot"].obs = no_sample_col["prot"].obs.drop(columns=["sample_id"])

    output = tmp_path / "cell_counts.tsv"
    run_component(
        [
            "--input",
            write_temp_h5mu(no_sample_col),
            "--modality",
            "rna",
            "--prot_modality",
            "prot",
            "--sample_id_column",
            "sample_id",
            "--sample_id",
            "my_workflow_id",
            "--rna_filter_columns",
            "filter_counts_rna;filter_quantile_rna",
            "--scrublet_filter_column",
            "filter_scrublet",
            "--prot_filter_columns",
            "filter_counts_prot;filter_quantile_prot",
            "--output",
            output,
        ]
    )
    report = pd.read_csv(output, sep="\t").set_index("sample_id")
    # All observations collapse into a single sample labelled with the provided id.
    assert set(report.index) == {"my_workflow_id"}
    # rna_keep across all 8 cells = [T,T,F,F,T,F,F,F] -> 3 ; all_keep = cell_0 only -> 1
    assert report.loc["my_workflow_id", "cell_count_after_rna_filter"] == 3
    assert report.loc["my_workflow_id", "cell_count_after_all_filter"] == 1


def test_cell_count_report_skip_scrublet_and_prot(
    run_component, h5mu, write_temp_h5mu, tmp_path
):
    output = tmp_path / "cell_counts.tsv"
    run_component(
        [
            "--input",
            write_temp_h5mu(h5mu),
            "--modality",
            "rna",
            "--sample_id_column",
            "sample_id",
            "--rna_filter_columns",
            "filter_counts_rna;filter_quantile_rna",
            "--output",
            output,
        ]
    )
    assert output.is_file()
    report = pd.read_csv(output, sep="\t").set_index("sample_id")
    # No scrublet, no protein modality: those columns are omitted.
    assert "cell_count_after_scrublet_filter" not in report.columns
    assert "cell_count_after_protein_filter" not in report.columns
    assert "cell_count_after_rna_filter" in report.columns
    assert "cell_count_after_all_filter" in report.columns
    # The RNA QC aggregations are still reported; the protein ones are omitted.
    assert "rna_median_total_counts_after_all_filter" in report.columns
    assert "prot_median_total_counts_after_all_filter" not in report.columns
    # With only the RNA filter, the combined count equals the RNA count.
    assert (
        report["cell_count_after_all_filter"] == report["cell_count_after_rna_filter"]
    ).all()
    assert report.loc["s1", "cell_count_after_rna_filter"] == 2
    assert report.loc["s2", "cell_count_after_rna_filter"] == 1
    # s1 RNA keep = [T,T,F,F] -> cells 0,1 survive ; total_counts median of {1000,2000} = 1500
    assert report.loc["s1", "rna_median_total_counts_after_all_filter"] == 1500.0
    assert report.loc["s1", "rna_max_total_counts_after_all_filter"] == 2000.0


if __name__ == "__main__":
    exit(pytest.main([__file__]))
