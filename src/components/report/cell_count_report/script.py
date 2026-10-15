import sys
import mudata as mu
import pandas as pd

### VIASH START
par = {
    "input": ["input.h5mu"],
    "modality": "rna",
    "prot_modality": "prot",
    "sample_id_column": "sample_id",
    "sample_id": None,
    "rna_filter_columns": [
        "filter_counts_rna",
        "filter_quantile_rna",
        "filter_mito_rna",
    ],
    "scrublet_filter_column": "filter_scrublet",
    "prot_filter_columns": ["filter_counts_prot", "filter_quantile_prot"],
    "rna_qc_metrics": ["total_counts", "pct_counts_mt", "n_genes_by_counts"],
    "prot_qc_metrics": ["total_counts", "n_genes_by_counts"],
    "qc_metric_aggregations": ["median", "mean", "max"],
    "output": "cell_counts.tsv",
    "output_column_rna": "cell_count_after_rna_filter",
    "output_column_scrublet": "cell_count_after_scrublet_filter",
    "output_column_prot": "cell_count_after_protein_filter",
    "output_column_all": "cell_count_after_all_filter",
}
### VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()

# QC aggregations are reported for the cells surviving every filter.
QC_SUFFIX = "_after_all_filter"
# Internal column prefixes used while building the per-cell table; stripped before output.
KEEP_PREFIX = "__keep__"
METRIC_SEP = "\x00"


def combine_keep_columns(obs, columns):
    """Logical AND of boolean keep-columns. Missing values are treated as not-kept."""
    keep = pd.Series(True, index=obs.index)
    for column in columns:
        if column not in obs:
            raise ValueError(f"Column '{column}' was not found in .obs.")
        keep &= obs[column].fillna(False).astype(bool)
    return keep


def metric_columns(obs, modality, metrics):
    """Per-cell QC metric values, renamed to '<modality><SEP><metric>' for later aggregation."""
    frame = pd.DataFrame(index=obs.index)
    for metric in metrics or []:
        if metric not in obs:
            raise ValueError(
                f"QC metric column '{metric}' was not found in the '{modality}' .obs."
            )
        frame[f"{modality}{METRIC_SEP}{metric}"] = obs[metric]
    return frame


def per_cell_table(input_file, sample_label):
    """Per-cell keep-masks, QC metrics and sample grouping for one h5mu."""
    logger.info("Reading modality '%s' from %s", par["modality"], input_file)
    rna_obs = mu.read_h5ad(input_file, mod=par["modality"]).obs

    # Determine the per-sample grouping. When the sample id column is absent, fall back to the
    # provided sample label (the workflow event id), or treat the whole input as one sample.
    if par["sample_id_column"] in rna_obs:
        groups = rna_obs[par["sample_id_column"]]
    else:
        label = sample_label if sample_label is not None else "sample"
        logger.info(
            "Column '%s' not found in .obs, labelling all observations as sample '%s'.",
            par["sample_id_column"],
            label,
        )
        groups = pd.Series(label, index=rna_obs.index)

    # RNA keep-mask (AND of the RNA filter columns).
    rna_keep = combine_keep_columns(rna_obs, par["rna_filter_columns"] or [])

    # Scrublet keep-mask (optional).
    if report_scrublet:
        scrublet_keep = combine_keep_columns(rna_obs, [par["scrublet_filter_column"]])
    else:
        scrublet_keep = pd.Series(True, index=rna_obs.index)

    # Per-cell RNA QC metrics, aggregated later over the surviving cells.
    metrics = metric_columns(rna_obs, par["modality"], par["rna_qc_metrics"])

    # Protein keep-mask and QC metrics (optional), aligned to the RNA observations.
    if report_prot:
        logger.info("Reading modality '%s' from %s", par["prot_modality"], input_file)
        prot_obs = mu.read_h5ad(input_file, mod=par["prot_modality"]).obs
        prot_keep = combine_keep_columns(prot_obs, par["prot_filter_columns"] or [])
        # Observations not present in the protein modality cannot pass the protein filter.
        prot_keep = prot_keep.reindex(rna_obs.index, fill_value=False)
        prot_metrics = metric_columns(
            prot_obs, par["prot_modality"], par["prot_qc_metrics"]
        ).reindex(rna_obs.index)
        metrics = pd.concat([metrics, prot_metrics], axis=1)
    else:
        prot_keep = pd.Series(True, index=rna_obs.index)

    # The combined mask: cells passing every filter (the intersection).
    all_keep = rna_keep & scrublet_keep & prot_keep

    table = pd.DataFrame({par["sample_id_column"]: groups})
    table[f"{KEEP_PREFIX}{par['output_column_rna']}"] = rna_keep
    if report_scrublet:
        table[f"{KEEP_PREFIX}{par['output_column_scrublet']}"] = scrublet_keep
    if report_prot:
        table[f"{KEEP_PREFIX}{par['output_column_prot']}"] = prot_keep
    table[f"{KEEP_PREFIX}{par['output_column_all']}"] = all_keep
    table["__all_keep__"] = all_keep
    return pd.concat([table, metrics], axis=1)


def aggregate_qc(survivors):
    """Aggregate the per-cell QC metrics over the surviving cells, per sample.

    Returns a per-sample DataFrame whose columns are named
    `<modality>_<aggregation>_<metric>_after_all_filter`, ordered metric-major.
    """
    metric_cols = [c for c in survivors.columns if METRIC_SEP in c]
    aggregations = par["qc_metric_aggregations"]
    grouped = survivors.groupby(par["sample_id_column"], sort=True)[metric_cols].agg(
        aggregations
    )
    # Flatten the MultiIndex columns to the reference column names, ordered metric-major (all
    # aggregations of one metric together). Each (metric, aggregation) column is selected
    # individually; indexing the MultiIndex frame with a list of tuples misaligns the columns.
    result = pd.DataFrame(index=grouped.index)
    for metric_col in metric_cols:
        modality, metric = metric_col.split(METRIC_SEP)
        for aggregation in aggregations:
            name = f"{modality}_{aggregation}_{metric}{QC_SUFFIX}"
            result[name] = grouped[(metric_col, aggregation)]
    return result.round(2)


# Whether the scrublet/protein stages are reported is the same across all inputs.
report_scrublet = par["scrublet_filter_column"] is not None
report_prot = par["prot_modality"] is not None

# Fallback sample labels, one per input file (in order). Absent entries fall back to "sample".
sample_ids = par.get("sample_id") or []

# Combine the per-cell keep-masks and QC metrics of every input file. Each file is a sample (or
# a set of samples distinguished by the sample id column); the report groups across them all.
per_cell = pd.concat(
    [
        per_cell_table(input_file, sample_ids[i] if i < len(sample_ids) else None)
        for i, input_file in enumerate(par["input"])
    ],
    ignore_index=True,
)

logger.info("Aggregating cell counts per sample")
# Summing booleans per group yields the number of cells passing each stage.
keep_cols = [c for c in per_cell.columns if c.startswith(KEEP_PREFIX)]
counts = (
    per_cell.groupby(par["sample_id_column"], sort=True)[keep_cols].sum().astype(int)
)
counts.columns = [c[len(KEEP_PREFIX) :] for c in counts.columns]

logger.info("Aggregating post-filter QC metrics over the surviving cells per sample")
qc = aggregate_qc(per_cell[per_cell["__all_keep__"]])

# Every sample appears in the report (left join); samples with no surviving cells keep their
# per-stage counts and get NaN QC aggregations.
report = counts.join(qc, how="left").reset_index()

logger.info("Writing cell count report to %s", par["output"])
report.to_csv(par["output"], sep="\t", index=False)
