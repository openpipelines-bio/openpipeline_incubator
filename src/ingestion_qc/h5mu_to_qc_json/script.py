import json
import pandas as pd
from pathlib import Path
import anndata as ad
import h5py
import sys
import os
import shutil

## VIASH START
# inputs = list(Path("data/sample_data/sample_data").glob("*.h5mu"))
# output = "data/sample-data.json"
inputs = list(Path("resources_test/qc_sample_data").glob("*.qc.cellbender.h5mu"))
output = "tmp.json"
par = {
    "input": sorted([str(x) for x in inputs]),
    "output": output,
    "output_reporting_json": "cr_struct.json",
    "modality": "rna",
    "input_obs_sample_id_key": "sample_id",
    "ingestion_method": "cellranger_multi",
    "input_obs_sample_id_key": "sample_id",
    "input_obs_total_counts_key": "total_counts",
    "input_obs_num_nonzero_vars_key": "num_nonzero_vars",
    "input_obs_fraction_mitochondrial_key": "fraction_mitochondrial",
    "input_obs_fraction_ribosomal_key": "fraction_ribosomal",
    "min_total_counts": 10,
    "min_num_nonzero_vars": 10,
    "obs_keys": [
        "total_counts",
        "num_nonzero_vars",
        "fraction_mitochondrial",
        "fraction_ribosomal",
    ],
    "cellbender_obs_keys": [
        "cellbender_background_fraction",
        "cellbender_cell_probability",
        "cellbender_cell_size",
        "cellbender_droplet_efficiency",
    ],
    "cellranger_metrics_uns_key": "metrics_cellranger",
    "metadata_obs_keys": []
}
meta = {
    "resources_dir": os.path.abspath("src/ingestion_qc/h5mu_to_qc_json/report_structure"),
}
i = 0
mudata_file = par["input"][i]

sys.path.append("src/utils")
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()

par["cellbender_obs_keys"] = {} if not par["cellbender_obs_keys"] else par["cellbender_obs_keys"]
par["metadata_obs_keys"] = {} if not par["metadata_obs_keys"] else par["metadata_obs_keys"]
# par["obs_keys"] = {} if not par["obs_keys"] else par["obs_keys"]


def transform_df(df):
    """Transform a DataFrame into the annotation object format."""
    columns = []
    for name in df.columns:
        data = df[name]

        # Determine dtype
        if pd.api.types.is_integer_dtype(data):
            dtype = "integer"
        elif pd.api.types.is_float_dtype(data):
            dtype = "numeric"
        elif pd.api.types.is_categorical_dtype(data):
            dtype = "categorical"
        else:
            raise ValueError(f"Unknown/unsupported data type for column {name}")

        column_info = {"name": name, "dtype": dtype}

        if dtype == "categorical":
            column_info["data"] = data.cat.codes.tolist()
            column_info["categories"] = data.cat.categories.tolist()
        else:
            column_info["data"] = data.tolist()

        columns.append(column_info)

    return {"num_rows": len(df), "num_cols": len(df.columns), "columns": columns}


def detect_categorical_keys(df):
    pass


def check_optional_obs_keys(obs, keys, message):
    missing_keys = [key for key in keys if key not in obs.columns]
    if missing_keys:
        logger.info(f"Missing keys in obs: {', '.join(missing_keys)}. {message}")


def transform_cellranger_metrics(uns, sample_id):
    if not par["cellranger_metrics_uns_key"] in uns:
        raise ValueError(f"Could not find cellranger metrics in uns: {par['cellranger_metrics_uns_key']}. Provide correct value for --cellranger_metrics_uns_key or make sure data was ingested using CellRanger multi.")

    cellranger_metrics = (
        uns[par["cellranger_metrics_uns_key"]]
        .pivot_table(
            index=[],
            columns="Metric Name",
            values="Metric Value",
            aggfunc="first",
        )
        .reset_index(drop=True)
    )

    cellranger_metrics.columns.name = None
    # Remove thousands separator and convert to numeric
    cellranger_metrics = cellranger_metrics.map(
        lambda x: (
            pd.to_numeric(x.replace(",", ""), errors="coerce")
            if isinstance(x, str)
            else x
        )
    )
    # Replace spaces with underscores in column names
    cellranger_metrics.columns = cellranger_metrics.columns.str.replace(" ", "_")
    for col in cellranger_metrics.columns:
        cellranger_metrics[col] = pd.to_numeric(cellranger_metrics[col], errors="coerce")
    cellranger_metrics["sample_id"] = [sample_id[0]]

    return cellranger_metrics


def generate_cellranger_stats(mod_obs, uns, sample_id, required_keys):

    # Check if celbender was run on the dataset
    if par["cellbender_obs_keys"]:
        check_optional_obs_keys(mod_obs, par["cellbender_obs_keys"], "Run cellbender first to include these metrics.")
    if par["metadata_obs_keys"]:
        check_optional_obs_keys(mod_obs, par["metadata_obs_keys"], "Make sure requested metadata colmuns are present in obs.")

    # Create cell RNA stats dataframe
    cell_rna_stats = pd.DataFrame(
        {
            "sample_id": pd.Categorical(sample_id),
            **{key: mod_obs[key] for key in required_keys},
            **{key: mod_obs[key] for key in par["cellbender_obs_keys"] if key in mod_obs.columns},
            **{key: mod_obs[key] for key in par["metadata_obs_keys"] if key in mod_obs.columns},
        }
    )

    cellranger_stats = transform_cellranger_metrics(uns, sample_id)

    return cell_rna_stats, cellranger_stats


def main(par):
    cell_stats_dfs = []
    sample_stats_dfs = []
    metrics_cellranger_dfs = []

    for i, mudata_file in enumerate(par["input"]):
        logger.info(f"Processing {mudata_file}")

        # read h5mu file
        file = h5py.File(mudata_file, "r")

        # read the necessary info
        grp_mod = file["mod"][par["modality"]]
        mod_obs = ad.experimental.read_elem(grp_mod["obs"])
        uns = ad.experimental.read_elem(file["uns"])

        # close the h5mu file
        file.close()

        barcodes_original_count = mod_obs.shape[0]

        # pre-filter cells
        if "min_total_counts" in par:
            mod_obs = mod_obs[mod_obs["total_counts"] >= par["min_total_counts"]]
        if "min_num_nonzero_vars" in par:
            mod_obs = mod_obs[mod_obs["num_nonzero_vars"] >= par["min_num_nonzero_vars"]]
        barcodes_filtered_count = mod_obs.shape[0]

        sample_id = (
            mod_obs[par["input_obs_sample_id_key"]].tolist()
            if par["input_obs_sample_id_key"] in mod_obs.columns
            else [f"sample_{i}"] * mod_obs.shape[0]
        )

        required_keys = [
            par["input_obs_total_counts_key"],
            par["input_obs_num_nonzero_vars_key"],
            par["input_obs_fraction_mitochondrial_key"],
            par["input_obs_fraction_ribosomal_key"]
            ]
        missing_keys = [key for key in required_keys if key not in mod_obs.columns]
        if missing_keys:
            raise ValueError(f"Missing keys in obs: {', '.join(missing_keys)}")

        sample_summary_stats = pd.DataFrame(
            {
                "sample_id": pd.Categorical([sample_id[0]]),
                "rna_num_barcodes": [barcodes_original_count],
                "rna_num_barcodes_filtered": [barcodes_filtered_count],
                "rna_sum_total_counts": [mod_obs[par["input_obs_total_counts_key"]].sum()],
                "rna_median_total_counts": [mod_obs[par["input_obs_total_counts_key"]].median()],
                "rna_overall_num_nonzero_vars": [mod_obs[par["input_obs_num_nonzero_vars_key"]].sum()],
                "rna_median_num_nonzero_vars": [mod_obs[par["input_obs_num_nonzero_vars_key"]].median()],
            }
        )

        if par["ingestion_method"] == "cellranger_multi":
            cell_rna_stats, cellranger_stats = generate_cellranger_stats(mod_obs, uns, sample_id, required_keys)
            metrics_cellranger_dfs.append(cellranger_stats)

        cell_stats_dfs.append(cell_rna_stats)
        sample_stats_dfs.append(sample_summary_stats)

    # Combine dataframes of all samples
    combined_cell_stats = pd.concat(cell_stats_dfs, ignore_index=True)
    combined_sample_stats = pd.concat(sample_stats_dfs, ignore_index=True)
    if par["ingestion_method"] == "cellranger_multi":
        combined_metrics_cellranger = pd.concat(metrics_cellranger_dfs, ignore_index=True)

    report_categories = [combined_cell_stats, combined_sample_stats]

    if par["ingestion_method"] == "cellranger_multi":
        report_categories.append(combined_metrics_cellranger)

    for df in report_categories:
        df["sample_id"] = pd.Categorical(df["sample_id"])

    output = {
        "cell_rna_stats": transform_df(combined_cell_stats),
        "sample_summary_stats": transform_df(combined_sample_stats)
    }

    if par["ingestion_method"] == "cellranger_multi":
        output["metrics_cellranger_stats"] = transform_df(combined_metrics_cellranger)

    logger.info(f"Writing output to {par['output']}")
    output_path = Path(par["output"])
    with open(output_path, "w") as f:
        json.dump(output, f, indent=2)

    report_structures = {
        "cellranger_multi": os.path.join(meta["resources_dir"], "report_structure/cellranger.json"),
        "xenium": os.path.join(meta["resources_dir"], "report_structure/xenium.json")
    }

    shutil.copy(report_structures[par["ingestion_method"]], par["output_reporting_json"])


if __name__ == "__main__":
    main(par)
