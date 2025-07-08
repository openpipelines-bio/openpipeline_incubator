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
inputs = list(Path("resources_test_after_running_script/qc_sample_data").glob("*.qc.h5mu"))
output = "tmp.json"
par = {
    "input": sorted([str(x) for x in inputs]),
    # "input": ["resources_test/spatial_qc_sample_data/xenium_tiny.qc.h5mu", "resources_test/spatial_qc_sample_data/xenium_tiny.qc.h5mu"],
    "output": "sc_data.json",
    "output_reporting_json": "sc_report_structure.json",
    "modality": "rna",
    "ingestion_method": "cellranger_multi",
    "obs_sample_id": "sample_id",
    "obs_total_counts": "total_counts",
    "obs_num_nonzero_vars": "num_nonzero_vars",
    "obs_fraction_mitochondrial": "fraction_mitochondrial",
    "obs_fraction_ribosomal": "fraction_ribosomal",
    "min_total_counts": 20,
    "min_num_nonzero_vars": 20,
    "obs_cellbender": [
        "cellbender_background_fraction",
        "cellbender_cell_probability",
        "cellbender_cell_size",
        "cellbender_droplet_efficiency",
    ],
    "uns_cellranger_metrics": "metrics_cellranger",
    "obs_metadata": ["cell_type"],
    "obs_nucleus_area": "nucleus_area",
    "obs_cell_area": "cell_area",
    "obs_x_coord": "x_coord",
    "obs_y_coord": "y_coord",
    "obs_control_probe_counts": "control_probe_counts",
    "obs_control_codeword_counts": "control_codeword_counts"
}
meta = {
    "resources_dir": os.path.abspath("src/ingestion_qc/h5mu_to_qc_json"),
}
i = 0
mudata_file = par["input"][i]

sys.path.append("src/utils")
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()

par["obs_cellbender"] = {} if not par["obs_cellbender"] else par["obs_cellbender"]


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


def check_optional_obs_keys(obs, keys, message):
    missing_keys = [key for key in keys if key not in obs.columns]
    if missing_keys:
        logger.info(f"Missing keys in obs: {', '.join(missing_keys)}. {message}")


def transform_cellranger_metrics(uns, sample_id):
    if not par["uns_cellranger_metrics"] in uns:
        raise ValueError(f"Could not find cellranger metrics in uns: {par['uns_cellranger_metrics']}. Provide correct value for --uns_cellranger_metrics or make sure data was ingested using CellRanger multi.")

    cellranger_metrics = (
        uns[par["uns_cellranger_metrics"]]
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


def format_cellbender_columns(mod_obs):
    # Check if celbender was run on the dataset
    if par["obs_cellbender"]:
        check_optional_obs_keys(mod_obs, par["obs_cellbender"], "Run cellbender first to include these metrics.")

    cellbender_obs_keys = [column for column in par["obs_cellbender"] if column in mod_obs]

    for key in cellbender_obs_keys:
        if not pd.api.types.is_float_dtype(mod_obs[key]):
            try:
                mod_obs[key] = mod_obs[key].astype("float16")
            except ValueError:
                raise ValueError(f"Could not convert column {key} to a float dtype. Please make sure all cellbender metrics are numeric.")

    return cellbender_obs_keys, mod_obs


def format_required_columns(required_keys, mod_obs):

    for key in required_keys:
        if not pd.api.types.is_numeric_dtype(mod_obs[key]):
            raise ValueError(f"Column {key} must be a numeric dtype.")

    if not pd.api.types.is_integer_dtype(mod_obs[par["obs_total_counts"]]):
        logger.info(f"Converting {par['obs_total_counts']} from {mod_obs[par['obs_total_counts']].dtype} to integer dtype...")
        mod_obs[par["obs_total_counts"]] = mod_obs[par["obs_total_counts"]].astype(int)

    if not pd.api.types.is_integer_dtype(mod_obs[par["obs_num_nonzero_vars"]]):
        logger.info(f"Converting {par['obs_num_nonzero_vars']} from {mod_obs[par['obs_num_nonzero_vars']].dtype} to integer dtype...")
        mod_obs[par["obs_num_nonzero_vars"]] = mod_obs[par["obs_num_nonzero_vars"]].astype(int)

    if not pd.api.types.is_float_dtype(mod_obs[par["obs_fraction_mitochondrial"]]):
        logger.info(f"Converting {par['obs_fraction_mitochondrial']} from {mod_obs[par['obs_fraction_mitochondrial']].dtype} to float dtype...")
        mod_obs[par["obs_fraction_mitochondrial"]] = mod_obs[par["obs_fraction_mitochondrial"]].astype("float16")

    if not pd.api.types.is_float_dtype(mod_obs[par["obs_fraction_ribosomal"]]):
        logger.info(f"Converting {par['obs_fraction_ribosomal']} from {mod_obs[par['obs_fraction_ribosomal']].dtype} to float dtype...")
        mod_obs[par["obs_fraction_ribosomal"]] = mod_obs[par["obs_fraction_ribosomal"]].astype("float16")

    return mod_obs


def format_categorical_columns(mod_obs):
    # Fetch all categorical columns for grouping if no columns are provided
    if not par["obs_metadata"]:
        metadata_obs_keys = mod_obs.select_dtypes(include=["object", "category"]).columns.tolist()
        if par["obs_sample_id"] in metadata_obs_keys:
            metadata_obs_keys.remove(par["obs_sample_id"])
    else:
        check_optional_obs_keys(mod_obs, par["obs_metadata"], "Make sure requested metadata colmuns are present in obs.")
        metadata_obs_keys = [key for key in par["obs_metadata"] if key in mod_obs]

    for key in metadata_obs_keys:
        if not isinstance(key, pd.CategoricalDtype):
            logger.info(f"{key} is not a categorical dtype. Converting {key} from {mod_obs[key].dtype} to categorical dtype...")
            mod_obs[key] = mod_obs[key].astype(str).astype("category")

    return metadata_obs_keys, mod_obs


def generate_cellranger_stats(mod_obs, uns, sample_id, required_keys):

    # Format required columns
    mod_obs = format_required_columns(required_keys, mod_obs)

    # Fetch and format  all categorical columns for grouping
    metadata_obs_keys, mod_obs = format_categorical_columns(mod_obs)

    # Fetch and format cellbender columns
    cellbender_obs_keys, mod_obs = format_cellbender_columns(mod_obs)

    # Create cell RNA stats dataframe
    cell_rna_stats = pd.DataFrame(
        {
            "sample_id": pd.Categorical(sample_id),
            **{key: mod_obs[key] for key in required_keys},
            **{key: mod_obs[key] for key in cellbender_obs_keys},
            **{key: mod_obs[key] for key in metadata_obs_keys},
        }
    )

    cellranger_stats = transform_cellranger_metrics(uns, sample_id)

    return cell_rna_stats, cellranger_stats


def format_xenium_columns(mod_obs):

    mod_obs["nucleus_ratio"] = mod_obs[par["obs_nucleus_area"]] / mod_obs[par["obs_cell_area"]]

    xenium_formatted_columns = [par["obs_cell_area"], "nucleus_ratio", "x_coord", "y_coord"]
    for key in xenium_formatted_columns:
        mod_obs[key] = mod_obs[key].astype("float16")

    return mod_obs, xenium_formatted_columns


def generate_xenium_stats(mod_obs, sample_id, required_keys):

    # Format required columns
    mod_obs = format_required_columns(required_keys, mod_obs)

    # Format xenium-specific columns
    mod_obs, xenium_formatted_columns = format_xenium_columns(mod_obs)

    # Fetch and format  all categorical columns for grouping
    metadata_obs_keys, mod_obs = format_categorical_columns(mod_obs)

    # Create cell RNA stats dataframe
    cell_rna_stats = pd.DataFrame(
        {
            "sample_id": pd.Categorical(sample_id),
            **{key: mod_obs[key] for key in required_keys},
            **{key: mod_obs[key] for key in xenium_formatted_columns},
            **{key: mod_obs[key] for key in metadata_obs_keys}
        }
    )

    return cell_rna_stats


def concatenate_dataframes(dfs):
    '''Concatenates a list of dataframes into a single dataframe, preserving categorical columns.'''
    df = pd.concat(dfs, ignore_index=True)

    # Find categorical columns that became object columms
    for col in df.columns:
        if any(df[col].dtype.name == 'category' for df in dfs if col in df.columns):
            # Get all categorical series for this column
            cat_series = [df[col] for df in dfs if col in df.columns and df[col].dtype.name == 'category']
            if cat_series:
                # Union the categories and apply to result
                unioned = pd.api.types.union_categoricals(cat_series)
                df[col] = pd.Categorical(df[col], categories=unioned.categories)
    return df


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
        mod_obsm = ad.experimental.read_elem(grp_mod["obsm"])
        uns = ad.experimental.read_elem(file["uns"])

        # close the h5mu file
        file.close()

        barcodes_original_count = mod_obs.shape[0]

        # Add coordinates to obs before filtering
        if par["ingestion_method"] == "xenium":
            mod_obs["x_coord"] = mod_obsm["spatial"][:, 0]
            mod_obs["y_coord"] = mod_obsm["spatial"][:, 1]

        # Pre-filter cells
        logger.info("Pre-filtering cells based on counts...")
        if "min_total_counts" in par:
            mod_obs = mod_obs[mod_obs["total_counts"] >= par["min_total_counts"]]
        if "min_num_nonzero_vars" in par:
            mod_obs = mod_obs[mod_obs["num_nonzero_vars"] >= par["min_num_nonzero_vars"]]
        barcodes_filtered_count = mod_obs.shape[0]

        # Detect sample id's
        logger.info("Detecting sample id's...")
        sample_id = (
            mod_obs[par["obs_sample_id"]].tolist()
            if par["obs_sample_id"] in mod_obs.columns
            else [f"sample_{i}"] * mod_obs.shape[0]
        )

        # Generating sample summary statistics
        logger.info("Generating sample summary statistics...")
        required_keys = [
            par["obs_total_counts"],
            par["obs_num_nonzero_vars"],
            par["obs_fraction_mitochondrial"],
            par["obs_fraction_ribosomal"]
            ]
        missing_keys = [key for key in required_keys if key not in mod_obs.columns]
        if missing_keys:
            raise ValueError(f"Missing keys in obs: {', '.join(missing_keys)}")

        sample_summary = {
            "sample_id": pd.Categorical([sample_id[0]]),
            "rna_num_barcodes": [barcodes_original_count],
            "rna_num_barcodes_filtered": [barcodes_filtered_count],
            "rna_sum_total_counts": [mod_obs[par["obs_total_counts"]].sum()],
            "rna_median_total_counts": [mod_obs[par["obs_total_counts"]].median()],
            "rna_overall_num_nonzero_vars": [mod_obs[par["obs_num_nonzero_vars"]].sum()],
            "rna_median_num_nonzero_vars": [mod_obs[par["obs_num_nonzero_vars"]].median()],
        }

        if par["ingestion_method"] == "xenium":
            sample_summary["control_probe_percentage"] = mod_obs[par["obs_control_probe_counts"]].sum() / mod_obs["total_counts"].sum() * 100
            sample_summary["negative_decoding_percentage"] = mod_obs[par["obs_control_codeword_counts"]].sum() / mod_obs["total_counts"].sum() * 100

        sample_summary_stats = pd.DataFrame(sample_summary)

        if par["ingestion_method"] == "cellranger_multi":
            cell_rna_stats, cellranger_stats = generate_cellranger_stats(mod_obs, uns, sample_id, required_keys)
            metrics_cellranger_dfs.append(cellranger_stats)

        if par["ingestion_method"] == "xenium":
            cell_rna_stats = generate_xenium_stats(mod_obs, sample_id, required_keys)

        cell_stats_dfs.append(cell_rna_stats)
        sample_stats_dfs.append(sample_summary_stats)

    # Combine dataframes of all samples
    logger.info("Combining data of all samples into single object...")
    combined_cell_stats = concatenate_dataframes(cell_stats_dfs)
    combined_sample_stats = concatenate_dataframes(sample_stats_dfs)
    if par["ingestion_method"] == "cellranger_multi":
        combined_metrics_cellranger = concatenate_dataframes(metrics_cellranger_dfs)

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

    logger.info(f"Writing output data json to {par['output']}")
    output_path = Path(par["output"])
    with open(output_path, "w") as f:
        json.dump(output, f, indent=2)

    report_structures = {
        "cellranger_multi": os.path.join(meta["resources_dir"], "report_structure/cellranger.json"),
        "xenium": os.path.join(meta["resources_dir"], "report_structure/xenium.json")
    }

    logger.info(f"Writing output report structure json to {par['output_reporting_json']}")
    shutil.copy(report_structures[par["ingestion_method"]], par["output_reporting_json"])


if __name__ == "__main__":
    main(par)
