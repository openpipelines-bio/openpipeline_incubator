import json
import pandas as pd
from pathlib import Path
import anndata as ad
import h5py
import sys

## VIASH START
inputs = ["resources_test/spatial_qc_sample_data/xenium/xenium_tiny_qc.h5mu", "resources_test/spatial_qc_sample_data/xenium/xenium_tiny_qc.h5mu"]
output = "tmp.json"
par = {
    "input": sorted([str(x) for x in inputs]),
    "output": output,
    "modality": "rna",
    "sample_id_key": "sample_id",
    "min_total_counts": 1,
    "min_num_nonzero_vars": 1,
    "obs_keys": [
        "total_counts",
        "num_nonzero_vars",
        "fraction_mitochondrial_genes",
        "fraction_ribosomal_genes",
        "control_probe_counts",
        "control_codeword_counts",
        "cell_area",
        "nucleus_area"
    ],
    "spatial_obsm_key": "spatial",
    "metadata_obs_keys": None
}
meta = {
        "resources_dir": "src/utils"
    }

# i = 0
# mudata_file = par["input"][i]
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()

par["spatial_obsm_key"] = {} if not par["spatial_obsm_key"] else par["spatial_obsm_key"]
par["metadata_obs_keys"] = {} if not par["metadata_obs_keys"] else par["metadata_obs_keys"]
par["obs_keys"] = {} if not par["obs_keys"] else par["obs_keys"]


def transform_df(df):
    """Transform a DataFrame into the annotation object format."""
    columns = []
    for name in df.columns:
        print(f"Processing column {name}")
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


def main(par):
    cell_stats_dfs = []
    sample_stats_dfs = []

    for i, mudata_file in enumerate(par["input"]):
        logger.info(f"Processing {mudata_file}")

        # read h5mu file
        file = h5py.File(mudata_file, "r")

        # read the necessary info
        grp_mod = file["mod"][par["modality"]]
        mod_obs = ad.experimental.read_elem(grp_mod["obs"])
        mod_obsm = ad.experimental.read_elem(grp_mod["obsm"])

        mod_obs["spatial_coord_x"] = mod_obsm[par["spatial_obsm_key"]][:, 0]
        mod_obs["spatial_coord_y"] = mod_obsm[par["spatial_obsm_key"]][:, 1]

        # close the h5mu file
        file.close()

        barcodes_original_count = mod_obs.shape[0]

        # pre-filter cells
        if "min_total_counts" in par:
            mod_obs = mod_obs[mod_obs["total_counts"] >= par["min_total_counts"]]
        if "min_num_nonzero_vars" in par:
            mod_obs = mod_obs[mod_obs["num_nonzero_vars"] >= par["min_num_nonzero_vars"]]
        barcodes_filtered_count = mod_obs.shape[0]

        missing_keys = [key for key in par["obs_keys"] if key not in mod_obs.columns]
        if missing_keys:
            raise ValueError(f"Missing keys in obs: {', '.join(missing_keys)}")

        if par["metadata_obs_keys"]:
            check_optional_obs_keys(mod_obs, par["metadata_obs_keys"], "Make sure requested metadata colmuns are present in obs.")

        sample_id = (
            mod_obs[par["sample_id_key"]].tolist()
            if par["sample_id_key"] in mod_obs.columns
            else [f"sample_{i}"] * mod_obs.shape[0]
        )

        cell_rna_stats = pd.DataFrame(
            {
                "sample_id": pd.Categorical(sample_id),
                **{key: mod_obs[key] for key in par["obs_keys"]},
                **{key: mod_obs[key] for key in ["spatial_coord_x", "spatial_coord_y"]},
                **{key: mod_obs[key] for key in par["metadata_obs_keys"] if key in mod_obs.columns},
            }
        )

        sample_summary_stats = pd.DataFrame(
            {
                "sample_id": pd.Categorical([sample_id[0]]),
                "rna_num_barcodes": [barcodes_original_count],
                "rna_num_barcodes_filtered": [barcodes_filtered_count],
                "rna_sum_total_counts": [mod_obs["total_counts"].sum()],
                "rna_median_total_counts": [mod_obs["total_counts"].median()],
                "rna_overall_num_nonzero_vars": [mod_obs["num_nonzero_vars"].sum()],
                "rna_median_num_nonzero_vars": [mod_obs["num_nonzero_vars"].median()],
            }
        )

        cell_stats_dfs.append(cell_rna_stats)
        sample_stats_dfs.append(sample_summary_stats)

    combined_cell_stats = pd.concat(cell_stats_dfs, ignore_index=True)
    combined_sample_stats = pd.concat(sample_stats_dfs, ignore_index=True)

    for df in [combined_cell_stats, combined_sample_stats]:
        df["sample_id"] = pd.Categorical(df["sample_id"])

    output = {
        "cell_rna_stats": transform_df(combined_cell_stats),
        "sample_summary_stats": transform_df(combined_sample_stats)
    }

    output_path = Path(par["output"])
    with open(output_path, "w") as f:
        json.dump(output, f, indent=2)


if __name__ == "__main__":
    main(par)
