import json
import pandas as pd
from pathlib import Path
import anndata as ad
import h5py

## VIASH START
# inputs = list(Path("data/sample_data/sample_data").glob("*.h5mu"))
# output = "data/sample-data.json"
inputs = list(Path("data/sample_data/various_cart").glob("*.h5mu"))
output = "data/various_cart.json"
par = {
    "input": sorted([str(x) for x in inputs]),
    "output": output,
    "modality": "rna",
    "sample_id_key": "sample_id",
    "min_total_counts": 10,
    "min_num_nonzero_vars": 10,
    "obs_keys": [
        "total_counts",
        "num_nonzero_vars",
        "fraction_mitochondrial_genes",
        "fraction_ribosomal_genes",
        "pct_of_counts_in_top_50_vars",
    ],
    "cellranger_metrics_uns_key": "metrics_cellranger",
}
i = 0
mudata_file = par["input"][i]
## VIASH END


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


def main(par):
    cell_stats_dfs = []
    sample_stats_dfs = []
    metrics_cellranger_dfs = []

    print(par["input"])

    for i, mudata_file in enumerate(par["input"]):
        print(f"Processing {mudata_file}")

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

        missing_keys = [key for key in par["obs_keys"] if key not in mod_obs.columns]
        if missing_keys:
            raise ValueError(f"Missing keys in obs: {', '.join(missing_keys)}")

        sample_id = (
            mod_obs[par["sample_id_key"]].tolist()
            if par["sample_id_key"] in mod_obs.columns
            else [f"sample_{i}"] * mod_obs.shape[0]
        )

        cell_rna_stats = pd.DataFrame(
            {
                "sample_id": pd.Categorical(sample_id),
                **{key: mod_obs[key] for key in par["obs_keys"]},
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

        if par["cellranger_metrics_uns_key"] in uns:
            metrics = (
                uns[par["cellranger_metrics_uns_key"]]
                .pivot_table(
                    index=[],
                    columns="Metric Name",
                    values="Metric Value",
                    aggfunc="first",
                )
                .reset_index(drop=True)
            )

            metrics.columns.name = None
            # Remove thousands separator and convert to numeric
            metrics = metrics.map(
                lambda x: (
                    pd.to_numeric(x.replace(",", ""), errors="coerce")
                    if isinstance(x, str)
                    else x
                )
            )
            # Replace spaces with underscores in column names
            metrics.columns = metrics.columns.str.replace(" ", "_")
            for col in metrics.columns:
                metrics[col] = pd.to_numeric(metrics[col], errors="coerce")
            metrics["sample_id"] = [sample_id[0]]
            metrics_cellranger_dfs.append(metrics)

        cell_stats_dfs.append(cell_rna_stats)
        sample_stats_dfs.append(sample_summary_stats)

    combined_cell_stats = pd.concat(cell_stats_dfs, ignore_index=True)
    combined_sample_stats = pd.concat(sample_stats_dfs, ignore_index=True)
    combined_metrics_cellranger = pd.concat(metrics_cellranger_dfs, ignore_index=True)

    for df in [combined_cell_stats, combined_sample_stats, combined_metrics_cellranger]:
        df["sample_id"] = pd.Categorical(df["sample_id"])

    output = {
        "cell_rna_stats": transform_df(combined_cell_stats),
        "sample_summary_stats": transform_df(combined_sample_stats),
        "metrics_cellranger_stats": transform_df(combined_metrics_cellranger),
    }

    output_path = Path(par["output"])
    with open(output_path, "w") as f:
        json.dump(output, f, indent=2)


if __name__ == "__main__":
    main(par)
