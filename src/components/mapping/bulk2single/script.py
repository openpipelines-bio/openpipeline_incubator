import sys
import anndata as ad
import pandas as pd
import omicverse as ov

## VIASH START
par = {
    "bulk_data": "input_bulk.tsv",
    "input_sc": "input_sc.h5ad",
    "gene_id_mapping": None,
    "celltype_key": "Cell_type",
    "bulk_group": None,
    "max_single_cells": 5000,
    "top_marker_num": 500,
    "ratio_num": 1,
    "gpu": "0",
    "batch_size": 512,
    "learning_rate": 1e-4,
    "hidden_size": 256,
    "epoch_num": 3500,
    "patience": 50,
    "leiden_size": 25,
    "n_comps": 100,
    "output": "output.h5ad",
    "output_cell_fractions": None,
    "output_correlation": None,
    "output_compression": None,
}
meta = {"temp_dir": "/tmp", "resources_dir": "."}
## VIASH END


sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
logger = setup_logger()

logger.info(f"Reading bulk data from '{par['bulk_data']}'...")
bulk_data = ov.read(par["bulk_data"], index_col=0)

if par.get("gene_id_mapping"):
    logger.info(f"Mapping bulk gene identifiers using '{par['gene_id_mapping']}'...")
    bulk_data = ov.bulk.Matrix_ID_mapping(bulk_data, par["gene_id_mapping"])

logger.info(f"Reading single-cell data from '{par['input_sc']}'...")
adata_sc = ad.read_h5ad(par["input_sc"])

if par["celltype_key"] not in adata_sc.obs:
    raise ValueError(
        f"Cell type column '{par['celltype_key']}' not found in .obs of the single-cell input."
    )

gpu_par = str(par.get("gpu", "0"))
if gpu_par == "mps":
    gpu = "mps"
elif gpu_par == "cpu":
    gpu = -1
else:
    gpu = int(gpu_par)

logger.info(f"Predicting cell fraction of '{par['bulk_data']}'...")

model = ov.bulk2single.Bulk2Single(
    bulk_data=bulk_data,
    single_data=adata_sc,
    celltype_key=par["celltype_key"],
    bulk_group=par["bulk_group"],
    max_single_cells=par["max_single_cells"],
    top_marker_num=par["top_marker_num"],
    ratio_num=par["ratio_num"],
    gpu=gpu,
)
cell_fraction_prediction = model.predicted_fraction()

if par.get("output_cell_fractions"):
    logger.info(f"Writing cell fraction predictions to '{par['output_cell_fractions']}'...")
    cell_fraction_prediction.to_csv(par["output_cell_fractions"])

logger.info("Training: preprocessing...")
model.bulk_preprocess_lazy()
model.single_preprocess_lazy()
model.prepare_input()

logger.info("Training: training the VAE...")
vae_net = model.train(
    batch_size=par["batch_size"],
    learning_rate=par["learning_rate"],
    hidden_size=par["hidden_size"],
    epoch_num=par["epoch_num"],
    patience=par["patience"],
    vae_save_dir=meta["temp_dir"],
    vae_save_name="bulk2single_vae",
    generate_save_dir=meta["temp_dir"],
    generate_save_name="bulk2single_generated",
    save=False,
)

logger.info("Generating: generating single-cell data...")
generate_adata = model.generate()

logger.info("Generating: filtering out noise...")
generate_adata = model.filtered(generate_adata, leiden_size=par["leiden_size"], n_comps=par["n_comps"])


logger.info("Evaluating output...")
if par.get("output_correlation"):
    logger.info(f"Writing reference-vs-generated correlation to '{par['output_correlation']}'...")
    correlation_table = ov.bulk2single.bulk2single_plot_correlation(
        adata_sc, generate_adata, celltype_key=par["celltype_key"], return_table=True
    )
    # bulk2single_plot_correlation (omicverse 2.2.3) returns a bare, unlabeled
    # ndarray; rows/columns follow each dataset's cell types in the same
    # alphabetically-sorted order a pandas groupby(...).mean() produces, which
    # is how the row/column labels are internally derived. Verify against
    # omicverse's source if this component is ever repinned to a newer version.
    reference_celltypes = sorted(adata_sc.obs[par["celltype_key"]].unique())
    generated_celltypes = sorted(generate_adata.obs[par["celltype_key"]].unique())
    pd.DataFrame(
        correlation_table, index=reference_celltypes, columns=generated_celltypes
    ).to_csv(par["output_correlation"])

logger.info(f"Writing output to '{par['output']}'...")
generate_adata.write_h5ad(par["output"], compression=par.get("output_compression"))

logger.info("Done!")