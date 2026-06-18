import mudata as mu
import scanpy as sc
import sys

### VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    "modality": "rna",
    "layer": None,
    "qc_vars": None,
    "top_n_vars": None,
    "log1p": True,
    "output": "output.h5mu",
    "output_compression": "gzip",
}
meta = {"resources_dir": "src/utils/"}
### VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from compress_h5mu import write_h5ad_to_h5mu_with_compression

logger = setup_logger()

logger.info("Reading modality %s from %s", par["modality"], par["input"])
adata = mu.read_h5ad(par["input"], mod=par["modality"]).copy()

qc_vars = par["qc_vars"] if par["qc_vars"] else []
for qc_var in qc_vars:
    if qc_var not in adata.var.columns:
        raise ValueError(
            f"Column '{qc_var}' provided to --qc_vars not found in .var "
            f"of modality '{par['modality']}'."
        )

logger.info("Calculating QC metrics with scanpy.pp.calculate_qc_metrics")
sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=qc_vars,
    percent_top=par["top_n_vars"],
    layer=par["layer"],
    log1p=par["log1p"],
    inplace=True,
)

logger.info("Writing output data to %s", par["output"])
write_h5ad_to_h5mu_with_compression(
    par["output"], par["input"], par["modality"], adata, par["output_compression"]
)

logger.info("Finished")
