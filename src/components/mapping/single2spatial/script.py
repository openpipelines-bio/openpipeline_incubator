import sys

import anndata as ad
import numpy as np
import omicverse as ov

## VIASH START
par = {
    "input_sc": "input_sc.h5ad",
    "input_spatial": "input_spatial.h5ad",
    "celltype_key": "Cell_type",
    "obsm_spatial_coordinates": "spatial",
    "spot_key": ["xcoord", "ycoord"],
    "spot_num": 500,
    "cell_num": 10,
    "top_marker_num": 500,
    "marker_used": True,
    "k": 10,
    "num_epochs": 1000,
    "batch_size": 1000,
    "predicted_size": 32,
    "layer": None,
    "device_type": "cpu",
    "output": "output.h5ad",
    "obsm_output": "spatial",
    "output_compression": None,
}
meta = {"temp_dir": "/tmp", "resources_dir": "."}
## VIASH END


sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()

logger.info(f"Reading single-cell data from '{par['input_sc']}'...")
adata_sc = ad.read_h5ad(par["input_sc"])

logger.info(f"Reading spatial data from '{par['input_spatial']}'...")
adata_sp = ad.read_h5ad(par["input_spatial"])

if par["celltype_key"] not in adata_sc.obs:
    raise ValueError(
        f"Cell type column '{par['celltype_key']}' not found in .obs of the single-cell input."
    )

if par.get("layer"):
    for label, adata in (("single-cell", adata_sc), ("spatial", adata_sp)):
        if par["layer"] not in adata.layers:
            raise ValueError(f"Layer '{par['layer']}' not found in the {label} input.")
        adata.X = adata.layers[par["layer"]]

spot_key = list(par["spot_key"])
if len(spot_key) != 2:
    raise ValueError("--spot_key must contain exactly two column names (x and y).")

obsm_key = par.get("obsm_spatial_coordinates")
if obsm_key and obsm_key in adata_sp.obsm:
    logger.info(f"Copying spatial coordinates from .obsm['{obsm_key}']...")
    coords = np.asarray(adata_sp.obsm[obsm_key])[:, :2]
    adata_sp.obs[spot_key[0]] = coords[:, 0]
    adata_sp.obs[spot_key[1]] = coords[:, 1]
elif not all(key in adata_sp.obs for key in spot_key):
    raise ValueError(
        f"Spatial coordinates not found: '{obsm_key}' not in .obsm and {spot_key} not in .obs of the spatial input."
    )

device_type_par = par.get("device_type", "cpu")
if device_type_par == "gpu":
    gpu = 0
else:
    gpu = -1

logger.info("Fitting Single2Spatial mapper...")
st_model = ov.bulk2single.Single2Spatial(
    single_data=adata_sc,
    spatial_data=adata_sp,
    celltype_key=par["celltype_key"],
    spot_key=spot_key,
    top_marker_num=par["top_marker_num"],
    marker_used=par["marker_used"],
    gpu=gpu,
)

if device_type_par != "cpu" and st_model.used_device.type == "cpu":
    logger.warning(
        f"--device_type {device_type_par!r} was requested but Single2Spatial fell "
        f"back to CPU (used_device={st_model.used_device}); this container's torch "
        "build may lack CUDA support, or no compatible device was found."
    )
else:
    logger.info(f"Single2Spatial is using device: {st_model.used_device}")

logger.info("Training mapping model and predicting spatial coordinates...")
sp_adata = st_model.train(
    spot_num=par["spot_num"],
    cell_num=par["cell_num"],
    df_save_dir=meta["temp_dir"],
    df_save_name="single2spatial_model",
    k=par["k"],
    num_epochs=par["num_epochs"],
    batch_size=par["batch_size"],
    predicted_size=par["predicted_size"],
    save=False,
)

obsm_output = par["obsm_output"]
if obsm_output != "X_spatial":
    sp_adata.obsm[obsm_output] = sp_adata.obsm.pop("X_spatial")

# Single2Spatial returns a pandas nullable "string"-dtype var index, which
# anndata refuses to write by default (opt-in only, for backwards compat
# with anndata < 0.11 readers).
sp_adata.var.index = sp_adata.var.index.astype(str)

logger.info(f"Writing output to '{par['output']}'...")
sp_adata.write_h5ad(par["output"], compression=par.get("output_compression"))

logger.info("Done!")
