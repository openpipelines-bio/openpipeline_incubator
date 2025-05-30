import shutil
import anndata as ad
import h5py
import sys

## VIASH START
par = {
    "input": "resources_test/qc_sample_data/sample_one.qc.h5mu",
    # "input": "resources_test/spatial_qc_sample_data/xenium/xenium_tiny_qc.h5mu",
    # "input": "/Users/dorienroosen/code/openpipeline_spatial/resources_test/cosmx/Lung5_Rep2_tiny.h5mu",
    "output": "output.h5mu",
    "output_uns_ingestion_method": "ingestion_method",
    "modality": "rna"
}
meta = {
    "resources_dir": "src/utils"
}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()


def main(par):

    # read h5mu file
    with h5py.File(par["input"], "r") as file:
        mod = file["mod"][par["modality"]]
        uns = ad.experimental.read_elem(file["uns"])
        mod_obs = ad.experimental.read_elem(mod["obs"])
        mod_uns = ad.experimental.read_elem(mod["uns"])

    # detect ingestion method
    ingestion_methods = {
            "cellranger_multi": "metrics_cellranger" in uns,
            "xenium": all(key in mod_obs for key in ["segmentation_method", "nucleus_area"]),
            "cosmx": "spatial" in mod_uns and "fov" in mod_obs
        }

    # make sure only one ingestion method is detected
    detected_methods = [method for method, detected in ingestion_methods.items() if detected]
    methods_count = len(detected_methods)

    if methods_count == 1:
        detected_method = detected_methods[0]
        logger.info(f"Detected ingestion method {detected_method}")

    elif methods_count == 0:
        raise ValueError("No ingestion method detected")

    else:
        raise ValueError(f"Multiple ingestion methods detected: {', '.join(detected_methods)}")

    # write ingestion method to .uns if not already present and save updated h5mu
    if par["output_uns_ingestion_method"] in mod_uns:
        if mod_uns[par["output_uns_ingestion_method"]] != detected_method:
            raise ValueError(f"Output uns key {par['output_uns_ingestion_method']} already exists and contains different value `{mod_uns[par['output_uns_ingestion_method']]}` then detected method `{detected_method}`.")
        shutil.copy(par["input"], par["output"])

    else:
        shutil.copy(par["input"], par["output"])
        with h5py.File (par["output"], "r+") as out_file:
            out_file["uns"][par["output_uns_ingestion_method"]] = detected_method


if __name__ == "__main__":
    main(par)
