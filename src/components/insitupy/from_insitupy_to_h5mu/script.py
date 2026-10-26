from insitupy import InSituData, InSituExperiment
from pathlib import Path
import sys
import mudata as mu
from anndata import AnnData

## VIASH START
par = {
    "input": "resources_test/insitupy/insitupy_experiment",  # Not existing yet
    "cells_layer": None,
    "obs_sample_id": "sample_id",
    "output": "insitupy_experiment.h5mu",
    "modality": "rna",
    "output_compression": None,
}
meta = {"temp_dir": "/tmp", "resources_dir": "."}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()
input_file = Path(par["input"])
output_file = Path(par["output"])


def convert_experiment(input_path: Path) -> AnnData:
    exp = InSituExperiment.read(input_path)
    exp.load_cells()
    return exp.to_anndata(
        cells_layer=par["cells_layer"],
        label_col=par["obs_sample_id"],
    )


def convert_insitudata(input_path: Path) -> AnnData:
    xd = InSituData.read(input_path, load_all=False)
    xd.load_cells()
    if par["cells_layer"]:
        return xd.cells[par["cells_layer"]].table
    return xd.cells.table


def main():
    logger.info("Loading InSituPy object...")
    if (input_file / "metadata.parquet").exists():
        adata = convert_experiment(input_file)
    else:
        adata = convert_insitudata(input_file)

    logger.info(
        "Writing %d cells x %d genes to modality '%s' of '%s'...",
        adata.n_obs,
        adata.n_vars,
        par["modality"],
        output_file,
    )
    mdata = mu.MuData({par["modality"]: adata})
    mdata.write_h5mu(output_file, compression=par["output_compression"])

    logger.info("Done!")


if __name__ == "__main__":
    main()
