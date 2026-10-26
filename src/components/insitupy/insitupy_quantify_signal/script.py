from insitupy import InSituData, InSituExperiment
from pathlib import Path
import sys

## VIASH START
par = {
    "input": "resources_test/insitupy/insitupy_experiment",  # Not existing yet
    "cells_layer": None,
    "image_name": "CD20",
    "cells_compartment": "cells",
    "method": "median",
    "downsample_factor": None,
    "tile_size": None,
    "output": "insitupy_experiment",
}
meta = {"temp_dir": "/tmp", "resources_dir": "."}
## VIASH END
sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()
input_file = Path(par["input"])
output_file = Path(par["output"])


def insitupy_quantify_signal(data: InSituData):
    data.quantify_signal(
        image_name=par["image_name"],
        cells_layer=par["cells_layer"],
        cells_compartment=par["cells_compartment"],
        method=par["method"],
        downsample_factor=par["downsample_factor"],
        tile_size=par["tile_size"],
        add_to_obs=True,
    )


def main():
    logger.info("Loading InSituPy object...")
    if (input_file / "metadata.parquet").exists():
        data = InSituExperiment.read(input_file)
        data.load_all()
        for _, xd in data.iterdata():
            insitupy_quantify_signal(xd)
    else:
        data = InSituData.read(input_file)  # load_all=True by default
        insitupy_quantify_signal(data)

    logger.info("Saving quantified project to '%s'...", output_file)
    data.saveas(output_file)
    logger.info("Done!")


if __name__ == "__main__":
    main()
