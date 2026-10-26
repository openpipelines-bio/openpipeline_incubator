from insitupy import InSituData, InSituExperiment
from insitupy.tools import calc_distance_of_cells_from
from pathlib import Path
import sys

## VIASH START
par = {
    "input": "resources_test/insitupy/insitupy_experiment",  # Not existing yet
    "cells_layer": None,
    "annotation_tuple": ["pathologist", "tumor"],
    "region_tuple": None,
    "key_to_save": None,
    "output": "insitupy_experiment",
}
meta = {"temp_dir": "/tmp", "resources_dir": "."}
## VIASH END
sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()
input_file = Path(par["input"])
output_file = Path(par["output"])


def to_tuple(values, arg):
    if values is None:
        return None
    if len(values) != 2:
        raise ValueError(f"--{arg} needs exactly two values (key, name), got {values}.")
    return tuple(values)


annotation_tuple = to_tuple(par["annotation_tuple"], "annotation_tuple")
region_tuple = to_tuple(par["region_tuple"], "region_tuple")


def calc_distance(data: InSituData):
    # The region only restricts which cells are included; the distances are always
    # calculated from the annotation.
    calc_distance_of_cells_from(
        data=data,
        annotation_tuple=annotation_tuple,
        region_tuple=region_tuple,
        cells_layer=par["cells_layer"],
        key_to_save=par["key_to_save"],
    )


def main():
    logger.info("Loading InSituPy object...")
    if (input_file / "metadata.parquet").exists():
        data = InSituExperiment.read(input_file)
        data.load_all()
        for _, xd in data.iterdata():
            calc_distance(xd)
    else:
        data = InSituData.read(input_file)  # load_all=True by default
        calc_distance(data)

    logger.info("Saving project with distances to '%s'...", output_file)
    data.saveas(output_file)
    logger.info("Done!")


if __name__ == "__main__":
    main()
