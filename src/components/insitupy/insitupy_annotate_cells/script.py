from insitupy import InSituData, InSituExperiment
from pathlib import Path
import sys

## VIASH START
par = {
    "input": "resources_test/insitupy/insitupy_experiment",  # Not existing yet
    "cells_layers": None,
    "input_regions": None,
    "region_keys": None,
    "region_scale_factor": None,
    "input_annotations": None,
    "annotation_keys": None,
    "annotation_scale_factor": None,
    "output": "insitupy_experiment",
}
meta = {"temp_dir": "/tmp", "resources_dir": "."}
## VIASH END
sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()
input_file = Path(par["input"])
output_file = Path(par["output"])


def import_region_annot(data: InSituData):
    if par["input_regions"]:
        data.import_regions(
            files=par["input_regions"],
            keys=par["region_keys"],
            scale_factor=par["region_scale_factor"],
        )
        data.assign_regions(cells_layers=par["cells_layers"], add_to_obs=True)
    if par["input_annotations"]:
        data.import_annotations(
            files=par["input_annotations"],
            keys=par["annotation_keys"],
            scale_factor=par["annotation_scale_factor"],
        )
        data.assign_annotations(cells_layers=par["cells_layers"], add_to_obs=True)


def main():
    logger.info("Loading InSituPy object...")
    if (input_file / "metadata.parquet").exists():
        data = InSituExperiment.read(input_file)
        data.load_all()
        for _, xd in data.iterdata():
            import_region_annot(xd)
    else:
        data = InSituData.read(input_file)  # load_all=True by default
        import_region_annot(data)

    logger.info("Saving annotated project to '%s'...", output_file)
    data.saveas(output_file)
    logger.info("Done!")


if __name__ == "__main__":
    main()
