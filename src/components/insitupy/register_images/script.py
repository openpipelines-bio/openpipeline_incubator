from insitupy import InSituData
from insitupy.tools import register_images
from insitupy.images.io import read_image
from pathlib import Path
import shutil
import sys

## VIASH START
par = {
    "input": "resources_test/insitupy/xenium_tiny_project", # Not existing yet
    "image_path": "resources_test/xenium_multichannel/xenium_multicellseg_tiny/morphology_focus/ch0000_dapi.ome.tif",
    "channel_names": ["DAPI"],
    "channel_name_for_registration": "DAPI",
    "template_image_name": "nuclei",
    "min_good_matches_per_area": 5,
    "test_flipping": True,
    "decon_scale_factor": 0.2,
    "deconvolve_template": False,
    "physicalsize": None,
    "debug": False,
    "rank_matches_for_qc": True,
    "identifier": None,
    "force_failure_qc": False,
    "raise_on_insufficient_matches": False,
    "output": "insitupy_project_registered",
}
meta = {"temp_dir": "/tmp", "resources_dir": "."}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()


def main():
    input_path = Path(par["input"])
    output_path = Path(par["output"])

    read_kwargs = {
        "image_path": par["image_path"],
        "channel_names": par["channel_names"],
        "channel_name_for_registration": par["channel_name_for_registration"],
        "template_image_name": par["template_image_name"],
        "save_registered_images": True,
        # written to <output_dir>/registered_images, moved into --output after saving
        "output_dir": meta["temp_dir"],
        "min_good_matches_per_area": par["min_good_matches_per_area"],
        "test_flipping": par["test_flipping"],
        "decon_scale_factor": par["decon_scale_factor"],
        "deconvolve_template": par["deconvolve_template"],
        # default unit is micrometers, written as an escape to keep the file ASCII-only
        "physicalsize": par["physicalsize"] or "\u00b5m",
        "debug": par["debug"],
        "rank_matches_for_qc": par["rank_matches_for_qc"],
        "identifier": par["identifier"],
        "force_failure_qc": par["force_failure_qc"],
        "raise_on_insufficient_matches": par["raise_on_insufficient_matches"],
    }

    _, _, image_axes, _ = read_image(par["image_path"])
    if image_axes in ["CYX", "YXC"]:
        assert par["channel_name_for_registration"], (
            f"'{par['image_path']}' is an IF image (axes '{image_axes}'), "
            "so --channel_name_for_registration must be provided."
        )

    logger.info("Reading InSituPy project '%s'...", input_path)
    xd = InSituData.read(input_path)

    logger.info("Registering images...")
    register_images(xd, **read_kwargs)

    logger.info("Saving output to '%s'...", output_path)
    xd.saveas(output_path, verbose=False)

    # saveas replaces the output folder, so the registered images are moved in afterwards
    registered_images_dir = Path(meta["temp_dir"]) / "registered_images"
    if registered_images_dir.exists():
        shutil.move(registered_images_dir, output_path / "registered_images")

    logger.info("Done!")


if __name__ == "__main__":
    main()