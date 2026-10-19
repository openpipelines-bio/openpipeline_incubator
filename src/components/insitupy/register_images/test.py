import sys
from pathlib import Path

import numpy as np
import pytest
import tifffile
from insitupy.io import read_xenium


## VIASH START
meta = {
    "name": "register_images",
    "resources_dir": "resources_test/xenium_multichannel",
}
## VIASH END

input_path = Path(meta["resources_dir"]) / "xenium_multicellseg_tiny"

# Xenium pixel size in micrometers (see experiment.xenium)
PIXEL_SIZE = 0.2125

# Channel names of the morphology_focus images
CHANNEL_NAMES = ["DAPI", "ATP1A1_CD45_E-Cadherin", "18S", "alphaSMA_Vimentin"]


def assert_registered_output(output: Path, registered_channels: list):
    assert output.is_dir(), "Output folder was not created."

    registered_images_dir = output / "registered_images"
    assert registered_images_dir.is_dir(), (
        "Registered images folder was not created in the output folder."
    )

    # registered images are saved as <slide>__<sample>__<channel>__registered.ome.tif
    for channel in registered_channels:
        images = list(registered_images_dir.glob(f"*__{channel}__registered.ome.tif"))
        assert images, f"No registered image found for channel '{channel}'."
        assert images[0].stat().st_size > 0, (
            f"Registered image for channel '{channel}' is empty."
        )


# While the xenium- (and visium-) converter components are still being made,
# this function converts the Xenium bundle itself (no insitupy resources on s3)
@pytest.fixture(scope="module")
def insitupy_project(tmp_path_factory):
    project = tmp_path_factory.mktemp("input") / "insitupy_project"
    xd = read_xenium(input_path, verbose=False)
    xd.saveas(project, verbose=False)
    return project


# Render stacked images to register
@pytest.fixture(scope="module")
def if_image(tmp_path_factory):
    focus_files = sorted((input_path / "morphology_focus").glob("ch000*.ome.tif"))
    stack = np.stack([tifffile.imread(f, is_ome=False) for f in focus_files])

    image = tmp_path_factory.mktemp("image") / "if_image.ome.tif"
    tifffile.imwrite(
        image,
        stack,
        ome=True,
        metadata={
            "axes": "CYX",
            "Channel": {"Name": CHANNEL_NAMES},
            "PhysicalSizeX": PIXEL_SIZE,
            "PhysicalSizeY": PIXEL_SIZE,
        },
    )
    return image


def test_basic_execution(run_component, insitupy_project, if_image, tmp_path):
    output = tmp_path / "insitupy_project_registered"
    run_component(
        [
            "--input",
            str(insitupy_project),
            "--image_path",
            str(if_image),
            "--channel_names",
            ";".join(CHANNEL_NAMES),
            "--channel_name_for_registration",
            "DAPI",
            "--raise_on_insufficient_matches",
            "--output",
            str(output),
        ]
    )

    # DAPI channel is only used for alignment and is not saved
    assert_registered_output(output, [n for n in CHANNEL_NAMES if n != "DAPI"])


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
