import sys
from pathlib import Path
import pandas as pd

import numpy as np
import pytest
import tifffile
from insitupy.io import read_xenium, read_visium
from PIL import Image


## VIASH START
meta = {
    "name": "register_images",
    "resources_dir": "resources_test",
}
## VIASH END


def input_spec(technology: str):
    input_specs = {}
    if technology == "xenium":
        input_path = (
            Path(meta["resources_dir"])
            / "xenium_multichannel"
            / "xenium_multicellseg_tiny"
        )
        # Xenium pixel size in micrometers (see experiment.xenium)
        pixel_size = 0.2125
        # Channel names of the morphology_focus images
        channels = ["DAPI", "ATP1A1_CD45_E-Cadherin", "18S", "alphaSMA_Vimentin"]
    else:
        input_path = (
            Path(meta["resources_dir"])
            / "visium"
            / "Visium_FFPE_Human_Ovarian_Cancer_tiny_spaceranger"
        )
        # Based on Space Ranger glossary and output documentation, PIXEL_SIZE is inferred
        # Because "--fullres_pixel_size" is a required argument for read_visium
        pos = pd.read_csv(input_path / "spatial" / "tissue_positions.csv")
        # neighbour two columns to the right
        right = pos.assign(array_col=pos.array_col - 2)
        pairs = pos.merge(right, on=["array_row", "array_col"], suffixes=("", "_nb"))
        dist = np.hypot(
            pairs.pxl_row_in_fullres_nb - pairs.pxl_row_in_fullres,
            pairs.pxl_col_in_fullres_nb - pairs.pxl_col_in_fullres,
        )
        pixel_size = 100 / dist.median()
        channels = ["HE"]
    input_specs["input_path"] = input_path
    input_specs["pixel_size"] = pixel_size
    input_specs["channels"] = channels
    return input_specs


# While the xenium- (and visium-) converter components are still being made,
# this fixture converts the Xenium/Visium bundle itself (no insitupy resources on s3).
# The technology is passed by the test with indirect parametrization.
@pytest.fixture(scope="module")
def insitupy_project(request, tmp_path_factory):
    technology = request.param
    spec = input_spec(technology)
    project = tmp_path_factory.mktemp(technology) / "insitupy_project"
    if technology == "xenium":
        xd = read_xenium(spec["input_path"], verbose=False)
    else:
        xd = read_visium(
            spec["input_path"], fullres_pixel_size=spec["pixel_size"], verbose=False
        )
    xd.saveas(project, verbose=False)
    return project


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


def assert_debug_output(output: Path, qc_suffixes: list):
    # with --debug, QC files are saved as <identifier>__<suffix> in registered_images/registration_qc
    qc_dir = output / "registered_images" / "registration_qc"
    assert qc_dir.is_dir(), "Registration QC folder was not created with --debug."

    for suffix in qc_suffixes:
        files = list(qc_dir.glob(f"*__{suffix}"))
        assert files, f"No registration QC file found ending in '__{suffix}'."
        assert files[0].stat().st_size > 0, (
            f"Registration QC file '{files[0].name}' is empty."
        )


# Render stacked images to register for Xenium (CYX) image to be registered
@pytest.fixture(scope="module")
def if_image(tmp_path_factory):
    spec = input_spec("xenium")
    focus_files = sorted(
        (spec["input_path"] / "morphology_focus").glob("ch000*.ome.tif")
    )
    stack = np.stack([tifffile.imread(f, is_ome=False) for f in focus_files])

    image = tmp_path_factory.mktemp("image") / "if_image.ome.tif"
    tifffile.imwrite(
        image,
        stack,
        ome=True,
        metadata={
            "axes": "CYX",
            "Channel": {"Name": spec["channels"]},
            "PhysicalSizeX": spec["pixel_size"],
            "PhysicalSizeY": spec["pixel_size"],
        },
    )
    return image


# Render flipped H&E image to register for Visium (YXS), so the registration has a
# transform to find. tissue_hires_scalef is 1.0 for this dataset, so the hires image
# is at full resolution.
@pytest.fixture(scope="module")
def he_image(tmp_path_factory):
    spec = input_spec("visium")
    rgb = np.asarray(
        Image.open(spec["input_path"] / "spatial" / "tissue_hires_image.png")
    )[..., :3]

    image = tmp_path_factory.mktemp("image") / "he_image.ome.tif"
    tifffile.imwrite(
        image,
        np.ascontiguousarray(np.fliplr(rgb)),
        photometric="rgb",
        ome=True,
        metadata={
            "axes": "YXS",
            "PhysicalSizeX": spec["pixel_size"],
            "PhysicalSizeY": spec["pixel_size"],
        },
    )
    return image


@pytest.mark.parametrize("insitupy_project", ["xenium"], indirect=True)
def test_xenium_if_image(run_component, insitupy_project, if_image, tmp_path):
    channels = input_spec("xenium")["channels"]
    output = tmp_path / "insitupy_project_registered"
    run_component(
        [
            "--input",
            str(insitupy_project),
            "--image_path",
            str(if_image),
            "--channel_names",
            ";".join(channels),
            "--channel_name_for_registration",
            "DAPI",
            "--raise_on_insufficient_matches",
            "--debug",
            "--output",
            str(output),
        ]
    )

    # DAPI channel is only used for alignment and is not saved
    assert_registered_output(output, [n for n in channels if n != "DAPI"])
    assert_debug_output(output, ["transform.csv", "matches_overview.png"])


@pytest.mark.parametrize("insitupy_project", ["visium"], indirect=True)
def test_visium_he_image(run_component, insitupy_project, he_image, tmp_path):
    channels = input_spec("visium")["channels"]
    output = tmp_path / "insitupy_project_registered"
    run_component(
        [
            "--input",
            str(insitupy_project),
            "--image_path",
            str(he_image),
            "--channel_names",
            ";".join(channels),
            "--template_image_name",
            "hires",
            "--deconvolve_template",
            "--raise_on_insufficient_matches",
            "--debug",
            "--output",
            str(output),
        ]
    )

    assert_registered_output(output, channels)
    assert_debug_output(
        output,
        [
            "transform.csv",
            "matches_overview.png",
            "deconvolved_target.png",
            "deconvolved_template.png",
        ],
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
