
import subprocess
import sys
from pathlib import Path

import pytest
from insitupy import InSituData


## VIASH START
meta = {
    "name": "from_xenium_to_insitupy",
    "resources_dir": "resources_test/xenium_multichannel/xenium_multicellseg_tiny",
}
## VIASH END

input_path = Path(meta["resources_dir"])


def test_basic_execution(run_component, random_path):
    output = random_path()
    run_component(
        "--input",
        input_path, 
        "--backend",
        backend,
        "nuclei_type", 
        ,
        "load_cell_segmentation_images": par["load_cell_segmentation_images"],
        "load_background_images": par["load_background_images"],
        "transcript_mode": par["transcript_mode"],
        "restructure_transcripts": par["restructure_transcripts"],
    )




if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))