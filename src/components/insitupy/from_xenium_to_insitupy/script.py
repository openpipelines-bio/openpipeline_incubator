from insitupy.io import read_xenium
from pathlib import Path

## VIASH START
par = {
    "input": "resources_test/xenium_multichannel/xenium_multicellseg_tiny",
    "output": "insitupy_project",
    "backend": "insitupy",
    "nuclei_type": "focus",
    "load_cell_segmentation_images": False,
    "load_background_images": False,
    "transcript_mode": "dask",
    "restructure_transcripts": False,
    "slide_id": "slide_id",
    "sample_id": "sample_id",
}
## VIASH END


def main():
    input_path = Path(par["input"])
    output_path = Path(par["output"])

    read_kwargs = {
        "backend": par["backend"],
        "nuclei_type": par["nuclei_type"],
        "load_cell_segmentation_images": par["load_cell_segmentation_images"],
        "load_background_images": par["load_background_images"],
        "transcript_mode": par["transcript_mode"],
        "restructure_transcripts": par["restructure_transcripts"],
        "verbose": False,
    }
    # slide_id and sample_id are only used by the spatialdata backend
    if par["backend"] == "spatialdata":
        read_kwargs["slide_id"] = par["slide_id"]
        read_kwargs["sample_id"] = par["sample_id"]

    print(
        f"Reading Xenium bundle '{input_path}' with the '{par['backend']}' backend...",
        flush=True,
    )
    xd = read_xenium(input_path, **read_kwargs)

    print(f"Saving InSituPy project to '{output_path}'...", flush=True)
    xd.saveas(output_path, verbose=False)

    print("Done!", flush=True)


if __name__ == "__main__":
    main()
