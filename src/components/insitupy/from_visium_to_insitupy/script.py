from insitupy.io import read_visium
from pathlib import Path

## VIASH START
par = {
    "input": "resources_test/visium/Visium_FFPE_Human_Ovarian_Cancer_tiny_spaceranger",
    "output": "insitupy_project",
    "dataset_name": None,
    "sample_name": None,
    "counts_file": "filtered_feature_bc_matrix.h5",
    "fullres_image_file": None,
    "fullres_pixel_size": 4.035,
}
## VIASH END


def main():
    input_path = Path(par["input"])
    output_path = Path(par["output"])

    read_kwargs = {
        "dataset_name": par["dataset_name"],
        "sample_name": par["sample_name"],
        "counts_file": par["counts_file"],
        "fullres_image_file": par["fullres_image_file"],
        "fullres_pixel_size": par["fullres_pixel_size"],
    }

    print(f"Reading Visium bundle '{input_path}'...", flush=True)
    xd = read_visium(input_path, **read_kwargs)

    print(f"Saving InSituPy project to '{output_path}'...", flush=True)
    xd.saveas(output_path, verbose=False)

    print("Done!", flush=True)


if __name__ == "__main__":
    main()
