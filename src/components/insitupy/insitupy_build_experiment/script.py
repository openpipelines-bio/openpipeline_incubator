from insitupy import InSituData, InSituExperiment
from pathlib import Path
import sys
import pandas as pd

## VIASH START
par = {
    "input": [
        "resources_test/insitupy/sample_1_insitupy",  # Not existing yet
        "resources_test/insitupy/sample_2_insitupy",  # Not existing yet
    ],
    "sample_id": ["sample_1", "sample_2"],
    "metadata": None,
    "crop": False,
    "xlim": [0, 2000],
    "ylim": [0, 2000],
    "output": "insitupy_experiment",
}
meta = {"temp_dir": "/tmp", "resources_dir": "."}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()


def check_arguments():
    assert len(par["input"]) == len(par["sample_id"]), (
        f"--input ({len(par['input'])} values) and --sample_id "
        f"({len(par['sample_id'])} values) must have the same length."
    )
    assert len(set(par["sample_id"])) == len(par["sample_id"]), (
        "--sample_id values must be unique."
    )
    if par["crop"]:
        for lim in ("xlim", "ylim"):
            assert par[lim] and len(par[lim]) == 2, (
                f"--crop requires --{lim} with exactly two values (min max)."
            )


def read_sample_metadata():
    # Extra per-sample metadata (e.g. patient, condition) beyond slide_id and sample_id
    if not par["metadata"]:
        return {}

    logger.info("Reading sample metadata from '%s'...", par["metadata"])
    df = pd.read_csv(par["metadata"], dtype=str)
    assert "sample_id" in df.columns, "--metadata must contain a 'sample_id' column."
    df = df.set_index("sample_id")

    missing = set(par["sample_id"]) - set(df.index)
    assert not missing, f"--metadata has no rows for sample IDs: {sorted(missing)}"

    return df.to_dict(orient="index")


def main():
    check_arguments()
    sample_meta = read_sample_metadata()

    exp = InSituExperiment()

    for input_path, sample_id in zip(par["input"], par["sample_id"]):
        logger.info("Reading InSituPy project '%s' (%s)...", input_path, sample_id)
        xd = InSituData.read(Path(input_path))

        # crop() returns a new, detached dataset, so it must happen before add()
        if par["crop"]:
            logger.info(
                "Cropping '%s' to xlim=%s, ylim=%s...",
                sample_id,
                par["xlim"],
                par["ylim"],
            )
            xd = xd.crop(xlim=tuple(par["xlim"]), ylim=tuple(par["ylim"]))

        exp.add(
            data=xd,
            metadata={
                "slide_id": xd.slide_id,
                "sample_id": sample_id,
                **sample_meta.get(sample_id, {}),
            },
        )

    logger.info(
        "Saving InSituExperiment with %d samples to '%s'...",
        len(par["input"]),
        par["output"],
    )
    exp.saveas(par["output"], verbose=False)

    logger.info("Done!")


if __name__ == "__main__":
    main()
