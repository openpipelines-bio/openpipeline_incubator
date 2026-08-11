import csv
import re
import sys
from pathlib import Path

import pandas as pd

### VIASH START
par = {
    "input": "fastq/",
    "sample_sheet": None,
    "include_undetermined": False,
    "output": "fastq_manifest.csv",
}
### VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()

# Matches bcl-convert and bases2fastq FASTQ naming. The lane group is optional
# because bcl-convert drops it under NoLaneSplitting. Both .fastq.gz and .fq.gz
# are accepted.
FASTQ_RE = re.compile(
    r"^(?P<sample>[A-Za-z0-9\-_.]+)_S(?P<snum>\d+)(_L(?P<lane>\d+))?_(?P<read>[RI]\d)_(?P<chunk>\d+)\.f(ast)?q\.gz$"
)

UNDETERMINED_SAMPLE = "Undetermined"
# Sort key used for the (missing) blank lane, so no-lane-splitting files sort before lane 1, 2, ...
NO_LANE_SORT_KEY = -1
# read_type sort order: R1, R2, I1, I2 (genomic reads before index reads), not lexical order.
READ_TYPE_ORDER = {"R1": 0, "R2": 1, "I1": 2, "I2": 3}


def find_fastq_files(input_dir):
    """Recursively find every file under input_dir."""
    return sorted(p for p in Path(input_dir).rglob("*") if p.is_file())


def parse_filename(path):
    """Parse a FASTQ filename against FASTQ_RE. Returns None if it doesn't match."""
    match = FASTQ_RE.match(path.name)
    if match is None:
        return None

    sample = match.group("sample")
    lane = match.group("lane")
    read = match.group("read")

    return {
        "sample_id": sample,
        # Leading zeros are dropped (L001 -> "1"). Blank when NoLaneSplitting
        # was used, i.e. there is no _L00N_ component in the filename.
        "lane": str(int(lane)) if lane is not None else "",
        "read_type": read,
        "filename": path.name,
        "path": str(path.resolve()),
    }


def parse_sample_sheet_samples(sample_sheet_path):
    """Parse an Illumina-format sample sheet's [*_Data] section."""
    with open(sample_sheet_path, newline="") as f:
        lines = [line.rstrip("\n").rstrip("\r") for line in f]

    data_section_lines = []
    in_data_section = False
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            section_name = stripped[1:-1]
            in_data_section = section_name.endswith("_Data") or section_name == "Data"
            continue
        if in_data_section and stripped != "":
            data_section_lines.append(line)

    if not data_section_lines:
        logger.warning(
            "No '[*_Data]' section found in sample sheet '%s'; skipping cross-validation.",
            sample_sheet_path,
        )
        return None

    reader = csv.DictReader(data_section_lines)
    header = reader.fieldnames or []
    id_column = None
    for candidate in ("Sample_ID", "Sample_Name"):
        if candidate in header:
            id_column = candidate
            break

    if id_column is None:
        logger.warning(
            "No 'Sample_ID' or 'Sample_Name' column found in sample sheet '%s'; skipping "
            "cross-validation.",
            sample_sheet_path,
        )
        return None

    samples = {row[id_column] for row in reader if row.get(id_column)}
    return samples


def main():
    input_dir = Path(par["input"])
    if not input_dir.is_dir():
        raise NotADirectoryError(f"Input directory '{input_dir}' does not exist.")

    files = find_fastq_files(input_dir)

    rows = []
    for path in files:
        parsed = parse_filename(path)
        if parsed is None:
            logger.debug(
                "Skipping file that does not match the FASTQ naming pattern: %s", path
            )
            continue

        is_undetermined = parsed[
            "sample_id"
        ] == "Undetermined" and path.name.startswith("Undetermined_S0_")
        if is_undetermined and not par["include_undetermined"]:
            logger.debug(
                "Excluding undetermined file (--include_undetermined is false): %s",
                path,
            )
            continue
        if is_undetermined:
            parsed["sample_id"] = UNDETERMINED_SAMPLE

        if path.stat().st_size == 0:
            logger.warning("FASTQ file is empty: %s", path)

        rows.append(parsed)

    if not rows:
        raise ValueError(f"No FASTQ files found under '{input_dir}'.")

    manifest = pd.DataFrame(
        rows, columns=["sample_id", "lane", "read_type", "filename", "path"]
    )

    # Validate: every sample needs an R1, R1/R2 counts must match, and no R2 may
    # be unpaired.
    # Lane is part of the pairing key so multi-lane R1/R2 counts are compared
    # per lane.
    errors = []
    for sample_id, sample_rows in manifest.groupby("sample_id"):
        r1_by_lane = sample_rows.loc[
            sample_rows["read_type"] == "R1", "lane"
        ].value_counts()
        r2_by_lane = sample_rows.loc[
            sample_rows["read_type"] == "R2", "lane"
        ].value_counts()

        if r1_by_lane.sum() == 0:
            errors.append(f"Sample '{sample_id}' has no R1 files.")
            continue

        if r2_by_lane.sum() > 0:
            if r1_by_lane.sum() != r2_by_lane.sum():
                errors.append(
                    f"Sample '{sample_id}' has {r1_by_lane.sum()} R1 file(s) but "
                    f"{r2_by_lane.sum()} R2 file(s)."
                )
            unpaired_lanes = set(r2_by_lane.index) - set(r1_by_lane.index)
            if unpaired_lanes:
                errors.append(
                    f"Sample '{sample_id}' has R2 file(s) in lane(s) "
                    f"{sorted(unpaired_lanes)} with no matching R1."
                )

    if errors:
        raise ValueError("FASTQ manifest validation failed:\n" + "\n".join(errors))

    if par.get("sample_sheet"):
        sheet_samples = parse_sample_sheet_samples(par["sample_sheet"])
        if sheet_samples is not None:
            found_samples = set(manifest["sample_id"]) - {UNDETERMINED_SAMPLE}
            if found_samples != sheet_samples:
                missing_from_sheet = found_samples - sheet_samples
                missing_from_fastq = sheet_samples - found_samples
                raise ValueError(
                    "Sample names found in the FASTQ directory do not match the sample "
                    "sheet.\n"
                    f"In FASTQ directory but not in sample sheet: {sorted(missing_from_sheet)}\n"
                    f"In sample sheet but not in FASTQ directory: {sorted(missing_from_fastq)}"
                )

    # Sort by sample_id, then lane (numerically, blank lane first), then
    # read_type in the fixed order R1, R2, I1, I2 (genomic reads before index
    # reads).
    manifest["__lane_sort__"] = manifest["lane"].apply(
        lambda x: int(x) if x != "" else NO_LANE_SORT_KEY
    )
    manifest["__read_type_sort__"] = manifest["read_type"].map(READ_TYPE_ORDER)
    manifest = (
        manifest.sort_values(
            by=["sample_id", "__lane_sort__", "__read_type_sort__"], kind="stable"
        )
        .drop(columns=["__lane_sort__", "__read_type_sort__"])
        .reset_index(drop=True)
    )

    logger.info("Writing manifest with %d row(s) to %s", len(manifest), par["output"])
    manifest.to_csv(par["output"], index=False)


if __name__ == "__main__":
    main()
