import sys
from pathlib import Path

import pandas as pd
import pytest

## VIASH START
par = {
    "output_fastq": "fastq/",
    "output_fastq_manifest": "fastq_manifest.csv",
}
## VIASH END

# Columns create_fastq_manifest documents and always writes.
EXPECTED_COLUMNS = {"sample_id", "lane", "read_type", "filename", "path"}


def test_run():
    manifest_path = Path(par["output_fastq_manifest"])
    assert manifest_path.is_file(), f"Manifest '{manifest_path}' should be a file."

    output_fastq_dir = Path(par["output_fastq"]).resolve()
    assert output_fastq_dir.is_dir(), (
        f"'{output_fastq_dir}' (--output_fastq) should be a directory."
    )

    manifest = pd.read_csv(manifest_path, dtype=str).fillna("")
    assert not manifest.empty, f"Manifest '{manifest_path}' should not be empty."

    # Documented columns must be present
    assert EXPECTED_COLUMNS.issubset(manifest.columns), (
        f"Manifest should contain columns {sorted(EXPECTED_COLUMNS)}, "
        f"found: {manifest.columns.to_list()}"
    )

    # Sample ids must be non-empty for every row
    blank_rows = manifest[manifest["sample_id"].str.strip() == ""]
    assert blank_rows.empty, (
        f"All manifest rows should have a non-empty sample_id. Blank rows:\n{blank_rows}"
    )

    # Every referenced FASTQ file must actually exist under --output_fastq, be
    # non-empty, and agree with its own 'filename' column.
    for row in manifest.to_dict("records"):
        fastq_file = Path(row["path"])
        assert fastq_file.is_file(), (
            f"Manifest row for sample '{row['sample_id']}' references '{fastq_file}', "
            "which does not exist."
        )
        assert fastq_file.stat().st_size > 0, (
            f"FASTQ file '{fastq_file}' for sample '{row['sample_id']}' should not be empty."
        )
        assert fastq_file.name == row["filename"], (
            f"Manifest 'filename' ({row['filename']}) should match the basename of "
            f"'path' ({fastq_file.name})."
        )
        assert fastq_file.resolve().is_relative_to(output_fastq_dir), (
            f"'{fastq_file}' should live under the reported --output_fastq directory "
            f"'{output_fastq_dir}'."
        )

    # R1/R2 pairing: every (sample_id, lane) group needs at least one R1 file, and
    # if R2 files are present their count must match R1's for that group.
    read_type_counts = manifest.pivot_table(
        index=["sample_id", "lane"], columns="read_type", aggfunc="size", fill_value=0
    )
    assert "R1" in read_type_counts.columns, (
        "Manifest should contain R1 rows for every sample."
    )
    r1_counts = read_type_counts["R1"]
    assert (r1_counts > 0).all(), (
        f"Every (sample_id, lane) group should have at least one R1 file:\n{r1_counts}"
    )
    if "R2" in read_type_counts.columns:
        r2_counts = read_type_counts["R2"]
        mismatched = read_type_counts[(r2_counts > 0) & (r1_counts != r2_counts)]
        assert mismatched.empty, (
            f"R1/R2 counts should match per (sample_id, lane) group:\n{mismatched}"
        )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "--import-mode=importlib"]))
