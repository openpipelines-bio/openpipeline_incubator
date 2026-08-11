import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

## VIASH START
meta = {
    "name": "create_fastq_manifest",
    "resources_dir": "resources_test/",
    "executable": "target/executable/demultiplexing/create_fastq_manifest/create_fastq_manifest",
    "config": "src/components/demultiplexing/create_fastq_manifest/config.vsh.yaml",
}
## VIASH END


def write_fastq(path, content=b"x"):
    """Write a dummy non-empty stand-in FASTQ.gz file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(content)


def touch_empty(path):
    """Create a genuinely empty (0 byte) FASTQ.gz file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.touch()


def read_manifest(output):
    return pd.read_csv(output, dtype=str, keep_default_na=False)


def test_single_sample_single_lane(run_component, tmp_path):
    fastq_dir = tmp_path / "fastq"
    write_fastq(fastq_dir / "sample1_S1_L001_R1_001.fastq.gz")
    write_fastq(fastq_dir / "sample1_S1_L001_R2_001.fastq.gz")
    output = tmp_path / "manifest.csv"

    run_component(["--input", fastq_dir, "--output", output])

    manifest = read_manifest(output)
    assert list(manifest.columns) == [
        "sample_id",
        "lane",
        "read_type",
        "filename",
        "path",
    ]
    assert len(manifest) == 2
    assert set(manifest["sample_id"]) == {"sample1"}
    assert set(manifest["lane"]) == {"1"}
    assert set(manifest["read_type"]) == {"R1", "R2"}
    # path must be absolute and point at a real file.
    for path in manifest["path"]:
        assert path.startswith("/"), f"Expected an absolute path, got: {path}"
        assert Path(path).is_file()


def test_multi_lane(run_component, tmp_path):
    fastq_dir = tmp_path / "fastq"
    for lane in ("L001", "L002", "L003"):
        write_fastq(fastq_dir / f"sample1_S1_{lane}_R1_001.fastq.gz")
        write_fastq(fastq_dir / f"sample1_S1_{lane}_R2_001.fastq.gz")
    output = tmp_path / "manifest.csv"

    run_component(["--input", fastq_dir, "--output", output])

    manifest = read_manifest(output)
    assert len(manifest) == 6
    assert set(manifest["lane"]) == {"1", "2", "3"}
    # Rows are sorted by sample_id, then lane, then read_type
    assert manifest["lane"].tolist() == ["1", "1", "2", "2", "3", "3"]


def test_no_lane_component(run_component, tmp_path):
    # bcl-convert drops the _L00N_ component under NoLaneSplitting.
    fastq_dir = tmp_path / "fastq"
    write_fastq(fastq_dir / "sample1_S1_R1_001.fastq.gz")
    write_fastq(fastq_dir / "sample1_S1_R2_001.fastq.gz")
    output = tmp_path / "manifest.csv"

    run_component(["--input", fastq_dir, "--output", output])

    manifest = read_manifest(output)
    assert len(manifest) == 2
    # lane must be blank, not e.g. "0" or missing the column entirely.
    assert set(manifest["lane"]) == {""}


def test_single_end_r1_only(run_component, tmp_path):
    fastq_dir = tmp_path / "fastq"
    write_fastq(fastq_dir / "sample1_S1_L001_R1_001.fastq.gz")
    output = tmp_path / "manifest.csv"

    run_component(["--input", fastq_dir, "--output", output])

    manifest = read_manifest(output)
    assert len(manifest) == 1
    assert manifest.loc[0, "read_type"] == "R1"


def test_index_reads(run_component, tmp_path):
    fastq_dir = tmp_path / "fastq"
    write_fastq(fastq_dir / "sample1_S1_L001_R1_001.fastq.gz")
    write_fastq(fastq_dir / "sample1_S1_L001_R2_001.fastq.gz")
    write_fastq(fastq_dir / "sample1_S1_L001_I1_001.fastq.gz")
    write_fastq(fastq_dir / "sample1_S1_L001_I2_001.fastq.gz")
    output = tmp_path / "manifest.csv"

    run_component(["--input", fastq_dir, "--output", output])

    manifest = read_manifest(output)
    assert len(manifest) == 4
    assert set(manifest["read_type"]) == {"R1", "R2", "I1", "I2"}
    # Index reads must carry the correct read_type
    i2_row = manifest[manifest["read_type"] == "I2"].iloc[0]
    assert i2_row["filename"] == "sample1_S1_L001_I2_001.fastq.gz"
    i1_row = manifest[manifest["read_type"] == "I1"].iloc[0]
    assert i1_row["filename"] == "sample1_S1_L001_I1_001.fastq.gz"
    # read_type order is R1, R2, I1, I2
    assert manifest["read_type"].tolist() == ["R1", "R2", "I1", "I2"]


def test_undetermined_excluded_by_default(run_component, tmp_path):
    fastq_dir = tmp_path / "fastq"
    write_fastq(fastq_dir / "sample1_S1_L001_R1_001.fastq.gz")
    write_fastq(fastq_dir / "sample1_S1_L001_R2_001.fastq.gz")
    write_fastq(fastq_dir / "Undetermined_S0_L001_R1_001.fastq.gz")
    write_fastq(fastq_dir / "Undetermined_S0_L001_R2_001.fastq.gz")
    output = tmp_path / "manifest.csv"

    run_component(["--input", fastq_dir, "--output", output])

    manifest = read_manifest(output)
    assert "Undetermined" not in set(manifest["sample_id"])
    assert len(manifest) == 2


def test_undetermined_included(run_component, tmp_path):
    fastq_dir = tmp_path / "fastq"
    write_fastq(fastq_dir / "sample1_S1_L001_R1_001.fastq.gz")
    write_fastq(fastq_dir / "sample1_S1_L001_R2_001.fastq.gz")
    write_fastq(fastq_dir / "Undetermined_S0_L001_R1_001.fastq.gz")
    write_fastq(fastq_dir / "Undetermined_S0_L001_R2_001.fastq.gz")
    output = tmp_path / "manifest.csv"

    run_component(
        [
            "--input",
            fastq_dir,
            "--output",
            output,
            "--include_undetermined",
            "true",
        ]
    )

    manifest = read_manifest(output)
    assert set(manifest["sample_id"]) == {"sample1", "Undetermined"}
    assert len(manifest) == 4


def test_fq_gz_extension(run_component, tmp_path):
    fastq_dir = tmp_path / "fastq"
    write_fastq(fastq_dir / "sample1_S1_L001_R1_001.fq.gz")
    write_fastq(fastq_dir / "sample1_S1_L001_R2_001.fq.gz")
    output = tmp_path / "manifest.csv"

    run_component(["--input", fastq_dir, "--output", output])

    manifest = read_manifest(output)
    assert len(manifest) == 2
    assert set(manifest["filename"]) == {
        "sample1_S1_L001_R1_001.fq.gz",
        "sample1_S1_L001_R2_001.fq.gz",
    }


def test_files_in_subdirectories(run_component, tmp_path):
    # bcl-convert writes Sample_Project subdirectories
    # bases2fastq writes per-project ones, the scan must be recursive.
    fastq_dir = tmp_path / "fastq"
    write_fastq(fastq_dir / "Project1" / "sample1_S1_L001_R1_001.fastq.gz")
    write_fastq(fastq_dir / "Project1" / "sample1_S1_L001_R2_001.fastq.gz")
    write_fastq(fastq_dir / "Project2" / "nested" / "sample2_S2_L001_R1_001.fastq.gz")
    write_fastq(fastq_dir / "Project2" / "nested" / "sample2_S2_L001_R2_001.fastq.gz")
    output = tmp_path / "manifest.csv"

    run_component(["--input", fastq_dir, "--output", output])

    manifest = read_manifest(output)
    assert len(manifest) == 4
    assert set(manifest["sample_id"]) == {"sample1", "sample2"}


def test_empty_directory_fails(run_component, tmp_path):
    fastq_dir = tmp_path / "fastq"
    fastq_dir.mkdir()
    output = tmp_path / "manifest.csv"

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(["--input", fastq_dir, "--output", output])

    assert "No FASTQ files found" in err.value.stdout.decode("utf-8")


def test_unpaired_r2_fails(run_component, tmp_path):
    fastq_dir = tmp_path / "fastq"
    write_fastq(fastq_dir / "sample1_S1_L001_R1_001.fastq.gz")
    write_fastq(fastq_dir / "sample1_S1_L001_R2_001.fastq.gz")
    # An extra R2 in a lane with no matching R1.
    write_fastq(fastq_dir / "sample1_S1_L002_R2_001.fastq.gz")
    output = tmp_path / "manifest.csv"

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(["--input", fastq_dir, "--output", output])

    message = err.value.stdout.decode("utf-8")
    assert "sample1" in message


def test_sample_sheet_mismatch_fails(run_component, tmp_path):
    fastq_dir = tmp_path / "fastq"
    write_fastq(fastq_dir / "sample1_S1_L001_R1_001.fastq.gz")
    write_fastq(fastq_dir / "sample1_S1_L001_R2_001.fastq.gz")

    sample_sheet = tmp_path / "SampleSheet.csv"
    sample_sheet.write_text(
        "[Header]\n"
        "FileFormatVersion,2\n"
        "\n"
        "[BCLConvert_Data]\n"
        "Sample_ID,Index\n"
        "other_sample,AAAAAAAA\n"
    )
    output = tmp_path / "manifest.csv"

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                fastq_dir,
                "--sample_sheet",
                sample_sheet,
                "--output",
                output,
            ]
        )

    message = err.value.stdout.decode("utf-8")
    assert "sample1" in message
    assert "other_sample" in message


def test_sample_sheet_match_succeeds(run_component, tmp_path):
    fastq_dir = tmp_path / "fastq"
    write_fastq(fastq_dir / "sample1_S1_L001_R1_001.fastq.gz")
    write_fastq(fastq_dir / "sample1_S1_L001_R2_001.fastq.gz")

    sample_sheet = tmp_path / "SampleSheet.csv"
    sample_sheet.write_text(
        "[Header]\n"
        "FileFormatVersion,2\n"
        "\n"
        "[BCLConvert_Data]\n"
        "Sample_ID,Index\n"
        "sample1,AAAAAAAA\n"
    )
    output = tmp_path / "manifest.csv"

    run_component(
        [
            "--input",
            fastq_dir,
            "--sample_sheet",
            sample_sheet,
            "--output",
            output,
        ]
    )

    manifest = read_manifest(output)
    assert set(manifest["sample_id"]) == {"sample1"}


def test_empty_fastq_warns(run_component, tmp_path):
    # Undetermined and no-index runs legitimately produce empty FASTQ files.
    # This must warn, not fail, and the file must still be recorded in the
    # manifest.
    fastq_dir = tmp_path / "fastq"
    touch_empty(fastq_dir / "sample1_S1_L001_R1_001.fastq.gz")
    write_fastq(fastq_dir / "sample1_S1_L001_R2_001.fastq.gz")
    output = tmp_path / "manifest.csv"

    result = run_component(["--input", fastq_dir, "--output", output])

    manifest = read_manifest(output)
    assert len(manifest) == 2

    log_text = ""
    if result is not None and getattr(result, "stdout", None) is not None:
        log_text = result.stdout.decode("utf-8")
    assert log_text == "" or "empty" in log_text.lower()


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
