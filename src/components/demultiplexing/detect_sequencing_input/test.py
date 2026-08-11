import json
import subprocess
import sys

import pytest

## VIASH START
meta = {
    "name": "detect_sequencing_input",
    "resources_dir": "resources_test/",
    "executable": "target/executable/components/demultiplexing/detect_sequencing_input/detect_sequencing_input",
    "config": "src/components/demultiplexing/detect_sequencing_input/config.vsh.yaml",
}
## VIASH END

RUN_INFO_READS = [
    {
        "number": 1,
        "num_cycles": 50,
        "is_indexed": False,
        "is_reverse_complement": False,
    },
    {"number": 2, "num_cycles": 8, "is_indexed": True, "is_reverse_complement": False},
    {"number": 3, "num_cycles": 24, "is_indexed": True, "is_reverse_complement": True},
    {
        "number": 4,
        "num_cycles": 49,
        "is_indexed": False,
        "is_reverse_complement": False,
    },
]


def run_info_xml(reads, instrument="A00834", shuffle=False):
    """Hand-written minimal RunInfo.xml content for the given read structure."""
    ordered_reads = list(reversed(reads)) if shuffle else reads
    read_lines = []
    for read in ordered_reads:
        attrs = (
            f'Number="{read["number"]}" NumCycles="{read["num_cycles"]}" '
            f'IsIndexedRead="{"Y" if read["is_indexed"] else "N"}"'
        )
        if read.get("is_reverse_complement"):
            attrs += ' IsReverseComplement="Y"'
        read_lines.append(f"      <Read {attrs} />")
    reads_block = "\n".join(read_lines)
    instrument_line = f"    <Instrument>{instrument}</Instrument>" if instrument else ""
    return f"""<?xml version="1.0"?>
<RunInfo>
  <Run Id="220101_A00834_0001_AHNCF7DSXX" Number="1">
    <Flowcell>HNCF7DSXX</Flowcell>
{instrument_line}
    <Date>220101</Date>
    <Reads>
{reads_block}
    </Reads>
  </Run>
</RunInfo>
"""


@pytest.fixture
def illumina_run_folder(tmp_path):
    run_dir = tmp_path / "illumina_run"
    run_dir.mkdir()
    (run_dir / "RunInfo.xml").write_text(run_info_xml(RUN_INFO_READS))
    return run_dir


@pytest.fixture
def element_run_folder(tmp_path):
    run_dir = tmp_path / "element_run"
    run_dir.mkdir()
    (run_dir / "RunParameters.json").write_text(json.dumps({"RunID": "20230404"}))
    return run_dir


def touch_fastq_files(directory, filenames):
    directory.mkdir(parents=True, exist_ok=True)
    for filename in filenames:
        (directory / filename).write_bytes(b"")


def read_output(output_path):
    with open(output_path) as f:
        return json.load(f)


def test_illumina_run_folder(run_component, illumina_run_folder, tmp_path):
    output = tmp_path / "detection.json"
    run_component(["--input", illumina_run_folder, "--output", output])

    result = read_output(output)
    assert result["input_type"] == "illumina"
    assert result["demultiplexer"] == "bclconvert"
    assert result["sample_sheet"] is None
    assert result["instrument"] == "A00834"
    assert result["read_structure"] == RUN_INFO_READS


def test_illumina_run_folder_read_structure_ordered_by_number(run_component, tmp_path):
    # Reads are written out of order in the XML, the output must still be ordered by number.
    run_dir = tmp_path / "illumina_run_shuffled"
    run_dir.mkdir()
    (run_dir / "RunInfo.xml").write_text(run_info_xml(RUN_INFO_READS, shuffle=True))

    output = tmp_path / "detection.json"
    run_component(["--input", run_dir, "--output", output])

    result = read_output(output)
    assert [r["number"] for r in result["read_structure"]] == [1, 2, 3, 4]
    assert result["read_structure"] == RUN_INFO_READS


def test_instrument_extraction(run_component, tmp_path):
    run_dir = tmp_path / "illumina_run_instrument"
    run_dir.mkdir()
    (run_dir / "RunInfo.xml").write_text(
        run_info_xml(RUN_INFO_READS, instrument="LH00123")
    )

    output = tmp_path / "detection.json"
    run_component(["--input", run_dir, "--output", output])

    result = read_output(output)
    assert result["instrument"] == "LH00123"


def test_instrument_absent_is_null(run_component, tmp_path):
    run_dir = tmp_path / "illumina_run_no_instrument"
    run_dir.mkdir()
    (run_dir / "RunInfo.xml").write_text(run_info_xml(RUN_INFO_READS, instrument=None))

    output = tmp_path / "detection.json"
    run_component(["--input", run_dir, "--output", output])

    result = read_output(output)
    assert result["instrument"] is None


@pytest.mark.parametrize("marker_filename", ["RunParameters.json", "RunUploaded.json"])
def test_element_run_folder(run_component, tmp_path, marker_filename):
    run_dir = tmp_path / f"element_run_{marker_filename}"
    run_dir.mkdir()
    (run_dir / marker_filename).write_text(json.dumps({"outcome": "OutcomeCompleted"}))

    output = tmp_path / "detection.json"
    run_component(["--input", run_dir, "--output", output])

    result = read_output(output)
    assert result["input_type"] == "element"
    assert result["demultiplexer"] == "bases2fastq"
    assert result["instrument"] is None
    assert result["read_structure"] == []


def test_fastq_only_directory(run_component, tmp_path):
    run_dir = tmp_path / "fastq_dir"
    touch_fastq_files(
        run_dir,
        [
            "test_sample_S1_L001_R1_001.fastq.gz",
            "test_sample_S1_L001_R2_001.fastq.gz",
        ],
    )

    output = tmp_path / "detection.json"
    run_component(["--input", run_dir, "--output", output])

    result = read_output(output)
    assert result["input_type"] == "fastq"
    assert result["demultiplexer"] is None
    assert result["instrument"] is None
    assert result["read_structure"] == []


def test_fastq_in_subdirectories(run_component, tmp_path):
    run_dir = tmp_path / "fastq_nested"
    touch_fastq_files(
        run_dir / "Project1" / "Sample1",
        ["Sample1_S1_L001_R1_001.fastq.gz", "Sample1_S1_L001_R2_001.fastq.gz"],
    )

    output = tmp_path / "detection.json"
    run_component(["--input", run_dir, "--output", output])

    result = read_output(output)
    assert result["input_type"] == "fastq"
    assert result["demultiplexer"] is None


def test_empty_directory_fails(run_component, tmp_path):
    run_dir = tmp_path / "empty_dir"
    run_dir.mkdir()

    output = tmp_path / "detection.json"
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(["--input", run_dir, "--output", output])

    stdout = err.value.stdout.decode("utf-8")
    assert "Could not detect the sequencing input type" in stdout
    assert "empty directory" in stdout


def test_both_illumina_and_element_markers_fails(run_component, tmp_path):
    run_dir = tmp_path / "conflicting_run"
    run_dir.mkdir()
    (run_dir / "RunInfo.xml").write_text(run_info_xml(RUN_INFO_READS))
    (run_dir / "RunParameters.json").write_text(json.dumps({"RunID": "20230404"}))

    output = tmp_path / "detection.json"
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(["--input", run_dir, "--output", output])

    stdout = err.value.stdout.decode("utf-8")
    assert "RunInfo.xml" in stdout
    assert "RunParameters.json" in stdout


def test_lowercase_run_parameters_xml_still_detected_as_illumina(
    run_component, tmp_path
):
    run_dir = tmp_path / "illumina_run_lowercase_params"
    run_dir.mkdir()
    (run_dir / "RunInfo.xml").write_text(run_info_xml(RUN_INFO_READS))
    (run_dir / "runParameters.xml").write_text("<RunParameters></RunParameters>")

    output = tmp_path / "detection.json"
    run_component(["--input", run_dir, "--output", output])

    result = read_output(output)
    assert result["input_type"] == "illumina"
    assert result["demultiplexer"] == "bclconvert"


def test_demultiplexer_override_agreeing_with_detection(
    run_component, illumina_run_folder, tmp_path
):
    output = tmp_path / "detection.json"
    run_component(
        [
            "--input",
            illumina_run_folder,
            "--demultiplexer",
            "bclconvert",
            "--output",
            output,
        ]
    )

    result = read_output(output)
    assert result["input_type"] == "illumina"
    assert result["demultiplexer"] == "bclconvert"


def test_demultiplexer_override_disagreeing_with_detection_fails(
    run_component, illumina_run_folder, tmp_path
):
    output = tmp_path / "detection.json"
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                illumina_run_folder,
                "--demultiplexer",
                "bases2fastq",
                "--output",
                output,
            ]
        )

    stdout = err.value.stdout.decode("utf-8")
    assert "bases2fastq" in stdout
    assert "bclconvert" in stdout


def test_skip_demultiplexing_forces_fastq_passthrough(
    run_component, illumina_run_folder, tmp_path
):
    output = tmp_path / "detection.json"
    run_component(
        ["--input", illumina_run_folder, "--skip_demultiplexing", "--output", output]
    )

    result = read_output(output)
    assert result["input_type"] == "fastq"
    assert result["demultiplexer"] is None
    assert result["read_structure"] == []
    assert result["instrument"] is None


def test_sample_sheet_found_in_input_directory(
    run_component, illumina_run_folder, tmp_path
):
    sample_sheet = illumina_run_folder / "SampleSheet.csv"
    sample_sheet.write_text("[Header]\n")

    output = tmp_path / "detection.json"
    run_component(["--input", illumina_run_folder, "--output", output])

    result = read_output(output)
    assert result["sample_sheet"] == str(sample_sheet.resolve())


def test_sample_sheet_override(run_component, illumina_run_folder, tmp_path):
    sample_sheet = tmp_path / "custom_sheet.csv"
    sample_sheet.write_text("[Header]\n")

    output = tmp_path / "detection.json"
    run_component(
        [
            "--input",
            illumina_run_folder,
            "--sample_sheet",
            sample_sheet,
            "--output",
            output,
        ]
    )

    result = read_output(output)
    assert result["sample_sheet"] == str(sample_sheet.resolve())


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
