import json
import re
import subprocess
import sys

import pytest

## VIASH START
meta = {
    "name": "apply_bclconvert_preset",
    "resources_dir": "resources_test/",
    "executable": "target/executable/demultiplexing/apply_bclconvert_preset/apply_bclconvert_preset",
    "config": "src/components/demultiplexing/apply_bclconvert_preset/config.vsh.yaml",
}
## VIASH END


# A minimal but realistic sample sheet:
BASE_SAMPLE_SHEET = """[Header]
FileFormatVersion,2
[Reads]
Read1Cycles,50
Read2Cycles,49
Index1Cycles,8
Index2Cycles,24
[BCLConvert_Settings]
SoftwareVersion,4.2.7
[BCLConvert_Data]
Sample_ID,index,index2
sample1,AAAAAAAA,CCCCCCCCCCCCCCCCCCCCCCCC
"""


def make_read(number, num_cycles, is_indexed, is_reverse_complement=False):
    return {
        "number": number,
        "num_cycles": num_cycles,
        "is_indexed": is_indexed,
        "is_reverse_complement": is_reverse_complement,
    }


def read_structure(r1=50, i1=8, i2=24, r2=49):
    """The read_structure of the schema example: R1, I1, I2 (reverse complement), R2."""
    return [
        make_read(1, r1, False),
        make_read(2, i1, True),
        make_read(3, i2, True, is_reverse_complement=True),
        make_read(4, r2, False),
    ]


def detection_dict(structure, **overrides):
    base = {
        "input_type": "illumina",
        "demultiplexer": "bclconvert",
        "sample_sheet": "/path/to/SampleSheet.csv",
        "instrument": "A00834",
        "read_structure": structure,
    }
    base.update(overrides)
    return base


@pytest.fixture
def write_sample_sheet(tmp_path):
    def wrapper(text, name="SampleSheet.csv"):
        path = tmp_path / name
        path.write_text(text)
        return path

    return wrapper


@pytest.fixture
def write_detection_json(tmp_path):
    def wrapper(structure, name="detection.json", **overrides):
        path = tmp_path / name
        path.write_text(json.dumps(detection_dict(structure, **overrides)))
        return path

    return wrapper


def section_names(text):
    """The ordered list of '[SectionName]' headers appearing in a sample sheet."""
    return re.findall(r"^\[(.+)\]$", text, flags=re.MULTILINE)


def section_lines(text, name):
    """The raw content lines of one section (excluding its '[Name]' header line)."""
    lines = text.splitlines()
    start = None
    for i, line in enumerate(lines):
        if line.strip() == f"[{name}]":
            start = i + 1
            break

    if start is None:
        return None

    end = len(lines)
    for i in range(start, len(lines)):
        if re.match(r"^\[.+\]$", lines[i].strip()):
            end = i
            break

    return lines[start:end]


def parse_settings(text):
    """Parse the [BCLConvert_Settings] section into a dict, for easy assertions."""
    lines = section_lines(text, "BCLConvert_Settings") or []
    settings = {}
    for line in lines:
        if not line.strip():
            continue
        key, _, value = line.partition(",")
        settings[key] = value

    return settings


def test_override_cycles_u16_single_cell_atac(
    run_component, write_sample_sheet, write_detection_json, tmp_path
):
    sample_sheet = write_sample_sheet(BASE_SAMPLE_SHEET)
    detection_json = write_detection_json(read_structure(r1=50, i1=8, i2=16, r2=49))
    output = tmp_path / "output.csv"

    run_component(
        [
            "--sample_sheet",
            sample_sheet,
            "--detection_json",
            detection_json,
            "--preset",
            "10x_atac",
            "--output",
            output,
        ]
    )

    settings = parse_settings(output.read_text())
    assert settings["OverrideCycles"] == "R1:Y50;I1:I8;I2:U16;R2:Y49"
    assert settings["CreateFastqForIndexReads"] == "1"
    assert settings["TrimUMI"] == "0"


def test_override_cycles_u24_multiome_atac(
    run_component, write_sample_sheet, write_detection_json, tmp_path
):
    sample_sheet = write_sample_sheet(BASE_SAMPLE_SHEET)
    detection_json = write_detection_json(read_structure(r1=50, i1=8, i2=24, r2=49))
    output = tmp_path / "output.csv"

    run_component(
        [
            "--sample_sheet",
            sample_sheet,
            "--detection_json",
            detection_json,
            "--preset",
            "10x_atac",
            "--output",
            output,
        ]
    )

    settings = parse_settings(output.read_text())
    assert settings["OverrideCycles"] == "R1:Y50;I1:I8;I2:U24;R2:Y49"


def test_override_cycles_read_structure_order_independent(
    run_component, write_sample_sheet, write_detection_json, tmp_path
):
    scrambled = [
        make_read(4, 49, False),
        make_read(2, 10, True),
        make_read(1, 51, False),
        make_read(3, 16, True, is_reverse_complement=True),
    ]
    sample_sheet = write_sample_sheet(BASE_SAMPLE_SHEET)
    detection_json = write_detection_json(scrambled)
    output = tmp_path / "output.csv"

    run_component(
        [
            "--sample_sheet",
            sample_sheet,
            "--detection_json",
            detection_json,
            "--preset",
            "10x_atac",
            "--output",
            output,
        ]
    )

    settings = parse_settings(output.read_text())
    assert settings["OverrideCycles"] == "R1:Y51;I1:I10;I2:U16;R2:Y49"


def test_creates_missing_bclconvert_settings_section(
    run_component, write_sample_sheet, write_detection_json, tmp_path
):
    sheet_without_settings = (
        "[Header]\n"
        "FileFormatVersion,2\n"
        "[BCLConvert_Data]\n"
        "Sample_ID,index,index2\n"
        "sample1,AAAAAAAA,CCCCCCCCCCCCCCCCCCCCCCCC\n"
    )
    sample_sheet = write_sample_sheet(sheet_without_settings)
    detection_json = write_detection_json(read_structure(i2=16))
    output = tmp_path / "output.csv"

    run_component(
        [
            "--sample_sheet",
            sample_sheet,
            "--detection_json",
            detection_json,
            "--preset",
            "10x_atac",
            "--output",
            output,
        ]
    )

    text = output.read_text()
    settings = parse_settings(text)
    assert settings["CreateFastqForIndexReads"] == "1"
    assert settings["TrimUMI"] == "0"
    assert settings["OverrideCycles"] == "R1:Y50;I1:I8;I2:U16;R2:Y49"
    names = section_names(text)
    assert names == ["Header", "BCLConvert_Settings", "BCLConvert_Data"]


def test_overwrites_conflicting_existing_values(
    run_component, write_sample_sheet, write_detection_json, tmp_path
):
    sheet_with_conflicts = (
        "[Header]\n"
        "FileFormatVersion,2\n"
        "[BCLConvert_Settings]\n"
        "CreateFastqForIndexReads,0\n"
        "TrimUMI,1\n"
        "OverrideCycles,R1:Y999;I1:I999;I2:I999;R2:Y999\n"
        "[BCLConvert_Data]\n"
        "Sample_ID,index,index2\n"
        "sample1,AAAAAAAA,CCCCCCCCCCCCCCCCCCCCCCCC\n"
    )
    sample_sheet = write_sample_sheet(sheet_with_conflicts)
    detection_json = write_detection_json(read_structure(i2=16))
    output = tmp_path / "output.csv"

    run_component(
        [
            "--sample_sheet",
            sample_sheet,
            "--detection_json",
            detection_json,
            "--preset",
            "10x_atac",
            "--output",
            output,
        ]
    )

    settings = parse_settings(output.read_text())
    assert settings["CreateFastqForIndexReads"] == "1"
    assert settings["TrimUMI"] == "0"
    assert settings["OverrideCycles"] == "R1:Y50;I1:I8;I2:U16;R2:Y49"


@pytest.mark.parametrize("preset_args", [["--preset", "10x_atac"], []])
def test_adapter_blanking_is_unconditional(
    run_component, write_sample_sheet, write_detection_json, tmp_path, preset_args
):
    sheet_with_adapters = (
        "[Header]\n"
        "FileFormatVersion,2\n"
        "[BCLConvert_Settings]\n"
        "SoftwareVersion,4.2.7\n"
        "Adapter,AGATCGGAAGAGC\n"
        "AdapterRead1,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA\n"
        "AdapterRead2,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT\n"
        "[BCLConvert_Data]\n"
        "Sample_ID,index,index2\n"
        "sample1,AAAAAAAA,CCCCCCCCCCCCCCCCCCCCCCCC\n"
    )
    sample_sheet = write_sample_sheet(sheet_with_adapters)
    detection_json = write_detection_json(read_structure(i2=16))
    output = tmp_path / "output.csv"

    run_component(
        [
            "--sample_sheet",
            sample_sheet,
            "--detection_json",
            detection_json,
            "--output",
            output,
        ]
        + preset_args
    )

    settings = parse_settings(output.read_text())
    assert settings["Adapter"] == ""
    assert settings["AdapterRead1"] == ""
    assert settings["AdapterRead2"] == ""
    assert settings["SoftwareVersion"] == "4.2.7"
    if not preset_args:
        assert "CreateFastqForIndexReads" not in settings
        assert "TrimUMI" not in settings
        assert "OverrideCycles" not in settings


def test_unknown_sections_preserved_byte_for_byte(
    run_component, write_sample_sheet, write_detection_json, tmp_path
):
    weird_section = [
        "Foo,   Bar space kept   ",
        "# a comment-like line, with a comma",
        ",,,",
        "",
        "Trailing,blank,line,above",
    ]
    sheet_with_unknown_section = (
        BASE_SAMPLE_SHEET + "[MyWeirdSection]\n" + "\n".join(weird_section) + "\n"
    )
    sample_sheet = write_sample_sheet(sheet_with_unknown_section)
    detection_json = write_detection_json(read_structure(i2=16))
    output = tmp_path / "output.csv"

    run_component(
        [
            "--sample_sheet",
            sample_sheet,
            "--detection_json",
            detection_json,
            "--preset",
            "10x_atac",
            "--output",
            output,
        ]
    )

    text = output.read_text()
    assert section_lines(text, "MyWeirdSection") == weird_section


def test_section_order_preserved(
    run_component, write_sample_sheet, write_detection_json, tmp_path
):
    sheet_with_extra_section = BASE_SAMPLE_SHEET + "[MyWeirdSection]\nFoo,Bar\n"
    sample_sheet = write_sample_sheet(sheet_with_extra_section)
    detection_json = write_detection_json(read_structure(i2=16))
    output = tmp_path / "output.csv"

    run_component(
        [
            "--sample_sheet",
            sample_sheet,
            "--detection_json",
            detection_json,
            "--preset",
            "10x_atac",
            "--output",
            output,
        ]
    )

    assert section_names(output.read_text()) == [
        "Header",
        "Reads",
        "BCLConvert_Settings",
        "BCLConvert_Data",
        "MyWeirdSection",
    ]


def test_refuses_existing_per_sample_override_cycles_column(
    run_component, write_sample_sheet, write_detection_json, tmp_path
):
    sheet_with_per_sample_override_cycles = (
        "[Header]\n"
        "FileFormatVersion,2\n"
        "[BCLConvert_Settings]\n"
        "SoftwareVersion,4.2.7\n"
        "[BCLConvert_Data]\n"
        "Sample_ID,index,index2,OverrideCycles\n"
        "sample1,AAAAAAAA,CCCCCCCCCCCCCCCCCCCCCCCC,R1:Y50;I1:I8;I2:U16;R2:Y49\n"
    )
    sample_sheet = write_sample_sheet(sheet_with_per_sample_override_cycles)
    detection_json = write_detection_json(read_structure(i2=16))
    output = tmp_path / "output.csv"

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--sample_sheet",
                sample_sheet,
                "--detection_json",
                detection_json,
                "--preset",
                "10x_atac",
                "--output",
                output,
            ]
        )

    message = err.value.stdout.decode("utf-8")
    assert "OverrideCycles" in message
    assert "per-sample" in message


@pytest.mark.parametrize(
    "structure",
    [
        pytest.param(
            [make_read(1, 50, False), make_read(2, 8, True), make_read(3, 49, False)],
            id="only_one_index_read",
        ),
        pytest.param(
            [
                make_read(1, 50, False),
                make_read(2, 8, True),
                make_read(3, 24, True, is_reverse_complement=True),
                make_read(4, 30, False),
                make_read(5, 19, False),
            ],
            id="three_genomic_reads",
        ),
    ],
)
def test_refuses_read_structure_not_two_genomic_two_indexed(
    run_component, write_sample_sheet, write_detection_json, tmp_path, structure
):
    sample_sheet = write_sample_sheet(BASE_SAMPLE_SHEET)
    detection_json = write_detection_json(structure)
    output = tmp_path / "output.csv"

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--sample_sheet",
                sample_sheet,
                "--detection_json",
                detection_json,
                "--preset",
                "10x_atac",
                "--output",
                output,
            ]
        )

    message = err.value.stdout.decode("utf-8")
    assert "genomic" in message
    assert "indexed" in message


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
