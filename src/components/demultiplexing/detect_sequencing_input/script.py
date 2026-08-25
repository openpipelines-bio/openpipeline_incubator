import json
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

### VIASH START
par = {
    "input": "run_dir",
    "sample_sheet": None,
    "demultiplexer": None,
    "skip_demultiplexing": False,
    "output": "detection.json",
}
meta = {"resources_dir": "src/utils"}
### VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()

# Run-folder markers, checked at the top level of --input only
ELEMENT_MARKERS = ["RunParameters.json", "RunUploaded.json"]
ILLUMINA_MARKER = "RunInfo.xml"

# Illumina FASTQ naming convention, matching create_fastq_manifest
FASTQ_PATTERN = re.compile(
    r"^(?P<sample>[A-Za-z0-9\-_.]+)_S(?P<snum>\d+)(_L(?P<lane>\d+))?_(?P<read>[RI]\d)_(?P<chunk>\d+)\.f(ast)?q\.gz$"
)


def find_top_level_markers(input_dir, names):
    """Names from `names` that are present as files at the top level of `input_dir`."""
    return [name for name in names if (input_dir / name).is_file()]


def find_fastq_files(input_dir):
    """FASTQ files matching the naming convention, searched recursively."""
    return [
        p for p in input_dir.rglob("*") if p.is_file() and FASTQ_PATTERN.match(p.name)
    ]


def parse_run_info(run_info_path):
    """Parse RunInfo.xml for the instrument name and the per-read structure."""
    root = ET.parse(run_info_path).getroot()

    instrument_el = root.find("./Run/Instrument")
    instrument = (
        instrument_el.text.strip()
        if instrument_el is not None and instrument_el.text
        else None
    )

    read_structure = []
    for read_el in root.findall("./Run/Reads/Read"):
        read_structure.append(
            {
                "number": int(read_el.attrib["Number"]),
                "num_cycles": int(read_el.attrib["NumCycles"]),
                "is_indexed": read_el.attrib.get("IsIndexedRead", "N") == "Y",
                "is_reverse_complement": read_el.attrib.get("IsReverseComplement", "N")
                == "Y",
            }
        )
    read_structure.sort(key=lambda read: read["number"])

    return instrument, read_structure


def resolve_sample_sheet(par, input_dir):
    """The resolved sample sheet path (string), or None when not found/given."""
    if par.get("sample_sheet"):
        return str(Path(par["sample_sheet"]).resolve())

    candidate = input_dir / "SampleSheet.csv"
    if candidate.is_file():
        return str(candidate.resolve())

    return None


def detect(par):
    input_dir = Path(par["input"])
    if not input_dir.is_dir():
        raise ValueError(f"--input '{input_dir}' is not a directory.")

    element_markers = find_top_level_markers(input_dir, ELEMENT_MARKERS)
    illumina_markers = find_top_level_markers(input_dir, [ILLUMINA_MARKER])

    if element_markers and illumina_markers:
        raise ValueError(
            f"Conflicting run-folder markers found in '{input_dir}': Illumina marker(s) "
            f"{illumina_markers} and Element AVITI marker(s) {element_markers} are both "
            "present at the top level. This looks like two run folders were merged into "
            "one directory. Provide an input directory for a single platform."
        )

    if element_markers:
        logger.info(
            "Detected Element AVITI run folder in '%s' (marker(s): %s).",
            input_dir,
            element_markers,
        )
        input_type = "element"
        demultiplexer = "bases2fastq"
        instrument = None
        read_structure = []
    elif illumina_markers:
        logger.info(
            "Detected Illumina run folder in '%s' (marker: %s).",
            input_dir,
            illumina_markers[0],
        )
        input_type = "illumina"
        demultiplexer = "bclconvert"
        instrument, read_structure = parse_run_info(input_dir / ILLUMINA_MARKER)
    else:
        fastq_files = find_fastq_files(input_dir)
        if fastq_files:
            logger.info(
                "No run-folder markers found in '%s', detected %d FASTQ file(s). "
                "Treating as already-demultiplexed input.",
                input_dir,
                len(fastq_files),
            )
            input_type = "fastq"
            demultiplexer = None
            instrument = None
            read_structure = []
        else:
            top_level = sorted(p.name for p in input_dir.iterdir())
            raise ValueError(
                f"Could not detect the sequencing input type of '{input_dir}'. Found no "
                "Illumina run-folder marker (RunInfo.xml), no Element AVITI run-folder "
                "marker (RunParameters.json or RunUploaded.json), and no FASTQ files "
                "matching the expected naming convention (searched recursively). "
                f"Top-level contents: {top_level if top_level else '(empty directory)'}."
            )

    requested_demultiplexer = par.get("demultiplexer")
    if requested_demultiplexer is not None and requested_demultiplexer != demultiplexer:
        raise ValueError(
            f"--demultiplexer '{requested_demultiplexer}' does not match the "
            f"demultiplexer detected from '{input_dir}' (detected: {demultiplexer!r}, "
            f"input_type: '{input_type}')."
        )

    if par.get("skip_demultiplexing"):
        logger.info(
            "--skip_demultiplexing set: forcing input_type='fastq', demultiplexer=null."
        )
        input_type = "fastq"
        demultiplexer = None
        instrument = None
        read_structure = []

    return {
        "input_type": input_type,
        "demultiplexer": demultiplexer,
        "sample_sheet": resolve_sample_sheet(par, input_dir),
        "instrument": instrument,
        "read_structure": read_structure,
    }


result = detect(par)

logger.info("Writing detection result to %s", par["output"])
with open(par["output"], "w") as f:
    json.dump(result, f, indent=2)
