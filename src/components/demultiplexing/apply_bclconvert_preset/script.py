import csv
import io
import json
import re
import sys

### VIASH START
par = {
    "sample_sheet": "SampleSheet.csv",
    "detection_json": "detection.json",
    "preset": None,
    "output": "SampleSheet.patched.csv",
}
meta = {"resources_dir": "src/utils"}
### VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()

# Sample sheets are a sequence of '[SectionName]' headers each followed by CSV-style lines
SECTION_HEADER_RE = re.compile(r"^\[(.+)\]$")

# The BCLConvert_Settings keys this component touches
ADAPTER_KEYS = ("Adapter", "AdapterRead1", "AdapterRead2")
BCLCONVERT_SETTINGS = "BCLConvert_Settings"
BCLCONVERT_DATA = "BCLConvert_Data"


def parse_sections(text):
    """Parse a sample sheet into an ordered list of sections.

    Each section is a dict with 'name', 'header' and 'lines'.
    """
    sections = []
    current = {"name": None, "header": None, "lines": []}
    for line in text.splitlines():
        match = SECTION_HEADER_RE.match(line.strip())
        if match:
            if current["header"] is not None or current["lines"]:
                sections.append(current)
            current = {"name": match.group(1), "header": line, "lines": []}
        else:
            current["lines"].append(line)
    if current["header"] is not None or current["lines"]:
        sections.append(current)
    return sections


def serialize_sections(sections):
    lines = []
    for section in sections:
        if section["header"] is not None:
            lines.append(section["header"])
        lines.extend(section["lines"])

    return "\n".join(lines) + "\n"


def find_section(sections, name):
    for section in sections:
        if section["name"] == name:
            return section

    return None


def get_or_create_bclconvert_settings(sections):
    """Return the [BCLConvert_Settings] section, creating and inserting one if absent."""
    section = find_section(sections, BCLCONVERT_SETTINGS)
    if section is not None:
        return section

    logger.info("No [%s] section found; creating one.", BCLCONVERT_SETTINGS)
    new_section = {
        "name": BCLCONVERT_SETTINGS,
        "header": f"[{BCLCONVERT_SETTINGS}]",
        "lines": [],
    }

    # Insert before the first "_Data" section
    insert_at = len(sections)
    for index, existing in enumerate(sections):
        if existing["name"] is not None and existing["name"].endswith("_Data"):
            insert_at = index
            break
    sections.insert(insert_at, new_section)

    return new_section


def parse_key_value_lines(lines):
    """Parse 'Key,Value' lines into an order-preserving dict, skipping blank lines."""
    settings = {}
    for line in lines:
        if not line.strip():
            continue
        row = next(csv.reader([line]))
        if not row:
            continue
        key = row[0]
        value = row[1] if len(row) > 1 else ""
        settings[key] = value

    return settings


def serialize_key_value(settings):
    lines = []
    for key, value in settings.items():
        buf = io.StringIO()
        csv.writer(buf, lineterminator="").writerow([key, value])
        lines.append(buf.getvalue())

    return lines


def parse_header_columns(lines):
    """The column names of the first non-blank line of a '*_Data' section."""
    for line in lines:
        if not line.strip():
            continue
        return next(csv.reader([line]))

    return []


def compute_override_cycles(read_structure):
    """Compute the OverrideCycles string for the 10x_atac preset.

    Every genomic (non-indexed) read passes through unchanged as 'Y<cycles>'. Of
    the two indexed reads, the first stays a demultiplexing index 'I<cycles>',
    the second (index2) is reclassified from index to UMI ('U<cycles>') because
    the 10x ATAC i5 is a cell barcode with too many possible sequences to be
    used for demultiplexing. The two published 10x examples differ only in this
    UMI length (U16 for single-cell ATAC, U24 for Multiome ATAC).

    Refuses (raises ValueError) unless there are exactly two genomic and two
    indexed reads, since the ATAC assumption this preset encodes does not hold
    for anything else.
    """
    ordered = sorted(read_structure, key=lambda read: read["number"])
    genomic = [read for read in ordered if not read["is_indexed"]]
    indexed = [read for read in ordered if read["is_indexed"]]

    if len(genomic) != 2 or len(indexed) != 2:
        raise ValueError(
            "The 10x_atac preset requires exactly two genomic (non-indexed) reads and "
            "two indexed reads in the detection JSON's read_structure, but found "
            f"{len(genomic)} genomic and {len(indexed)} indexed read(s)."
        )

    r1, r2 = genomic
    i1, i2 = indexed
    return (
        f"R1:Y{r1['num_cycles']};I1:I{i1['num_cycles']};I2:U{i2['num_cycles']};"
        f"R2:Y{r2['num_cycles']}"
    )


def set_setting(settings, key, value):
    """Set settings[key] = value, warning if this overwrites."""
    if key in settings and settings[key] != value:
        logger.warning(
            "Overwriting existing '%s' value '%s' with '%s'.",
            key,
            settings[key],
            value,
        )

    settings[key] = value


def blank_adapters(settings):
    """Blank any Adapter/AdapterRead1/AdapterRead2 entries, keeping the keys present."""
    for key in ADAPTER_KEYS:
        if settings.get(key):
            logger.info("Blanking existing '%s' value '%s'.", key, settings[key])
            settings[key] = ""


def apply_10x_atac_preset(settings, sections, read_structure):
    # BCL Convert (from 4.1) forbids a setting specified both globally and
    # per-sample. A mixed GEX+ATAC (Multiome) flowcell commonly has a per-sample
    # OverrideCycles column, so writing a global one here would silently break
    # it.
    data_section = find_section(sections, BCLCONVERT_DATA)
    if data_section is not None:
        columns = parse_header_columns(data_section["lines"])
        if "OverrideCycles" in columns:
            raise ValueError(
                f"The sample sheet's [{BCLCONVERT_DATA}] section already has a per-sample "
                "'OverrideCycles' column. BCL Convert forbids specifying a setting both "
                "globally and per-sample, so the 10x_atac preset refuses to also write a "
                "global OverrideCycles setting."
            )

    override_cycles = compute_override_cycles(read_structure)
    set_setting(settings, "CreateFastqForIndexReads", "1")
    set_setting(settings, "TrimUMI", "0")
    set_setting(settings, "OverrideCycles", override_cycles)


logger.info("Reading detection JSON from %s", par["detection_json"])
with open(par["detection_json"]) as f:
    detection = json.load(f)

read_structure = detection.get("read_structure") or []

logger.info("Reading sample sheet from %s", par["sample_sheet"])
with open(par["sample_sheet"]) as f:
    sections = parse_sections(f.read())

bclconvert_settings_section = get_or_create_bclconvert_settings(sections)
settings = parse_key_value_lines(bclconvert_settings_section["lines"])

# Always blank adapter settings
blank_adapters(settings)

if par["preset"] == "10x_atac":
    logger.info("Applying preset '10x_atac'")
    apply_10x_atac_preset(settings, sections, read_structure)
elif par["preset"]:
    raise ValueError(f"Unknown preset '{par['preset']}'.")

bclconvert_settings_section["lines"] = serialize_key_value(settings)

logger.info("Writing patched sample sheet to %s", par["output"])
with open(par["output"], "w") as f:
    f.write(serialize_sections(sections))

logger.info("Finished")
