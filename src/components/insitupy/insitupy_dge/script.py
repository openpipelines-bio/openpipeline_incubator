from insitupy import InSituData, InSituExperiment
from insitupy.tools import dge
from pathlib import Path
import sys

## VIASH START
par = {
    "input": "resources_test/insitupy/insitupy_project",  # Not existing yet
    "cells_layer": None,
    "target_id": 0,
    "ref_id": None,
    "target_annotation_tuple": ["pathologist", "tumor"],
    "target_cell_type_tuple": ["cell_type", "T cell"],
    "target_region_tuple": ["tissue_regions", "region_1"],
    "target_name": "tumor",
    "target_metadata": [ "condition=treated", "patient=P1" ],
    "ref": None,
    "ref_annotation_tuple": [ "pathologist", "healthy" ],
    "ref_cell_type_tuple": ["same"],
    "ref_region_tuple": ["same"],
    "ref_name": "healthy",
    "ref_metadata": [ "condition=control", "patient=P2" ],
    "method": "t-test",
    "consider_neighbors": False,
    "exclude_ambiguous_assignments": False,
    "force_assignment": False,
    "skip_log1p_check": False,
    "output": "dge_results",
}
meta = {"temp_dir": "/tmp", "resources_dir": "."}
## VIASH END
sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()
input_file = Path(par["input"])
output_file = Path(par["output"])

def to_tuple(values, arg, keywords=()):
    if not values:
        return None
    if len(values) == 1 and values[0] in keywords:
        return values[0]
    if len(values) != 2:
        allowed = f" or one of {list(keywords)}" if keywords else ""
        raise ValueError(f"--{arg} needs exactly two values (key, name){allowed}, got {values}.")
    return tuple(values)

def to_dict(input_pairs: list):
    dict_pairs = {}
    for pair in input_pairs: 
        pair_splitted = pair.split("=")
        dict_pairs[pair_splitted[0]] = pair_splitted[1]
    return dict_pairs

# Arguments accepted by both insitupy.tools.dge and InSituExperiment.dge
kwargs = {
    "target_annotation_tuple": to_tuple(par["target_annotation_tuple"], "target_annotation_tuple"),
    "target_cell_type_tuple": to_tuple(par["target_cell_type_tuple"], "target_cell_type_tuple"),
    "target_region_tuple": to_tuple(par["target_region_tuple"], "target_region_tuple"),
    "ref_annotation_tuple": to_tuple(par["ref_annotation_tuple"], "ref_annotation_tuple", ("same", "rest")),
    "ref_cell_type_tuple": to_tuple(par["ref_cell_type_tuple"], "ref_cell_type_tuple", ("same", "rest")),
    "ref_region_tuple": to_tuple(par["ref_region_tuple"], "ref_region_tuple", ("same",)),
    "method": par["method"],
    "exclude_ambiguous_assignments": par["exclude_ambiguous_assignments"],
    "force_assignment": par["force_assignment"],
    "assert_log1p": not par["skip_log1p_check"],
}

# Arguments only accepted by insitupy.tools.dge (single-sample input)
insitudata_kwargs = {
    "cells_layer": par["cells_layer"],
    "target_name": par["target_name"],
    "ref_name": par["ref_name"],
    "target_metadata": to_dict(par["target_metadata"]) if par["target_metadata"] else None,
    "ref_metadata": to_dict(par["ref_metadata"]) if par["ref_metadata"] else None,
    "consider_neighbors": par["consider_neighbors"],
}

def convert_ref_id(ref_input: list[str] | None):
    if not ref_input:
        return None
    if len(ref_input) == 1:
        if ref_input[0] == "rest":
            return "rest"
        return int(ref_input[0])
    return [int(x) for x in ref_input]

def dge_insitudata(xd_obj: InSituData, ref: list[InSituData] | None = None):
    dge_results = dge(
        target = xd_obj, 
        ref = ref,
        **kwargs,
        **insitudata_kwargs,
    )
    return dge_results

def dge_experiment(exp: InSituExperiment, target_id: int, ref_id: list[int] | str | None ):
    dge_results = exp.dge(
        target_id = target_id, 
        ref_id = convert_ref_id(ref_id),
        **kwargs,
    )
    return dge_results

def main():
    logger.info("Loading InSituPy object...")
    if (input_file / "metadata.parquet").exists():
        exp = InSituExperiment.read(input_file)
        exp.load_all()
        dge_results = dge_experiment(exp, par["target_id"], par["ref_id"])

    else:
        xd = InSituData.read(input_file)
        xd.load_all()
        ref = None
        if par["ref"]:
            ref = [InSituData.read(path) for path in par["ref"]]
            for ref_xd in ref:
                ref_xd.load_all()
        dge_results = dge_insitudata(xd, ref)

    logger.info("Saving dge results to '%s'...", output_file)
    dge_results.save(output_file)
    logger.info("Done!")


if __name__ == "__main__":
    main()

