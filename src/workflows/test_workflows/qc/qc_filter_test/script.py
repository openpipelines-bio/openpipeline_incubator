import sys
import pytest
import pandas as pd
from mudata import read_h5mu

## VIASH START
par = {
    "input": "output.h5mu",
    "og_input": "input.h5mu",
    "modality": "rna",
    "keep_flag_columns": [
        "filter_counts_rna",
        "filter_quantile_rna",
        "filter_mito_rna",
    ],
}
## VIASH END


def test_run():
    og_mudata = read_h5mu(par["og_input"])
    output_mudata = read_h5mu(par["input"])

    modality = par["modality"]

    # The workflow operates in place on a single modality and emits the whole object: every
    # modality present in the input must still be present in the output.
    assert set(og_mudata.mod.keys()) == set(output_mudata.mod.keys()), (
        f"Modalities should be unchanged. Input: {list(og_mudata.mod.keys())}, "
        f"output: {list(output_mudata.mod.keys())}"
    )
    assert modality in output_mudata.mod, (
        f"Tested modality '{modality}' should be present in the output."
    )

    og_mod = og_mudata.mod[modality]
    output_mod = output_mudata.mod[modality]

    # No filtering: the filters run in flag-only mode, so no cells may be dropped.
    assert output_mod.n_obs == og_mod.n_obs, (
        f"No cells should be filtered out. Input: {og_mod.n_obs} cells, "
        f"output: {output_mod.n_obs} cells."
    )
    pd.testing.assert_index_equal(
        og_mod.obs_names, output_mod.obs_names, check_order=False
    )

    # No gene filtering either: the sub-workflows never subset .var.
    assert output_mod.n_vars == og_mod.n_vars, (
        f"No genes should be filtered out. Input: {og_mod.n_vars} genes, "
        f"output: {output_mod.n_vars} genes."
    )
    pd.testing.assert_index_equal(
        og_mod.var_names, output_mod.var_names, check_order=False
    )

    # The expected boolean keep-flag columns must have been written to .obs.
    for column in par["keep_flag_columns"]:
        assert column in output_mod.obs.columns, (
            f"Expected keep-flag column '{column}' to be present in .obs of modality "
            f"'{modality}'. Found: {output_mod.obs.columns.to_list()}"
        )
        # Accepts both numpy bool (count/quantile filters) and pandas nullable boolean
        # (scrublet writes the latter).
        assert pd.api.types.is_bool_dtype(output_mod.obs[column]), (
            f"Keep-flag column '{column}' should be boolean. "
            f"Found dtype: {output_mod.obs[column].dtype}"
        )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "--import-mode=importlib"]))
