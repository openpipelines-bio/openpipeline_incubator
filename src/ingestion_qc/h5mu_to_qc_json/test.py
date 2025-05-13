import pytest
import os
import json
import sys
import numpy as np

## VIASH START
meta = {
    "resources_dir": "resources_test",
    "executable": "./target/executable/ingestion_qc/h5mu_to_qc_json/h5mu_to_qc_json"
}
## VIASH END


def test_simple_execution(run_component, tmp_path):
    output_json_path = tmp_path / "output.json"
    output_reporting_json_path = tmp_path / "output_reporting.json"

    run_component(
        [
            "--input", meta["resources_dir"] + "/resources_test/qc_sample_data/sample_one.qc.cellbender.h5mu",
            "--input", meta["resources_dir"] + "/resources_test/qc_sample_data/sample_two.qc.cellbender.h5mu",
            "--ingestion_method", "cellranger_multi",
            "--output", output_json_path,
            "--output_reporting_json", output_reporting_json_path
        ]
    )

    assert os.path.exists(output_json_path), "Output file was not created"

    with open(output_json_path, "r") as f:
        output_json_dict = json.load(f)

    assert output_json_dict.keys() == {"cell_rna_stats", "sample_summary_stats", "metrics_cellranger_stats"}

    column_names_cell = [col["name"] for col in output_json_dict["cell_rna_stats"]["columns"]]
    expected_column_names = [
        "sample_id", "total_counts", "num_nonzero_vars",
        "fraction_mitochondrial",  "fraction_ribosomal",
        "cellbender_background_fraction", "cellbender_cell_probability",
        "cellbender_cell_size", "cellbender_droplet_efficiency"
        ]
    assert np.all([column in column_names_cell for column in expected_column_names])

    for key in output_json_dict.keys():
        assert output_json_dict[key].keys() == {"num_rows", "num_cols", "columns"}
        for col in output_json_dict[key]["columns"]:
            assert {"name", "dtype", "data"}.issubset(col.keys())


def test_set_filters(run_component, tmp_path):
    output_json_path = tmp_path / "output.json"
    output_reporting_json_path = tmp_path / "output_reporting.json"

    run_component(
        [
            "--input", meta["resources_dir"] + "/resources_test/qc_sample_data/sample_one.qc.cellbender.h5mu",
            "--input", meta["resources_dir"] + "/resources_test/qc_sample_data/sample_two.qc.cellbender.h5mu",
            "--ingestion_method", "cellranger_multi",
            "--output", output_json_path,
            "--output_reporting_json", output_reporting_json_path,
            "--input_obs_sample_id_key", "sample_id",
            "--input_obs_total_counts_key", "total_counts",
            "--input_obs_num_nonzero_vars_key", "num_nonzero_vars",
            "--input_obs_fraction_mitochondrial_key", "fraction_mitochondrial",
            "--input_obs_fraction_ribosomal_key", "fraction_ribosomal",
            "--min_total_counts", "20",
            "--min_num_nonzero_vars", "20"
        ]
    )

    assert os.path.exists(output_json_path), "Output file was not created"

    with open(output_json_path, "r") as f:
        output_json_dict = json.load(f)

    assert output_json_dict.keys() == {"cell_rna_stats", "sample_summary_stats", "metrics_cellranger_stats"}

    column_names = [col["name"] for col in output_json_dict["cell_rna_stats"]["columns"]]
    expected_column_names = [
        "sample_id", "total_counts", "num_nonzero_vars",
        "fraction_mitochondrial",  "fraction_ribosomal",
        "cellbender_background_fraction", "cellbender_cell_probability",
        "cellbender_cell_size", "cellbender_droplet_efficiency"
        ]
    assert np.all([column in column_names for column in expected_column_names])
    for key in output_json_dict.keys():
        assert output_json_dict[key].keys() == {"num_rows", "num_cols", "columns"}
        for col in output_json_dict[key]["columns"]:
            assert {"name", "dtype", "data"}.issubset(col.keys())

    total_counts = next(col for col in output_json_dict["cell_rna_stats"]["columns"] if col["name"] == "total_counts")
    assert min(total_counts["data"]) >= 20

    num_nonzero_vars = next(col for col in output_json_dict["cell_rna_stats"]["columns"] if col["name"] == "num_nonzero_vars")
    assert min(num_nonzero_vars["data"]) >= 20


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
