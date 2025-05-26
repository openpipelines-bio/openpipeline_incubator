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


def test_cellranger_execution(run_component, tmp_path):
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
            "--obs_sample_id", "sample_id",
            "--obs_total_counts", "total_counts",
            "--obs_num_nonzero_vars", "num_nonzero_vars",
            "--obs_fraction_mitochondrial", "fraction_mitochondrial",
            "--obs_fraction_ribosomal", "fraction_ribosomal",
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


def test_xenium_execution(run_component, tmp_path):
    output_json_path = tmp_path / "output.json"
    output_reporting_json_path = tmp_path / "output_reporting.json"

    run_component(
        [
            "--input", meta["resources_dir"] + "/resources_test/spatial_qc_sample_data/xenium_tiny.qc.h5mu",
            "--input", meta["resources_dir"] + "/resources_test/spatial_qc_sample_data/xenium_tiny.qc.h5mu",
            "--ingestion_method", "xenium",
            "--min_num_nonzero_vars", "1",
            "--output", output_json_path,
            "--output_reporting_json", output_reporting_json_path
        ]
    )

    assert os.path.exists(output_json_path), "Output file was not created"

    with open(output_json_path, "r") as f:
        output_json_dict = json.load(f)

    assert output_json_dict.keys() == {"cell_rna_stats", "sample_summary_stats"}
    assert "metrics_cellranger_stats" not in output_json_dict.keys()

    column_names_cell = [col["name"] for col in output_json_dict["cell_rna_stats"]["columns"]]
    expected_column_names = [
        "sample_id", "total_counts", "num_nonzero_vars",
        "fraction_mitochondrial",  "fraction_ribosomal",
        "cell_area", "nucleus_ratio",
        "x_coord", "y_coord"
        ]
    assert np.all([column in column_names_cell for column in expected_column_names])

    for key in output_json_dict.keys():
        assert output_json_dict[key].keys() == {"num_rows", "num_cols", "columns"}
        for col in output_json_dict[key]["columns"]:
            assert {"name", "dtype", "data"}.issubset(col.keys())


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
