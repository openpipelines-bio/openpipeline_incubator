import pytest
import mudata as mu
import os
import json
import sys

## VIASH START
meta = {
    "resources_dir": "resources_test",
    "executable": "./target/executable/ingestion_qc/h5mu_to_qc_json/h5mu_to_qc_json"
}
## VIASH END

def test_simple_execution(run_component, tmp_path):
    output_json_path = tmp_path / "output.json"
    
    run_component(
        [
            "--input", meta["resources_dir"] + "/resources_test/qc_sample_data/sample_one.qc.h5mu",
            "--input", meta["resources_dir"] + "/resources_test/qc_sample_data/sample_two.qc.h5mu",
            "--output", output_json_path,
        ]
    )
    
    assert os.path.exists(output_json_path), "Output file was not created"
    
    with open(output_json_path, "r") as f:
        output_json_dict = json.load(f)
    
    assert output_json_dict.keys() == {"cell_rna_stats", "sample_summary_stats", "metrics_cellranger_stats"}
    
    column_names = [col["name"] for col in output_json_dict["cell_rna_stats"]["columns"]]
    assert column_names == ["sample_id", "total_counts", "num_nonzero_vars", "fraction_mitochondrial", "fraction_ribosomal"]
        
    for key in output_json_dict.keys():
        assert output_json_dict[key].keys() == {"num_rows", "num_cols", "columns"}
        for col in output_json_dict[key]["columns"]:
            assert {"name", "dtype", "data"}.issubset(col.keys())
        

def test_set_filters(run_component, tmp_path):
    output_json_path = tmp_path / "output.json"
    
    run_component(
        [
            "--input", meta["resources_dir"] + "/resources_test/qc_sample_data/sample_one.qc.h5mu",
            "--input", meta["resources_dir"] + "/resources_test/qc_sample_data/sample_two.qc.h5mu",
            "--output", output_json_path,
            "--sample_id_key", "sample_id",
            "--min_total_counts", "10",
            "--min_num_nonzero_vars", "10",
            "--obs_keys", "total_counts", 
            "--obs_keys", "num_nonzero_vars"
        ]
    )
        
    assert os.path.exists(output_json_path), "Output file was not created"
    
    with open(output_json_path, "r") as f:
        output_json_dict = json.load(f)
    
    assert output_json_dict.keys() == {"cell_rna_stats", "sample_summary_stats", "metrics_cellranger_stats"}
    
    column_names = [col["name"] for col in output_json_dict["cell_rna_stats"]["columns"]]
    assert column_names == ["sample_id", "total_counts", "num_nonzero_vars"]
        
    for key in output_json_dict.keys():
        assert output_json_dict[key].keys() == {"num_rows", "num_cols", "columns"}
        for col in output_json_dict[key]["columns"]:
            assert {"name", "dtype", "data"}.issubset(col.keys())
            
    total_counts = next(col for col in output_json_dict["cell_rna_stats"]["columns"] if col["name"] == "total_counts")
    assert min(total_counts["data"]) >= 10
    
    num_nonzero_vars = next(col for col in output_json_dict["cell_rna_stats"]["columns"] if col["name"] == "num_nonzero_vars")
    assert min(num_nonzero_vars["data"]) >= 10
    
if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))