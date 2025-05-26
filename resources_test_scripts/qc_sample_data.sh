#/bin/bash

OUT_DIR=resources_test/qc_sample_data
OUT_DIR_SPATIAL=resources_test/spatial_qc_sample_data

[ ! -d "$OUT_DIR" ] && mkdir -p "$OUT_DIR"
[ ! -d "$OUT_DIR_SPATIAL" ] && mkdir -p "$OUT_DIR_SPATIAL"

# fetch/create h5mu from somewhere
cat > /tmp/params_create_h5mu.yaml <<EOF
param_list:
  - id: sample_one
    input_id: sample_one
    input: s3://openpipelines-data/10x_5k_anticmv/5k_human_antiCMV_T_TBNK_connect_qc.h5mu
  - id: sample_two
    input_id: sample_two
    input: s3://openpipelines-data/10x_5k_anticmv/5k_human_antiCMV_T_TBNK_connect_qc.h5mu
output: '\$id.qc.h5mu'
output_compression: gzip
publish_dir: "$OUT_DIR"
EOF

# add the sample ID to the mudata object
nextflow run openpipelines-bio/openpipeline \
  -latest \
  -r 2.0.0 \
  -main-script target/nextflow/metadata/add_id/main.nf \
  -profile docker \
  -params-file /tmp/params_create_h5mu.yaml \
  -resume

cat > /tmp/params_subset.yaml <<EOF
param_list:
  - id: sample_one
    input: resources_test/qc_sample_data/sample_one.qc.h5mu
  - id: sample_two
    input: resources_test/qc_sample_data/sample_two.qc.h5mu
output: '\$id.qc.h5mu'
number_of_observations: 10000
output_compression: gzip
publish_dir: "$OUT_DIR"
EOF

# subset h5mus
nextflow run openpipelines-bio/openpipeline \
  -latest \
  -r 2.0.0 \
  -main-script target/nextflow/filter/subset_h5mu/main.nf \
  -profile docker \
  -params-file /tmp/params_subset.yaml \
  -resume

cat > /tmp/add_metadata_obs.py <<EOF
import mudata as mu
import glob
import numpy as np
import pandas as pd
import os

# Directory containing the h5mu files
out_dir = "$(pwd)/resources_test/qc_sample_data"

# List of h5mu files
h5mu_files = glob.glob(os.path.join(out_dir, "*.h5mu"))
print(f"Found {len(h5mu_files)} h5mu files: {h5mu_files}")

# Metadata values to randomly assign
donor_ids = ["donor_1", "donor_2", "donor_3"]
cell_types = ["CD4+ T cell", "CD8+ T cell", "B cell", "NK cell", "Monocyte"]
batches = ["batch_A", "batch_B"]
conditions = ["treated", "control"]

for h5mu_file in h5mu_files:
    print(f"Processing {h5mu_file}...")
    
    # Load MuData object
    mdata = mu.read_h5mu(h5mu_file)
    rna = mdata.mod["rna"]
    n_obs = rna.n_obs
    
    # Generate random metadata
    np.random.seed(42 + hash(h5mu_file) % 100)  # Different seed for each file but reproducible
    
    # Create metadata
    rna.obs["donor_id"] = np.random.choice(donor_ids, size=n_obs)
    rna.obs["cell_type"] = np.random.choice(cell_types, size=n_obs)
    rna.obs["batch"] = np.random.choice(batches, size=n_obs)
    rna.obs["condition"] = np.random.choice(conditions, size=n_obs)
    
    # Add a continuous variable too
    rna.obs["quality_score"] = np.random.uniform(0, 1, size=n_obs)
    
    # Save the modified MuData object
    mu.write_h5mu(h5mu_file, mdata)
    print(f"Added metadata to {h5mu_file}")

print("All files processed successfully!")
EOF

# Execute the Python script
python /tmp/add_metadata_obs.py

# generate cellbender out for testing
cat > /tmp/params_cellbender.yaml <<EOF
param_list:
  - id: sample_one
    input: resources_test/qc_sample_data/sample_one.qc.h5mu
  - id: sample_two
    input: resources_test/qc_sample_data/sample_two.qc.h5mu
output: '\$id.qc.cellbender.h5mu'
epochs: 5
output_compression: gzip
publish_dir: "$OUT_DIR"
EOF

nextflow run openpipelines-bio/openpipeline \
  -latest \
  -r 2.0.0 \
  -main-script target/nextflow/correction/cellbender_remove_background/main.nf \
  -profile docker \
  -params-file /tmp/params_cellbender.yaml \
  -resume

# fetch spatial sample data from s3
aws s3 sync \
  --profile di \
  s3://openpipelines-bio/openpipeline_incubator/resources_test/spatial_qc_sample_data \
  "$OUT_DIR_SPATIAL"

# generate json for testing
viash run src/ingestion_qc/h5mu_to_qc_json/config.vsh.yaml --engine docker -- \
  --input "$OUT_DIR"/sample_one.qc.cellbender.h5mu \
  --input "$OUT_DIR"/sample_two.qc.cellbender.h5mu \
  --ingestion_method cellranger_multi \
  --obs_metadata "donor_id;cell_type;batch;condition" \
  --output "$OUT_DIR"/sc_dataset.json \
  --output_reporting_json "$OUT_DIR"/sc_report_structure.json

viash run src/ingestion_qc/h5mu_to_qc_json/config.vsh.yaml --engine docker -- \
  --input "$OUT_DIR_SPATIAL"/xenium_tiny.qc.h5mu \
  --input "$OUT_DIR_SPATIAL"/xenium_tiny.qc.h5mu \
  --ingestion_method xenium \
  --min_num_nonzero_vars 1 \
  --output "$OUT_DIR_SPATIAL"/xenium_dataset.json \
  --output_reporting_json "$OUT_DIR_SPATIAL"/xenium_report_structure.json


# copy to s3
aws s3 sync \
  --profile di \
  resources_test/qc_sample_data \
  s3://openpipelines-bio/openpipeline_incubator/resources_test/qc_sample_data \
  --delete --dryrun \
  --exclude "*.yaml" --include "*.h5mu"  --include "*.json"


aws s3 sync \
  --profile di \
  resources_test/spatial_qc_sample_data \
  s3://openpipelines-bio/openpipeline_incubator/resources_test/spatial_qc_sample_data \
  --delete --dryrun \
  --exclude "*.yaml" --include "*.h5mu"  --include "*.json"
