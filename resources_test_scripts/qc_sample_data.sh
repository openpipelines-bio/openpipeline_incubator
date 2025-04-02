#/bin/bash

OUT_DIR=resources_test/qc_sample_data

[ ! -d "$OUT_DIR" ] && mkdir -p "$OUT_DIR"

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

# generate json for testing
viash run src/ingestion_qc/h5mu_to_qc_json/config.vsh.yaml --engine docker -- \
  --input "$OUT_DIR"/sample_one.qc.cellbender.h5mu \
  --input "$OUT_DIR"/sample_two.qc.cellbender.h5mu \
  --output "$OUT_DIR"/dataset.json

# copy to s3
aws s3 sync \
  --profile di \
  resources_test/qc_sample_data \
  s3://openpipelines-bio/openpipeline_incubator/resources_test/qc_sample_data \
  --delete --dryrun
