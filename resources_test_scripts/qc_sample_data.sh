#/bin/bash

OUT_DIR=resources_test/qc_sample_data

[ ! -d "$OUT_DIR" ] && mkdir -p "$OUT_DIR"

# fetch/create h5mu from somewhere
cat > /tmp/params.yaml <<EOF
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
  -params-file /tmp/params.yaml \
  -resume

# copy to s3
aws s3 sync \
  --profile di \
  resources_test/qc_sample_data \
  s3://openpipelines-bio/openpipeline_incubator/resources_test/resources_test/qc_sample_data \
  --delete --dryrun
