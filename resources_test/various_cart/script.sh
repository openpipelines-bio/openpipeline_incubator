#!/bin/bash

cat > /tmp/params.yaml << HERE
param_list:
  - id: various_cart_sample_7
    input: resources_test/various_cart/various_cart_sample_7.from_cellranger_multi_to_h5mu.output_0.h5mu
  - id: various_cart_sample_11
    input: resources_test/various_cart/various_cart_sample_11.from_cellranger_multi_to_h5mu.output_0.h5mu
  - id: various_cart_sample_15
    input: resources_test/various_cart/various_cart_sample_15.from_cellranger_multi_to_h5mu.output_0.h5mu
publish_dir: "data/sample_data/various_cart"
HERE

# nextflow run openpipelines-bio/incubator -latest \
nextflow run . \
  -main-script target/nextflow/ingestion_qc/generate_report/main.nf \
  -resume \
  -profile docker,mount_temp \
  -params-file /tmp/params.yaml \
  -c src/labels.config
