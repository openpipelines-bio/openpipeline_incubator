#/bin/bash

OUT_DIR=resources_test/spatial_qc_sample_data

[ ! -d "$OUT_DIR" ] && mkdir -p "$OUT_DIR"

# fetch/create h5mu from somewhere
cat > /tmp/qc.yaml <<EOF
param_list:
  - id: xenium_tiny
    input: s3://openpipelines-bio/openpipeline_spatial/resources_test/xenium/xenium_tiny.h5mu
  - id: Lung5_Rep2_tiny
    input: s3://openpipelines-bio/openpipeline_spatial/resources_test/cosmx/Lung5_Rep2_tiny.h5mu
var_name_mitochondrial_genes: mitochondrial
var_name_ribosomal_genes: ribosomal
output: '\$id.qc.h5mu'
output_compression: gzip
publish_dir: "$OUT_DIR"
EOF

nextflow run openpipelines-bio/openpipeline \
  -latest \
  -r 2.1.0 \
  -main-script target/nextflow/workflows/qc/qc/main.nf \
  -profile docker \
  -params-file /tmp/qc.yaml \
  -resume \
  -config src/configs/labels_ci.config

# copy to s3
aws s3 sync \
  --profile di \
  resources_test/spatial_qc_sample_data \
  s3://openpipelines-bio/openpipeline_incubator/resources_test/spatial_qc_sample_data \
  --delete --dryrun \
  --exclude "*" --include "*.h5mu" \
  
