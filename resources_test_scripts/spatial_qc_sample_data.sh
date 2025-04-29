#/bin/bash

OUT_DIR=resources_test/spatial_qc_sample_data

[ ! -d "$OUT_DIR" ] && mkdir -p "$OUT_DIR"

# Fetch Xenium H5MU data
cat > /tmp/xenium_tiny.yaml <<EOF
id: xenium_tiny
input: s3://openpipelines-bio/openpipeline_spatial/resources_test/xenium/xenium_tiny.h5mu
output: xenium_tiny_qc.h5mu
var_name_mitochondrial_genes: mitochondrial_genes
var_name_ribosomal_genes: ribosomal_genes
output_compression: gzip
publish_dir: resources_test/spatial_qc_sample_data/xenium/
EOF

# Run QC workflow on xenium data
nextflow run openpipelines-bio/openpipeline \
  -latest \
  -r 2.1.0 \
  -main-script target/nextflow/workflows/qc/qc/main.nf \
  -profile docker \
  -params-file /tmp/xenium_tiny.yaml \
  -resume \
  -config src/configs/labels_ci.config

# generate json for testing
viash run src/spatial_ingestion_qc/h5mu_to_qc_json/config.vsh.yaml --engine docker -- \
  --input "$OUT_DIR"/xenium/xenium_tiny_qc.h5mu \
  --min_total_counts "1" \
  --min_num_nonzero_vars "1" \
  --input "$OUT_DIR"/xenium/xenium_tiny_qc.h5mu \
  --output "$OUT_DIR"/xenium/xenium_tiny.json