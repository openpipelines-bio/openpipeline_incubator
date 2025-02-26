#!/bin/bash

viash ns build --setup cb --parallel

cat > /tmp/params.yaml <<EOF
param_list:
  - input: resources_test/sample_data/sample_1.qc.output.h5mu
    id: sample_one
  - input: resources_test/sample_data/sample_2.qc.output.h5mu
    id: sample_two
output_qc_json: output_qc.json
output_html: output_report.html
EOF

nextflow run . \
  -main-script target/nextflow/ingestion_qc/generate_report/main.nf \
  -params-file /tmp/params.yaml \
  -profile docker \
  --publish_dir test_results \
  --resume
