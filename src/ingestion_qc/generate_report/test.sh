#!/bin/bash

viash ns build --setup cb --parallel

nextflow run . \
  -main-script target/nextflow/ingestion_qc/generate_report/main.nf \
  -params-file src/ingestion_qc/generate_report/example.yaml \
  -profile docker \
  --publish_dir test_results \
  --resume 