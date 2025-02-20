#!/bin/bash

viash ns build --setup cb --parallel -q generate_report

nextflow run . \
  -main-script src/ingestion_qc/generate_report/main.nf \
  -params-file src/ingestion_qc/generate_report/example.yaml \
  -entry run_wf \
  --publish_dir test_results \
  --resume 