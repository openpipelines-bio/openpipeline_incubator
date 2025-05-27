#!/usr/bin/env bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

viash ns build --setup cb -q generate_report

nextflow run . \
  -main-script src/ingestion_qc/generate_report/test.nf \
  -profile docker,no_publish,local \
  -entry test_no_cellbender \
  -c src/configs/labels_ci.config \
  -resume

nextflow run . \
  -main-script src/ingestion_qc/generate_report/test.nf \
  -profile docker,no_publish,local \
  -entry test_with_cellbender \
  -c src/configs/labels_ci.config \
  -resume

nextflow run . \
  -main-script src/ingestion_qc/generate_report/test.nf \
  -profile docker,no_publish,local \
  -entry test_multiple_reports \
  -c src/configs/labels_ci.config \
  -resume
