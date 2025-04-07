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
  -c src/config/labels.config \
  --resources_test s3://openpipelines-bio/openpipeline_incubator/resources_test/qc_sample_data/ \
  -resume

nextflow run . \
  -main-script src/ingestion_qc/generate_report/test.nf \
  -profile docker,no_publish,local \
  -entry test_with_cellbender \
  -c src/config/labels.config \
  --resources_test s3://openpipelines-bio/openpipeline_incubator/resources_test/qc_sample_data/ \
  -resume