#!/bin/bash

set -eo pipefail

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

nextflow \
  run . \
  -main-script src/atlas_service/test.nf \
  -entry test_wf \
  -resume \
  -profile docker \
  -c src/configs/labels_ci.config \
  -c src/configs/integration_tests.config \
  --publish_dir test

nextflow \
  run . \
  -main-script src/atlas_service/test.nf \
  -profile docker,no_publish \
  -resume \
  -entry test_wf_2 \
  -c src/configs/labels_ci.config \
  -c src/configs/integration_tests.config

  nextflow \
  run . \
  -main-script src/atlas_service/test.nf \
  -profile docker,no_publish \
  -resume \
  -entry test_wf_3 \
  -c src/configs/labels_ci.config \
  -c src/configs/integration_tests.config

  nextflow \
  run . \
  -main-script src/atlas_service/test.nf \
  -profile docker,no_publish \
  -resume \
  -entry test_wf_4 \
  -c src/configs/labels_ci.config \
  -c src/configs/integration_tests.config