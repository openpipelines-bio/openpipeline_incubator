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
  -profile docker,no_publish \
  -c src/utils/labels_ci.config \
  -c src/utils/integration_tests.config
