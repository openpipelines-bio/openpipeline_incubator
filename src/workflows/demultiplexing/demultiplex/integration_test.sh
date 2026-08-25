#!/bin/bash

set -eo pipefail

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

# NOTE: no `-stub -profile no_publish` smoke run here. `main.nf` parses
# `detect_sequencing_input`'s JSON output inline right after that component
# runs (`new groovy.json.JsonSlurper().parseText(state.detection_json.text)`)
# to decide which demultiplexer branch to take. Viash's auto-generated Nextflow
# stub just touches an empty placeholder for that output, so `-stub` fails
# before any real wiring is exercised ("Text must not be null or empty"),
# regardless of which entrypoint is used. No integration_test.sh in
# openpipeline v4.0.3/v4.1.0/v4.1.1 (nor openpipeline_qc/openpipeline_rapids/
# openpipeline_spatial) actually uses `-stub` either, so there is no working
# precedent to follow here.

nextflow \
  run . \
  -main-script src/workflows/demultiplexing/demultiplex/test.nf \
  -entry test_wf_fastq_passthrough \
  -profile docker,no_publish \
  -c src/workflows/utils/labels_ci.config \
  -c src/workflows/utils/integration_tests.config

nextflow \
  run . \
  -main-script src/workflows/demultiplexing/demultiplex/test.nf \
  -entry test_wf \
  -profile docker,no_publish \
  -c src/workflows/utils/labels_ci.config \
  -c src/workflows/utils/integration_tests.config

nextflow \
  run . \
  -main-script src/workflows/demultiplexing/demultiplex/test.nf \
  -entry test_wf_preset \
  -profile docker,no_publish \
  -c src/workflows/utils/labels_ci.config \
  -c src/workflows/utils/integration_tests.config

nextflow \
  run . \
  -main-script src/workflows/demultiplexing/demultiplex/test.nf \
  -entry test_wf_bases2fastq \
  -profile docker,no_publish \
  -c src/workflows/utils/labels_ci.config \
  -c src/workflows/utils/integration_tests.config

nextflow \
  run . \
  -main-script src/workflows/demultiplexing/demultiplex/test.nf \
  -entry test_wf_qc \
  -profile docker,no_publish \
  -c src/workflows/utils/labels_ci.config \
  -c src/workflows/utils/integration_tests.config
