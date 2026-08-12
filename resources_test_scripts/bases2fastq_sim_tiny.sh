#!/bin/bash

set -eo pipefail

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

# ---------------------------------------------------------------------------
# A Element Biosciences simulated run
# ---------------------------------------------------------------------------
ID=bases2fastq_sim_tiny
OUT="resources_test/$ID"

if [ ! -f "$OUT/RunManifest.csv" ]; then
  mkdir -p "$OUT"

  TMPDIR=$(mktemp -d)
  trap 'rm -rf "$TMPDIR"' EXIT

  wget -q http://element-public-data.s3.amazonaws.com/bases2fastq-share/bases2fastq-v2/20230404-bases2fastq-sim-151-151-9-9.tar.gz \
    -O "$TMPDIR/20230404-bases2fastq-sim-151-151-9-9.tar.gz"
  tar -xf "$TMPDIR/20230404-bases2fastq-sim-151-151-9-9.tar.gz" --strip-components=1 -C "$OUT"
fi

aws s3 sync \
  "$OUT" \
  s3://openpipelines-bio/openpipeline_incubator/resources_test/"$ID" \
  --delete \
  --dryrun
