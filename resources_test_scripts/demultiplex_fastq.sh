#!/bin/bash

set -eo pipefail

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

ID=demultiplex_fastq
OUT=resources_test/$ID
DIR="$OUT"

[ -d "$DIR" ] || mkdir -p "$DIR"

# Synthesize a tiny FASTQ-only fixture for the demultiplex workflow's passthrough path.
# create_fastq_manifest only inspects filenames (and file size, to warn on empty files), never
# FASTQ content, so a short dummy record is enough; real gzip compression is not required.
# Covers:
# - sample1: multi-lane (L001, L002), paired-end
# - sample2: single-end (R1 only, no R2)
# - Undetermined_S0: excluded from the manifest by default (--include_undetermined not set)
write_fastq () {
  cat > "$1" <<FASTQ
@read1
ACGTACGTAC
+
IIIIIIIIII
FASTQ
}

write_fastq "$DIR/sample1_S1_L001_R1_001.fastq.gz"
write_fastq "$DIR/sample1_S1_L001_R2_001.fastq.gz"
write_fastq "$DIR/sample1_S1_L002_R1_001.fastq.gz"
write_fastq "$DIR/sample1_S1_L002_R2_001.fastq.gz"
write_fastq "$DIR/sample2_S2_L001_R1_001.fastq.gz"
write_fastq "$DIR/Undetermined_S0_L001_R1_001.fastq.gz"
write_fastq "$DIR/Undetermined_S0_L001_R2_001.fastq.gz"

aws s3 sync \
  "$DIR" \
  s3://openpipelines-bio/openpipeline_incubator/resources_test/"$ID" \
  --delete \
  --dryrun
