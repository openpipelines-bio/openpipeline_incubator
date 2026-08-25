#!/bin/bash

set -eo pipefail

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

# ---------------------------------------------------------------------------
# ATAC: a dual-index Illumina run (two genomic + two indexed reads),
# for the demultiplex workflow's test_wf_preset (bcl_convert with
# --preset 10x_atac).
#
# This is 10x's BCL-Convert-native ARC ATAC tiny-bcl dataset. RunInfo.xml
# declares two genomic reads (50bp/49bp) and two indexed reads (8bp/24bp).
# ---------------------------------------------------------------------------
ID=cellranger_arc_tiny_bcl_atac
OUT="resources_test/$ID"

if [ ! -f "$OUT/SampleSheet.csv" ]; then
  mkdir -p "$OUT"

  TMPDIR=$(mktemp -d)
  trap 'rm -rf "$TMPDIR"' EXIT

  wget -q https://cf.10xgenomics.com/supp/cell-arc/cellranger-arc-tiny-bcl-atac-2.0.0.tar.gz \
    -O "$TMPDIR/cellranger-arc-tiny-bcl-atac-2.0.0.tar.gz"
  tar -xf "$TMPDIR/cellranger-arc-tiny-bcl-atac-2.0.0.tar.gz" --strip-components=1 -C "$OUT"

  rm -f "$OUT/IlluminaSampleSheet.csv" \
    "$OUT/Data/Intensities/BaseCalls/SampleSheet.csv" \
    "$OUT/Data/Intensities/BaseCalls/IlluminaSampleSheet.csv"

  # Create a v2 SampleSheet, the dataset comes with v1
  # Don't include [BCLConvert_Settings] to test adding this in the workflow
  cat > "$OUT/SampleSheet.csv" << SAMPLESHEET
[Header]
FileFormatVersion,2

[BCLConvert_Data]
Lane,Sample_ID,index
1,test_sample_atac,TTGTAAGA
1,test_sample_atac,GGCGTTTC
1,test_sample_atac,CCTACCAT
1,test_sample_atac,AAACGGCG
SAMPLESHEET
fi

aws s3 sync \
  "$OUT" \
  s3://openpipelines-bio/openpipeline_incubator/resources_test/"$ID" \
  --delete \
  --dryrun
