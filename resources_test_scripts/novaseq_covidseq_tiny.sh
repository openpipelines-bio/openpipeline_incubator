#!/bin/bash

set -eo pipefail

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

# ---------------------------------------------------------------------------
# GEX fixture: a single-index Illumina run, for the demultiplex workflow's
# test_wf (bcl_convert, no preset).
#
# This is nf-core/demultiplex's NovaSeq6000 COVIDSeq test dataset. RunInfo.xml=
# declares one 36bp genomic read and two 10bp indexed reads (single-end).
# ---------------------------------------------------------------------------
ID=novaseq_covidseq_tiny
OUT="resources_test/$ID"

if [ ! -f "$OUT/SampleSheet.csv" ]; then
  mkdir -p "$OUT"

  TMPDIR=$(mktemp -d)
  trap 'rm -rf "$TMPDIR"' EXIT

  wget -q https://github.com/nf-core/test-datasets/raw/demultiplex/testdata/NovaSeq6000/200624_A00834_0183_BHMTFYDRXX.tar.gz \
    -O "$TMPDIR/200624_A00834_0183_BHMTFYDRXX.tar.gz"
  tar -xf "$TMPDIR/200624_A00834_0183_BHMTFYDRXX.tar.gz" --strip-components=1 -C "$OUT"

  # Create a v2 SampleSheet, the dataset comes with v1
  cat > "$OUT/SampleSheet.csv" << SAMPLESHEET
[Header]
FileFormatVersion,2

[Reads]
Read1Cycles,36
Index1Cycles,10
Index2Cycles,10

[BCLConvert_Settings]
AdapterRead1,CTGTCTCTTATACACATCT

[BCLConvert_Data]
Lane,Sample_ID,index,index2
1,Sample1,GAACTGAGCG,TCGTGGAGCG
1,SampleA,AGGTCAGATA,CTACAAGATA
1,Sample23,CGTCTCATAT,TATAGTAGCT
1,sampletest,ATTCCATAAG,TGCCTGGTGG
SAMPLESHEET
fi

aws s3 sync \
  "$OUT" \
  s3://openpipelines-bio/openpipeline_incubator/resources_test/"$ID" \
  --delete \
  --dryrun
