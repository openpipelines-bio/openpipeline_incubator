#!/usr/bin/env bash

# Fetch the bulk2single tutorial's own datasets (dentate gyrus bulk + single-cell
# RNA-seq, from https://omicverse.readthedocs.io/en/latest/Tutorials-bulk2single/t_bulk2single/)
# so the component can be tested against the same real data the algorithm was
# demonstrated on, instead of synthetic/pseudobulk data.
#
# Run from the root of the repository.
set -euo pipefail

DATA_DIR="./resources_test/bulk2single"
mkdir -p "$DATA_DIR"
cd "$DATA_DIR"

# Bulk RNA-seq count matrix (GEO GSE74985): genes (Ensembl mouse gene IDs) in
# rows, samples in columns. Downloaded gzipped, then decompressed to a plain
# .tsv since --bulk_data reads it directly with no decompression step.
BULK_URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE74nnn/GSE74985/suppl/GSE74985_mergedCount.txt.gz"
if [[ -f "input_bulk.tsv" ]]; then
  echo "Skipping input_bulk.tsv (already present)"
else
  echo "Downloading bulk count matrix ..."
  curl -fsSL "$BULK_URL" -o "GSE74985_mergedCount.txt.gz"
  gunzip -c "GSE74985_mergedCount.txt.gz" > "input_bulk.tsv"
  rm "GSE74985_mergedCount.txt.gz"
fi

# Single-cell reference (dentate gyrus neurogenesis, adjusted/preprocessed
# variant of scvelo's scv.datasets.dentategyrus()) as AnnData - gene symbols in
# .var, cluster/cell-type labels in .obs['clusters'].
SINGLE_CELL_URL="https://github.com/theislab/scvelo_notebooks/raw/master/data/DentateGyrus/10X43_1.h5ad"
if [[ -f "input_sc.h5ad" ]]; then
  echo "Skipping input_sc.h5ad (already present)"
else
  echo "Downloading single-cell reference ..."
  curl -fsSL "$SINGLE_CELL_URL" -o "input_sc.h5ad"
fi

# Ensembl gene ID <-> gene symbol mapping pair (GRCm39/mouse), needed because
# the bulk matrix is indexed by Ensembl IDs while the single-cell reference
# uses gene symbols. Source: omicverse's own
# ov.utils.download_geneid_annotation_pair() (Stanford mirror).
GENE_ID_MAPPING_URL="https://stacks.stanford.edu/file/cv694yk7414/pair_GRCm39.tsv"
if [[ -f "pair_GRCm39.tsv" ]]; then
  echo "Skipping pair_GRCm39.tsv (already present)"
else
  echo "Downloading gene ID mapping pair ..."
  curl -fsSL "$GENE_ID_MAPPING_URL" -o "pair_GRCm39.tsv"
fi

echo "Test data in $(pwd)"
echo "Done."