openpipeline_incubator x.x.x (unreleased)

## New functionality

### Components

* `components/qc/calculate_qc_metrics_scanpy`: add per-cell/per-gene QC metrics to an `.h5mu`
  modality using `scanpy.pp.calculate_qc_metrics` (PR #1).

* `components/report/cell_count_report`: read the per-modality boolean keep-columns and emit a
  per-sample TSV of the cells surviving each filter stage (PR #1).

* `components/mapping/single2spatial`: deep-learning mapper (based on the Bulk2Space deep-forest
  algorithm) that projects single-cell RNA-seq profiles onto spatial coordinates, taking a
  single-cell and a spatial AnnData as input (PR #4).

### Workflows

* `workflows/qc/qc_filter_rna`: single-sample, unimodal RNA sub-workflow that computes QC and
  scrublet doublet keep-flags in flag-only mode (no subsetting) (PR #1).

* `workflows/qc/qc_filter_prot`: single-sample, unimodal protein sub-workflow that computes QC
  keep-flags in flag-only mode (PR #1).

* `workflows/processing/process_single_sample`: multimodal (RNA + protein) orchestrator that
  splits modalities, runs each through its `qc_filter_<modality>` sub-workflow, merges them back,
  writes the cell-count report, applies the keep-flags, filters rare RNA genes, and intersects
  observations across modalities (PR #1, #3).
