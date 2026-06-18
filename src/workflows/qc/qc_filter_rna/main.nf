workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    // Flag mitochondrial genes in .var with a boolean `mt` column, used by
    // calculate_qc_metrics below to compute per-cell mitochondrial metrics.
    | grep_annotation_column.run(
      key: "grep_mitochondrial_genes",
      fromState: {id, state ->
        [
          "input": state.input,
          "modality": "rna",
          "input_layer": state.layer,
          "input_column": state.var_gene_names,
          "matrix": "var",
          "regex_pattern": state.mitochondrial_gene_regex,
          "output_match_column": state.output_match_column,
        ]
      },
      toState: ["input": "output"]
    )
    // Calculate total counts, log1p, top-gene proportions and the mitochondrial metrics
    // (pct_counts_mt, total_counts_mt, log1p_total_counts_mt) from the `mt` var column.
    | calculate_qc_metrics.run(
      key: "rna_calculate_qc_metrics",
      fromState: {id, state ->
        [
          "input": state.input,
          "modality": "rna",
          "layer": state.layer,
          "qc_vars": [state.output_match_column],
          "top_n_vars": state.top_n_vars,
          "log1p": true,
        ]
      },
      toState: ["input": "output"]
    )
    // Flag cells with too few total counts (empty droplets / low quality).
    | delimit_counts.run(
      key: "rna_delimit_total_counts",
      fromState: {id, state ->
        [
          "input": state.input,
          "modality": "rna",
          "obs_count_column": "total_counts",
          "obs_name_filter": state.obs_name_filter_total_counts,
          "min_count": state.min_total_counts,
          "max_count": state.max_total_counts,
        ]
      },
      toState: ["input": "output"]
    )
    // Flag cells outside the quantile range of the log1p total counts. Only run when at
    // least one quantile bound is provided. The `log1p_total_counts` column was already
    // log1p-transformed by calculate_qc_metrics, so disable the filter's own log1p step.
    | filter_with_quantile.run(
      key: "rna_quantile_filter_total_counts",
      runIf: {id, state ->
        state.max_log1p_total_counts_quantile != null || state.min_log1p_total_counts_quantile != null
      },
      fromState: {id, state ->
        [
          "input": state.input,
          "modality": "rna",
          "obs_column": "log1p_total_counts",
          "obs_log1p_transform": false,
          "obs_min_quantile": state.min_log1p_total_counts_quantile,
          "obs_max_quantile": state.max_log1p_total_counts_quantile,
          "obs_name_filter": state.obs_name_filter_quantile,
        ]
      },
      toState: ["input": "output"]
    )
    // Flag cells with a mitochondrial percentage above the threshold.
    | delimit_counts.run(
      key: "rna_delimit_pct_counts_mt",
      fromState: {id, state ->
        [
          "input": state.input,
          "modality": "rna",
          "obs_count_column": "pct_counts_${state.output_match_column}",
          "min_count": state.min_pct_counts_mitochondrial,
          "max_count": state.max_pct_counts_mitochondrial,
          "obs_name_filter": state.obs_name_filter_mitochondrial,
        ]
      },
      toState: ["input": "output"]
    )
    // Doublet detection (flag only). Skipped when --skip_scrublet is set.
    | filter_with_scrublet.run(
      key: "rna_scrublet_doublet_detection",
      runIf: {id, state -> !state.skip_scrublet},
      fromState: {id, state ->
        [
          "input": state.input,
          "modality": "rna",
          "layer": state.layer,
          "obs_name_filter": state.obs_name_filter_scrublet,
          "obs_name_doublet_score": state.obs_name_doublet_score,
          "obs_name_predicted_doublets": state.obs_name_predicted_doublets,
          "scrublet_score_threshold": state.scrublet_score_threshold,
          "expected_doublet_rate": state.scrublet_expected_doublet_rate,
          "min_counts": state.scrublet_min_counts,
          "min_cells": state.scrublet_min_cells,
          "min_gene_variablity_percent": state.scrublet_min_gene_variability_percent,
          "num_pca_components": state.scrublet_num_pca_components,
        ]
      },
      args: [output_compression: "gzip"],
      toState: ["input": "output"]
    )
    // Emit the flagged (but not subsetted) data. Subsetting and cross-modality
    // intersection are handled by the top-level qc_filter workflow.
    | setState(["output": "input"])

  emit:
  output_ch
}
