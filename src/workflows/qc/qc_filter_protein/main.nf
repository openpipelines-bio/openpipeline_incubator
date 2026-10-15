workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    // Calculate the total counts (and its log1p) per cell, used by the count and
    // quantile filters below.
    | calculate_qc_metrics.run(
      key: "prot_calculate_qc_metrics",
      fromState: {id, state ->
        [
          "input": state.input,
          "modality": "prot",
          "layer": state.layer,
          "top_n_vars": state.top_n_vars,
          "log1p": true,
        ]
      },
      toState: ["input": "output"]
    )
    // Flag cells with too few total counts.
    | delimit_counts.run(
      key: "prot_delimit_total_counts",
      fromState: {id, state ->
        [
          "input": state.input,
          "modality": "prot",
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
      key: "prot_quantile_filter_total_counts",
      runIf: {id, state ->
        state.max_log1p_total_counts_quantile != null || state.min_log1p_total_counts_quantile != null
      },
      fromState: {id, state ->
        [
          "input": state.input,
          "modality": "prot",
          "obs_column": "log1p_total_counts",
          "obs_log1p_transform": false,
          "obs_min_quantile": state.min_log1p_total_counts_quantile,
          "obs_max_quantile": state.max_log1p_total_counts_quantile,
          "obs_name_filter": state.obs_name_filter_quantile,
        ]
      },
      toState: ["input": "output"]
    )
    // Emit the flagged (but not subsetted) data. Subsetting and cross-modality
    // intersection are handled by the top-level qc_filter workflow.
    | setState(["output": "input"])

  emit:
  output_ch
}
