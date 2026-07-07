workflow run_wf {
  take:
  input_ch

  main:

  // Maps each modality to the sub-workflow arguments it needs, keyed by the sub-workflow
  // argument name and pointing to the top-level state key holding the value. Adding a new
  // modality is a matter of adding a qc_filter_<modality> sub-workflow and an entry here.
  def modality_arguments = [
    "rna": [
      "layer": "rna_layer",
      "min_total_counts": "rna_min_total_counts",
      "max_total_counts": "rna_max_total_counts",
      "min_log1p_total_counts_quantile": "rna_min_log1p_total_counts_quantile",
      "max_log1p_total_counts_quantile": "rna_max_log1p_total_counts_quantile",
      "min_pct_counts_mitochondrial": "min_pct_counts_mitochondrial",
      "max_pct_counts_mitochondrial": "max_pct_counts_mitochondrial",
      "top_n_vars": "rna_top_n_vars",
      "mitochondrial_gene_regex": "mitochondrial_gene_regex",
      "var_gene_names": "var_gene_names",
      "skip_scrublet": "skip_scrublet",
      "scrublet_score_threshold": "scrublet_score_threshold",
      "scrublet_expected_doublet_rate": "scrublet_expected_doublet_rate",
      "scrublet_min_counts": "scrublet_min_counts",
      "scrublet_min_cells": "scrublet_min_cells",
      "scrublet_min_gene_variability_percent": "scrublet_min_gene_variability_percent",
      "scrublet_num_pca_components": "scrublet_num_pca_components",
    ],
    "prot": [
      "layer": "prot_layer",
      "min_total_counts": "prot_min_total_counts",
      "max_total_counts": "prot_max_total_counts",
      "min_log1p_total_counts_quantile": "prot_min_log1p_total_counts_quantile",
      "max_log1p_total_counts_quantile": "prot_max_log1p_total_counts_quantile",
      "top_n_vars": "prot_top_n_vars",
    ],
  ].asImmutable()

  split_ch = input_ch
    // Stamp the sample id onto a .obs column so downstream grouping (e.g. the cell-count
    // report) can identify which sample each cell belongs to.
    | add_id.run(
      fromState: {id, state ->
        [
          "input": state.input,
          "input_id": id,
          "obs_output": state.sample_id_column,
        ]
      },
      toState: ["input": "output"]
    )
    // Split the multimodal input into one file per modality for separate processing.
    | split_modalities.run(
      fromState: {id, state -> ["input": state.input]},
      toState: ["split_output": "output", "split_types": "output_types"]
    )
    // Expand the split directory + types CSV into one event per modality.
    | flatMap {id, state ->
      def outputDir = state.split_output
      def csv = state.split_types.splitCsv(strip: true, sep: ",").findAll{!it[0].startsWith("#")}
      def header = csv.head()
      def types = csv.tail().collect { row -> [header, row].transpose().collectEntries() }
      types.collect{ dat ->
        [id, state + ["input": outputDir.resolve(dat.filename), "modality": dat.name, "_meta": ["join_id": id]]]
      }
    }
    | view {"After splitting modalities: $it"}

  // Known modalities (rna, prot) get their QC + doublet keep-flags computed per modality.
  known_ch = split_ch
    | runEach(
      components: [qc_filter_rna, qc_filter_prot],
      filter: {id, state, component -> "qc_filter_" + state.modality == component.config.name},
      fromState: {id, state, component ->
        def args = modality_arguments.get(state.modality).collectEntries{key_, value_ -> [key_, state[value_]]}
        return args + ["id": id, "input": state.input]
      },
      toState: ["input": "output"]
    )

  // Modalities without a dedicated sub-workflow pass through unchanged.
  unknown_ch = split_ch
    | filter {id, state -> !modality_arguments.containsKey(state.modality)}

  filtered_ch = known_ch.mix(unknown_ch)
    // Regroup the per-modality files of each sample.
    | map {id, state -> [state._meta.join_id, state]}
    | groupTuple(by: 0, sort: "hash")
    | map { id, states ->
        def new_input = states.collect{it.input}
        def modalities = states.collect{it.modality}.unique()
        def other_state_keys = states.inject([].toSet()){ current_keys, state ->
          current_keys + state.keySet()
        }.minus(["input", "modality", "_meta"])
        def new_state = other_state_keys.inject([:]){ old_state, argument_name ->
          def argument_values = states.collect{it.get(argument_name)}.unique()
          assert argument_values.size() == 1, "Arguments should be the same across modalities. \
                                               Argument name: $argument_name, values: $argument_values"
          old_state + [(argument_name): argument_values[0]]
        }
        [id, new_state + ["input": new_input, "modalities": modalities]]
    }
    // Merge the per-modality files back into a single multimodal object. All keep-flags are
    // preserved in each modality's .obs; no cells have been dropped yet.
    | merge.run(
      fromState: ["input": "input"],
      args: [output_compression: "gzip"],
      toState: ["input": "output"]
    )
    // Stash the merged, still-unsubset file (carrying every keep-flag) so the cell-count
    // report can read per-stage counts after the chain below has subset `input`.
    | map {id, state -> [id, state + ["flagged_input": state.input]]}
    // Apply the RNA cell keep-flags.
    | do_filter.run(
      key: "rna_filter_cells",
      runIf: {id, state -> state.modalities.contains("rna")},
      fromState: {id, state ->
        def obs_filter = ["filter_counts_rna", "filter_quantile_rna", "filter_mito_rna"]
        if (!state.skip_scrublet) {
          obs_filter += ["filter_scrublet"]
        }
        [
          "input": state.input,
          "modality": "rna",
          "obs_filter": obs_filter,
        ]
      },
      toState: ["input": "output"]
    )
    // Flag and remove RNA genes expressed in too few cells (evaluated after cell filtering).
    | filter_with_counts.run(
      key: "rna_gene_filter",
      runIf: {id, state -> state.modalities.contains("rna") && state.filter_genes_min_cells != null},
      fromState: {id, state ->
        [
          "input": state.input,
          "modality": "rna",
          "layer": state.rna_layer,
          "var_name_filter": "filter_genes_rna",
          "min_cells_per_gene": state.filter_genes_min_cells,
        ]
      },
      toState: ["input": "output"]
    )
    | do_filter.run(
      key: "rna_filter_genes",
      runIf: {id, state -> state.modalities.contains("rna") && state.filter_genes_min_cells != null},
      fromState: {id, state ->
        [
          "input": state.input,
          "modality": "rna",
          "var_filter": ["filter_genes_rna"],
        ]
      },
      toState: ["input": "output"]
    )
    // Apply the protein cell keep-flags.
    | do_filter.run(
      key: "prot_filter_cells",
      runIf: {id, state -> state.modalities.contains("prot")},
      fromState: {id, state ->
        [
          "input": state.input,
          "modality": "prot",
          "obs_filter": ["filter_counts_prot", "filter_quantile_prot"],
        ]
      },
      toState: ["input": "output"]
    )
    // Keep only cells passing every filter in every QC-filtered modality. Only modalities with
    // a qc_filter_<modality> sub-workflow participate; pass-through modalities (e.g. vdj) are
    // excluded so their barcode set never narrows the result.
    | intersect_obs.run(
      runIf: {id, state -> state.intersect_obs && state.modalities.count{modality_arguments.containsKey(it)} > 1},
      fromState: {id, state ->
        [
          "input": state.input,
          "modalities": state.modalities.findAll{modality_arguments.containsKey(it)},
        ]
      },
      toState: ["input": "output"]
    )

  // Barrier: gather every sample into a single event holding the list of filtered files
  // (`filtered_inputs`) and the list of still-unsubset, flag-carrying files (`flagged_inputs`).
  // The workflow output arguments and report knobs are identical across samples.
  combined_ch = filtered_ch
    | joinStates { ids, states ->
      def combined_state = [
        "filtered_inputs": states.collect{it.input},
        "flagged_inputs": states.collect{it.flagged_input},
        // Per-sample event ids, aligned with flagged_inputs, used as the cell-count report's
        // sample id when the data carries no sample-id .obs column.
        "sample_ids": ids,
        "modalities": states.collect{it.modalities}.flatten().unique(),
        "sample_id_column": states[0].sample_id_column,
        "skip_scrublet": states[0].skip_scrublet,
        "output_processed_h5mu": states[0].output_processed_h5mu,
        "cell_count_report": states[0].cell_count_report,
        "_meta": ["join_id": ids[0]],
      ]
      ["combined", combined_state]
    }

  // Collect the filtered files into a single output directory, then write a single cell-count
  // report across all samples from the unsubset keep-flags. Both run off the one joinStates
  // barrier above, chained linearly (disjoint output keys, no fork to rejoin).
  output_ch = combined_ch
    | move_files_to_directory.run(
      fromState: {id, state ->
        [
          "input": state.filtered_inputs,
          "output": state.output_processed_h5mu,
        ]
      },
      toState: ["output_processed_h5mu": "output"]
    )
    | cell_count_report.run(
      fromState: {id, state ->
        def args = [
          "input": state.flagged_inputs,
          "sample_id": state.sample_ids,
          "sample_id_column": state.sample_id_column,
          "rna_filter_columns": ["filter_counts_rna", "filter_quantile_rna", "filter_mito_rna"],
          "output": state.cell_count_report,
          "modality": "rna"
        ]
        if (state.modalities.contains("prot")) {
          args.prot_modality = "prot"
          args.prot_filter_columns = ["filter_counts_prot", "filter_quantile_prot"]
        }
        if (!state.skip_scrublet) {
          args.scrublet_filter_column = "filter_scrublet"
        }
        args
      },
      toState: ["cell_count_report": "output"]
    )
    | setState(["output_processed_h5mu", "cell_count_report", "_meta"])

  emit:
  output_ch
}
