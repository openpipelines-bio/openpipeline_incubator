workflow run_wf {
  take:
    input_ch

  main:
    output_ch = input_ch
    | map { id, state ->
      def new_state = state + [ "query_processed": state.output, "_meta": ["join_id": id] ]
      [id, new_state]
    }
    | process_samples_workflow.run(
      fromState: {id, state ->
        def newState = [
          "input": state.input, 
          "id": id,
          "rna_layer": state.input_layer,
          "rna_min_counts": state.rna_min_counts,
          "rna_max_counts": state.rna_max_counts,
          "rna_min_genes_per_cell": state.rna_min_genes_per_cell,
          "rna_max_genes_per_cell": state.rna_max_genes_per_cell,
          "rna_min_cells_per_gene": state.rna_min_cells_per_gene,
          "rna_min_fraction_mito": state.rna_min_fraction_mito,
          "rna_max_fraction_mito": state.rna_max_fraction_mito,
          "highly_variable_features_var_output": state.highly_variable_features_var_output,
          "highly_variable_features_obs_batch_key": state.highly_variable_features_obs_batch_key,
          "var_name_mitochondrial_genes": state.var_name_mitochondrial_genes,
          "var_gene_names": state.input_var_gene_names,
          "mitochondrial_gene_regex": state.mitochondrial_gene_regex,
          "var_qc_metrics": state.var_qc_metrics,
          "top_n_vars": state.top_n_vars,
        ]  
      },
      args: [
        "pca_overwrite": "true",
        "add_id_obs_output": "sample_id"
      ],
      toState: ["query_processed": "output"], 
    )

    | view {"After processing query: $it"}

    | scgpt_annotation_workflow.run(
      runIf: { id, state -> state.annotation_methods.contains("scgpt_annotation") },
      fromState: { id, state ->
        [ 
          "id": id,
          "input": state.query_processed,
          "modality": state.modality,
          "input_layer": state.input_layer,
          "input_var_gene_names": state.input_var_gene_names,
          "model": state.scgpt_model,
          "model_config": state.scgpt_model_config,
          "model_vocab": state.scgpt_model_vocab,
          "finetuned_checkpoints_key": state.scgpt_finetuned_checkpoints_key,
          "label_mapper_key": state.scgpt_label_mapper_key,
          "pad_token": state.scgpt_pad_token,
          "pad_value": state.scgpt_pad_value,
          "n_hvg": state.scgpt_n_hvg,
          "dsbn": state.scgpt_dsbn,
          "batch_size": state.scgpt_batch_size,
          "n_input_bins": state.scgpt_n_input_bins,
          "seed": state.scgpt_seed
        ]
      },
      args: [
        "input_obs_batch_label": "sample_id",
        "output_obs_predictions": "scgpt_pred",
        "output_obs_probability": "scgpt_proba"
      ],
      toState: [ "query_processed": "output" ]
    )

    | view {"After scgpt: $it"}

    // | celltypist.run(
    //   runIf: { id, state -> state.annotation_methods.contains("celltypist") && state.celltypist_model },
    //   fromState: [ 
    //     "input": "query_processed",
    //     "modality": "modality",
    //     "input_layer": "input_layer",
    //     "input_var_gene_names": "input_var_gene_names",
    //     "input_reference_gene_overlap": "input_reference_gene_overlap",
    //     "model": "celltypist_model",
    //     "majority_voting": "celltypist_majority_voting"
    //   ],
    //   args: [
    //     "output_obs_predictions": "celltypist_pred",
    //     "output_obs_probability": "celltypist_proba"
    //   ],
    //   toState: [ "query_processed": "output" ]
    // )

    // | view {"After celltypist: $it"}
    | map {id, state ->
      def new_state = state + ["output": state.query_processed]
      [id, new_state]
    }
    | view {"After mapping: $it"}
    | setState(["output", "_meta"])
    | view {"After setstate: $it"}

  emit:
    output_ch
}