workflow run_wf {
  take:
    input_ch

  main:
    output_ch = input_ch
    | map { id, state ->
      def new_state = state + [ "query_processed": state.output, "_meta": ["join_id": id] ]
      [id, new_state]
    }
    // Enforce annotation method-specific required arguments
    | map { id, state ->
      def new_state = [:]
      // Check scGPT arguments
      if (state.annotation_methods.contains("scgpt_annotation") && 
        (!state.scgpt_model || !state.scgpt_model_config || !state.scgpt_model_vocab)) {
        throw new RuntimeException("Using scgpt_annotation requires --scgpt_model, --scgpt_model_config and --scgp_model_vocab parameters.")
        }
      // Check CellTypist arguments
      if (state.annotation_methods.contains("celltypist") && 
        (!state.celltypist_model && !state.reference)) {
        throw new RuntimeException("Celltypist was selected as an annotation method. Either --celltypist_model or --reference must be provided.")
        }
      if (state.annotation_methods.contains("celltypist") && state.celltypist_model && state.reference )  {
        System.err.println(
          "Warning: --celltypist_model is set and a --reference was provided. \
          The pre-trained Celltypist model will be used for annotation, the reference will be ignored."
          )
        }
      // Check Harmony KNN arguments
      if ((state.annotation_methods.contains("harmony_knn") || state.annotation_methods.contains("scvi_knn"))  && !state.reference ) {
        throw new RuntimeException("When `harmony_knn` or `scvi_knn` are selected as an annotation method, a --reference dataset must be provided.")

      [id, state + new_state]
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
          "var_name_mitochondrial_genes": state.var_name_mitochondrial_genes,
          "var_gene_names": state.input_var_gene_names,
          "mitochondrial_gene_regex": state.mitochondrial_gene_regex,
          "var_qc_metrics": state.var_qc_metrics
        ]  
      },
      args: [
        "pca_overwrite": "true",
        "add_id_obs_output": "sample_id"
      ],
      toState: ["query_processed": "output"], 
    )

    | scgpt_annotation.run(
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

    | celltypist_annotation.run(
      runIf: { id, state -> state.annotation_methods.contains("celltypist") && state.celltypist_model },
      fromState: [ 
        "input": "query_processed",
        "modality": "modality",
        "input_var_gene_names": "input_var_gene_names",
        "input_reference_gene_overlap": "input_reference_gene_overlap",
        "model": "celltypist_model",
        "majority_voting": "celltypist_majority_voting"
      ],
      args: [
        // log normalized counts are expected for celltypist
        "input_layer": "log_normalized",
        "output_obs_predictions": "celltypist_pred",
        "output_obs_probability": "celltypist_proba"
      ],
      toState: [ "query_processed": "output" ]
    )

    | celltypist_annotation.run(
      runIf: { id, state -> state.annotation_methods.contains("celltypist") && !state.celltypist_model },
      fromState: [
        "input": "query_processed",
        "modality": "modality",
        "input_var_gene_names": "input_var_gene_names",
        "input_reference_gene_overlap": "input_reference_gene_overlap",
        "reference": "reference",
        "reference_layer": "reference_layer_lognormalized_counts",
        "reference_obs_target": "reference_obs_label",
        "reference_var_gene_names": "reference_var_gene_names",
        "reference_obs_batch": "reference_obs_batch",
        "reference_var_input": "reference_var_input",
        "feature_selection": "celltypist_feature_selection",
        "C": "celltypist_C",
        "max_iter": "celltypist_max_iter",
        "use_SGD": "celltypist_use_SGD",
        "min_prop": "celltypist_min_prop",
        "majority_voting": "celltypist_majority_voting"
      ],
      args: [
        // log normalized counts are expected for celltypist
        "input_layer": "log_normalized",
        "output_obs_predictions": "celltypist_pred",
        "output_obs_probability": "celltypist_proba"
      ],
      toState: [ "query_processed": "output" ]
    )

    | harmony_knn_annotation.run(
      runIf: { id, state -> state.annotation_methods.contains("harmony_knn") },
      fromState: { id, state ->
        [ 
          "id": id,
          "input": state.query_processed,
          "modality": state.modality,
          "input_var_gene_names": state.input_var_gene_names,
          "input_reference_gene_overlap": state.input_reference_gene_overlap,
          "reference": state.reference,
          "reference_layer": state.reference_layer_lognormalized_counts,
          "reference_obs_target": state.reference_obs_label,
          "reference_var_gene_names": state.reference_var_gene_names,
          "reference_obs_batch_label": state.reference_obs_batch,
          "n_hvg": state.n_hvg,
          "harmony_theta": state.harmony_theta,
        ]
      },
      args: [
        "input_layer": "log_normalized",
        "input_obs_batch_label": "sample_id",
        "output_obs_predictions": "harmony_knn_pred",
        "output_obs_probability": "harmony_knn_proba",
        "output_obsm_integrated": "X_integrated_harmony",
        "overwrite_existing_key": "true"
      ],
      toState: [ "query_processed": "output" ]
    )

    | scvi_knn_annotation.run(
      runIf: { id, state -> state.annotation_methods.contains("harmony_knn") },
      fromState: { id, state ->
        [ 
          "id": id,
          "input": state.query_processed,
          "modality": state.modality,
          "input_layer": state.input_layer,
          "input_var_gene_names": state.input_var_gene_names,
          "input_reference_gene_overlap": state.input_reference_gene_overlap,
          "reference": state.reference,
          "reference_layer": state.reference_layer_raw_counts,
          "reference_layer_lognormalized": state.reference_layer_lognormalized_counts,
          "reference_obs_target": state.reference_obs_label,
          "reference_var_gene_names": state.reference_var_gene_names,
          "reference_obs_batch_label": state.reference_obs_batch,
          "n_hvg": state.n_hvg,
          "scvi_early_stopping": state.scvi_early_stopping,
          "scvi_early_stopping_patience": state.scvi_early_stopping_patience,
          "scvi_early_stopping_min_delta": state.scvi_early_stopping_min_delta,
          "scvi_max_epochs": state.scvi_max_epochs,
          "scvi_reduce_lr_on_plateau": state.scvi_reduce_lr_on_plateau,
          "scvi_lr_factor": state.scvi_lr_factor,
          "scvi_lr_patience": state.scvi_lr_patience
        ]
      },
      args: [
        "input_layer_lognormalized": "log_normalized",
        "input_obs_batch_label": "sample_id",
        "output_obs_predictions": "harmony_knn_pred",
        "output_obs_probability": "harmony_knn_proba",
        "output_obsm_integrated": "X_integrated_harmony",
        "overwrite_existing_key": "true"
      ],
      toState: [ "query_processed": "output" ]
    )

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