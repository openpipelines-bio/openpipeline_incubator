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
      }

      [id, state + new_state]
    }
    | process_samples_workflow.run(
      fromState: [
        "input": "input", 
        "id": "id",
        "rna_layer": "input_layer",
        "rna_min_counts": "rna_min_counts",
        "rna_max_counts": "rna_max_counts",
        "rna_min_genes_per_cell": "rna_min_genes_per_cell",
        "rna_max_genes_per_cell": "rna_max_genes_per_cell",
        "rna_min_cells_per_gene": "rna_min_cells_per_gene",
        "rna_min_fraction_mito": "rna_min_fraction_mito",
        "rna_max_fraction_mito": "rna_max_fraction_mito",
        "rna_min_fraction_ribo": "rna_min_fraction_ribo",
        "rna_max_fraction_ribo": "rna_max_fraction_ribo",
        "var_name_mitochondrial_genes": "var_name_mitochondrial_genes",
        "var_name_ribosomal_genes": "var_name_ribosomal_genes",
        "var_gene_names": "input_var_gene_names",
        "mitochondrial_gene_regex": "mitochondrial_gene_regex",
        "ribosomal_gene_regex": "ribosomal_gene_regex",
        "var_qc_metrics": "var_qc_metrics"
      ],
      args: [
        "pca_overwrite": "true",
        "add_id_obs_output": "sample_id"
      ],
      toState: ["query_processed": "output"], 
    )

    | scgpt_annotation.run(
      runIf: { id, state -> state.annotation_methods.contains("scgpt_annotation") },
      fromState: [ 
        "id": "id",
        "input": "query_processed",
        "modality": "modality",
        "input_var_gene_names": "input_var_gene_names",
        "model": "scgpt_model",
        "model_config": "scgpt_model_config",
        "model_vocab": "scgpt_model_vocab",
        "finetuned_checkpoints_key": "scgpt_finetuned_checkpoints_key",
        "label_mapper_key": "scgpt_label_mapper_key",
        "pad_token": "scgpt_pad_token",
        "pad_value": "scgpt_pad_value",
        "n_hvg": "scgpt_n_hvg",
        "dsbn": "scgpt_dsbn",
        "batch_size": "scgpt_batch_size",
        "n_input_bins": "scgpt_n_input_bins",
        "seed": "scgpt_seed"
      ],
      args: [
        "input_layer": "log_normalized",
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
      fromState: [ 
        "id": "id",
        "input": "query_processed",
        "modality": "modality",
        "input_var_gene_names": "input_var_gene_names",
        "input_reference_gene_overlap": "input_reference_gene_overlap",
        "reference": "reference",
        "reference_layer": "reference_layer_lognormalized_counts",
        "reference_obs_target": "reference_obs_label",
        "reference_var_gene_names": "reference_var_gene_names",
        "reference_obs_batch_label": "reference_obs_batch",
        "n_hvg": "n_hvg",
        "harmony_theta": "harmony_theta",
        "leiden_resolution": "leiden_resolution",
        "knn_weights": "knn_weights",
        "knn_n_neighbors": "knn_n_neighbors"
      ],
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
      fromState: [ 
        "id": "id",
        "input": "query_processed",
        "modality": "modality",
        "input_layer": "input_layer",
        "input_var_gene_names": "input_var_gene_names",
        "input_reference_gene_overlap": "input_reference_gene_overlap",
        "reference": "reference",
        "reference_layer": "reference_layer_raw_counts",
        "reference_layer_lognormalized": "reference_layer_lognormalized_counts",
        "reference_obs_target": "reference_obs_label",
        "reference_var_gene_names": "reference_var_gene_names",
        "reference_obs_batch_label": "reference_obs_batch",
        "n_hvg": "n_hvg",
        "early_stopping": "early_stopping",
        "early_stopping_patience": "early_stopping_patience",
        "early_stopping_min_delta": "early_stopping_min_delta",
        "max_epochs": "max_epochs",
        "reduce_lr_on_plateau": "reduce_lr_on_plateau",
        "lr_factor": "lr_factor",
        "lr_patience": "lr_patience",
        "leiden_resolution": "leiden_resolution",
        "knn_weights": "knn_weights",
        "knn_n_neighbors": "knn_n_neighbors"
      ],
      args: [
        "input_layer_lognormalized": "log_normalized",
        "input_obs_batch_label": "sample_id",
        "output_obs_predictions": "scvi_knn_pred",
        "output_obs_probability": "scvi_knn_proba",
        "output_obsm_integrated": "X_integrated_scvi",
        "overwrite_existing_key": "true"
      ],
      toState: [ "query_processed": "output" ]
    )

    | scanvi_scarches_annotation.run(
      runIf: { id, state -> state.annotation_methods.contains("scanvi_scarches")},
      fromState: [
        "id": "id",
        "input": "query_processed",
        "modality": "modality",
        "layer": "input_layer",
        "input_var_gene_names": "input_var_gene_names",
        "reference": "reference",
        "reference_obs_target": "reference_obs_label",
        "reference_obs_batch_label": "reference_obs_batch",
        "reference_var_hvg": "reference_var_input",
        "reference_var_gene_names": "reference_var_gene_names",
        "unlabeled_category": "reference_obs_label_unlabeled_category",
        "early_stopping": "early_stopping",
        "early_stopping_monitor": "early_stopping_monitor",
        "early_stopping_patience": "early_stopping_patience",
        "early_stopping_min_delta": "early_stopping_min_delta",
        "max_epochs": "max_epochs",
        "reduce_lr_on_plateau": "reduce_lr_on_plateau",
        "lr_factor": "lr_factor",
        "lr_patience": "lr_patience",
        "leiden_resolution": "leiden_resolution",
        "knn_weights": "knn_weights",
        "knn_n_neighbors": "knn_n_neighbors"
      ],
      args: [
        "input_obs_batch_label": "sample_id",
        "output_obs_predictions": "scanvi_knn_pred",
        "output_obs_probability": "scanvi_knn_proba"
      ],
      toState: [ "query_processed": "output" ]
    )

    | map {id, state ->
      def new_state = state + ["output": state.query_processed]
      [id, new_state]
    }

    | setState(["output", "_meta"])

  emit:
    output_ch
}