nextflow.enable.dsl=2

include { atlas_service } from params.rootDir + "/target/nextflow/atlas_service/main.nf"
params.resources_test = params.rootDir + "/resources_test"

workflow test_wf {
  resources_test = file(params.resources_test)

  output_ch = Channel.fromList(
    [
      [
        id: "simple_execution_test",
        input: resources_test.resolve("pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"),
        reference: resources_test.resolve("annotation_test_data/TS_Blood_filtered.h5mu"),
        reference_var_gene_names: "ensemblid",
        reference_layer_lognormalized_counts: "log_normalized",
        reference_obs_batch: "donor_assay",
        reference_obs_label: "cell_type",
        max_epochs: "5",
        annotation_methods: "celltypist;scvi_knn;harmony_knn;scanvi_scarches"
      ]
    ])
    | view {"State at start: $it"}
    | map{ state -> [state.id, state] }
    | atlas_service 
    | view {"After AaaS: $it"}
    | view { output ->
      assert output.size() == 2 : "Outputs should contain two elements; [id, state]"

      // check id
      def id = output[0]
      assert id == "merged" : "Output ID should be `merged`"

      // check output
      def state = output[1]
      assert state instanceof Map : "State should be a map. Found: ${state}"
      assert state.containsKey("output") : "Output should contain key 'output'."
      assert state.output.isFile() : "'output' should be a file."
      assert state.output.toString().endsWith(".h5mu") : "Output file should end with '.h5mu'. Found: ${state.output}"
    
    "Output: $output"
  }
}

workflow test_wf_2 {
  resources_test = file(params.resources_test)

  output_ch = Channel.fromList(
    [
      [
        id: "pbmc",
        input: resources_test.resolve("pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"),
        var_name_mitochondrial_genes: 'mitochondrial',
        rna_min_counts: 2,
        prot_min_counts: 3,
        add_id_to_obs: true,
        add_id_make_observation_keys_unique: true,
        add_id_obs_output: "sample_id",
        reference: resources_test.resolve("annotation_test_data/TS_Blood_filtered.h5mu"),
        reference_var_gene_names: "ensemblid",
        reference_layer_lognormalized_counts: "log_normalized",
        reference_obs_batch: "donor_assay",
        reference_obs_label: "cell_type",
        annotation_methods: "celltypist"
      ],
      [
        id: "pbmc_with_more_params",
        input: resources_test.resolve("pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"),
        rna_min_counts: 2,
        rna_max_counts: 1000000,
        rna_min_genes_per_cell: 1,
        rna_max_genes_per_cell: 1000000,
        rna_min_cells_per_gene: 1,
        rna_min_fraction_mito: 0.0,
        rna_max_fraction_mito: 1.0,
        prot_min_counts: 3,
        prot_max_counts: 1000000,
        prot_min_proteins_per_cell: 1,
        prot_max_proteins_per_cell: 1000000,
        prot_min_cells_per_protein: 1,
        var_name_mitochondrial_genes: 'mitochondrial',
        obs_name_mitochondrial_fraction: 'fraction_mitochondrial',
        add_id_to_obs: true,
        add_id_make_observation_keys_unique: true,
        add_id_obs_output: "sample_id",
        reference: resources_test.resolve("annotation_test_data/TS_Blood_filtered.h5mu"),
        reference_var_gene_names: "ensemblid",
        reference_layer_lognormalized_counts: "log_normalized",
        reference_obs_batch: "donor_assay",
        reference_obs_label: "cell_type",
        annotation_methods: "celltypist"
      ]
    ])
    | view {"State at start: $it"}
    | map { state -> [state.id, state] }
    | atlas_service 
    | view {"After AaaS: $it"}
    | view { output ->
      assert output.size() == 2 : "Outputs should contain two elements; [id, state]"

      // check id
      def id = output[0]
      assert id == "merged" : "Output ID should be `merged`"

      // check output
      def state = output[1]
      assert state instanceof Map : "State should be a map. Found: ${state}"
      assert state.containsKey("output") : "Output should contain key 'output'."
      assert state.output.isFile() : "'output' should be a file."
      assert state.output.toString().endsWith(".h5mu") : "Output file should end with '.h5mu'. Found: ${state.output}"
    
      "Output: $output"
    }
  }

workflow test_wf_3 {
  resources_test = file(params.resources_test)

  output_ch = Channel.fromList(
    [
      [
        id: "celltypist_model",
        input: resources_test.resolve("pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"),
        celltypist_model: resources_test.resolve("annotation_test_data/celltypist_model_Immune_All_Low.pkl"),
        annotation_methods: "celltypist",
        input_var_gene_names: "gene_symbol"
      ]
    ])
    | view {"State at start: $it"}
    | map{ state -> [state.id, state] }
    | atlas_service 
    | view {"After AaaS: $it"}
    | view { output ->
      assert output.size() == 2 : "Outputs should contain two elements; [id, state]"

      // check id
      def id = output[0]
      assert id == "merged" : "Output ID should be `merged`"

      // check output
      def state = output[1]
      assert state instanceof Map : "State should be a map. Found: ${state}"
      assert state.containsKey("output") : "Output should contain key 'output'."
      assert state.output.isFile() : "'output' should be a file."
      assert state.output.toString().endsWith(".h5mu") : "Output file should end with '.h5mu'. Found: ${state.output}"
    
    "Output: $output"
  }
}

workflow test_wf_4 {
  resources_test = file(params.resources_test)

  output_ch = Channel.fromList(
    [
      [
        id: "scgpt",
        input: resources_test.resolve("pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"),
        annotation_methods: "scgpt_annotation",
        input_var_gene_names: "gene_symbol",
        scgpt_model: resources_test.resolve("scgpt/finetuned_model/best_model.pt"),
        scgpt_model_config: resources_test.resolve("scgpt/source/args.json"),
        scgpt_model_vocab: resources_test.resolve("scgpt/source/vocab.json"),
        annotation_methods: "scgpt_annotation"
      ]
    ])
    | view {"State at start: $it"}
    | map{ state -> [state.id, state] }
    | atlas_service 
    | view {"After AaaS: $it"}
    | view { output ->
      assert output.size() == 2 : "Outputs should contain two elements; [id, state]"

      // check id
      def id = output[0]
      assert id == "merged" : "Output ID should be `merged`"

      // check output
      def state = output[1]
      assert state instanceof Map : "State should be a map. Found: ${state}"
      assert state.containsKey("output") : "Output should contain key 'output'."
      assert state.output.isFile() : "'output' should be a file."
      assert state.output.toString().endsWith(".h5mu") : "Output file should end with '.h5mu'. Found: ${state.output}"
    
    "Output: $output"
  }
}
