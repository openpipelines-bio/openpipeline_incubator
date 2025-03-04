nextflow.enable.dsl=2

include { atlas_service } from params.rootDir + "/target/nextflow/atlas_service/main.nf"
params.resources_test = params.rootDir + "/resources_test"

workflow test_wf {
  // allow changing the resources_test dir
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
        annotation_methods: "celltypist"
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
      assert id.endsWith("_test") : "Output ID should be same as input ID"

      // check output
      def state = output[1]
      assert state instanceof Map : "State should be a map. Found: ${state}"
      assert state.containsKey("output") : "Output should contain key 'output'."
      assert state.output.isFile() : "'output' should be a file."
      assert state.output.toString().endsWith(".h5mu") : "Output file should end with '.h5mu'. Found: ${state.output}"
    
    "Output: $output"
    }
    | toSortedList({a, b -> a[0] <=> b[0]})
    | map { output_list ->
      assert output_list.size() == 2 : "output channel should contain 2 events"
      assert output_list.collect{it[0]} == ["no_leiden_resolutions_test", "simple_execution_test"]
    }
    }