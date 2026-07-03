nextflow.enable.dsl=2

include { qc_filter_prot } from params.rootDir + "/target/nextflow/workflows/qc/qc_filter_prot/main.nf"
include { qc_filter_test } from params.rootDir + "/target/_test/nextflow/test_workflows/qc/qc_filter_test/main.nf"

params.resources_test = params.rootDir + "/resources_test"

workflow test_wf {

  resources_test = file(params.resources_test)

  output_ch = Channel.fromList([
      [
        id: "qc_filter_prot_test",
        input: resources_test.resolve("pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"),
        min_total_counts: 5,
        max_log1p_total_counts_quantile: 0.99,
        output: "qc_filter_prot_test.output.h5mu",
      ],
    ])
    | map { state -> [state.id, state] }
    | qc_filter_prot.run(
      toState: { id, output, state -> output + [og_input: state.input] }
    )
    // Smoke test: the channel event and the emitted output file have the expected shape.
    | view { output ->
      assert output.size() == 2 : "outputs should contain two elements; [id, state]"
      def state = output[1]
      assert state.containsKey("output") : "Output should contain key 'output'."
      assert state.output.isFile() : "'output' should be a file."
      assert state.output.toString().endsWith(".h5mu") : "Output should be a h5mu file. Found: ${state.output}"
      "Output: $output"
    }
    // Data check: filters ran in flag-only mode (no cells/genes dropped) and the protein
    // keep-flag columns were written to .obs.
    | qc_filter_test.run(
      fromState: { id, state ->
        [
          "input": state.output,
          "og_input": state.og_input,
          "modality": "prot",
          "keep_flag_columns": ["filter_counts_prot", "filter_quantile_prot"],
        ]
      }
    )
    | toSortedList()
    | map { output_list ->
      assert output_list.size() == 1 : "output channel should contain one event"
      assert output_list[0][0] == "qc_filter_prot_test" : "Output ID should be 'qc_filter_prot_test'"
    }
}
