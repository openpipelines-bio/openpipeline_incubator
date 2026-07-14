nextflow.enable.dsl=2

include { process_single_sample } from params.rootDir + "/target/nextflow/workflows/processing/process_single_sample/main.nf"
include { process_single_sample_test } from params.rootDir + "/target/_test/nextflow/test_workflows/processing/process_single_sample_test/main.nf"

params.resources_test = "s3://openpipelines-bio/openpipeline_incubator/resources_test/"

workflow test_wf {

  resources_test = file(params.resources_test)

  output_ch = Channel.fromList([
      [
        id: "qc_filter_test",
        input: resources_test.resolve("pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"),
        rna_min_total_counts: 100,
        rna_max_log1p_total_counts_quantile: 0.99,
        max_pct_counts_mitochondrial: 50,
        var_gene_names: "gene_symbol",
        prot_min_total_counts: 5,
        prot_max_log1p_total_counts_quantile: 0.99,
        filter_genes_min_cells: 3,
        output_processed_h5mu: "qc_filter_test.processed_h5mu",
        cell_count_report: "qc_filter_test.cell_counts.tsv"
      ],
    ])
    | map{ state -> [state.id, state] }
    | process_single_sample.run(
      // Capture the original per-sample id before the workflow collapses every sample into the
      // single "combined" event, so the stamped sample-id column and report grouping can be
      // checked against it.
      toState: { id, output, state -> output + [og_input: state.input, og_id: state.id] }
    )
    // Resolve the single per-sample h5mu inside the published output directory.
    | map { id, state ->
      def h5mu_files = state.output_processed_h5mu.listFiles().findAll { it.name.endsWith(".h5mu") }
      assert h5mu_files.size() == 1 : "Expected one processed h5mu in ${state.output_processed_h5mu}, found: ${h5mu_files}"
      [id, state + [processed_file: h5mu_files[0]]]
    }
    // Smoke test: the channel event and the emitted output files have the expected shape.
    | view { output ->
      assert output.size() == 2 : "outputs should contain two elements; [id, state]"
      assert output[1].processed_file.toString().endsWith(".h5mu") : "Processed output should be a h5mu file. Found: ${output[1].processed_file}"
      assert output[1].cell_count_report.toString().endsWith(".tsv") : "Cell count report should be a tsv file. Found: ${output[1].cell_count_report}"
      "Output: $output"
    }
    // Data check: the workflow subset the data, intersected the modalities and produced a
    // cell-count report consistent with the final output.
    | process_single_sample_test.run(
      fromState: { id, state ->
        [
          "input": state.processed_file,
          "og_input": state.og_input,
          "cell_count_report": state.cell_count_report,
          "rna_modality": "rna",
          "prot_modality": "prot",
          "expected_sample_id": state.og_id,
          "expect_scrublet": true,
        ]
      }
    )
    | toSortedList()
    | map { output_list ->
      assert output_list.size() == 1 : "output channel should contain one event"
      assert output_list[0][0] == "combined" : "Output ID should be 'combined' (all samples collapsed into one event)"
    }
}

workflow test_wf_skip_scrublet {

  resources_test = file(params.resources_test)

  output_ch = Channel.fromList([
      [
        id: "qc_filter_skip_scrublet_test",
        input: resources_test.resolve("pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"),
        rna_min_total_counts: 100,
        rna_max_log1p_total_counts_quantile: 0.99,
        max_pct_counts_mitochondrial: 50,
        var_gene_names: "gene_symbol",
        prot_min_total_counts: 5,
        prot_max_log1p_total_counts_quantile: 0.99,
        filter_genes_min_cells: 3,
        skip_scrublet: true,
        output_processed_h5mu: "qc_filter_skip_scrublet_test.processed_h5mu",
        cell_count_report: "qc_filter_skip_scrublet_test.cell_counts.tsv"
      ],
    ])
    | map{ state -> [state.id, state] }
    | process_single_sample.run(
      // Capture the original per-sample id before the workflow collapses every sample into the
      // single "combined" event, so the stamped sample-id column and report grouping can be
      // checked against it.
      toState: { id, output, state -> output + [og_input: state.input, og_id: state.id] }
    )
    // Resolve the single per-sample h5mu inside the published output directory.
    | map { id, state ->
      def h5mu_files = state.output_processed_h5mu.listFiles().findAll { it.name.endsWith(".h5mu") }
      assert h5mu_files.size() == 1 : "Expected one processed h5mu in ${state.output_processed_h5mu}, found: ${h5mu_files}"
      [id, state + [processed_file: h5mu_files[0]]]
    }
    // Smoke test: the channel event and the emitted output files have the expected shape.
    | view { output ->
      assert output.size() == 2 : "outputs should contain two elements; [id, state]"
      assert output[1].processed_file.toString().endsWith(".h5mu") : "Processed output should be a h5mu file. Found: ${output[1].processed_file}"
      assert output[1].cell_count_report.toString().endsWith(".tsv") : "Cell count report should be a tsv file. Found: ${output[1].cell_count_report}"
      "Output: $output"
    }
    // Data check: the workflow subset the data and produced a consistent cell-count report.
    // Scrublet was skipped, so no scrublet stage is expected in the report.
    | process_single_sample_test.run(
      fromState: { id, state ->
        [
          "input": state.processed_file,
          "og_input": state.og_input,
          "cell_count_report": state.cell_count_report,
          "rna_modality": "rna",
          "prot_modality": "prot",
          "expected_sample_id": state.og_id,
          "expect_scrublet": false,
        ]
      }
    )
    | toSortedList()
    | map { output_list ->
      assert output_list.size() == 1 : "output channel should contain one event"
      assert output_list[0][0] == "combined" : "Output ID should be 'combined' (all samples collapsed into one event)"
    }
}
