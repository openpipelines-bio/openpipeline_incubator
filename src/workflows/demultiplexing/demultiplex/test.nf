nextflow.enable.dsl=2

include { demultiplex } from params.rootDir + "/target/nextflow/workflows/demultiplexing/demultiplex/main.nf"

params.resources_test = "s3://openpipelines-bio/openpipeline_incubator/resources_test/"

workflow test_wf_fastq_passthrough {

  resources_test = file(params.resources_test)

  output_ch = Channel.fromList([
      [
        id: "demultiplex_fastq_passthrough_test",
        input: resources_test.resolve("demultiplex_fastq"),
        output_fastq: "demultiplex_fastq_passthrough_test.fastq",
        output_fastq_manifest: "demultiplex_fastq_passthrough_test.manifest.csv",
      ],
    ])
    | map { state -> [state.id, state] }
    | demultiplex.run(
      toState: { id, output, state -> output + [og_input: state.input] }
    )
    // Check expected output exists
    | view { output ->
      assert output.size() == 2 : "outputs should contain two elements; [id, state]"
      def state = output[1]
      assert state.containsKey("output_fastq") : "Output should contain key 'output_fastq'."
      assert state.output_fastq.isDirectory() : "'output_fastq' should be a directory."
      assert state.containsKey("output_fastq_manifest") : "Output should contain key 'output_fastq_manifest'."
      assert state.output_fastq_manifest.isFile() : "'output_fastq_manifest' should be a file."
      "Output: $output"
    }
    // Check output FASTQ file match
    | map { id, state ->
      def og_files = state.og_input.listFiles().collect { it.name }.toSet()
      def out_files = state.output_fastq.listFiles().collect { it.name }.toSet()
      assert out_files == og_files : \
        "output_fastq should still contain the original FASTQ files. Expected: $og_files, found: $out_files"

      def rows = state.output_fastq_manifest.splitCsv(header: true)
      assert rows.size() == 5 : "Expected 5 manifest rows, found ${rows.size()}: $rows"
      assert rows*.sample_id.toSet() == ["sample1", "sample2"].toSet() : \
        "Unexpected sample_ids: ${rows*.sample_id}"
      assert rows.findAll { it.sample_id == "sample1" }*.lane.toSet() == ["1", "2"].toSet() : \
        "sample1 should span lanes 1 and 2."
      assert rows.findAll { it.sample_id == "sample2" }*.lane == ["1"] : \
        "sample2 should be in lane 1 only."
      assert rows.findAll { it.sample_id == "sample2" }*.read_type == ["R1"] : \
        "sample2 should be single-end (R1 only)."
      assert !rows*.sample_id.contains("Undetermined") : \
        "Undetermined reads should be excluded by default (--include_undetermined not set)."
      [id, state]
    }
    | toSortedList()
    | map { output_list ->
      assert output_list.size() == 1 : "output channel should contain one event"
      assert output_list[0][0] == "demultiplex_fastq_passthrough_test" : \
        "Output ID should be 'demultiplex_fastq_passthrough_test'"
    }
}
