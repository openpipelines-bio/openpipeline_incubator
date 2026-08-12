nextflow.enable.dsl=2

include { demultiplex } from params.rootDir + "/target/nextflow/workflows/demultiplexing/demultiplex/main.nf"
include { demultiplex_test } from params.rootDir + "/target/_test/nextflow/test_workflows/demultiplexing/demultiplex_test/main.nf"

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
      assert state.output_sample_qc == null : "'output_sample_qc' should be null when --run_qc is not set."
      assert state.output_multiqc == null : "'output_multiqc' should be null when --run_qc is not set."
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
    | demultiplex_test.run(
      fromState: { id, state ->
        [
          "output_fastq": state.output_fastq,
          "output_fastq_manifest": state.output_fastq_manifest,
        ]
      }
    )
    | toSortedList()
    | map { output_list ->
      assert output_list.size() == 1 : "output channel should contain one event"
      assert output_list[0][0] == "demultiplex_fastq_passthrough_test" : \
        "Output ID should be 'demultiplex_fastq_passthrough_test'"
    }
}


workflow test_wf_qc {

  resources_test = file(params.resources_test)

  output_ch = Channel.fromList([
      [
        id: "demultiplex_qc_test",
        input: resources_test.resolve("demultiplex_fastq"),
        run_qc: true,
        output_fastq: "demultiplex_qc_test.fastq",
        output_fastq_manifest: "demultiplex_qc_test.manifest.csv",
      ],
    ])
    | map { state -> [state.id, state] }
    | demultiplex.run(
      toState: { id, output, state -> output + [og_input: state.input] }
    )
    // Check falco + multiqc actually ran and produced non-empty output
    | view { output ->
      assert output.size() == 2 : "outputs should contain two elements; [id, state]"
      def state = output[1]
      assert state.containsKey("output_sample_qc") : "Output should contain key 'output_sample_qc'."
      assert state.output_sample_qc.isDirectory() : "'output_sample_qc' should be a directory."
      assert state.output_sample_qc.list().size() > 0 : "'output_sample_qc' should not be empty."
      assert state.containsKey("output_multiqc") : "Output should contain key 'output_multiqc'."
      assert state.output_multiqc.isFile() : "'output_multiqc' should be a file."
      assert state.output_multiqc.size() > 0 : "'output_multiqc' should not be empty."
      assert state.output_multiqc.name.endsWith(".html") : "'output_multiqc' should be an HTML report."
      "Output: $output"
    }
    | demultiplex_test.run(
      fromState: { id, state ->
        [
          "output_fastq": state.output_fastq,
          "output_fastq_manifest": state.output_fastq_manifest,
        ]
      }
    )
    | toSortedList()
    | map { output_list ->
      assert output_list.size() == 1 : "output channel should contain one event"
      assert output_list[0][0] == "demultiplex_qc_test" : \
        "Output ID should be 'demultiplex_qc_test'"
    }
}


// Count the FASTQ records in a gzipped FASTQ file (4 lines per record).
def countFastqReads(file) {
  def n_lines = 0
  new java.util.zip.GZIPInputStream(new FileInputStream(file)).eachLine { n_lines++ }
  return n_lines.intdiv(4)
}

workflow test_wf {

  resources_test = file(params.resources_test)

  // Real read counts per sample
  expected_counts = [
    "Sample1": 9920,
    "SampleA": 8560,
    "Sample23": 10111,
    "sampletest": 8925,
  ]

  output_ch = Channel.fromList([
      [
        id: "demultiplex_bcl_convert_test",
        input: resources_test.resolve("novaseq_covidseq_tiny"),
        output_fastq: "demultiplex_bcl_convert_test.fastq",
        output_fastq_manifest: "demultiplex_bcl_convert_test.manifest.csv",
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
    // Check the manifest reflects real bcl-convert output
    | map { id, state ->
      def rows = state.output_fastq_manifest.splitCsv(header: true)
      // RunInfo.xml declares a single genomic read (36bp)
      assert rows.size() == 4 : "Expected 4 manifest rows (one per sample, single-end), found ${rows.size()}: $rows"
      assert rows*.read_type.toSet() == ["R1"].toSet() : \
        "Expected only R1 rows (single-end run), found read types: ${rows*.read_type}"
      assert rows*.sample_id.toSet() == expected_counts.keySet() : \
        "Unexpected sample_ids: ${rows*.sample_id}"
      assert !rows*.sample_id.contains("Undetermined") : \
        "Undetermined reads should be excluded by default (--include_undetermined not set)."
      rows.each { row ->
        def fastq_file = new File(row.path)
        assert fastq_file.length() > 0 : "${row.sample_id} FASTQ file should not be empty."
        def n_reads = countFastqReads(fastq_file)
        assert n_reads == expected_counts[row.sample_id] : \
          "${row.sample_id} should have ${expected_counts[row.sample_id]} reads, found ${n_reads}."
      }
      [id, state]
    }
    | demultiplex_test.run(
      fromState: { id, state ->
        [
          "output_fastq": state.output_fastq,
          "output_fastq_manifest": state.output_fastq_manifest,
        ]
      }
    )
    | toSortedList()
    | map { output_list ->
      assert output_list.size() == 1 : "output channel should contain one event"
      assert output_list[0][0] == "demultiplex_bcl_convert_test" : \
        "Output ID should be 'demultiplex_bcl_convert_test'"
    }
}

workflow test_wf_preset {

  resources_test = file(params.resources_test)

  output_ch = Channel.fromList([
      [
        id: "demultiplex_preset_test",
        input: resources_test.resolve("cellranger_arc_tiny_bcl_atac"),
        preset: "10x_atac",
        output_fastq: "demultiplex_preset_test.fastq",
        output_fastq_manifest: "demultiplex_preset_test.manifest.csv",
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
    // Check index reads (I1/I2) reached the manifest: without --preset 10x_atac,
    // CreateFastqForIndexReads stays off and they never would.
    | map { id, state ->
      def rows = state.output_fastq_manifest.splitCsv(header: true)
      def read_types = rows*.read_type.toSet()
      assert read_types.containsAll(["R1", "R2", "I1", "I2"]) : \
        "Expected R1/R2/I1/I2 read types in the manifest, found: $read_types"
      def i1_rows = rows.findAll { it.read_type == "I1" }
      def i2_rows = rows.findAll { it.read_type == "I2" }
      assert new File(i1_rows[0].path).length() > 0 : "I1 FASTQ file should not be empty."
      assert new File(i2_rows[0].path).length() > 0 : "I2 FASTQ file should not be empty."
      [id, state]
    }
    | demultiplex_test.run(
      fromState: { id, state ->
        [
          "output_fastq": state.output_fastq,
          "output_fastq_manifest": state.output_fastq_manifest,
        ]
      }
    )
    | toSortedList()
    | map { output_list ->
      assert output_list.size() == 1 : "output channel should contain one event"
      assert output_list[0][0] == "demultiplex_preset_test" : \
        "Output ID should be 'demultiplex_preset_test'"
    }
}

workflow test_wf_bases2fastq {

  resources_test = file(params.resources_test)

  output_ch = Channel.fromList([
      [
        id: "demultiplex_bases2fastq_test",
        input: resources_test.resolve("bases2fastq_sim_tiny"),
        output_fastq: "demultiplex_bases2fastq_test.fastq",
        output_fastq_manifest: "demultiplex_bases2fastq_test.manifest.csv",
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
    // Check the manifest reflects real bases2fastq output: 5 real samples, each
    // paired-end (R1/R2) and split across the fixture's 2 lanes.
    | map { id, state ->
      def rows = state.output_fastq_manifest.splitCsv(header: true)
      def expected_samples = (0..4).collect { "sample_${it}" }.toSet()
      assert rows*.sample_id.toSet() == expected_samples : \
        "Unexpected sample_ids: ${rows*.sample_id.unique()}"
      assert !rows*.sample_id.contains("Undetermined") : \
        "Undetermined reads should be excluded by default (--include_undetermined not set)."
      assert rows*.read_type.toSet() == ["R1", "R2"].toSet() : \
        "Expected only R1/R2 rows (paired-end run), found read types: ${rows*.read_type.unique()}"
      assert rows*.lane.toSet() == ["1", "2"].toSet() : \
        "Expected reads split across lanes 1 and 2, found lanes: ${rows*.lane.unique()}"
      expected_samples.each { sample ->
        ["1", "2"].each { lane ->
          ["R1", "R2"].each { read_type ->
            def matches = rows.findAll { it.sample_id == sample && it.lane == lane && it.read_type == read_type }
            assert matches.size() == 1 : \
              "Expected exactly one ${read_type} file for ${sample} in lane ${lane}, found ${matches.size()}."
            def fastq_file = new File(matches[0].path)
            assert fastq_file.length() > 0 : "${sample} lane ${lane} ${read_type} FASTQ file should not be empty."
          }
        }
      }
      [id, state]
    }
    | demultiplex_test.run(
      fromState: { id, state ->
        [
          "output_fastq": state.output_fastq,
          "output_fastq_manifest": state.output_fastq_manifest,
        ]
      }
    )
    | toSortedList()
    | map { output_list ->
      assert output_list.size() == 1 : "output channel should contain one event"
      assert output_list[0][0] == "demultiplex_bases2fastq_test" : \
        "Output ID should be 'demultiplex_bases2fastq_test'"
    }
}
