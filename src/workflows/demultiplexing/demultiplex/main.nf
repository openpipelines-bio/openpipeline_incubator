workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch
    // Unpack the input if it looks like a tar archive
    | untar.run(
      runIf: {id, state -> state.input.toString() ==~ /.*\.(tar|tar\.gz|tgz)$/},
      fromState: ["input": "input"],
      toState: ["input": "output"]
    )
    // Detect whether the input is an Illumina/Element run folder or already FASTQ
    | detect_sequencing_input.run(
      fromState: {id, state ->
        [
          "input": state.input,
          "sample_sheet": state.sample_sheet,
          "demultiplexer": state.demultiplexer,
          "skip_demultiplexing": state.skip_demultiplexing,
        ]
      },
      toState: ["detection_json": "output"]
    )
    // Read the detection results
    | map {id, state ->
      def detection = new groovy.json.JsonSlurper().parseText(state.detection_json.text)
      [id, state + ["detection": detection]]
    }
    // Build a manifest from the FASTQ directory
    | create_fastq_manifest.run(
      fromState: {id, state ->
        [
          "input": state.input,
          "sample_sheet": state.detection.sample_sheet ? file(state.detection.sample_sheet) : null,
          "include_undetermined": state.include_undetermined,
        ]
      },
      toState: ["output_fastq_manifest": "output"]
    )
    | map {id, state -> [id, state + ["output_fastq": state.input]]}
    | setState(
      [
        "output_fastq", "output_fastq_manifest"
      ]
    )

  emit:
  output_ch
}
