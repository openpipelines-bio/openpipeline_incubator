workflow run_wf {
  take:
  input_ch

  main:
  ch = input_ch
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
    // Read the detection results.
    | map {id, state ->
      def detection = new groovy.json.JsonSlurper().parseText(state.detection_json.text)
      def logsRequested = state.output_demultiplexer_logs \
        && !state.output_demultiplexer_logs.toString().endsWith('$id.$key.output_demultiplexer_logs')

      [id, state + [
        "detection": detection,
        // Temporarily set optional outputs to null, they will be set later
        "output_demultiplexer_logs": logsRequested ? state.output_demultiplexer_logs : null,
        "output_sample_qc": null,
        "output_multiqc": null,
      ]]
    }

  // Split into two literal channels rather than using runEach's filter, which
  // silently drops non-matching events (chPassthrough would stay Channel.empty()
  // unless runIf is set).
  chDemux = ch
    | filter {id, state -> state.detection.demultiplexer != null}
    // Track the sample sheet to actually feed bcl-convert; patched below when
    // a preset applies, left as the detected sheet otherwise.
    | map {id, state ->
      def sheet = state.detection.sample_sheet ? file(state.detection.sample_sheet) : null
      [id, state + ["bclconvert_sample_sheet": sheet]]
    }
    // Patch the sample sheet for a 10x chemistry. Only for bcl-convert: there is
    // no GEX preset (needs none of this) and bases2fastq presets don't exist.
    | apply_bclconvert_preset.run(
      runIf: {id, state -> state.preset && state.detection.demultiplexer == "bclconvert"},
      fromState: {id, state ->
        [
          "sample_sheet": state.bclconvert_sample_sheet,
          "detection_json": state.detection_json,
          "preset": state.preset,
        ]
      },
      toState: ["bclconvert_sample_sheet": "output"]
    )
    // Two separate .run() calls rather than runEach: detect_sequencing_input
    // guarantees demultiplexer is exactly one of bclconvert/bases2fastq/null
    // here, so whichever runIf doesn't match passes its event through
    // untouched, and runEach doesn't support the per-call args/key used below.
    | bcl_convert.run(
      runIf: {id, state -> state.detection.demultiplexer == "bclconvert"},
      fromState: {id, state ->
        [
          "bcl_input_directory": state.input,
          "sample_sheet": state.bclconvert_sample_sheet,
          "no_lane_splitting": state.bclconvert_no_lane_splitting,
          "tiles": state.bclconvert_tiles,
          "exclude_tiles": state.bclconvert_exclude_tiles,
          "first_tile_only": state.bclconvert_first_tile_only,
          "strict_mode": state.bclconvert_strict_mode,
          "fastq_gzip_compression_level": state.bclconvert_fastq_gzip_compression_level,
          "num_unknown_barcodes_reported": state.bclconvert_num_unknown_barcodes_reported,
          // Explicit nulls (rather than omitting the keys) so Viash's
          // auto-generated default output filename doesn't kick in: an
          // unrequested optional output must not get produced at all.
          "reports": null,
          "logs": state.output_demultiplexer_logs,
        ]
      },
      toState: {id, output, state ->
        def newState = state + ["input": output.output_directory]
        state.output_demultiplexer_logs ? newState + ["output_demultiplexer_logs": output.logs] : newState
      }
    )
    | bases2fastq.run(
      runIf: {id, state -> state.detection.demultiplexer == "bases2fastq"},
      fromState: {id, state ->
        [
          "analysis_directory": state.input,
          "legacy_fastq": state.bases2fastq_legacy_fastq,
          "split_lanes": state.bases2fastq_split_lanes,
          "detect_adapters": state.bases2fastq_detect_adapters,
          "r1_cycles": state.bases2fastq_r1_cycles,
          "r2_cycles": state.bases2fastq_r2_cycles,
          "i1_cycles": state.bases2fastq_i1_cycles,
          "i2_cycles": state.bases2fastq_i2_cycles,
          "report": null,
          "logs": state.output_demultiplexer_logs,
        ]
      },
      toState: {id, output, state ->
        def newState = state + ["input": output.output_directory]
        state.output_demultiplexer_logs ? newState + ["output_demultiplexer_logs": output.logs] : newState
      }
    )

  chPassthrough = ch
    | filter {id, state -> state.detection.demultiplexer == null}

  chMixed = chDemux.mix(chPassthrough)

  output_ch = chMixed
    // Build a manifest from the FASTQ directory: the demultiplexed output for
    // chDemux events, the original directory as-is for chPassthrough events.
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
    // Run falco once over every FASTQ path in the manifest
    | falco.run(
      runIf: {id, state -> state.run_qc},
      fromState: {id, state ->
        def fastq_paths = state.output_fastq_manifest.splitCsv(header: true)*.path.collect{file(it)}.unique()
        [
          "input": fastq_paths,
          "allow_empty_input": true,
          // Disable falco reports, they don't work with multiple input files
          "data_filename": null,
          "report_filename": null,
          "summary_filename": null,
        ]
      },
      toState: ["output_sample_qc": "outdir"]
    )
    | multiqc.run(
      runIf: {id, state -> state.run_qc},
      fromState: {id, state ->
        [
          "input": [state.output_sample_qc],
          // Set the falco file paths
          "cl_config": "sp: {fastqc/data: {fn: '*_fastqc_data.txt'}}",
        ]
      },
      toState: ["output_multiqc": "output_report"]
    )
    | setState(
      [
        "output_fastq", "output_fastq_manifest", "output_demultiplexer_logs",
        "output_sample_qc", "output_multiqc"
      ]
    )

  emit:
  output_ch
}
