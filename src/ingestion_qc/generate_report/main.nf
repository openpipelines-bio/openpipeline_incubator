workflow run_wf {
  take: input_ch
  main:
  output_ch = input_ch

    // store join id
    | map { id, state ->
      [id, state + [_meta: [join_id: id]]]
    }

    // add sample ids to each state
    | add_id.run(
      fromState: [
        input_id: "id", 
        input: "input"
      ],
      toState: [ "input": "output" ]
    )

    // run cellbender
    | cellbender.run(
      runIf: {id, state -> state.run_cellbender},
      fromState: [
        id: "id",
        input: "input",
        epochs: "cellbender_epochs",
      ],
      toState: ["input": "output"]
    )

    // run qc on each sample
    | qc.run(
      fromState: [
        id: "id",
        input: "input",
        var_gene_names: "var_gene_names",
        var_name_mitochondrial_genes: "var_name_mitochondrial_genes",
        var_name_ribosomal_genes: "var_name_ribosomal_genes"
      ],
      toState: [ "output": "output" ]
    )

    // combine files into one state
    | joinStates { ids, states ->
      def newId = "combined"
      def newState = [
        input: states.collect{it.output},
        _meta: states[0]._meta,
        output_html: states[0].output_html,
      ]
      [ newId, newState ]
    }

    // move all processed h5mu files to the same folder
    | move_files_to_directory.run(
      fromState: [
        input: "input",
        output: "output_processed_h5mu"
      ],
      toState: [ "output_processed_h5mu": "output" ]
    )

    // generate qc json
    | h5mu_to_qc_json.run(
      fromState: ["input"],
      args: [
        sample_id_key: "sample_id",
        metadata_obs_keys: "metadata_obs_keys",
      ],
      toState: [
        output_qc_json: "output"
      ]
    )

    | generate_html.run(
      fromState: [ input: "output_qc_json" ],
      toState: [
        output_qc_report: "output_qc_report"
      ]
    )

    | setState([ "_meta", "output_qc_report", "output_processed_h5mu" ])

  emit: output_ch
}