workflow run_wf {
  take: input_ch
  main:
  h5mu_ch = input_ch

    // run cellbender
    | cellbender.run(
      runIf: {id, state -> state.run_cellbender},
      fromState: [
        id: "id",
        input: "input",
        epochs: "cellbender_epochs",
      ],
      toState: { id, output, state -> 
        state + ["input": output.output,
                 "metadata_obs_keys": state.metadata_obs_keys]
      }
    )

    // run qc on each sample
    | qc_wf.run(
      fromState: [
        id: "id",
        input: "input",
        var_gene_names: "var_gene_names",
        var_name_mitochondrial_genes: "var_name_mitochondrial_genes",
        var_name_ribosomal_genes: "var_name_ribosomal_genes"
      ],
      toState: ["output_processed_h5mu": "output", 
                "metadata_obs_keys": "metadata_obs_keys"]
    )

    qc_ch = h5mu_ch

    // add sample ids to each state
    | add_id.run(
      fromState: [input_id: "id", input: "output_processed_h5mu"],
      toState: ["output": "output",
                "output_processed_h5mu": "output_processed_h5mu"]
    )

    // combine files into one state
    | joinStates { ids, states ->
      def newId = "combined"
      def newState = [
        input: states.collect{it.output},
        _meta: states[0]._meta,
        output_html: states[0].output_html,
        metadata_obs_keys: states[0].metadata_obs_keys,
        output_processed_h5mu: states.collect{it.output_processed_h5mu}
      ]
      [newId, newState]
    }

    // | view

    // generate qc json
    | h5mu_to_qc_json.run(
      fromState: ["input"],
      args: [
        sample_id_key: "sample_id",
        metadata_obs_keys: "metadata_obs_keys",
      ],
      toState: [output_qc_json: "output",
                output_processed_h5mu: "output_processed_h5mu"]
    )

    | generate_html.run(
      fromState: [input: "output_qc_json"],
      toState: [output_qc_report: "output_qc_report",
                output_processed_h5mu: "output_processed_h5mu"]
    )

    | view

    output_ch = h5mu_ch.combine(qc_ch)
      | map {sample_id, sample_state, combined_id, combined_state ->
          // Create new state by adding the QC report to the sample state
          def new_state = sample_state + ["output_qc_report": combined_state.output_qc_report]
          // Return the tuple with sample_id and the new state
          return [sample_id, new_state]
      }
      // view for debugging
      | view
      // emit output
      | setState(["output_qc_report", "output_processed_h5mu"])

  emit: output_ch
}