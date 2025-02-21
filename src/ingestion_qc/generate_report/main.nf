workflow run_wf {
  take: input_ch
  main:
  output_ch = input_ch

    // store join id
    | map { id, state ->
      [id, state + [_meta: [join_id: id]]]
    }

    // run qc on each sample
    | qc_wf.run(
      fromState: ["id", "input"],
      args: [
        var_gene_names: "gene_symbol",
        var_name_mitochondrial_genes: "mitochondrial_genes",
        var_name_ribosomal_genes: "ribosomal_genes"
      ],
      toState: ["output"]
    )

    // add sample ids to each state
    | add_id.run(
      fromState: [input_id: "id", input: "output"],
      toState: ["output"]
    )

    // combine files into one state
    | joinStates { ids, states ->
      def newId = "combined"
      def newState = [
        input: states.collect{it.output},
        _meta: states[0]._meta,
        output_html: states[0].output_html
      ]
      [newId, newState]
    }

    // generate qc json
    | h5mu_to_qc_json.run(
      fromState: ["input"],
      args: [sample_id_key: "sample_id"],
      toState: [output_qc_json: "output"]
    )

    | generate_html.run(
      fromState: [input: "output_qc_json"],
      toState: [output_html: "output"]
    )

    // emit output
    | setState(["output_qc_json", "output_html", "_meta"])

  emit: output_ch
}