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
      toState: [ "input": "output" ]
    )

    | joinStates { ids, states ->
      def newId = "combined"
      // gather keys with unique values across states that should be combined
      def new_state_non_unique_values = [
        input: states.collect{it.input},
        _meta: states.collect{it._meta.join_id}
      ]
      // gather keys from different states
      def all_state_keys = states.inject([].toSet()){ current_keys, state ->
            def new_keys = current_keys + state.keySet()
            return new_keys
        }.minus(["output", "id", "input", "_meta"])
      // Create the new state from the keys, values should be the same across samples
      def new_state = all_state_keys.inject([:]){ old_state, argument_name ->
            argument_values = states.collect{it.get(argument_name)}.unique()
            assert argument_values.size() == 1, "Arguments should be the same across samples. Argument name: $argument_name, \
                                                 argument value: $argument_values"
            // take the unique value from the set (there is only one)
            def argument_value
            argument_values.each { argument_value = it }
            def current_state = old_state + [(argument_name): argument_value]
            return current_state
        }
      def final_state = new_state_non_unique_values + new_state
      [ newId, final_state ]
    }

    | view {"After combining states: $it"}

    // move all processed h5mu files to the same folder
    | move_files_to_directory.run(
      fromState: [
        input: "input",
        output: "output_processed_h5mu"
      ],
      toState: [ "output_processed_h5mu": "output" ]
    )

    | flatMap { id, state ->

        def totalInputs = state.input.size()
        def maxSamplesPerGroup =2
        def numGroups = Math.max(1, Math.ceil(totalInputs / maxSamplesPerGroup) as Integer)
        def baseSamplesPerGroup = totalInputs.intdiv(numGroups)
        def remainder = totalInputs % numGroups

        println "Splitting ${totalInputs} samples into ${numGroups} groups (max ${maxSamplesPerGroup} per group)"
        
        def inputs = []
        for (int i = 0; i < state.input.size(); i++) {
            inputs << [input: state.input[i], _meta: [join_id: state._meta[i]]]
        }

        println "inputs: ${inputs}"

        def sortedInputs = inputs.sort { it._meta.join_id }

        println "sorted inputs: ${sortedInputs}"
        
        def groups = []
        def itemIndex = 0

        (0..<numGroups).each { groupNum ->
            def samplesInGroup = baseSamplesPerGroup + (groupNum < remainder ? 1 : 0)
            def groupItems = sortedInputs[itemIndex..<(itemIndex + samplesInGroup)]
            
            def newId = "combined_${groupNum + 1}_of_${numGroups}"
            def newState = state.clone()  // Copy all the original state

            println "new state: ${newState}"
            
            // Override the input and _meta with the grouped items
            newState.input = groupItems.collect { it.input }
            newState._meta = groupItems[0]._meta
            
            println "Group ${groupNum + 1}: ${samplesInGroup} samples - ${newState._meta}"
            
            groups << [newId, newState]
            itemIndex += samplesInGroup
        }
        
        return groups
    }

    | view {"After grouping: $it"}

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