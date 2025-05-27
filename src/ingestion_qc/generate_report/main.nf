workflow run_wf {
  take: input_ch
  main:
  qc_ch = input_ch
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
      toState: { id, output, state ->
        def keysToRemove = ["var_gene_names", "var_name_mitochondrial_genes", "var_name_ribosomal_genes", "run_cellbender", "cellbender_epochs"]
        def newState = state.findAll{it.key !in keysToRemove}
        newState + ["input": output.output]
      }
    )

    | joinStates { ids, states ->
      def newId = "qc_data"
      // gather keys with unique values across states that should be combined
      def new_state_non_unique_values = [
        input: states.collect{it.input},
        join_ids: states.collect{it._meta.join_id},
        _meta: [join_id: ids[0]]
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
      def data_state = new_state_non_unique_values + new_state
      [ newId, data_state ]
    }

  processed_files_ch = qc_ch

    // move all processed h5mu files to the same folder
    | move_files_to_directory.run(
      fromState: [
        input: "input",
        output: "output_processed_h5mu"
      ],
      toState: [ "output_processed_h5mu": "output" ]
    )
    | setState(["output_processed_h5mu"])

  report_ch = qc_ch
    // group the processed samples to generate one or multiple reports
    | flatMap { id, state ->

        // calculate number of reports to be generated and number of samples per report
        def totalInputs = state.input.size()
        def maxSamplesPerGroup =2
        def numGroups = Math.max(1, Math.ceil(totalInputs / maxSamplesPerGroup) as Integer)
        def baseSamplesPerGroup = totalInputs.intdiv(numGroups)
        def remainder = totalInputs % numGroups

        println "Splitting ${totalInputs} samples into ${numGroups} groups (max ${maxSamplesPerGroup} per group)"
        
        // sort inputs to make grouping deterministic
        def inputs = []
        for (int i = 0; i < state.input.size(); i++) {
            inputs << [input: state.input[i], _meta: [join_id: state.join_ids[i]]]
        }

        def sortedInputs = inputs.sort { it._meta.join_id }
  
        def groups = []
        def itemIndex = 0

        // create one channel per report
        (0..<numGroups).each { groupNum ->
            def samplesInGroup = baseSamplesPerGroup + (groupNum < remainder ? 1 : 0)
            def groupItems = sortedInputs[itemIndex..<(itemIndex + samplesInGroup)]
            
            def newId = "combined_${groupNum + 1}_of_${numGroups}"
            def newState = state.clone()  // Copy all the original state
            
            // Override the input and _meta with the grouped items
            newState.input = groupItems.collect { it.input }
            newState._meta = groupItems[0]._meta
            
            println "Group ${groupNum + 1}: ${samplesInGroup} samples - ${newState._meta}"
            
            groups << [newId, newState]
            itemIndex += samplesInGroup
        }
        
        return groups
    }

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

    // generate html report
    | generate_html.run(
      fromState: [ input: "output_qc_json" ],
      toState: [
        output_qc_report: "output_qc_report"
      ]
    )

    // collect the reports into a single channel
    | joinStates { ids, states ->
      def newId = "qc_report"
      def report_state = [
        output_qc_report: states.collect{it.output_qc_report},
        _meta: states[0]._meta
      ]
      [ newId, report_state ]
    }
    
  output_ch = report_ch.mix(processed_files_ch)

    | joinStates { ids, states ->

      assert states.size() == 2, "Expected 2 states, but got ${states.size()}"
      assert ids.contains('qc_report'), "Expected one channel to have the id `qc_report`, but got ${ids}"
      assert ids.contains('qc_data'), "Expected one channel to have the id `qc_data`, but got ${ids}"

      def newId = "combined"
      def combined_state = states[0] + states [1]

      [ newId, combined_state ]
    }

  emit: output_ch
}