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
      args: [
        obs_output: "sample_id"
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
      args: [
        obs_background_fraction: "cellbender_background_fraction",
        obs_cell_probability: "cellbender_cell_probability",
        obs_droplet_efficiency: "cellbender_droplet_efficiency",
        obs_cell_size: "cellbender_cell_size",
      ],
      toState: ["input": "output"]
    )

    // run qc on each sample
    | qc.run(
      fromState: [
        id: "id",
        input: "input",
        var_gene_names: "var_gene_names"
      ],
      args: [
        var_name_mitochondrial_genes: "mitochondrial",
        var_name_ribosomal_genes: "ribosomal",
        output_obs_num_nonzero_vars: "num_nonzero_vars",
        output_obs_total_counts_vars: "total_counts"
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
        def maxSamplesPerGroup = state.max_samples_per_report
        def numGroups = Math.max(1, Math.ceil(totalInputs / maxSamplesPerGroup) as Integer)
        def baseSamplesPerGroup = totalInputs.intdiv(numGroups)
        def remainder = totalInputs % numGroups

        println "Dividing ${totalInputs} sample(s) over ${numGroups} report(s) (max ${maxSamplesPerGroup} per report)"
        
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

    // Set aside output for QC report instructions
    | map { id, state -> 
      def new_state = state + ["output_reporting_json": "reporting_json.json"]
      [id, new_state]
    }

    // Set report filter settings
    | map { id, state -> 
      def conditionalValues = [
        "cellranger_multi": [
          "min_total_counts": 10,
          "min_num_nonzero_vars": 10
        ],
        "xenium": [
          "min_total_counts": 10, 
          "min_num_nonzero_vars": 1
        ]
      ]
      
      def method = state.ingestion_method
      def additionalParams = conditionalValues[method]
      
      [ id, state + additionalParams ]
    }

    | view {"After setting filters: $it"}

    // generate qc json
    | h5mu_to_qc_json.run(
      fromState: [
        input: "input",
        ingestion_method: "ingestion_method",
        obs_metadata: "obs_metadata",
        min_total_counts: "min_total_counts",
        min_num_nonzero_vars: "min_num_nonzero_vars"
      ],
      args: [
        obs_sample_id: "sample_id",
        obs_total_counts: "total_counts",
        obs_num_nonzero_vars: "num_nonzero_vars",
        obs_fraction_mitochondrial: "fraction_mitochondrial",
        obs_fraction_ribosomal: "fraction_ribosomal",
      ],
      toState: [
        output: "output",
        output_reporting_json: "output_reporting_json"
      ]
    )

    // generate html report
    | generate_html.run(
      fromState: [ 
        input_data: "output",
        input_structure: "output_reporting_json"
      ],
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