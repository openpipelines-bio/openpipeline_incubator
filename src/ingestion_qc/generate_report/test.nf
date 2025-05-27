nextflow.enable.dsl=2
targetDir = params.rootDir + "/target/nextflow/ingestion_qc"

include { generate_report } from targetDir + "/generate_report/main.nf"

params.resources_test = "s3://openpipelines-bio/openpipeline_incubator/resources_test/"

workflow test_no_cellbender {

  resources_test_file = file(params.resources_test)

  output_ch = Channel.fromList([
      [
        id: "sample_1",
        input: resources_test_file.resolve("qc_sample_data/sample_one.qc.h5mu"),
        run_cellbender: false,
        metadata_obs_keys: ["donor_id", "cell_type", "batch", "condition"],
        output_html: "report.html",
        publish_dir: "test_out"
      ],
      [
        id: "sample_2",
        input: resources_test_file.resolve("qc_sample_data/sample_two.qc.h5mu"),
        run_cellbender: false,
        metadata_obs_keys: ["donor_id", "cell_type", "batch", "condition"],
        output_html: "report.html",
        publish_dir: "test_out"
      ]
    ])

    | map{ state -> [state.id, state] }
    | generate_report

    | view { output ->
        assert output.size() == 2 : "Outputs should contain two elements; [id, state]"
        def id = output[0]
        def state = output [1]
        assert id == "combined": "Output ID should be `combined`"
        assert state instanceof Map : "State should be a map. Found: ${state}"
        assert state.containsKey("output_qc_report"): "Output should contain key `output_qc_report`"
        assert state.containsKey("output_processed_h5mu"): "Output should contain key `output_processed_h5mu`"
        assert state.output_qc_report.size() == 1 : "Expected exactly one output HTML file to be generated"
        assert state.output_qc_report.every { it.isFile()} : "All output HTML report file should exist"
        assert state.output_processed_h5mu.isDirectory() : "Output directory should exist"
        def files = state.output_processed_h5mu.listFiles().findAll { it.isFile() }
        assert files.size() == 2 : "Output directory should contain exactly 2 files, but found ${files.size()} files"
        "Output: $output"
    }
}

workflow test_with_cellbender {

  resources_test_file = file(params.resources_test)

  output_ch = Channel.fromList([
      [
        id: "sample_one",
        input: resources_test_file.resolve("qc_sample_data/sample_one.qc.h5mu"),
        run_cellbender: true,
        cellbender_epochs: 1,
        output_html: "report.html",
        publish_dir: "test_out"
      ],
      [
        id: "sample_two",
        input: resources_test_file.resolve("qc_sample_data/sample_two.qc.h5mu"),
        run_cellbender: true,
        cellbender_epochs: 1,
        output_html: "report.html",
        publish_dir: "test_out"
      ]
    ])

    | map{ state -> [state.id, state] }
    | generate_report

    | view { output ->
        assert output.size() == 2 : "Outputs should contain two elements; [id, state]"
        def id = output[0]
        def state = output [1]
        assert id == "combined": "Output ID should be `combined`"
        assert state instanceof Map : "State should be a map. Found: ${state}"
        assert state.containsKey("output_qc_report"): "Output should contain key `output_qc_report`"
        assert state.containsKey("output_processed_h5mu"): "Output should contain key `output_processed_h5mu`"
        assert state.output_qc_report.size() == 1 : "Expected exactly one output HTML file to be generated"
        assert state.output_qc_report.every { it.isFile()} : "All output HTML report file should exist"
        assert state.output_processed_h5mu.isDirectory() : "Output directory should exist"
        def files = state.output_processed_h5mu.listFiles().findAll { it.isFile() }
        assert files.size() == 2 : "Output directory should contain exactly 2 files, but found ${files.size()} files"
        "Output: $output"
    }
}

workflow test_multiple_reports {

  resources_test_file = file(params.resources_test)

  output_ch = Channel.fromList([
      [
        id: "sample_1",
        input: resources_test_file.resolve("qc_sample_data/sample_one.qc.h5mu"),
        run_cellbender: false,
        metadata_obs_keys: ["donor_id", "cell_type", "batch", "condition"],
        output_html: "report.html",
        max_samples_per_report: 2,
        publish_dir: "test_out"
      ],
      [
        id: "sample_2",
        input: resources_test_file.resolve("qc_sample_data/sample_two.qc.h5mu"),
        run_cellbender: false,
        metadata_obs_keys: ["donor_id", "cell_type", "batch", "condition"],
        output_html: "report.html",
        max_samples_per_report: 2,
        publish_dir: "test_out"
      ],
      [
        id: "sample_3",
        input: resources_test_file.resolve("qc_sample_data/sample_one.qc.h5mu"),
        run_cellbender: false,
        metadata_obs_keys: ["donor_id", "cell_type", "batch", "condition"],
        output_html: "report.html",
        max_samples_per_report: 2,
        publish_dir: "test_out"
      ],
      [
        id: "sample_4",
        input: resources_test_file.resolve("qc_sample_data/sample_two.qc.h5mu"),
        run_cellbender: false,
        metadata_obs_keys: ["donor_id", "cell_type", "batch", "condition"],
        output_html: "report.html",
        max_samples_per_report: 2,
        publish_dir: "test_out"
      ],
      [
        id: "sample_5",
        input: resources_test_file.resolve("qc_sample_data/sample_one.qc.h5mu"),
        run_cellbender: false,
        metadata_obs_keys: ["donor_id", "cell_type", "batch", "condition"],
        output_html: "report.html",
        max_samples_per_report: 2,
        publish_dir: "test_out"
      ]
    ])

    | map{ state -> [state.id, state] }
    | generate_report

    | view { output ->
        assert output.size() == 2 : "Outputs should contain two elements; [id, state]"
        def id = output[0]
        def state = output [1]
        assert id == "combined": "Output ID should be `combined`"
        assert state instanceof Map : "State should be a map. Found: ${state}"
        assert state.containsKey("output_qc_report"): "Output should contain key `output_qc_report`"
        assert state.containsKey("output_processed_h5mu"): "Output should contain key `output_processed_h5mu`"
        assert state.output_qc_report.size() == 3 : "Expected exactly one output HTML file to be generated"
        assert state.output_qc_report.every { it.isFile()} : "All output HTML report file should exist"
        assert state.output_processed_h5mu.isDirectory() : "Output directory should exist"
        def files = state.output_processed_h5mu.listFiles().findAll { it.isFile() }
        assert files.size() == 5 : "Output directory should contain exactly 5 files, but found ${files.size()} files"
        "Output: $output"
    }
}