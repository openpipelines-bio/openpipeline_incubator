nextflow.enable.dsl=2
targetDir = params.rootDir + "/target/nextflow/ingestion_qc"

include { generate_report } from targetDir + "/generate_report/main.nf"

params.resources_test = "s3://openpipelines-bio/openpipeline_incubator/resources_test/qc_sample_data/"

workflow test_no_cellbender {

  resources_test_file = file(params.resources_test)

  output_ch = Channel.fromList([
      [
        id: "sample_one",
        input: resources_test_file.resolve("sample_one.qc.h5mu"),
        run_cellbender: false,
        metadata_obs_keys: ["donor_id", "cell_type", "batch", "condition"],
        output_html: "report.html",
        publish_dir: "test_out"
      ],
      [
        id: "sample_two",
        input: resources_test_file.resolve("sample_two.qc.h5mu"),
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
        assert output[1].output.isFile() : "Output HTML report file should exist"
        "Output: $output"
    }
}


workflow test_with_cellbender {

  resources_test_file = file(params.resources_test)

  output_ch = Channel.fromList([
      [
        id: "sample_one",
        input: resources_test_file.resolve("sample_one.qc.h5mu"),
        run_cellbender: true,
        cellbender_epochs: 1,
        output_html: "report.html",
        publish_dir: "test_out"
      ],
      [
        id: "sample_two",
        input: resources_test_file.resolve("sample_two.qc.h5mu"),
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
        assert output[1].output.isFile() : "Output HTML report file should exist"
        "Output: $output"
    }
}