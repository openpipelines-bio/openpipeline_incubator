echo ">> Generating report"
"$meta_executable" \
    --input_data "$meta_resources_dir/sc_dataset.json" \
    --input_structure "$meta_resources_dir/sc_report_structure.json" \
    --output_qc_report "index.html" \

echo ">> Checking output"
[ ! -f "index.html" ] && echo "Error: Output report does not exist." && exit 1

echo ">> Test succesful" && exit 0