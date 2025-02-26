echo ">> Generating report"
"$meta_executable" \
    --input "$meta_resources_dir/dataset.json" \
    --output "$meta_resources_dir/index.html" \

echo ">> Checking output"
[ ! -f "$meta_resources_dir/index.html" ] && echo "Error: Output report does not exist." && exit 1

echo ">> Test succesful" && exit 0