ABSOLUTE_INPUT_DATA=$(realpath $par_input_data)
ABSOLUTE_INPUT_STRUCTURE=$(realpath $par_input_structure)
ABSOLUTE_OUTPUT=$(realpath $par_output_qc_report)

cd /opt/incubator_ingestion_qc
mkdir src/data

echo "Compressing input data..."
pnpm run compress_data "$ABSOLUTE_INPUT_DATA" "src/data/dataset.ts"

echo "Compressing report structure..."
pnpm run compress_data "$ABSOLUTE_INPUT_STRUCTURE" "src/data/report_structure.ts"

echo "Generating HTML..."
pnpm run build

echo "Copying HTML to output directory..."
cp dist/index.html "$ABSOLUTE_OUTPUT"
