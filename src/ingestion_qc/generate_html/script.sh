ABSOLUTE_INPUT=$(realpath $par_input)
ABSOLUTE_OUTPUT=$(realpath $par_output)

cd /opt/incubator_ingestion_qc
mkdir src/data

echo "Absolute input path: $ABSOLUTE_INPUT"
echo "Absolute output path: $ABSOLUTE_OUTPUT"

echo "Compressing input data..."
pnpm run compress_data "$ABSOLUTE_INPUT" "src/data/dataset.ts"

echo "Generating HTML..."
pnpm run build

cp dist/index.html "$ABSOLUTE_OUTPUT"