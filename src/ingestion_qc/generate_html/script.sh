cd /opt/incubator_ingestion_qc
cp $par_input /opt/incubator_ingestion_qc/data/dataset.json

echo "Compressing input data..."
pnpm run compress_data

echo "Generating HTML..."
pnpm run build

echo "Copying HTML to output..."
cp dist/index.html $par_output