# Test required arguments scGPT
cat > params.yaml << HERE
id: scgpt
input: /Users/dorienroosen/code/openpipeline/resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu
modality: rna
input_var_gene_names: gene_symbol
annotation_methods: scgpt_annotation
scgpt_model_config: /Users/dorienroosen/code/openpipeline/resources_test/scgpt/source/args.json
scgpt_model_vocab: /Users/dorienroosen/code/openpipeline/resources_test/scgpt/source/vocab.json
HERE

nextflow run . \
-main-script target/nextflow/atlas_service/main.nf \
-params-file params.yaml \
-resume \
-profile docker,no_publish \
-c target/nextflow/atlas_service/nextflow.config \
-c /Users/dorienroosen/code/openpipeline/src/workflows/utils/labels_ci.config

# Test required arguments CellTypist
cat > params.yaml << HERE
id: celltypist
input: /Users/dorienroosen/code/openpipeline/resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu
modality: rna
input_var_gene_names: gene_symbol
annotation_methods: scgpt,celltypist
HERE

nextflow run . \
-main-script target/nextflow/atlas_service/main.nf \
-params-file params.yaml \
-resume \
-profile docker,no_publish \
-c target/nextflow/atlas_service/nextflow.config \
-c /Users/dorienroosen/code/openpipeline/src/workflows/utils/labels_ci.config

# Test required arguments CellTypist
cat > params.yaml << HERE
id: celltypist
input: /Users/dorienroosen/code/openpipeline/resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu
modality: rna
input_var_gene_names: gene_symbol
annotation_methods: harmony_knn
HERE

nextflow run . \
-main-script target/nextflow/atlas_service/main.nf \
-params-file params.yaml \
-resume \
-profile docker,no_publish \
-c target/nextflow/atlas_service/nextflow.config \
-c /Users/dorienroosen/code/openpipeline/src/workflows/utils/labels_ci.config
