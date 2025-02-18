# Test required arguments scGPT
cat > params.yaml << HERE
id: scgpt_no_params
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
id: celltypist_overlapping_params
input: /Users/dorienroosen/code/openpipeline/resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu
modality: rna
input_var_gene_names: gene_symbol
annotation_methods: celltypist
HERE

nextflow run . \
-main-script target/nextflow/atlas_service/main.nf \
-params-file params.yaml \
-resume \
-profile docker,no_publish \
-c target/nextflow/atlas_service/nextflow.config \
-c /Users/dorienroosen/code/openpipeline/src/workflows/utils/labels_ci.config

cat > params.yaml << HERE
id: celltypist_no_params
input: /Users/dorienroosen/code/openpipeline/resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu
modality: rna
input_var_gene_names: gene_symbol
reference: /Users/dorienroosen/code/openpipeline/resources_test/annotation_test_data/TS_Blood_filtered.h5mu
annotation_methods: scgpt,celltypist
celltypist_model: /Users/dorienroosen/code/openpipeline/resources_test/annotation_test_data/celltypist_model_Immune_All_Low.pkl
HERE

nextflow run . \
-main-script target/nextflow/atlas_service/main.nf \
-params-file params.yaml \
-resume \
-profile docker,no_publish \
-c target/nextflow/atlas_service/nextflow.config \
-c /Users/dorienroosen/code/openpipeline/src/workflows/utils/labels_ci.config
