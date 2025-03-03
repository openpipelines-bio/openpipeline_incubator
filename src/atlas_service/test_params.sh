# Test required arguments scGPT
cat > params.yaml << HERE
id: scgpt
input: /home/dorienroosen/openpipeline/resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu
modality: rna
input_var_gene_names: gene_symbol
annotation_methods: scgpt_annotation
scgpt_model_config: /home/dorienroosen/openpipeline/resources_test/scgpt/source/args.json
scgpt_model_vocab: /home/dorienroosen/openpipeline/resources_test/scgpt/source/vocab.json
HERE

nextflow run . \
-main-script target/nextflow/atlas_service/main.nf \
-params-file params.yaml \
-resume \
-profile docker,no_publish \
-c target/nextflow/atlas_service/nextflow.config \
-c /home/dorienroosen/openpipeline/src/workflows/utils/labels_ci.config

# Test required arguments CellTypist
cat > params.yaml << HERE
id: celltypist_1
input: /home/dorienroosen/openpipeline/resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu
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
-c /home/dorienroosen/openpipeline/src/workflows/utils/labels_ci.config

cat > params.yaml << HERE
id: celltypist_2
input: /home/dorienroosen/openpipeline/resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu
celltypist_model: /home/dorienroosen/openpipeline/resources_test/annotation_test_data/celltypist_model_Immune_All_Low.pkl
reference: /home/dorienroosen/openpipeline/resources_test/annotation_test_data/TS_Blood_filtered.h5mu
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
-c /home/dorienroosen/openpipeline/src/workflows/utils/labels_ci.config


# Test required arguments Harmony
cat > params.yaml << HERE
id: harmony
input: /home/dorienroosen/openpipeline/resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu
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
-c /home/dorienroosen/openpipeline/src/workflows/utils/labels_ci.config

# Test required arguments SCVI
cat > params.yaml << HERE
id: scvi
input: /home/dorienroosen/openpipeline/resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu
modality: rna
input_var_gene_names: gene_symbol
annotation_methods: scvi_knn
HERE

nextflow run . \
-main-script target/nextflow/atlas_service/main.nf \
-params-file params.yaml \
-resume \
-profile docker,no_publish \
-c target/nextflow/atlas_service/nextflow.config \
-c /home/dorienroosen/openpipeline/src/workflows/utils/labels_ci.config