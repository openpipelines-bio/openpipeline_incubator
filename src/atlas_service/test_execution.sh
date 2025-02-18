cat > params.yaml << HERE
id: harmony
input: /Users/dorienroosen/code/openpipeline/resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu
modality: rna
reference: /Users/dorienroosen/code/openpipeline/resources_test/annotation_test_data/TS_Blood_filtered.h5mu
reference_var_gene_names: ensemblid
reference_obs_batch: donor_assay
reference_obs_label: cell_type
celltypist_model: /Users/dorienroosen/code/openpipeline/resources_test/annotation_test_data/celltypist_model_Immune_All_Low.pkl
annotation_methods: harmony_knn
HERE

nextflow run . \
-main-script target/nextflow/atlas_service/main.nf \
-params-file params.yaml \
-resume \
-profile docker,no_publish \
-c target/nextflow/atlas_service/nextflow.config \
-c /Users/dorienroosen/code/openpipeline/src/workflows/utils/labels_ci.config

cat > params.yaml << HERE
id: celltypist
input: /Users/dorienroosen/code/openpipeline/resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu
modality: rna
input_var_gene_names: gene_symbol
annotation_methods: celltypist
publish_dir: output
HERE

nextflow run . \
-main-script target/nextflow/atlas_service/main.nf \
-params-file params.yaml \
-resume \
-profile docker,no_publish \
-c target/nextflow/atlas_service/nextflow.config \
-c /Users/dorienroosen/code/openpipeline/src/workflows/utils/labels_ci.config

cat > params.yaml << HERE
id: scgpt
input: /Users/dorienroosen/code/openpipeline/resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu
modality: rna
annotation_methods: scgpt_annotation
input_var_gene_names: gene_symbol
scgpt_model: /Users/dorienroosen/code/openpipeline/resources_test/scgpt/finetuned_model/best_model.pt
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
