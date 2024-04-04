# scgog

Repository dedicated to optimize use of common Data Processing, Feature Exraction, and Feature Engineering techniques of Single Cell data mainly for Machine Learning purposes.

- Input: scRNA-seq, ATAC-seq
- User input: choice of Data Analysis techniques and ML models
- Output: Dataset and Model checkpoint

Inspired by BIOINF579 class.


Tutorial:
!wget -O pbmc_tutorial/data/pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix.h5 https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_3k/pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix.h5
!wget -O pbmc_tutorial/data/pbmc_granulocyte_sorted_3k_atac_peak_annotation.tsv.gz https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_3k/pbmc_granulocyte_sorted_3k_atac_peak_annotation.tsv.gz

python -m build
pip install dist/scgog-0.0.1.tar.gz

commands:
preprocess 
Launches CLI to guide the used through preprocessing
visualize
Launches CLI to guide the used through visualizization
ml_benchmarking
Launches CLI to guide the used through ML benchmarking

