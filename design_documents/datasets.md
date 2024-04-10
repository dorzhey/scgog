# scRNA-seq and ATAC-seq
in .h5 format of 10x Genomics

# Peak annotations
in .tsv format

# Real dataset for answering a biological question using the tool

PBMC from a healthy donor - granulocytes removed through cell sorting (3k) [link](https://support.10xgenomics.com/single-cell-multiome-atac-gex/datasets/1.0.0/pbmc_granulocyte_sorted_3k)

Single Cell Multiome ATAC + Gene Exp. Dataset by Cell Ranger ARC 1.0.0

The PBMC (Peripheral Blood Mononuclear Cells) dataset from a healthy donor, especially one that has been processed to remove granulocytes and then used for single-cell Multiome ATAC + Gene Expression analysis by 10x Genomics, is popular in tutorials for several reasons:

* Comprehensive Data Type: Offers both ATAC-seq and gene expression data from the same cells, showcasing integrative genomic analyses.
* Standardized and High-Quality: Generated using well-documented protocols by 10x Genomics, ensuring consistent data quality ideal for educational examples.
* Broad Research Relevance: Relevant to a wide range of fields like immunology and cancer research, making it an engaging and versatile teaching dataset.
* Open Access: Licensed under Creative Commons, allowing free and widespread use in tutorials and educational material.


## Research Question:

How can we identify and characterize the distinct cell types present in a healthy human PBMC sample using integrative analysis of single-cell ATAC-seq and gene expression data?

## Expected Results:

* Cluster Identification: The preprocessing, normalization, and visualization steps are expected to reveal distinct clusters in the data, each corresponding to different cell types or states within the PBMC sample. These clusters should be distinguishable based on their gene expression and chromatin accessibility profiles.
* Machine Learning Validation: Using clustering results as labels for machine learning models should validate the initial cluster identification, demonstrating that the data contains sufficient information to accurately classify cells into their respective types and states without prior labeling.
* Expected Answer: The integrative analysis of single-cell ATAC-seq and gene expression data enables effective identification and characterization of cell types within healthy human PBMC samples, providing insights into the cellular composition and functional states of immune cells in a healthy individual.
