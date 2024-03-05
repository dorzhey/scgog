import pandas as pd
import muon as mu
import scanpy as sc


def preprocess_omics_data(file_path: str) -> mu.MuData:
    """
    Preprocess omics data by reading from an .h5 file and making variable names unique.

    Parameters:
    - file_path (str): Path to the .h5 file.

    Returns:
    - mu.MuData: Preprocessed MuData object.
    """
    mdata = mu.read_10x_h5(file_path)
    mdata.var_names_make_unique()
    
    return mdata


def read_ann_data(file_path: str) -> pd.DataFrame:
    """
    Read annotation data from a .tsv.gz file.

    Parameters:
    - file_path (str): Path to the .tsv.gz file.

    Returns:
    - pd.DataFrame: Annotation data.
    """
    ann_data = pd.read_csv(file_path, sep='\t', compression='gzip')
    return ann_data


def quality_control(mdata: mu.MuData) -> mu.MuData:
    """
    Apply quality control measures to multimodal (RNA and ATAC) data.

    Parameters:
    - mdata (mu.MuData): MuData object containing RNA and ATAC data.

    Returns:
    - mu.MuData: MuData object after applying quality control.
    """
    # Compute QC metrics for RNA
    mdata['rna'].var['mt'] = mdata['rna'].var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(mdata['rna'], qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    
    # Compute QC metrics for ATAC
    sc.pp.calculate_qc_metrics(mdata['atac'], percent_top=None, log1p=False, inplace=True)
    
    mu.pp.filter_obs(mdata['rna'], 'n_genes_by_counts', lambda x: (x >= 200) & (x < 5000))
    mu.pp.filter_obs(mdata['rna'], 'total_counts', lambda x: x < 15000)
    mu.pp.filter_obs(mdata['rna'], 'pct_counts_mt', lambda x: x < 20)
    # Filter ATAC based on quality metrics
    mu.pp.filter_obs(mdata['atac'], 'n_genes_by_counts', lambda x: (x >= 2000) & (x <= 15000))
    mu.pp.filter_obs(mdata['atac'], 'total_counts', lambda x: (x >= 4000) & (x <= 40000))
    # Intersect observations to keep only cells present in both modalities
    mu.pp.intersect_obs(mdata)
    
    return mdata


def merge_data(mdata: mu.MuData) -> mu.MuData:
    """
    Merge scRNA-seq and scATAC-seq data into a single MuData object.

    Parameters:
    - mdata (mu.MuData): MuData object containing RNA and ATAC data.

    Returns:
    - mu.MuData: A MuData object containing the merged datasets.
    """
    return mu.MuData({'rna': mdata['rna'], 'atac': mdata['atac']})
