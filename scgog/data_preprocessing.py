import pandas as pd
import muon as mu
import scanpy as sc


def preprocess_omics_data(file_path):
    mdata = mu.read_10x_h5(file_path)
    mdata.var_names_make_unique()
    
    return mdata

def read_ann_data(file_path):
    ann_data = pd.read_csv(file_path, sep='\t', compression='gzip')
    return ann_data

# QUALITY CONTROL
def quality_control(mdata):
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

def merge_data(mdata):
    """
    Merges scRNA-seq and scATAC-seq data into a single MuData object.
    :param adata_rna: An AnnData object containing the scRNA-seq data.
    :param adata_atac: An AnnData object containing the scATAC-seq data.
    :return: A MuData object containing the merged datasets.
    """
    return mu.MuData({'rna': mdata['rna'], 'atac': mdata['atac']})

