import pandas as pd
import muon as mu
import scanpy as sc
from muon import atac as ac
import numpy as np
import os
import warnings

class DataPreprocessor():
    def __init__(self, h5_file_path, ann_file_path) -> None:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.mdata = mu.read_10x_h5(h5_file_path)
        print("After reading object has ", self.mdata.shape[0], "observations", self.mdata.shape[1], "variables")
        self.mdata.var_names_make_unique()
        print("After make_unique object has ", self.mdata.shape[0], "observations", self.mdata.shape[1], "variables")
        self.quality_control()
        print("After quality control object has ", self.mdata.shape[0], "observations", self.mdata.shape[1], "variables")
        self.normalize_data()
        print("After normalization object has ", self.mdata.shape[0], "observations", self.mdata.shape[1], "variables")
        self.merge_ann_data(ann_file_path)
        print("After merging with annotation object has ", self.mdata.shape[0], "observations", self.mdata.shape[1], "variables")

    def quality_control(self) -> mu.MuData:
        """
        Apply quality control measures to multimodal (RNA and ATAC) data.

        Parameters:
        - mdata (mu.MuData): MuData object containing RNA and ATAC data.

        Returns:
        - mu.MuData: MuData object after applying quality control.
        """
        mdata = self.mdata

        # Compute QC metrics for RNA
        mdata['rna'].var['mt'] = mdata['rna'].var_names.str.startswith('MT-')
        sc.pp.calculate_qc_metrics(mdata['rna'], qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

        mu.pp.filter_var(mdata['rna'], 'n_cells_by_counts', lambda x: x >= 3)
        mu.pp.filter_obs(mdata['rna'], 'n_genes_by_counts', lambda x: (x >= 200) & (x < 5000))
        mu.pp.filter_obs(mdata['rna'], 'total_counts', lambda    x: x < 15000)
        mu.pp.filter_obs(mdata['rna'], 'pct_counts_mt', lambda x: x < 20)
        
        # Compute QC metrics for ATAC
        sc.pp.calculate_qc_metrics(mdata['atac'], percent_top=None, log1p=False, inplace=True)
        # Filter ATAC based on quality metrics
        mu.pp.filter_var(mdata['atac'], 'n_cells_by_counts', lambda x: x >= 50)
        mu.pp.filter_obs(mdata['atac'], 'n_genes_by_counts', lambda x: (x >= 2000) & (x <= 15000))
        mu.pp.filter_obs(mdata['atac'], 'total_counts', lambda x: (x >= 4000) & (x <= 40000))
        # Intersect observations to keep only cells present in both modalities
        mu.pp.intersect_obs(mdata)


    def normalize_data(self) -> mu.MuData:
        """
        Normalize and Scale both scRNA-seq and scATAC-seq data.
        
        Parameters:
        - mdata (mu.MuData): MuData object containing RNA and ATAC data.
        
        Returns:
        - mu.MuData: A MuData object containing normalized and scaled data.
        """
        mdata = self.mdata

        mdata['rna'].layers["counts"] = mdata['rna'].X.copy()
        sc.pp.normalize_total(mdata['rna'], target_sum=1e4)
        sc.pp.log1p(mdata['rna'])

        sc.pp.highly_variable_genes(mdata['rna'], min_mean=0.02, max_mean=4, min_disp=0.5)
        sc.pp.scale(mdata['rna'], max_value=10)
        sc.tl.pca(mdata['rna'], svd_solver='arpack')

        mdata['atac'].layers["counts"] = mdata['atac'].X
        ac.pp.tfidf(mdata['atac'], scale_factor=None)
        ac.tl.lsi(mdata['atac'])
        # mdata['atac'].obsm['X_lsi'] = mdata['atac'].obsm['X_lsi'][:,1:]
        # mdata['atac'].varm["LSI"] = mdata['atac'].varm["LSI"][:,1:]
        # mdata['atac'].uns["lsi"]["stdev"] = mdata['atac'].uns["lsi"]["stdev"][1:]
        

    def merge_ann_data(self, file_path: str) -> mu.MuData:
        """
        Merge CellRanger peak annotation file to the existing MuData file
        
        Parameters:
        - file_patgh (str): Path to CellRanger peak annotation file
        - mdata (mu.MuData): MuData object containing RNA and ATAC data.
        
        Returns:
        - mu.MuData: A MuData object containing normalized and scaled data.
        """

        # read cellranger peak annotation with vectorized interval creation
        peak_annotation = pd.read_csv(file_path, sep='\t')
        peak_annotation['interval'] = peak_annotation['chrom'] + ':' + peak_annotation['start'].astype(str) + '-' + peak_annotation['end'].astype(str)

        # Using groupby to handle duplicates and aggregate values
        aggregated = peak_annotation.groupby('interval').agg({
        'gene': lambda x: ';'.join(x.dropna().astype(str)),
        'distance': lambda x: ';'.join(map(str, x.dropna())),
        'peak_type': lambda x: ';'.join(x.dropna().astype(str))
        }).reset_index()

        # Merge directly without creating extra DataFrames
        mdata = self.mdata
        mdata['atac'].var = pd.merge(mdata['atac'].var, aggregated, on='interval', how='left')

        # Update the var_names if necessary
        # for col in ['gene', 'distance', 'peak_type']:
        #    if col in mdata['atac'].var.columns:
        #        mdata['atac'].var[col] = mdata['atac'].var[col].apply(lambda x: str(x) if not pd.isna(x) else '')
    
    # def save_data(self):
    #     self.mdata.write('mudata.h5mu')
    #     file_dir = os.path.dirname(__file__)
    #     save_path = os.path.join(file_dir, 'mudata.h5mu')
    #     return save_path



def preprocess_omics_data(h5_file_path, ann_file_path):
    processor = DataPreprocessor(h5_file_path, ann_file_path)
    return processor.mdata
    

#preprocess_omics_data("\\Users\Dorzhey\\OneDrive\\Desktop\\lab\\projects\\re_design\\10x_data\\pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix.h5",
  #                    "\\Users\Dorzhey\\OneDrive\\Desktop\\lab\\projects\\re_design\\10x_data\\pbmc_granulocyte_sorted_3k_atac_peak_annotation.tsv")