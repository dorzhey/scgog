import pandas as pd
import muon as mu
import scanpy as sc
from muon import atac as ac
import numpy as np
import os

class DataPreprocessor():
    def __init__(self, h5_file_path, ann_file_path) -> None:
        self.mdata = mu.read_10x_h5(h5_file_path)
        self.mdata.var_names_make_unique()
        self.quality_control()
        print(self.mdata.shape)
        self.normalize_data()
        self.merge_ann_data(ann_file_path)

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
        mdata['atac'].obsm['X_lsi'] = mdata['atac'].obsm['X_lsi'][:,1:]
        mdata['atac'].varm["LSI"] = mdata['atac'].varm["LSI"][:,1:]
        mdata['atac'].uns["lsi"]["stdev"] = mdata['atac'].uns["lsi"]["stdev"][1:]
        self.mdata = mdata
        

    def merge_ann_data(self, file_path: str) -> mu.MuData:
        """
        Merge CellRanger peak annotation file to the existing MuData file
        
        Parameters:
        - file_patgh (str): Path to CellRanger peak annotation file
        - mdata (mu.MuData): MuData object containing RNA and ATAC data.
        
        Returns:
        - mu.MuData: A MuData object containing normalized and scaled data.
        """
        # read cellranger peak annotation
        peak_annotation = pd.read_csv(file_path, sep='\t')

        # add the interval to match atac var
        peak_annotation['interval'] = peak_annotation['chrom'] + ':' + peak_annotation['start'].astype(str) + '-' + peak_annotation['end'].astype(str)

        # merge dupes in cellranger peak annotation
        unique_interval = np.unique(peak_annotation['interval'])
        ind_col = []
        g_col = []
        d_col = []
        t_col = []
        for iu in unique_interval:
            ind = peak_annotation.index[peak_annotation['interval'] == iu].tolist()
            ind_col.append(iu)
            g_col.append(peak_annotation.iloc[ind]['gene'].tolist())
            d_col.append(peak_annotation.iloc[ind]['distance'].tolist())
            t_col.append(peak_annotation.iloc[ind]['peak_type'].tolist())

        pivot_annotation = {'interval': ind_col, 
                            'gene': g_col, 
                            'distance': d_col, 
                            'peak_type': t_col}
        pivot_annotation = pd.DataFrame(pivot_annotation)
        
        mdata = self.mdata
        # merge the cellranger peak annotation to the muon/anndata data
        mdata['atac'].var = mdata['atac'].var.merge(pivot_annotation, on='interval', how='left')
        # put back the interval var names for the atac modality
        mdata.mod['atac'].var_names = mdata.mod['atac'].var.interval
        self.mdata = mdata
    
    def save_data(self):
        self.mdata.write('mudata.h5mu')
        file_dir = os.path.dirname(__file__)
        save_path = os.path.join(file_dir, 'mudata.h5mu')
        return save_path



def preprocess_omics_data(h5_file_path, ann_file_path):
    processor = DataPreprocessor(h5_file_path, ann_file_path)
    return processor.save_data()
    

#preprocess_omics_data("\\Users\Dorzhey\\OneDrive\\Desktop\\lab\\projects\\re_design\\10x_data\\pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix.h5",
  #                    "\\Users\Dorzhey\\OneDrive\\Desktop\\lab\\projects\\re_design\\10x_data\\pbmc_granulocyte_sorted_3k_atac_peak_annotation.tsv")