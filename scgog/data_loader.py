import torch
from torch.utils.data import Dataset, DataLoader
import numpy as np
import scanpy as sc
import muon as mu
import pandas as pd

class MultiModalDataset(Dataset):
    """
    A PyTorch Dataset for loading multimodal (scRNA-seq and scATAC-seq) data with Leiden clusters as labels,
    integrating RNA and ATAC data using weighted nearest neighbors (WNN).
    """
    def __init__(self, mdata):
        """
        Initializes the dataset with multimodal data and computes Leiden clusters using WNN integration as labels.

        Parameters:
        - mdata (mu.MuData): The multimodal data object processed by DataPreprocessor.
        """

        if "leiden_wnn" not in mdata.obs:
            self.mdata = self.compute_clusters(mdata)
        else:
            self.mdata = mdata
        # Prepare features and labels
        self.features, self.labels = self.prepare_data()

    def compute_clusters(self, mdata):
        """
        Computes Leiden clusters for the dataset using WNN integration of RNA and ATAC data.
        """
        # Compute neighbors for RNA
        sc.pp.neighbors(mdata['rna'], n_neighbors=10, n_pcs=20)
        # Compute neighbors for ATAC using LSI representation
        sc.pp.neighbors(mdata['atac'], use_rep="X_lsi", n_neighbors=10, n_pcs=20)
        # Run multimodal nearest neighbors on the graph
        mu.pp.neighbors(mdata, key_added='wnn')

        # Run multimodal UMAP on the WNN
        mu.tl.umap(mdata, neighbors_key='wnn', random_state=10)
        # Save the UMAP coordinates
        mdata.obsm["X_wnn_umap"] = mdata.obsm["X_umap"]
        # Cluster the cells with Leiden algorithm
        sc.tl.leiden(mdata, resolution=.3, neighbors_key='wnn', key_added='leiden_wnn')
        
        return mdata

    def prepare_data(self):
        """
        Prepares the features and labels for the dataset based on the computed clusters.
        """
        # Features can be derived from concatenated RNA and ATAC data or other representations
        # Here, we assume UMAP coordinates or original modalities can be used as features
        # This is a placeholder and should be tailored based on the actual data analysis needs
        
        # Placeholder for feature extraction, modify as needed
        features = np.concatenate([self.mdata['rna'].obsm['X_pca'], self.mdata['atac'].obsm['X_lsi']], axis=1)

        # Extract cluster labels
        labels = self.mdata.obs['leiden_wnn'].to_numpy().astype(int)
        
        
        
        return features, labels

    def __len__(self):
        """
        Returns the size of the dataset.
        """
        return self.features.shape[0]

    def __getitem__(self, idx):
        """
        Returns a single sample from the dataset.
        """
        return torch.tensor(self.features[idx], dtype=torch.float), torch.tensor(self.labels[idx], dtype=torch.long)
    

class PeakTFDataset(Dataset):
    def __init__(
        self,
        datafile,
        seqlen=300,
        subset=None,
        label = "cluster"
    ):

        data = pd.read_csv(datafile, nrows=subset) if subset is not None else pd.read_csv(datafile)
        x_train_seq = np.array([self.one_hot_encode_dna_sequence(x, seqlen) for x in data.sequence if "N" not in x])
        self.sequence = np.array([x.tolist() for x in x_train_seq])

        if label not in data.columns:
            raise Exception("No such label")
        self.label = np.array(data[label])
        if label == "peaktype":
            self.encoding, self.label = np.unique(self.label, return_inverse=True)

        print("loaded")
        
    def one_hot_encode_dna_sequence(self,seq, length):
        # Symmetrically trim the sequence to the specified length
        if len(seq) > length:
            trim_start = (len(seq) - length) // 2
            trim_end = trim_start + length
            seq = seq[trim_start:trim_end]
        
        # Define the mapping of nucleotides to one-hot encoding
        nucleotide_map = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        
        # Initialize the one-hot encoded sequence
        one_hot_seq = np.zeros(length, dtype=int)
        
        # Fill in the one-hot encoded sequence
        for i, nucleotide in enumerate(seq):
            if nucleotide in nucleotide_map:
                one_hot_seq[i] = nucleotide_map[nucleotide]
        
        return one_hot_seq

    def __len__(self):
        return(self.sequence.shape[0])
        
    def __getitem__(self,idx):
        return self.sequence[idx], self.label[idx]
    
def get_loader(data, use_sequence_data= False, label="cluster", batch_size=128, shuffle=True):
    """
    Creates a DataLoader for the multimodal dataset with Leiden clusters computed using WNN as labels.

    Parameters:
    - mdata (mu.MuData): The multimodal data object prepared for clustering.
    - batch_size (int): Batch size for the DataLoader.
    - shuffle (bool): Whether to shuffle the dataset.

    Returns:
    - DataLoader: The DataLoader for the dataset.
    """
    dataset = PeakTFDataset(data, label=label) if use_sequence_data else MultiModalDataset(data)

    return DataLoader(dataset, batch_size=batch_size, shuffle=shuffle)

def get_dataset(data, use_sequence_data= False, label="cluster"):
    return PeakTFDataset(data, label=label) if use_sequence_data else MultiModalDataset(data)
