import torch
from torch.utils.data import Dataset, DataLoader
import numpy as np
from muon import MuData

class MultiModalDataset(Dataset):
    """
    A PyTorch Dataset for loading multimodal (scRNA-seq and scATAC-seq) data.
    """
    def __init__(self, mdata, label_type='gene'):
        """
        Initializes the dataset with multimodal data.

        Parameters:
        - mdata (mu.MuData): The multimodal data object processed by DataPreprocessor.
        - label_type (str): The type of label to use ('gene' or 'peak_type').
        """
        super(MultiModalDataset, self).__init__()
        self.mdata = mdata
        self.label_type = label_type
        # Prepare features and labels
        self.features, self.labels = self.prepare_data()

    def prepare_data(self):
        """
        Prepares the features and labels for the dataset.
        """
        # Extracting RNA and ATAC data as features
        rna_data = self.mdata['rna'].X.toarray()  # Convert sparse matrix to dense
        atac_data = self.mdata['atac'].X.toarray()
        features = np.concatenate([rna_data, atac_data], axis=1)

        # Extracting labels
        if self.label_type == 'gene':
            labels = np.array(self.mdata['atac'].var['gene'].values.tolist())
        else:  # peak_type
            labels = np.array(self.mdata['atac'].var['peak_type'].values.tolist())

        # Encoding labels as integers for classification
        from sklearn.preprocessing import LabelEncoder
        le = LabelEncoder()
        labels_encoded = le.fit_transform(labels)
        
        return features, labels_encoded

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
    
def get_loader(mdata:MuData, label:str, batch_size=32, shuffle=True):
    """
    Creates a DataLoader with clustering-based labels.

    Parameters:
    - mdata (str): Path to the .h5 file.
    - batch_size (int): Batch size for the DataLoader.
    - shuffle (bool): Whether to shuffle the dataset.

    Returns:
    - DataLoader: The DataLoader for the dataset.
    """
    dataset = MultiModalDataset(mdata, label)
    return DataLoader(dataset, batch_size=batch_size, shuffle=shuffle)

# Usage example
# Assuming 'processed_mdata' is your processed MuData object from DataPreprocessor
