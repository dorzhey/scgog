import pytest
import numpy as np
import muon as mu
import pandas as pd
import torch
from torch.utils.data import DataLoader, Dataset
from scgog.data_loader import MultiModalDataset, PeakTFDataset  # Adjust the import path based on your project structure

# Create a fixture for a mock MuData object
@pytest.fixture
def mock_mudata():
    # Simulate minimal RNA and ATAC data with required structure
    rna_data = np.random.rand(10, 5)
    atac_data = np.random.rand(10, 5)
    mdata = mu.MuData({
        'rna': mu.AnnData(rna_data, obs={'leiden_wnn': np.random.randint(0, 3, 10)}),
        'atac': mu.AnnData(atac_data)
    })
    return mdata

# Create a fixture for a mock CSV file path for sequence data
@pytest.fixture
def mock_csv(tmp_path):
    df = pd.DataFrame({
        'sequence': ['ACTG' * 75, 'TACG' * 75],  # Example sequences
        'cluster': [1, 2]
    })
    csv_path = tmp_path / "data.csv"
    df.to_csv(csv_path, index=False)
    return str(csv_path)

### Testing MultiModalDataset
def test_multimodal_dataset_initialization(mock_mudata):
    """Test initialization of MultiModalDataset."""
    dataset = MultiModalDataset(mock_mudata)
    assert len(dataset) == 10  # Assuming 10 samples
    assert isinstance(dataset, Dataset)

def test_multimodal_dataset_getitem(mock_mudata):
    """Test the retrieval of a single data item."""
    dataset = MultiModalDataset(mock_mudata)
    features, labels = dataset[0]
    assert isinstance(features, torch.Tensor)
    assert isinstance(labels, torch.Tensor)
    assert labels.dtype == torch.int64  # Labels should be integers

### Testing PeakTFDataset
def test_peaktf_dataset_initialization(mock_csv):
    """Test initialization of PeakTFDataset with a mock CSV file."""
    dataset = PeakTFDataset(datafile=mock_csv)
    assert len(dataset) == 2  # Based on the provided sequences
    assert isinstance(dataset, Dataset)

def test_peaktf_dataset_getitem(mock_csv):
    """Test the retrieval of a single sequence item."""
    dataset = PeakTFDataset(datafile=mock_csv)
    sequence, label = dataset[0]
    assert isinstance(sequence, np.ndarray)
    assert label in [1, 2]  # Cluster labels

### Running Tests
# If using pytest to execute, include this in your test script
if __name__ == "__main__":
    pytest.main([__file__])