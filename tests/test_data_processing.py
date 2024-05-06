import pytest
from unittest.mock import patch, MagicMock
import pandas as pd
import numpy as np
from scgog.data_preprocessing import DataPreprocessor

# Setup mock data that reflects all necessary attributes
def create_mock_mudata():
    mock = MagicMock()
    mock.obs = {'leiden_wnn': np.random.randint(0, 3, size=10)}
    mock['rna'] = MagicMock()
    mock['atac'] = MagicMock()
    mock['rna'].var_names = ['gene1', 'gene2', 'gene3']
    mock['rna'].var = pd.DataFrame(index=mock['rna'].var_names)
    mock['rna'].obsm = {'X_pca': np.random.rand(10, 5)}
    mock['atac'].obsm = {'X_lsi': np.random.rand(10, 5)}
    return mock

@patch('scgog.data_preprocessing.mu.read_10x_h5', return_value=create_mock_mudata())
@patch('scgog.data_preprocessing.pd.read_csv', return_value=pd.DataFrame({'gene': ['gene1'], 'peak': [1]}))
def test_data_preprocessor_initialization(read_10x_h5_mock, read_csv_mock):
    processor = DataPreprocessor('/fake/path/h5file.h5', '/fake/path/annotation.tsv')
    assert processor is not None
    assert 'X_pca' in processor.mdata['rna'].obsm

def test_quality_control():
    processor = DataPreprocessor()
    processor.mdata = create_mock_mudata()
    processor.quality_control()
    assert 'n_cells_by_counts' in processor.mdata['rna'].var

# Run the test
if __name__ == "__main__":
    pytest.main([__file__])
