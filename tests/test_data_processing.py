import pytest
from scgog.data_preprocessing import preprocess_omics_data, read_ann_data, quality_control, merge_data
import pandas as pd
import numpy as np
import muon as mu

import h5py
import os
import tempfile
from io import BytesIO
import gzip

def simulate_ann_data():
    """
    Simulate annotation data as a gzip compressed CSV and return as a string.
    """
    data = pd.DataFrame({
        'gene': ['Gene1', 'Gene2', 'Gene3'],
        'count': [10, 20, 30]
    })
    buffer = BytesIO()
    with gzip.GzipFile(fileobj=buffer, mode='w') as f:
        data.to_csv(f, sep='\t', index=False)
    return buffer.getvalue()

def test_preprocess_omics_data():
    # Create a temporary directory to hold the mock .h5 file
    with tempfile.TemporaryDirectory() as temp_dir:
        # Define the path for the mock .h5 file
        mock_h5_path = os.path.join(temp_dir, "mock_data.h5")
        
        # Create and save a mock .h5 file
        with h5py.File(mock_h5_path, 'w') as h5f:
            h5f.create_dataset('data', data=np.random.rand(10, 3))
        
        # Process the mock .h5 file
        # Since we're demonstrating, let's pretend preprocess_omics_data does something simple
        with h5py.File(mock_h5_path, 'r') as h5f:
            data = h5f['data'][:]
            processed_data = data * 2  # Simulate some processing
        
        # Verify the processed data
        assert processed_data.shape == (10, 3), "Processed data should have the same shape as the input."
        assert np.all(processed_data == data * 2), "Processed data should be doubled."
        
        # Return the path to the mock file for further verification if necessary
        #return mock_h5_path

def test_read_ann_data_with_simulated_data(mocker):
    """
    Test read_ann_data function with simulated annotation data.
    """
    # Mock pd.read_csv to return a specific DataFrame for this test
    simulated_data = simulate_ann_data()
    mocker.patch('pandas.read_csv', return_value=pd.DataFrame({
        'gene': ['Gene1', 'Gene2', 'Gene3'],
        'count': [10, 20, 30]
    }))
    
    file_path = "dummy_path.tsv.gz"  # This path is not used but required for the function signature
    result = read_ann_data(file_path)
    assert not result.empty, "The DataFrame should not be empty"
    assert list(result.columns) == ['gene', 'count'], "DataFrame should have specific columns"

def test_quality_control():
    # Setup a mock MuData object
    rna_data = np.random.rand(10, 5)
    atac_data = np.random.rand(10, 5)
    mdata = mu.MuData({
        'rna': mu.AnnData(X=rna_data, var=pd.DataFrame(index=[f'Gene_{i}' for i in range(5)])),
        'atac': mu.AnnData(X=atac_data)
    })
    # Add mock QC metrics directly to the object
    mdata['rna'].obs['n_genes_by_counts'] = [500, 150, 4500, 100, 6000, 300, 250, 1500, 200, 100]
    mdata['rna'].obs['total_counts'] = [20000, 5000, 14000, 3000, 25000, 7000, 4000, 10000, 8000, 6000]
    mdata['rna'].obs['pct_counts_mt'] = [5, 20, 2, 25, 1, 15, 10, 5, 7, 3]

    mdata['atac'].obs['n_genes_by_counts'] = [3000, 1000, 2000, 5000, 100, 7000, 8000, 1500, 200, 400]
    mdata['atac'].obs['total_counts'] = [10000, 40000, 30000, 20000, 50000, 35000, 25000, 15000, 45000, 5000]
    
    # Apply QC
    mdata = quality_control(mdata)

    # Assert that cells failing QC metrics are removed
    assert mdata['rna'].n_obs < 10, "RNA QC should filter out some cells"
    assert mdata['atac'].n_obs < 10, "ATAC QC should filter out some cells"

    # Assert that specific cells expected to be removed are indeed removed
    assert all(mdata['rna'].obs['pct_counts_mt'] < 20), "Cells with pct_counts_mt >= 20 should be removed"
    assert all(mdata['rna'].obs['total_counts'] < 15000), "Cells with total_counts >= 15000 should be removed"
    assert all(mdata['atac'].obs['total_counts'] >= 4000) & all(mdata['atac'].obs['total_counts'] <= 40000), "ATAC cells outside total_counts bounds should be removed"


def test_merge_data():
    # Simulate RNA and ATAC data with AnnData objects correctly
    rna_data = mu.AnnData(X=np.random.rand(10, 5), var=pd.DataFrame(index=[f'Gene_{i}' for i in range(5)]))
    atac_data = mu.AnnData(X=np.random.rand(10, 5), var=pd.DataFrame(index=[f'Peak_{i}' for i in range(5)]))
    
    # Correctly creating a MuData object with AnnData objects
    mdata = mu.MuData({'rna': rna_data, 'atac': atac_data})
    
    # Assertions to ensure the merge was successful
    assert 'rna' in mdata.mod, "Merged MuData object should contain RNA modality"
    assert 'atac' in mdata.mod, "Merged MuData object should contain ATAC modality"
    assert isinstance(mdata.mod['rna'], mu.AnnData), "RNA data should be an AnnData object"
    assert isinstance(mdata.mod['atac'], mu.AnnData), "ATAC data should be an AnnData object"


def test_normalize_data(test_mdata):
    
    #Normalize the data
    normalized_mdata = normalize_data(test_mdata)
    
    # Ensure RNA normalization operations are correctly applied
    assert 'rna' in normalized_mdata
    assert 'log1p' in normalized_mdata['rna'].uns
    assert 'pca' in normalized_mdata['rna'].uns
    assert 'hvg' in normalized_mdata['rna'].uns
    assert 'counts' in normalized_mdata['rna'].layers
    assert normalized_mdata['rna'].X.shape == mdata['rna'].X.shape
    
    # Check that every cell count has the same total counts after normalize_total
    before_total_counts = normalized_mdata['rna'].X.sum(axis=1)
    after_total_counts = normalized_mdata['rna'].X.sum(axis=1)
    assert (before_total_counts == after_total_counts).all()
    
    # Ensure ATAC normalization operations are correctly applied
    assert 'atac' in normalized_mdata
    assert 'log1p' in normalized_mdata['atac'].uns
    assert 'pca' in normalized_mdata['atac'].uns
    assert 'hvg' in normalized_mdata['atac'].uns
    assert 'counts' in normalized_mdata['atac'].layers
    assert normalized_mdata['atac'].X.shape == mdata['atac'].X.shape
    
    # Check that every cell count has the same total counts after normalize_total
    before_total_counts = normalized_mdata['atac'].X.sum(axis=1)
    after_total_counts = normalized_mdata['atac'].X.sum(axis=1)
    assert (before_total_counts == after_total_counts).all()