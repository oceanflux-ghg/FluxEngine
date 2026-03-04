#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 27 19:02:31 2026

"""


import pytest
import numpy as np
from os import path

from fluxengine.core.datalayer import DataLayerMetaData, DataLayer
import fluxengine.core.data_preprocessing as fe_preprocessing
from fluxengine.core.fe_setup_tools import get_fluxengine_root



def get_metadata_init_args():
    return {
        "name":"internal_name", "netCDFName":"nc_name", "units":"mol",
        "minBound":0, "maxBound":1000000,
        "standardName":"standard name",
        "longName":"long from description of the variable",
        "fillValue":DataLayer.missing_value
    }




def test_DataLayer_construction_basic():
    """
    Tests correct initialisation:
    1) metadata attributes are copied over
    2) data is copied over correctly, and has the correct dimensions
    3) fdata is created correctly
    4) missing values copied correctly, even if not using the DataLayer missing_value value
    
    
    """
    metadata = DataLayerMetaData(**get_metadata_init_args())
    data = np.array(([1, 2, 3], [4, 5, 6]), dtype=float)
    
    dataLayer = DataLayer(data.copy(), metadata, fillValue=DataLayer.fill_value)
    
    #1) Metadata correctly copied over to DataLayer
    for attrName in vars(metadata):
        assert getattr(dataLayer, attrName) == getattr(metadata, attrName)
    
    #2) Data is copied correctly, has the correct dimensions
    assert dataLayer.ny == data.shape[0]
    assert dataLayer.nx == data.shape[1]
    assert np.allclose(dataLayer.data, data)
    
    #3) fdata is copied correctly
    assert np.allclose(dataLayer.fdata, data.ravel())
    
    #4) Missing values are correctly copied and transformed to the DataLayer's missing_value value
    #   e.g. using np.nan
    np.array(([np.nan, 2, 3], [4, 5, np.nan]), dtype=float)
    dataLayer = DataLayer(data.copy(), metadata, fillValue=np.nan)
    expectedOutput = np.where(np.isnan(data), DataLayer.missing_value, data)
    assert np.allclose(dataLayer.data, expectedOutput)
    assert np.allclose(dataLayer.fdata, expectedOutput.ravel())



def test_DataLayer_construction_min_max_masking():
    """
    Tests DataLayers are constructed and mask/copies mssing_value where data is outside min/max range
    """
    
    #When both min and max are provided
    data = np.array(([0.0, -0.0, 1.0], [-0.1, 1.1, 0.5]), dtype=float)
    expectedOutput = np.where(np.logical_or(data<0, data>1), DataLayer.missing_value, data)
    
    metadata = DataLayerMetaData(**get_metadata_init_args())
    metadata.minBound = 0.0 #inclusive
    metadata.maxBound = 1.0 #also inclusive
    dataLayer = DataLayer(data.copy(), metadata, fillValue=DataLayer.fill_value)
    assert np.allclose(dataLayer.data, expectedOutput)
    assert np.allclose(dataLayer.fdata, expectedOutput.ravel())
    
    #Only min bound is provided
    metadata.minBound = None
    data = np.array(([0.0, -0.0, 1.0], [-0.1, 1.1, 0.5]), dtype=float)
    expectedOutput = np.where(data>1, DataLayer.missing_value, data)
    dataLayer = DataLayer(data.copy(), metadata, fillValue=DataLayer.fill_value)
    assert np.allclose(dataLayer.data, expectedOutput)
    assert np.allclose(dataLayer.fdata, expectedOutput.ravel())
    
    #Only max bound is provided
    metadata.minBound = 0.0
    metadata.maxBound = None
    data = np.array(([0.0, -0.0, 1.0], [-0.1, 1.1, 0.5]), dtype=float)
    expectedOutput = np.where(data<0, DataLayer.missing_value, data)
    dataLayer = DataLayer(data.copy(), metadata, fillValue=DataLayer.fill_value)
    assert np.allclose(dataLayer.data, expectedOutput)
    assert np.allclose(dataLayer.fdata, expectedOutput.ravel())



def test_DataLayer_construction_applies_preprocessing():
    """
    Tests data preprocessing functions are correctly applied (when supplied)
    """
    preprocessingFunctions = [fe_preprocessing.pow2, fe_preprocessing.pow3]
    
    data = np.array(([0.0, DataLayer.missing_value, 2.0], [-2.0, 0.5, DataLayer.missing_value]), dtype=float)
    expectedOutput = np.where(data!=DataLayer.missing_value, (data**2)**3, data)
    
    metadata = DataLayerMetaData(**get_metadata_init_args())
    dataLayer = DataLayer(data.copy(), metadata, fillValue=DataLayer.fill_value, preprocessing=preprocessingFunctions)
    
    assert np.allclose(dataLayer.data, expectedOutput)
    assert np.allclose(dataLayer.fdata, expectedOutput.ravel())
    


def test_create_empty_datalayer():
    """
    Tests the creation of empty DataLayer:
    1) Correctly filled with missing_value
    2) Correct name and dimensions
    3) Correct metadata
    """
    
    metadata = DataLayerMetaData(**get_metadata_init_args())
    dataLayer = DataLayer.create_empty_datalayer(5, 3, metadata)
    
    #1) Correct fill value
    assert np.all(dataLayer.data == DataLayer.missing_value)
    assert np.all(dataLayer.fdata == DataLayer.missing_value)
    
    #2) Correct dimensions
    assert (dataLayer.nx, dataLayer.ny) == (5, 3)
    
    #3) Correct metadata
    for attrName in vars(metadata):
        assert getattr(dataLayer, attrName) == getattr(metadata, attrName)

    

def test_create_from_file_basic():
    """
    Test correct creation of DataLayer objects from netCDF files
    """
    dataNcPath = path.join(get_fluxengine_root(), "data", "verification_data", "SST","2010","20100101_OCF-SST-GLO-1M-100-ATS-ARC.nc")
    
    metadata = DataLayerMetaData(name="test_data_from_file")
    dataLayer = DataLayer.create_from_file(dataNcPath, "sst_skin_mean", metadata, 0)
    
    #Correct number of dimensions and name
    assert dataLayer.data.shape == (180, 360)
    assert dataLayer.name == "test_data_from_file"
    #Don't check every value, just mean and number of non-missing values match expectations
    assert np.allclose(np.mean(dataLayer.data[dataLayer.data != DataLayer.missing_value]), 290.0087384972042)
    assert np.sum(dataLayer.data != DataLayer.missing_value) == 36512
    
    #Check that invalid time dimension throws
    metadata = DataLayerMetaData(name="test_data_from_file")
    metadata.timeDimensionName = "invalid_time"
    with pytest.raises(RuntimeError):
        dataLayer = DataLayer.create_from_file(dataNcPath, "sst_skin_mean", metadata, 0)




if __name__ == "__main__":
    pytest.main(["-v", "--pdb"])
