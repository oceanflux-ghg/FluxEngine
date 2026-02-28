#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 27 19:02:31 2026

"""


import pytest
import numpy as np
import fluxengine.core.data_preprocessing as fe_preprocessing
from fluxengine.core.datalayer import DataLayer
#from test_tools.mock import mock_datalayer_metadata


def test_transpose_basic(mockDataLayerMetaData):
    """
    Check that data is correctly transposed. Fdata is updated. Dimensions are correct
    """
    testData = np.arange(1, 7).reshape(2, 3)
    testDataLayer = DataLayer("test_data", testData.copy(), mockDataLayerMetaData, -999.9)
    fe_preprocessing.transpose(testDataLayer)
    #dimensions transposed
    assert(testDataLayer.nx == 2) #Note: FluxEngine uses column first indexing (possibly to match netCDF4?)
    assert(testDataLayer.ny == 3)
    assert(np.all(testDataLayer.data == np.array([[1,4],[2,5],[3,6]]))) #transpose works as expected
    assert(np.all(testDataLayer.fdata == np.array([[1,4],[2,5],[3,6]]).ravel())) #fdata view is updated





# @pytest.fixture
# def mock_datalayer_simple():
#     class MockDataLayer:
#         def __init__(self, fdata:np.ndarray, name:str="test_datalayer", missing_value:float=DataLayer.missing_value):
#             self.name = "test_datalayer"
#             self.missing_value = missing_value
#             self.fdata = np.array(fdata, copy=True)
#     return MockDataLayer




#### Testing common 'contract' for unit conversions
#Many of the preprocessing functions have the same required behaviour, so test these using
# parameterised tests:
#   1) input is transformed to expected output
#   2) conversion is reflected in the shaped 'data' field, and not just the flat view of the data
#   3) missing_values are not changed
#conversion_func: function being tested
#inputData, expectedOutput: the input and expected output for the function being tested
#hasStableElements: True, if elements do not move (e.g. no flip or roll). When False, missing values must be tracked with transformation.
@pytest.mark.parametrize(
    "conversion_func, hasStableElements, inputData, expectedOutput",
        [
            (
                fe_preprocessing.kelvin_to_celsius, True,
                np.array([[0.0, 273.15], [280.0, 300.0]]),
                np.array([[-273.15, 0.0], [6.85, 26.85]])
            ),
            (
                fe_preprocessing.celsius_to_kelvin, True,
                np.array([[0.0, -273.15], [6.85, 26.85]]),
                np.array([[273.15, 0.0], [280.0, 300.0]])
            ),
            (
                fe_preprocessing.pascal_to_millibar, True,
                np.array([[101000, 101325], [102020, 98011]], dtype=float),
                np.array([[1010.0, 1013.25], [1020.2, 980.11]])
            ),
            (
                fe_preprocessing.percent_to_proportion, True,
                np.array([[0.0, 100.0], [100.0/3, 75.0]]),
                np.array([[0.0, 1.0], [1.0/3, 0.75]])
            ),
            (
                fe_preprocessing.pow2, True,
                np.array([[-1, 2.0, 3.0, 10.0], [-5.0, 2.5, 0.42, -0.09]], dtype=float),
                np.array([[1, 4.0, 9.0, 100.0], [25.0, 6.25, 0.42**2, (-0.09)**2]], dtype=float)
            ),
            (
                fe_preprocessing.pow3, True,
                np.array([[-1, 2.0, 3.0, 10.0], [-5.0, 2.5, 0.42, -0.09]], dtype=float),
                np.array([[-1, 8.0, 27.0, 1000.0], [-125.0, 15.625, 0.42**3, (-0.09)**3]], dtype=float)
            ),
            (
                #Note: daytohour is named poorly. It does the opposite, converting horus to days. TODO: make github issue...
                fe_preprocessing.daytohour, True,
                np.array([[24, 12, 1], [48, 168, 732]], dtype=float),
                np.array([[1, 0.5, 1/24], [2, 7, 30.5]], dtype=float)
            ),
            (
                fe_preprocessing.flip_longitude, False,
                np.array([[1, 2, 3], [10, 20, 30]], dtype=float),
                np.array([[10, 20, 30], [1, 2, 3]], dtype=float)
            ),
            (
                fe_preprocessing.flip_latitude, False,
                np.array([[1.0, 2.0, 3.0], [1.5, 2.5, 3.5]], dtype=float),
                np.array([[3.0, 2.0, 1.0], [3.5, 2.5, 1.5]], dtype=float)
            ),
            
            (
                fe_preprocessing.longitude_roll_180, False,
                np.array([[1, 2, 3, 4], [5, 6, 7, 8]], dtype=float),
                np.array([[3, 4, 1, 2], [7, 8, 5, 6]], dtype=float)
            ),
            
        ]
    )

def test_unit_conversion_common_contract(conversion_func, hasStableElements, inputData, expectedOutput, mockDataLayerMetaData):
    """
    Tests common behavioural requirements for unit conversion functions:
        1) input is transformed to expected output
        2) conversion is reflected in the shaped 'data' field, and not just the flat view of the data
        3) missing_values are not changed
    """
    
    expectedFlatOutput = expectedOutput.flatten()
    
    #1) correct conversion (in shaped 'data' field)
    testData = DataLayer("test_data", inputData.copy(), mockDataLayerMetaData, DataLayer.missing_value)
    #testData = mock_datalayer_simple(fdata=inputData.copy())
    conversion_func(testData)
    assert np.allclose(testData.data, expectedOutput)
    
    #2) updates propagated to the flat data
    assert np.allclose(testData.fdata, expectedFlatOutput)
    
    #3) missing values remain unchanged
    inputDataWithMissing = inputData.copy()
    inputDataWithMissing[0, 0] = DataLayer.missing_value
    inputDataWithMissing[-1, -1] = DataLayer.missing_value
    testData = DataLayer("test_data", inputDataWithMissing, mockDataLayerMetaData, DataLayer.missing_value)
    #testData = mock_datalayer_simple(fdata=inputDataWithMissing)
    wMissing = testData.data == testData.missing_value
    conversion_func(testData)
    if hasStableElements == False: #If the positions aren't stable, missing data values will move. Put the missing data locations through the same transformation to track the where they end up.
        missingLocs = DataLayer("track_missing", wMissing, mockDataLayerMetaData, DataLayer.missing_value)
        conversion_func(missingLocs)
        wMissing = missingLocs.data
    assert np.all(testData.data[wMissing] == testData.missing_value)
    assert np.allclose(testData.data[wMissing==False], expectedOutput[wMissing==False])
    
        

def test_lat_grid_lines_to_centre_of_cells_basic(mockDataLayerMetaData):
    """
    Tests expected output is correct (values should be means across a two element wide lattitude window)
    Tests that the dimensions are correct (n-1, m) where (n, m) are the input dimensions
    Tests that changes are propagated to fdata
    """
    inputData = np.array([[2, 20], [4, 40], [6, 60], [8, 80], [10, 100]], dtype=float)
    expectedOutput = np.array([[3, 30], [5, 50], [7, 70], [9, 90]], dtype=float)
    
    metadata = mockDataLayerMetaData
    testData = DataLayer("test_data", inputData.copy(), metadata, DataLayer.missing_value)
    fe_preprocessing.lat_grid_lines_to_centre_of_cells(testData)
    
    assert(testData.ny == expectedOutput.shape[0])
    assert(testData.nx == expectedOutput.shape[1])
    assert(np.allclose(testData.data, expectedOutput))
    assert(np.allclose(testData.fdata, expectedOutput.ravel()))





#def test_foc_to_epsilon
#def test_foc_to_epsilon_craig1994


if __name__ == "__main__":
    pytest.main(["-v", "--pdb"])
    #pytest.main()

    
    
    