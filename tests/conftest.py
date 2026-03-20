#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 28 14:45:21 2026

Contains common fixtures / other reused testing tools

"""

import pytest
import numpy as np

#Minimal mocked DataLayerMetaData object
@pytest.fixture(name="mockDataLayerMetaData")
def mock_datalayer_metadata():
    class MockDataLayerMetaData:
        def __init__(self):
            self.name = "mocked_DataLayerMetaData"
            self.minBound = None
            self.maxBound = None
    return MockDataLayerMetaData()




#Mock FluxEngine object which mocks FluxEngine's finished state, i.e. containing output data
@pytest.fixture(name="mockFluxEngineObject_outputs")
def mock_fluxengine_outputs():
    class MockRunConfig_outputs:
        def __init__(self):
            self.run_count = 0
            self.exclude_outputs = ""
            self.output_temporal_chunking = 0
            self.temporal_resolution = None
            self.year = 2010
            self.month = 1
            self.day = 1
            self.hour = 0
            self.minute = 0
            self.second = 0
            self.run_count = 0
            self.output_temporal_chunking = 1
            #self.time_data = ???
            self.output_path = None #Should be overwritten as needed!
    class MockFluxEngineObject_outputs:
        def __init__(self, runParams):
            self.time_data = 0
            self.data = np.arange(4*5, dtype=float)
            self.runParams = runParams
            self.nx = 4
            self.ny = 5
            self.longitude_data = np.arange(self.nx)
            self.latitude_data = np.arange(self.ny)
            self.longitude_grid = np.repeat(self.longitude_data, self.ny).reshape(self.ny, self.nx)
            self.latitude_grid = np.tile(self.latitude_data, self.nx).reshape(self.ny, self.nx)
            
    return MockFluxEngineObject_outputs(MockRunConfig_outputs());