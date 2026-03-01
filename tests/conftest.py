#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 28 14:45:21 2026

Contains common fixtures / other reused testing tools

"""

import pytest

#Minimal mocked DataLayerMetaData object
@pytest.fixture(name="mockDataLayerMetaData")
def mock_datalayer_metadata():
    class MockDataLayerMetaData:
        def __init__(self):
            self.name = "mocked_DataLayerMetaData"
            self.minBound = None
            self.maxBound = None
    return MockDataLayerMetaData()

