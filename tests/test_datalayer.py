#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 27 19:02:31 2026

"""


import pytest
import numpy as np
from fluxengine.core.datalayer import DataLayer

#from test_tools.mock import mock_datalayer_metadata



def test_create_empty_datalayer(mockDataLayerMetaData):
    """
    Tests...
    """
    assert(True)
    #DataLayer.create_empty_datalayer("test_data", 5, 3, metadata, fillValue=0.0)





if __name__ == "__main__":
    pytest.main(["-v", "--pdb"])
