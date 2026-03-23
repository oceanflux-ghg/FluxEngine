#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 06:24:22 2026

"""

import pytest
from os import path

import fluxengine.core.fe_core as fe_core


def test_RunParameters():
    """
    Run parameter objects are created and updated correctly
    """
    #RunParameters are constructed with no parameters
    params = fe_core.RunParameters()
    assert len(vars(params)) == 0
    
    #Adding parameters adds them as attributes
    params.set_parameters({"testParam1": 42, "testParam2": 3.141})
    assert len(vars(params)) == 2
    assert params.testParam1 == 42
    assert params.testParam2 == 3.141
    
    #Setting parameters clears any previous parameters
    params.set_parameters({"exampleParam1": "some value"})
    assert hasattr(params, "testParam1") == False
    assert hasattr(params, "testParam2") == False
    assert len(vars(params)) == 1
    assert params.exampleParam1 == "some value"



# def test_write_netcdf(mockFluxEngineObject_outputs):
#     """
#     ...
#     """
#     # import tempfile #for auto-deleting temporary directories
#     # tmpDir = tempfile.TemporaryDirectory()
    
#     # mockFluxEngineObject_outputs.runParams.output_path = path.join(tmpDir.name, "test_output.nc")
    
#     # fe_core.write_netcdf(mockFluxEngineObject_outputs)
    
#     pass
    
    
    


# def test_calculate_solubility_distilled():
#     """
#     ...
#     """
#     pass


# def test_calculate_whitecapping():
#     """
#     ...
#     """
#     pass


# def test_add_noise():
#     """
#     ...
#     """
#     pass




#def test_add_noise_and_bias_wind
#def test_add_bias_k_biology_wind
#def test_ass_sst_rain_bias
#def test_median_filter2D

#Note: These tests check against function output generated from 2026.03.02 and will only identify if output has changed from this baseline.
#      It would be better to compare to published empirical data, e.g. from the original publications?
#def test_schmidt_Wanninkhot1992
#def test_schmidt_Wanninkhot2014
#def test_solubility_Wanninkhot1992
#def test_solubility_Wanninkhot2014


#def test_output_unit_molecular_mass
#def test_calculate_concw
#def test_calculate_conca
#def test_copy_missing_values
#def test_check_output_dataset
#def test_average_pixels
#def test_check_dimensions


#def test_ ... Next test FluxEngine class itself.





if __name__ == "__main__":
    pytest.main(["-v", "--pdb"])
