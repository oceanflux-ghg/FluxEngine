#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 06:24:22 2026

"""

import pytest

import fluxengine.core.fe_core


def test_RunParameters():
    """
    Run parameters are copied and updated correctly
    """
    pass


def test_write_netcdf():
    """
    ...
    """
    pass


def test_calculate_solubility_distilled():
    """
    ...
    """
    pass


def test_calculate_whitecapping():
    """
    ...
    """
    pass


def test_add_noise():
    """
    ...
    """
    pass




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
