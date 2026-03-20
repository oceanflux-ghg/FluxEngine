#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 06:13:09 2026

"""

import pytest
from os import path
import tempfile #for auto-deleting temporary directories
import subprocess
import sys

from netCDF4 import Dataset

from fluxengine.core.fe_setup_tools import get_fluxengine_root


def test_script_fe_run_basic():
    """
    Minimal test calling fe_run.py to run fluxengine. Checks that there are no errors and an output file with data variables is written
    """    
    with tempfile.TemporaryDirectory() as tmpDir:
        scriptPath = path.abspath(path.join(get_fluxengine_root(), "scripts", "fe_run.py"))
        configPath = path.abspath(path.join(get_fluxengine_root(), "test_data", "test_config_valid_full.conf"))
        
        result = subprocess.run([sys.executable, scriptPath, configPath,
                        "-l", #process indicator layers off
                        "-S1", #only run a single time point
                        "-o", tmpDir, #overwrite output directory
                        ],
                        check=True)
        
        #Script exited correctly without issue. Note: if exception occurred, test will automatically fail.
        assert result.returncode == 0
        
        #Output file exists
        expectedOutputPath = path.join(tmpDir, "2010", "01", "OceanFluxGHG-month01-jan-2010-v0.nc")
        assert path.exists(expectedOutputPath)
        
        #Ocean gas flux output variable exists
        dataset = Dataset(expectedOutputPath, "r")
        assert "OF" in dataset.variables.keys()


def test_script_fe_run_custom_gtv():
    """
    Minimal test calling fe_run.py to run fluxengine using a custom gas transfer velocity functor
    """
    
    pass

    with tempfile.TemporaryDirectory() as tmpDir:
        scriptPath = path.abspath(path.join(get_fluxengine_root(), "scripts", "fe_run.py"))
        configPath = path.abspath(path.join(get_fluxengine_root(), "test_data", "test_config_valid_custom_gtv.conf"))
        customGtvPath = path.abspath(path.join(get_fluxengine_root(), "test_data", "example_custom_gtv.py"))
        
        
        result = subprocess.run([sys.executable, scriptPath, configPath,
                        "-l", #process indicator layers off
                        "-S1", #only run a single time point
                        "-o", tmpDir, #overwrite output directory
                        "-custom_gas_transfer_parameterisation", customGtvPath, #path containing the custom gas transfer velocity implementation file
                        ],
                        check=True)
        
        #Script exited correctly without issue. Note: if exception occurred, test will automatically fail.
        assert result.returncode == 0
        
        #Output file exists
        expectedOutputPath = path.join(tmpDir, "2010", "01", "OceanFluxGHG-month01-jan-2010-v0.nc")
        assert path.exists(expectedOutputPath)
        
        #Ocean gas flux output variable exists
        dataset = Dataset(expectedOutputPath, "r")
        assert "OF" in dataset.variables.keys()
        
    
    
    
    


    


#TODO: These functions are not (yet?) tested
# match_filenames
 #fe_obj_from_run_parameters
# generate_datetime_points
# run_fluxengine
# process_timestep



if __name__ == "__main__":
    pytest.main(["-v", "--pdb"])