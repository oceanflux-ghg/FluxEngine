#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 06:13:09 2026

"""

import pytest
from os import path
import datetime
import numpy as np

import fluxengine.core.fe_setup_tools as fe_setup_tools
from fluxengine.core.datalayer import DataLayer


#Test importing of configuration files


def test_get_fluxengine_root():
    """
    Checks that get_fluxengine_root returns a path from which the 'data', 'configs' and 'scripts' directories can be found
    """
    rootPath = fe_setup_tools.get_fluxengine_root()
    
    #data directory exists, and at least one expected data file can be found
    assert path.exists(path.join(rootPath, "data"))
    assert path.exists(path.join(rootPath, "data", "Longhurst-provinces-mask.nc"))
    
    #configs directory exists, and at least one expected data file can be found
    assert path.exists(path.join(rootPath, "configs"))
    assert path.exists(path.join(rootPath, "configs", "socatv4_sst_salinity_gradients-N00.conf"))
    
    #scripts directory exists, and at least one expected data file can be found
    assert path.exists(path.join(rootPath, "scripts"))
    assert path.exists(path.join(rootPath, "scripts", "fe_run.py"))
    
    

def test_read_config_file():
    """
    Tests that config files can be read and a dictionary containing the appropriate variables is returned
    """
    
    #Raises exception if config file doesn't exist
    with pytest.raises(ValueError):
        fe_setup_tools.read_config_file(path.join("none/existent/config.conf"))
    
    ##Malformed config files raise an exception
    #config variable is missing the assignment operator
    with pytest.raises(ValueError):
        fe_setup_tools.read_config_file(path.join(fe_setup_tools.get_fluxengine_root(), "test_data", "test_config_malformed_no_assignment.conf"))
    #config variable is missing a value after the assignment
    with pytest.raises(ValueError):
        fe_setup_tools.read_config_file(path.join(fe_setup_tools.get_fluxengine_root(), "test_data", "test_config_malformed_no_value.conf"))
    #config has an unexpected line which doesn't correspond to a known entry
    with pytest.raises(ValueError):
        fe_setup_tools.read_config_file(path.join(fe_setup_tools.get_fluxengine_root(), "test_data", "test_config_malformed_unexpected_line.conf"))
    #config variable uses an incorrect assignment operator
    with pytest.raises(ValueError):
        fe_setup_tools.read_config_file(path.join(fe_setup_tools.get_fluxengine_root(), "test_data", "test_config_malformed_wrong_assignment_operator.conf"))
    
    ######
    #Valid config
    configPath = path.join(fe_setup_tools.get_fluxengine_root(), "test_data", "test_config_valid_partial.conf")
    config = fe_setup_tools.read_config_file(configPath)
    assert (isinstance(config, dict)) #A dictionary was returned
    #All entries are keys
    for key in config.keys():
        assert isinstance(config[key], str)
    #All expected entries are present
    expectedConfig = {"flux_calc": "rapid",
                      "sst_gradients": "yes",
                      "cool_skin_difference": "0.17",
                      "saline_skin_value": "0.1",
                      "axes_data_layer": "sstskin",
                      "latitude_prod": "lat",
                      "longitude_prod": "lon",
                      "time_prod": "time",
                      "k_parameterisation": "k_Nightingale2000",
                      "output_dir": "test_output/test_config_output",
                      "sstskin_path": "<FEROOT>/data/verification_data/SST/<YYYY>/<YYYY><MM>01_OCF-SST-GLO-1M-100-ATS-ARC.nc",
                      "sstskin_prod": "sst_skin_mean",
                      "sstskinC_maxBound": "100",
                      "sstskinC_minBound": "-1.8",
                      "sstskin_maxBound": "400",
                      "sstskin_minBound": "271.35",
                      "sstfndC_maxBound": "100",
                      "sstfndC_minBound": "-1.8",
                      "pco2_sst_maxBound": "30.5",
                      "pco2_sst_minBound": "-1.8"}
    assert config == expectedConfig
    


def test_parse_config_version_tag():
    """
    Tests the configuration file version can be correctly parsed
    """
    
    #Version is parsed correctly
    assert fe_setup_tools.parse_config_version_tag("#?FluxEngineConfigVersion:1") == 1.0
    assert fe_setup_tools.parse_config_version_tag("#?FluxEngineConfigVersion:2") == 2.0
    #non-integer versions parsed correctly
    assert fe_setup_tools.parse_config_version_tag("#?FluxEngineConfigVersion:3.141") == 3.141
    #tolerant to whitespace
    assert fe_setup_tools.parse_config_version_tag("#?FluxEngineConfigVersion: 2.718") == 2.718
    assert fe_setup_tools.parse_config_version_tag("#?FluxEngineConfigVersion : 2.718") == 2.718
    assert fe_setup_tools.parse_config_version_tag("#?FluxEngineConfigVersion :2.718") == 2.718
    assert fe_setup_tools.parse_config_version_tag(" #?FluxEngineConfigVersion: 2.718") == 2.718
    assert fe_setup_tools.parse_config_version_tag("#?FluxEngineConfigVersion:\t2.718") == 2.718
    
    
    
def test_read_config_metadata():
    """
    Some minimal testing that the settings.xml metadata file is parsed correctly. Checks the existance of metadata for key configuration file and data layer variables
    """
    
    metadata = fe_setup_tools.read_config_metadata(path.join(fe_setup_tools.get_fluxengine_root(), "core", "settings.xml"))
    
    #Config variable metadata
    assert metadata["output_dir"] == {"name": "output_dir", "required": "true", "type": "path"}
    assert metadata["latitude_prod"] == {"name": "latitude_prod", "required": "true", "type": "string"}
    assert metadata["longitude_prod"] == {"name": "longitude_prod", "required": "true", "type": "string"}
    assert metadata["flux_calc"] == {"name": "flux_calc", "required": "true", "type": "multioption", "options": {"rapid": 1, "equilibrium": 2, "bulk": 3}}
    assert metadata["exclude_outputs"] == {"name": "exclude_outputs", "required": "false", "type": "string", "default": ""}
    
    
    #Check for existence of some data layers, and the metadata for their associated variables
    for dataLayerName in ["pco2_sw", "windu10", "ice", "pressure", "salinity"]:
        assert metadata[dataLayerName+"_path"] == {"name": dataLayerName, "required": "true", "type": "DataLayerPath"}
        assert metadata[dataLayerName+"_prod"] == {"name": dataLayerName, "required": "true", "type": "string"}
        assert metadata[dataLayerName+"_netCDFName"] == {"name": dataLayerName, "required": "false", "type": "string"}
        assert metadata[dataLayerName+"_units"] == {"name": dataLayerName, "required": "false", "type": "string"}
        assert metadata[dataLayerName+"_minBound"] == {"name": dataLayerName, "required": "false", "type": "float"}
        assert metadata[dataLayerName+"_maxBound"] == {"name": dataLayerName, "required": "true", "type": "float"}
        assert metadata[dataLayerName+"_standardName"] == {"name": dataLayerName, "required": "false", "type": "string"}
        assert metadata[dataLayerName+"_longName"] == {"name": dataLayerName, "required": "false", "type": "string"}
        assert metadata[dataLayerName+"_temporalChunking"] == {"name": dataLayerName, "required": "false", "type": "integer"}
        assert metadata[dataLayerName+"_temporalSkipInterval"] == {"name": dataLayerName, "required": "false", "type": "integer"}
        assert metadata[dataLayerName+"_timeDimensionName"] == {"name": dataLayerName, "required": "false", "type": "string"}




def test_verify_config_variables():
    """
    Test that loaded configurations are verified and converted to the appropriate data types:
    1) Missing variables added if they have default values
    2) Values are converted to their correct types
    3) Correct type inference for custom variables (for which there is no metadata)
    4) Exception thrown if a DataLayer doesn't have all of the required variables to define it
    """
    configPath = path.join(fe_setup_tools.get_fluxengine_root(), "test_data", "test_config_valid_full.conf")
    baseConfig = fe_setup_tools.read_config_file(configPath)
    metadata = fe_setup_tools.read_config_metadata(path.join(fe_setup_tools.get_fluxengine_root(), "core", "settings.xml"))
    
    #1) Missing variables added if they have default values
    config = baseConfig.copy()
    del config["bias_k"]
    fe_setup_tools.verify_config_variables(config, metadata)
    assert config["bias_k"] == False
    
    #Values are converted to their correct types
    assert type(config["output_structure"]) == str
    assert type(config["bias_k_value"]) == float
    assert type(config["output_temporal_chunking"]) == int
    
    #3) Correct type inference for custom variables (for which there is no metadata)
    assert type(config["custom_expected_string"]) == str
    assert type(config["custom_expected_float"]) == float
    assert type(config["custom_expected_int"]) == int
    assert type(config["custom_expected_bool"]) == bool
    
    #4) Exception thrown if a DataLayer doesn't have all of the required variables to define it
    config = baseConfig.copy()
    del config["sstskin_prod"]
    with pytest.raises(ValueError):
        fe_setup_tools.verify_config_variables(config, metadata)
    


def test_substitute_tokens():
    """
    Tests that date and filepath tokens are correctly substituted into the input string
    """
    
    #Test date and time related substitutions
    curTimeMock = datetime.datetime(1961, 4, 12, 6, 7)
    expectedOutput = "1961 61 04 4 APR Apr apr 12 102 06 07"
    assert fe_setup_tools.substitute_tokens("<YYYY> <YY> <MM> <M> <MMM> <Mmm> <mmm> <DD> <DDD> <hh> <mm>", curTimeMock) == expectedOutput
    
    #Test leap year cumulative days in a year figure is correct
    curTimeMock = datetime.datetime(2000, 3, 1, 0, 0) #leap year
    assert fe_setup_tools.substitute_tokens("<DDD>", curTimeMock) == "061"
    curTimeMock = datetime.datetime(2001, 3, 1, 0, 0) #non-leap year
    assert fe_setup_tools.substitute_tokens("<DDD>", curTimeMock) == "060"
    
    #Test file path substitution of root directory
    assert fe_setup_tools.substitute_tokens("<FEROOT>", curTimeMock) == fe_setup_tools.get_fluxengine_root()
    


def test_build_k_functor():
    """
    The correct gas transfer velocity parameterisation function can be retrieved from configuration files
    1) Build-in functor is retrieved from string config specification (without arguments)
    2) Build-in functor is retrieved from string config specification (with arguments)
    3) Custom functor specified by a different python file are correctly retrieved (without arguments)
    4) Custom functor specified by a different python file are correctly retrieved (with arguments)
    """
    import fluxengine.core.rate_parameterisation
    
    
    
    #1) Build-in functor is retrieved from string config specification (without arguments))
    runParameters = {"k_parameterisation": "k_Ho2006"};
    kFunctor = fe_setup_tools.build_k_functor(runParameters)
    assert isinstance(kFunctor, fluxengine.core.rate_parameterisation.k_Ho2006)
    #Test another parameterisation
    runParameters = {"k_parameterisation": "k_Nightingale2000"};
    kFunctor = fe_setup_tools.build_k_functor(runParameters)
    assert isinstance(kFunctor, fluxengine.core.rate_parameterisation.k_Nightingale2000)
    
    #2) Build-in functor is retrieved from string config specification (with arguments)
    runParameters = {"k_parameterisation": "kt_OceanFluxGHG", "kb_weighting": 1.5, "kd_weighting": 0.5};
    kFunctor = fe_setup_tools.build_k_functor(runParameters)
    assert isinstance(kFunctor, fluxengine.core.rate_parameterisation.kt_OceanFluxGHG)
    
    #3) Custom functor specified by a different python file are correctly retrieved (without arguments))
    runParameters = {"k_parameterisation": "example_custom_k_parameter_without_args"};
    kFunctor = fe_setup_tools.build_k_functor(runParameters, customGTVPath=path.join(fe_setup_tools.get_fluxengine_root(), "test_data", "test_custom_k_parameterisation.py"))
    assert str(type(kFunctor)) == "<class 'example_custom_k_parameter_without_args'>"
    
    #4) Custom functor specified by a different python file are correctly retrieved (with arguments)
    runParameters = {"k_parameterisation": "example_custom_k_parameter_with_args", "customConfigVar": 101.5};
    kFunctor = fe_setup_tools.build_k_functor(runParameters, customGTVPath=path.join(fe_setup_tools.get_fluxengine_root(), "test_data", "test_custom_k_parameterisation.py"))
    assert str(type(kFunctor)) == "<class 'example_custom_k_parameter_with_args'>"
    assert kFunctor.customConfigVar == 101.5



def test_get_preprocessing_funcs():
    """
    Preprocessing functions are correctly retrieved using their string names
    1) empty list or comma-only when no preprocessing functions are specified
    2) correct ordering when multiple preprocessing functions are specified
    3) whitespace is ignored when specifying functions
    4) duplicates are allowed (and included in duplicate)
    5) trailing commas allowed (and ignored)
    6) exception thrown when non-existing functions are specified
    """
    
    import fluxengine.core.data_preprocessing as dpp
    
    #1) empty list when no preprocessing functions are specified
    funcs = fe_setup_tools.get_preprocessing_funcs("")
    assert isinstance(funcs, list)
    assert len(funcs) == 0
    funcs = fe_setup_tools.get_preprocessing_funcs(",")
    assert isinstance(funcs, list)
    assert len(funcs) == 0
    
    #2) correct ordering when multiple preprocessing functions are specified
    funcs = fe_setup_tools.get_preprocessing_funcs("kelvin_to_celsius,pascal_to_millibar,nano_to_micro")
    assert funcs == [dpp.kelvin_to_celsius, dpp.pascal_to_millibar, dpp.nano_to_micro]
    
    #3) whitespace is ignored when specifying functions
    funcs = fe_setup_tools.get_preprocessing_funcs("\t\tkelvin_to_celsius,  pascal_to_millibar ,nano_to_micro, flip_latitude ,pow2  ")
    assert funcs == [dpp.kelvin_to_celsius, dpp.pascal_to_millibar, dpp.nano_to_micro, dpp.flip_latitude, dpp.pow2]
    
    #4) duplicates are allowed (and included in duplicate)
    funcs = fe_setup_tools.get_preprocessing_funcs("kelvin_to_celsius, pow2, pow2")
    assert funcs == [dpp.kelvin_to_celsius, dpp.pow2, dpp.pow2]
    
    #5) trailing commas allowed (and ignored)
    funcs = fe_setup_tools.get_preprocessing_funcs("kelvin_to_celsius, pow2, pow3,")
    assert funcs == [dpp.kelvin_to_celsius, dpp.pow2, dpp.pow3]
    
    #6) exception thrown when non-existing functions are specified
    with pytest.raises(ValueError):
        funcs = fe_setup_tools.get_preprocessing_funcs("kelvin_to_celsius, made_up_function")
        assert funcs == [dpp.kelvin_to_celsius, dpp.pow2, dpp.pow3]
    with pytest.raises(ValueError):
        funcs = fe_setup_tools.get_preprocessing_funcs("made_up_function")
        assert funcs == [dpp.kelvin_to_celsius, dpp.pow2, dpp.pow3]
    

def test_create_run_parameters():
    """
    Converts '_infile' variables for data layers with correct datetime substitutions
    """
    configVars = fe_setup_tools.read_config_file(path.join(fe_setup_tools.get_fluxengine_root(), "test_data", "test_config_valid_partial.conf"))
    metadata = fe_setup_tools.read_config_metadata(path.join(fe_setup_tools.get_fluxengine_root(), "core", "settings.xml"))
    fe_setup_tools.verify_config_variables(configVars, metadata)
    
    timePoint = datetime.datetime(2010, 1, 1, 0, 0, 0)
    runParams = fe_setup_tools.create_run_parameters(configVars, metadata, timePoint, 0, None, None , False, None)
    assert "sstskin_infile" in runParams
    assert "20100101" in runParams["sstskin_infile"] #<YYYY><MM><DD> substition
    
    #Increment one month
    timePoint = datetime.datetime(2010, 2, 1, 0, 0, 0)
    runParams = fe_setup_tools.create_run_parameters(configVars, metadata, timePoint, 0, None, None , False, None)
    assert "20100201" in runParams["sstskin_infile"] #<YYYY><MM><DD> substition

    

def test_run_fluxengine_basic():
    """
    Minimal test demonstrating that fluxengine can be ran
    """
    
    import tempfile
    from netCDF4 import Dataset
    
    
    with tempfile.TemporaryDirectory() as tempOutputDir:
        configPath = path.abspath(path.join(fe_setup_tools.get_fluxengine_root(), "test_data", "test_config_valid_full.conf"))
        
        (status, feObj) = fe_setup_tools.run_fluxengine(configPath, 2010, 2010, singleRun=True,
                                                        outputDirOverride=tempOutputDir);
        
        #Script exited correctly without issue. Note: if exception occurred, test will automatically fail.
        assert status == 0
        
        #Output file exists
        expectedOutputPath = path.join(tempOutputDir, "2010", "01", "OceanFluxGHG-month01-jan-2010-v0.nc")
        assert path.exists(expectedOutputPath)
        
        #Ocean gas flux output variable exists
        dataset = Dataset(expectedOutputPath, "r")
        assert "OF" in dataset.variables.keys()

    
    
def test_run_fluxengine_custom_gtv():
    """
    Minimal test demonstrating that fluxengine can be ran using a custom gas transfer velocity parameterisation implementation
    """
    
    import tempfile
    from netCDF4 import Dataset
    
    with tempfile.TemporaryDirectory() as tempOutputDir:
        configPath = path.abspath(path.join(fe_setup_tools.get_fluxengine_root(), "test_data", "test_config_valid_custom_gtv.conf"))
        customGtvPath = path.abspath(path.join(fe_setup_tools.get_fluxengine_root(), "test_data", "example_custom_gtv.py"))
        
        (status, feObj) = fe_setup_tools.run_fluxengine(configPath, 2010, 2010, singleRun=True,
                                                        outputDirOverride=tempOutputDir,
                                                        customGTVPath = customGtvPath);
        
        #Script exited correctly without issue. Note: if exception occurred, test will automatically fail.
        assert status == 0
        
        #Output file exists
        expectedOutputPath = path.join(tempOutputDir, "2010", "01", "OceanFluxGHG-month01-jan-2010-v0.nc")
        assert path.exists(expectedOutputPath)
        
        #Ocean gas flux output variable exists
        dataset = Dataset(expectedOutputPath, "r")
        assert "OF" in dataset.variables.keys()
        
        #Output selected gas transfer velocity data exiss, and has the expected value
        assert "OK3" in dataset.variables.keys()
        assert "SC" in dataset.variables.keys() #'scskin' name in NetCDF output is 'SC'
        kData = dataset.variables["OK3"][:]
        scskinData = dataset.variables["SC"][:]
        #Checking 
        missingMaskK = kData == DataLayer.missing_value
        missingMaskScskin = scskinData == DataLayer.missing_value
        assert np.all(missingMaskK == missingMaskScskin) #missing values match
        assert np.all(kData[missingMaskK==False] == scskinData[missingMaskScskin==False]*2.5) #2.5 is the test GTV's parameter, and scskin*2.5 is the example nonsense k calculation    


#TODO: These functions are not (yet?) tested
# match_filenames
# fe_obj_from_run_parameters
# generate_datetime_points
# process_timestep



if __name__ == "__main__":
    pytest.main(["-v", "--pdb"])