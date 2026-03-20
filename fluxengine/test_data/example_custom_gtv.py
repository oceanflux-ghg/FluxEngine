#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 12:44:43 2026

@author: tom holding
"""

from fluxengine.core.rate_parameterisation import KCalculationBase
from fluxengine.core.datalayer import DataLayer


#A nonsense custom gas transfer velocity parameterisation, implemented in a separate file
#Used to test third-party execution of custom GTVs
class k_CustomGtvTest(KCalculationBase):
    # Optional initialiser arguments can be used, but their names must correspond to names in the config file.
    # For example some_parameter is used here.
    #Note: In the input arguments, it's best to append the parameterisation's
    # name or some other unique namespace to the arg names. This avoids polluting
    # the config/runParameter's namespace (minimises risk of name clashes)
    def __init__(self, k_CustomGtvTest_some_parameter):
        self.name = self.__class__.__name__
        self.someParameter = k_CustomGtvTest_some_parameter  # 'example_init_parameter' would need to be defined in the configuration file

    # Must return a list of strings corresponding to the input data layers required by the k calculations. These must already exist.
    def input_names(self):
        return ["scskin"]

    # Must return a list of strings corresponding to the names of data layers which it will write to. These may or may not already exist.
    # If they don't exist, they'll be created by FluxEngine as float DataLayer objects, and initialised to missing_value before the GTV is executed.
    def output_names(self):
        return ["k"]

    # Main k calculation. Input and output datalayers can be extracted from 'data'.
    # Should modify output layers in place (i.e. without copying), and return True or False to indicate successful execution.
    def __call__(self, data):
        print("Using", self.name, "with some_parameter =", self.someParameter)
        
        validMask = data["scskin"].fdata != DataLayer.missing_value
        
        #Example nonsense gas transfer velocity calculation
        data["k"].fdata[validMask] = data["scskin"].fdata[validMask] * self.someParameter

        #Return true, indicating everything went correctly
        return True