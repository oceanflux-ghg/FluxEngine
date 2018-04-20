#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 11:32:08 2018

Allows FluxEngine to be run from the commandline.
Parses command line arguments, verifies configuration file and runs FluxEngine for each month/year

@author: tomholding
"""

#TODO: systematic method to support non-mandatory vs mandatory variables
#TODO: process indicator layers (and parameters) should be passed in config as k_parameterisation is.
#TODO: fully remove -t options in here and fe_core.py (TAKAHASHI_DRIVER)

#Useful cl args
#output/SOCATv4_WoolfRuns/sst_salinity_gradients-N00/sst_salinity_gradients-N00.conf -l -m1
#configs/socatv4_sst_salinity_gradients-N00.conf -l -m1
#configs/takahashi09_validation.conf -l -t -s 2000 -e 2000 -m1


import argparse; #parsing command line arguments
import inspect; #stack

import fluxengine_src.fe_setup_tools as setup;
                

if __name__ == "__main__":
    function = inspect.stack()[0][1]+", main";
    
    #Setup command line parser, and parse arguments.
    description = """Util for generating the ESA OceanFlux Greenhouse Gases global climatology and time series fluxes.
        (development, use with care)
        
        AUTHOR
        Originally by Jamie Shutler <jams@pml.ac.uk> 2012-2013
        Transcribed to python and updated by Tom Holding <t.m.holding@exeter.ac.uk> 2018
        """;
    clParser = argparse.ArgumentParser(description=description, epilog="Both this script and the FluxEngine are in continual development. Use with care.");
    clParser.add_argument("config", help="path to the config file");
    clParser.add_argument("-pco2_dir_override", "-p", metavar='\b', help="pCO2 directory (overrides config specification)");
    clParser.add_argument("-output_dir_override", "-o", metavar='\b', help="output directory (overrides config specification)");
    clParser.add_argument("-process_layers_off", "-l", help="turns process indicator layers off (reduces runtime)",
                        action="store_true");
    clParser.add_argument("-use_takahashi_driver", "-t", help="indicates that this is a Takahashi validation run. This should only be used in conjunction with the supplied takahashi09_validation/conf file.",
                          action="store_true", default=False);
    #clParser.add_argument("-use_takahashi_validation", "-t", help="DEPRECIATED: This is non-functional and should not be used. Instead use the takahashi validation configuration file provided in the 'configs' folder.",
    #                    action="store_true");
    clParser.add_argument("-verbose", "-v", help="verbose: increases the amount of information sent to stdout.",
                        action="store_true");
    clParser.add_argument("-year_start", "-s", metavar='\b', help="first year to evaluate (default=2010)", type=int, default=2010);
    clParser.add_argument("-year_end", "-e", metavar='\b', help="final year to evaluate (default=2010)", type=int, default=2010);
    #clParser.add_argument("-list_k, help="...", action="store_true", default=False);
    #clParser.add_argument("-list_preprocessing", help="...", action="store_true", default=False);
    clParser.add_argument("-m1", help="only run a single month (useful for testing).", action="store_true", default=False);
    clArgs = clParser.parse_args();
    
    
    #Months and years to run the model over
    yearsToRun = range(clArgs.year_start, clArgs.year_end+1);
    monthsToRun = range(0,12);
    if clArgs.m1 == True: #If m1 flag is set, just run the first month.
        monthsToRun = [monthsToRun[0]];
        yearsToRun = [yearsToRun[0]];
    
    returnCode, fe = setup.run_fluxengine(clArgs.config, yearsToRun, monthsToRun, verbose=clArgs.verbose,
                   processLayersOff=clArgs.process_layers_off,
                   takahashiDriver=clArgs.use_takahashi_driver,
                   pco2DirOverride=clArgs.pco2_dir_override,
                   outputDirOverride=clArgs.output_dir_override);


