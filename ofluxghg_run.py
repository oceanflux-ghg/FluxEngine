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


from os import getcwd, path, makedirs;
import argparse; #parsing command line arguments
import socket; #for hostname
import time;
import calendar; #string abbreviations of months
import inspect; #stack

import fluxengine_src.fe_setup_tools as setup;


#Takes a config file, a list of years and months, and runs the flux engine for each month/year combination.
def run_fluxengine(configFilePath, yearsToRun, monthsToRun, verbose=False, processLayersOff=True,
                   takahashiDriver=False, pco2DirOverride=None, outputDirOverride=None):
    function = inspect.stack()[0][1]+", "+inspect.stack()[0][3];
    hostname = socket.gethostname();
    workingDirectory = getcwd();
    if verbose:
        print "Hostname identified as: ", hostname;
        print "Working directory is: ", workingDirectory;
    
    global fe; #For debugging purposes this makes it easy to query the FluxEngine object after execution has finished.
    
    
    #Parse config file
    configPath = path.join(workingDirectory, configFilePath);
    configVariables = setup.read_config_file(configPath, verbose=verbose);
    
    #Parse settings file for default metadata about the config variables
    settingsPath = path.join(workingDirectory, configVariables["src_home"], "settings.xml");
    metadata = setup.read_config_metadata(settingsPath, verbose=verbose);
    
    
    #Substitute commandline override arguments (-pco2_dir_override, -output_dir_override)
    if (pco2DirOverride != None):
        print "Using optional override for pCO2w data directory. '%s' will be set to '%s' (ie overriding both directory of pco2 data and selection in the configuration file)." % (configVariables["pco2"], clArgs.pco2_dir_override);
        configVariables["pco2_path"] = path.abspath(pco2DirOverride);
    if (outputDirOverride != None):
        configVariables["output_dir"] = path.abspath(outputDirOverride);
        print "Using optional override for output directory, output_dir will be set to % (overriding the configuration file value)." % (outputDirOverride);
    
    #Printing some feedback...
    if processLayersOff == True and verbose:
        print "Switching off generation of processing indicator layers. This reduces the processing time by appx. 50% (switch: process_layers_off).";
        #configVariables["process_indicator_layers"] = None;
    if takahashiDriver == True and verbose:
        print "This is a takahashi validation run. Ensure config is the configs/takahashi09_validation.conf file supplied with FluxEngine.";
    
    #Using the metadata, process the config variables appropriately
    #Checks types are valid, converts strings to the data types required by fluxengine
    #Checks data layers contain at least a prod and path, etc.
    setup.verify_config_variables(configVariables, metadata, verbose=verbose);
    
    ############
    ##TODO:
    #Alert user to any unused variables (can we infer all required variables yet?).    
    #Error if required variable hasn't been specified.
    ############
    
    
    #Begin main execution logic
    processTimeStr = time.strftime("%d/%m/%Y %H:%M:%S");
    print "Executing on '%s' at %s" % (hostname, processTimeStr);
    
    for year in yearsToRun:
        for monthNum in monthsToRun:
            #Run parameters can vary depending on the current month and year (e.g. paths and filenames,
            #So these must be generated on a per-month/year basis.
            try:
                runParameters = setup.create_run_parameters(configVariables, metadata, year, monthNum, processTimeStr, configFilePath, processLayersOff);
                
                #TODO: temporary stop-gap. Takahashi driver switch will be removed from future releases and moved to the configuration file.
                if takahashiDriver == True:
                    runParameters["TAKAHASHI_DRIVER"] = True;
                else:
                    runParameters["TAKAHASHI_DRIVER"] = False;
            except ValueError as e:
                print e.args;
                return;
            except OSError as e:
                print e.args;
                return;
            
            #Create output file path
            try:
                if path.exists(runParameters["output_dir"]) == False:
                    makedirs(runParameters["output_dir"]);
            except OSError as e:
                print "Couldn't create output directory '%s'. Do you have write access?" % runParameters["output_dir"];
                print type(e), e.args;
            
            #Create fluxengine object setup according to runParameters
            fe = setup.fe_obj_from_run_parameters(runParameters, metadata, processLayersOff, verbose=False);
            
            #Run fluxengine            
            if fe != None:
                #try:
                    returnCode = fe.run();
                #except Exception as e:
                #    print "\n\n%s: Exception caught while running FluxEngine:" % function;
                #    print type(e), e.args;
                #    return 1;
            
            #Check for successful run, if one fails don't run the rest.
            if returnCode != 0:
                print ("%s: There was an error running flux engine:\n\n"%function), e.args[0];
                print "Exiting...";
                return returnCode;
            else:
                print "Flux engine exited with exit code:", returnCode;
                print calendar.month_abbr[monthNum+1], year, "completed successfully.\n";
            
    #runStatus = {};
    #runStatus["return_code"] = returnCode;
    #runStatus["output_dir"] = runParameters["output_dir"];
    #runStatus["config_used"] = configFilePath;
    return 0;
                

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
    
    run_fluxengine(clArgs.config, yearsToRun, monthsToRun, verbose=clArgs.verbose,
                   processLayersOff=clArgs.process_layers_off,
                   takahashiDriver=clArgs.use_takahashi_driver,
                   pco2DirOverride=clArgs.pco2_dir_override,
                   outputDirOverride=clArgs.output_dir_override);


