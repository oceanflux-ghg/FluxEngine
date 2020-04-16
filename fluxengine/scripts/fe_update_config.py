#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 09:12:03 2020

Takes a configuration file designed for FluxEngine v3.x and writes a copy of
the configuration file updated for FluxEngine 4.0

@author: Tom H
"""

import argparse;
from sys import exit;
import fluxengine;
import fluxengine.core.fe_setup_tools as setup;


#Create a string specifying a valid configuration file version tag
def generate_version_tag(version):
    return "#?FluxEngineConfigVersion:"+str(version)+"\n";


if __name__ == "__main__":
    description = """Utility for updating FluxEngine configuration files to be 
                     compatible with FluxEngine v4.0.""";
    
    clParser = argparse.ArgumentParser(description=description, epilog="Both this script and the FluxEngine are in continual development. Use with care.");
    clParser.add_argument("inputPath", help="path to a FluxEngine configuration file to be updated");
    clParser.add_argument("outputPath", help="path to write the updated FluxEngine configuration file to");
    clParser.add_argument("-alwayswrite", help="Forces the output file to be written even if no changes were made to the input configuration file.", action="store_true", default=False);
    clParser.add_argument("-oldversion", help="FluxEngine version that the old configuration file was designed for (e.g. 3.1). This is only used if the version cannot be determined from the configuration file itself.", type=float, default=None);
    clArgs = clParser.parse_args();
    
    #Read input
    file = open(clArgs.inputPath, 'r');
    config = file.readlines();
    file.close();
    
    #extract config version
    try:
        inputVersion = setup.parse_config_version_tag(config[0]);
    except ValueError as e:
        print(e);
        if clArgs.oldversion is not None:
            inputVersion = clArgs.oldversion;
        else:
            print("Error extracting input configuration file version. Either add this information to the old configuration file or specify it using the -oldversion flag.");
            exit(0);
    
    #Check that the input file is actually an old version
    if inputVersion == fluxengine.__version__:
        print("Given configuration file is already conformant to the current version.");
        if clArgs.alwayswrite == False:
            print("No output file has been written.");
            exit(0);
    
    
    #Update 'co2' config variables to  'gas' config variables
    for i in range(len(config)):
        line = config[i];
        if "vco2_air" in line:
            config[i] = line.replace("vco2_air", "vgas_air");
        if "pco2_air" in line:
            config[i] = line.replace("pco2_air", "pgas_air");
        if "pco2_sw" in line:
            config[i] = line.replace("pco2_sw", "pgas_sw");
    
    #set new config version
    if inputVersion >= 4.0:
        config[0] = generate_version_tag(4.0);
    else: #In the case of v3 and earlier, there is no config version tag
        config.insert(0, generate_version_tag(4.0));
    
    #write to file
    outputFile = open(clArgs.outputPath, 'w');
    outputFile.writelines(config);
    outputFile.close();
    print("Updated config file written to: "+clArgs.outputPath); 