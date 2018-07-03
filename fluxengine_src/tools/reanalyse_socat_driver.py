#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  2 13:44:02 2018

@author: Tom Holding
"""

import inspect;
import argparse;
from os import path, makedirs;
import shutil; ##move desired output from reanalyse_socat to the specified output folder and delete undesired output

#TODO: 

import reanalyse_socat.reanalyse_socat_v2 as rs;


if __name__ == "__main__":
    function = inspect.stack()[0][1]+", main";
    temporaryOutputPath = path.join(path.dirname(__file__), "reanalyse_socat/output/");
    
    #Setup command line parser, and parse arguments.
    description = """Utility which allows reanalysis of SOCAT data to a consistent temperature and depth.
        (under continual development, use with care)
        
        This tool uses previously developed code by ???.
        
        AUTHOR
        Python driver script by Tom Holding <t.m.holding@exeter.ac.uk> 2018
        Everything else is from the reanalyse_socat open source project (see the reanalyse_socat directory in the FE tools folder).
        """;
    clParser = argparse.ArgumentParser(description=description, epilog="Both this script and the FluxEngine are in continual development. Use with care.");
    
    argParsePathGroup = clParser.add_argument_group(title='Data paths');
    argParsePathGroup.add_argument("-socat_dir", help="Path to the directory containing tsv formatted SOCAT data.", default="");
    argParsePathGroup.add_argument("-sst_dir", help="Path to the directory containing netCDF formatted sea sureface temperature (either Reynolds or AATSR).", default="");
    argParsePathGroup.add_argument("-sst_tail", help="The SST file name which follows yyyymm. Do not include the year or month digits.", default="");
    argParsePathGroup.add_argument("-output_dir", help="Path to the output directory", default="~/socat_reanalysis_output/");
    
    argParseTimeGroup = clParser.add_argument_group(title='Run period');
    argParseTimeGroup.add_argument("-startyr", type=int, help='The year to start importing SOCAT data for.', default=2010);
    argParseTimeGroup.add_argument("-endyr", type=int, help='The last year to importing SOCAT data for (includes this year).', default=2010);
    
    argParseSettingsGroup = clParser.add_argument_group(title='Optional settings');
    argParseSettingsGroup.add_argument("-regions", metavar='<keyword_list>', type=str, nargs='*', help='A list of region codes to process (e.g. NA for North Atlantic).', default=None);
    argParseSettingsGroup.add_argument('-socatversion', type=int, dest='socatversion', help="The version of the SOCAT data files to read",default=2);
    argParseSettingsGroup.add_argument('-asciioutput', dest='asciioutput', action='store_true', help="To output data as ascii lists rather than gridded netcdf.", default=True);
    argParseSettingsGroup.add_argument('-useaatsr', action='store_true', help="To use the AATSR SST data.", default=False);
    argParseSettingsGroup.add_argument('-usereynolds', action='store_true', help="To use the Reynolds SST data.", default=True);
    
    clArgs = clParser.parse_args();
    
    #Check no temporary output was left in place from a previous run (if it is, delete it)
    if path.exists(temporaryOutputPath) == True:
        shutil.rmtree(temporaryOutputPath);
    
    exitCode = rs.RunReanalyseSocat(socatdir=clArgs.socat_dir,
                         sstdir=clArgs.sst_dir,
                         ssttail=clArgs.sst_tail,
                         socatversion=clArgs.socatversion,
                         regions=clArgs.regions,
                         usereynolds=clArgs.usereynolds,
                         startyr=clArgs.startyr,
                         endyr=clArgs.endyr,
                         asciioutput=clArgs.asciioutput,
                         output="reanalyse_socat/output/");
        
    #copy output files to the output_dir folder.
    if exitCode == 0:
        print "SOCAT reanalysis completed successfully."
        outputDir = path.abspath(path.expanduser(clArgs.output_dir));
        
        #Create output directory
        if path.exists(outputDir) == False:
            makedirs(outputDir);
            
        try:
            #move relevant reanalyse_socat output files
            shutil.move(path.join(temporaryOutputPath, "reanalysed_socat"), outputDir);
            
            #remove irrelevant reanalyse_socat output files
            shutil.rmtree(temporaryOutputPath);
        except shutil.Error as e: #If the output directory already exists don't delete anything and tell the user.
            if "already exists" in e.args[0]:
                print e.args[0];
                print "Output can be viewed at:", path.join(temporaryOutputPath, "reanalysed_socat");
            else: #Otherwise, propagate the exception
                raise e;
        
        
        
























