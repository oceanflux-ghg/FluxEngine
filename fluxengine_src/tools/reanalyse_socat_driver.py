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


#Example options:
#Single global input file, gridded output, no coastal file:
#   -socat_dir ~/data/SOCAT_ascii/SOCATv4/ -socat_files SOCATv4.tsv -sst_dir ~/data/ocean_flux_ftp/SST_reynolds_avhrr/ -sst_tail 01_OCF-SST-GLO-1M-100-REYNOLDS.nc -output_dir ~/Files/fluxengine_v3/FluxEngine/output/reanalysis_socat_output_gridded -socatversion 4 -usereynolds -startyr 2009 -endyr 2009 -keepduplicates -keeptempfiles
#Multiple regions, gridded output, no coastal file:
#   -socat_dir ~/data/SOCAT_ascii/SOCATv4/ -socat_files SOCATv4_NorthAtlantic.tsv SOCATv4_Indian.tsv -regions NA IN -sst_dir ~/data/ocean_flux_ftp/SST_reynolds_avhrr/ -sst_tail 01_OCF-SST-GLO-1M-100-REYNOLDS.nc -output_dir ~/Files/fluxengine_v3/FluxEngine/output/reanalysis_socat_output_gridded -socatversion 4 -usereynolds -startyr 2009 -endyr 2009 -keepduplicates -keeptempfiles
#Multiple regions, gridded output, with coastal file:
#   -socat_dir ~/data/SOCAT_ascii/SOCATv4/ -socat_files SOCATv4_NorthAtlantic.tsv SOCATv4_Indian.tsv SOCATv4_Coastal.tsv -regions NA IN CO -sst_dir ~/data/ocean_flux_ftp/SST_reynolds_avhrr/ -sst_tail 01_OCF-SST-GLO-1M-100-REYNOLDS.nc -output_dir ~/Files/fluxengine_v3/FluxEngine/output/reanalysis_socat_output_gridded -socatversion 4 -usereynolds -startyr 2009 -endyr 2009 -keepduplicates -keeptempfiles -withcoastal CO
#Multiple regions, ascii output, no coastal file:
#   -socat_dir ~/data/SOCAT_ascii/SOCATv4/ -socat_files SOCATv4_NorthAtlantic.tsv SOCATv4_Indian.tsv -regions NA IN -sst_dir ~/data/ocean_flux_ftp/SST_reynolds_avhrr/ -sst_tail 01_OCF-SST-GLO-1M-100-REYNOLDS.nc -output_dir ~/Files/fluxengine_v3/FluxEngine/output/reanalysis_socat_output_gridded -socatversion 4 -usereynolds -startyr 2009 -endyr 2009 -keepduplicates -keeptempfiles -asciioutput
#
#Using SOCATv5:
#SOCATv5 ASCII:
#   -socat_dir ~/data/SOCAT_ascii/SOCATv5/ -socat_files SOCATv5.tsv -sst_dir ~/data/ocean_flux_ftp/SST_reynolds_avhrr/ -sst_tail 01_OCF-SST-GLO-1M-100-REYNOLDS.nc -output_dir ~/Files/fluxengine_v3/FluxEngine/output/reanalysis_socat_output_gridded -socatversion 4 -usereynolds -startyr 2009 -endyr 2009 -keepduplicates -keeptempfiles -asciioutput
#SOCATv5 gridded:
#   -socat_dir ~/data/SOCAT_ascii/SOCATv5/ -socat_files SOCATv5.tsv -sst_dir ~/data/ocean_flux_ftp/SST_reynolds_avhrr/ -sst_tail 01_OCF-SST-GLO-1M-100-REYNOLDS.nc -output_dir ~/Files/fluxengine_v3/FluxEngine/output/reanalysis_socat_output_gridded -socatversion 4 -usereynolds -startyr 2009 -endyr 2009 -keepduplicates -keeptempfiles
#
#Using SOCATv6:
#SOCATv6 ASCII:
#   -socat_dir ~/data/SOCAT_ascii/SOCATv6/ -socat_files SOCATv6.tsv -sst_dir ~/data/ocean_flux_ftp/SST_reynolds_avhrr/ -sst_tail 01_OCF-SST-GLO-1M-100-REYNOLDS.nc -output_dir ~/Files/fluxengine_v3/FluxEngine/output/reanalysis_socat_output_gridded -socatversion 4 -usereynolds -startyr 2009 -endyr 2009 -keepduplicates -keeptempfiles -asciioutput
#SOCATv6 gridded, 2009:
#   -socat_dir ~/data/SOCAT_ascii/SOCATv6/ -socat_files SOCATv6.tsv -sst_dir ~/data/ocean_flux_ftp/SST_reynolds_avhrr/ -sst_tail 01_OCF-SST-GLO-1M-100-REYNOLDS.nc -output_dir ~/Files/fluxengine_v3/FluxEngine/output/reanalysis_socat_output_gridded -socatversion 4 -usereynolds -startyr 2009 -endyr 2009 -keepduplicates -keeptempfiles
#SOCATv6 gridded, 2016:
#   -socat_dir ~/data/SOCAT_ascii/SOCATv6/ -socat_files SOCATv6.tsv -sst_dir ~/data/ocean_flux_ftp/SST_reynolds_avhrr/ -sst_tail 01_OCF-SST-GLO-1M-100-REYNOLDS.nc -output_dir ~/Files/fluxengine_v3/FluxEngine/output/reanalysis_socat_output_gridded -socatversion 4 -usereynolds -startyr 2016 -endyr 2016 -keepduplicates -keeptempfiles

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
    argParsePathGroup.add_argument("-socat_files", type=str, nargs='*', help='A list of region codes to process (e.g. NA for North Atlantic). If no region codes are specified the global database is assumed.', default=None);
    argParsePathGroup.add_argument("-sst_dir", help="Path to the directory containing netCDF formatted sea sureface temperature (either Reynolds or AATSR).", default="");
    argParsePathGroup.add_argument("-sst_tail", help="The SST file name which follows yyyymm. Do not include the year or month digits.", default="");
    argParsePathGroup.add_argument("-output_dir", help="Path to the output directory", default="~/socat_reanalysis_output/");
    
    argParseTimeGroup = clParser.add_argument_group(title='Run period');
    argParseTimeGroup.add_argument("-startyr", type=int, help='The year to start importing SOCAT data for.', default=2010);
    argParseTimeGroup.add_argument("-endyr", type=int, help='The last year to importing SOCAT data for (includes this year).', default=2010);
    
    argParseSettingsGroup = clParser.add_argument_group(title='Optional settings');
    argParseSettingsGroup.add_argument("-regions", metavar='<keyword_list>', type=str, nargs='*', help='A list of region codes to process. These will be used to name output files and should correspond to the list of file names given in -socat_files. If no region code is specified it is assumed a single file containing global data is provided.', default=["GL"]);
    argParseSettingsGroup.add_argument('-withcoastal',dest='withcoastal',type=str,help="The region code (defined using -regions) which corresponds to coastal data. Coastal data is appended to other regions where gridcells overlap and any remaining data is analysed seperately. Do not specify if no coastal data is used. Default is None (no coastal data used).",default=None)
    argParseSettingsGroup.add_argument('-socatversion', type=int, dest='socatversion', help="The version of the SOCAT data files to read",default=2);
    argParseSettingsGroup.add_argument('-asciioutput', dest='asciioutput', action='store_true', help="To output data as ascii lists rather than gridded netcdf.", default=False);
    argParseSettingsGroup.add_argument('-useaatsr', action='store_true', help="To use the AATSR SST data.", default=False);
    argParseSettingsGroup.add_argument('-usereynolds', action='store_true', help="To use the Reynolds SST data.", default=False);
    argParseSettingsGroup.add_argument('-keepduplicates', action='store_true', help="Do not detect and remove duplicate data points.", default=False);
    argParseSettingsGroup.add_argument('-keeptempfiles', action='store_true', help="Do not delete the temporary (full output) files produced by the reanalyse_socat.", default=False);
    
    clArgs = clParser.parse_args();
    
    if clArgs.socat_dir == None or clArgs.socat_files == None or clArgs.sst_dir == None or clArgs.sst_tail == None:
        raise SystemExit("Error: -socat_dir, -socat_files, -sst_dir and -sst_tail are all mandatory arguments. At least one of these has not been specified.");
    
    #Check no temporary output was left in place from a previous run (if it is, delete it)
    if path.exists(temporaryOutputPath) == True:
        shutil.rmtree(temporaryOutputPath);
    
    exitCode = rs.RunReanalyseSocat(socatdir=clArgs.socat_dir,
                         socatfiles=clArgs.socat_files,
                         sstdir=clArgs.sst_dir,
                         ssttail=clArgs.sst_tail,
                         socatversion=clArgs.socatversion,
                         regions=clArgs.regions,
                         usereynolds=clArgs.usereynolds,
                         startyr=clArgs.startyr,
                         endyr=clArgs.endyr,
                         asciioutput=clArgs.asciioutput,
                         keepduplicates=clArgs.keepduplicates,
                         withcoastal=clArgs.withcoastal,
                         output="reanalyse_socat/output/");
        
    #copy output files to the output_dir folder.
    if exitCode == 0:
        print "Reanalysis completed successfully."
        outputDir = path.abspath(path.expanduser(clArgs.output_dir));
        
        #Create output directory
        #if path.exists(outputDir) == False:
        #    makedirs(outputDir);
        try:
            #copy relevant reanalyse_socat output files
            #shutil.move(path.join(temporaryOutputPath, "reanalysed_data"), outputDir);
            shutil.copytree(path.join(temporaryOutputPath, "reanalysed_data"), outputDir);
            
            #remove irrelevant reanalyse_socat output files
            if clArgs.keeptempfiles == False:
                shutil.rmtree(temporaryOutputPath);
            else:
                print "-keeptempfiles was set so full output can be viewed at:", "Output can be viewed at:", path.abspath(temporaryOutputPath);
        except shutil.Error as e: #If the output directory already exists don't delete anything and tell the user.
            if "already exists" in e.args[0]:
                print e.args[0];
                print "Output can be viewed at:", path.join(temporaryOutputPath, "reanalysed_data");
            else: #Otherwise, propagate the exception
                raise e;
        
        
        

shutil.copytree(temporaryOutputPath, path.expanduser("~/Desktop/test/directory/output/folder"))






















