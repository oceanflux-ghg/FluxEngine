#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  2 13:44:02 2018
Drives the reanalyse_socat tool and provides access via commandline.
@author: Tom Holding
"""

#Example options:
#Single small region test run, ascii output, no coastal file.
#   -input_dir ~/data/SOCAT_ascii/SOCATv6/ -input_files SOCATv6_Indian.tsv -sst_dir ~/data/ocean_flux_ftp/SST_reynolds_avhrr/ -sst_tail 01_OCF-SST-GLO-1M-100-REYNOLDS.nc -output_dir ~/Files/fluxengine_v3/FluxEngine/output/reanalysis_socat_output_gridded -socatversion 6 -usereynolds -startyr 2008 -endyr 2015 -keepduplicates -keeptempfiles -asciioutput
#   -input_dir ~/Desktop/PostDoc/Backup/data/SOCAT_ascii/SOCATv6/ -socat_files SOCATv6_Indian.tsv -sst_dir ~/Desktop/PostDoc/Backup/data/uptodate/reynolds_avhrr_only_monthly_calculated_tmh -sst_tail 01_OCF-SST-GLO-1M-100-REYNOLDS_TMH.nc -output_dir ~/Files/fluxengine_v3/FluxEngine/output/reanalysis_socat_output_gridded -socatversion 6 -usereynolds -startyr 2008 -endyr 2015 -keepduplicates -keeptempfiles -asciioutput

#Single global input file, gridded output, no coastal file:
#   -input_dir ~/data/SOCAT_ascii/SOCATv4/ -input_files SOCATv4.tsv -sst_dir ~/data/ocean_flux_ftp/SST_reynolds_avhrr/ -sst_tail 01_OCF-SST-GLO-1M-100-REYNOLDS.nc -output_dir ~/Files/fluxengine_v3/FluxEngine/output/reanalysis_socat_output_gridded -socatversion 4 -usereynolds -startyr 2009 -endyr 2009 -keepduplicates -keeptempfiles
#Multiple regions, gridded output, no coastal file:
#   -input_dir ~/data/SOCAT_ascii/SOCATv4/ -input_files SOCATv4_NorthAtlantic.tsv SOCATv4_Indian.tsv -regions NA IN -sst_dir ~/data/ocean_flux_ftp/SST_reynolds_avhrr/ -sst_tail 01_OCF-SST-GLO-1M-100-REYNOLDS.nc -output_dir ~/Files/fluxengine_v3/FluxEngine/output/reanalysis_socat_output_gridded -socatversion 4 -usereynolds -startyr 2009 -endyr 2009 -keepduplicates -keeptempfiles
#Multiple regions, gridded output, with coastal file:
#   -input_dir ~/data/SOCAT_ascii/SOCATv4/ -input_files SOCATv4_NorthAtlantic.tsv SOCATv4_Indian.tsv SOCATv4_Coastal.tsv -regions NA IN CO -sst_dir ~/data/ocean_flux_ftp/SST_reynolds_avhrr/ -sst_tail 01_OCF-SST-GLO-1M-100-REYNOLDS.nc -output_dir ~/Files/fluxengine_v3/FluxEngine/output/reanalysis_socat_output_gridded -socatversion 4 -usereynolds -startyr 2009 -endyr 2009 -keepduplicates -keeptempfiles -withcoastal CO
#Multiple regions, ascii output, no coastal file:
#   -input_dir ~/data/SOCAT_ascii/SOCATv4/ -input_files SOCATv4_NorthAtlantic.tsv SOCATv4_Indian.tsv -regions NA IN -sst_dir ~/data/ocean_flux_ftp/SST_reynolds_avhrr/ -sst_tail 01_OCF-SST-GLO-1M-100-REYNOLDS.nc -output_dir ~/Files/fluxengine_v3/FluxEngine/output/reanalysis_socat_output_gridded -socatversion 4 -usereynolds -startyr 2009 -endyr 2009 -keepduplicates -keeptempfiles -asciioutput
#
#Using SOCATv5:
#SOCATv5 ASCII:
#   -input_dir ~/data/SOCAT_ascii/SOCATv5/ -input_files SOCATv5.tsv -sst_dir ~/data/ocean_flux_ftp/SST_reynolds_avhrr/ -sst_tail 01_OCF-SST-GLO-1M-100-REYNOLDS.nc -output_dir ~/Files/fluxengine_v3/FluxEngine/output/reanalysis_socat_output_gridded -socatversion 4 -usereynolds -startyr 2009 -endyr 2009 -keepduplicates -keeptempfiles -asciioutput
#SOCATv5 gridded:
#   -input_dir ~/data/SOCAT_ascii/SOCATv5/ -input_files SOCATv5.tsv -sst_dir ~/data/ocean_flux_ftp/SST_reynolds_avhrr/ -sst_tail 01_OCF-SST-GLO-1M-100-REYNOLDS.nc -output_dir ~/Files/fluxengine_v3/FluxEngine/output/reanalysis_socat_output_gridded -socatversion 4 -usereynolds -startyr 2009 -endyr 2009 -keepduplicates -keeptempfiles
#
#Using SOCATv6:
#SOCATv6 ASCII:
#   -input_dir ~/data/SOCAT_ascii/SOCATv6/ -input_files SOCATv6.tsv -sst_dir ~/data/ocean_flux_ftp/SST_reynolds_avhrr/ -sst_tail 01_OCF-SST-GLO-1M-100-REYNOLDS.nc -output_dir ~/Files/fluxengine_v3/FluxEngine/output/reanalysis_socat_output_gridded -socatversion 4 -usereynolds -startyr 2009 -endyr 2009 -keepduplicates -keeptempfiles -asciioutput
#SOCATv6 gridded, 2009:
#   -input_dir ~/data/SOCAT_ascii/SOCATv6/ -input_files SOCATv6.tsv -sst_dir ~/data/ocean_flux_ftp/SST_reynolds_avhrr/ -sst_tail 01_OCF-SST-GLO-1M-100-REYNOLDS.nc -output_dir ~/Files/fluxengine_v3/FluxEngine/output/reanalysis_socat_output_gridded -socatversion 4 -usereynolds -startyr 2009 -endyr 2009 -keepduplicates -keeptempfiles
#SOCATv6 gridded, 2016:
#   -input_dir ~/data/SOCAT_ascii/SOCATv6/ -input_files SOCATv6.tsv -sst_dir ~/data/ocean_flux_ftp/SST_reynolds_avhrr/ -sst_tail 01_OCF-SST-GLO-1M-100-REYNOLDS.nc -output_dir ~/Files/fluxengine_v3/FluxEngine/output/reanalysis_socat_output_gridded -socatversion 4 -usereynolds -startyr 2016 -endyr 2016 -keepduplicates -keeptempfiles

#Whole global SOCATv6 reanalysis (gridded)
#   -input_dir ~/data/SOCAT_ascii/SOCATv6/ -input_files SOCATv6.tsv -sst_dir ~/data/ocean_flux_ftp/SST_reynolds_avhrr/ -sst_tail 01_OCF-SST-GLO-1M-100-REYNOLDS.nc -output_dir ~/Files/fluxengine_v3/FluxEngine/output/reanalysis_socat_output_gridded -socatversion 4 -usereynolds -startyr 1957 -endyr 2017 -keepduplicates -keeptempfiles
#Whole global SOCATv6 reanalysis (ASCII)
#   -input_dir ~/data/SOCAT_ascii/SOCATv6/ -input_files SOCATv6.tsv -sst_dir ~/data/ocean_flux_ftp/SST_reynolds_avhrr/ -sst_tail 01_OCF-SST-GLO-1M-100-REYNOLDS.nc -output_dir ~/Files/fluxengine_v3/FluxEngine/output/reanalysis_socat_output_gridded -socatversion 4 -usereynolds -startyr 1957 -endyr 2017 -keepduplicates -keeptempfiles -asciioutput

import inspect;
import argparse;
from os import path, makedirs;
import shutil; #move desired output from reanalyse_socat to the specified output folder and delete undesired output

import fluxengine.tools.reanalyse_fco2.reanalyse_socat_v2 as rs;


if __name__ == "__main__":
    function = inspect.stack()[0][1]+", main";
    temporaryOutputPath = path.join(path.dirname(__file__), "reanalyse_socat", "output");
    
    
    #Setup command line parser, and parse arguments.
    description = """Utility which allows reanalysis of in situ data to a consistent temperature and depth.
        (under continual development, use with care)
        
        AUTHORS:
        Python driver script by Tom Holding <t.m.holding@exeter.ac.uk> 2018
        Everything else is directly modified from the reanalyse_socat open source project (see the reanalyse_fco2 directory in the FE tools folder).
        """;
    clParser = argparse.ArgumentParser(description=description, epilog="Both this script and the FluxEngine are in continual development. Use with care.");
    
    argParsePathGroup = clParser.add_argument_group(title='Data paths');
    argParsePathGroup.add_argument("-input_dir", help="Path to the directory containing tsv formatted input data (for example SOCAT data).", default="");
    argParsePathGroup.add_argument("-input_files", type=str, nargs='*', help='A list of input file names (for example "SOCATv5.tsv")', default=None);
    argParsePathGroup.add_argument("-sst_dir", help="Path to the directory containing netCDF formatted sea sureface temperature (either Reynolds or AATSR).", default="");
    argParsePathGroup.add_argument("-sst_tail", help="The SST file name which follows yyyymm. Do not include the year or month digits.", default="");
    argParsePathGroup.add_argument("-output_dir", help="Path to the output directory", default="~/socat_reanalysis_output/");
    
    argParseTimeGroup = clParser.add_argument_group(title='Run period');
    argParseTimeGroup.add_argument("-startyr", type=int, help='The year to start importing SOCAT data for.', default=2010);
    argParseTimeGroup.add_argument("-endyr", type=int, help='The last year to importing SOCAT data for (includes this year).', default=2010);
    
    argParseSettingsGroup = clParser.add_argument_group(title='Optional settings');
    argParseSettingsGroup.add_argument("-regions", metavar='<keyword_list>', type=str, nargs='*', help='A list of region codes to process. These will be used to name output files and should correspond to the list of file names given in -socat_files. If no region code is specified it is assumed a single file containing global data is provided.', default=["GL"]);
    argParseSettingsGroup.add_argument('-withcoastal',dest='withcoastal',type=str,help="The region code (defined using -regions) which corresponds to coastal data. Coastal data is appended to other regions where gridcells overlap and any remaining data is analysed seperately. Do not specify if no coastal data is used. Default is None (no coastal data used).",default=None)
    argParseSettingsGroup.add_argument('-asciioutput', dest='asciioutput', action='store_true', help="To output data as ascii lists rather than gridded netcdf.", default=False);
    argParseSettingsGroup.add_argument('-useaatsr', action='store_true', help="To use the AATSR SST data.", default=False);
    argParseSettingsGroup.add_argument('-usereynolds', action='store_true', help="To use the Reynolds SST data.", default=False);
    argParseSettingsGroup.add_argument('-keepduplicates', action='store_true', help="Do not detect and remove duplicate data points.", default=False);
    argParseSettingsGroup.add_argument('-keeptempfiles', action='store_true', help="Do not delete the temporary (full output) files produced by the reanalyse_socat.", default=False);
    
    argParseRunModeGroup = clParser.add_argument_group(title='Mode (socat vs non-socat)');
    argParseRunModeGroup.add_argument('-socatversion', type=int, dest='socatversion', help="If running in SOCAT mode this determines the version of the SOCAT data files to read",default=None);
    argParseRunModeGroup.add_argument('-notsocatformat', action='store_true', help="Set if running the tool on data which is not formatted exactly as socat files are. This allows custom column names/values to be set using the other options in this group (see below). Note that no filtering based on expedition code or quality control flags will be conducted in this run mode.", default=False);
    argParseRunModeGroup.add_argument('-year_col', type=str, help="Column name or number (indexed from 0) of the column which specifies year.", default=None);
    argParseRunModeGroup.add_argument('-month_col', type=str, help="Column name or number (indexed from 0) of the column which specifies month.", default=None);
    argParseRunModeGroup.add_argument('-day_col', type=str, help="Column name or number (indexed from 0) of the column which specifies day of the month.", default=None);
    argParseRunModeGroup.add_argument('-hour_col', type=str, help="Column name or number (indexed from 0) of the column which specifies hours.", default=None);
    argParseRunModeGroup.add_argument('-minute_col', type=str, help="Column name or number (indexed from 0) of the column which specifies minutes.", default=None);
    argParseRunModeGroup.add_argument('-second_col', type=str, help="Column name or number (indexed from 0) of the column which specifies seconds.", default=None);
    argParseRunModeGroup.add_argument('-longitude_col', type=str, help="Column name or number (indexed from 0) of the column which specifies longitude [0, 360).", default=None);
    argParseRunModeGroup.add_argument('-latitude_col', type=str, help="Column name or number (indexed from 0) of the column which specifies latitude (-90, 90).", default=None);
    argParseRunModeGroup.add_argument('-salinity_col', type=str, help="Column name or number (indexed from 0) of the column which specifies salinity.", default=None);
    argParseRunModeGroup.add_argument('-salinity_sub_col', type=str, help="Column name or number (indexed from 0) of the column which specifies values which will be used to substitute missing salinity values (e.g. modelled salinity field). If not specified then the value 35 will be substituted.", default=None);
    argParseRunModeGroup.add_argument('-SST_C_col', type=str, help="Column name or number (indexed from 0) of the column which specifies sea surface temperature (sub-skin) in Censius.", default=None);
    argParseRunModeGroup.add_argument('-Tequ_col', type=str, help="Column name or number (indexed from 0) of the column which specifies water temperature at equlibrium in Celsius.", default=None);
    argParseRunModeGroup.add_argument('-air_pressure_col', type=str, help="Column name or number (indexed from 0) of the column which specifies air pressure.", default=None);
    argParseRunModeGroup.add_argument('-air_pressure_sub_col', type=str, help="Column name or number (indexed from 0) of the column which specifies values which will be used to substitue missing air pressure values (e.g. modelled air pressure).", default=None);
    argParseRunModeGroup.add_argument('-air_pressure_equ_col', type=str, help="Column name or number (indexed from 0) of the column which specifies air pressure at equilibrium.", default=None);
    argParseRunModeGroup.add_argument('-fCO2_col', type=str, help="Column name or number (indexed from 0) of the column which specifies fugacity of carbon dioxide.", default=None);
    argParseRunModeGroup.add_argument('-expocode_col', type=str, help="Column name or number (indexed from 0) of the column which specifies a cruise or dataset specific code (e.g. 'Expocode' in SOCAT datasets).", default=None);
    
    
    
    
    clArgs = clParser.parse_args();
    
    if clArgs.input_dir == None or clArgs.input_files == None or clArgs.sst_dir == None or clArgs.sst_tail == None:
        raise SystemExit("Error: -input_dir, -input_files, -sst_dir and -sst_tail are all mandatory arguments. At least one of these has not been specified.");
    
    #Check no temporary output was left in place from a previous run (if it is, delete it)
    if path.exists(temporaryOutputPath) == True:
        shutil.rmtree(temporaryOutputPath);
    
    
    exitCode = rs.RunReanalyseSocat(socatdir=clArgs.input_dir,
                         socatfiles=clArgs.input_files,
                         sstdir=clArgs.sst_dir,
                         ssttail=clArgs.sst_tail,
                         regions=clArgs.regions,
                         usereynolds=clArgs.usereynolds,
                         startyr=clArgs.startyr,
                         endyr=clArgs.endyr,
                         asciioutput=clArgs.asciioutput,
                         keepduplicates=clArgs.keepduplicates,
                         
                         socatversion=clArgs.socatversion,
                         notsocatformat=clArgs.notsocatformat,
                         year_col=clArgs.year_col,
                         month_col=clArgs.month_col,
                         day_col=clArgs.day_col,
                         hour_col=clArgs.hour_col,
                         minute_col=clArgs.minute_col,
                         second_col=clArgs.second_col,
                         longitude_col=clArgs.longitude_col,
                         latitude_col=clArgs.latitude_col,
                         salinity_col=clArgs.salinity_col,
                         salinity_sub_col=clArgs.salinity_sub_col,
                         SST_C_col=clArgs.SST_C_col,
                         Tequ_col=clArgs.Tequ_col,
                         air_pressure_col=clArgs.air_pressure_col,
                         air_pressure_sub_col=clArgs.air_pressure_sub_col,
                         air_pressure_equ_col=clArgs.air_pressure_equ_col,
                         fCO2_col=clArgs.fCO2_col,
                         expocode_col=clArgs.expocode_col,
                         
                         withcoastal=clArgs.withcoastal,
                         #output=clArgs.output_dir);
                         output=temporaryOutputPath);#"reanalyse_socat/output/");
    
    #copy output files to the output_dir folder.
    if exitCode == 0:
        print("Reanalysis completed successfully.")
        outputDir = path.abspath(path.expanduser(clArgs.output_dir));
        
        #Create output directory
        if path.exists(outputDir) == False:
            makedirs(outputDir);
        try:
            #copy relevant reanalyse_socat output files
            shutil.move(path.join(temporaryOutputPath), outputDir);
            #shutil.copytree(path.join(temporaryOutputPath), outputDir);
            
            #remove irrelevant reanalyse_socat output files
            if clArgs.keeptempfiles == False:
                try:
                    shutil.rmtree(temporaryOutputPath);
                except:
                    pass;
            else:
                print("-keeptempfiles was set so full output can be viewed at:", "Output can be viewed at:", path.abspath(temporaryOutputPath));
        except shutil.Error as e: #If the output directory already exists don't delete anything and tell the user.
            if "already exists" in e.args[0]:
                print(e.args[0]);
                print("Output can be viewed at:", path.join(temporaryOutputPath, "reanalysed_data"));
            else: #Otherwise, propagate the exception
                raise e;
        
        




















