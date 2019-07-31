#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Converts a netCDF output file produced by FluxEngine into a text file with one
column for each data layer and one row for each non-missing value grid cell.

Created on Fri May 25 11:04:03 2018

@author: tom holding
"""

import pathlib;
from netCDF4 import Dataset;
import numpy as np;
import pandas as pd;

#Takes a list of file globs and returns a list of every file that matches any of the file globs.
#If a file matches more than one glob it will be included multiple times.
#   inFiles: List of file globs (uses standard Unix match patterns)
#def match_input_file_globs(inFiles):
#    expandedGlobs = [];
#    for inFile in inFiles:
#        expandedGlobs += glob(inFile);
#    return expandedGlobs;
def match_input_file_globs(inFiles):
    expandedGlobs = [];
    for inFile in inFiles:
        #pathlib does not support globbing with absolute paths, must make a relative path from the path anchor
        pathAnchor = pathlib.Path(pathlib.Path(inFile).anchor);
        relativePath = pathlib.Path(inFile).relative_to(pathAnchor);
        expandedGlobs += pathAnchor.glob(str(relativePath)); #Find files matching glob
        
    expandedGlobs = [str(p) for p in expandedGlobs]; #Convert from pathlib.Path to string
    if len(expandedGlobs) != 0:
        return expandedGlobs;
    else:
        return inFiles;

def convert_ncdf_to_text(inFiles, textOutPath, delim="\t", commentChar="#", missingValue="nan", encoding="uft-8"):
    #Expand and match any globs to generate the full list of input files.
    inFiles = match_input_file_globs(inFiles);
    
    #Keep track of failures.
    failedFiles = []; #Files which could not be opened / processed.
    
    #Process each input netCDF file and produce a text file
    for inFile in inFiles:
        print("Processing file at:", inFile);
        try:
            ncFile = Dataset(inFile, 'r');
        except IOError:
            print("Failed to open", inFile);
            failedFiles.append(inFile);
            continue; #Move to the next file
            
        varNames = [varName for varName in list(ncFile.variables.keys()) if varName not in ["time", "latitude", "longitude"]];
        #varNames = ["OF"];
        unitNames = [ncFile.variables[var].units for var in varNames];
        colNames = [varNames[i]+" ["+unitNames[i]+"]" if unitNames[i] != "" else varNames[i] for i in range(0, len(varNames))];#  else varNames[i]];
        
        #store the lon and lat data
        latArray = ncFile.variables["latitude"][:];
        lonArray = ncFile.variables["longitude"][:];
        
        #construct a row for each grid cell:
        rowList = []; #Store rows in a list and just perform a single append to avoid repeatedly reallocating memory for a growing dataframe
        for ilat, lat in enumerate(latArray):
            for ilon, lon in enumerate(lonArray):
                #For each grid cell with a gas flux value, create a row.
                if np.ma.is_masked(ncFile.variables["OF"][0,lat,lon]) == False:
                    row = np.empty([len(varNames)+2]);
                    row[0] = lat;
                    row[1] = lon;
                    for i, varName in enumerate(varNames):
                        idx = i+2; #+2 for lat and lon
                        val = ncFile.variables[varNames[i]][0,lat,lon];
                        if np.ma.is_masked(val) == False:
                            row[idx] = val;
                        else:
                            row[idx] = np.nan;
                    
                    #Add row to the row list
                    rowList.append(row);
        
        #Create dataframe and write as tsv
        df = pd.DataFrame(rowList, columns = ["Latitude", "Longitude"]+colNames);
        outFile = inFile.rsplit('.',1)[0] + ".tsv";
        df.to_csv(outFile, sep="\t", index=False, encoding='utf-8');
        print("outputFile written to:", outFile);
