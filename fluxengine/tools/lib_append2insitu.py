#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Converts a netCDF output file produced by FluxEngine into a text file with one
column for each data layer and one row for each non-missing value grid cell.

Created on Fri May 25 11:04:03 2018

@author: tom holding
"""

from netCDF4 import Dataset;
import numpy as np;
import pandas as pd;


#Return the index of the nearest element in an array
#   val: Value to be matched
#   array: array to match against
def find_nearest(val, arr):
    return (np.abs(arr - val)).argmin()


def append_to_in_situ(feOutputPath, insituDataPath, outputPath, varsToAppend = ["OF", "OK1"], delim="\t", latCol="Latitude", lonCol="Longitude", dateIndex=0, rowsToSkip=[], missingValue='nan', encoding='utf-8', verbose=False):
    
    #Keep track of failures.
    #failedFiles = []; #Files which could not be opened / processed.
    
    #Process each input netCDF file and produce a text file
    if verbose:
        print("Combining files at ", feOutputPath, "and", insituDataPath);
    
    #Read in situ data file
    insituData = pd.read_table(insituDataPath, sep=delim, skiprows=rowsToSkip, parse_dates=[dateIndex], encoding=encoding);
    
    #Read feOutput file
    ncFile = Dataset(feOutputPath, 'r');
    
    #Create empty vectors for each netCDF variable
    newVectors = {};
    for variable in varsToAppend:
        newVectors[variable] = np.zeros([len(insituData)]);
    
    #Store lat and lon arrays
    latArr =  ncFile.variables["latitude"][:];
    lonArr = ncFile.variables["longitude"][:];
    
    #Iterate through each line in the in situ data file and append FluxEngine output variables to it.
    for i, row in insituData.iterrows():
        lat = row[latCol];
        lon = row[lonCol];
        x = find_nearest(lat, latArr);
        y = find_nearest(lon, lonArr);
        
        for variable in varsToAppend:
            value = ncFile.variables[variable][0,x,y];
            if np.ma.is_masked(value):
                newVectors[variable][i] = missingValue;
                if verbose:
                    print("Row", i, "Warning: No data found for lat", lat, "lon", str(lon)+".", missingValue, "inserted instead");
            else:
                newVectors[variable][i] = float(value);
    
    #Add new vectors to dataframe and export
    for variable in varsToAppend:
        colName = variable+" ["+ncFile.variables[variable].units+"]";
        insituData[colName] = newVectors[variable];
    
    insituData.to_csv(outputPath, sep=delim, encoding=encoding, index=False);
