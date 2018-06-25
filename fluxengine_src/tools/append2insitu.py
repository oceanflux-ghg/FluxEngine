#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Converts a netCDF output file produced by FluxEngine into a text file with one
column for each data layer and one row for each non-missing value grid cell.

Created on Fri May 25 11:04:03 2018

@author: tom holding
"""

import argparse;
from glob import glob;
from netCDF4 import Dataset;
import numpy as np;
#from datetime import datetime, timedelta;
import pandas as pd;


#Return the index of the nearest element in an array
#   val: Value to be matched
#   array: array to match against
def find_nearest(val, arr):
    return (np.abs(arr - val)).argmin()


def parse_cl_arguments():
    description = unicode("""Converts netCDF3 data produced as output from FluxEngine into text data (e.g. csv, tsv).
    The text output will contain one column for each data layer, and one row for each cell of the grid (if there is no missing data).
    One text file will be created for each netCDF3 file.
    """, 'utf-8');
    
    parser = argparse.ArgumentParser(description=description);
    parser.add_argument("feOutput",
                        help="path to a FluxEngine output file");
    parser.add_argument("insituData",
                        help="path to a plain text in situ data file");
    parser.add_argument("outPath",
                        help="Path to store the produced text file.");
    parser.add_argument("--varsToAppend", nargs="+", default=["OF", "OK1"],
                        help="List of netCDF variable names which will be appended to the in situ data file. Default is the air-sea gas flux (OF) and the gas transfer velocity (K)");
    parser.add_argument("--delim", default="\t",
                        help="delimiter token used to seperate data in the in situ data file. Default is a tab.");
    parser.add_argument("--latCol", default="Latitude",
                        help="Name or column number (starting from 0) of the latitude column in the in situ input file. Default is a 'latitude'.");
    parser.add_argument("--lonCol", default="Longitude",
                        help="Name or column number (starting from 0) of the longitude column in the in situ input file. Default is 'longitude'.");
    parser.add_argument("-c", "--commentChar", default="#",
                        help="Character prefix which indicates a comment. Default is a '#'. Set to an empty string to turn comments off.");
    parser.add_argument("-m", "--missingValue", default="nan",
                        help="Value used to indicate a missing value in the input text file. Default is 'nan'.");
    parser.add_argument("--encoding", default="utf-8",
                        help="Encoding of the input text file. Default it 'utf-8'");
    parser.add_argument("-d", "--dateIndex", type=int, default=0,
        help="The column number (starting from 0) which contains the date/time field in the in situ input data file. Default is a 0.");
    parser.add_argument("-n", "--numCommentLines", nargs="+", type=int, default=[],
        help="Number of comment lines before the header in your in situ input data file. These lines will be ignored). Default is a 0.");
    clArgs = parser.parse_args();
    
    return clArgs



if __name__ == "__main__":
    #Parse commandline arguments
    clArgs = parse_cl_arguments();
    feOutputPath = clArgs.feOutput;
    insituDataPath = clArgs.insituData;
    outputPath = clArgs.outPath;
    varsToAppend = clArgs.varsToAppend;
    delim = clArgs.delim;
    dateIndex = clArgs.dateIndex;
    rowsToSkip = clArgs.numCommentLines;
    missingValue = clArgs.missingValue;
    encoding = clArgs.encoding;
    latCol = clArgs.latCol;
    lonCol = clArgs.lonCol;
    
    #Keep track of failures.
    #failedFiles = []; #Files which could not be opened / processed.
    
    #Process each input netCDF file and produce a text file
    print "Combining files at ", feOutputPath, "and", insituDataPath;
    
    #Read in situ data file
    insituData = pd.read_table(insituDataPath, sep=delim, skiprows=rowsToSkip, parse_dates=[dateIndex], encoding=encoding);
    
    #Read feOutput file
    try:
        ncFile = Dataset(feOutputPath, 'r');
    except IOError:
        print "Failed to open", feOutputPath;
    
    #Create empty vectors for each netCDF variable
    newVectors = {};
    for variable in varsToAppend:
        newVectors[variable] = np.zeros([len(insituData)]);
    
    #Store lat and lon arrays
    latArr =  ncFile.variables["latitude"][:];
    lonArr = ncFile.variables["longitude"][:];
    
    #Iterate through each line in the in situ data file and append FluxEngine output variables to it.
    for i, row in insituData.iterrows():
        print i;
        lat = row[latCol];
        lon = row[lonCol];
        x = find_nearest(lat, latArr);
        y = find_nearest(lon, lonArr);
        
        for variable in varsToAppend:
            value = ncFile.variables[variable][0,x,y];
            if np.ma.is_masked(value):
                newVectors[variable][i] = missingValue;
                print "Warning: No data found for lat", lat, "lon", lon;
            else:
                newVectors[variable][i] = float(value);
    
    #Add new vectors to dataframe and export
    for variable in varsToAppend:
        insituData[variable] = newVectors[variable];
    
    insituData.to_csv(outputPath, sep=delim, encoding=encoding, index=False);










