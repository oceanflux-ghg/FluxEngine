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
from datetime import datetime, timedelta;
import pandas as pd;


#Return the index of the nearest element in an array
#   val: Value to be matched
#   array: array to match against
def find_nearest(val, arr):
    return (np.abs(arr - val)).argmin()

#Return the time step this point is in
#   currentTime: time of the current val/sample in seconds from startTime
#   temporalResolution: time step between time points in seconds
#   startTime: time in seconds of the first time step
def find_timestep(currentTime, temporalResolution, startTime):
    dt = currentTime - startTime;
    index = int(np.floor(dt/temporalResolution));
    return index;


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
    timeArr = ncFile.variables["time"][:];
    startTime = timeArr[0];
    
    #If the temporal dimension is used in the FluxEngine output files
    if len(timeArr) > 1:
        #get temporal resolution
        dtStr = ncFile.temporal_resolution;
        dtDays = int(dtStr.split(" ")[0]);
        dtHours = int(dtStr.split(" ")[2].split(":")[0]);
        dtMinutes = int(dtStr.split(" ")[2].split(":")[1]);
        dtSeconds = int(dtStr.split(" ")[2].split(":")[2]);

        temporalResolutionSecs = dtSeconds + 60*dtMinutes + 60*60*dtHours + 60*60*24*dtDays;


    
    #Iterate through each line in the in situ data file and append FluxEngine output variables to it.
    for i, row in insituData.iterrows():
        #Match temporal dimension firect
        #calculate time in seconds since 1970-01-01 for this row
        timeSecs = (datetime(row[dateIndex].year, row[dateIndex].month, row[dateIndex].day, row[dateIndex].hour, row[dateIndex].minute, row[dateIndex].second) - datetime(1970, 1, 1)).total_seconds();

        if len(timeArr) > 1:
            t = find_timestep(timeSecs, temporalResolutionSecs, startTime); #indices of time step containing this row
        else:
            t = 0;

        #Match spatial location after time
        lat = row[latCol];
        lon = row[lonCol];
        x = find_nearest(lat, latArr); #indices of nearest lat value in the netCDF file
        y = find_nearest(lon, lonArr); #indices of nearest lon value in the netCDF file
        
        for variable in varsToAppend:
            value = ncFile.variables[variable][t,x,y];
            if np.ma.is_masked(value):
                newVectors[variable][i] = missingValue;
                if verbose:
                    print("Row", i, "Warning: No data found for timesetp", t, "lat", lat, "lon", str(lon)+".", missingValue, "inserted instead");
            else:
                newVectors[variable][i] = float(value);
    
    #Add new vectors to dataframe and export
    for variable in varsToAppend:
        colName = variable+" ["+ncFile.variables[variable].units+"]";
        insituData[colName] = newVectors[variable];
    
    insituData.to_csv(outputPath, sep=delim, encoding=encoding, index=False);
