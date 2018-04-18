#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 19 07:33:12 2018
Contains functions to fascilitate comparison of Datalayers between two netCDF files.

@author: tomholding
"""

from netCDF4 import Dataset;
import numpy as np;
import os;

#calculates a matrix which contains the difference between two 2D matrices.
#if missing_value is specified it will ignore these values (difference is set to 0) if either input matrix contain this value.
def calc_diff(old, new, missingValue=None):
    if old.shape != old.shape:
        raise ValueError("old and new shape do not match.");
    
    output = np.zeros(old.shape);
    
    #Ensure we're working with np.array objects, but don't change the originals.    
    old2 = np.array(old);
    new2 = np.array(new);
    
    for y in range(0, old2.shape[1]):
        for x in range(0, old2.shape[0]):
            if missingValue == None or (old2[x,y] != missingValue and new2[x,y] != missingValue):
                curDiff =  old2[x,y] - new2[x,y];
                output[x,y] = curDiff;

    return output;

#Converts nans to missingValue. Returns modified version and does not modify original matrix.
def nans_to_missing_value(data, missingValue=-999.0):
    fdata = np.ravel(data.copy());
    if np.ma.is_masked(fdata):
        fdata.unshare_mask();
    
    for i in range(len(fdata)):
        if np.isnan(fdata[i]):
            fdata[i] = missingValue;
    fdata.shape = data.shape;
    return fdata;

#Compares two FluxEngine output files, calculates the sum of absolute differences (old-new) for each
#grid point, as well as differences in the number of masked values.
#returns ...
def compare_fe_output_files(pathToOld, pathToNew, missingValue=-999.0, verbose=False):
    oldNetCDF = Dataset(pathToOld, 'r');
    newNetCDF = Dataset(pathToNew, 'r');
    varsToCompare = old.variables.keys();
    
    diffs = {};
    for varToCompare in varsToCompare:
        #Ignore "time", "latitude" and "longitude"
        if varToCompare in [u"time", u"latitude", u"longitude"]:
            continue;
    
    oldVar = oldNetCDF.variables[varToCompare][:];
    #oldVar = np.nan_to_num(oldVar);
    
    try:
        newVar = newNetCDF.variables[varToCompare][:];
        #newVar = np.nan_to_num(newVar);
        
        #remove time dimensions
        if len(oldVar.shape) == 3:
            oldVar = oldVar[0,:,:];
        if len(newVar.shape) == 3:
            newVar = newVar[0,:,:];
        #else:
        #    oldVar = oldVar[:,:];
        #    newVar = newVar[:,:];
            
        if np.any(np.isnan(oldVar)):
            oldVar = nans_to_missing_value(oldVar);
        if np.any(np.isnan(newVar)):
            newVar = nans_to_missing_value(newVar);
        
        if np.ma.is_masked(oldVar) and np.ma.is_masked(newVar):
            oldMasked = np.sum(oldVar.mask);
            newMasked = np.sum(newVar.mask);
            maskedDiff = oldMasked-newMasked;
        else:
            maskedDiff = "at least one matrix is not masked.";

        diff = calc_diff(oldVar, newVar, missingValue=missingValue);
        if verbose:
            print varToCompare+":", np.sum(np.abs(diff)), "   ", maskedDiff;
        diffs[varToCompare] = (varToCompare, np.sum(np.abs(diff)), maskedDiff);
        
    #KeyError if there was no new value when there was an old value
    except KeyError:
        print "no new value for: "+varToCompare;

    #Finally check there are no variables in newNetCDF which are not in the old netCDF.
    for v in newNetCDF.variables.keys():
        if v not in varsToCompare: #if variable in the new output but not in the old reference output
            print "no OLD value for: "+v;

if __name__ == "__main__":
    from sys import argv;
    if len(argv) != 3 or len(argv) != 4:
        print "Two or three arguments are required; <path to reference/old output> <path to new output>. %d arguments detected." % len(argv)-1;
    if len(argv) == 4:
        missingValue = float(argv[3]);
    else:
        missingValue = -999.0;
    
    comparison = compare_fe_output_files(argv[1], argv[2], missingValue);
    
    print "Comparison of old and new output files:"
    print "Variable name | sum(abs(old-new)) | difference in number of masked values"
    for element in comparison:
        print element[0]+":", element[1], "   ", element[2];
