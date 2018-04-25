#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 18 12:48:03 2018
Provides validation tools which can validate FluxEngine output against a known reference data set.

@author: tomholding
"""

#Use SET_PATHS to ensure that we are running the development version of the FluxEngine and not any previously installed version.
#import SET_PATHS;
#rootPath, srcPath = SET_PATHS.set_paths(verbose=True);


from fluxengine_src.fe_setup_tools import run_fluxengine, read_config_file;
from fluxengine_src.tools.ofluxghg_flux_budgets import run_flux_budgets;
from fluxengine_src.tools.compare_net_budgets import calc_net_budget_percentages;
from argparse import Namespace;
from os import path;
from netCDF4 import Dataset;
import calendar;

import numpy as np;

##requires 2D matrices as inputs (old and new). returns a 2D matrix which contains the difference between old and new (old-new).
#def calc_diff(old, new):
#    if old.shape != old.shape:
#        raise ValueError("old and new shape do not match.");
#    
#    #convert to standard np.arrays, but don't modify original type.
#    old2 = np.array(old);
#    new2 = np.array(new);
#    
#    output = np.zeros(old.shape);
#    
#    #
#    for y in range(0, old.shape[1]):
#        for x in range(0, old.shape[0]):
#            #if old[x,y] != -999.0 and new[x,y] != -999.0:
#            curDiff =  old2[x,y] - new2[x,y];
#            output[x,y] = curDiff;
#
#    return output;

#returns a copy of data with all nans replaced with missingValues
def nans_to_missing_value(data, missingValue=-999.0):
    fdata = np.ravel(data.copy());
    if np.ma.is_masked(fdata):
        fdata.unshare_mask();
    
    for i in range(len(fdata)):
        if np.isnan(fdata[i]):
            fdata[i] = missingValue;
    fdata.shape = data.shape;
    return fdata;

#Do everything needed to run a validation run.
def validate_run(name, year, configFilePath, referencePath, referenceFluxBudgetsFilename="_global.txt",
                 runFluxEngine=True, runFluxBudgets=True, #these two arguments force the fluxengine / flux-budgets to run overwriting any previous data.
                 validateNetCDFOutput=True, failThresholdNetCDFOutput=0.000001, #Compare netCDF output variables? Threshold is for sum of abs difference
                  validateFluxBudgets=True, failThresholdFluxBudgetsOutput=1.0, #Compare flux budget outputs? Threshold in percent
                 takahashiRun=False, verbose=False):

    validationSuccessful = True; #Assume validation was successful until something goes wrong.
    rootPath = path.abspath(path.join(__file__, "../../..")); #TODO: Should be replaced with a proper resource manager
    print "\n\nRunning validation run (year "+str(year)+") for:", name;  
    
    
    #Does the flux engine need running? if so run it
    if runFluxEngine:
        print "Running FluxEngine... (%s)" % name;
        runStatus, fe = run_fluxengine(configFilePath, [year], range(0,12), processLayersOff=True, takahashiDriver=takahashiRun, verbose=False);
    
    #extract output directory from config file
    validationOutputPath = path.join(read_config_file(configFilePath, verbose=False)["output_dir"], ''); #Add '/' to the end (required for flux_budgets)
    print validationOutputPath;
    
    
    #Do the flux-budgets need to be calculated?
    if runFluxBudgets:
        fluxBudgetsArgs = Namespace(LooseIce=False, cidataset='OIC1', cwdataset='OSFC', dir=validationOutputPath,
                                fluxdataset='OF', gridarea=0, gridareadataset='area', gridareafile='no_file',
                                icePercent=False, icedataset='P1', kwdataset='OK3', landdataset='land_proportion',
                                landfile=path.join(rootPath, 'data/onedeg_land.nc'),
                                maskdatasets=[], maskfile=path.join(rootPath, 'data/World_Seas-final-complete_IGA.nc'),
                                outroot=validationOutputPath, places=10, ref=None,
                                regions=[], verbosity=1, window=None);
        run_flux_budgets(fluxBudgetsArgs);
    
    #Compare flux budgets between validation and reference?
    if validateFluxBudgets:
        print "Comparing flux-budgets between validation run and reference data set (%s)." % name;
        newBudgetsOutputFilePath = path.join(validationOutputPath, "_global.txt");
        referenceBudgetsOutputFile = path.join(referencePath, referenceFluxBudgetsFilename);
        failedVariables = validate_flux_budgets_output(newBudgetsOutputFilePath, referenceBudgetsOutputFile, failThreshold=failThresholdFluxBudgetsOutput, verbose=verbose);
        if len(failedVariables) != 0:
            validationSuccessful = False;
    
    
    #Compare netCDF variables
    if validateNetCDFOutput:
        print "Validating netCDF output for '%s' against reference:" % name;
        failedMonths = validate_netCDF_output(validationOutputPath, referencePath, year, failThreshold=failThresholdNetCDFOutput, verbose=verbose);
        
        if len(failedMonths) != 0:
            validationSuccessful = False;
        
    return validationSuccessful;


#Compares the flux-budgets output between a new validation run and a known reference run.
#Returns a list of variables which didn't meet +/- the threshold.
def validate_flux_budgets_output(newBudgetsOutputFilePath, referenceBudgetsOutputFile, failThreshold=1.0, verbose=False):
    diffsFEv1 = calc_net_budget_percentages(newBudgetsOutputFilePath, referenceBudgetsOutputFile, verbose=False);
    
    #Check each flux budgets entry.
    failedVariables = [];
    for key in diffsFEv1:
        if diffsFEv1[key] > (100.0+failThreshold) or diffsFEv1[key] < (100.0-failThreshold):
            print key, "percentage difference from reference:", diffsFEv1[key];
            failedVariables.append(key);
    
    if len(failedVariables) == 0:
        print "All flux-budget values are within threshold limits:";
        for key in diffsFEv1:
            print "\t"+key+": "+str(diffsFEv1[key])+"%";
    else:
        print "Validation failed because %d flux-budget values were outside threshold range compared to the reference provided." % len(failedVariables);
    
    return failedVariables;


#Compares each variable in the netCDF output files for one year between a new run and a reference run.  
#Returns a list of months which fail the validation.        
def validate_netCDF_output(newOutputPath, referenceOutputPath, year, failThreshold=0.000001, verbose=False):
    #Test3: Compare individual data layers with known FE reference output
    #TODO: Use compare_datalayers.py...
    failedMonths = []; #Which months failed this validation test?
    diffsByMonth = {}; #Stores the sum of abs difference for each variable for each month
    for i in range(12):
        monthStr = calendar.month_abbr[i+1].lower();
        monthNum = "%02d"%(i+1);
        valiData = Dataset(path.join(newOutputPath, str(year), monthNum, "OceanFluxGHG-month"+monthNum+"-"+monthStr+"-"+str(year)+"-v0.nc"), 'r');
        refData = Dataset(path.join(referenceOutputPath, str(year), monthNum, "OceanFluxGHG-month"+monthNum+"-"+monthStr+"-"+str(year)+"-v0.nc"), 'r');
        
        diffs = {}; #Stores the sum of abs difference for each variable
        failedDiffs = []; #Which data layers failed the sum(abs(diffs)) threshold?
        failedMaskedDiffs = []; #Which data layers differ in the number of masked/missing values?
        diffsMissingValiData = []; #Which data layers exist in the reference dataset but don't exist in the validation dataset?
        diffsMissingRefData = []; #Which data layers exist in the validation dataset but don't exist in the reference dataset?
    
        #Compare the new and old value of each variable
        varsToCompare = refData.variables.keys();
        for varToCompare in varsToCompare:
            if varToCompare in [u"time", u"latitude", u"longitude"]:
                continue;
            
            oldVar = refData.variables[varToCompare][:];
            oldVar = np.nan_to_num(oldVar);
    
            try:
                newVar = valiData.variables[varToCompare][:];
                newVar = np.nan_to_num(newVar);
                
                if len(oldVar.shape) == 3:
                    oldVar = oldVar[0,:,:];
                    newVar = newVar[0,:,:];
                else:
                    oldVar = oldVar[:,:];
                    newVar = newVar[:,:];
                if np.any(np.isnan(oldVar)):
                    oldVar = nans_to_missing_value(oldVar);
                if np.any(np.isnan(newVar)):
                    newVar = nans_to_missing_value(newVar);
                
                #Calculate differences.
                diff = oldVar-newVar;
                sumAbsDiff = np.sum(np.abs(diff));
                diffs[varToCompare] = np.sum(np.abs(diff));
                if np.ma.is_masked(oldVar) and np.ma.is_masked(newVar):
                    oldMasked = np.sum(oldVar.mask);
                    newMasked = np.sum(newVar.mask);
                    maskedDiff = oldMasked-newMasked;
                else:
                    maskedDiff = "na";

                if verbose >= 2:
                    print varToCompare+":", sumAbsDiff, "    ", maskedDiff;
                if sumAbsDiff > failThreshold:
                    failedDiffs.append(varToCompare);
                    if verbose == 1:
                        print varToCompare+":", sumAbsDiff, "    ", maskedDiff;
                
            except KeyError:
                diffsMissingValiData.append(varToCompare);
                #print "no new value for: "+varToCompare;
        
        for v in valiData.variables.keys():
            if v not in varsToCompare: #if variable in the new output but not in the old reference output
                #print "no OLD value for: "+v;
                diffsMissingRefData.append(varToCompare);
    
        #Print number of missing or extra variables.
        if verbose:
            if len(diffsMissingValiData) != 0 or len (diffsMissingRefData) != 0:
                print "\n"+monthStr+":";
                print "%d variables exist in reference dataset but are not present in the validation dataset:" % len(diffsMissingValiData);
                print "\t", diffsMissingValiData;
                print "%d variables exist in reference dataset but are not present in the validation dataset:" % len(diffsMissingRefData);
                print "\t", diffsMissingRefData;
    
        #Did this month pass or fail?
        if len(failedDiffs) != 0 or len(failedMaskedDiffs) != 0:
            failedMonths.append(monthStr);
            print "\tMonth '%s' failed. %d variables differed in their sum of absolute differences and %d differed in the number of masked values." % (monthStr, len(failedDiffs), len(failedMaskedDiffs));
        else:
            print "\tMonth '%s' passed." % monthStr;
        
        #Store the differences for all variables (useful for interactive debugging)
        diffsByMonth[monthStr] = diffs;
    
    
    #Check for any months which failed:
    if len(failedMonths) != 0:
        print "\n\n%d months contained datalayers which deviated from the reference values:" % len(failedMonths);
        print "\t", failedMonths;
    else:
        print "Verification netCDF variables match reference run for all months.";
    
    return failedMonths; #Return the list of failed months


#def compare_netCDF_files(newFile, referenceFile, failThreshold=0.001, verbose=False):
#    pass;



