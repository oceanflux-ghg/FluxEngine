#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 18 12:48:03 2018
Runs a more indepth version of the validation which:
    * Generates output using all inputs from Takahashi2009
    * Compared this validation run with the results from the Takahashi2009 paper
    * Compares the flux-budgets output of this validation run with known reference
        output from versions 1.0, 2.0 and 3.0 of FluxEngine
    * Compares each datalayer in the validation run with a known reference for
        version 3.0 of FluxEngine
@author: tomholding
"""

#Use SET_PATHS to ensure that we are running the development version of the FluxEngine and not any previously installed version.
import SET_PATHS;
rootPath, srcPath = SET_PATHS.set_paths(verbose=True);


#from ofluxghg_run import run_fluxengine;
from fe_setup_tools import run_fluxengine;
from tools.ofluxghg_flux_budgets import run_flux_budgets;
from tools.compare_net_budgets import read_global_core_budgets, calc_net_budget_percentages;
from argparse import Namespace;
from os import path;
from netCDF4 import Dataset;
import calendar;

import numpy as np;

#requires 2D matrices as inputs (old and new). returns a 2D matrix which contains the difference between old and new (old-new).
def calc_diff(old, new):
    if old.shape != old.shape:
        raise ValueError("old and new shape do not match.");
    
    #convert to standard np.arrays, but don't modify original type.
    old2 = np.array(old);
    new2 = np.array(new);
    
    output = np.zeros(old.shape);
    
    #
    for y in range(0, old.shape[1]):
        for x in range(0, old.shape[0]):
            #if old[x,y] != -999.0 and new[x,y] != -999.0:
            curDiff =  old2[x,y] - new2[x,y];
            output[x,y] = curDiff;

    return output;

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


#Run flux engine with takahashi2009 inputs
print "Running FluxEngine for year 2000 using takahashi09 validation";
year = 2000;
configFilePath = path.join(rootPath, "configs/takahashi09_validation.conf");
runStatus, fe = run_fluxengine(configFilePath, [year], range(0,12), processLayersOff=True, takahashiDriver=True, verbose=False);

#run net budgets
print "\n\nNow calculating flux budgets...";
outputPath = path.dirname(path.dirname(fe.runParams.output_dir));


fluxBudgetsArgs = Namespace(LooseIce=False, cidataset='OIC1', cwdataset='OSFC', dir=outputPath,
                            fluxdataset='OF', gridarea=0,
                            gridareadataset='area', gridareafile='no_file', icePercent=False, icedataset='P1',
                            kwdataset='OK3', landdataset='land_proportion', landfile=path.join(rootPath, 'data/onedeg_land.nc'),
                            maskdatasets=[], maskfile=path.join(rootPath, 'data/World_Seas-final-complete_IGA.nc'),
                            outroot=outputPath, places=10, ref=None,
                            regions=[], verbosity=0, window=None);

run_flux_budgets(fluxBudgetsArgs);


#Keep track of tests passed
N_TESTS = 5;
nTestsPassed = 0;


#Test 1: compare similarity flux budgets output to results in Takahashi 2009 paper
print "\n\nValidation Test 1) Comparing validation output to Takahashi 2009 new flux results.";
budgetsOutputFilePath = path.join(outputPath, "_global.txt");
(colNames, vals) = read_global_core_budgets(budgetsOutputFilePath);

takahashi09GlobalFlux = -1.42*1000.0; #TgC y^-1
percentDiff = takahashi09GlobalFlux / vals[0] * 100.0;
print "Calculated (validation) global net CO2 flux:", vals[0];
print "Takahashi global net CO2 flux:", takahashi09GlobalFlux;
print "Output difference from Takahashi09 global net flux is: "+str(percentDiff)+"%";

if percentDiff > 105.5 or percentDiff < 94.5:
    print "Validation run did not meet similarity threshold of 5.5%.";
else:
    print "Validation run is sufficiently similar to Takahashi 2009 (<5.5% difference).";
    nTestsPassed += 1;


#Test 2: Compare similarity of flux budgets to reference data for FEv1:
print "\n\nValidation Test 2) Comparing validation output to FluxEngine v1.0 output:";
refFEv1Path = path.join(rootPath, "data/validation_data/validation_reference_output/takahashi09_FEv1/", "_global.txt");
diffsFEv1 = calc_net_budget_percentages(budgetsOutputFilePath, refFEv1Path, verbose=False);
numFailed = 0;

for key in diffsFEv1:
    if diffsFEv1[key] > 100.01 or diffsFEv1[key] < 99.99:
        print key, "percentage difference from reference:", diffsFEv1[key];
        numFailed += 1;

if numFailed == 0:
    print "Validation run is sufficiently similar to FluxEngine v1.0 output! All values are within threshold limits:";
    for key in diffsFEv1:
        print "\t"+key+": "+str(diffsFEv1[key])+"%";
    nTestsPassed += 1;
else:
    print "Validation failed because %d values were outside threshold compared to FluxEngine v1.0." % numFailed;

#Test 3: Compare similarity of flux budgets to reference data for FEv2:
print "\n\nValidation Test 3) Comparing validation output to FluxEngine v2.0 output:";
refFEv2Path = path.join(rootPath, "data/validation_data/validation_reference_output/takahashi09_FEv2/", "_global.txt");
diffsFEv2 = calc_net_budget_percentages(budgetsOutputFilePath, refFEv2Path, verbose=False);
numFailed = 0;

for key in diffsFEv1:
    if diffsFEv2[key] > 100.01 or diffsFEv2[key] < 99.99:
        print key, "percentage difference from reference:", diffsFEv2[key];
        numFailed += 1;

if numFailed == 0:
    print "Validation run is sufficiently similar to FluxEngine v2.0 output! All values are within threshold limits:";
    for key in diffsFEv2:
        print "\t"+key+": "+str(diffsFEv2[key])+"%";
    nTestsPassed += 1;
else:
    print "Validation failed because %d values were outside threshold compared to FluxEngine v2.0." % numFailed;


#Test 4: Compare similarity of flux budgets to reference data for FEv2:
print "\n\nValidation Test 4) Comparing validation output to FluxEngine v3.0 reference output:";
refFEv2Path = path.join(rootPath, "data/validation_data/validation_reference_output/takahashi09_FEv3/", "_global.txt");
diffsFEv2 = calc_net_budget_percentages(budgetsOutputFilePath, refFEv2Path, verbose=False);
numFailed = 0;

for key in diffsFEv1:
    if diffsFEv2[key] > 100.01 or diffsFEv2[key] < 99.99:
        print key, "percentage difference from reference:", diffsFEv2[key];
        numFailed += 1;

if numFailed == 0:
    print "Validation run is sufficiently similar to FluxEngine v3.0 reference output! All values are within threshold limits:";
    for key in diffsFEv2:
        print "\t"+key+": "+str(diffsFEv2[key])+"%";
    nTestsPassed += 1;
else:
    print "Validation failed because %d values were outside threshold compared to FluxEngine v3.0 reference." % numFailed;


#Test5: Compare individual data layers with known FEv3.0 reference output
print "\n\nValidation Test 5) Comparing netCDF4 variables in validation output to FluxEngine v3.0 reference output:";
failedMonths = []; #Which months failed this validation test?
diffsByMonth = {}; #Stores the sum of abs difference for each variable for each month
for i in range(12):
    monthStr = calendar.month_abbr[i+1].lower();
    monthNum = "%02d"%(i+1);
    valiData = Dataset(path.join(outputPath, str(year), monthNum, "OceanFluxGHG-month"+monthNum+"-"+monthStr+"-"+str(year)+"-v0.nc"), 'r');
    refData = Dataset(path.join(rootPath, "data/validation_data/validation_reference_output/takahashi09_FEv3", str(year), monthNum, "OceanFluxGHG-month"+monthNum+"-"+monthStr+"-"+str(year)+"-v0.nc"), 'r');
    
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
            #newVar = np.array(new.variables[varToCompare][:]);
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
            
            #Calculate differences in the number of masked values between validation and reference variables.
            if np.ma.is_masked(oldVar) and np.ma.is_masked(newVar):
                oldMasked = np.sum(oldVar.mask);
                newMasked = np.sum(newVar.mask);
                maskedDiff = oldMasked-newMasked;
            else:
                maskedDiff = "na";
            
            #Is there a difference in the number of masked values?
            if maskedDiff != "na" and maskedDiff > 0:
                failedMaskedDiffs.append(varToCompare);
            
            #Calculate and check sum of absolute difference is within threshold value
            diff = calc_diff(oldVar, newVar);
            diffs[varToCompare] = np.sum(np.abs(diff));
            sumAbsDiff = np.sum(np.abs(diff));
            if np.sum(np.abs(diff)) > 0.000001:
                failedDiffs.append(varToCompare);
            
        except KeyError:
            diffsMissingValiData.append(varToCompare);
            #print "no new value for: "+varToCompare;
    
    for v in valiData.variables.keys():
        if v not in varsToCompare: #if variable in the new output but not in the old reference output
            #print "no OLD value for: "+v;
            diffsMissingRefData.append(varToCompare);

    #Print number of missing or extra variables.
    if len(diffsMissingValiData) != 0 or len (diffsMissingRefData) != 0:
        print "\n"+monthStr+":";
        print "%d variables exist in reference dataset but are not present in the validation dataset:" % diffsMissingValiData;
        print "\t", diffsMissingValiData;
        print "%d variables exist in reference dataset but are not present in the validation dataset:" % diffsMissingRefData;
        print "\t", diffsMissingRefData;

    #Did this month pass or fail?
    if len(failedDiffs) != 0 or len(failedMaskedDiffs) != 0:
        failedMonths.append(monthStr);
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
    nTestsPassed += 1;
    

#Summary:
if nTestsPassed == N_TESTS:
    print "\n\n***** Validation successful! *****";
    validationSuccessful = True;
else:
    print "\n\n***** Validation failed! *****";
    validationSuccessful = False;

print "%d of %d validation tests passed." % (nTestsPassed, N_TESTS);

#return {"run status":runStatus, "tests passed":nTestsPassed, "total tests":N_TESTS,
#        "total net flux percentage difference from Takahashi 2009 result":percentDiff,
#        "net budgets comparison with FluxEngine v1.0":diffsFEv1,
#        "net budgets comparison with FluxEngine v2.0":diffsFEv2,
#        "validation successful?":validationSuccessful};
