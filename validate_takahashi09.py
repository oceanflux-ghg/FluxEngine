#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 15 13:55:38 2018

Runs the flux engine for using SOCATv4 pco2 data for 2010.
Compares output to a known reference output for verification.

@author: verwirrt
"""

from ofluxghg_run import run_fluxengine;
from fluxengine_src.tools.ofluxghg_flux_budgets import run_flux_budgets;
from fluxengine_src.tools.compare_net_budgets import read_global_core_budgets, calc_net_budget_percentages
from argparse import Namespace;
from os import getcwd, path;

#Run flux engine
print "Running FluxEngine for year 2000 using takahashi09 validation";
configFilePath = "configs/takahashi09_validation.conf";
runStatus = run_fluxengine(configFilePath, [2000], range(0,12), processLayersOff=True, takahashiDriver=True, verbose=False);


#run net budgets
print "\n\nNow calculating flux budgets...";
outputPath = path.join(getcwd(), "output/validate_takahashi09/")

fluxBudgetsArgs = Namespace(LooseIce=False, cidataset='OIC1', cwdataset='OSFC', dir=path.join(getcwd(), outputPath),
                            fluxdataset='OF', gridarea=0,
                            gridareadataset='area', gridareafile='no_file', icePercent=False, icedataset='P1',
                            kwdataset='OK3', landdataset='land_proportion', landfile=path.join(getcwd(), 'data/onedeg_land.nc'),
                            maskdatasets=[], maskfile=path.join(getcwd(), 'data/World_Seas-final-complete_IGA.nc'),
                            outroot=path.join(getcwd(), outputPath), places=10, ref=None,
                            regions=[], verbosity=0, window=None);

run_flux_budgets(fluxBudgetsArgs);


#Keep track of tests passed
nTestsPassed = 0;
N_TESTS = 3;


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
refFEv1Path = path.join(getcwd(), "data/validation_data/validation_reference_output/takahashi09_FEv1/", "_global.txt");
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
refFEv2Path = path.join(getcwd(), "data/validation_data/validation_reference_output/takahashi09_FEv2/", "_global.txt");
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


#Summary:
if nTestsPassed == N_TESTS:
    print "\n\n***** Validation successful! *****";
else:
    print "\n\n***** Validation failed! *****";

print "%d of %d validation tests passed." % (nTestsPassed, N_TESTS);






