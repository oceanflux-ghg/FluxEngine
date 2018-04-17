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
from fluxengine_src.tools.compare_net_budgets import read_global_core_budgets
from argparse import Namespace;
from os import getcwd, path;

#Run flux engine
print "Running FluxEngine for year 2000 using takahashi09 validation";
configFilePath = "configs/takahashi09_validation.conf";
#runStatus = run_fluxengine(configFilePath, [2000], range(0,12), processLayersOff=True, takahashiDriver=True, verbose=False);


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


#compare similarity flux budgets output between new and ref runs
print "\n\nComparing output to reference data...";
budgetsOutputFilePath = path.join(outputPath, "_global.txt");
#(colNames, vals) = read_global_core_budgets(budgetsOutputFilePath);
(colNames, vals) = read_global_core_budgets("/Users/tomholding/Documents/Files/fluxengine_v3/fluxengine_v3_ref_runs/data/validation_data/validation_reference_output/takahashiv2009_allinputs_no_gradients/_global.txt")

takahashi09GlobalFlux = -1.42*1000.0; #TgC y^-1
percentDiff = takahashi09GlobalFlux / vals[0] * 100.0;
print "Calculated global net CO2 flux:", vals[0];
print "Takahashi global net CO2 flux:", takahashi09GlobalFlux;
print "Output difference from Takahashi09 global net flux is: "+str(percentDiff)+"%";

if percentDiff > 105.0 or percentDiff < 95.0:
    print "Validation failed.";
else:
    print "Validation successful!";

for i, field in enumerate(colNames):
    print vals[i]/1000.0, "\t("+field+")";
