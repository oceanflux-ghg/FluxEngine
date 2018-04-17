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
from fluxengine_src.tools.compare_net_budgets import calc_net_budget_percentages
from argparse import Namespace;
from os import getcwd, path;

#Run flux engine
print "Running FluxEngine for year 2010...";
configFilePath = "configs/socatv4_sst_salinity_gradients-N00.conf";
runStatus = run_fluxengine(configFilePath, [2010], range(0,12), processLayersOff=True, verbose=False);


#run net budgets
print "\n\nNow calculating flux budgets...";
fluxBudgetsArgs = Namespace(LooseIce=False, cidataset='OIC1', cwdataset='OSFC', dir=path.join(getcwd(),
                            'output/validate_socatv4_sst_salinity_N00/'), fluxdataset='OF', gridarea=0,
                            gridareadataset='area', gridareafile='no_file', icePercent=False, icedataset='P1',
                            kwdataset='OK3', landdataset='land_proportion', landfile=path.join(getcwd(), 'data/onedeg_land.nc'),
                            maskdatasets=[], maskfile=path.join(getcwd(), 'data/World_Seas-IHO-mask.nc'),
                            outroot=path.join(getcwd(), 'output/validate_socatv4_sst_salinity_N00/'), places=10, ref=None,
                            regions=[], verbosity=0, window=None);
run_flux_budgets(fluxBudgetsArgs);


#compare similarity flux budgets output between new and ref runs
print "\n\nComparing output to reference data...";
newPath = "output/validate_socatv4_sst_salinity_N00/_global.txt";
refPath = "data/validation_data/validation_reference_output/socatv4_sst_salinity_N00_reference/SST_Salinity_gradients-N00_global.txt";
diffs = calc_net_budget_percentages(newPath, refPath, verbose=False);

numFailed = 0;
for key in diffs:
    if diffs[key] > 101.0 or diffs[key] < 99.0:
        print key, "percentage difference from reference:", diffs[key];
        numFailed += 1;

if numFailed == 0:
    print "\nValidation successful! All values are within threshold limits:";
    for key in diffs:
        print "\t"+key+": "+str(diffs[key])+"%";
else:
    print "\nValidation failed because %d values were outside threshold limits." % numFailed;

