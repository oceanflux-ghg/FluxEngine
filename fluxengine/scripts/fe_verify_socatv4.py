#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 15 13:55:38 2018

Runs the flux engine for using SOCATv4 pco2 data for 2010.
Compares output to a known reference output for verification.

@author: verwirrt
"""

#from ofluxghg_run import run_fluxengine;
from fluxengine.core.fe_setup_tools import run_fluxengine, get_fluxengine_root;
from fluxengine.tools.lib_ofluxghg_flux_budgets import run_flux_budgets;
from fluxengine.tools.lib_compare_net_budgets import calc_net_budget_percentages
from argparse import Namespace;
from os import path;

#Runs the verification proceedure for socat using sst salinity gradients and Nightinggale 2000 k parameterisation.
def run_socat_sst_salinity_gradients_N00_verification(verbose=True):
    #Get the path of the FluxEngine root directory.
    feRoot = get_fluxengine_root();
    
    #Run flux engine
    if verbose:
        print("Running FluxEngine for year 2010...");
    configFilePath = path.join(feRoot, "configs/socatv4_sst_salinity_gradients-N00.conf");
    runStatus, fe = run_fluxengine(configFilePath, 2010, 2010, processLayersOff=True, verbose=False);
    
    
    #run net budgets
    if verbose:
        print("\n\nNow calculating flux budgets...");
    outputFilePath = path.dirname(path.dirname(fe.runParams.output_dir));
    fluxBudgetsArgs = Namespace(LooseIce=False, cidataset='OIC1', cwdataset='OSFC',
                                dir=path.join(outputFilePath, ''), fluxdataset='OF', gridarea=0,
                                gridareadataset='area', gridareafile='no_file', icePercent=False, icedataset='P1',
                                kwdataset='OK3', landdataset='land_proportion', landfile=path.join(feRoot, 'data/onedeg_land.nc'),
                                maskdatasets=[], maskfile=path.join(feRoot, 'data/World_Seas-IHO-mask.nc'),
                                outroot=path.join(outputFilePath, ''), places=10, ref=None,
                                regions=[], verbosity=0, window=None);
    run_flux_budgets(fluxBudgetsArgs);
    
    
    #compare similarity flux budgets output between new and ref runs
    if verbose:
        print("\n\nComparing output to reference data...");
    newPath = path.join(outputFilePath, "_global.txt");
    refPath = path.join(feRoot, "data", "verification_data", "verification_reference_netflux", "socatv4_sst_salinity_N00_reference_FEv3", "SST_Salinity_gradients-N00_global.txt");
    diffs = calc_net_budget_percentages(newPath, refPath, verbose=False);
    
    numFailed = 0;
    for key in diffs:
        if diffs[key] > 100.0001 or diffs[key] < 99.9999:
            if verbose:
                print(key, "percentage difference from reference:", diffs[key]);
            numFailed += 1;
    
    if numFailed == 0:
        if verbose:
            print("Verification successful! All values are within threshold limits:");
            verificationSuccessful = True;
        for key in diffs:
            if verbose:
                print("\t"+key+": "+str(diffs[key])+"%");
    else:
        if verbose:
            print("\Verification failed because %d values were outside threshold limits." % numFailed);
            verificationSuccessful = True;
    
    return {"run status":runStatus, "number of net budgets exceeding threshold":numFailed,
            "verification successful?":verificationSuccessful, "percentage difference from reference (dictionary)":diffs};
    
if __name__ == "__main__":
    run_socat_sst_salinity_gradients_N00_verification(verbose=True);

