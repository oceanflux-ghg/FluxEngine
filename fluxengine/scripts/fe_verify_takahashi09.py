#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 15 13:55:38 2018

Runs the flux engine for using SOCATv4 pco2 data for 2010.
Compares output to a known reference output for verification.

@author: Tom Holding
"""

from os import path;

from fluxengine.core.fe_setup_tools import run_fluxengine, get_fluxengine_root;
from fluxengine.tools.lib_ofluxghg_flux_budgets import run_flux_budgets;
from fluxengine.tools.lib_compare_net_budgets import read_global_core_budgets, calc_net_budget_percentages;
from argparse import Namespace;

#Runs the verification proceedure for socat using sst salinity gradients and Nightinggale 2000 k parameterisation.
#verify_socatv4_sst_salinity_gradients_N00
def run_takahashi09_verification(verbose=True):
    #Get the path of the FluxEngine root directory.
    feRoot = get_fluxengine_root();
    
    #Run flux engine
    if verbose:
        print("Running FluxEngine for year 2000 using takahashi09 verification");
    configFilePath = path.join(feRoot, "configs", "takahashi09_verification.conf");
    runStatus, fe = run_fluxengine(configFilePath, 2000, 2000, processLayersOff=True, takahashiDriver=True, verbose=False);
    
    
    #run net budgets
    if verbose:
        print("\n\nNow calculating flux budgets...");
    outputPath = path.dirname(path.dirname(fe.runParams.output_dir));
    
    fluxBudgetsArgs = Namespace(LooseIce=False, cidataset='OIC1', cwdataset='OSFC', dir=outputPath,
                                fluxdataset='OF', gridarea=0,
                                gridareadataset='area', gridareafile='no_file', icePercent=False, icedataset='P1',
                                kwdataset='OK3', landdataset='land_proportion', landfile=path.join(feRoot, "data", "onedeg_land.nc"),
                                maskdatasets=[], maskfile=path.join(feRoot, "data", "World_Seas-final-complete_IGA.nc"),
                                outroot=path.join(outputPath, ''), places=10, ref=None,
                                regions=[], verbosity=0, window=None);
    
    run_flux_budgets(fluxBudgetsArgs);
    
    
    #Keep track of tests passed
    nTestsPassed = 0;
    N_TESTS = 4;
    
    
    #Test 1: compare similarity flux budgets output to results in Takahashi 2009 paper
    if verbose:
        print("\n\nVerification Test 1) Comparing verification output to Takahashi 2009 new flux results.");
    budgetsOutputFilePath = path.join(outputPath, "_global.txt");
    (colNames, vals) = read_global_core_budgets(budgetsOutputFilePath);
    
    takahashi09GlobalFlux = -1.42*1000.0; #TgC y^-1
    percentDiff = takahashi09GlobalFlux / vals[0] * 100.0;
    if verbose:
        print("Calculated (verification) global net CO2 flux:", vals[0]);
    if verbose:
        print("Takahashi global net CO2 flux:", takahashi09GlobalFlux);
    if verbose:
        print("Output difference from Takahashi09 global net flux is: "+str(percentDiff)+"%");
    
    if percentDiff > 105.5 or percentDiff < 94.5:
        if verbose:
            print("Verification run did not meet similarity threshold of 5.5%.");
    else:
        if verbose:
            print("Verification run is sufficiently similar to Takahashi 2009 (<5.5% difference).");
        nTestsPassed += 1;
    
    
    #Test 2: Compare similarity of flux budgets to reference data for FEv1:
    if verbose:
        print("\n\nVerification Test 2) Comparing verification output to FluxEngine v1.0 output:");
    refFEv1Path = path.join(feRoot, "data","verification_data","verification_reference_netflux","takahashi09_FEv1", "_global.txt");
    diffsFEv1 = calc_net_budget_percentages(budgetsOutputFilePath, refFEv1Path, verbose=False);
    numFailed = 0;
    
    for key in diffsFEv1:
        if diffsFEv1[key] > 100.01 or diffsFEv1[key] < 99.99:
            if verbose:
                print(key, "percentage difference from reference:", diffsFEv1[key]);
            numFailed += 1;
    
    if numFailed == 0:
        if verbose:
            print("Verification run is sufficiently similar to FluxEngine v1.0 output! All values are within threshold limits:");
        for key in diffsFEv1:
            if verbose:
                print("\t"+key+": "+str(diffsFEv1[key])+"%");
        nTestsPassed += 1;
    else:
        if verbose:
            print("Verification failed because %d values were outside threshold compared to FluxEngine v1.0." % numFailed);
    
    
    #Test 3: Compare similarity of flux budgets to reference data for FEv2:
    if verbose:
        print("\n\nVerification Test 3) Comparing verification output to FluxEngine v2.0 output:");
    refFEv2Path = path.join(feRoot, "data","verification_data","verification_reference_netflux","takahashi09_FEv2", "_global.txt");
    diffsFEv2 = calc_net_budget_percentages(budgetsOutputFilePath, refFEv2Path, verbose=False);
    numFailed = 0;
    
    for key in diffsFEv1:
        if diffsFEv2[key] > 100.01 or diffsFEv2[key] < 99.99:
            if verbose:
                print(key, "percentage difference from reference:", diffsFEv2[key]);
            numFailed += 1;
    
    if numFailed == 0:
        if verbose:
            print("Verification run is sufficiently similar to FluxEngine v2.0 output! All values are within threshold limits:");
        for key in diffsFEv2:
            if verbose:
                print("\t"+key+": "+str(diffsFEv2[key])+"%");
        nTestsPassed += 1;
    else:
        print("Verification failed because %d values were outside threshold compared to FluxEngine v2.0." % numFailed);
    
    
    #Test 4: Compare similarity of flux budgets to reference data for FEv3:
    if verbose:
        print("\n\nVerification Test 4) Comparing verification output to FluxEngine v3.0 output:");
    refFEv3Path = path.join(feRoot, "data","verification_data","verification_reference_netflux","takahashi09_FEv3", "_global.txt");
    diffsFEv3 = calc_net_budget_percentages(budgetsOutputFilePath, refFEv3Path, verbose=False);
    numFailed = 0;
    
    for key in diffsFEv1:
        if diffsFEv3[key] > 100.01 or diffsFEv2[key] < 99.99:
            if verbose:
                print(key, "percentage difference from reference:", diffsFEv3[key]);
            numFailed += 1;
    
    if numFailed == 0:
        if verbose:
            print("Verification run is sufficiently similar to FluxEngine v3.0 output! All values are within threshold limits:");
        for key in diffsFEv3:
            if verbose:
                print("\t"+key+": "+str(diffsFEv3[key])+"%");
        nTestsPassed += 1;
    else:
        print("Verification failed because %d values were outside threshold compared to FluxEngine v3.0." % numFailed);
    
    #Summary:
    if nTestsPassed == N_TESTS:
        if verbose:
            print("\n\n***** Verification successful! *****");
        verificationSuccessful = True;
    else:
        if verbose:
            print("\n\n***** Verification failed! *****");
        verificationSuccessful = False;
    
    print("%d of %d verification tests passed." % (nTestsPassed, N_TESTS));

    return {"run status":runStatus, "tests passed":nTestsPassed, "total tests":N_TESTS,
            "total net flux percentage difference from Takahashi 2009 result":percentDiff,
            "net budgets comparison with FluxEngine v1.0":diffsFEv1,
            "net budgets comparison with FluxEngine v2.0":diffsFEv2,
            "verification successful?":verificationSuccessful};

if __name__ == "__main__":
    run_takahashi09_verification(verbose=True);



