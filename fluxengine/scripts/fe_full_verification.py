#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 23 14:12:13 2018

@author: tomholding
"""

import fluxengine as fluxengine_src;
from fluxengine.tools.lib_verification_tools import verify_run;
from os import path;

rootPath = "full_validation";
configPath = path.join(rootPath, "configs");
referencePath = path.join(rootPath, "reference_data");

regenerateData = True; #When true this will force the flux engine and flux-budgets calculations to run, regenerating/overwriting any previous output.

##To specify a particular version of the fluxengine (e.g. an in-development version instead of the installed version)
##   set 'feRootPath' in the code below.
##To use the installed version of FluxEngine comment out the below code.
if __name__ == "__main__":
#    import sys;
#    if len(sys.argv) != 2:
#        print "Error: Requires one input argument which specifies the root directory of the version of FluxEngine you want to validate. Found %d arguments." % (len(sys.argv)-1);
#        exit(1);
#    
#    feRootPath = sys.argv[1];
#    if feRootPath not in sys.path:
#        sys.path.insert(1, feRootPath); #needs to be at the start of the path list to preceed any other versions of FluxEngine.
#    

    print("\nThe full validation suite will be run for the FluxEngine installation located at:\n", path.dirname(fluxengine_src.__file__)); 
    ans = input("\nThis take several hours (longer on slower computers). Continue? (y/n): ");
    if ans.lower() == "y":
        print("Running validation suite...");
    else:
        print("Cancelled.");
        exit(0);
    
    
    validationOutcome = {}; #A place to store validation outcomes...
    
    name = "SOCATv4 no gradients Nightingale2000 FEv2";
    validationOutcome[name] = verify_run(name, 2010, configFilePath=path.abspath(path.join(configPath,"socatv4_no_gradients-N00.conf")),
                 referencePath=path.join(referencePath, "socatv4_no_gradients_N00_reference_FEv2"),
                 referenceFluxBudgetsFilename="no_gradients-N00_global.txt",
                 runFluxEngine=regenerateData, runFluxBudgets=regenerateData,
                 validateNetCDFOutput=True, failThresholdNetCDFOutput=0.000001,
                 validateFluxBudgets=True, failThresholdFluxBudgetsOutput=1.2, ###Should double check this against handover version.
                 takahashiRun=False, verbose=False);
    
    name = "SOCATv4 no gradients Nightingale2000 FEv3";
    validationOutcome[name] = verify_run(name, 2010, configFilePath=path.abspath(path.join(configPath,"socatv4_no_gradients-N00.conf")),
                 referencePath=path.join(referencePath, "socatv4_no_gradients_N00_reference_FEv3"),
                 referenceFluxBudgetsFilename="no_gradients-N00_global.txt",
                 runFluxEngine=False, runFluxBudgets=False, #Don't bother regenerating data, this was just done above
                 validateNetCDFOutput=True, failThresholdNetCDFOutput=0.000001,
                 validateFluxBudgets=True, failThresholdFluxBudgetsOutput=1.0,
                 takahashiRun=False, verbose=False);
    
    
    name = "SOCATv4 SST gradients Nightingale2000 FEv3";
    validationOutcome[name] = verify_run(name, 2010, configFilePath=path.abspath(path.join(configPath,"socatv4_sst_gradients-N00.conf")),
                 referencePath=path.join(referencePath,"socatv4_sst_gradients_N00_reference_FEv3"),
                 referenceFluxBudgetsFilename="_global.txt",
                 runFluxEngine=regenerateData, runFluxBudgets=regenerateData, #Don't bother regenerating data, this was just done above
                 validateNetCDFOutput=True, failThresholdNetCDFOutput=0.000001,
                 validateFluxBudgets=True, failThresholdFluxBudgetsOutput=1.0,
                 takahashiRun=False, verbose=False);
    
    
    name = "SOCATv4 SST salinity gradients Nightingale2000 FEv2";
    validationOutcome[name] = verify_run(name, 2010, configFilePath=path.join(configPath,"socatv4_sst_salinity_gradients-N00.conf"),
                 referencePath=path.join(referencePath,"socatv4_sst_salinity_N00_reference_FEv2"),
                 referenceFluxBudgetsFilename="SST_Salinity_gradients-N00_global.txt",
                 runFluxEngine=regenerateData, runFluxBudgets=regenerateData,
                 validateNetCDFOutput=True, failThresholdNetCDFOutput=0.000001,
                 validateFluxBudgets=True, failThresholdFluxBudgetsOutput=1.0,
                 takahashiRun=False, verbose=False);
    
    name = "SOCATv4 SST salinity gradients Nightingale2000 FEv3";
    validationOutcome[name] = verify_run(name, 2010, configFilePath=path.abspath(path.join(configPath,"socatv4_sst_salinity_gradients-N00.conf")),
                 referencePath=path.join(referencePath,"socatv4_sst_salinity_N00_reference_FEv3"),
                 referenceFluxBudgetsFilename="SST_Salinity_gradients-N00_global.txt",
                 runFluxEngine=False, runFluxBudgets=False, #Don't bother regenerating data, this was just done above
                 validateNetCDFOutput=True, failThresholdNetCDFOutput=0.000001,
                 validateFluxBudgets=True, failThresholdFluxBudgetsOutput=1.0,
                 takahashiRun=False, verbose=False);
    
    
    #Note that SMOS data only starts May 2010, so the year 2011 is used instead
    name = "SOCATv4 SST salinity gradients SMOS2011 Nightingale2000 FEv3";
    validationOutcome[name] = verify_run(name, 2010, configFilePath=path.abspath(path.join(configPath,"socatv4_sst_salinity_gradients_smos2011-N00.conf")),
                 referencePath=path.join(referencePath,"socatv4_sst_salinity_gradients_smos2011_N00_reference_FEv3"),
                 referenceFluxBudgetsFilename="_global.txt",
                 runFluxEngine=regenerateData, runFluxBudgets=regenerateData, #Don't bother regenerating data, this was just done above
                 validateNetCDFOutput=True, failThresholdNetCDFOutput=0.000001,
                 validateFluxBudgets=True, failThresholdFluxBudgetsOutput=1.0,
                 takahashiRun=False, verbose=False);
    
    
    name = "SOCATv4 SSTfnd, SST salinity gradients, K parameterisation: Nightingale2000, FEv3";
    validationOutcome[name] = verify_run(name, 2010, configFilePath=path.abspath(path.join(configPath,"socatv4_sst_salinity_gradients_sstfnd-N00.conf")),
                 referencePath=path.join(referencePath,"socatv4_sst_salinity_gradients_sstfnd_N00_reference_FEv3"),
                 referenceFluxBudgetsFilename="_global.txt",
                 runFluxEngine=regenerateData, runFluxBudgets=regenerateData, #Don't bother regenerating data, this was just done above
                 validateNetCDFOutput=True, failThresholdNetCDFOutput=0.000001,
                 validateFluxBudgets=True, failThresholdFluxBudgetsOutput=1.0,
                 takahashiRun=False, verbose=False);
    
    
    
    
    name = "SOCATv4 SST salinity gradients K_generic a0 FEv2";
    validationOutcome[name] = verify_run(name, 2010, configFilePath=path.abspath(path.join(configPath,"socatv4_sst_salinity_gradients-K0.conf")),
                 referencePath=path.join(referencePath,"socatv4_sst_salinity_gradients_K0_reference_FEv2"),
                 referenceFluxBudgetsFilename="SST_Salinity_gradients-K0_global.txt",
                 runFluxEngine=regenerateData, runFluxBudgets=regenerateData,
                 validateNetCDFOutput=True, failThresholdNetCDFOutput=0.000001,
                 validateFluxBudgets=True, failThresholdFluxBudgetsOutput=1.0,
                 takahashiRun=False, verbose=False);
    
    name = "SOCATv4 SST salinity gradients K_generic a0 FEv3";
    validationOutcome[name] = verify_run(name, 2010, configFilePath=path.abspath(path.join(configPath,"socatv4_sst_salinity_gradients-K0.conf")),
                 referencePath=path.join(referencePath,"socatv4_sst_salinity_gradients_K0_reference_FEv3"),
                 referenceFluxBudgetsFilename="SST_Salinity_gradients-K0_global.txt",
                 runFluxEngine=False, runFluxBudgets=False, #Don't bother regenerating data, this was just done above
                 validateNetCDFOutput=True, failThresholdNetCDFOutput=0.000001,
                 validateFluxBudgets=True, failThresholdFluxBudgetsOutput=1.0,
                 takahashiRun=False, verbose=False);
    
    
    
    name = "SOCATv4 SST salinity gradients K_generic a1 FEv2";
    validationOutcome[name] = verify_run(name, 2010, configFilePath=path.abspath(path.join(configPath,"socatv4_sst_salinity_gradients-K1.conf")),
                 referencePath=path.join(referencePath,"socatv4_sst_salinity_gradients_K1_reference_FEv2"),
                 referenceFluxBudgetsFilename="SST_Salinity_gradients-K1_global.txt",
                 runFluxEngine=regenerateData, runFluxBudgets=regenerateData,
                 validateNetCDFOutput=True, failThresholdNetCDFOutput=0.000001,
                 validateFluxBudgets=True, failThresholdFluxBudgetsOutput=1.0,
                 takahashiRun=False, verbose=False);
    
    name = "SOCATv4 SST salinity gradients K_generic a1 FEv3";
    validationOutcome[name] = verify_run(name, 2010, configFilePath=path.abspath(path.join(configPath,"socatv4_sst_salinity_gradients-K1.conf")),
                 referencePath=path.join(referencePath,"socatv4_sst_salinity_gradients_K1_reference_FEv3"),
                 referenceFluxBudgetsFilename="SST_Salinity_gradients-K1_global.txt",
                 runFluxEngine=False, runFluxBudgets=False, #Don't bother regenerating data, this was just done above
                 validateNetCDFOutput=True, failThresholdNetCDFOutput=0.000001,
                 validateFluxBudgets=True, failThresholdFluxBudgetsOutput=1.0,
                 takahashiRun=False, verbose=False);
   
    
    
    name = "SOCATv4 SST salinity gradients K_generic a2 FEv2";
    validationOutcome[name] = verify_run(name, 2010, configFilePath=path.abspath(path.join(configPath,"socatv4_sst_salinity_gradients-K2.conf")),
                 referencePath=path.join(referencePath,"socatv4_sst_salinity_gradients_K2_reference_FEv2"),
                 referenceFluxBudgetsFilename="SST_Salinity_gradients-K2_global.txt",
                 runFluxEngine=regenerateData, runFluxBudgets=regenerateData,
                 validateNetCDFOutput=True, failThresholdNetCDFOutput=0.000001,
                 validateFluxBudgets=True, failThresholdFluxBudgetsOutput=1.0,
                 takahashiRun=False, verbose=False);
    
    name = "SOCATv4 SST salinity gradients K_generic a2 FEv3";
    validationOutcome[name] = verify_run(name, 2010, configFilePath=path.abspath(path.join(configPath,"socatv4_sst_salinity_gradients-K2.conf")),
                 referencePath=path.join(referencePath,"socatv4_sst_salinity_gradients_K2_reference_FEv3"),
                 referenceFluxBudgetsFilename="SST_Salinity_gradients-K2_global.txt",
                 runFluxEngine=False, runFluxBudgets=False, #Don't bother regenerating data, this was just done above
                 validateNetCDFOutput=True, failThresholdNetCDFOutput=0.00000,
                 validateFluxBudgets=True, failThresholdFluxBudgetsOutput=1.0,
                 takahashiRun=False, verbose=False);
    
    
    
    name = "SOCATv4 SST salinity gradients K_generic a3 FEv2";
    validationOutcome[name] = verify_run(name, 2010, configFilePath=path.abspath(path.join(configPath,"socatv4_sst_salinity_gradients-K3.conf")),
                 referencePath=path.join(referencePath,"socatv4_sst_salinity_gradients_K3_reference_FEv2"),
                 referenceFluxBudgetsFilename="SST_Salinity_gradients-K3_global.txt",
                 runFluxEngine=regenerateData, runFluxBudgets=regenerateData,
                 validateNetCDFOutput=True, failThresholdNetCDFOutput=0.000001,
                 validateFluxBudgets=True, failThresholdFluxBudgetsOutput=1.0,
                 takahashiRun=False, verbose=False);
    
    name = "SOCATv4 SST salinity gradients K_generic a3 FEv3";
    validationOutcome[name] = verify_run(name, 2010, configFilePath=path.abspath(path.join(configPath,"socatv4_sst_salinity_gradients-K3.conf")),
                 referencePath=path.join(referencePath,"socatv4_sst_salinity_gradients_K3_reference_FEv3"),
                 referenceFluxBudgetsFilename="SST_Salinity_gradients-K3_global.txt",
                 runFluxEngine=False, runFluxBudgets=False, #Don't bother regenerating data, this was just done above
                 validateNetCDFOutput=True, failThresholdNetCDFOutput=0.000001,
                 validateFluxBudgets=True, failThresholdFluxBudgetsOutput=1.0,
                 takahashiRun=False, verbose=False);
    
    
                     
    
    name = "Takahashi09 all inputs FEv1";
    validationOutcome[name] = verify_run(name, 2000, configFilePath=path.abspath(path.join(configPath,"takahashi09_all_inputs.conf")),
                 referencePath=path.join(referencePath,"takahashi09_all_inputs_reference_FEv1"),
                 referenceFluxBudgetsFilename="_global.txt",
                 runFluxEngine=regenerateData, runFluxBudgets=regenerateData,
                 validateNetCDFOutput=True, failThresholdNetCDFOutput=17.0, #SC output deviates by ~16 (sum of abs difference)
                 validateFluxBudgets=True, failThresholdFluxBudgetsOutput=0.1,
                 takahashiRun=True, verbose=False);
    
    name = "Takahashi09 all inputs FEv2";
    validationOutcome[name] = verify_run(name, 2000, configFilePath=path.abspath(path.join(configPath,"takahashi09_all_inputs.conf")),
                 referencePath=path.join(referencePath,"takahashi09_all_inputs_reference_FEv2"),
                 referenceFluxBudgetsFilename="_global.txt",
                 runFluxEngine=False, runFluxBudgets=False, #Don't bother regenerating data, this was just done above
                 validateNetCDFOutput=True, failThresholdNetCDFOutput=17.0, #SC output deviates by ~16 (sum of abs difference)
                 validateFluxBudgets=True, failThresholdFluxBudgetsOutput=0.1,
                 takahashiRun=True, verbose=False);
    
    name = "Takahashi09 all inputs FEv3";
    validationOutcome[name] = verify_run(name, 2000, configFilePath=path.abspath(path.join(configPath,"takahashi09_all_inputs.conf")),
                 referencePath=path.join(referencePath,"takahashi09_all_inputs_reference_FEv3"),
                 referenceFluxBudgetsFilename="_global.txt",
                 runFluxEngine=False, runFluxBudgets=False, #Don't bother regenerating data, this was just done above
                 validateNetCDFOutput=True, failThresholdNetCDFOutput=0.000001,
                 validateFluxBudgets=True, failThresholdFluxBudgetsOutput=0.01,
                 takahashiRun=True, verbose=False);
    
    
    
    name = "Takahashi09 pCO2 no gradients FEv3 Nightingale2000";
    validationOutcome[name] = verify_run(name, 2010, configFilePath=path.abspath(path.join(configPath,"takahashi09_pco2_no_gradients-N00.conf")),
                 referencePath=path.join(referencePath,"takahashi09_pco2_no_gradients_N00_reference_FEv3"),
                 referenceFluxBudgetsFilename="_global.txt",
                 runFluxEngine=regenerateData, runFluxBudgets=regenerateData,
                 validateNetCDFOutput=True, failThresholdNetCDFOutput=0.000001,
                 validateFluxBudgets=True, failThresholdFluxBudgetsOutput=0.1,
                 takahashiRun=False, verbose=False);
    
    
    name = "Takahashi09 pCO2 SST salinity gradients FEv3 Nightingale2000";
    validationOutcome[name] = verify_run(name, 2010, configFilePath=path.abspath(path.join(configPath,"takahashi09_pco2_sst_salinity_gradients-N00.conf")),
                 referencePath=path.join(referencePath,"takahashi09_pco2_sst_salinity_gradients_N00_reference_FEv3"),
                 referenceFluxBudgetsFilename="_global.txt",
                 runFluxEngine=regenerateData, runFluxBudgets=regenerateData,
                 validateNetCDFOutput=True, failThresholdNetCDFOutput=0.000001,
                 validateFluxBudgets=True, failThresholdFluxBudgetsOutput=0.1,
                 takahashiRun=False, verbose=False);
    
    #Come brief output.
    for key in sorted(validationOutcome.keys()):
        print(key+": ", validationOutcome[key]);




#