#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 19 07:33:12 2018

@author: tomholding
"""

from netCDF4 import Dataset;
import numpy as np;
import os;

#requires 2D data.
def calc_diff(old, new):
    output = np.zeros(old.shape);
    if old.shape != old.shape:
        raise ValueError("old and new shape do not match.");
    
    old2 = np.array(old);
    new2 = np.array(new);
    
    for y in range(0, old.shape[1]):
        for x in range(0, old.shape[0]):
            #if old[x,y] != -999.0 and new[x,y] != -999.0:
            curDiff =  old2[x,y] - new2[x,y];
            output[x,y] = curDiff;

    return output;

def nans_to_missing_value(data, missingValue=-999.0):
    fdata = np.ravel(data);
    fdata.unshare_mask();
    
    for i in range(len(fdata)):
        if np.isnan(fdata[i]):
            fdata[i] = missingValue;
    fdata.shape = data.shape;
    return fdata;



path = "output/validate_takahashi09/"
#path = "output/SOCATv4_WoolfRuns/no_gradients-N00"
month = "01"; monthStr = "jan";
#old = Dataset(os.path.join(path, "reference_output/2010/"+month+"/OceanFluxGHG-month"+month+"-"+monthStr+"-2010-v0.nc"), 'r');

oldPath = "/Users/tomholding/Documents/Files/fluxengine_v3/FluxEngine/output/validate_takahashi09/2000/REF_OceanFluxGHG-month01-jan-2000-v0.nc";
newPath = "/Users/tomholding/Documents/Files/fluxengine_v3/FluxEngine/output/validate_takahashi09/2000/01/OceanFluxGHG-month01-jan-2000-v0.nc";
old = Dataset(oldPath, 'r');
new = Dataset(newPath, 'r');
#old = Dataset(os.path.join(path, "2000/"+month+"/REF_OceanFluxGHG-month"+month+"-"+monthStr+"-2010-v0.nc"), 'r');
#new = Dataset(os.path.join(path, "2010/"+month+"/OceanFluxGHG-month"+month+"-"+monthStr+"-2010-v0.nc"), 'r');

varsToCompare = old.variables.keys();

diffs = {};
for varToCompare in varsToCompare:
    if varToCompare in [u"time", u"latitude", u"longitude"]:
        continue;
    
    oldVar = old.variables[varToCompare][:];
    #oldVar = np.array(old.variables[varToCompare][:]);
    oldVar = np.nan_to_num(oldVar);
    
    #if varToCompare == "ST1_mean":
    #    varToCompare = "ST1_Kelvin_mean";
    
    try:
        newVar = new.variables[varToCompare][:];
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
        
        #oldVar.mask = False;
        #newVar.mask = False;
        diff = oldVar-newVar;
        if np.ma.is_masked(oldVar) and np.ma.is_masked(newVar):
            oldMasked = np.sum(oldVar.mask);
            newMasked = np.sum(newVar.mask);
            maskedDiff = oldMasked-newMasked;
        else:
            maskedDiff = "na";
        print varToCompare+":", np.sum(np.abs(diff)), "   ", maskedDiff;
        
        #diff = calc_diff(oldVar, newVar);
        #print np.sum(np.abs(diff));
        diffs[varToCompare] = np.sum(np.abs(diff));
        
        
        #print np.all(np.isnan(oldVar)), np.all(np.isnan(newVar));
        #print np.all(oldVar == -999.0), np.all(newVar == -999.0);
        #print "";
    except KeyError:
        print "no new value for: "+varToCompare;

for v in new.variables.keys():
    if v not in varsToCompare: #if variable in the new output but not in the old reference output
        print "no OLD value for: "+v;

import matplotlib.pyplot as plt;


#for outputVar in ["FH06", "kt", "k", "kd", "kb", "salinity", "sstskinC", "concw", "conca", "dpco2_cor", "sstfndC", "pco2_sw_cor"]:
#    if outputVar in fe.d:
#        plt.figure();
#        plt.imshow(fe.d[outputVar]);
#        plt.colorbar(); plt.clim(-3,3);
#        plt.title(outputVar);


ok3 = old.variables["OK3"][0,:,:];
nk3 = new.variables["OK3"][0,:,:];
dk3 = ok3-nk3;

ok3.mask=False;
nk3.mask=False;
dk3 = ok3-nk3;

plt.figure(); plt.imshow(dk3); plt.colorbar(); plt.title("diff OK3");
#plt.figure(); plt.imshow(ok3); plt.colorbar(); plt.title("old OK3");
#plt.figure(); plt.imshow(nk3); plt.colorbar(); plt.title("new OK3");

of = old.variables["OF"][0,:,:];
nf = new.variables["OF"][0,:,:];
df = of-nf;
plt.figure(); plt.imshow(nf); plt.colorbar(); plt.title("diff OF");


ost1 = old.variables["ST1_mean"][0,:,:];
nst1 = new.variables["ST1_mean"][0,:,:];
dst1 = ost1-nst1;
plt.figure(); plt.imshow(ost1); plt.colorbar(); plt.title("ref sstskin");
plt.figure(); plt.imshow(nst1); plt.colorbar(); plt.title("new sstskin");
plt.figure(); plt.imshow(dst1); plt.colorbar(); plt.title("ref-new sstskin");
