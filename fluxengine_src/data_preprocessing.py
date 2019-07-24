#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 26 14:42:46 2018

@author: tomholding
"""


#Used as default pre-processing function. Do not modify.
def no_preprocessing(datalayer):
    pass;

#transpose data
def transpose(datalayer):
    from numpy import transpose as nptranspose;
    datalayer.data = nptranspose(datalayer.data);
    datalayer.calculate_fdata(); #must updata fdata after changing data.
    datalayer.ny = datalayer.data.shape[0];
    datalayer.nx = datalayer.data.shape[1];

#converts kelvin to celsius
def kelvin_to_celsius(datalayer):
    print("Converting %s from kelvin to celsius." % datalayer.name);
    
    for i in range(len(datalayer.fdata)):
        if datalayer.fdata[i] != datalayer.missing_value:
            datalayer.fdata[i] -= 273.15;

#converts celsius to kelvin
def celsius_to_kelvin(datalayer):
    print("Converting %s from celsius to kelvin." % datalayer.name);
    
    for i in range(len(datalayer.fdata)):
        if datalayer.fdata[i] != datalayer.missing_value:
            datalayer.fdata[i] += 273.15;

def pascal_to_millibar(datalayer):
    print("Converting %s from Pa to mbar" % datalayer.name);
    for i in range(0, len(datalayer.fdata)):
        if (datalayer.fdata[i] != datalayer.missing_value):
            datalayer.fdata[i] = datalayer.fdata[i] * 0.01;

#convert percentages to proportions
def percent_to_proportion(datalayer):
    print("Converting %s from percent to fraction" % datalayer.name);
    
    for i in range(len(datalayer.fdata)):
        if datalayer.fdata[i] != datalayer.missing_value:
            datalayer.fdata[i] /= 100.0;

def nano_to_micro(datalayer):
    print("Converting %s from nano<units> to micro<units>" % datalayer.name);

    for i in range(len(datalayer.fdata)):
        if datalayer.fdata[i] != datalayer.missing_value:
            datalayer.fdata[i] /= 1000.0;

##Reorder axes. The Required order by FluxEngine is lat, lon.
##Reorder (lon, lat) to (lat, long)
##TODO: Same as transpose, so remove!
#def reorder_axes_lon_lat_t(datalayer):
#    print "Reordering axes from (lon, lat, time) to (time, lat, lon)";
#    from numpy import swapaxes;
#    datalayer.data = swapaxes(datalayer.data, 0, 1);
#    datalayer.calculate_fdata(); #must update fdata after changing data.


def flip_longitude(datalayer):
    print("Preprocessing %s: Flipping longitude orientation." % datalayer.name);
    from numpy import flipud;
    datalayer.data = flipud(datalayer.data);
    datalayer.calculate_fdata(); #must update fdata after changing data.

def flip_latitude(datalayer):
    print("Preprocessing %s: Flipping latitude orientation." % datalayer.name);
    from numpy import fliplr;
    datalayer.data = fliplr(datalayer.data);
    datalayer.calculate_fdata(); #must update fdata after changing data.


#E.g. as rough way to approximate second moment (not recommented)
def pow2(datalayer):
    datalayer.fdata = datalayer.fdata**2;

#E.g. as rough way to approximate third moment (not recommended)
def pow3(datalayer):
    datalayer.fdata = datalayer.fdata**3;

#Resamples the datalayer which used latitude grid lines starting on the 'line'
#to be at the centre of grid points. This results in a latitude dimension with
#a size one smaller than the original. E.g. for a 1x1 global grid, on the line 
#latitudes would have 181 points (-90 and 90, inclusive). This function will convert
#to 180 points (-89.5 to 89.5, inclusive).
#It is and is achieved by using a rolling mean of two values across the latitude timension.
def lat_grid_lines_to_centre_of_cells(datalayer):
    from numpy import empty, mean;
    newData = empty((datalayer.data.shape[0]-1, datalayer.data.shape[1]), dtype=float);
    
    for i in range(datalayer.data.shape[0]-1):
        newData[i,:] = mean(datalayer.data[[i,i+1],:], axis=0);
    datalayer.data = newData;
    datalayer.calculate_fdata(); #Update the 'fdata' after changing 'data'
    datalayer.ny, datalayer.nx = datalayer.data.shape; #Shape has changed so this must be updated too

#Rolls the dataset 180 degree in the east-west direction. This is useful for converting between
#longitude conventions of -180 to 180 and 0 to 360.
def longitude_roll_180(datalayer):
    from numpy import roll;
    datalayer.data = roll(datalayer.data, 180, axis=1);
    datalayer.calculate_fdata(); #recalculate the 'fdata' after changing 'data'


#Converts 'wave to ocean energy' (foc in WaveWatch) to dissipation rate of turbulent kinetic energy (epsilon)
#Calculates dissipation rate of turbulent energy in the top 10m (mean)
def foc_to_epsilon(datalayer):
    print("Converting datalayer '%s' from 'wave to ocean energy' to 'dissipation rate of turbulent kinetic energy'." % datalayer.name);
    
    waterDensity = 1026.0;
    for i in range(len(datalayer.fdata)):
        if datalayer.fdata[i] != datalayer.missing_value:
            datalayer.fdata[i] = datalayer.fdata[i] / (10.0 * waterDensity); #Assumes wave energy (input at surface) is dissipated through over a 10m layer
            

#Converts 'wave to ocean energy' (foc in WaveWatch) to dissipation rate of turbulent kinetic energy (epsilon)
#Calculates dissipation rate of turbulent energy at 0.2m by fitting mean dissipation rate over the top 10m using the relationship in Craig and Banner 1994 (figure 8)
def foc_to_epsilon_craig1994(datalayer):
    print("Converting datalayer '%s' from 'wave to ocean energy' to 'dissipation rate of turbulent kinetic energy'." % datalayer.name);
    import numpy as np;
    from scipy import optimize;
    
    #Calculate mean dissipation from FOC
    focMeanEps = datalayer.fdata.copy();
    waterDensity = 1026.0;
    for i in range(len(datalayer.fdata)):
        if datalayer.fdata[i] != datalayer.missing_value:
            focMeanEps[i] = datalayer.fdata[i] / (10.0 * waterDensity); #Assumes wave energy (input at surface) is dissipated through over a 10m layer
    
    #Fit beta for craig depth curve for each grid cell
    def calc_eps_mean(beta, depth0, depth1): #calculate epsilon mean given beta
        b, c = (-2.00270875e+01, -7.75494213e-01);
        integral = (beta/b) * np.exp(b*depth1 + c) - (beta/b) * np.exp(b*depth0 + c);    
        return integral / (depth1-depth0);
    
    targetEpsMean = None;
    def calib_error(calibrationFactor): #Error function for the fitting algorithm
        testEpsMean = calc_eps_mean(calibrationFactor, 0.0, 10.0);
        return np.abs(testEpsMean-targetEpsMean);
    
    #given an epsilon-depth relationship, calculate epsilon at a given depth
    #depth in metres
    def depth_function(depth, a=1.85748581e-03, b=-2.00270875e+01, c=-7.75494213e-01):
        dissipationRate = beta*np.exp(b*depth + c);
        return dissipationRate;


    for i in range(len(datalayer.fdata)):
        if datalayer.fdata[i] != datalayer.missing_value:
            targetEpsMean = focMeanEps[i];
            
            #find beta that matches means.
            result = optimize.minimize_scalar(calib_error);
            
            print(i, targetEpsMean, result.success, result.x); # check if solver was successful
            beta = result.x;
            #Calculate epsilon at 2cm and write to data layer.
            datalayer.fdata[i] = depth_function(0.2, a=beta);


