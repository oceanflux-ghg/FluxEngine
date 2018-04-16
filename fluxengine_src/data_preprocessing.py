#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 26 14:42:46 2018

@author: tomholding
"""


def no_preprocessing(datalayer):
    pass;
    #return datalayer.fdata;

#transpose data
def transpose(datalayer):
    from numpy import transpose as nptranspose;
    datalayer.data = nptranspose(datalayer.data);
    datalayer.calculate_fdata(); #must updata fdata after changing data.
    datalayer.ny = datalayer.data.shape[0];
    datalayer.nx = datalayer.data.shape[1];

#converts kelvin to celsius
def kelvin_to_celsius(datalayer):
    print "Converting %s from kelvin to celsius." % datalayer.name;
    
    for i in range(len(datalayer.fdata)):
        if datalayer.fdata[i] != datalayer.missing_value:
            datalayer.fdata[i] -= 273.15;

#converts celsius to kelvin
def celsius_to_kelvin(datalayer):
    print "Converting %s from celsius to kelvin." % datalayer.name;
    
    for i in range(len(datalayer.fdata)):
        if datalayer.fdata[i] != datalayer.missing_value:
            datalayer.fdata[i] += 273.15;

def pascal_to_millibar(datalayer):
    print "Converting %s from Pa to mbar" % datalayer.name;
    for i in range(0, len(datalayer.fdata)):
        if (datalayer.fdata[i] != datalayer.missing_value):
            datalayer.fdata[i] = datalayer.fdata[i] * 0.01;

#convert percentages to proportions
def percent_to_proportion(datalayer):
    print "Converting %s from percent to fraction" % datalayer.name;
    
    for i in range(len(datalayer.fdata)):
        if datalayer.fdata[i] != datalayer.missing_value:
            datalayer.fdata[i] /= 100.0;


##Reorder axes. The Required order by FluxEngine is lat, lon./Users/tomholding/Documents/Files/standalone_fluxengine/output/SOCATv4_WoolfRuns/no_gradients-N00/netFlux_output.log
##Reorder (lon, lat) to (lat, long)
##TODO: Same as transpose, so remove!
#def reorder_axes_lon_lat_t(datalayer):
#    print "Reordering axes from (lon, lat, time) to (time, lat, lon)";
#    from numpy import swapaxes;
#    datalayer.data = swapaxes(datalayer.data, 0, 1);
#    datalayer.calculate_fdata(); #must update fdata after changing data.


def flip_longitude(datalayer):
    print "Preprocessing %s: Flipping longitude orientation." % datalayer.name;
    from numpy import flipud;
    datalayer.data = flipud(datalayer.data);
    datalayer.calculate_fdata(); #must update fdata after changing data.

def flip_latitude(datalayer):
    print "Preprocessing %s: Flipping latitude orientation." % datalayer.name;
    from numpy import fliplr;
    datalayer.data = fliplr(datalayer.data);
    datalayer.calculate_fdata(); #must update fdata after changing data.

