#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 27 08:29:33 2018

@author: tomholding
"""

def calc_mean(data, missing_value = -999.0):
    
    if hasattr(data, "fdata"): #if we have a datalayer or datalayer like object
        data = data.fdata;
    
    t=0.0; n=0;
    for i in range(len(data)):
        if data[i] != missing_value:
            t+= data[i];
            n+=1;
    try:
        return t/n;
    except ZeroDivisionError:
        return "Divide by zero - All missing_value.";