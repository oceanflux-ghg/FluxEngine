#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 16 14:54:46 2018

This file defines process indicator layer functors which allow process indicator layers to be
computed in a modular fashion.

@author: tomholding
"""

from .datalayer import DataLayer;
from numpy import arange;

#Base class from which process indicator layers are derived.
class ProcessIndicatorBase:
    def __init__(self):
        self.name = self.__class__.__name__;
    
    def input_names(self):
        raise NotImplementedError("input_names must be implemented in all derived classes.");
    
    def output_names(self):
        raise NotImplementedError("output_names must be implementedin all derived classes.");
    
    #This should take 
    def __call__(self, data):
        raise NotImplementedError("__call__ must be implemented in all derived classes.");



#Example showing required methods and format.
#Class name should correspond to a valid value of 'k_parameterisation' in the config file, e.g. "k_Wanninkhof1992".
class process_indicator_example(ProcessIndicatorBase):
    #Optional parameters can be passed as arguments, but their names must correspond to names in the config file.
    def __init__(self, someParameter=0.0):
        self.name = self.__class__.__name__;
        self.someParameter = someParameter;
    
    #Must return a list of strings corresponding to the input data layers required by the calculation. These must already exist.
    def input_names(self):
        return ["list", "of", "DataLayer", "names"];
    
    #Must return a list of strings corresponding to the names of data layers which it will write to.
    #These may or may not already exist becauce empty datalayers will be created before execution if neccessary.
    def output_names(self):
        return ["these", "are", "written", "to"];
    
    #Main calculation is performed here. Input and output datalayers can be accessed from 'data'.
    #Output layers should be modified in place (i.e. without copying), and return True or False to indicate successful execution.
    def __call__(self, data):
        print(self.long_name, "with 'someParameter' as", self.someParameter);
        return True;


#Calculate low_wind indicator layer
class low_wind_indicator(ProcessIndicatorBase):
    def __init__(self, lowWindThreshold=5.0):
        self.name = self.__class__.__name__;
        self.lowWindThreshold = lowWindThreshold;
        
    def input_names(self):
        return ["windu10"];
    
    def output_names(self):
        return ["low_wind"];
        
    def __call__(self, data):
        function = "(process_indicator_layers.py: "+self.name+".__call__)";
        try:
            windu10 = data["windu10"].fdata;
            low_wind = data["low_wind"].fdata;
        except KeyError as e:
            print("%s: Required data layer for process indicator layer was not found." % function);
            print(type(e), e.args);
            return False;
       
        #Main logic; generate low_wind mask
        for i in arange(len(windu10)):
            if (windu10[i] != DataLayer.missing_value):
                if (windu10[i] <= self.lowWindThreshold):
                    low_wind[i] = 1;
                elif (windu10[i] > self.lowWindThreshold):
                    low_wind[i] = 0;
            else:
                low_wind[i] = DataLayer.missing_value_int;
        
        return True;


#biological activity layer (classifies level of biological activity into 1 (low) 2 (mid) or 3 (high))
class bioclass_indicator(ProcessIndicatorBase):
    def __init__(self, lowThreshold=0.5, midThreshold=1.0):
        self.name = self.__class__.__name__;
        self.lowThreshold = lowThreshold;
        self.midThreshold = midThreshold;
        
    def input_names(self):
        return ["biology"];
    def output_names(self):
        return ["bioclass"];
    
    def __call__(self, data):
        function = "(process_indicator_layers.py: "+self.name+".__call__)";
        try:
            biology = data["biology"].fdata;
            bioclass = data["bioclass"].fdata;
        except KeyError as e:
            print("%s: Required data layer for process indicator layer was not found." % function);
            print(type(e), e.args);
            return False;
        
        for i in arange(len(biology)):
            #adding in the > 0.0 condition as preliminary ESA CCI data are missing a load of data attributes
            if biology[i] != DataLayer.missing_value and biology[i] > 0.0:
                if  biology[i] <= self.lowThreshold:
                    bioclass[i] = 1;
                elif biology[i] <= self.midThreshold:
                    bioclass[i] = 2;
                elif biology[i] > self.midThreshold:
                    bioclass[i] = 3;
            else:
                bioclass[i] = DataLayer.missing_value_int;
        
        return True;


#diurnal warming when (sstskin - sstfnd) > differenceThreshold
class diurnal_warming_indicator(ProcessIndicatorBase):
    #differenceThreshold is the degree warmer sstskin must be than sstfnd to be considered a region of diurnal warming
    def __init__(self, differenceThreshold=0.5):
        self.name = self.__class__.__name__;
        self.differenceThreshold = differenceThreshold;
    
    def input_names(self):
        return ["sstskin", "sstfnd"];
    def output_names(self):
        return ["diurnal_warming"];
    
    def __call__(self, data):
        function = "(process_indicator_layers.py: "+self.name+".__call__)";
        try:
            sstskin = data["sstskin"].fdata;
            sstfnd = data["sstfnd"].fdata;
            diurnal_warming = data["diurnal_warming"].fdata;
        except KeyError as e:
            print("%s: Required data layer for process indicator layer was not found." % function);
            print(type(e), e.args);
            return False;
        
        #diurnal warming when (sstskin - sstfnd) > differenceThreshold
        for i in arange(len(sstskin)):
            if (sstskin[i] != DataLayer.missing_value) and (sstfnd[i] != DataLayer.missing_value):
                sstDiff = sstskin[i] - sstfnd[i]
                if ( (sstskin[i] > sstfnd[i]) and (sstDiff > self.differenceThreshold) ): #differenceThreshold condition allows focus on strong gradients
                    diurnal_warming[i] = 1;
                else:
                    diurnal_warming[i] = 0;
            else:
                diurnal_warming[i] = DataLayer.missing_value_int;
        
        return True;


#oceanic basins
class oceanic_basins_indicator(ProcessIndicatorBase):
    #differenceThreshold is the degree warmer sstskin must be than sstfnd to be considered a region of diurnal warming
    def __init__(self, differenceThreshold=0.5):
        self.name = self.__class__.__name__;
        self.differenceThreshold = differenceThreshold;
    
    def input_names(self):
        return ["atlantic_ocean_mask", "pacific_ocean_mask", "southern_ocean_mask", "indian_ocean_mask"];
    def output_names(self):
        return ["atlantic_ocean_mask", "pacific_ocean_mask", "southern_ocean_mask", "indian_ocean_mask"];
    
    def __call__(self, data):
        function = "(process_indicator_layers.py: "+self.name+".__call__)";
        try:
            atlantic_ocean_mask = data["atlantic_ocean_mask"].fdata;
            pacific_ocean_mask = data["pacific_ocean_mask"].fdata;
            southern_ocean_mask = data["southern_ocean_mask"].fdata;
            indian_ocean_mask = data["indian_ocean_mask"].fdata;
        except KeyError as e:
            print("%s: Required data layer for process indicator layer was not found." % function);
            print(type(e), e.args);
            return False;
        
        #re-assinging values and adding in missing_value entries
        for i in arange(len(atlantic_ocean_mask)):
            if atlantic_ocean_mask[i] != DataLayer.missing_value:
                if atlantic_ocean_mask[i] == 30.0:     
                    atlantic_ocean_mask[i] = 1.0
                else:
                    atlantic_ocean_mask[i] = DataLayer.missing_value
            
            if pacific_ocean_mask[i] != DataLayer.missing_value:
                if pacific_ocean_mask[i] == 70.0:  
                    pacific_ocean_mask[i] = 1.0
                else:
                    pacific_ocean_mask[i] = DataLayer.missing_value
                  
            if southern_ocean_mask[i] != DataLayer.missing_value:
                if southern_ocean_mask[i] == 90.0:     
                    southern_ocean_mask[i] = 1.0
                else:
                    southern_ocean_mask[i] = DataLayer.missing_value
                  
            if indian_ocean_mask[i] != DataLayer.missing_value:
                if indian_ocean_mask[i] == 50.0:   
                    indian_ocean_mask[i] = 1.0
                else:
                    indian_ocean_mask[i] = DataLayer.missing_value
        
        return True;
    

#longhurst provinces
#TODO: need the full longhurst mask file, currently only 3 provinces
class longhurst_provinces_indicator(ProcessIndicatorBase):
    #differenceThreshold is the degree warmer sstskin must be than sstfnd to be considered a region of diurnal warming
    def __init__(self, differenceThreshold=0.5):
        self.name = self.__class__.__name__;
        self.differenceThreshold = differenceThreshold;
    
    def input_names(self):
        return ["longhurst_mask"];
    def output_names(self):
        return ["longhurst_mask"];
    
    def __call__(self, data):
        function = "(process_indicator_layers.py: "+self.name+".__call__)";
        try:
            longhurst_mask = data["longhurst_mask"].fdata;
        except KeyError as e:
            print("%s: Required data layer for process indicator layer was not found." % function);
            print(type(e), e.args);
            return False;

            #reassign longhurst provinces
            for i in arange(len(longhurst_mask)):
               if longhurst_mask[i] != DataLayer.missing_value:
                   if  longhurst_mask[i] == 0.0:
                       longhurst_mask[i] = DataLayer.missing_value_int
                   elif longhurst_mask[i] == 30.0:
                       longhurst_mask[i] = 1
               else:
                   longhurst_mask[i] = DataLayer.missing_value_int

        return True;
