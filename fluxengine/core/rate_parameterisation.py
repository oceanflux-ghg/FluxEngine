#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 16 14:54:46 2018

@author: tomholding
"""

#TODO: in fe, select and pack/unpack the appropriate data layers so that k_parameterisation module does not have full access to all data.

from .datalayer import DataLayer;
from numpy import arange, array;
from math import sqrt, exp;
import inspect;

###Methods
def GM12_kd_wind(windu10_fdata, windu10_moment2_fdata, windu10_moment3_fdata, scskin_fdata, nx, ny):
    #kd gas transfer from Goddijn-Murphy et al., JGR 2012
   
    kdwind_fdata = array([DataLayer.missing_value] * nx * ny)
    
    for i in arange(nx * ny):   
       kdwind_fdata[i] = DataLayer.missing_value
       if ( (windu10_fdata[i] != DataLayer.missing_value) and (windu10_moment2_fdata[i] != DataLayer.missing_value) and (windu10_moment3_fdata[i] != DataLayer.missing_value) and (scskin_fdata[i] != DataLayer.missing_value) and (scskin_fdata[i] > 0.0) ):
          kdwind_fdata[i] = (2.2 * windu10_fdata[i]) - 3.4
          kdwind_fdata[i] = kdwind_fdata[i] * sqrt(660.0/scskin_fdata[i])              
       else:
          kdwind_fdata[i] = DataLayer.missing_value
      
    return kdwind_fdata

def OceanFluxGHG_k(sigma0_fdata, sig_wv_ht_fdata, windu10_fdata, windu10_moment2_fdata, sstskinC_fdata, pco2_sw_fdata, scskin_fdata):
   dataLength = len(sigma0_fdata);
    # determine the combined Goddijn-Murphy 2012 and Fangohr and Woolf parameterisation
   kinematic_fdata = array([DataLayer.missing_value] * dataLength)
   CD_fdata = array([DataLayer.missing_value] * dataLength)
   friction_fdata = array([DataLayer.missing_value] * dataLength)
   
   #kt_fdata = array([missing_value] * nx*ny) #TMH: This doens't appear to be used...
   kd_fdata = array([DataLayer.missing_value] * dataLength)
   kb_fdata = array([DataLayer.missing_value] * dataLength)
   
   for i in arange(dataLength):

      # kinematic viscosity
     if ( (sstskinC_fdata[i] != DataLayer.missing_value) ):        
         # kinematic viscosity
        kinematic_fdata[i] = 0.00000183 * exp( (-(sstskinC_fdata[i])) / 36.0)
     else:
        kinematic_fdata[i] = DataLayer.missing_value
     pco2_sw_fdata.shape = (dataLength)

      # wind drag coefficient
      # algorithm is only value for a wind speed of up to 26 ms^-1
     if ( (windu10_fdata[i] != DataLayer.missing_value) and (windu10_moment2_fdata[i] != DataLayer.missing_value) and (windu10_fdata[i] < 26.0) and (windu10_fdata[i] > 0.0)):
        if (windu10_fdata[i] >= 6.0):
           CD_fdata[i] = 0.60 + (0.070 * windu10_fdata[i])
        elif (windu10_fdata[i] < 6.0):
           CD_fdata[i] = 0.29 + (3.1 / windu10_fdata[i]) + ( 7.7 /(windu10_moment2_fdata[i]))
        else:
           CD_fdata[i] = DataLayer.missing_value
     else:
        CD_fdata[i] = DataLayer.missing_value
     
      # friction velocity in units of ?
     if ( (windu10_fdata[i] != DataLayer.missing_value) and (CD_fdata[i] != DataLayer.missing_value) ):
        friction_fdata[i] = sqrt(CD_fdata[i] * 0.001) * windu10_fdata[i]
     else:
        friction_fdata[i] = DataLayer.missing_value
     
      # kb component     
     if ( (sig_wv_ht_fdata[i] != DataLayer.missing_value) and (friction_fdata[i] != DataLayer.missing_value) and (kinematic_fdata[i] != DataLayer.missing_value) and  (friction_fdata[i] > 0.0) and (sig_wv_ht_fdata[i] > 0.0) ):  
         
        kb_fdata[i] = 0.00002 * ( ((sig_wv_ht_fdata[i] * 100.0) * (friction_fdata[i]) * 100 * 3600) / (kinematic_fdata[i] * 36000000.0) )
        kb_fdata[i] = kb_fdata[i]#/36.0 conversion from cm/h to 10^-4 m/s (100/3600) = 1/36
     else:
        kb_fdata[i] = DataLayer.missing_value
      
      # kd component                 
     if ( (scskin_fdata[i] != DataLayer.missing_value) and (sigma0_fdata[i] != DataLayer.missing_value) and (scskin_fdata[i] > 0.0) ):
        
         # conversion from sigma0 dB units to linear units
        sigma0_fdata[i] = sigma0_fdata[i]/10.0
        sigma0_fdata[i] = pow(10.0,sigma0_fdata[i])
        
        kd_fdata[i] = (( (2100/(sigma0_fdata[i]*sigma0_fdata[i]) ) + 0.1 ) * (pow( (scskin_fdata[i] / 600.0),-0.5))) 
         
        kd_fdata[i] = kd_fdata[i]# /36.0 # unit conversion from cm/h to 10^-4 m/s (100/3600) = 1/36
     else:
        kd_fdata[i] = DataLayer.missing_value
                
   return kd_fdata, kb_fdata

# kb_weighting and kd_weighting: Weighting for kb and kd components of k_GoddijnMurphy_Fangohr2012 k parameterisation
# Setting both equal to 1.0 means that the total k will simply be a linear combination
# These need to both be valid real numbers
def OceanFluxGHG_kt(kd_fdata, kb_fdata, kb_weighting, kd_weighting):
   #combining the Oceanflux kd and kb components
   dataLength = len(kd_fdata);
   ktotal_fdata = array([DataLayer.missing_value] * dataLength)
   for i in arange(dataLength):  
       # summing the results
       # units are in 10^-4 m/s
      if ( (kd_fdata[i] != DataLayer.missing_value) and (kb_fdata[i] != DataLayer.missing_value) ):
          # optional weighting in here
         ktotal_fdata[i] = (kb_weighting*kb_fdata[i]) + (kd_weighting*kd_fdata[i])
      elif ( (kd_fdata[i] != DataLayer.missing_value) and (kb_fdata[i] == DataLayer.missing_value) ): 
         ktotal_fdata[i] = kd_fdata[i]
      elif ( (kd_fdata[i] == DataLayer.missing_value) and (kb_fdata[i] != DataLayer.missing_value) ): 
         ktotal_fdata[i] = DataLayer.missing_value
      elif ( (kd_fdata[i] == DataLayer.missing_value) and (kb_fdata[i] == DataLayer.missing_value) ):
         ktotal_fdata[i] = DataLayer.missing_value
      else: #not possible to get into this state, but included for completeness
         ktotal_fdata[i] = DataLayer.missing_value

   return ktotal_fdata


class KCalculationBase:
    def __init__(self):
        self.name = self.__class__.__name__;
    
    def input_names(self):
        raise NotImplementedError("input_names must be implemented in all derived classes.");
    
    def output_names(self):
        raise NotImplementedError("output_names must be implementedin all derived classes.");
    
    #This should take 
    def __call__(self, data):
        raise NotImplementedError("__call__ must be implemented in all derived classes.");


###K parameterisation functors (these performs the k calculation and encapsulate prerequisites etc.)
#Example showing required methods and format.
#Class name should correspond to a valid value of 'k_parameterisation' in the config file, e.g. "k_Wanninkhof1992".
class k_example(KCalculationBase):
    #Optional initialiser arguments can be used, but their names must correspond to names in the config file.
    #For example example_init_parameter is used here.
    def __init__(self, example_init_parameter):
        self.name = self.__class__.__name__;
        self.parameter = example_init_parameter; #'example_init_parameter' would need to be defined in the configuration file
    
    #Must return a list of strings corresponding to the input data layers required by the k calculations. These must already exist.
    def input_names(self):
        return ["list", "of", "DataLayer", "names"];
    
    #Must return a list of strings corresponding to the names of data layers which it will write to. These may or may not already exist.
    def output_names(self):
        return ["these", "are", "written", "to"];
    
    #Main k calculation. Input and output datalayers can be extracted from 'data'.
    #Should modify output layers in place (i.e. without copying), and return True or False to indicate successful execution.
    def __call__(self, data):
        print(self.name, "with example_init_parameter as", self.parameter);
        return True;


#rain wet deposition
class rain_wet_deposition(KCalculationBase):
    def __init__(self):
        self.name = self.__class__.__name__;
        
    def input_names(self):
        return [];
    
    def output_names(self):
        return ["k"];
        
    def __call__(self, data):
        function = "(rate_parameterisation.py: rain_wet_deposition.__call__)"
        print("%s Using the rain_wet_deposition k parameterisation" % (function));
        try:
            k = data["k"].fdata;
            data["k"].standardName = "wet_deposition so no k parameterisation";
            data["k"].longName = "wet deposition of DIC by rain so no k parameterisation";
        except KeyError as e:
            print("%s: Required data layer for selected k parameterisation was not found." % function);
            print(type(e), e.args);
            return False;

        #setting all k values to 0.0
        k[:] = 0.0;
        return 0;


#Ho et al., 2006
#Ho, D.T., Law, C.S., Smith, M.J., Schlosser, P., Harvey, M. and Hill, P., 2006. Measurements of air‐sea gas exchange at high wind speeds in the Southern Ocean: Implications for global parameterizations. Geophysical Research Letters, 33(16).
class k_Ho2006(KCalculationBase):
    def __init__(self):
        self.name = self.__class__.__name__;
    
    def input_names(self):
        return ["windu10", "windu10_moment2", "windu10_moment3", "scskin"];
    
    def output_names(self):
        return ["k"];
    
    def __call__(self, data):
        function = "(rate_parameterisation.py: k_Ho2006.__call__)"
        print("%s Using the Ho et al., 2006 (H06) k parameterisation (default option)" % (function));
        try:
            #for ease of access, simply assign attributes to each input/output.
            for name in self.input_names() + self.output_names():
                setattr(self, name, data[name].fdata);
            data["k"].standardName="gas_transfer_velocity_of_carbon_dioxide";
            data["k"].longName="Ho et al., 2006 (H06) gas transfer velocity";
            
        except KeyError as e:
            print("%s: Required data layer for selected k parameterisation was not found." % function);
            print(type(e), e.args);
            return False;
        
        #determine the Ho et al 2006 k relationship
        for i in arange(len(self.k)):
            self.k[i] = DataLayer.missing_value;
            if ((self.windu10[i] != DataLayer.missing_value) and (self.windu10_moment2[i] != DataLayer.missing_value) and (self.windu10_moment3[i] != DataLayer.missing_value) and (self.scskin[i] != DataLayer.missing_value) and (self.scskin[i] > 0.0) ):
                self.k[i] = 0.266 * self.windu10_moment2[i];
                self.k[i] = self.k[i] * sqrt(600.0/self.scskin[i]);
                #self.k[i] = self.k[i]; #/36.0 # conversion from cm/h to 10^-4 m/s (100/3600) = 1/36
            else:
                self.k[i] = DataLayer.missing_value
        return True;


#Nightingale, P. D., et al. 2000. In situ evaluation of air-sea gas exchange parameterizations using novel conserva-tive  and  volatile  tracers.  Global  Biogeochem.  Cycles14:373-387 [doi:10. 1029/ 1999GB900091].
class k_Nightingale2000(KCalculationBase):
    def __init__(self):
        self.name = self.__class__.__name__;
    
    def input_names(self):
        return ["windu10", "windu10_moment2", "scskin"];
        
    def output_names(self):
        return ["k"];
    
    def __call__(self, data):
        function = "(rate_parameterisation.py: k_Nightingale2000.__call__)";
        print("%s Using the Nightingale et al., 2000 (N00) k parameterisation" % (function))
        try:
            #for ease of access, simply assign attributes to each input/output.
            for name in self.input_names() + self.output_names():
                setattr(self, name, data[name].fdata);
            data["k"].standardName="gas_transfer_velocity_of_carbon_dioxide" ;
            data["k"].longName="Nightingale et al., 2000 (N00) gas transfer velocity";
        except KeyError as e:
            print("%s: Required data layer for selected k parameterisation was not found." % function);
            print(type(e), e.args);
            return False;
        
        #determine the Nightingale et al 2000 k relationship
        for i in arange(len(self.k)):   
            self.k[i] = DataLayer.missing_value
            if ( (self.windu10[i] != DataLayer.missing_value) and (self.windu10_moment2[i] != DataLayer.missing_value) and (self.scskin[i] != DataLayer.missing_value) and (self.scskin[i] > 0.0) ):#SOCATv4 - No need for wind moment3
                self.k[i] = 0.222 * self.windu10_moment2[i] + (0.333 * self.windu10[i])  
                self.k[i] = self.k[i] * sqrt(600.0/self.scskin[i])
            else:
                self.k[i] = DataLayer.missing_value
        
        return True;



class kt_OceanFluxGHG(KCalculationBase):
    def __init__(self, kb_weighting, kd_weighting):
        self.name = self.__class__.__name__;
        self.kb_weighting = kb_weighting;
        self.kd_weighting = kd_weighting;
    
    def input_names(self):
        return ["sigma0", "sig_wv_ht", "windu10", "windu10_moment2", "sstskinC", "pco2_sw", "scskin"];
    
    def output_names(self):
        return ["k", "kb", "kd", "kt"];
    
    def __call__(self, data):
        #using OceanFlux GHG kt approach
        function = "(rate_parameterisation.py: kt_OceanFluxGHG.__call__)";
        print("%s Using the OceanFluxGHG kt parameterisation (kt = kd_backscatter + kb)" % (function))
        
        try:
            #for ease of access, simply assign attributes to each input/output.
            for name in self.input_names() + self.output_names():
                setattr(self, name, data[name].fdata);
            data["k"].standardName = "total (direct_from_backscatter plus bubble mediated) component of gas transfer velocity of carbon dioxide";
            data["k"].longName = "total (direct_from_backscatter plus bubble mediated) component of gas transfer velocity of carbon dioxide";
        except KeyError as e:
            print("%s: Required data layer for selected k parameterisation was not found." % function);
            print(type(e), e.args);
            return False;
        
        self.kd, self.kb = OceanFluxGHG_k(self.sigma0, self.sig_wv_ht, self.windu10, self.windu10_moment2, self.sstskinC, self.pco2_sw, self.scskin);
        # calculate the total kt
        self.kt = OceanFluxGHG_kt(self.kd, self.kb, self.kb_weighting, self.kd_weighting);
        self.k[:] = self.kt;
        
        return True;


class k_Wanninkhof1992(KCalculationBase):
    def __init__(self):
        self.name = self.__class__.__name__;
    
    def input_names(self):
        return ["windu10", "windu10_moment2", "scskin"];
    
    def output_names(self):
        return ["k"];
    
    def __call__(self, data):
        function = "(rate_parameterisation.py: k_Wanninkhof1992.__call__)";
        print("%s Using the Wanninkhof 1992 (W92) k parameterisation" % (function))
        try:
            #for ease of access, simply assign attributes to each input/output.
            for name in self.input_names() + self.output_names():
                setattr(self, name, data[name].fdata);
            data["k"].standardName="gas_transfer_velocity_of_carbon_dioxide";
            data["k"].longName="Wanninkhof, 1992 (W92) gas transfer velocity";
        except KeyError as e:
            print("%s: Required data layer for selected k parameterisation was not found." % function);
            print(type(e), e.args);
            return False;
        
        #determine the Wanninkhof 1992 k relationship
        for i in arange(len(self.k)):   
           self.k[i] = DataLayer.missing_value
          #if ( (self.windu10[i] != DataLayer.missing_value) and (self.windu10_moment2[i] != DataLayer.missing_value) and (self.windu10_moment3[i] != DataLayer.missing_value) and (self.scskin[i] != DataLayer.missing_value) and (self.scskin[i] > 0.0) ):
           if ( (self.windu10[i] != DataLayer.missing_value) and (self.windu10_moment2[i] != DataLayer.missing_value) and (self.scskin[i] != DataLayer.missing_value) and (self.scskin[i] > 0.0) ):#SOCATv4 - no need for wind moment 3
              self.k[i] = 0.31 * self.windu10_moment2[i]    
              self.k[i] = self.k[i] * sqrt(660.0/self.scskin[i])
               
              self.k[i] = self.k[i]# /36.0 # conversion from cm/h to 10^-4 m/s (100/3600) = 1/36
           else:
              self.k[i] = DataLayer.missing_value
        return True;


class k_McGillis2001(KCalculationBase):
    def __init__(self):
        self.name = self.__class__.__name__;
    
    def input_names(self):
        return ["windu10", "windu10_moment2", "windu10_moment3", "scskin"];
    
    def output_names(self):
        return ["k"];
    
    def __call__(self, data):
        function = "(rate_parameterisation.py: k_McGillis2001.__call__)";
        print("%s Using the McGillis et al., 2001 (M01) k parameterisation" % (function))
        try:
            #for ease of access, simply assign attributes to each input/output.
            for name in self.input_names() + self.output_names():
                setattr(self, name, data[name].fdata);
            data["k"].standardName="gas_transfer_velocity_of_carbon_dioxide";
            data["k"].longName="McGillis et al., 2001 (M01) gas transfer velocity";
        except KeyError as e:
            print("%s: Required data layer for selected k parameterisation was not found." % function);
            print(type(e), e.args);
            return False;
        
        #determine the Wanninkhof and McGillis 1999 k relationship
        for i in arange(len(self.k)):   
            self.k[i] = DataLayer.missing_value
            if ( (self.windu10[i] != DataLayer.missing_value) and (self.windu10_moment2[i] != DataLayer.missing_value) and (self.windu10_moment3[i] != DataLayer.missing_value) and (self.scskin[i] != DataLayer.missing_value) and (self.scskin[i] > 0.0) ):
                self.k[i] = 3.3 + (0.026 * self.windu10_moment3[i])
                self.k[i] = self.k[i] * sqrt((660.0/self.scskin[i]))
                     
                self.k[i] = self.k[i]# /36.0 # conversion from cm/h to 10^-4 m/s (100/3600) = 1/36.0
            else:
                self.k[i] = DataLayer.missing_value 
        return True;


class k_Ho1997(KCalculationBase):
    def __init__(self):
        self.name = self.__class__.__name__;
    
    def input_names(self):
        return ["windu10", "windu10_moment2", "windu10_moment3", "scskin", "rain"];
        
    def output_names(self):
        return ["k", "krain"];
    
    def __call__(self, data):
        function = "(rate_parameterisation.py: k_Ho1997.__call__)";
        print("%s Using the Ho1997 et al., 2000 (N00) k parameterisation" % (function))
        try:
            #for ease of access, simply assign attributes to each input/output.
            for name in self.input_names() + self.output_names():
                setattr(self, name, data[name].fdata);
            data["k"].standardName="gas_transfer_velocity_of_carbon_dioxide_due_to_rain";
            data["k"].longName="Ho et al., 1997 (H97) gas transfer velocity";
        except KeyError as e:
            print("%s: Required data layer for selected k parameterisation was not found." % function);
            print(type(e), e.args);
            return False;

        print("%s Using the Ho et al., 1997 (H97) k parameterisation (rain driven k)" % (function))
        #Ho et al, Tellus, 1997 rain component of k
        #note based on rain rate (Rn), but filtered based on wind component so still looking at the same datapoints
        for i in arange(len(self.k)):   
            self.k[i] = DataLayer.missing_value
            if ( (self.windu10[i] != DataLayer.missing_value) and (self.windu10_moment2[i] != DataLayer.missing_value) and (self.windu10_moment3[i] != DataLayer.missing_value) and (self.scskin[i] != DataLayer.missing_value) and (self.scskin[i] > 0.0) and (self.rain[i] != DataLayer.missing_value) ):
                self.k[i] = 0.929 + (0.679 * self.rain[i]) - (0.0015*pow(self.rain[i], 2.0))
                self.k[i] = self.k[i] * sqrt(600.0/self.scskin[i])

                #copying the result into the krain array
                #k_fdata[i] = k_fdata[i]# /36.0 # conversion from cm/h to 10^-4 m/s (100/3600) = 1/36.0
                self.krain[i] = self.k[i]
            else:
                self.k[i] = DataLayer.missing_value
                self.krain[i] = DataLayer.missing_value
        return True;


#Goddijn-Murphy, L., D.K. Woolf, B. Chapron, P. Queffeulou (2013) Improvements to estimating the air–sea gas transfer velocity by using dual-frequency, altimeter backscatter, Remote Sensing of Environment,  139, 1-5, doi:10.1016/j.rse.2013.07.026
class kd_OceanFluxGHG_backscatter(KCalculationBase):
    def __init__(self):
        self.name = self.__class__.__name__;
    
    def input_names(self):
        return ["sigma0", "sig_wv_ht", "windu10", "windu10_moment2", "sstskinC", "pco2_sw", "scskin"];
        
    def output_names(self):
        return ["k", "kd"];
    
    def __call__(self, data):
        function = "(rate_parameterisation.py: kd_OceanFluxGHG_backscatter.__call__)";
        print("%s Using the OceanFluxGHG kd-backscatter (direct component) parameterisation" % (function));
        
        try:
            #for ease of access, simply assign attributes to each input/output.
            for name in self.input_names() + self.output_names():
                setattr(self, name, data[name].fdata);
            data["k"].standardName="direct_gas_transfer_velocity_of_carbon_dioxide_from_backscatter";
            data["k"].longName="direct component of gas transfer velocity of carbon dioxide derived from radar backscatter";
        except KeyError as e:
            print("%s: Required data layer for selected k parameterisation was not found." % function);
            print(type(e), e.args);
            return False;
        
        #direct component of k using OceanFluxGHG radar backscatter parameterisation
        self.kd, kb_discard = OceanFluxGHG_k(self.sigma0, self.sig_wv_ht, self.windu10, self.windu10_moment2, self.sstskinC, self.pco2_sw, self.scskin);
        self.k[:] = self.kd #disregards kb
        return True;
       

class kb_OceanFluxGHG(KCalculationBase):
    def __init__(self):
        self.name = self.__class__.__name__;
    
    def input_names(self):
        return ["sigma0", "sig_wv_ht", "windu10", "windu10_moment2", "sstskinC", "pco2_sw", "scskin"];
    
    def output_names(self):
        return ["k", "kb"];
    
    def __call__(self, data):
        function = "(rate_parameterisation.py: kb_OceanFluxGHG.__call__)";
        print("%s Using the OceanFluxGHG kb (bubble mediated) parameterisation" % (function));
        
        try:
            #for ease of access, simply assign attributes to each input/output.
            for name in self.input_names() + self.output_names():
                setattr(self, name, data[name].fdata);
            data["k"].standardName="bubble_gas_transfer_velocity_of_carbon_dioxide";
            data["k"].longName="bubble mediated component of gas transfer velocity of carbon dioxide";
        except KeyError as e:
            print("%s: Required data layer for selected k parameterisation was not found." % function);
            print(type(e), e.args);
            return False;

        #bubble mediated component of k using OceanFluxGHG parameterisation   
        self.kd, self.kb = OceanFluxGHG_k(self.sigma0, self.sig_wv_ht, self.windu10, self.windu10_moment2, self.sstskinC, self.pco2_sw, self.scskin);
        self.k[:] = self.kb #disregards kd
        return True;


class k_generic(KCalculationBase):
    def __init__(self, k_generic_a0, k_generic_a1, k_generic_a2, k_generic_a3, k_generic_sc):
        function = inspect.stack()[0][1]+", "+inspect.stack()[0][3];
        self.name = self.__class__.__name__;
        self.k_generic_a0 = k_generic_a0;
        self.k_generic_a1 = k_generic_a1;
        self.k_generic_a2 = k_generic_a2;
        self.k_generic_a3 = k_generic_a3;
        self.k_generic_sc = k_generic_sc;
        
        if self.k_generic_sc != 600 and self.k_generic_sc != 660:
            raise ValueError("%s: Could not initialise k_generic rate parameterisation. Only schmidt numbers (k_generic_sc) of 600 and 660 are allowed, but a schmidt number of %d was provided." % (function, k_generic_sc));
    
    def input_names(self):
        return ["windu10", "windu10_moment2", "windu10_moment3", "scskin"];
        
    def output_names(self):
        return ["k"];
    
    def __call__(self, data):
        function = inspect.stack()[0][1]+", "+inspect.stack()[0][3];
        print("%s Using the general cubic form of k with user specified values (k_generic_a0:%lf k_generic_a1:%lf k_generic_a2:%lf k_generic_a3:%lf k_generic_sc:%lf)" % (function, self.k_generic_a0, self.k_generic_a1, self.k_generic_a2, self.k_generic_a3, self.k_generic_sc))
        
        try:
            #for ease of access, simply assign attributes to each input/output.
            for name in self.input_names() + self.output_names():
                setattr(self, name, data[name].fdata);
            data["k"].standardName="generic_gas_transfer_velocity_formulation";
            data["k"].longName="User defined generic formulation of k (a0:%lf a1:%lf a2:%lf a3:%lf sc:%lf)" % (self.k_generic_a0, self.k_generic_a1, self.k_generic_a2, self.k_generic_a3, self.k_generic_sc);
        except KeyError as e:
            print("%s: Required data layer for selected k parameterisation was not found." % function);
            print(type(e), e.args);
            return False;
        
        # general form
        # kw = (600/Sc)^0.5 [a0 + a1*U + a2*U2 + a3*U3]        
        for i in arange(len(self.k)):
            self.k[i] = DataLayer.missing_value
            #if ( (self.windu10[i] != DataLayer.missing_value) and (self.windu10_moment2[i] != DataLayer.missing_value) and (self.windu10_moment3[i] != DataLayer.missing_value) and (self.scskin[i] != DataLayer.missing_value) and (self.scskin[i] > 0.0) ):
            if ( (self.windu10[i] != DataLayer.missing_value) and (self.windu10_moment2[i] != DataLayer.missing_value) and (self.scskin[i] != DataLayer.missing_value) and (self.scskin[i] > 0.0) ):
                
                #general form kw = (600/Sc)^0.5 [a0 + a1*U + a2*U2 + a3*U3]
                self.k[i] = self.k_generic_a0 + (self.k_generic_a1 * self.windu10[i]) +  (self.k_generic_a2 * self.windu10_moment2[i]) + (self.k_generic_a3 * self.windu10_moment3[i])
                self.k[i] = self.k[i] * sqrt((self.k_generic_sc/self.scskin[i]))
                #k_fdata[i] = k_fdata[i]# /36.0 # conversion from cm/h to 10^-4 m/s (100/3600) = 1/36.0
            else:
                self.k[i] = DataLayer.missing_value  
        
        return True;


#direct component using the Goddijn-Murphy et al., 2012 wind parameterisation for kd
class kd_OceanFluxGHG_wind(KCalculationBase):
    def __init__(self):
        self.name = self.__class__.__name__;
    
    def input_names(self):
        return ["windu10", "windu10_moment2", "windu10_moment3", "scskin"];
    
    def output_names(self):
        return ["k", "kd"];
    
    def __call__(self, data):
        function = "(rate_parameterisation.py: kd_OceanFluxGHG_wind.__call__)";
        print("%s Using the Goddijn-Murphy et al., 2012 kd-wind (direct-component) parameterisation" % (function));
        
        try:
             #for ease of access, simply assign attributes to each input/output.
            for name in self.input_names() + self.output_names():
                setattr(self, name, data[name].fdata);
            data["k"].standardName="direct_gas_transfer_velocity_of_carbon_dioxide_from_wind_speed";
            data["k"].longName="direct component of gas transfer velocity of carbon dioxide derived from wind speed (Goddijn-Murphy et al., 2012)";
        except KeyError as e:
            print("%s: Required data layer for selected k parameterisation was not found." % function);
            print(type(e), e.args);
            return False;
    
        #direct component of k using OceanFluxGHG kd-wind parameterisation
        self.kd = GM12_kd_wind(self.windu10, self.windu10_moment2, self.windu10_moment3, self.scskin, self.nx, self.ny);
        self.k[:] = self.kd

        return True;


#Goddijn-Murphy et. al. 2015
class kt_OceanFluxGHG_kd_wind(KCalculationBase):
    def __init__(self, kb_weighting, kd_weighting):
        self.name = self.__class__.__name__;
        self.kb_weighting = kb_weighting;
        self.kd_weighting = kd_weighting;
    
    def input_names(self):
        return ["sigma0", "sig_wv_ht", "windu10", "windu10_moment2", "windu10_moment3", "sstskinC", "pco2_sw", "scskin"];
    
    def output_names(self):
        return ["k", "kb", "kd", "kt"];
    
    def __call__(self, data):
        # using OceanFlux GHG kt approach with kd based on the wind parameterisation of Goddijn-Murphy et al., 2015
        function = "(rate_parameterisation.py: kt_OceanFluxGHG_kd_wind.__call__)";
        print("%s Using the OceanFluxGHG kt parameterisation (kt = kd_wind_speed + kb)" % (function));
        
        try:
            #for ease of access, simply assign attributes to each input/output.
            for name in self.input_names() + self.output_names():
                setattr(self, name, data[name].fdata);
            data["k"].standardName = "total_gas_transfer_velocity_of_carbon_dioxide";
            data["k"].longName="total (direct_from_wind_speed plus bubble mediated) component of gas transfer velocity of carbon dioxide";
        except KeyError as e:
            print("%s: Required data layer for selected k parameterisation was not found." % function);
            print(type(e), e.args);
            return False;
        
        #kb portion will be discarded
        self.kd, self.kb = OceanFluxGHG_k(self.sigma0, self.sig_wv_ht, self.windu10, self.windu10_moment2, self.sstskinC, self.pco2_sw, self.scskin);
        
        #overwrite kd with Goddijn-Murphy et al., JGR 2012 gas transfer...
        self.kd = GM12_kd_wind(self.windu10, self.windu10_moment2, self.windu10_moment3, self.scskin, self.nx, self.ny);
        self.kt = OceanFluxGHG_kt(self.kd, self.kb, self.kb_weighting, self.kd_weighting);
        self.k[:] = self.kt
        return True;



class k_Wanninkhof2014(KCalculationBase):
    #def __init__(self, kb_weighting, kd_weighting):
    def __init__(self):
        self.name = self.__class__.__name__;
    
    def input_names(self):
        return ["windu10", "windu10_moment2", "scskin"];
    
    def output_names(self):
        return ["k"];
    
    def __call__(self, data):
        # using OceanFlux GHG kt approach with kd based on Wanninkhof2014
        # Wanninkhof, Rik. "Relationship between wind speed and gas exchange over the ocean revisited." Limnology and Oceanography: Methods 12.6 (2014): 351-362.
        function = "(rate_parameterisation.py: k_Wanninkhof2014.__call__)";
        print("%s Using the Wanninkhof 2014 k parameterisation" % (function));
        
        try:
            #for ease of access, simply assign attributes to each input/output.
            for name in self.input_names() + self.output_names():
                setattr(self, name, data[name].fdata);
            data["k"].standardName="W14" 
            data["k"].longName="Wanninkhof 2014, Limnol. Oceanogr.: Methods 12, 2014, 351-362"#IGA
        except KeyError as e:
           print("%s: Required data layer for selected k parameterisation was not found." % function);
           print(type(e), e.args);
           return False;
        
        #determine the k relationship
        for i in arange(len(self.k)):   
            self.k[i] = DataLayer.missing_value
            #if ( (self.windu10[i] != DataLayer.missing_value) and (self.windu10_moment2[i] != DataLayer.missing_value) and (self.windu10_moment3[i] != DataLayer.missing_value) and (self.scskin[i] != DataLayer.missing_value) and (self.scskin[i] > 0.0) ):
            if ( (self.windu10[i] != DataLayer.missing_value) and (self.windu10_moment2[i] != DataLayer.missing_value) and (self.scskin[i] != DataLayer.missing_value) and (self.scskin[i] > 0.0) ):#SOCATv4 - no need for wind moment 3
               self.k[i] = 0.251 * self.windu10_moment2[i]
               self.k[i] = self.k[i] * sqrt(660.0/self.scskin[i])
               
               #self.k[i] = self.k[i]# /36.0 # conversion from cm/h to 10^-4 m/s (100/3600) = 1/36
            else:
               self.k[i] = DataLayer.missing_value
        return True;




#Biological surfactant dampening with Nightingale2000 parameterisation. First applies Nightingale parameterisation then applied surfactant suppression
#Implements method described in:
#   Pereira, Ryan, Ian Ashton, Bita Sabbaghzadeh, Jamie D. Shutler, and Robert C. Upstill-Goddard. "Reduced air–sea CO 2 exchange in the Atlantic Ocean due to biological surfactants." Nature Geoscience (2018): 1.
class k_Nightingale2000_with_surfactant_suppression(KCalculationBase):
    def __init__(self):
        self.name = self.__class__.__name__;
    
    def input_names(self):
        return ["windu10", "windu10_moment2", "windu10_moment3", "scskin", "sstskin"];
        
    def output_names(self):
        return ["k"];
    
    def __call__(self, data):
        function = "(rate_parameterisation.py: k_Nightingale2000_with_surfactant_suppression.__call__)";
        print("%s Using the Nightingale 2000 with surfactant suppression k parameterisation" % (function));
        
        #Calculate Nightingale2000 k
        nightingaleFunctor = k_Nightingale2000();
        nightingaleFunctor(data);
        
        #Apply surfactant suppression using methodology in Pereira et al. 2018.
        try:
            #for ease of access, simply assign attributes to each input/output.
            self.SSTSkin = data["sstskin"].fdata;
            data["k"].standardName="Nightingale2000 with surfactant surpression" 
            data["k"].longName="Nightingale 2000 with surfactant surpression: Pereira, R. et al., Nature Geoscience";
            self.k = data["k"].fdata; #Contains Nightingale parameterisation
        except KeyError as e:
           print("%s: Required data layer for selected k parameterisation was not found." % function);
           print(type(e), e.args);
           return False;
        
        #Apply surfactant suppression
        for i in arange(len(self.k)):   
            if ( (self.SSTSkin[i] != DataLayer.missing_value) and (self.k[i] != DataLayer.missing_value) ):
               self.k[i] = self.k[i] * (1.0 - (0.0046 * (self.SSTSkin[i]-273.15)**2.5673)/100.0);
               #percentage surpression, so (1.0 - ((0.0046 * (self.SSTSkin[i]-273.15)**2.5673))/100.0))
        return True;


#Zappa et al 2007. k-parameterisation using dissipation rate of turbulent kinetic energy (epsilon)
#Zappa, Christopher J., Wade R. McGillis, Peter A. Raymond, James B. Edson, Eric J. Hintsa, Hendrik J. Zemmelink, John WH Dacey, and David T. Ho. "Environmental turbulent mixing controls on air‐water gas exchange in marine and aquatic systems." Geophysical Research Letters 34, no. 10 (2007).
#Uses a custom configuration file parameter 'k_Zappa2007_epsilon_calibration' which 
class k_Zappa2007(KCalculationBase):
    def __init__(self, k_Zappa2007_epsilon_calibration=1.0):
        self.name = self.__class__.__name__;
        self.k_Zappa2007_epsilon_calibration = k_Zappa2007_epsilon_calibration; #Allows calibration e.g. to different depths.
    
    def input_names(self):
        return ["tke_dissipation", "scskin", "sstskin"];
    
    def output_names(self):
        return ["k"];
    
    def si_viscosity_to_cSt(self, viscosity):
        return viscosity*1000000.0;
    
    def __call__(self, data):
        # using OceanFlux GHG kt approach with kd based on Wanninkhof2014
        # Wanninkhof, Rik. "Relationship between wind speed and gas exchange over the ocean revisited." Limnology and Oceanography: Methods 12.6 (2014): 351-362.
        function = "(rate_parameterisation.py: k_Zappa2007.__call__)";
        print("%s Using the Zappa 2007 k parameterisation" % (function));
        
        try:
            #for ease of access, simply assign attributes to each input/output.
            for name in self.input_names() + self.output_names():
                setattr(self, name, data[name].fdata);
            data["k"].standardName="Zappa2007" 
            data["k"].longName="Zappa 2007: using the dissipation rate of turbulent kinetic energy. Zappa, Christopher J., Wade R. McGillis, Peter A. Raymond, James B. Edson, Eric J. Hintsa, Hendrik J. Zemmelink, John WH Dacey, and David T. Ho. Environmental turbulent mixing controls on air‐water gas exchange in marine and aquatic systems. Geophysical Research Letters 34, no. 10 (2007). The dissipation rate of turbulent kinetic energy was scaled by calibration factor: "+str(self.k_Zappa2007_epsilon_calibration);
        except KeyError as e:
           print("%s: Required data layer for selected k parameterisation was not found." % function);
           print(type(e), e.args);
           return False;
        
        #determine the k relationship
        for i in arange(len(self.k)):   
            if ( (self.tke_dissipation[i] != DataLayer.missing_value) and (self.scskin[i] != DataLayer.missing_value) and (self.scskin[i] > 0.0) ):#SOCATv4 - no need for wind moment 3
                kinematicViscosity = 0.00000183 * exp( -(self.sstskin[i]-273.15) / 36.0); #Parameterisation from previously used study (GHGOceanFlux, see above k parameterisations)
                #kinematicViscosity = self.si_viscosity_to_cSt(kinematicViscosity); #Convert from m^2 s^-1 to cSt (centistokes)
                
                self.k[i] = self.k_Zappa2007_epsilon_calibration * (0.419*(self.scskin[i]**-0.5)) * ((self.tke_dissipation[i]*kinematicViscosity)**0.25); #new
                self.k[i] *= 100.0 * 3600.0; #Convert from m s^-1 to cm hr^-1
            else:
               self.k[i] = DataLayer.missing_value
        return True;






#Base class for extensions to the core k-calculation, e.g. 
class KCalculationExtension:
    def __init__(self):
        self.name = self.__class__.__name__;
    
    def input_names(self):
        raise NotImplementedError("input_names must be implemented in all derived classes.");
    
    def output_names(self):
        raise NotImplementedError("output_names must be implementedin all derived classes.");
    
    #This should take 
    def __call__(self, data):
        raise NotImplementedError("__call__ must be implemented in all derived classes.");


#linear additive k/krain for rain case
#add ho1997 value to data from chosen k parameterisation (adds rain component to the results from the existing choice of parameterisation)
#see Ashton2016
#note: data are filtered based on windu10 to ensure that they cover the same data as the wind based paramterisations 
class AddKRainLinearHo1997(KCalculationExtension):
    def __init__(self):
        self.name = self.__class__.__name__;
    
    def input_names(self):
        return ["k", "scskin", "rain", "windu10", "windu10_moment2", "windu10_moment3"];
    
    def output_names(self):
        return ["k", "krain"];
    
    def __call__(self, data):
        function = "(rate_parameterisation.py, "+self.name+".__call__)";
        print("Adding linear rain component to k parameterisation (Ashton2016)");
        try:
            #for ease of access, simply assign attributes to each input/output.
            for name in self.input_names() + self.output_names():
                setattr(self, name, data[name].fdata);
        except KeyError as e:
           print("%s: Required data layer for selected k parameterisation was not found." % function);
           print(type(e), e.args);
           return False;
        
        for i in arange(len(self.k)):
            if ( (self.k[i] != DataLayer.missing_value) and (self.scskin[i] != DataLayer.missing_value) and (self.rain[i] != DataLayer.missing_value) and (self.scskin[i] > 0.0) and (self.windu10[i] != DataLayer.missing_value) and (self.windu10_moment2[i] != DataLayer.missing_value) and (self.windu10_moment3[i] != DataLayer.missing_value)):
                self.krain[i] = 0.929 + (0.679 * self.rain[i]) - (0.0015*pow(self.rain[i], 2.0));
                self.krain[i] = self.krain[i] * sqrt(600.0/self.scskin[i]);
                #self.krain[i] = self.krain[i]# /36.0 # conversion from cm/h to 10^-4 m/s (100/3600) = 1/36.0
                self.k[i] = self.k[i] + self.krain[i];
            else:
                self.krain[i] = DataLayer.missing_value;
        return True;



#nonlinear changes to k for rain and wind.
#Parameterisation from Harrison et al., JGR 2012, equations 11, 12, 13, 14 (see Ashton2016)
#and air_density (in kg/m3) = air_pressure (in pascals) / R(gas constant for dry air) * T(in K)
#R gas constant for dry air = 287.058 J/(kg K)
#typical air pressure is 100 000 pascals
#typical SST in kelvin in 300 K
#typical density at sea level = 1.25 kg m^-3
#need to select non-linear parameterisation and the generic k (runParams.k_parameterisation == 9)
#untested 06/01/2014
class AddKRainNonlinearHarrison2012(KCalculationExtension):
    def __init__(self):
        self.name = self.__class__.__name__;
    def input_names(self):
        return ["k", "windu10", "windu10_moment2", "windu10_moment3", "scskin", "pressure", "sstskin", "rain"];
    def output_names(self):
        return ["k", "krain"];
    
    def __call__(self, data):
        function = "(rate_parameterisation.py, "+self.name+".__call__)";
        print("Adding non-linear rain component to k parameterisation (Ashton2016)");
        try:
            #for ease of access, simply assign attributes to each input/output.
            for name in self.input_names() + self.output_names():
                setattr(self, name, data[name].fdata);
        except KeyError as e:
           print("%s: Required data layer for selected k parameterisation was not found." % function);
           print(type(e), e.args);
           return False;

        alpha_r = 0.3677 # from Harrison et al., 2012, equation 12
        beta_r = 0.0
        rho_a = 0.0
        # gas constant
        R = 287.058 # in J/(kg K)
        length = len(self.k);
        
        windstress_fdata = array([DataLayer.missing_value] * length)
        winddrag_fdata = array([DataLayer.missing_value] * length)   
         # need to calculate wind drag and then the wind stress
        for i in arange(length):
           if ( (self.windu10[i] != DataLayer.missing_value) and (self.windu10_moment2[i] != DataLayer.missing_value) ):
            # wind drag
            # Yelland and Taylor 1996
              if ( (self.windu10[i] >= 6.0) and (self.windu10[i] < 26.0) ):
                 winddrag_fdata[i] = 0.60 + (0.070 * self.windu10[i])
              elif ( (self.windu10[i] < 6.0) and (self.windu10[i] > 0.0) ):
                 winddrag_fdata[i] = 0.29 + (3.1 / self.windu10[i]) + ( 7.7 /(self.windu10_moment2[i]) )
              else:
                 winddrag_fdata[i] = -999.0
              
               # wind stress
              if (winddrag_fdata[i] == -999.0 or self.windu10_moment2[i] == -999.0):
                 windstress_fdata[i] = -999.0
              elif (winddrag_fdata[i] != 0.0 and self.windu10_moment2[i] != 0.0):
                 windstress_fdata[i] = sqrt( (winddrag_fdata[i]/1000.0) * self.windu10_moment2[i] )
                 #print "CD:%lf U^2:%lf U10:%lf stress: %lf" % (winddrag_fdata[i], windu10_moment2[i], windu10[i], windstress_fdata[i])
              else:
                 windstress_fdata[i] = 0.0
           else:
              windstress_fdata[i] = DataLayer.missing_value
              winddrag_fdata[i] = DataLayer.missing_value
        
        for i in arange(length):
            if ( (self.k[i] != DataLayer.missing_value) and (self.scskin[i] != DataLayer.missing_value) and (self.rain[i] != DataLayer.missing_value) and (self.scskin[i] > 0.0) and (self.windu10[i] != DataLayer.missing_value) and (self.windu10_moment2[i] != DataLayer.missing_value) and (self.windu10_moment3[i] != DataLayer.missing_value) and (self.pressure[i] != DataLayer.missing_value) and (self.sstskin[i] != DataLayer.missing_value) and (self.windstress_fdata[i] != DataLayer.missing_value)):
                #guts in here
                # data["pressure"].fdata are in mb, need it in Pascals, so *100 to get Pascals (ie convert from mb to P)
                rho_a = (self.pressure[i]*100.0) /(R*self.sstskin[i])    
                beta_r = (0.0112*self.rain[i]) / (rho_a*pow(windstress_fdata[i],3)) 
                #print "beta_r:%lf stress^3:%lf rho_a:%lf pres:%lf R:%lf sstskinK:%lf" % (beta_r, pow(windstress_fdata[i],3), rho_a, self.data["pressure"].fdata[i], R, self.data["sstskin"].fdata[i])
                self.krain[i] = (1-exp(-(alpha_r*beta_r))) * (63.02 * pow(0.0112*self.rain[i], 0.6242))
                self.krain[i] = self.krain[i] * sqrt(600.0/self.scskin[i])
                self.k[i] = self.k[i] + self.krain[i]
            else:
                self.krain[i] = DataLayer.missing_value
                self.k[i] = DataLayer.missing_value
        return True;