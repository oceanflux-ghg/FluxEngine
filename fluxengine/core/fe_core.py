#!/usr/bin/env python2

# ofluxghg-flux-calc.py IGA working version - now using GIT at OCEANFLUX-SHARED/workspace/i.g.c.ashton_exeter.ac.uk/FE_IGA_new
#test
# utility to load various input netcdf datasets and determine
# air-sea flux of CO2 based on parameterisations in the ESA STSE OceanFlux Greenhouse Gases Technical Specification (TS) 


# History
# date, decsription, author, email, institution
# v0 10/2012  Basic structure and output, Jamie Shutler, jams@pml.ac.uk, Plymouth Marine Laboratory.
# v1 24/10/2012  improved error handling and some method defs, Jamie Shutler, jams@pml.ac.uk, Plymouth Marine Laboratory.
# v3 03/12/2012  additional functionality for uncertainty analyses, improved error handling, jams@pml.ac.uk, Plymouth Marine Laboratory.
# v4 06/08/2013 approaching final form and almost compliant with TS, jams@pml.ac.uk, Plymouth Marine Laboratory
# v5 06/08/2014 version used for final runs and FluxEngine publication
#Further changes marked by #IGA within the text. Version control within a GIT repo

 # netcdf bits
from netCDF4 import Dataset
import sys
from math import log, exp, pow, isnan;
from numpy import size, flipud, mean, zeros, nonzero, array, resize, ma, arange, dtype, ones, meshgrid, where;
from numpy import any as npany;
from random import normalvariate
import logging;
from os import path;
from string import Template;

from .datalayer import DataLayer, DataLayerMetaData;
from .settings import Settings;
#from .debug_tools import calc_mean; #calculate mean ignoring missing values.

from datetime import timedelta, datetime;

 # debug mode switches
DEBUG = False
DEBUG_PRODUCTS = True
VERIFICATION_RUNS = False # forces flux calculations to use Takahashi SST_t as the SST dataset for the pco2 data
DEBUG_LOGGING = True; # TH: Added to make debugging easier. Sets up a logger object and writes a to the specified path (or cwd)
#workingDirectory = os.getcwd(); #Only used to make a full path for logs, every other path is suppled as an absolute filepath.
                                #Should really be passed as an argument.


 # missing value set for intermediate data sets and output dataset
missing_value = -999.0
fill_value = -999.0
missing_value_int = -999
fill_value_int = -999

# # valid data ranges for testing output products and flagging violations
# # orignally set in TS, some have been modified based on the actual input datasets
# # these are used in the NetCDF metadata and to check the valid data ranges to create the
# # OFA11 process indicator layer
#OF_min = -0.75
#OF_max = 0.6
#OK1_min = 0.0
#OK1_max = 100.0
#OK3_min = 0.0
#OK3_max = 100.0
#OKD_min = 0.0
#OKD_max = 100.0
#OKR_min = 0.0
#OKR_max = 100.0
#OKB1_min = 0.0
#OKB1_max = 100.0
#OKS1_min = 15.0
#OKS1_max = 40.0
#OKT1_min = -1.8
#OKT1_max = 30.5
#OFWR_min = -0.75
#OFWR_max = 0.6
#
# # these thresholds are in g-C m^-3 (rather than ppm as in the TS)
#OSC1_min = 0.095
#OSC1_max = 0.27
#OIC1_min = 0.11
#OIC1_max = 0.22
#
# # these thresholds are in uatm
#OBPC_min = 200
#OBPC_max = 800

# function name definition for stderr and stdout messages
function = "(ofluxghg-flux-calc, main)"



#Stores a set of parameters required to run the FluxEngine, allowing easy centralised access.
class RunParameters:
    def set_parameters(self, parameterDict):
        for key in list(vars(self).keys()): #remove all previously stored values
            delattr(self, key);
        for key in list(parameterDict.keys()): #Create instance variables for every key:value pair
            setattr(self, key, parameterDict[key]);


#
# method definitions
#
# writing the final netcdf output
def write_netcdf(fluxEngineObject, verbose=False):
    timeData = fluxEngineObject.time_data;
    dataLayers = fluxEngineObject.data;
    runParams = fluxEngineObject.runParams;

    outputChunk = int(runParams.run_count % runParams.output_temporal_chunking);
    #fluxEngineObject.logger.debug("Writing netCDF output with output_chunk = %d", outputChunk);
    
    if outputChunk == 0: #Create a new file and write output to it.
        nx = fluxEngineObject.nx;
        ny = fluxEngineObject.ny;
        latitudeData = fluxEngineObject.latitude_data;
        longitudeData = fluxEngineObject.longitude_data;
        latitudeGrid = fluxEngineObject.latitude_grid;
        longitudeGrid = fluxEngineObject.longitude_grid;
        
    
        function="write_netcdf";
        
        #open a new netCDF file for writing.
        #need to set format type, defaults to NetCDF4
        ncfile = Dataset(runParams.output_path, 'w', format="NETCDF3_64BIT_OFFSET");#,format='NETCDF3_CLASSIC');
    
        #Assign units attributes to coordinate var data. This attaches a
        #text attribute to each of the coordinate variables, containing the
        #units.
    
        #create the lat and lon dimensions.
        if len(latitudeData.shape)<2:#IGA - If the initial latitude data was a vector, make latitude and longitude as dimensions
            ncfile.createDimension('latitude',ny)
            ncfile.createDimension('longitude',nx)
            ncfile.createDimension('time',int(fluxEngineObject.runParams.output_temporal_chunking));
            dims = tuple(('time','latitude','longitude'))
        else:
            ncfile.createDimension('y',ny)
            ncfile.createDimension('x',nx)
            ncfile.createDimension('time',)
            dims = tuple(('time','y','x'))
    
        secs = ncfile.createVariable('time',dtype('float64').char,('time',))
        secs.units = 'seconds since 1970-01-01 00:00:00'
        secs.axis = "T"
        secs.long_name = "Time - seconds since 1970-01-01 00:00:00"
        secs.standard_name = "time"
        if (runParams.temporal_resolution != "monthly") and (runParams.temporal_resolution != None): #Fill in all the time points. Can't really do this with monthly temporal resolution.
            for t in range(0, len(secs[:])):
                secs[t] = timeData + runParams.temporal_resolution.total_seconds()*t;
        else:
            secs[outputChunk] = timeData;
        secs.valid_min = 0.0
        secs.valid_max = 1.79769313486232e+308
    
        # Define the coordinate variables. They will hold the coordinate
        # information, that is, the latitudes and longitudes.
        if len(latitudeData.shape)<2:#IGA - If the initial latitude data was a vector, write data as vectors
            lats2 = ncfile.createVariable('latitude',dtype('float64').char,('latitude'))
            lons2 = ncfile.createVariable('longitude',dtype('float64').char,('longitude'))
            # Assign units attributes to coordinate var data. This attaches a
            # text attribute to each of the coordinate variables, containing the
            # units.
            lats2.units = 'degrees_north'
            lats2.axis = "Y"
            lats2.long_name = "Latitude North"
            lats2.standard_name = "latitude"
    
            lons2.units = 'degrees_east'
            lons2.axis = "X"
            lons2.long_name = "Longitude East"
            lons2.standard_name = "longitude"
           
            lats2[:] = latitudeData
            lons2[:] = longitudeData
    
            lats2.valid_min = -90.0
            lats2.valid_max = 180.0
    
            lons2.valid_min = -180.0
            lons2.valid_max = 360.0
        else:# if the input lat/long was a grid, write output only as a grid.
            lats = ncfile.createVariable('latitude', dtype('float64').char, dims);
            lons = ncfile.createVariable('longitude', dtype('float64').char, dims);
            # Assign units attributes to coordinate var data. This attaches a
            # text attribute to each of the coordinate variables, containing the
            # units.
            lats.units = 'degrees_north'
            lats.axis = "Y"
            lats.long_name = "Latitude North"
            lats.standard_name = "latitude"
    
            lons.units = 'degrees_east'
            lons.axis = "X"
            lons.long_name = "Longitude East"
            lons.standard_name = "longitude"
           
            lats[:] = latitudeGrid
            lons[:] = longitudeGrid
    
            lats.valid_min = -90.0
            lats.valid_max = 180.0
    
            lons.valid_min = -180.0
            lons.valid_max = 360.0
    
        #
        #data layers
        #
        for dataLayerName in dataLayers:
            try:
                if verbose:
                    print("Writing datalayer '"+dataLayerName+"' to netCDF file as "+dataLayers[dataLayerName].netCDFName);
                variable = ncfile.createVariable(dataLayers[dataLayerName].netCDFName, dtype('float64').char, dims, fill_value=DataLayer.fill_value)
                data = dataLayers[dataLayerName].fdata; #fdata is usually a view by sometimes a copy so it has to be done this way. There is probably a better way to do this.
                data.shape = (dataLayers[dataLayerName].nx, dataLayers[dataLayerName].ny);
                variable[outputChunk, :, :] = data;
            except AttributeError as e:
                print("%s:No netCDFName or data attribute found in DataLayer '%s'." % (function, dataLayerName));
                raise e;
            except ValueError as e:
                print(type(e), e.args);
                print("%s: Cannot resize datalayer '%s'" % (function, dataLayers[dataLayerName].name));
                raise e;
            
            variable.missing_value = missing_value;
            variable.scale_factor = 1.0;
            variable.add_offset = 0.0;
            
            try:
                if dataLayers[dataLayerName].units != None:
                    variable.units = dataLayers[dataLayerName].units;
            except AttributeError:
                print("%s: No units found for datalayer named '%s'." % (function, dataLayerName));
            
            try:
                if dataLayers[dataLayerName].minBound != None:
                    variable.valid_min = dataLayers[dataLayerName].minBound;
            except AttributeError:
                print("%s: No minBound found for datalayer named '%s'." % (function, dataLayerName));
            
            try:
                if dataLayers[dataLayerName].maxBound != None:
                    variable.valid_max = dataLayers[dataLayerName].maxBound;
            except AttributeError:
                print("%s: No maxBound found for datalayer named '%s'." % (function, dataLayerName));
            
            try:
                if dataLayers[dataLayerName].standardName != None:
                    variable.standard_name = dataLayers[dataLayerName].standardName;
            except AttributeError:
                print("%s: No standardName found for datalayer named '%s'." % (function, dataLayerName));
            
            try:
                if dataLayers[dataLayerName].longName != None:
                    variable.long_name = dataLayers[dataLayerName].longName;
            except AttributeError:
                print("%s: No longName found for datalayer named '%s'." % (function, dataLayerName));
        
        #set some global attributes
        setattr(ncfile, 'Conventions', 'CF-1.6') 
        setattr(ncfile, 'Institution', 'Originally developed by the partners of the ESA OceanFlux GHG and OceanFlux GHG Evolution projects. Now continued by the CarbonLab team at the University of Exeter.') 
        setattr(ncfile, 'Contact', 'email: j.d.shutler@exeter.ac.uk')
        
        #Output all the parameters used in this run.
        for paramName in list(vars(runParams).keys()):
            paramValue = getattr(runParams, paramName);
            if paramValue is not None:
                if type(paramValue) is bool: #netCDF does not support bool types.
                    paramValue = int(paramValue);
                elif isinstance(paramValue, timedelta): #netCDF does not support object instances.
                    paramValue = str(paramValue);
                elif paramValue == None: #netCDF does not support None type.
                    paramValue = "None";
                
                setattr(ncfile, paramName, paramValue);
        
        if int(fluxEngineObject.runParams.output_temporal_chunking) != 1:
            setattr(ncfile, "start_year", fluxEngineObject.runParams.year);
            setattr(ncfile, "start_month", fluxEngineObject.runParams.month);
            setattr(ncfile, "start_day", fluxEngineObject.runParams.day);
            setattr(ncfile, "start_hour", fluxEngineObject.runParams.hour);
            setattr(ncfile, "start_minute", fluxEngineObject.runParams.minute);
            setattr(ncfile, "start_second", fluxEngineObject.runParams.second);
     
        ncfile.close();
    
    else: #Not the first temporal point, so open existing file and write to it
        ncfile = Dataset(fluxEngineObject.runParams.output_path, 'r+');
        
        #Write time
        ncfile.variables["time"][outputChunk] = fluxEngineObject.time_data;
                
        #Update data layers
        dataLayers = fluxEngineObject.data;
        for dataLayerName in dataLayers:
            try:
                if verbose:
                    print("Writing datalayer '"+dataLayerName+"' to netCDF file as "+dataLayers[dataLayerName].netCDFName);
                #variable = ncfile.createVariable(dataLayers[dataLayerName].netCDFName, dtype('float64').char, dims, fill_value=DataLayer.fill_value)
                data = dataLayers[dataLayerName].fdata; #fdata is usually a view by sometimes a copy so it has to be done this way. There is probably a better way to do this.
                data.shape = (dataLayers[dataLayerName].nx, dataLayers[dataLayerName].ny);
                ncfile.variables[dataLayers[dataLayerName].netCDFName][outputChunk,:,:] = data;
            except AttributeError:
                print("%s:No netCDFName or data attribute found in DataLayer '%s'." % (function, dataLayerName));
        
        #Update data range
        setattr(ncfile, "end_year", fluxEngineObject.runParams.year);
        setattr(ncfile, "end_month", fluxEngineObject.runParams.month);
        setattr(ncfile, "end_day", fluxEngineObject.runParams.day);
        setattr(ncfile, "end_hour", fluxEngineObject.runParams.hour);
        setattr(ncfile, "end_minute", fluxEngineObject.runParams.minute);
        setattr(ncfile, "end_second", fluxEngineObject.runParams.second);
        
        
        ncfile.close();
        
    return 0;


#Calculates solubility of distilled water (i.e. assuming 0 salinity) given global temperature.
#overwrites solubilityDistilled with the calculated value.
#TODO: no need to pass nx, ny into all these functions.
#TODO: rain_wet_deposition_switch isn't used!
def calculate_solubility_distilled(solubilityDistilled, salinity, rain_wet_deposition_switch, sstskin, deltaT, nx, ny, schmidtParameterisation):
    #First create a 0 salinity dataset
    salDistil = array([missing_value] * len(salinity))
    for i in arange(nx * ny):
        if (salinity[i] != missing_value):
            salDistil[i] = 0.0
        else:
            salDistil[i] = missing_value
    
    #Next calculate the solubility using the zero salinity 'distilled water' dataset
    if schmidtParameterisation == "schmidt_Wanninkhof2014":
        solubilityDistilled = solubility_Wanninkhof2014(sstskin, salDistil, deltaT, nx, ny, True);
    elif schmidtParameterisation == "schmidt_Wanninkhof1992":
        solubilityDistilled = solubility_Wanninkhof1992(sstskin, salDistil, deltaT, nx, ny, True);
    return solubilityDistilled;

#whitecapping data using relationship from TS and parameters from table 1 of Goddijn-Murphy et al., 2010, equation r1
#Takes two DataLayers as input. whitecap is written in place.
def calculate_whitecapping(windu10, whitecap):
    if (not isinstance(windu10, DataLayer)) or (not isinstance(whitecap, DataLayer)):
        raise ValueError("ofluxghg_flux_calc.calculate_whitecapping: Invalid arguments. Arguments must be DataLayer type.");
    if windu10.nx != whitecap.nx or windu10.ny != whitecap.ny:
        raise ValueError("ofluxghg_flux_calc.calculate_whitecapping: Invalid arguments. windu10 and whitecap dimensions do not match (%d, %d vs %d, %d)." % (windu10.nx, windu10.ny, whitecap.nx, whitecap.ny));
    
    for i in arange(windu10.nx * windu10.ny):
        if (windu10.fdata[i] != missing_value):
            whitecap.fdata[i] = 0.00159 * pow(windu10.fdata[i], 2.7)
        else:
            whitecap.fdata[i] = missing_value
    return whitecap.fdata;


def add_noise(data, rmse_value, mean_value, no_elements, clipNegative=False):
# randomly adding noise to data array, based log normal distribution
# rmse value is used as the standard deviation of the noise function
    #Add noise
    numClipped = 0;
    totalClipped = 0.0;
    for i in arange(no_elements):
        if ( (data[i] != missing_value) and (data[i] != 0.0) ):
            newVal = data[i] + normalvariate(0, rmse_value);
            data[i] = newVal;
            if clipNegative and newVal < 0: #Clip negative values to 0.
                data[i] = 0.0;
                numClipped += 1;
                totalClipped -= newVal;
    
    #If noise was greater than 20% of the mean value, print warning and report clipping.
    if clipNegative:
        meanVal = mean(data[where(data != missing_value)]);
        if rmse_value/abs(meanVal) > 0.2:
            print("WARNING: Adding noise to input datalayer but RMSE value is > 20% of the mean value for this input. This may result in negative values which will be clipped at 0, and therefore indirectly add bias to the input data.");
        print("INFO: A total of %d grid cells were clipped resulting in an approximate bias of %f." % (numClipped, totalClipped/no_elements));
    
    return data;


def add_noise_and_bias_wind(winduData, moment2, moment3, rmseValue, meanValue, biasValue, numElements):
    # randomly adding noise to wind data layer, based log normal distribution
    # rmse value is used as the standard deviation of the noise function
    # noise is added to wind moment2 and moment3 scaled by power 2 or power 3
    
    print("INFO: Adding noise to wind overwrites windu10_moment2 and windu10_moment3 with (windu10+epsilon)^2 and (windu10+epsilon)^3.");
    # intialise the random number generator
    numClipped = 0;
    totalClipped = 0.0;
    for i in arange(numElements):
        if ( (winduData[i] != missing_value) and (winduData[i] != 0.0) ):
            newVal = winduData[i] + normalvariate(0, rmseValue) + biasValue;
            winduData[i] = newVal;
            if newVal < 0: #Clip negative values to 0
                winduData[i] = 0.0;
                numClipped += 1;
                totalClipped -= newVal;

            moment2[i] = winduData[i]**2;
            moment3[i] = winduData[i]**3;
    
    #If noise was greater than 20% of the mean value, print warning and report clipping.
    meanVal = mean(winduData[where(winduData != missing_value)]);
    if rmseValue/meanVal > 0.2:
        print("WARNING: Adding noise to input datalayer but RMSE value is > 20% of the mean value for this input. This may result in negative values which will be clipped at 0, and therefore indirectly add bias to the input data.");
        print("WARNING: A total of %d grid cells were clipped resulting in an approximate additional bias of %f." % (numClipped, totalClipped/numElements));
    
    return (winduData, moment2, moment3);


######Ians alternative version (thanks PLand) - Untested#############

# def add_noise(data, err_value, mean_value, no_elements):
#  # randomly adding noise to data array, based log normal distribution
#  # rmse value is used as the standard deviation of the noise function

#     if type(err_value) is not float:#If not a float, its an array of spatially varying error variances (rain)
#         sizetup = data.shape
#         data = (data + np.random.normal(size = sizetup) * err_value).clip(0)
#     else:
#         # intialise the random number generator
#         for i in arange(no_elements):
#             if ( (data[i] != missing_value) and (data[i] != 0.0) ):         
#                 orig = data[i]
#                 value = log(data[i])
#                 stddev = float(err_value/orig)
#                 noise = normalvariate(0,stddev) # determines the random noise value based on the rain input data and the standard deviation of the uncertainty
#                 value = value + noise
#                 data[i] = exp(value)
#         return data
######Ians altered version############END

def add_bias_k_biology_wind(data, bias_k_value, biology_data, biology_value, wind_data, wind_value, no_elements, percent_switch):
 # adding bias offset to data array based on biology and wind conditions

   for i in arange(no_elements):
      if ( (data[i] != missing_value) and (biology_data[i] != missing_value) and (wind_data[i] != missing_value) ):         
         if ( (biology_data[i] >= biology_value) and (wind_data[i] <= wind_value) ):
            if percent_switch == 1:
               data[i] -= (data[i]*(bias_k_value/100.0))
            else:
               data[i] += bias_k_value
   return data
   
def add_bias(data, bias_value, no_elements):
 # adding bias offset to data array

   for i in arange(no_elements):
      if ( (data[i] != missing_value) ):         
         orig = data[i]
         data[i] = orig + bias_value
   return data

def add_sst_rain_bias(data, bias_value, rain_intensity, rain_data, wind_speed, wind_data, no_elements):
 # adding bias offset to data array based on rain intensity and wind speed

   for i in arange(no_elements):
      if ( (data[i] != missing_value) ):
         if ((rain_data[i] >= rain_intensity) and (wind_data[i] <= wind_speed)):        
            orig = data[i]
            data[i] = orig + bias_value
   return data


def median_filter2D(datain, nx, ny):
   # median filter, converted from C routine

    # calculating patch dimensions
   size = 5;  # patch size (square)
   w = (size-1)/2

   data_temp = zeros([nx, ny])
   vector = array([0.0] * size*size)

    # calculating patches and storing in 'vector', ready for sorting
   for x in range(w,nx-w):
      for y in range(w,ny-w):
         for xx in range(0,size):
            for yy in range(0,size):
               pic_y = y - w + yy
               pic_x = x - w + xx
               pixel = datain[pic_x, pic_y]               
               vector[xx+(yy*size)] = pixel
      
         masked_data = ma.masked_array(vector, vector == missing_value)
         if (datain[x,y] > missing_value):
            data_temp[x,y] = ma.median(masked_data)
         else:
            data_temp[x,y] = missing_value      

     # storing result back into original
   datain = data_temp

   return datain

#determine the schmidt number
#based on Schmid relationship from Wanninkhof1992 - Relationship between wind speed and gas exchange over the ocean, JGR Oceans
def schmidt_Wanninkhof1992(sstC_fdata, nx, ny, gas):
#calculating the schmidt data

   sc_fdata = array([missing_value] * nx*ny)
   if 'o2' in gas.lower():
       for i in arange(nx * ny):
          if (sstC_fdata[i] != missing_value):
             sc_fdata[i] = 1953.4 - (128.0 * sstC_fdata[i]) + (3.9918 * (sstC_fdata[i] * sstC_fdata[i])) - (0.050091 * (sstC_fdata[i] * sstC_fdata[i] * sstC_fdata[i]))#IGA-O2
          else:
        # assigning invalid values
             sc_fdata[i] = missing_value

   if 'n2o' in gas.lower():
       for i in arange(nx * ny):
          if (sstC_fdata[i] != missing_value):            
             sc_fdata[i] = 2301.1 - (151.1 * sstC_fdata[i]) + (4.7364 * (sstC_fdata[i] * sstC_fdata[i])) - (0.059431 * (sstC_fdata[i] * sstC_fdata[i] * sstC_fdata[i]))#IGA-N2O
          else:
        # assigning invalid values
             sc_fdata[i] = missing_value

   if 'ch4' in gas.lower():
       for i in arange(nx * ny):
          if (sstC_fdata[i] != missing_value):
              sc_fdata[i] = 2039.2 - (120.31 * sstC_fdata[i]) + (3.4209 * (sstC_fdata[i] * sstC_fdata[i])) - (0.040437 * (sstC_fdata[i] * sstC_fdata[i] * sstC_fdata[i]))#IGA-CH4
          else:
        # assigning invalid values
             sc_fdata[i] = missing_value
   if 'co2' in gas.lower():
       for i in arange(nx * ny):
          if (sstC_fdata[i] != missing_value):
              # relationship is only valid for temperatures <=30.0 oC
              sc_fdata[i] = 2073.1 - (125.62 * sstC_fdata[i]) + (3.6276 * (sstC_fdata[i] * sstC_fdata[i])) - (0.043219 * (sstC_fdata[i] * sstC_fdata[i] * sstC_fdata[i]))#IGA-CO2
          else:
        # assigning invalid values
             sc_fdata[i] = missing_value
   return sc_fdata

#based on Schmid relationship from Wanninkhof2014 - Relationship between wind speed and gas exchange over the ocean revisited, Limnology and Oceanography
def schmidt_Wanninkhof2014(sstC_fdata, nx, ny, gas):
    sc_fdata = array([missing_value] * nx*ny)
    if 'o2' in gas.lower():
        for i in arange(nx * ny):
            if (sstC_fdata[i] != missing_value):
                sc_fdata[i] = 1920.4 + (-135.6 * sstC_fdata[i]) + (5.2122 * sstC_fdata[i]**2) + (0.10939 * sstC_fdata[i]**3) + (0.00093777 * sstC_fdata[i]**4);
            else:
            # assigning invalid values
                sc_fdata[i] = missing_value

    if 'n2o' in gas.lower():
        for i in arange(nx * ny):
            if (sstC_fdata[i] != missing_value):            
                sc_fdata[i] = 2356.2 + (-166.38 * sstC_fdata[i]) + (6.3952 * sstC_fdata[i]**2) + (-0.13422 *sstC_fdata[i]**3) + (0.0011506 * sstC_fdata[i]**4);
            else:
            # assigning invalid values
                sc_fdata[i] = missing_value

    if 'ch4' in gas.lower():
        for i in arange(nx * ny):
            if (sstC_fdata[i] != missing_value):
                sc_fdata[i] = 2101.2 + (-131.54 * sstC_fdata[i]) + (4.4931 * sstC_fdata[i]**2) + (-0.08676 * sstC_fdata[i]**3) + (0.00070663 * sstC_fdata[i]**4);
            else:
            # assigning invalid values
                sc_fdata[i] = missing_value
    if 'co2' in gas.lower():
        for i in arange(nx * ny):
            if (sstC_fdata[i] != missing_value):
                # relationship is only valid for temperatures <=30.0 oC
                sc_fdata[i] = 2116.8 + (-136.25 * sstC_fdata[i]) + (4.7353 * sstC_fdata[i]**2) + (-0.092307 * sstC_fdata[i]**3) + (0.0007555 * sstC_fdata[i]**4);
            else:
            # assigning invalid values
                sc_fdata[i] = missing_value
    return sc_fdata
    

#solubility calculation equation from Table A2 of Wanninkkhof, JGR, 1992
def solubility_Wanninkhof1992(sstK, sal, deltaT, nx, ny, flux_calc, gas):
    #solubility calculation
    #equation from Table A2 of Wanninkkhof, JGR, 1992
    sol = array([missing_value] * nx*ny)
    if gas == 'co2':
        for i in arange(nx * ny):
            if ( (sstK[i] != missing_value) and (sal[i] != missing_value) and (sstK[i] > 0.0) ):
                sol[i] = -60.2409 + ( 93.4517*(100.0 / sstK[i]) ) + (23.3585 * (log(sstK[i]/100.0))) + (sal[i] * (0.023517 + ( (-0.023656)*(sstK[i]/100.0)) + (0.0047036*( (sstK[i]/100.0)*(sstK[i]/100.0) ) ) ) );
                sol[i] = exp(sol[i])
                #runParams.flux_calc is a switch to remove Delta_T-Sb component - ie selects use of RAPID or EQUILIBRIUM flux models from Woolf et al., 2016
                if flux_calc != 2:
                    deltaT[i] = 0.0
                sol[i] = sol[i] * (1 - (0.015*deltaT[i]))
            else:
                sol[i] = missing_value
    elif gas == 'o2':
        for i in arange(nx * ny):
            if ( (sstK[i] != missing_value) and (sal[i] != missing_value) and (sstK[i] > 0.0) ):
                sol[i] = -58.3877 + ( 85.8079*(100.0 / sstK[i]) ) + (23.8439 * (log(sstK[i]/100.0))) + (sal[i] * (-0.034892 + ( (0.015568)*(sstK[i]/100.0)) + (-0.0019387*( (sstK[i]/100.0)*(sstK[i]/100.0) ) ) ) );
                sol[i] = exp(sol[i])
                #runParams.flux_calc is a switch to remove Delta_T-Sb component - ie selects use of RAPID or EQUILIBRIUM flux models from Woolf et al., 2016
                if flux_calc != 2:
                    deltaT[i] = 0.0
                sol[i] = sol[i] * (1 - (0.015*deltaT[i]))
            else:
                sol[i] = missing_value
        
    elif gas == 'n2o':
        for i in arange(nx * ny):
            if ( (sstK[i] != missing_value) and (sal[i] != missing_value) and (sstK[i] > 0.0) ):
                sol[i] = -64.8539 + ( 100.2520*(100.0 / sstK[i]) ) + (25.2049 * (log(sstK[i]/100.0))) + (sal[i] * (-0.062544 + ( (0.035337)*(sstK[i]/100.0)) + (-0.0054699*( (sstK[i]/100.0)*(sstK[i]/100.0) ) ) ) );
                sol[i] = exp(sol[i])
                #runParams.flux_calc is a switch to remove Delta_T-Sb component - ie selects use of RAPID or EQUILIBRIUM flux models from Woolf et al., 2016
                if flux_calc != 2:
                    deltaT[i] = 0.0
                sol[i] = sol[i] * (1 - (0.015*deltaT[i]))
            else:
                sol[i] = missing_value
        
    elif gas == 'ch4':
        for i in arange(nx * ny):
            if ( (sstK[i] != missing_value) and (sal[i] != missing_value) and (sstK[i] > 0.0) ):
                sol[i] = -68.8862 + ( 101.4956*(100.0 / sstK[i]) ) + (28.7314 * (log(sstK[i]/100.0))) + (sal[i] * (-0.076146 + ( (0.043970)*(sstK[i]/100.0)) + (-0.0068672*( (sstK[i]/100.0)*(sstK[i]/100.0) ) ) ) );
                sol[i] = exp(sol[i])
                #runParams.flux_calc is a switch to remove Delta_T-Sb component - ie selects use of RAPID or EQUILIBRIUM flux models from Woolf et al., 2016
                if flux_calc != 2:
                    deltaT[i] = 0.0
                sol[i] = sol[i] * (1 - (0.015*deltaT[i]))
            else:
                sol[i] = missing_value
    
    return sol

#solubility calculation equation Table 2 of Wanninkhof, Rik. "Relationship between wind speed and gas exchange over the ocean revisited." Limnology and Oceanography: Methods 12.6 (2014): 351-362.
def solubility_Wanninkhof2014(sstK, sal, deltaT, nx, ny, flux_calc, gas):
    #solubility calculation
    #equation from Table A2 of Wanninkkhof, JGR, 1992
    sol = array([missing_value] * nx*ny)
    if gas == 'co2':
        for i in arange(nx * ny):
            if ( (sstK[i] != missing_value) and (sal[i] != missing_value) and (sstK[i] > 0.0) ):
                sol[i] = -58.0931 + ( 90.5069*(100.0 / sstK[i]) ) + (22.2940 * (log(sstK[i]/100.0))) + (sal[i] * (0.027766 + ( (-0.025888)*(sstK[i]/100.0)) + (0.0050578*( (sstK[i]/100.0)*(sstK[i]/100.0) ) ) ) );
                sol[i] = exp(sol[i])
                #runParams.flux_calc is a switch to remove Delta_T-Sb component - ie selects use of RAPID or EQUILIBRIUM flux models from Woolf et al., 2016
                if flux_calc != 2:
                    deltaT[i] = 0.0
                sol[i] = sol[i] * (1 - (0.015*deltaT[i]))
            else:
                sol[i] = missing_value
    elif gas == 'o2':
        for i in arange(nx * ny):
            if ( (sstK[i] != missing_value) and (sal[i] != missing_value) and (sstK[i] > 0.0) ):
                sol[i] = -58.3877 + ( 85.8079*(100.0 / sstK[i]) ) + (23.8439 * (log(sstK[i]/100.0))) + (sal[i] * (-0.034892 + ( (0.015568)*(sstK[i]/100.0)) + (-0.0019387*( (sstK[i]/100.0)*(sstK[i]/100.0) ) ) ) );
                sol[i] = exp(sol[i])
                #runParams.flux_calc is a switch to remove Delta_T-Sb component - ie selects use of RAPID or EQUILIBRIUM flux models from Woolf et al., 2016
                if flux_calc != 2:
                    deltaT[i] = 0.0
                sol[i] = sol[i] * (1 - (0.015*deltaT[i]))
            else:
                sol[i] = missing_value
        
    elif gas == 'n2o':
        for i in arange(nx * ny):
            if ( (sstK[i] != missing_value) and (sal[i] != missing_value) and (sstK[i] > 0.0) ):
                sol[i] = -62.7062 + ( 97.3066*(100.0 / sstK[i]) ) + (24.1406 * (log(sstK[i]/100.0))) + (sal[i] * (-0.058420 + ( (0.033193)*(sstK[i]/100.0)) + (-0.0051313*( (sstK[i]/100.0)*(sstK[i]/100.0) ) ) ) );
                sol[i] = exp(sol[i])
                #runParams.flux_calc is a switch to remove Delta_T-Sb component - ie selects use of RAPID or EQUILIBRIUM flux models from Woolf et al., 2016
                if flux_calc != 2:
                    deltaT[i] = 0.0
                sol[i] = sol[i] * (1 - (0.015*deltaT[i]))
            else:
                sol[i] = missing_value
        
    elif gas == 'ch4':
        for i in arange(nx * ny):
            if ( (sstK[i] != missing_value) and (sal[i] != missing_value) and (sstK[i] > 0.0) ):
                sol[i] = -68.8862 + ( 101.4956*(100.0 / sstK[i]) ) + (28.7314 * (log(sstK[i]/100.0))) + (sal[i] * (-0.076146 + ( (0.043970)*(sstK[i]/100.0)) + (-0.0068672*( (sstK[i]/100.0)*(sstK[i]/100.0) ) ) ) );
                sol[i] = exp(sol[i])
                #runParams.flux_calc is a switch to remove Delta_T-Sb component - ie selects use of RAPID or EQUILIBRIUM flux models from Woolf et al., 2016
                if flux_calc != 2:
                    deltaT[i] = 0.0
                sol[i] = sol[i] * (1 - (0.015*deltaT[i]))
            else:
                sol[i] = missing_value
    
    return sol


#Returns the molecular mass of the compound used in the output units. This depends on the gas for which fluxes are being calculated.
def get_output_unit_molecular_mass(gas):
    if gas.lower() == "co2":
        return 12.0108; #Output units will be in carbon.
    elif gas.lower() == "n2o":
        return 44.013; #N2O output units
    elif gas.lower() == "ch4":
        return 16.04; #CH4 output units
    else:
        raise ValueError("Unrecognised gas specification in configuration file: ", gas);


#Calculate the mass boundary layer concentration (ie concentration in the water)
#Each argument should be supplied as fdata matrix (flattened matrix)
def calculate_concw(concFactor, foundationSolubility, fCO2water, concw):
    for i in range(len(concw)):
        if ( (foundationSolubility[i] != missing_value) and (fCO2water[i] != missing_value)):
            concw[i] = foundationSolubility[i] * concFactor * fCO2water[i];
        else:
            concw[i] = missing_value;

#Calculate the interfacial concentration (i.e. at the interface between the ocean and the atmosphere)
#Each argument should be supplied as fdata matrix (flattened matrix)
def calculate_conca(concFactor, skinSolubility, fCO2air, conca):
    for i in range(len(conca)):
        if ( (skinSolubility[i] != missing_value) and (fCO2air[i] != missing_value)):
            conca[i] = skinSolubility[i] * concFactor * fCO2air[i];
        else:
            conca[i] = missing_value;


#Copies missing_values from 'master' to 'derived' e.g. for making mean and stddev datasets consistent.
#master and derived should be DataLayers
def copy_missing_values(master, derived, missingValue=DataLayer.missing_value):
    if master.nx != derived.nx and master.ny != derived.ny:
        raise ValueError("copy_missing_values: master and derived DataLayers do not have the same dimensions: (%d, %d) versus (%d, %d)." % (master.nx, master.ny, derived.nx, derived.ny));
    for i in arange(master.nx * master.ny):
        if (master.fdata[i] == missingValue):
            derived.fdata[i] = missing_value;

#Creates / adds to the failed_quality_fdata layer for each element of 'data' which exceeds the specified min/max range.
def check_output_dataset(datalayer, failed_quality_fdata):
    #function = "(check_output, main)"
    #checking the contents of a dataset
    for i in arange(len(datalayer.fdata)):
        if (datalayer.fdata[i] != DataLayer.missing_value):
            if datalayer.maxBound != None and datalayer.fdata[i] > datalayer.maxBound:
                #print "\n%s Dataset %s fails OceanFluxGHG TS table 7 valid limits, (defined min/max are: %lf/%lf, found %lf at grid point %d (%lf) exiting" % (function, name, min_range, max_range, data[i], i, i/nx)
                #print "\n%s Coincident data values sstskinC_fdata:%lf sstfndC_fdata:%lf windu10_fdata:%lf sig_wv_ht_fdata:%lf solfnd_fdata:%lf solskin_fdata:%lf k_fdata:%lf concw_fdata:%lf conca_fdata:%lf sal_fdata:%lf" % (function, sstskinC_fdata[i], sstfndC_fdata[i], windu10_fdata[i], sig_wv_ht_fdata[i], solfnd_fdata[i], solskin_fdata[i], k_fdata[i], concw_fdata[i], conca_fdata[i], sal_fdata[i])
                if failed_quality_fdata[i] == DataLayer.missing_value:
                    failed_quality_fdata[i] = 1; # first entry so need to initiase
                else:
                    failed_quality_fdata[i] += 1;
            elif datalayer.minBound != None and datalayer.fdata[i] < datalayer.minBound:
                #print "\n%s Dataset %s fails OceanFluxGHG TS table 7 valid limits, (defined min/max are: %lf/%lf, found %lf at grid point %d (%lf) exiting" % (function, name, min_range, max_range, data[i], i, i/nx)
                #print "\n%s Coincident data values sstskinC_fdata:%lf sstfndC_fdata:%lf windu10_fdata:%lf sig_wv_ht_fdata:%lf solfnd_fdata:%lf solskin_fdata:%lf k_fdata:%lf concw_fdata:%lf conca_fdata:%lf sal_fdata:%lf" % (function, sstskinC_fdata[i], sstfndC_fdata[i], windu10_fdata[i], sig_wv_ht_fdata[i], solfnd_fdata[i], solskin_fdata[i], k_fdata[i], concw_fdata[i], conca_fdata[i], sal_fdata[i])
                if failed_quality_fdata[i] == DataLayer.missing_value:
                    failed_quality_fdata[i] = 1; # first entry so need to initiase
                else:
                    failed_quality_fdata[i] += 1;


#datalayer is a DataLayer object, nx and ny are the output dimensions (must be a factor of the matrix dimensions).
def average_pixels(datalayer, nx, ny, missing_value):
    '''averages arr into superpixels each consisting of the mean of a n x n
    window in arr. Seems to be intended to go from a 0.5 x 0.5 degree grid to a
    1 x 1 degree grid, in which case n must be set to 2. Checks for and ignores
    values equal to missing_value.'''
    function = "(average_pixels, main)"
    
    #copy fdata, because data may be out of data if not a view.
    data = datalayer.fdata;
    data.shape = datalayer.data.shape; #Reshape to resemble fdata.
    
    #determin n
    if len(data.shape) == 2:
        datany = data.shape[0];
        datanx = data.shape[1];
    elif len(data.shape) == 3: #Sometimes there is a time dimension with only one entry in
        datany = data.shape[1];
        datanx = data.shape[2];
    else:
        raise ValueError("%s: Unexpected number of dimensions (%d) when trying to resize data layer." % (function, len(data.shape)));
    
    if datany%ny != 0 and datanx%nx != 0:
        print(datanx, datany)
        raise ValueError("%s: Cannot rescale data layer because global data layer dimensions are not a whole multiple of the data layer's dimensions." % function);
    n = (datanx/nx);
    n2 = (datany/ny);
    if n!=n2:
        raise ValueError("%s: Scaling data layers by irregular scaling factors is not supported (e.g. %d != %d)." % (function, n, n2));
    
    
    #print "%s Averaging sstgrad_fdata into 1x1 degree grid (N=%d)" % (function, n)
    nj0, ni0 = data.shape
    nj, ni = [nj0 / n, ni0 / n]
    if nj * n != nj0 or ni * n != ni0:
       print("Dimensions ", nj0, ni0, " indivisible by ", n);
    out = resize(missing_value, [nj, ni])
    for j in range(nj):
       j0 = j * n
       j1 = j0 + n
       for i in range(ni):
          i0 = i * n
          a = data[j0:j1, i0:i0 + n]
          w = nonzero(a != missing_value)
          if size(w) > 0:
             out[j, i] = mean(a[w])
    
    datalayer.data = data;
    datalayer.calculate_fdata(); #Resampled data, so now we must recalculate fdata.
    return out

#Returns false if dataLayer dimensions do not match the reference dimensions.
def check_dimensions(dataLayer, ref_nx, ref_ny, DEBUG=False):
   function = "(check_dimensions, main)"
   if dataLayer.data.shape[1] == ref_nx and dataLayer.data.shape[0] == ref_ny:
      if DEBUG:
         print("\n%s Input data (%s) have identical dimensions to reference values (%s, %s) "% (function, dataLayer.name, dataLayer.nx, dataLayer.ny))
         return True;
   else:
      print("\n%s Input data ('%s') dimensions are non-identical to reference (new: %s, %s is not equal to: %s, %s)." % (function, dataLayer.name, dataLayer.nx, dataLayer.ny, ref_nx, ref_ny))
      return False;


class FluxEngine:
    def __init__(self, parameterDict):
        self.runParams = RunParameters();
        self.runParams.set_parameters(parameterDict);
        
        #self.rootDirectory = path.abspath(path.dirname(inspect.stack()[0][1]));
        #TODO: Shouldn't used src_home, this is going to be removed from config file soon. Use rootDirectory instead.
        self.defaultSettings = Settings(path.join(path.dirname(__file__), "settings.xml")); #Load default settings metadata for datalayers.
        self.data = {};
        self.kParameterisationFunctors = [] #List of functor objects which encapsulates the k rate calculation.
                                            #The k rate calculation can be extended in a modular fashion.
                                            #Each functor is called in order.
        self.processIndicatorFunctors = []; #List of functor objects which calculate process indicator layers.
                                            #Each functor is called in order, hence functors later in the list can use outputs from those earlier as inputs.
        
        self._load_lon_lat_time();
    
    #Adds a data layer by reading a netCDF file. Optionally adds stddev and count data from the same file.
    #preprocessing is an optional list of functions to transform the data before storing.
    #TODO: transposeData should now be added as a preprocessing function
    def add_data_layer(self, name, infile, prod, stddevProd=None, countProd=None, transposeData=False, preprocessing=None):
        function = "(fe_core, FluxEngine.add_data_layer)";
        
        try:
            self._add_single_data_layer(name, infile, prod, transposeData, preprocessing=preprocessing);
            
        except KeyError as e:
            print("\n%s: Data variable '%s' is missing from %s input (%s)" % (function, prod, name, infile));
            print(e, e.args);
            raise e;
           # return False;
#        except ValueError as e: #E.g. incorrect number of dimensions
#            print "\n%s: %s" % (function, e.args);
#            raise e;
#            #return False;
        
        #TODO: stddev and count should be turned into separate datalayers in the config file processing stage, rather
        #      than testing for _prods and adding them here. Then _add_single_data_layer can be merged with this function.
        #If they have been specified, search for stddev.
        if stddevProd != None:
            if stddevProd == "auto": #automatically detect stddev product
                stddevProd = prod.rsplit("_", 1)[0]+"_stddev";
            try:
                self._add_single_data_layer(name+"_stddev", infile, stddevProd, transposeData=transposeData);
                copy_missing_values(self.data[name], self.data[name+"_stddev"]); #Any missing values in 'name' should be copied to the stddev dataset for consistency.
            #except TypeError:
            except KeyError as e:
                print("here", type(e), e.args);
                print("%s: Did not find stddev variable ('%s') in netCDF file for %s." % (function, stddevProd, name));
            except ValueError as e: #E.g. incorrect number of dimensions
                print("\n%s: %s" % (function, e.args));
                return False;
        
        if countProd != None:
            if countProd == "auto": #automatically detect count product
                countProd = prod.rsplit("_", 1)[0]+"_count";
            try:
                self._add_single_data_layer(name+"_count", infile, countProd, transposeData=transposeData);
                copy_missing_values(self.data[name], self.data[name+"_count"]); #Any missing values in 'name' should be copied to the count dataset for consistency.
            except KeyError:
                print("%s: Did not find stddev variable ('%s') in netCDF file for %s." % (function, countProd, name));
            except ValueError as e: #E.g. incorrect number of dimensions
                print("\n%s: %s" % (function, e.args));
                return False;
            
        #Datalayer was successfully added.
        return True;
    
    #Adds a single datalayer and acquires metadata
    #TODO: transposeData is now be handled by preprocessing, so this cna be removed...
    def _add_single_data_layer(self, name, infile, prod, transposeData=False, preprocessing=None):
        #function = "(ofluxghg_flux_calc, FluxEngine._add_single_data_layer)";
        
        metaData = self._extract_data_layer_meta_data(name);
        
        inputChunk = int(self.runParams.run_count % metaData.temporalChunking);
        timeIndex = inputChunk * (metaData.temporalSkipInterval+1);
        dl = DataLayer.create_from_file(name, infile, prod, metaData, timeIndex, transposeData=transposeData, preprocessing=preprocessing);
        self.data[name] = dl;
    
    #Creates a DataLayer which is filled (by default) with DataLayer.missing_value.
    #The DataLayer will be added to self.data and is accessable by it's 'name', e.g. self.data["new_datalayer"].
    #If no dimensions (nx, ny) are specified, then the default self.nx and self.ny dimensions are used.
    #If there is known metadata for the datalayer (e.g. min, max, units etc.) these will be loaded from settings.xml
    #Metadata can be overwritten in the config file using the the datalayer name and the xml attribute from the settings.xml file, e.g.:
    #       datalayername_units = C m^2s^-1
    #       datalayername_maxBound = 100.0
    def add_empty_data_layer(self, name, nx=None, ny=None, fillValue=DataLayer.missing_value):
        if nx==None: nx=self.nx;
        if ny==None: ny=self.ny;
        
        metaData = self._extract_data_layer_meta_data(name);        
        dl = DataLayer.create_empty_datalayer(name, nx, ny, metaData, fillValue=fillValue);
        
        self.data[name] = dl;
    
    #Returns a DataLayerMetaData instance for the named DataLayer.
    #First it will search for default values specified in the settings xml file
    #Failing that the default default values are used
    #Finally any metadata values which are specified in the config file are used to overwrite default values
    def _extract_data_layer_meta_data(self, name):
        #Extract default metadata
        if name in self.defaultSettings.allDataLayers: #If default values are specified in the settings.xml file
            metaData = self.defaultSettings.allDataLayers[name];
        else: #Otherwise create metadata with 'default' default values
            metaData = DataLayerMetaData(name);
        
        #overwrite default metadata with any found in the run parameters that correspond to the specified name.
        for attribute in vars(metaData):
            if attribute != name:
                if name+"_"+attribute in vars(self.runParams):
                    if attribute in ["temporalChunking", "temporalSkipInterval"]: #Integers
                        setattr(metaData, attribute, int(getattr(self.runParams, name+"_"+attribute))); #TODO: should be handled in setup_tools::create_run_parameters really.
                    else: #floats
                        setattr(metaData, attribute, getattr(self.runParams, name+"_"+attribute));
        
        return metaData;
    
    #Defines the k parameterisation to be used. Requires an callable object which:
    #   implements input_names() and output_names()
    #   and contains standard_name and long_name string attributes.
    def add_k_parameterisation_component(self, kObject):
        function = "(FluxEngine.add_k_parameterisation_component)";
        if kObject != None:
            self.kParameterisationFunctors.append(kObject);
        else:
            raise ValueError("%s: Trying to add 'None' k parameterisation component." % function);
    
    def add_process_indicator_functor(self, piObject):
        function = "(FluxEngine.add_process_indicator_functor)";
        if piObject != None:
            self.processIndicatorFunctors.append(piObject);
        else:
            raise ValueError("%s: Trying to add 'None' process indicator layer component." % function);

    def run(self):
        status = self._check_datalayers(); #Check for consistency of datalayers (e.g. dimensions all match one another).
        if status == True:
            return self._run_fluxengine(self.runParams);
        else:
            print("Not all datalayers are consistent (check dimensions of data).");
            return status;
    
    #Reads longitude, latitude, time and dimension sizes (nx, ny).
    def _load_lon_lat_time(self):
        function = "(ofluxghg_flux_calc, FluxEngine._load_lon_lat_time)";
        
        #Find the path of the relevent datalayer.
        try:
            axesDatalayerInfile = getattr(self.runParams, self.runParams.axes_data_layer+"_infile");
        except ValueError as e:
            print("Couldn't find file path for axes_data_layer (%s). Check that this is correctly defined in the config file (it must match the name of another datalayer)." % self.runParams.axes_data_layer);
            print(e.args);
        
        #Read longitude, latitude and time data from the sstskin infile, so open this file.
        try:
            dataset = Dataset(axesDatalayerInfile);
        except IOError as e:
            print("\n%s: axes_data_layer (%s) inputfile %s does not exist" % (function, self.runParams.axes_data_layer, axesDatalayerInfile));
            print(type(e), e.args);
        
        #Read lat and lat
        try:
            self.latitude_data = dataset.variables[self.runParams.latitude_prod][:];
            self.longitude_data = dataset.variables[self.runParams.longitude_prod][:];
            if self.latitude_data[0]<0: #IGA - it is a vector that is in opposite orientation to 'taka'
                self.latitude_data = flipud(self.latitude_data);
        except KeyError:
            raise ValueError ("%s: Couldn't find longitude (%s) and/or latitude (%s) variables in %s. Have you set longitude_prod and latitude_prod correctly in your configuration file?" % (function, self.runParams.longitude_prod, self.runParams.latitude_prod, axesDatalayerInfile));

        #Determine if already a grid, if not calculate lon and lat grids.
        if len(self.latitude_data.shape) == 1: #not already a grid
            self.latitude_grid = meshgrid(ones((len(self.longitude_data))), self.latitude_data)[1]#IGA - gridding vectored geo data so that they can be treated the sam way as non-regular grids
            self.longitude_grid = meshgrid(self.longitude_data, ones((len(self.latitude_data))))[1]#IGA
        else: #IGA - geo data are already a grid
            self.latitude_grid = self.latitude_data;
            self.longitude_grid = self.longitude_data;
            
        #set time (since 1st Jan 1970)
        try:
            curDatetime = datetime(self.runParams.year, self.runParams.month, self.runParams.day, self.runParams.hour, self.runParams.minute, self.runParams.second);
            self.time_data = (curDatetime - datetime(1970, 1, 1)).total_seconds();
            #self.time_data = dataset.variables[self.runParams.time_prod][:];
        except KeyError:
            raise ValueError("%s: Couldn't find time (%s%) variables in %s. Have you set time_prod correctly in your configuration file?" % (function, self.runParams.time_prod, self.runParams.sstskin_infile));

        #set dimensions
        self.ny, self.nx = self.latitude_grid.shape;
        print("%s: Grid dimensions set to: (%d, %d)" % (function, self.ny, self.nx));

    #Ran after each datalayer is read in. Checks for consistency between datalayers.
    #Rescales data layers which can be rescaled.
    #returns true if all data layers are successfully validated.
    def _check_datalayers(self):
        #Check that dimensions match
        for key in self.data:
            if check_dimensions(self.data[key], self.nx, self.ny, DEBUG) == False:
                #Dimensions don't match, so try to rescale it.
                try:
                    print("Attempting to rescale datalayer '%s'."%key);
                    self.data[key].data = average_pixels(self.data[key], self.nx, self.ny, DataLayer.missing_value);
                    self.data[key].ny, self.data[key].nx = self.data[key].data.shape;
                    self.data[key].calculate_fdata();
                    print("Successfully rescaled datalayer '%s' to"%key, self.data[key].data.shape);
                except ValueError as e:
                    print(e.args);
                    return False;
        return True;
    
    #Helper function which returns True if a data layer was read from file as input
    #   as opposed to being generated by FluxEngine. Otherwise returns False.
    def _input_data_provided(self, dataLayerName):
        return (dataLayerName+"_path" in dir(self.runParams)) and (dataLayerName+"_prod" in dir(self.runParams));
    
    
    #Applies the mask to all data layers
    def _apply_mask(self):
        if "mask" in self.data:
            toIgnore = where(self.data["mask"].fdata == 0);
            
            #Apply mask to each data layer
            for dataLayerName in self.data:
                if dataLayerName != "mask":
                    self.data[dataLayerName].fdata[toIgnore] = self.data[dataLayerName].missing_value = True;


    def _run_fluxengine(self, runParams):
        function = "(ofluxghg_flux_calc, FluxEngine._run_fluxengine_)";
        
        #Set up logging object
        try:
            self.logger = logging.getLogger('FluxEngine_debug_log');
            #hdlr = logging.FileHandler(os.path.join(workingDirectory, runParams.LOG_PATH), filemode='w')
            hdlr = logging.FileHandler(runParams.LOG_PATH);
            formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s');
            hdlr.setFormatter(formatter);
            self.logger.addHandler(hdlr);
            if DEBUG_LOGGING == True:
                self.logger.setLevel(logging.DEBUG);
            else:
                self.logger.setLevel(logging.DEBUG);
        except:
            print("\n%s Couldn't initialise logger at path %s" % (function, runParams.LOG_PATH));
        
        #TODO: Replace directly with self.nx, self.ny, no need to use local variables here.
        nx = self.nx;
        ny = self.ny;
        
        #Apply mask to all data layers, if applicable.
        self._apply_mask(); #Checks for existance of mask internally.

        ### Adding empty data layers for data that are computed later.
        #If there isn't any sstfnd data create empty arrays to fill later. (i.e. sstfnd = sstskin + runParams.cool_skin_difference);
        #Temporarily add stddev and count data for any sst which is missing these values.
        if ("sstfnd" not in self.data) and (self._input_data_provided("sstfnd")==False):
            self.add_empty_data_layer("sstfnd");
            self.add_empty_data_layer("sstfnd_stddev");
            self.add_empty_data_layer("sstfnd_count");
        if "sstskin_stddev" not in self.data and "sstskin_count" not in self.data: #Remove this, stddev and count aren't used except when they are!
            self.add_empty_data_layer("sstskin_stddev");
            self.add_empty_data_layer("sstskin_count");
        if "sstfnd_stddev" not in self.data and "sstfnd_count" not in self.data: #Remove this, stddev and count aren't used except when they are!
            self.add_empty_data_layer("sstfnd_stddev");
            self.add_empty_data_layer("sstfnd_count");

        
        if runParams.TAKAHASHI_DRIVER == True:
           #need to generate the moment2 and moment3 data
           self.add_empty_data_layer("windu10_moment2");
           self.add_empty_data_layer("windu10_moment3");
           for i in arange(self.nx * self.ny):
              if self.data["windu10"].fdata[i] != missing_value:
                 self.data["windu10_moment2"].fdata[i] = self.data["windu10"].fdata[i]*self.data["windu10"].fdata[i]
                 self.data["windu10_moment3"].fdata[i] = self.data["windu10"].fdata[i]*self.data["windu10"].fdata[i]*self.data["windu10"].fdata[i]
        
        #some specific pco2 conditions
        #TODO: This shouldn't be randomly here.
        if "pgas_sw" in self.data:
            self.data["pgas_sw"].fdata[abs(self.data["pgas_sw"].fdata)<0.1]=DataLayer.missing_value#IGA_SOCATv4

        if (runParams.pco2_data_selection != 1):
            #signifies that we're NOT using SOCAT data, which means there is no stddev data for pgas_sw
            self.add_empty_data_layer("pgas_sw_stddev");
            
        #initialising some data structures for the calculations
        self.add_empty_data_layer("scskin");
        self.add_empty_data_layer("scfnd");
        self.add_empty_data_layer("solubility_skin");
        self.add_empty_data_layer("solubility_fnd");
        self.add_empty_data_layer("pH2O");
        
        #apply additive saline skin value, but why selectively apply it?
        #TODO: Why are these hard-coded values. Should be using minBound and maxBound?
        self.add_empty_data_layer("salinity_skin");
        for i in arange(self.nx * self.ny):
            if (self.data["salinity"].fdata[i] >= 0.0) and (self.data["salinity"].fdata[i] <= 50.0):
                if ( (self.data["salinity"].fdata[i] + runParams.saline_skin_value) <= 50.0):
                    self.data["salinity_skin"].fdata[i] = self.data["salinity"].fdata[i] + runParams.saline_skin_value
            else:
                self.data["salinity"].fdata[i] = missing_value
                self.data["salinity_skin"].fdata[i] = missing_value
        if (runParams.saline_skin_value != 0.0):
            print("%s Using the saline skin model (%lf psu added to skin salinities)" % (function, runParams.saline_skin_value))
        
        #conversion of rain data from mm day-1 to mm hr^-1
        if "rain" in self.data:
            for i in arange(self.nx * self.ny):   
               if (self.data["rain"].fdata[i] != DataLayer.missing_value):
                   self.data["rain"].fdata[i] /= 24.0;
        
        
		#ability to randomly perturb the input datasets
        #needed for the ensemble analyses
        #stddev of noise is using published RMSE for each dataset
        #all datasets are considered to have bias=0, hence mean of noise=0
        if (runParams.random_noise_windu10_switch == 1):
           #add_noise(self.data["windu10"].fdata, 0.44, 0.0, nx*ny)
           #For the 0.8 value, see: https://podaac.jpl.nasa.gov/Cross-Calibrated_Multi-Platform_OceanSurfaceWindVectorAnalyses
           add_noise_and_bias_wind(self.data["windu10"].fdata, self.data["windu10_moment2"].fdata, self.data["windu10_moment3"].fdata, runParams.windu10_noise, 0.0, runParams.windu10_bias, nx*ny);
           print("%s Adding random noise to windu10_mean, windu10_moment2 and windu10_moment3 (mean 0.0, stddev 0.44 ms^-1 - assuming using ESA GlobWave data)" % (function))
        
        if (runParams.random_noise_sstskin_switch == 1):
           add_noise(self.data["sstskin"].fdata, runParams.sstskin_noise, 0.0, nx*ny)
           print("%s Adding random noise to sstskin (mean 0.0, stddev 0.14 ^oC - assuming using ESA CCI ARC data)" % (function))
        
        if (runParams.random_noise_sstfnd_switch == 1):
           add_noise(self.data["sstfnd"].fdata, runParams.sstfnd_noise, 0.0, nx*ny)
           print("%s Adding random noise to sstfnd (mean 0.0, stddev 0.6 ^oC - assuming using OSTIA data)" % (function))
        
        if (runParams.random_noise_pco2_switch == 1):
           print("/n%s Shape of pco2 data", self.data["pgas_sw"].fdata.shape)
           add_noise(self.data["pgas_sw"].fdata, runParams.pco2_noise, 0.0, nx*ny, clipNegative=True) #SOCAT uncertainty
           print("%s Adding random noise to pco2/fco2 data (mean 0.0, stddev 6 uatm - Using Candyfloss data, value provided by J Shutler)" % (function))
        
        
        #interpreting fnd_data option
        ####Derive sstskin as necessary.
        #Two possible datasets: sstskin and sstfnd, from different parts of the water column.
        #One option: sst gradients
        #If using only one dataset then copy that over the other. Otherwise keep both.
        #runParams.cool_skin_difference is the assumed temperature difference between the skin and foundation layer
        
        #using sstfnd, so copy sstfnd into sstskin
        #sstskin = sstfnd
        if runParams.sst_gradients_switch == 0 and (self._input_data_provided("sstskin")==False) and (self._input_data_provided("sstfnd")==True):
           print("%s SST gradient handling is off, using SSTfnd data selection in configuration file for all components of the flux calculation (this will ignore any SSTskin data in configuration file)." % (function))
           #copy sstfnd data into the sstskin dataset to make sure
           if "sstskin" not in self.data: #Must add the sstskin layer first!
                self.add_empty_data_layer("sstskin");
                if "sstfnd_stddev" in self.data:
                    self.add_empty_data_layer("sstskin_stddev");
                if "sstfnd_count" in self.data:    
                    self.add_empty_data_layer("sstskin_count");
           for i in arange(nx * ny):
               if self.data["sstfnd"].fdata[i] != missing_value:
                   self.data["sstskin"].fdata[i] = self.data["sstfnd"].fdata[i]
                   #if "sstskin_stddev" in self.data and "sstskin_count" in self.data:
                   self.data["sstskin_stddev"].fdata[i] = self.data["sstfnd_stddev"].fdata[i]
                   self.data["sstskin_count"].fdata[i] = self.data["sstfnd_count"].fdata[i]
               else:
                   self.data["sstskin"].fdata[i] = missing_value;
        
        
        #IGA added for the case where only foundation is provided and gradients are on------------------------------
        #must estimate sstskin (= sstfnd - runParams.cool_skin_difference)
        elif runParams.sst_gradients_switch == 1 and (self._input_data_provided("sstskin")==False) and (self._input_data_provided("sstfnd")==True):
            print("%s Using SSTfnd data selection with correction for skin temperature (SSTskin = SSTfnd - %f)(ignoring SSTskin data in configuration file)." % (function, runParams.cool_skin_difference));
            #actually copy sstfnd data into the sstskin dataset to make sure
            if "sstskin" not in self.data: #Must add the sstskin layer first!
                self.add_empty_data_layer("sstskin");
                if "sstfnd_stddev" in self.data:
                    self.add_empty_data_layer("sstskin_stddev");
                if "sstfnd_count" in self.data:    
                    self.add_empty_data_layer("sstskin_count");
            for i in arange(nx * ny): #sstdkin = sstfnd - runParams.cool_skin_difference
                if self.data["sstfnd"].fdata[i] != missing_value:
                    self.data["sstskin"].fdata[i] = self.data["sstfnd"].fdata[i]-runParams.cool_skin_difference
                    self.data["sstskin_stddev"].fdata[i] =  self.data["sstfnd_stddev"].fdata[i]
                    self.data["sstskin_count"].fdata[i] = self.data["sstfnd_count"].fdata[i]
                else:
                    self.data["sstskin"].fdata[i] = missing_value;
        
        #Using sstskin, so calculate it from sstfnd.
        elif runParams.sst_gradients_switch == 0 and (self._input_data_provided("sstskin")==True) and (self._input_data_provided("sstfnd")==False):
           print("%s SST gradient handling is off, using SSTskin to derive SSTfnd (SSTfnd = SSTskin + %f) for flux calculation (ignoring SSTfnd data in configuration file)." % (function, runParams.cool_skin_difference));
           #In this case SST gradients are not used and only sstskin is provided. sstfnd is needed for the flux calculation, so:
           #    Calculate sstfnd from sstskin
           #    Copy sstfnd over sstskin to prevent accidental use of sst gradients
           #    TODO: Infer sst_gradients from flux equation, infer use_sstskin and use_sstfnd from data input sources.
           
           #setting sstfnd_ data fields to skin values
           for i in arange(nx * ny):
              if self.data["sstskin"].fdata[i] != missing_value:
                 
                 self.data["sstfnd"].fdata[i] = self.data["sstskin"].fdata[i] + runParams.cool_skin_difference #calculate sstfnd from sstskin for use in bulk equation
                 self.data["sstskin"].fdata[i] = self.data["sstfnd"].fdata[i] #Copy sstfnd over sstskin to prevent accidental use of sst_gradients
                 self.data["sstfnd_stddev"].fdata[i] = self.data["sstskin_stddev"].fdata[i]
                 self.data["sstfnd_count"].fdata[i] = self.data["sstskin_count"].fdata[i]
              else:
                 self.data["sstfnd"].fdata[i] = missing_value
                 self.data["sstskin"].fdata[i] = missing_value
                 
        elif runParams.sst_gradients_switch == 1 and (self._input_data_provided("sstskin")==True) and (self._input_data_provided("sstfnd")==False):
           print("%s SST gradient handling is on, using SSTskin and SSTfnd = SSTskin + %f for flux calculation (ignoring SSTfnd data in configuration file)." % (function, runParams.cool_skin_difference));
            #setting sstfnd_ data fields to skin values  
           for i in arange(nx * ny):
              if self.data["sstskin"].fdata[i] != missing_value:
                 self.data["sstfnd"].fdata[i] = self.data["sstskin"].fdata[i] + runParams.cool_skin_difference
                 self.data["sstfnd_stddev"].fdata[i] =  self.data["sstskin_stddev"].fdata[i]
                 self.data["sstfnd_count"].fdata[i] = self.data["sstskin_count"].fdata[i]
              else:
                 self.data["sstfnd"].fdata[i] = missing_value
        elif runParams.sst_gradients_switch == 0 and (self._input_data_provided("sstskin")==True) and (self._input_data_provided("sstfnd")==True):
           print("%s SST gradient handling is off, using SSTfnd and SSTskin from the configuration file." % (function))
        elif runParams.sst_gradients_switch == 1 and (self._input_data_provided("sstskin")==True) and (self._input_data_provided("sstfnd")==True):
           print("%s SST gradient handling is on, using SSTfnd and SSTskin from the configuration file." % (function))
        else:
           print("\n%s sst_gradients_switch (%s), sstskin provided (%s), sstfnd provided (%s) combination not recognised, exiting." % (function, runParams.sst_gradients_switch==1, self._input_data_provided("sstskin"), self._input_data_provided("sstfnd")))
           return 1;
       
        
        #If there is no pco2_sst data, we need to get it / generate it.
        if "pco2_sst" not in self.data: #SOCATv4
            try:
                print("No pco2_sst data was supplied. sstfnd will be used instead.")
                self.add_empty_data_layer("pco2_sst");
                self.data["pco2_sst"].fdata = self.data["sstfnd"].fdata-273.15; #copy/convert sstfnd
            except (IOError, KeyError, ValueError) as e:
                print("pco2_sst data not available and could read sstfnd so cannot proceed.");
                print(type(e), "\n"+e.args);
                return 1;
       

        #quality filtering and conversion of SST datasets
        #Calculate sstskinC and sstfndC
        self.add_empty_data_layer("sstskinC");
        self.add_empty_data_layer("sstfndC");
        for i in arange(nx * ny):
            if self.data["sstskin"].fdata[i] != missing_value:
                self.data["sstskinC"].fdata[i] = self.data["sstskin"].fdata[i] - 273.15;
                self.data["sstfndC"].fdata[i] = self.data["sstfnd"].fdata[i] - 273.15;
        
        
        
        #Note: wind bias added with noise now before calculating new first and second moments
#         # bias values to be added here
#        if (runParams.bias_windu10_switch == 1):
#           add_bias(self.data["windu10"].fdata, runParams.bias_windu10_value, nx*ny)
#            # makes no sense to add bias to second and third order moments, as any bias in the system 
#            # would only impact on the mean (which is the 1st order moment)
#           print("%s Adding bias to windu10_mean (not to second and third order moments) (value %lf ms^-1)" % (function, runParams.bias_windu10_value))
        
        if (runParams.bias_sstskin_switch == 1):
           add_bias(self.data["sstskin"].fdata, runParams.bias_sstskin_value, nx*ny)
           print("%s Adding bias noise to sstskin (value %lf ^oC)" % (function, runParams.bias_sstskin_value))
        
        if (runParams.bias_sstfnd_switch == 1):
           add_bias(self.data["sstfnd"].fdata, runParams.bias_sstfnd_value, nx*ny)
           print("%s Adding bias noise to sstfnd (value %lf ^oC)" % (function, runParams.bias_sstfnd_value))
        
        if (runParams.bias_pco2_switch == 1):
           add_bias(self.data["pgas_sw"].fdata, runParams.bias_pco2_value, nx*ny)
           print("%s Adding bias noise to pco2w (value %lf uatm)" % (function, runParams.bias_pco2_value))
        
         # bias based on change in sstskin due to rain
        if (runParams.bias_sstskin_due_rain_switch == 1):
           add_sst_rain_bias(self.data["sstskin"].fdata, runParams.bias_sstskin_due_rain_value, runParams.bias_sstskin_due_rain_intensity, self.data["rain"].fdata, runParams.bias_sstskin_due_rain_wind, self.data["windu10"].fdata, nx*ny)
        
         # quality filtering of wind and Hs data
        for i in arange(nx * ny):
           if (self.data["windu10"].fdata[i] != missing_value):
               # valid range taken from GlobWave Product User Guide Phase 3 (PUG3) doucment
              if (self.data["windu10"].fdata[i] > 50.0) or (self.data["windu10"].fdata[i] < 0.0):
                 self.data["windu10"].fdata[i] = missing_value
           if "sig_wv_ht" in self.data and self.data["sig_wv_ht"].fdata[i] != missing_value:
               # valid range taken from GlobWave Product User Guide Phase 3 (PUG3) doucment
              if (self.data["sig_wv_ht"].fdata[i] > 20.0) or (self.data["sig_wv_ht"].fdata[i] < 0.0):
                 self.data["sig_wv_ht"].fdata[i] = missing_value
           #if (self.data["sigma0"].fdata[i] != missing_value):
           #   if (self.data["sigma0"].fdata[i] > ??.0):
            #     self.data["sigma0"].fdata[i] = missing_value
        
        if (runParams.pco2_data_selection != 0 and VERIFICATION_RUNS != True):
           # signifies that we're using SOCAT data or in-situ data
           #print "%s Using the SOCAT data " % (function)
           #if self.data["pco2_sst"].fdata[i] != missing_value:
           for i in arange(nx * ny):
              if isnan(self.data["pco2_sst"].fdata[i]) != True:# and self.data["pco2_sst"].fdata[i] > 0.0 ): #SOCATv4_IGA
                 if self.data["pco2_sst"].fdata[i] > 260:
                   self.data["pco2_sst"].fdata[i] = self.data["pco2_sst"].fdata[i] - 273.15#IGA - If statement added because in-situ SST data may not be in K!
              
                 if self.data["pco2_sst"].fdata[i] > 30.5 or self.data["pco2_sst"].fdata[i] < -1.8:
                    self.data["pco2_sst"].fdata[i] = missing_value
              else:
                 self.data["pco2_sst"].fdata[i] = missing_value
        
        # quality control/contrain all SST data
        # check all SST data are within -1.8 - 30.5^oC (or 271.35 - 303.65K)
        for i in arange(nx * ny):
           if ((self.data["sstskin"].fdata[i] != missing_value) and (self.data["sstskinC"].fdata[i] != missing_value) and (self.data["sstfnd"].fdata[i] != missing_value) and (self.data["sstfndC"].fdata[i] != missing_value)):
              if ( (self.data["sstskinC"].fdata[i] > 30.5) or (self.data["sstfndC"].fdata[i] > 30.5) or (self.data["sstskinC"].fdata[i] < -1.8) or (self.data["sstfndC"].fdata[i] < -1.8)):
                 self.data["sstfnd"].fdata[i] = missing_value
                 self.data["sstskin"].fdata[i] = missing_value
                 self.data["sstfndC"].fdata[i] = missing_value
                 self.data["sstskinC"].fdata[i] = missing_value
                 
        # ensure that the different SST data all cover the same spatial regions
        # convert all missing values into standard value, rather than variations that seem to exist in some of these data
        #note this is an intersect operation (compared to the above)
        for i in arange(nx * ny):
           if ((self.data["sstfnd"].fdata[i] == self.data["sstfnd"].fillValue) or (self.data["sstskin"].fdata[i] == self.data["sstskin"].fillValue)):
              self.data["sstfnd"].fdata[i] = missing_value
              self.data["sstskin"].fdata[i] = missing_value
        
        # ---------------------------------------------------------------------------------------------------------------------------------------
        # -                                                                                                                                     -
        # converting pressure data from Pascals to millibar #IGA Need to formalise pressure data units - this converts from Pa.-------------START
        if runParams.TAKAHASHI_DRIVER != True:
            if npany(self.data["pressure"].fdata > 10000.0):
                for i in arange(nx * ny):
                    if (self.data["pressure"].fdata[i] != missing_value):
                        self.data["pressure"].fdata[i] = self.data["pressure"].fdata[i] * 0.01
                print("Converted pressure data to mbar from Pa")
        # ------------------------------------------------------IGA_SOCATv4 - Was removed as using pressure data in mbar.-------------END
        # -                                                                                                                                     -
        # ---------------------------------------------------------------------------------------------------------------------------------------
        
        #########################################################
        # calculations for components of the flux calculation
        #########################################################
         
        # pCO2/fCO2 extrapolation (if turned on) from reference year
        pco2_increment = (runParams.year - runParams.pco2_reference_year) * runParams.pco2_annual_correction;
        pco2_increment_air = pco2_increment;
            
        DeltaT_fdata = array([missing_value] * nx*ny)
        
        if runParams.flux_calc == 1:
           print("%s Using the RAPID model (from Woolf et al., 2016)" % (function))
        elif runParams.flux_calc == 2:
           print("%s Using the EQUILIBRIUM model (from Woolf et al., 2016)" % (function))
           for i in arange(nx * ny):
              if ( (self.data["sstfndC"].fdata[i] != missing_value) & (self.data["sstskinC"].fdata[i] != missing_value) ):
                 DeltaT_fdata[i] = self.data["sstfndC"].fdata[i] - self.data["sstskinC"].fdata[i]
              else:
                 DeltaT_fdata[i] = missing_value
        elif runParams.flux_calc == 3:
            print("%s Using the BULK model (ie Flux = k*S*delta_pCO2)" % (function))
            print("%s Note: this assumes that the SSTskin dataset is the only temperature dataset, and so this is what will be used to calculate k and the solubility" % (function))
        else:
           print("\n%s runParams.flux_calc from configuration not recognised, exiting." % (function))
           return 1;
        
        #Calculating the schmidt number at the skin and fnd
        if runParams.schmidt_parameterisation == "schmidt_Wanninkhof2014":
            self.data["scskin"].fdata = schmidt_Wanninkhof2014(self.data["sstskinC"].fdata, nx, ny, runParams.GAS)
            self.data["scfnd"].fdata = schmidt_Wanninkhof2014(self.data["sstfndC"].fdata, nx, ny, runParams.GAS)
        elif runParams.schmidt_parameterisation == "schmidt_Wanninkhof1992":
            self.data["scskin"].fdata = schmidt_Wanninkhof1992(self.data["sstskinC"].fdata, nx, ny, runParams.GAS)
            self.data["scfnd"].fdata = schmidt_Wanninkhof1992(self.data["sstfndC"].fdata, nx, ny, runParams.GAS)
        else:
            raise ValueError("Unrecognised schmidt/solubility parameterisation selected: "+runParams.schmidtParameterisation);
        
        #Calculate solubility
        if runParams.schmidt_parameterisation == "schmidt_Wanninkhof2014":
            #calculating the skin solubility, using skin sst and salinity
            self.data["solubility_skin"].fdata = solubility_Wanninkhof2014(self.data["sstskin"].fdata, self.data["salinity_skin"].fdata, DeltaT_fdata, nx, ny, True, runParams.GAS.lower());
            #calculating the interfacial solubility
            self.data["solubility_fnd"].fdata = solubility_Wanninkhof2014(self.data["sstfnd"].fdata, self.data["salinity"].fdata, DeltaT_fdata, nx, ny, runParams.flux_calc, runParams.GAS.lower());
        elif runParams.schmidt_parameterisation == "schmidt_Wanninkhof1992":
            #calculating the skin solubility, using skin sst and salinity
            self.data["solubility_skin"].fdata = solubility_Wanninkhof1992(self.data["sstskin"].fdata, self.data["salinity_skin"].fdata, DeltaT_fdata, nx, ny, True, runParams.GAS.lower());
            #calculating the interfacial solubility
            self.data["solubility_fnd"].fdata = solubility_Wanninkhof1992(self.data["sstfnd"].fdata, self.data["salinity"].fdata, DeltaT_fdata, nx, ny, runParams.flux_calc, runParams.GAS.lower());
        else:
            raise ValueError("Unrecognised schmidt/solubility parameterisation selected: "+runParams.schmidtParameterisation);
    
         # calculate pCO2 data using mean sea level pressure data
         # equation 26, Kettle et al, 2009, ACP
         # pCO2_air = X[CO2] ( P(t) - pH2O(t) )
         # where X[CO2] is from Takahashi climatology
         # pH2O(t) = 1013.25 exp (24.4543 - 67.4509 (100/Tk(t) - 4.8489 ln (Tk(t) / 100) - 0.000544 S)
         # S = salinity, Tk = temperature in Kelvin
         # pCO2_water = pCO2_water_tak exp (0.0423 (T_foundation - T_tak)
         # T_tak = Takahashi temperature
        
        sys.stdout.flush();

        
        ###############################
        # Calculate:                  #
        #   pH2O                      #
        ###############################
        for i in arange(nx * ny):
            #if ( (self.data["salinity_skin"].fdata[i] != missing_value) and (self.data["sstskin"].fdata[i] != missing_value) and (self.data["pressure"].fdata[i] != missing_value) and (self.data["vgas_air"].fdata[i] != missing_value) and (self.data["sstfnd"].fdata[i] != missing_value) and (self.data["pco2_sst"].fdata[i] != missing_value) and (self.data["pgas_sw"].fdata[i] != missing_value) and (self.data["sstskin"].fdata[i] !=0.0) ):
            if (self.data["salinity_skin"].fdata[i] != missing_value) and (self.data["sstskin"].fdata[i] != missing_value):
                #Equation A1 in McGillis, Wade R., and Rik Wanninkhof. "Aqueous CO2 gradients for air-sea flux estimates." Marine Chemistry 98.1 (2006): 100-108.
                self.data["pH2O"].fdata[i] = 1013.25 * exp(24.4543 - (67.4509 * (100.0/self.data["sstskin"].fdata[i])) - (4.8489 * log(self.data["sstskin"].fdata[i]/100.0)) - 0.000544 * self.data["salinity_skin"].fdata[i])
            else:
                self.data["pH2O"].fdata[i] = missing_value
                    
        #######################################################
        # Calculate corrected values for pco2_sw and pco2_air #
        # Only do this is pco2 data is provided               #
        #######################################################  
        # this may be needed when using SMOS salinity data
        # To-DO: awaiting info from Lonneke and David before implementing fully
        pCO2_salinity_term = 0;
        #dSalinity = salinity_rmse
        #if (salinity_option == 1):
        #  # need to determine dS/S using original salinity (its been modified above to add the salinity_rmse)
        #  # so (self.data["salinity"].fdata[i] - salinity_rmse) ensures that we are dealing with the original value of salinity
          #  # using Ys=1 as a global correction following Sarmiento and Gruber, 2006)
        #pCO2_salinity_term = 1.0*(dSalinity/(self.data["salinity"].fdata[i] - salinity_rmse) ) #TH: Commented this out, but should not be removed as may be used in the future.
        #else:
        # pCO2_salinity_term = 0.0
        
        #if pCO2 data in sea water is provided calculated corrected values.
        if "pgas_sw" in self.data: #Only calculate partial pressure data is available
            self.add_empty_data_layer("pgas_sw_cor");
            for i in range(len(self.data["pgas_sw"].fdata)):
                #If statement below to maintain a consistent calculation with previous versions. Perhaps not needed but would invalidate reference data otherwise.
                if ( (self.data["salinity_skin"].fdata[i] != missing_value) and (self.data["sstskin"].fdata[i] != missing_value) and (self.data["pressure"].fdata[i] != missing_value) and (self.data["sstfnd"].fdata[i] != missing_value) and (self.data["pco2_sst"].fdata[i] != missing_value) and (self.data["pgas_sw"].fdata[i] != missing_value) and (self.data["sstskin"].fdata[i] !=0.0) ):
                    if runParams.GAS == 'CO2' and runParams.pco2_data_selection != 3:
                        # correction to different years, correction is data and year specific.
                        # note for 2010, correction for SOCAT isn't strictly required. However the contents of the exponential will collapse
                        # to 1 (with some rounding error expected), so effectively no correction will be applied
                        self.data["pgas_sw_cor"].fdata[i] = pco2_increment + (self.data["pgas_sw"].fdata[i] *exp( (0.0423*(self.data["sstfndC"].fdata[i] - self.data["pco2_sst"].fdata[i])) - (0.0000435*((self.data["sstfndC"].fdata[i]*self.data["sstfndC"].fdata[i]) - (self.data["pco2_sst"].fdata[i]*self.data["pco2_sst"].fdata[i]) )) + pCO2_salinity_term) );
                    else:
                        self.data["pgas_sw_cor"].fdata[i] = self.data["pgas_sw"].fdata[i];
        
        
        #if interface/air CO2 data is provided as molar fraction and there is no partial pressure then calculate partial pressure.
        if ("vgas_air" in self.data) and ("pgas_air" not in self.data):
            # vco2 in ppm * 1000000 = atm
            # result /1000 to then convert from atm to uatm
            # hence * 0.001 factor
            self.add_empty_data_layer("pgas_air");
            for i in range(len(self.data["vgas_air"].fdata)):
                #If statement below to maintain a consistent calculation with previous versions. Perhaps not needed but would invalidate reference data otherwise.
                if ( (self.data["salinity_skin"].fdata[i] != missing_value) and (self.data["sstskin"].fdata[i] != missing_value) and (self.data["pressure"].fdata[i] != missing_value) and (self.data["vgas_air"].fdata[i] != missing_value) and (self.data["sstfnd"].fdata[i] != missing_value) and (self.data["pco2_sst"].fdata[i] != missing_value) and (self.data["pgas_sw"].fdata[i] != missing_value) and (self.data["sstskin"].fdata[i] !=0.0) ):
                    if runParams.GAS == 'CO2' and runParams.ATMGAS == 'V':
                        self.data["pgas_air"].fdata[i] = (self.data["vgas_air"].fdata[i] * (self.data["pressure"].fdata[i] - self.data["pH2O"].fdata[i])) / 1013.25;
                    else:
                        self.data["pgas_air"].fdata[i] = self.data["vgas_air"].fdata[i]
        
        #Now calculate corrected values for pCO2 at the interface/air
        ###Converts from ppm to microatm TH
        self.add_empty_data_layer("pgas_air_cor");
        #If statement added below to maintain a consistent calculation with previous versions. Perhaps not needed but would invalidate reference data otherwise.
        if runParams.TAKAHASHI_DRIVER == False: #Different for takahashi run to maintain compatability with verification run. This will be updated when verification runs are updated
            for i in range(len(self.data["pgas_air"].fdata)):
                if ( (self.data["salinity_skin"].fdata[i] != missing_value) and (self.data["sstskin"].fdata[i] != missing_value) and (self.data["pressure"].fdata[i] != missing_value) and (self.data["sstfnd"].fdata[i] != missing_value) and (self.data["pco2_sst"].fdata[i] != missing_value) and (self.data["pgas_sw"].fdata[i] != missing_value) and (self.data["sstskin"].fdata[i] !=0.0) ):
                    if runParams.GAS == 'CO2' and runParams.ATMGAS == 'V':
                        ##THtodo: 1e-6 can be removed...
                        self.data["pgas_air_cor"].fdata[i] = self.data["pgas_air"].fdata[i] + pco2_increment_air;
                    else:
                        self.data["pgas_air_cor"].fdata[i] = self.data["pgas_air"].fdata[i]
        else: #runParams.TAKAHASHI_DRIVER==True #Added to maintain compatability with takahashi verification. Should be removed when the verification runs are updated.
            for i in range(len(self.data["vgas_air"].fdata)):
                #If statement below to maintain a consistent calculation with previous versions. Perhaps not needed but would invalidate reference data otherwise.
                if ( (self.data["salinity_skin"].fdata[i] != missing_value) and (self.data["sstskin"].fdata[i] != missing_value) and (self.data["pressure"].fdata[i] != missing_value) and (self.data["vgas_air"].fdata[i] != missing_value) and (self.data["sstfnd"].fdata[i] != missing_value) and (self.data["pco2_sst"].fdata[i] != missing_value) and (self.data["pgas_sw"].fdata[i] != missing_value) and (self.data["sstskin"].fdata[i] !=0.0) ):
                    if runParams.GAS == 'CO2' and runParams.ATMGAS == 'V':
                        #THtodo: 1e-6 can be removed...
                        self.data["pgas_air_cor"].fdata[i] = (self.data["vgas_air"].fdata[i] * 1e-6 * (self.data["pressure"].fdata[i] - self.data["pH2O"].fdata[i]) / (1e-6 * 1013.25)) + (pco2_increment_air)
                    else:
                        self.data["pgas_air_cor"].fdata[i] = self.data["vgas_air"].fdata[i]

                    
            
            #SOCAT, so: conversion of pCO2 to fCO2 from McGillis and Wanninkhof 2006, Marine chemistry with correction from Weiss 1974 (as the equation in 2006 paper has a set of brackets missing)
            #runParams.pco2_data_selection ==2 signifies SOCAT fCO2 data, so converting pCO2_air_cor_fdata to fCO2_air_cor_fdata      
            if runParams.pco2_data_selection == 2 or runParams.pco2_data_selection == 4 or runParams.pco2_data_selection == 45:
                b11_fdata = array([missing_value] * nx*ny);
                d12_fdata = array([missing_value] * nx*ny);
                for i in range(len(self.data["pgas_air_cor"].fdata)):
                    #If statement below to maintain a consistent calculation with previous versions. Perhaps not needed but would invalidate reference data otherwise.
                    if ( (self.data["salinity_skin"].fdata[i] != missing_value) and (self.data["sstskin"].fdata[i] != missing_value) and (self.data["pressure"].fdata[i] != missing_value) and (self.data["vgas_air"].fdata[i] != missing_value) and (self.data["sstfnd"].fdata[i] != missing_value) and (self.data["pco2_sst"].fdata[i] != missing_value) and (self.data["pgas_sw"].fdata[i] != missing_value) and (self.data["sstskin"].fdata[i] !=0.0) ):
                        # conversion of pCO2 to fCO2 from McGillis and Wanninkhof 2006, Marine chemistry with correction from Weiss 1974 (as the equation in 2006 paper has a set of brackets missing)
                        b11_fdata[i] = -1636.75 + (12.0408*self.data["sstskin"].fdata[i]) - (0.0327957*self.data["sstskin"].fdata[i]*self.data["sstskin"].fdata[i]) + (3.16528e-5 * self.data["sstskin"].fdata[i]*self.data["sstskin"].fdata[i]*self.data["sstskin"].fdata[i])
                        d12_fdata[i] = 57.7 - (0.118*self.data["sstskin"].fdata[i])
                         # gas constant
                        R = 82.0578 # in [cm^3 atm/(mol K)]
                        # 1/0.987 = 1.0131712 - conversion between bar and atm, so *1013.25 is the conversion from millibar to atm.
                        # the combination of the B11 and d12 terms are in cm^3/mol and so these cancel with the P/RT term (in mol/cm^3) so the whole of the exp term is dimensionless
                        self.data["pgas_air_cor"].fdata[i] = self.data["pgas_air_cor"].fdata[i] * exp((b11_fdata[i] + (2*d12_fdata[i]) ) * 1e-6 * ((self.data["pressure"].fdata[i] * 1013.25)/(R*self.data["sstskin"].fdata[i]) ))

        
        ######################################
        # Takahashi verification information #  ##TODO: CHECK IF THIS CAN BE REMOVED?
        ######################################
        if runParams.TAKAHASHI_DRIVER: #Assumes CO2 data input is not suppled as concentrations
            # debuggin differences in pH20 values
            pCO2a_diff_fdata = array([missing_value] * nx*ny)
            dpCO2_diff_fdata = array([missing_value] * nx*ny)
            for i in arange(nx * ny):
                #Additional pCO2 outputs for Takahashi verification
                if self.data["pgas_air"].fdata[i] != missing_value:
                    pCO2a_diff_fdata[i] = self.data["pgas_air_cor"].fdata[i] - self.data["pgas_air"].fdata[i]
                    dpCO2_diff_fdata[i] = (self.data["pgas_sw_cor"].fdata[i] - self.data["pgas_air_cor"].fdata[i]) - (self.data["pgas_sw_cor"].fdata[i] - self.data["pgas_air"].fdata[i])
            
            pH2O_takahashi_fdata = array([missing_value] * nx*ny)
            humidity_fdata = array([missing_value] * nx*ny)
            pH2O_diff_fdata = array([missing_value] * nx*ny)
            for i in arange(nx * ny):
                #Additional humidity outputs for Takahashi verification
                if self.data["pgas_air"].fdata[i] != missing_value and self.data["pH2O"].fdata[i] != missing_value and self.data["pressure"].fdata[i] != missing_value and self.data["vgas_air"].fdata[i] != missing_value:
                    pH2O_takahashi_fdata[i] = self.data["pressure"].fdata[i] -  (self.data["pgas_air"].fdata[i] *1e-6 * 1013.25) / (self.data["vgas_air"].fdata[i] * 1e-6)
                    humidity_fdata[i] = pH2O_takahashi_fdata[i]/self.data["pH2O"].fdata[i];     
                    pH2O_diff_fdata[i] = ((humidity_fdata[i])-1.0) * 100.0;    
                else:
                    pH2O_takahashi_fdata[i] = missing_value;
                    humidity_fdata[i] = missing_value;
                    pH2O_diff_fdata[i] = missing_value;



        #######################################
        # Calculating gas transfer velocity k #
        #######################################
        for kParameterisationFunctor in self.kParameterisationFunctors:
            #Check each input exists #TODO: This should go in the pre-run checks!
            for inputDataName in kParameterisationFunctor.input_names():
                if inputDataName not in self.data: #This additional check isn't really needed as it is done in the functor and in the driver script.
                    raise KeyError("Selected kParameterisation ("+kParameterisationFunctor.name+") requires input data layers which have not been provided. Required the following DataLayers:\n"+str(kParameterisationFunctor.input_names()));
            
            #Before running it is necessary to create any non-existing output layers that are required by the k calculation
            for outputDataName in kParameterisationFunctor.output_names():
                if outputDataName not in self.data:
                    self.add_empty_data_layer(outputDataName);
            
            #Call the functor to calculate k
            kParameterisationOutput = kParameterisationFunctor(self.data);
            if kParameterisationOutput == False:
                raise RuntimeError("%s: k parameterisation component (%s) returned False indicating k has not been calculated successfully." % (function, kParameterisationFunctor.name));
        

         # ability to investigate bias on k due to surface biology/slicks
         # assumes that bias values are realistic and that they won't cause the k_fdata to become unrealistic
         # also assumes that self.data["biology"].fdata exist (ie it won't if they indicator layers are off)
        if (runParams.bias_k_switch == 1):
           add_bias_k_biology_wind(self.data["k"].fdata, runParams.bias_k_value, self.data["biology"].fdata, runParams.bias_k_biology_value, self.data["windu10"].fdata, runParams.bias_k_wind_value, nx*ny, runParams.bias_k_percent_switch)
           if runParams.bias_k_percent_switch==0:
              print("\n%s Adding bias to chosen k (k_fdata) parameterisation data (bias value of %lf ms^-1 added, where biology (biology fdata) is > %lf mg m^-3 and wind speed (windu10) is < %lf m s^-1)" % (function, runParams.bias_k_value, runParams.bias_k_biology_value, runParams.bias_k_wind_value))
           else:
              print("\n%s Adding bias to chosen k (k_fdata) parameterisation data (bias value of - %lf percent being used, where biology (biology fdata) is > %lf mg m^-3 and wind speed (windu10) is < %lf m s^-1)" % (function, runParams.bias_k_value, runParams.bias_k_biology_value, runParams.bias_k_wind_value))


         ############################
         # actual flux calculation  #
         ############################
         # determine flux based on kHO6_fdata derivation
         # CO2 flux = k * s * Delta_PCO2
         # s = salinity, Delta_PCO2 calculated above
          # solubility conversion between mol kg^-1 atm^-1 to g-C m^-3 uatm^-1
          # mol -> grams for CO2 x 12
          # kg=liter -> m^-3 = x 1000
          # atm -> uatm = /1000000
          # result = x12.0108/1000 #12.0108 is atomic mass of carbon
          # provides solubility in g-C m^-3 uatm^-1
          
          # multiplying solubility (g-C m^-3 uatm^-1) by pCO2 (uatm) = concentration in g-C m^-3
          
          # kg to m^-3 = 1 liter = 1/1000 m^-3
          # 1 liter = 1 kg
          
          # k conversion from m/s to to m/day
          # k are already in cm/h)
          # x 24 (converts from hours to days)
          # x 10^-2 (converts from cm to metres)
          # 24 * 10**{-2} = *24/100
        
          # this wil give a daily flux in g-C m^-2 day^-1
          # expanding equations
        
        #k_factor = 36.0 * 24.0 / 100.0 # used if using 10^-4 ms units for k, whereas we are now using cm/h for k
        k_factor = 24.0 / 100.0
        
        #Somewhere to store the net flux
        self.add_empty_data_layer("FH06");
        #Update FH06 (OF/gas flux data layer) description to reflect the gas transfer velocity calculation used in the calculation.
        self.data["FH06"].longName = self.data["FH06"].longName % self.data["k"].name;

        #If using rain wet deposition, calculate the gas solubility in distilled water
        if runParams.rain_wet_deposition_switch:
            self.add_empty_data_layer("FKo07");
            self.add_empty_data_layer("solubility_distilled");
            calculate_solubility_distilled(self.data["solubility_distilled"].fdata, self.data["salinity"].fdata,
                                           runParams.rain_wet_deposition_switch, self.data["sstskin"], DeltaT_fdata, self.nx, self.ny, runParams.schmidt_parameterisation);
        
        if ((runParams.kb_asymmetry != 1.0) and (runParams.k_parameterisation == 3)):
           print("%s kb asymetry has been enabled (runParams.kb_asymmetry:%lf and runParams.k_parameterisation:%d)" % (function, runParams.kb_asymmetry, runParams.k_parameterisation))
           if runParams.flux_calc == 3:
               raise ValueError("kb_asymmetry is not supported for the 'bulk' calculation. Try using 'rapid' or 'equilibrium' instead.");

        ###################################################
        # If concentration data are not provided as input #
        #    calculate them from corrected pco2 data      #
        ###################################################
        outputUnitMolecularMass = get_output_unit_molecular_mass(runParams.GAS); #Different gases mean the output units will be different. Thus different molecular masses are needed to correctly calculate moles.
        concFactor = (outputUnitMolecularMass/1000.0);
        if "concw" not in self.data:
            self.add_empty_data_layer("concw");
            if runParams.flux_calc == 3: #Bulk calculation, so should use the same solubility as conca
                calculate_concw(concFactor, self.data["solubility_skin"].fdata, self.data["pgas_sw_cor"].fdata, self.data["concw"].fdata); #calculate concw
            else: #Not bulk, so use the solubility at foundation layer
                calculate_concw(concFactor, self.data["solubility_fnd"].fdata, self.data["pgas_sw_cor"].fdata, self.data["concw"].fdata); #calculate concw
        if "conca" not in self.data:
            self.add_empty_data_layer("conca");
            calculate_conca(concFactor, self.data["solubility_skin"].fdata, self.data["pgas_air_cor"].fdata, self.data["conca"].fdata); #calculate conca
        
        
        ##############################
        # Main flux calculation loop # #assume corrected pco2 data at the moment #################################
        ##############################
        for i in arange(nx * ny):
            if ( (self.data["k"].fdata[i] != missing_value) and (self.data["concw"].fdata[i] != missing_value) and (self.data["conca"].fdata[i] != missing_value) ):
                #flux calculation
                #TODO: This is partially K-parameterisation dependent. Need to decouple this...
                
                #Rapid and equilibrium flux calculations
                if runParams.flux_calc == 1 or runParams.flux_calc == 2:
                    if ((runParams.kb_asymmetry != 1.0) and (runParams.k_parameterisation == 3)):
                        kd_component = self.data["kd"].fdata[i] * k_factor
                        kb_component = self.data["kb"].fdata[i] * k_factor
                        self.data["FH06"].fdata[i] = (kd_component * (self.data["concw"].fdata[i] - self.data["conca"].fdata[i])) + (kb_component * (self.data["concw"].fdata[i] - (runParams.kb_asymmetry *self.data["conca"].fdata[i]) ) )
                    else:
                        self.data["FH06"].fdata[i] = (self.data["k"].fdata[i] * k_factor) * (self.data["concw"].fdata[i] - self.data["conca"].fdata[i])

                # using simplified flux calculation with no separation of temperature between airside and waterside CO2
                # assumes that the skin temperature dataset is the only temperature dataset
                #Bulk model (F = k*solubility*(pCO2_water - pCO2_air)
                elif runParams.flux_calc == 3:
                    #Note that that concw is calculated differently if flux_calc ==3 vs !=3, uses skin solubility (same as conca) if flux_calc==3.
                    #The below calculation is therefore not identical to RAPID, and is the correct BULK forumulation.
                    self.data["FH06"].fdata[i] = self.data["k"].fdata[i] * k_factor * (self.data["concw"].fdata[i] - self.data["conca"].fdata[i]);
                
                #Some unexpected flux calculation was specified. Should never get this far as fe_setup.py will check this.
                else:
                    raise ValueError("Unrecognised flux calculation. Currently supported options are 'rapid' (1), 'equilibrium' (2) and 'bulk' (3). Recieved: "+str(runParams.flux_calc));
                
                
                # calculating and adding in flux component for wet deposition due to rain
                if runParams.rain_wet_deposition_switch:
                    if "pgas_air_cor" not in self.data: #Check that pCO2 / vCO2 data were supplied
                        raise ValueError("Cannot use wet deposition (rain_wet_deposition_switch) without specifying pCO2 or vCO2 data.");
                    else:
                        # relationship from Komori et al.m 2007
                        # need solubility of CO2 in fresh water
                        # solubility calculated using sstkin
                        # 24.0/1000.0 = conversion from mm hr^-1 to m day^-1
                        # flux is then in g C m^-2 day^-1
                        # flux is always negative, ie going into the ocean
                        self.data["FKo07"].fdata[i] = -(self.data["rain"].fdata[i] * (24.0/1000.0)) * (concFactor * self.data["solubility_distilled"].fdata[i]) * self.data["pgas_air_cor"].fdata[i]
                        if self.data["FH06"].fdata[i] != missing_value:
                            self.data["FH06"].fdata[i] += self.data["FKo07"].fdata[i]
                        else:
                            self.data["FH06"].fdata[i] = self.data["FKo07"].fdata[i]
            else:
                self.data["FH06"].fdata[i] = missing_value
        
        
        #Adding verification data outputs at same units as T09
        if runParams.TAKAHASHI_DRIVER == True:
          solskin_takadata = array([missing_value] * nx * ny)
          FH06_takadata = array([missing_value] * nx * ny)
          for i in arange(nx * ny):
           if self.data["solubility_skin"].fdata[i] != missing_value:
              solskin_takadata[i] = self.data["solubility_skin"].fdata[i]*1000 #from 'mol kg-1 atm-1' to 'mmol kg-1 atm-1'
           if self.data["FH06"].fdata[i] != missing_value:
              FH06_takadata[i] = self.data["FH06"].fdata[i]/30.5 #from flux per day to flux per month
        else:
          solskin_takadata = array([missing_value] * nx * ny)
          FH06_takadata = array([missing_value] * nx * ny)
        
        
        #
        # quality control of datasets, following TS (OceanFluxGHG_TS_D2-9_v1.8-signed.pdf) table 7
        #
        #calculate takahashi style DpCO2 and check range
        if "pgas_sw_cor" in self.data and "pgas_air_cor" in self.data:
            self.add_empty_data_layer("dpco2_cor");
            for i in arange(self.nx * self.ny):
               if ( (self.data["pgas_sw_cor"].fdata[i] != missing_value) and (self.data["pgas_air_cor"].fdata[i] != missing_value) ):
                  self.data["dpco2_cor"].fdata[i] = (self.data["pgas_sw_cor"].fdata[i] - self.data["pgas_air_cor"].fdata[i]) 
               else:
                  self.data["dpco2_cor"].fdata[i] = missing_value

        self.add_empty_data_layer("dpconc_cor");        
        for i in arange(self.nx * self.ny):
           if ( (self.data["concw"].fdata[i] != missing_value) and (self.data["conca"].fdata[i] != missing_value) ):
              self.data["dpconc_cor"].fdata[i] = (self.data["concw"].fdata[i] - self.data["conca"].fdata[i]) 
           else:
              self.data["dpconc_cor"].fdata[i] = missing_value

        #Calculate the total number of quality violations per grid location
        self.add_empty_data_layer("failed_quality")
        #self.d = {};
        for outputVar in ["FH06", "kt", "k", "kd", "kb", "salinity", "sstskinC", "concw", "conca", "dpco2_cor", "sstfndC", "pgas_sw_cor"]:
            if outputVar in self.data: #Only check outputs that we've actually created...
                check_output_dataset(self.data[outputVar], self.data["failed_quality"].fdata);
                #self.d[outputVar] = self.data["failed_quality"].fdata.copy();
                #self.d[outputVar].shape = self.data["failed_quality"].data.shape;

        #whitecapping data using relationship from TS and parameters from table 1 of Goddijn-Murphy et al., 2010, equation r1
        self.add_empty_data_layer("whitecap");
        self.data["whitecap"].fdata = calculate_whitecapping(self.data["windu10"], self.data["whitecap"]);
        
        
        #runParams.pco2_data_selection == 1 signifies that we're using SOCAT data 
        #this dataset is filled with missing_values prior to output, otherwise non-socat data will appear in the SFUG field in the netcdf
        #which will confuse the user
        #enabled for TAKAHASHI_DRIVER to enable checking
        if runParams.TAKAHASHI_DRIVER != True:
           if (self.runParams.pco2_data_selection == 0):
              self.data["pgas_sw"].fdata = array([missing_value] * self.nx*self.ny)
              self.data["pgas_sw_stddev"].fdata = array([missing_value] * self.nx*self.ny)
        
        #
        #procesing indictor attribute layers
        #
        #temp if clause using the -l flag for now #TODO: remove this, config should just specify required process indicator layers
        if (self._input_data_provided("sstfnd")==True):
            for piFunctor in self.processIndicatorFunctors:
                try:
                    for inputDataLayer in piFunctor.input_names(): #Check all required inputs exist.
                        if inputDataLayer not in self.data:
                            raise ValueError("%s: Missing input DataLayer (%s) for process indicator functor %s." % (function, inputDataLayer, piFunctor.name));
    
                    for outputDataLayer in piFunctor.output_names(): #Add any output DataLayers that don't already exist
                        if outputDataLayer not in self.data:
                            self.add_empty_data_layer(outputDataLayer);
                    #Execute the process indicator layer functor.
                    piFunctor(self.data);
                except ValueError as e:
                    print(e.args);
                    print("Exiting...");
                    return 1;
        
        #Substitute gas name into meta data / human-readable descriptions
        if "vgas_air" in self.data:
            self.data["vgas_air"].standardName = Template(self.data["vgas_air"].standardName ).safe_substitute(GAS=self.runParams.GAS);
            self.data["vgas_air"].longName = Template(self.data["vgas_air"].longName).safe_substitute(GAS=self.runParams.GAS);
        if "pgas_air" in self.data:
            self.data["pgas_air"].standardName = Template(self.data["pgas_air"].standardName ).safe_substitute(GAS=self.runParams.GAS);
            self.data["pgas_air"].longName = Template(self.data["pgas_air"].longName).safe_substitute(GAS=self.runParams.GAS);
        if "pgas_sw" in self.data:
            self.data["pgas_sw"].standardName = Template(self.data["pgas_sw"].standardName ).safe_substitute(GAS=self.runParams.GAS);
            self.data["pgas_sw"].longName = Template(self.data["pgas_sw"].longName).safe_substitute(GAS=self.runParams.GAS);
        if "pgas_sw_stddev" in self.data:
            self.data["pgas_sw_stddev"].standardName = Template(self.data["pgas_sw_stddev"].standardName ).safe_substitute(GAS=self.runParams.GAS);
            self.data["pgas_sw_stddev"].longName = Template(self.data["pgas_sw_stddev"].longName).safe_substitute(GAS=self.runParams.GAS);
        if "pgas_air_cor" in self.data:
            self.data["pgas_air_cor"].standardName = Template(self.data["pgas_air_cor"].standardName ).safe_substitute(GAS=self.runParams.GAS);
            self.data["pgas_air_cor"].longName = Template(self.data["pgas_air_cor"].longName).safe_substitute(GAS=self.runParams.GAS);
        if "pgas_sw_cor" in self.data:
            self.data["pgas_sw_cor"].standardName = Template(self.data["pgas_sw_cor"].standardName ).safe_substitute(GAS=self.runParams.GAS);
            self.data["pgas_sw_cor"].longName = Template(self.data["pgas_sw_cor"].longName).safe_substitute(GAS=self.runParams.GAS);
        if "FH06" in self.data:
            self.data["FH06"].standardName = Template(self.data["FH06"].standardName ).safe_substitute(GAS=self.runParams.GAS);
            self.data["FH06"].longName = Template(self.data["FH06"].longName).safe_substitute(GAS=self.runParams.GAS);
        if "conca" in self.data:
            self.data["conca"].standardName = Template(self.data["conca"].standardName ).safe_substitute(GAS=self.runParams.GAS);
            self.data["conca"].longName = Template(self.data["conca"].longName).safe_substitute(GAS=self.runParams.GAS);
        if "concw" in self.data:
            self.data["concw"].standardName = Template(self.data["concw"].standardName ).safe_substitute(GAS=self.runParams.GAS);
            self.data["concw"].longName = Template(self.data["concw"].longName).safe_substitute(GAS=self.runParams.GAS);
        if "kt" in self.data:
            self.data["kt"].longName = Template(self.data["kt"].standardName).safe_substitute(GAS=self.runParams.GAS);
        if "kd" in self.data:
            self.data["kd"].longName = Template(self.data["kt"].standardName).safe_substitute(GAS=self.runParams.GAS);
        if "kt" in self.data:
            self.data["kb"].longName = Template(self.data["kb"].standardName).safe_substitute(GAS=self.runParams.GAS);
        
        ### write out the final ouput to netcdf
        write_netcdf(self);        
        print("%s SUCCESS writing file %s" % (function, runParams.output_path))
#        
#        #Finally, close the logger
        handlers = self.logger.handlers[:];
        for handler in handlers:
            handler.close();
        
        sys.stdout.flush();
        return 0;

