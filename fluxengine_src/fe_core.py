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

# ToDo (wishList)
# Add uncertainty values as an option through configuration file
# Incorporate rain noise through configuration file
# Pressure units often vary - account for this through config?
# Line 1736 (re.match) returns an error if runParams.pco2_prod is 'pco2'
# 

 # netcdf bits
from netCDF4 import Dataset
import sys
from math import log, exp, pow, isnan;
from numpy import size, flipud, mean, zeros, nonzero, array, resize, ma, arange, dtype, ones, meshgrid;
from numpy import any as npany;
from random import normalvariate
import logging;
from os import path;

from datalayer import DataLayer, DataLayerMetaData;
from settings import Settings;
from debug_tools import calc_mean; #calculate mean ignoring missing values.

 # debug mode switches
DEBUG = False
DEBUG_PRODUCTS = True
TAKAHASHI_DRIVER = False # enables Takahashi data to drive code - used for verifying calculations
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
        for key in vars(self).keys(): #remove all previously stored values
            delattr(self, key);
        for key in parameterDict.keys(): #Create instance variables for every key:value pair
            setattr(self, key, parameterDict[key]);


#
# method definitions
#
# writing the final netcdf output
def write_netcdf(fluxEngineObject, verbose=False):
    latitudeData = fluxEngineObject.latitude_data;
    longitudeData = fluxEngineObject.longitude_data;
    timeData = fluxEngineObject.time_data;
    nx = fluxEngineObject.nx;
    ny = fluxEngineObject.ny;
    dataLayers = fluxEngineObject.data;
    runParams = fluxEngineObject.runParams;
    latitudeGrid = fluxEngineObject.latitude_grid;
    longitudeGrid = fluxEngineObject.longitude_grid;
    

    function="write_netcdf";
    
    #open a new netCDF file for writing.
    #need to set format type, defaults to NetCDF4
    ncfile = Dataset(runParams.output_path,'w',format='NETCDF3_CLASSIC');

    #Assign units attributes to coordinate var data. This attaches a
    #text attribute to each of the coordinate variables, containing the
    #units.

    #create the lat and lon dimensions.
    if len(latitudeData.shape)<2:#IGA - If the initial latitude data was a vector, make latitude and longitude as dimensions
        ncfile.createDimension('latitude',ny)
        ncfile.createDimension('longitude',nx)
        ncfile.createDimension('time',1)
        dims = tuple(('time','latitude','longitude'))
    else:
        ncfile.createDimension('y',ny)
        ncfile.createDimension('x',nx)
        ncfile.createDimension('time',1)
        dims = tuple(('time','y','x'))

    secs = ncfile.createVariable('time',dtype('float64').char,('time',))
    secs.units = 'seconds since 1970-01-01 00:00:00'
    secs.axis = "T"
    secs.long_name = "Time - seconds since 1970-01-01 00:00:00"
    secs.standard_name = "time"
    secs[:] = timeData
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
        lats2.valid_max = 90.0

        lons2.valid_min = -180.0
        lons2.valid_max = 180.0
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
        lats.valid_max = 90.0

        lons.valid_min = -180.0
        lons.valid_max = 180.0

    #
    #data layers
    #
    for dataLayerName in dataLayers:
        try:
            if verbose:
                print "Writing datalayer '"+dataLayerName+"' to netCDF file as "+dataLayers[dataLayerName].netCDFName;
            variable = ncfile.createVariable(dataLayers[dataLayerName].netCDFName, dtype('float64').char, dims, fill_value=DataLayer.fill_value)
            data = dataLayers[dataLayerName].fdata; #fdata is usually a view by sometimes a copy so it has to be done this way. There is probably a better way to do this.
            data.shape = (dataLayers[dataLayerName].nx, dataLayers[dataLayerName].ny);
            variable[:] = data;
        except AttributeError:
            print "%s:No netCDFName or data attribute found in DataLayer '%s'." % (function, dataLayerName);
        except RuntimeError as e:
            print "\n%s: Error when writing datalayer '%s' the netCDF variable ('%s') couldn't be created because it already exists. Check for netCDFName clashes." % (function, dataLayerName, dataLayers[dataLayerName].netCDFName);
            print e.args;
        except ValueError as e:
            print "%s: Cannot reside datalayer '%s'" % (function, dataLayers[dataLayerName].name);
            print type(e), e.args;
        
        variable.missing_value = missing_value;
        variable.scale_factor = 1.0;
        variable.add_offset = 0.0;
        
        try:
            if dataLayers[dataLayerName].units != None:
                variable.units = dataLayers[dataLayerName].units;
        except AttributeError:
            print "%s: No units found for datalayer named '%s'." % (function, dataLayerName);
        
        try:
            if dataLayers[dataLayerName].minBound != None:
                variable.valid_min = dataLayers[dataLayerName].minBound;
        except AttributeError:
            print "%s: No minBound found for datalayer named '%s'." % (function, dataLayerName);
        
        try:
            if dataLayers[dataLayerName].maxBound != None:
                variable.valid_max = dataLayers[dataLayerName].maxBound;
        except AttributeError:
            print "%s: No maxBound found for datalayer named '%s'." % (function, dataLayerName);
        
        try:
            if dataLayers[dataLayerName].standardName != None:
                variable.standard_name = dataLayers[dataLayerName].standardName;
        except AttributeError:
            print "%s: No standardName found for datalayer named '%s'." % (function, dataLayerName);
        
        try:
            if dataLayers[dataLayerName].longName != None:
                variable.long_name = dataLayers[dataLayerName].longName;
        except AttributeError:
            print "%s: No longName found for datalayer named '%s'." % (function, dataLayerName);
   
    #set some global attributes
    setattr(ncfile, 'Conventions', 'CF-1.6') 
    setattr(ncfile, 'Institution', 'Originally developed by the partners of the ESA OceanFlux GHG and OceanFlux GHG Evolution projects. Now continued by the CarbonLab team at the University of Exeter.') 
    setattr(ncfile, 'Contact', 'email: j.d.shutler@exeter.ac.uk')
    
    #Output all the parameters used in this run.
    for paramName in vars(runParams).keys():
        paramValue = getattr(runParams, paramName);
        if paramValue is not None:
            if type(paramValue) is bool: #netCDF does not support bool types.
                paramValue = int(paramValue);
            setattr(ncfile, paramName, paramValue);
 
    ncfile.close();
 
    return 0;


#Calculates solubility of distilled water (i.e. assuming 0 salinity) given global temperature.
#overwrites solubilityDistilled with the calculated value.
#TODO: no need to pass nx, ny into all these functions.
def calculate_solubility_distilled(solubilityDistilled, salinity, rain_wet_deposition_switch, sstskin, deltaT, nx, ny):
    #First create a 0 salinity dataset
    salDistil = array([missing_value] * len(salinity))
    for i in arange(nx * ny):
        if (salinity[i] != missing_value):
            salDistil[i] = 0.0
        else:
            salDistil[i] = missing_value
    
    #Next calculate the solubility using the zero salinity 'distilled water' dataset
    solubilityDistilled = solubility(sstskin, salDistil, deltaT, nx, ny, True);
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


def add_noise(data, rmse_value, mean_value, no_elements):
# randomly adding noise to data array, based log normal distribution
# rmse value is used as the standard deviation of the noise function

   # intialise the random number generator
  for i in arange(no_elements):
     if ( (data[i] != missing_value) and (data[i] != 0.0) ):
      orig = data[i]
      value = log(data[i])
      stddev = float(rmse_value/orig)
      noise = normalvariate(0,stddev) # determines the random noise value based on the input data and the standard deviation of the uncertainty
      value = value + noise
      data[i] = exp(value)
  return data

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

             # determine the schmidt number

def schmidt(sstC_fdata, nx, ny,gas):
 # calculating the schmidt data

   sc_fdata = array([missing_value] * nx*ny)
   if 'o2' in gas.lower():
       for i in arange(nx * ny):
          if (sstC_fdata[i] != missing_value):
             sc_fdata[i] = 1953.4 - (128.0 * sstC_fdata[i]) + (3.9918 * (sstC_fdata[i] * sstC_fdata[i])) - (0.050091 * (sstC_fdata[i] * sstC_fdata[i] * sstC_fdata[i])) + (0.00093777 * (sstC_fdata[i] * sstC_fdata[i] * sstC_fdata[i] * sstC_fdata[i]))#IGA-O2
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


def solubility(sstK, sal, deltaT, nx, ny, flux_calc):
 # solubility calculation
 # equation from Table A2 of Wanninkkhof, JGR, 1992
   sol = array([missing_value] * nx*ny)
   for i in arange(nx * ny):
      if ( (sstK[i] != missing_value) and (sal[i] != missing_value) and (sstK[i] > 0.0) ):
         sol[i] = -60.2409 + ( 93.4517*(100.0 / sstK[i]) ) + (23.3585 * (log(sstK[i]/100.0))) + (sal[i] * (0.023517 + ( (-0.023656)*(sstK[i]/100.0)) + (0.0047036*( (sstK[i]/100.0)*(sstK[i]/100.0) ) ) ) )
         sol[i] = exp(sol[i])
             #runParams.flux_calc is a switch to remove Delta_T-Sb component - ie selects use of RAPID or EQUILIBRIUM flux models from Woolf et al., 2016
         if flux_calc != 2:
            deltaT[i] = 0.0
         sol[i] = sol[i] * (1 - (0.015*deltaT[i]))
      else:
         sol[i] = missing_value
   return sol


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


#data is a matrix, nx and ny are the output dimensions (must be a factor of the matrix dimensions).
def average_pixels(data, nx, ny, missing_value):
    '''averages arr into superpixels each consisting of the mean of a n x n
    window in arr. Seems to be intended to go from a 0.5 x 0.5 degree grid to a
    1 x 1 degree grid, in which case n must be set to 2. Checks for and ignores
    values equal to missing_value.'''
    function = "(average_pixels, main)"
    
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
        print datanx, datany
        print "here:", datanx%nx, datany%ny
        raise ValueError("%s: Cannot rescale data layer because global data layer dimensions are not a whole multiple of the data layer's dimensions." % function);
    n = (datanx/nx);
    n2 = (datany/ny);
    if n!=n2:
        raise ValueError("%s: Scaling data layers by irregular scaling factors is not supported (e.g. %d != %d)." % (function, n, n2));
    
    
    print "%s Averaging sstgrad_fdata into 1x1 degree grid (N=%d)" % (function, n)
    nj0, ni0 = data.shape
    nj, ni = [nj0 / n, ni0 / n]
    if nj * n != nj0 or ni * n != ni0:
       print "Dimensions ", nj0, ni0, " indivisible by ", n;
    out = resize(missing_value, [nj, ni])
    for j in xrange(nj):
       j0 = j * n
       j1 = j0 + n
       for i in xrange(ni):
          i0 = i * n
          a = data[j0:j1, i0:i0 + n]
          w = nonzero(a != missing_value)
          if size(w) > 0:
             out[j, i] = mean(a[w])
    #print out;
    return out

#Returns false if dataLayer dimensions do not match the reference dimensions.
def check_dimensions(dataLayer, ref_nx, ref_ny, DEBUG=False):
   function = "(check_dimensions, main)"
   if dataLayer.nx == ref_nx and dataLayer.ny == ref_ny:
      if DEBUG:
         print "\n%s Input data (%s) have identical dimensions to reference values (%s, %s) "% (function, dataLayer.name, dataLayer.nx, dataLayer.ny)
         return True;
   else:
      print "\n%s Input data ('%s') dimensions are non-identical to reference (new: %s, %s is not equal to: %s, %s)." % (function, dataLayer.name, dataLayer.nx, dataLayer.ny, ref_nx, ref_ny)
      return False;


class FluxEngine: 
    def __init__(self, parameterDict):
        self.runParams = RunParameters();
        self.runParams.set_parameters(parameterDict);
        self.defaultSettings = Settings(path.join(self.runParams.src_home, "settings.xml")); #Load default settings metadata for datalayers.
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
        function = "(ofluxghg_flux_calc, FluxEngine.add_data_layer)";
        
        try:
            self._add_single_data_layer(name, infile, prod, transposeData, preprocessing=preprocessing);
            
        except KeyError:
            print "\n%s: Data variable '%s' is missing from %s input (%s)" % (function, prod, name, infile);
            return False;
        except ValueError as e: #E.g. incorrect number of dimensions
            print "\n%s: %s" % (function, e.args);
            return False;
        
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
                print "here", type(e), e.args;
                print "%s: Did not find stddev variable ('%s') in netCDF file for %s." % (function, stddevProd, name);
            except ValueError as e: #E.g. incorrect number of dimensions
                print "\n%s: %s" % (function, e.args);
                return False;
        
        if countProd != None:
            if countProd == "auto": #automatically detect count product
                countProd = prod.rsplit("_", 1)[0]+"_count";
            try:
                self._add_single_data_layer(name+"_count", infile, countProd, transposeData=transposeData);
                copy_missing_values(self.data[name], self.data[name+"_count"]); #Any missing values in 'name' should be copied to the count dataset for consistency.
            except KeyError:
                print "%s: Did not find stddev variable ('%s') in netCDF file for %s." % (function, countProd, name);
            except ValueError as e: #E.g. incorrect number of dimensions
                print "\n%s: %s" % (function, e.args);
                return False;
            
        #Datalayer was successfully added.
        return True;
    
    #Adds a single datalayer and acquires metadata
    #TODO: transposeData is now be handled by preprocessing, so this cna be removed...
    def _add_single_data_layer(self, name, infile, prod, transposeData=False, preprocessing=None):
        function = "(ofluxghg_flux_calc, FluxEngine._add_single_data_layer)";
        
        metaData = self._extract_data_layer_meta_data(name);
        
        try:
            dl = DataLayer.create_from_file(name, infile, prod, metaData, transposeData=transposeData, preprocessing=preprocessing);
            self.data[name] = dl;
            
            #If this is the first datalayer to be added, use this to set the nx and ny dimensions
            #   (all other datalayers will be checked against this to ensure conformity)
            #if len(self.data) == 1:
            #if name == "sstskin":
            #    self.nx = dl.nx;
            #    self.ny = dl.ny;
            #    print "Using %s to infer grid dimensions (%d, %d)." % (name, self.ny, self.nx);
            
        except IOError as e:
            print "\n%s: %s inputfile %s does not exist" % (function, name, infile)
            print e.args;
            return False;
    
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
        #for key in vars(self.runParams).keys():
        #    print key;
        
        #function = "(ofluxghg_flux_calc, FluxEngine.run)";
        status = self._check_datalayers(); #Check for consistency of datalayers (e.g. dimensions all match one another).
        if status == True:
            return self._run_fluxengine(self.runParams);
        else:
            print "Not all datalayers are consistent (check dimensions of data).";
            return status;
    
    #Reads longitude, latitude, time and dimension sizes (nx, ny).
    def _load_lon_lat_time(self):
        function = "(ofluxghg_flux_calc, FluxEngine._load_lon_lat_time)";
        
        #Find the path of the relevent datalayer.
        try:
            axesDatalayerInfile = getattr(self.runParams, self.runParams.axes_data_layer+"_infile");
        except ValueError as e:
            print "Couldn't find file path for axes_data_layer (%s). Check that this is correctly defined in the config file (it must match the name of another datalayer)." % self.runParams.axes_data_layer;
            print e.args;
        
        #Read longitude, latitude and time data from the sstskin infile, so open this file.
        try:
            dataset = Dataset(axesDatalayerInfile);
        except IOError as e:
            print "\n%s: axes_data_layer (%s) inputfile %s does not exist" % (function, self.runParams.axes_data_layer, axesDatalayerInfile);
            print type(e), e.args;
        
        #Read lat and lat
        try:
            self.latitude_data = dataset.variables[self.runParams.latitude_prod][:];
            self.longitude_data = dataset.variables[self.runParams.longitude_prod][:];
            if self.latitude_data[0]<0: #IGA - it is a vector that is in opposite orientation to 'taka'
                self.latitude_data = flipud(self.latitude_data);
        except KeyError as e:
            print "%s: Couldn't find longitude (%s%) and/or latitude (%s) variables in %s." % (function, self.runParams.latitude_prod, self.runParams.longitude_prod, self.runParams.sstskin_infile);

        #Determine if already a grid, if not calculate lon and lat grids.
        if len(self.latitude_data.shape) == 1: #not already a grid
            self.latitude_grid = meshgrid(ones((len(self.longitude_data))), self.latitude_data)[1]#IGA - gridding vectored geo data so that they can be treated the sam way as non-regular grids
            self.longitude_grid = meshgrid(self.longitude_data, ones((len(self.latitude_data))))[1]#IGA
        else: #IGA - geo data are already a grid
            self.latitude_grid = self.latitude_data;
            self.longitude_grid = self.longitude_data;
            
        #read time
        try:
            self.time_data = dataset.variables[self.runParams.time_prod][:];
        except KeyError as e:
            print "%s: Couldn't find time (%s%) variables in %s." % (function, self.runParams.time_prod, self.runParams.sstskin_infile);

        #set dimensions
        self.ny, self.nx = self.latitude_grid.shape;
        print "%s: Grid dimensions set to: (%d, %d)" % (function, self.ny, self.nx);

    #Ran after each datalayer is read in. Checks for consistency between datalayers.
    #Rescales data layers which can be rescaled.
    #returns true if all data layers are successfully validated.
    def _check_datalayers(self):
        #Check that dimensions match
        for key in self.data:
            if check_dimensions(self.data[key], self.nx, self.ny, DEBUG) == False:
                #Dimensions don't match, so try to rescale it.
                try:
                    print "Attempting to rescale.";
                    self.data[key].data = average_pixels(self.data[key].data, self.nx, self.ny, DataLayer.missing_value);
                    self.data[key].ny, self.data[key].nx = self.data[key].data.shape;
                    self.data[key].calculate_fdata();
                except ValueError as e:
                    print e.args;
                    return False;
        return True;


    def _run_fluxengine(self, runParams):
        function = "(ofluxghg_flux_calc, FluxEngine._run_fluxengine_)";
        #return 0;
        
        #Set up logging object
        try:
            logger = logging.getLogger('FluxEngine_debug_log');
            #hdlr = logging.FileHandler(os.path.join(workingDirectory, runParams.LOG_PATH), filemode='w')
            hdlr = logging.FileHandler(runParams.LOG_PATH);
            formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s');
            hdlr.setFormatter(formatter);
            logger.addHandler(hdlr);
            if DEBUG_LOGGING == True:
                logger.setLevel(logging.DEBUG);
            else:
                logger.setLevel(logging.DEBUG);
        except:
            print "\n%s Couldn't initialise logger at path %s" % (function, runParams.LOG_PATH);
    
        if DEBUG:
           print "\n%s Flux model: %d, k_parmaterisation: %d" % (function, runParams.flux_calc, runParams.k_parameterisation);
           print "\n%s Using the following files %s (sstskin), %s (sstfnd) %s (windu10) %s (sal) %s (rain) %s (bio) %s (sstgrad)" % (function, runParams.sstskin_infile, runParams.sstfnd_infile, runParams.windu10_infile, runParams.salinity_infile, runParams.rain_infile, runParams.biology_infile, runParams.sstgrad_infile);
           print "\n%s Using the following products %s (sstskin), %s (sstfnd) %s (windu10) %s (msl) %s (sal) %s (rain) %s (bio) %s (sstgrad)" % (function, runParams.sstskin_prod, runParams.sstfnd_prod, runParams.windu10_prod, runParams.pressure_prod, runParams.salinity_prod, runParams.rain_prod, runParams.biology_prod, runParams.sstgrad_prod);
        
        #TODO: Replace directly with self.nx, self.ny, no need to use local variables here.
        nx = self.nx;
        ny = self.ny;

        ### Adding empty data layers for data that are computed later.
        #If there isn't any sstfnd data create empty arrays to fill later. (i.e. sstfnd = sstskin + cool_skin_difference);
        if ("sstfnd" not in self.data) and (runParams.use_sstfnd_switch == 0):
            self.add_empty_data_layer("sstfnd");
            self.add_empty_data_layer("sstfnd_stddev");
            self.add_empty_data_layer("sstfnd_count");
        if "sstskin_stddev" not in self.data and "sstskin_stdfnd" not in self.data: #Remove this, stddev and count aren't used except when they are!
            self.add_empty_data_layer("sstskin_stddev");
            self.add_empty_data_layer("sstskin_count");

        
        if runParams.TAKAHASHI_DRIVER == True:
           #need to generate the moment2 and moment3 data
           self.add_empty_data_layer("windu10_moment2");
           self.add_empty_data_layer("windu10_moment3");
           for i in arange(self.nx * self.ny):
              if self.data["windu10"].fdata[i] != missing_value:
                 self.data["windu10_moment2"].fdata[i] = self.data["windu10"].fdata[i]*self.data["windu10"].fdata[i]
                 self.data["windu10_moment3"].fdata[i] = self.data["windu10"].fdata[i]*self.data["windu10"].fdata[i]*self.data["windu10"].fdata[i]
        
        #If not using takahashi then we need to generate pco2_air, so create an empty data layer in preparation.
            #No longer needed as this isn't automatically outputted.
        #if TAKAHASHI_DRIVER == False:
        #   self.add_empty_data_layer("pco2_air");
        
        #SOCATv4 - using input foundation temperature as the SST temp-------------START
        #If there is no pco2_sst data, we need to get it / generate it.
        if "pco2_sst" not in self.data: #SOCATv4
            try:
                print "No pco2_sst data was supplied."
                self.add_empty_data_layer("pco2_sst");
                self.data["pco2_sst"].fdata = self.data["sstfnd"].fdata-273.15; #copy/convert sstfnd
            except (IOError, KeyError, ValueError) as e:
                print "pco2_sst data not available and could read sstfnd so cannot proceed.";
                print type(e), "\n"+e.args;
                return 1;   
        
        #some specific pco2 conditions
        #TODO: This shouldn't be randomly here.
        self.data["pco2_sw"].fdata[abs(self.data["pco2_sw"].fdata)<0.1]=DataLayer.missing_value#IGA_SOCATv4

        if (runParams.pco2_data_selection != 1):
            #signifies that we're NOT using SOCAT data, which means there is no stddev data for pco2_sw
            self.add_empty_data_layer("pco2_sw_stddev");
            
        #initialising some data structures for the calculations
        self.add_empty_data_layer("scskin");
        self.add_empty_data_layer("scfnd");
        self.add_empty_data_layer("solubility_skin");
        self.add_empty_data_layer("solubility_fnd");
        self.add_empty_data_layer("pH2O");
        self.add_empty_data_layer("pco2_air_cor");
        self.add_empty_data_layer("pco2_sw_cor");
        
        #apply additive saline skin value, but why selectively apply it?
        self.add_empty_data_layer("salinity_skin");
        for i in arange(self.nx * self.ny):
            if (self.data["salinity"].fdata[i] >= 0.0) and (self.data["salinity"].fdata[i] <= 50.0):
                if ( (self.data["salinity"].fdata[i] + runParams.saline_skin_value) <= 50.0):
                    self.data["salinity_skin"].fdata[i] = self.data["salinity"].fdata[i] + runParams.saline_skin_value
            else:
                self.data["salinity"].fdata[i] = missing_value
                self.data["salinity_skin"].fdata[i] = missing_value
        if (runParams.saline_skin_value != 0.0):
            print "%s Using the saline skin model (%lf psu added to skin salinities)" % (function, runParams.saline_skin_value)
        
        #conversion of rain data from mm day-1 to mm hr^-1
        if "rain" in self.data:
            for i in arange(self.nx * self.ny):   
               if (self.data["rain"].fdata[i] != DataLayer.missing_value):
                   self.data["rain"].fdata[i] /= 24.0;
        
        
        #interpreting fnd_data option
        cool_skin_difference = 0.17; #0.16; #Relationship and value from Donlon et al., 2002; page 358, first paragraph
        if runParams.TAKAHASHI_DRIVER == True:
           cool_skin_difference = 0.0;
        
        ####Derive sstskin as necessary.
        #Two possible datasets: sstskin and sstfnd, from different parts of the water column.
        #One option: sst gradients
        #If using only one dataset then copy that over the other. Otherwise keep both.
        #cool_skin_difference is the assumed temperature difference between the skin and foundation layer
        
        #using sstfnd, so copy sstfnd into sstskin
        #sstskin = sstfnd
        if runParams.sst_gradients_switch == 0 and runParams.use_sstskin_switch == 0 and runParams.use_sstfnd_switch == 1:
           print "%s SST gradient handling is off, using SSTfnd data selection in configuration file for all components of the flux calculation (this will ignore any SSTskin data in configuration file)." % (function)
           #actually copy sstfnd data into the sstskin dataset to make sure
           for i in arange(nx * ny):
               if self.data["sstfnd"].fdata[i] != missing_value:
                   self.data["sstskin"].fdata[i] = self.data["sstfnd"].fdata[i]
                   #if "sstskin_stddev" in self.data and "sstskin_count" in self.data:
                   self.data["sstskin_stddev"].fdata[i] = self.data["sstfnd_stddev"].fdata[i]
                   self.data["sstskin_count"].fdata[i] = self.data["sstfnd_count"].fdata[i]
               else:
                   self.data["sstskin"].fdata[i] = missing_value;
        
        
        #IGA added for the case where only foundation is provided and gradients are on------------------------------
        #must estimate sstskin (= sstfnd - cool_skin_difference)
        elif runParams.sst_gradients_switch == 1 and runParams.use_sstskin_switch == 0 and runParams.use_sstfnd_switch == 1:
            print "%s Using SSTfnd data selection with correction for skin temperature (SSTskin = SSTfnd - 0.16)(ignoring SSTskin data in configuration file)." % (function)
            #actually copy sstfnd data into the sstskin dataset to make sure
            if "sstskin" not in self.data: #Must add the sstskin layer first!
                self.add_empty_data_layer("sstskin");
                if "sstfnd_stddev" in self.data:
                    self.add_empty_data_layer("sstskin_stddev");
                if "sstfnd_count" in self.data:    
                    self.add_empty_data_layer("sstskin_count");
            for i in arange(nx * ny): #sstdkin = sstfnd - cool_skin_difference
                if self.data["sstfnd"].fdata[i] != missing_value:
                    self.data["sstskin"].fdata[i] = self.data["sstfnd"].fdata[i]-cool_skin_difference
                    self.data["sstskin_stddev"].fdata[i] =  self.data["sstfnd_stddev"].fdata[i]
                    self.data["sstskin_count"].fdata[i] = self.data["sstfnd_count"].fdata[i]
                else:
                    self.data["sstskin"].fdata[i] = missing_value;
        
        #Using sstskin, so calculate it from sstfnd.
        elif runParams.sst_gradients_switch == 0 and runParams.use_sstskin_switch == 1 and runParams.use_sstfnd_switch == 0:
           print "%s SST gradient handling is off, using SSTskin to derive SSTfnd (SSTfnd = SSTskin + 0.17) for flux calculation (ignoring SSTfnd data in configuration file)." % (function)
           #setting sstfnd_ data fields to skin values
           for i in arange(nx * ny):
              if self.data["sstskin"].fdata[i] != missing_value:
                 
                 self.data["sstfnd"].fdata[i] = self.data["sstskin"].fdata[i] + cool_skin_difference
                 self.data["sstskin"].fdata[i] = self.data["sstskin"].fdata[i] + cool_skin_difference
                 self.data["sstfnd_stddev"].fdata[i] = self.data["sstskin_stddev"].fdata[i]
                 self.data["sstfnd_count"].fdata[i] = self.data["sstskin_count"].fdata[i]
              else:
                 self.data["sstfnd"].fdata[i] = missing_value
                 self.data["sstskin"].fdata[i] = missing_value
                 
        elif runParams.sst_gradients_switch == 1 and runParams.use_sstskin_switch == 1 and runParams.use_sstfnd_switch == 0:
           print "%s SST gradient handling is on, using SSTskin and SSTfnd = SSTskin + 0.16 for flux calculation (ignoring SSTfnd data in configuration file)." % (function)
            #setting sstfnd_ data fields to skin values  
           for i in arange(nx * ny):
              if self.data["sstskin"].fdata[i] != missing_value:
                 self.data["sstfnd"].fdata[i] = self.data["sstskin"].fdata[i] + cool_skin_difference
                 self.data["sstfnd_stddev"].fdata[i] =  self.data["sstskin_stddev"].fdata[i]
                 self.data["sstfnd_count"].fdata[i] = self.data["sstskin_count"].fdata[i]
              else:
                 self.data["sstfnd"].fdata[i] = missing_value
        elif runParams.sst_gradients_switch == 0 and runParams.use_sstskin_switch == 1 and runParams.use_sstfnd_switch == 1:
           print "%s SST gradient handling is off, using SSTfnd and SSTskin from the configuration file." % (function)        
        elif runParams.sst_gradients_switch == 1 and runParams.use_sstskin_switch == 1 and runParams.use_sstfnd_switch == 1:
           print "%s SST gradient handling is on, using SSTfnd and SSTskin from the configuration file." % (function)
        else:
           print "\n%s sst_gradients_switch (%d), use_sstskin_switch (%d) and use_sstfnd_switch (%d) combination in configuration not recognised, exiting." % (function, runParams.sst_gradients_switch, runParams.use_sstskin_switch, runParams.use_sstfnd_switch)
           return 1;
       

        #quality filtering and conversion of SST datasets
        #Calculate sstskinC and sstfndC
        self.add_empty_data_layer("sstskinC");
        self.add_empty_data_layer("sstfndC");
        #if TAKAHASHI_DRIVER != True:
        for i in arange(nx * ny):
            if self.data["sstskin"].fdata[i] != missing_value:
                self.data["sstskinC"].fdata[i] = self.data["sstskin"].fdata[i] - 273.15;
                self.data["sstfndC"].fdata[i] = self.data["sstfnd"].fdata[i] - 273.15;
        
#           for i in arange(nx * ny):
#              # convert SSTfnd to centrigrade and store the result
#              if (self.data["sstfnd"].fdata[i] != missing_value):
#                 self.data["sstfndC"].fdata[i] = self.data["sstfnd"].fdata[i] - 273.15
#              else:
#                 self.data["sstfndC"].fdata[i] = missing_value
#              # convert SSTskin to centrigrade and store the result
#              if (self.data["sstskin"].fdata[i] != missing_value):
#                 self.data["sstskinC"].fdata[i] = self.data["sstskin"].fdata[i] - 273.15
#              else:
#                 self.data["sstskinC"].fdata[i] = missing_value
        
         # quality filtering if using TAKAHASHI_DRIVER
         # TAKAHASHI data are in oC (and not Kelvin), so we have to apply the reverse
#        if TAKAHASHI_DRIVER == True:
#           for i in arange(nx * ny):
#              if (self.data["sstfnd"].fdata[i] != missing_value):
#                 self.data["sstfndC"].fdata[i] = self.data["sstfnd"].fdata[i]
#                 self.data["sstfnd"].fdata[i] = self.data["sstfndC"].fdata[i] + 273.15
#              else:
#                 self.data["sstfndC"].fdata[i] = missing_value
#        
#              if (self.data["sstskin"].fdata[i] != missing_value):
#                 self.data["sstskinC"].fdata[i] = self.data["sstskin"].fdata[i]
#                 self.data["sstskin"].fdata[i] = self.data["sstskinC"].fdata[i] + 273.15
#              else:
#                 self.data["sstskinC"].fdata[i] = missing_value
        
         # ability to randomly perturb the input datasets
         # needed for the ensemble analyses
         # stddev of noise is using published RMSE for each dataset
         # all datasets are considered to have bias=0, hence mean of noise=0
        if (runParams.random_noise_windu10_switch == 1):
           add_noise(self.data["windu10"].fdata, 0.44, 0.0, nx*ny)
           add_noise(self.data["windu10_moment2"].fdata, 0.44, 0.0, nx*ny)
           add_noise(self.data["windu10_moment3"].fdata, 0.44, 0.0, nx*ny)
           print "%s Adding random noise to windu10_mean, windu10_moment2 and windu10_moment3 (mean 0.0, stddev 0.44 ms^-1 - assuming using ESA GlobWave data)" % (function)
        
        if (runParams.random_noise_sstskin_switch == 1):
           add_noise(self.data["sstskin"].fdata, 0.14, 0.0, nx*ny)
           print "%s Adding random noise to sstskin (mean 0.0, stddev 0.14 ^oC - assuming using ESA CCI ARC data)" % (function)
        
        if (runParams.random_noise_sstfnd_switch == 1):
           add_noise(self.data["sstfnd"].fdata, 0.6, 0.0, nx*ny)
           print "%s Adding random noise to sstfnd (mean 0.0, stddev 0.6 ^oC - assuming using OSTIA data)" % (function)
        
        if (runParams.random_noise_pco2_switch == 1):
           print "/n%s Shape of pco2 data",self.data["pco2_sw"].fdata.shape
           add_noise(self.data["pco2_sw"].fdata, 6, 0.0, nx*ny)
           print "%s Adding random noise to pco2/fco2 data (mean 0.0, stddev 6 uatm - Using Candyfloss data, value provided by J Shutler)" % (function)
           #print "%s Adding random noise to pco2/fco2 data (mean 0.0, stddev 2.0 uatm - assuming using SOCAT flag A and flag B data)" % (function)
        
        #Ians rain noise test
        #if (random_noise_rain == 1):
        # self.data["rain"].data = reshape(self.data["rain"].data,self.data["rain"].data.shape[0]*self.data["rain"].data.shape[1])
        # rain_err_data = reshape(rain_err_data,rain_err_data.shape[0]*rain_err_data.shape[1])
        # add_noise(self.data["rain"].data, rain_err_data, 0.0, rain_nx*rain_ny)
        # print "%s Adding random noise to rain data using variance field" % (function)
        
        #Ians rain noise test end
        
         # bias values to be added here
        if (runParams.bias_windu10_switch == 1):
           add_bias(self.data["windu10"].fdata, runParams.bias_windu10_value, nx*ny)
            # makes no sense to add bias to second and third order moments, as any bias in the system 
            # would only impact on the mean (which is the 1st order moment)
           print "%s Adding bias to windu10_mean (not to second and third order moments) (value %lf ms^-1)" % (function, runParams.bias_windu10_value)
        
        if (runParams.bias_sstskin_switch == 1):
           add_bias(self.data["sstskin"].fdata, runParams.bias_sstskin_value, nx*ny)
           print "%s Adding bias noise to sstskin (value %lf ^oC)" % (function, runParams.bias_sstskin_value)
        
        if (runParams.bias_sstfnd_switch == 1):
           add_bias(self.data["sstfnd"].fdata, runParams.bias_sstfnd_value, nx*ny)
           print "%s Adding bias noise to sstfnd (value %lf ^oC)" % (function, runParams.bias_sstfnd_value)
        
        if (runParams.bias_pco2_switch == 1):
           add_bias(self.data["pco2_sw"].fdata, runParams.bias_pco2_value, nx*ny)
           print "%s Adding bias noise to pco2w (value %lf uatm)" % (function, runParams.bias_pco2_value)
        
        
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
                print "Converted pressure data to mbar from Pa"
        # ------------------------------------------------------IGA_SOCATv4 - Was removed as using pressure data in mbar.-------------END
        # -                                                                                                                                     -
        # ---------------------------------------------------------------------------------------------------------------------------------------
        
        #########################################################
        # calculations for components of the flux calculation
        #########################################################
         
        # year correction for pCO2/fCO2 data
        # if runParams.pco2_data_selection == 1, then using SOCAT data so no increment is added as SOCAT are normalised to 2010
        # if runParams.pco2_data_selection == 0, then using Takahashi et al data which are normalised to 2000
        if runParams.pco2_data_selection == 0:
        # increment to be added to pCO2w and pCO2a in Takahashi climatology
        # 1.5uatm per year since 2000 as set by Takahashi et al 2009.
        # Takahashi data are normalised to the year 2000
           pco2_increment = (runParams.year - 2000.0) * 1.5
           print "\n%s year: %d Takahashi pCO2_increment: %lf (uatm) (runParams.pco2_data_selection: %d)" % (function, runParams.year, pco2_increment,runParams.pco2_data_selection)
        elif (runParams.pco2_data_selection == 1): # signifies that SOCAT pco2 data are being used
           # need handling for SOCAT data for years before 2010
           # 1.5uatm per year since 2000 as set by Takahashi et al 2009.
           # so assuming -1.5uatm for each year prior to 2010
           pco2_increment = (runParams.year - 2010.0) * 1.5
           print "\n%s year: %d SOCAT pCO2_increment: %lf (uatm) (runParams.pco2_data_selection: %d)" % (function, runParams.year, pco2_increment,runParams.pco2_data_selection)
        elif (runParams.pco2_data_selection == 2): # signifies that SOCAT fco2 data are being used
           # need handling for SOCAT data for years before 2010
           # 1.5uatm per year since 2000 as set by Takahashi et al 2009.
           # so assuming -1.5uatm for each year prior to 2010
           pco2_increment = (runParams.year - 2010.0) * 1.5
           print "\n%s year: %d SOCAT fCO2_increment: %lf (uatm) (runParams.pco2_data_selection: %d)" % (function, runParams.year, pco2_increment,runParams.pco2_data_selection)
        elif (runParams.pco2_data_selection == 4): # signifies that in situ or time-series data are being used - no increment
           pco2_increment = 0.0
        elif (runParams.pco2_data_selection == 45): # signifies that SOCATv4 climatology is being used - increments apply
           pco2_increment = (runParams.year - 2010.0) * 1.5
           print "\n%s year: %d SOCATv4 climatology fCO2_increment: %lf (uatm) (runParams.pco2_data_selection: %d)" % (function, runParams.year, pco2_increment,runParams.pco2_data_selection)
        elif (runParams.pco2_data_selection == 3): # signifies that in situ or time-series data are being used - no increment
           pco2_increment = 0.0
           print "\ng%s year: %d insitu fCO2_increment: %lf (uatm) (runParams.pco2_data_selection: %d)" % (function, runParams.year, pco2_increment,runParams.pco2_data_selection)
        else:
           pco2_increment = 0.0
            
        DeltaT_fdata = array([missing_value] * nx*ny)
        
        if runParams.flux_calc == 1:
           print "%s Using the RAPID model (from Woolf et al., 2016)" % (function)
        elif runParams.flux_calc == 2:
           print "%s Using the EQUILIBRIUM model (from Woolf et al., 2016)" % (function)
           for i in arange(nx * ny):
              if ( (self.data["sstfndC"].fdata[i] != missing_value) & (self.data["sstskinC"].fdata[i] != missing_value) ):
                 DeltaT_fdata[i] = self.data["sstfndC"].fdata[i] - self.data["sstskinC"].fdata[i]
              else:
                 DeltaT_fdata[i] = missing_value
        elif runParams.flux_calc == 3:
            print "%s Using the BULK model (ie Flux = k*S*delta_pCO2)" % (function)
            print "%s Note: this assumes that the SSTskin dataset is the only temperature dataset, and so this is what will be used to calculate k and the solubility" % (function)
        else:
           print "\n%s runParams.flux_calc from configuration not recognised, exiting." % (function)
           return 1;
        
        # calculating the schmidt number at the skin and fnd
        self.data["scskin"].fdata = schmidt(self.data["sstskinC"].fdata, nx, ny, runParams.GAS)
        self.data["scfnd"].fdata = schmidt(self.data["sstfndC"].fdata, nx, ny, runParams.GAS)
        
         # calculating the skin solubility, using skin sst and salinity
        self.data["solubility_skin"].fdata = solubility(self.data["sstskin"].fdata, self.data["salinity_skin"].fdata, DeltaT_fdata, nx, ny, True)
        
         # calculating the interfacial solubility
        self.data["solubility_fnd"].fdata = solubility(self.data["sstfnd"].fdata, self.data["salinity"].fdata, DeltaT_fdata, nx, ny, runParams.flux_calc)

         # calculate pCO2 data using mean sea level pressure data
         # equation 26, Kettle et al, 2009, ACP
         # pCO2_air = X[CO2] ( P(t) - pH2O(t) )
         # where X[CO2] is from Takahashi climatology
         # pH2O(t) = 1013.25 exp (24.4543 - 67.4509 (100/Tk(t) - 4.8489 ln (Tk(t) / 100) - 0.000544 S)
         # S = salinity, Tk = temperature in Kelvin
         # pCO2_water = pCO2_water_tak exp (0.0423 (T_foundation - T_tak)
         # T_tak = Takahashi temperature
        
         # debuggin differences in pH20 values
        pH2O_takahashi_fdata = array([missing_value] * nx*ny)
        humidity_fdata = array([missing_value] * nx*ny)
        pH2O_diff_fdata = array([missing_value] * nx*ny)
        pCO2a_diff_fdata = array([missing_value] * nx*ny)
        #dpCO2_diff_fdata = array([missing_value] * nx*ny) #TODO: Should this be calculated/outputted?
        b11_fdata = array([missing_value] * nx*ny)
        d12_fdata = array([missing_value] * nx*ny)

        #TODO: This should be handled with options in the config file        
        # always using XCO2 from 2000 for SOCAT data
        if runParams.pco2_data_selection == 1 or runParams.pco2_data_selection == 2 or runParams.pco2_data_selection == 45:
           #pco2_increment_air = (runParams.year - 2000.0) * 1.5 #TODO: Changed.
           pco2_increment_air = pco2_increment;
           print "%s year: %d fCO2/pCO2_increment_air: %lf (uatm)" % (function, runParams.year, pco2_increment_air)
        else:
           pco2_increment_air = pco2_increment
           print "%s year: %d fCO2/pCO2_increment_air: %lf (uatm)" % (function, runParams.year, pco2_increment_air)
        
        sys.stdout.flush();
        
#        salskin_fc = 0;
#        sstskinK_fc = 0;
#        pres_fc = 0;
#        vco2_air_fc = 0;
#        sstfndK_fc = 0;
#        sstpco2_fc = 0;
#        pco2_sw_fc = 0;
#        sstskinK_fc0 = 0;
        
        for i in arange(nx * ny):
#            salskin_fc += (ma.is_masked(self.data["salinity_skin"].fdata[i]) == False) and (self.data["salinity_skin"].fdata[i] != missing_value);
#            sstskinK_fc += (ma.is_masked(self.data["sstskin"].fdata[i]) == False) and (self.data["sstskin"].fdata[i] != missing_value);
#            pres_fc += (ma.is_masked(self.data["pressure"].fdata[i]) == False) and (self.data["pressure"].fdata[i] != missing_value);
#            vco2_air_fc += (ma.is_masked(self.data["vco2_air"].fdata[i]) == False) and (self.data["vco2_air"].fdata[i] != missing_value);
#            sstfndK_fc += (ma.is_masked(self.data["sstfnd"].fdata[i]) == False) and (self.data["sstfnd"].fdata[i] != missing_value);
#            sstpco2_fc += (ma.is_masked(self.data["pco2_sst"].fdata[i]) == False) and (self.data["pco2_sst"].fdata[i] != missing_value);
#            pco2_sw_fc += (ma.is_masked(self.data["pco2_sw"].fdata[i]) == False) and (self.data["pco2_sw"].fdata[i] != missing_value);
#            sstskinK_fc0 += (ma.is_masked(self.data["sstskin"].fdata[i]) == False) and (self.data["sstskin"].fdata[i] != 0.0);
        
            
            if ( (self.data["salinity_skin"].fdata[i] != missing_value) and (self.data["sstskin"].fdata[i] != missing_value) and (self.data["pressure"].fdata[i] != missing_value) and (self.data["vco2_air"].fdata[i] != missing_value) and (self.data["sstfnd"].fdata[i] != missing_value) and (self.data["pco2_sst"].fdata[i] != missing_value) and (self.data["pco2_sw"].fdata[i] != missing_value) and (self.data["sstskin"].fdata[i] !=0.0) ):
            #if ( (self.data["salinity_skin"].fdata[i] != missing_value) and (self.data["sstskin"].fdata[i] != missing_value) and (self.data["pressure"].fdata[i] != missing_value) and (self.data["vco2_air"].fdata[i] != missing_value) and (self.data["sstfnd"].fdata[i] != missing_value) and (self.data["pco2_sw"].fdata[i] != missing_value) and (self.data["sstskin"].fdata[i] !=0.0) ):
                #print "1self.data["sstskin"].fdata: (%d,%d) %d %f log:%f\n" %(nx, ny,i,self.data["sstskin"].fdata[i]/100.0,log(self.data["sstskin"].fdata[i]/100.0))
                self.data["pH2O"].fdata[i] = 1013.25 * exp(24.4543 - (67.4509 * (100.0/self.data["sstskin"].fdata[i])) - (4.8489 * log(self.data["sstskin"].fdata[i]/100.0)) - 0.000544 * self.data["salinity_skin"].fdata[i])
                #print "2self.data["sstskin"].fdata: %d %f log:%f\n" %(i,self.data["sstskin"].fdata[i]/100.0,log(self.data["sstskin"].fdata[i]/100.0))
                  
                # this may be needed when using SMOS salinity data
                # To-DO: awaiting info from Lonneke and David before implementing fully
                pCO2_salinity_term = 0
                #dSalinity = salinity_rmse
                #if (salinity_option == 1):
                #  # need to determine dS/S using original salinity (its been modified above to add the salinity_rmse)
                #  # so (self.data["salinity"].fdata[i] - salinity_rmse) ensures that we are dealing with the original value of salinity
                  #  # using Ys=1 as a global correction following Sarmiento and Gruber, 2006)
                #pCO2_salinity_term = 1.0*(dSalinity/(self.data["salinity"].fdata[i] - salinity_rmse) ) #TH: Commented this out!
                #else:
                # pCO2_salinity_term = 0.0
                  
                # correction to different years, correction is data and year specific.
                # note for 2010, correction for SOCAT isn't strictly required. However the contents of the exponential will collapse
                # to 1 (with some rounding error expected), so effectively no correction will be applied
                if runParams.GAS == 'CO2' and runParams.pco2_data_selection != 3:
                    self.data["pco2_sw_cor"].fdata[i] = pco2_increment + (self.data["pco2_sw"].fdata[i] *exp( (0.0423*(self.data["sstfndC"].fdata[i] - self.data["pco2_sst"].fdata[i])) - (0.0000435*((self.data["sstfndC"].fdata[i]*self.data["sstfndC"].fdata[i]) - (self.data["pco2_sst"].fdata[i]*self.data["pco2_sst"].fdata[i]) )) + pCO2_salinity_term) )
                    # if pco2_increment + (self.data["sstfndC"].fdata[i] - self.data["pco2_sst"].fdata[i]) + ((self.data["sstfndC"].fdata[i]*self.data["sstfndC"].fdata[i]) - (self.data["pco2_sst"].fdata[i]*self.data["pco2_sst"].fdata[i]) ) >0:#SOCATv4
                    #   print pco2_increment, self.data["sstfndC"].fdata[i] - self.data["pco2_sst"].fdata[i], ((self.data["sstfndC"].fdata[i]*self.data["sstfndC"].fdata[i]) - (self.data["pco2_sst"].fdata[i]*self.data["pco2_sst"].fdata[i]) )#SOCATv4
                else:
                    self.data["pco2_sw_cor"].fdata[i] = self.data["pco2_sw"].fdata[i]
                # following disables temp correction when using SOCAT data
                #if runParams.pco2_data_selection != 0:
                    # check added for testing; only valid for 2010 data, otherwise test is meaningless
                    # e.g. (corrected value:15.115202, uncorrected value:350.351990); exponential should collapse to zero
            
                    # if runParams.year == 2010 and self.data["pco2_sw_cor"].fdata[i] != self.data["pco2_sw"].fdata[i]:
                #         print "%s Using SOCAT data, seawater correction miscalculation, corrected version is not identical to uncorrected version (corrected value:%lf, uncorrected value:%lf, SSTfnd %lf, SSTco2 %lf); exponential should collapse to zero (pco2_increment %lf)" % (function, self.data["pco2_sw_cor"].fdata[i],self.data["pco2_sw"].fdata[i], self.data["sstfndC"].fdata[i], self.data["pco2_sst"].fdata[i], pco2_increment)
                #         self.data["pco2_sw_cor"].fdata[i] = self.data["pco2_sw"].fdata[i]
                
                ###Converts from ppm to microatm THc
                # vco2 in ppm * 1000000 = atm
                # result /1000 to then convert from atm to uatm
                # hence * 0.001 factor
                if runParams.GAS == 'CO2' and runParams.ATMGAS == 'V':
                    #THtodo: 1e-6 can be removed...
                    self.data["pco2_air_cor"].fdata[i] = (self.data["vco2_air"].fdata[i] * 1e-6 * (self.data["pressure"].fdata[i] - self.data["pH2O"].fdata[i]) / (1e-6 * 1013.25)) + (pco2_increment_air)
                else:
                    self.data["pco2_air_cor"].fdata[i] = self.data["vco2_air"].fdata[i]
            
                #runParams.pco2_data_selection ==2 signifies SOCAT fCO2 data, so converting pCO2_air_cor_fdata to fCO2_air_cor_fdata      
                if runParams.pco2_data_selection == 2 or runParams.pco2_data_selection == 4 or runParams.pco2_data_selection == 45:
                    # conversion of pCO2 to fCO2 from McGillis and Wanninkhof 2006, Marine chemistry with correction from Weiss 1974 (as the equation in 2006 paper has a set of brackets missing)
                    b11_fdata[i] = -1636.75 + (12.0408*self.data["sstskin"].fdata[i]) - (0.0327957*self.data["sstskin"].fdata[i]*self.data["sstskin"].fdata[i]) + (3.16528e-5 * self.data["sstskin"].fdata[i]*self.data["sstskin"].fdata[i]*self.data["sstskin"].fdata[i])
                    d12_fdata[i] = 57.7 - (0.118*self.data["sstskin"].fdata[i])
                     # gas constant
                    R = 82.0578 # in [cm^3 atm/(mol K)]
                    # 1/0.987 = 1.0131712 - conversion between bar and atm, so *1013.25 is the conversion from millibar to atm.
                    # the combination of the B11 and d12 terms are in cm^3/mol and so these cancel with the P/RT term (in mol/cm^3) so the whole of the exp term is dimensionless
                    self.data["pco2_air_cor"].fdata[i] = self.data["pco2_air_cor"].fdata[i] * exp((b11_fdata[i] + (2*d12_fdata[i]) ) * 1e-6 * ((self.data["pressure"].fdata[i] * 1013.25)/(R*self.data["sstskin"].fdata[i]) ))
                    #self.data["pco2_air_cor"].fdata[i] = self.data["pco2_air_cor"].fdata[i] * exp((b11_fdata[i] + (2*d12_fdata[i]) ) * ((self.data["pressure"].fdata[i] / 1013.25)/(R*self.data["sstskin"].fdata[i]) ))
                    #print exp((b11_fdata[i] + (2*d12_fdata[i]) ) * 1e-6 * ((self.data["pressure"].fdata[i] * 1013.25)/(R*self.data["sstskin"].fdata[i]) ))
                    #print self.data["pressure"].fdata[i]
                  
                  #((self.data["vco2_air"].fdata[i] * (self.data["pressure"].fdata[i] - self.data["pH2O"].fdata[i]))*0.001) + (pco2_increment)
                  
                  #>>> V=373.769989/1e6
                  #>>> pco2a = 374.220001*1e-6*1013.25
                  
                if runParams.TAKAHASHI_DRIVER:           
                    if self.data["pco2_air"].fdata[i] != missing_value:
                     
                        pCO2a_diff_fdata[i] = self.data["pco2_air_cor"].fdata[i] - self.data["pco2_air"].fdata[i]
                        
                        #dpCO2_diff_fdata[i] = (self.data["pco2_sw_cor"].fdata[i] - self.data["pco2_air_cor"].fdata[i]) - (self.data["pco2_sw_cor"].fdata[i] - self.data["pco2_air"].fdata[i])
                        
                        pH2O_takahashi_fdata[i] = self.data["pressure"].fdata[i] -  (self.data["pco2_air"].fdata[i] *1e-6 * 1013.25) / (self.data["vco2_air"].fdata[i] * 1e-6)
                        #-( 1013.25*((self.data["pco2_air"].fdata[i]*1000000.0)/(self.data["vco2_air"].fdata[i]/1000000.0)) - self.data["pressure"].fdata[i])
                            
                        humidity_fdata[i] = pH2O_takahashi_fdata[i]/self.data["pH2O"].fdata[i]       
                        
                        pH2O_diff_fdata[i] = ((humidity_fdata[i])-1.0) * 100.0
                        
                    else:
                        pH2O_takahashi_fdata[i] = missing_value
                        humidity_fdata[i] = missing_value
                        pH2O_diff_fdata[i] = missing_value
                        # using Takahashi pCO2a to actual flux calculation
                        #self.data["pco2_air_cor"].fdata[i] = self.data["pco2_air"].fdata[i] #RECOMMENTED
                      
            else:
                self.data["pH2O"].fdata[i] = missing_value
                self.data["pco2_air_cor"].fdata[i] = missing_value
                self.data["pco2_sw_cor"].fdata[i] = missing_value
        
#        logger.debug("*** New month ***");
#        logger.debug("salskin: %d", salskin_fc);
#        logger.debug("sstskinK: %d", sstskinK_fc);
#        logger.debug("pres: %d", pres_fc);
#        logger.debug("vco2_air: %d", vco2_air_fc);
#        logger.debug("sstfndK: %d", sstfndK_fc);
#        logger.debug("pco2_sst: %d", sstpco2_fc);
#        logger.debug("pco2_sw: %d", pco2_sw_fc);
#        logger.debug("sstskinK0: %d", sstskinK_fc0);

        ########################################
        # Calculating gas transfer velocity k
        for kParameterisationFunctor in self.kParameterisationFunctors:
            #Check each input exists #TODO: This should go in the pre-run checks!
            for inputDataName in kParameterisationFunctor.input_names():
                if inputDataName not in self.data: #This additional check isn't really needed as it is done in the functor and in the driver script.
                    raise KeyError("Selected kParameterisation (%s) requires input data layers which have not been provided. Required the following DataLayers:\n"+str(kParameterisationFunctor.input_names()));
            
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
              print "\n%s Adding bias to chosen k (k_fdata) parameterisation data (bias value of %lf ms^-1 added, where biology (biology fdata) is > %lf mg m^-3 and wind speed (windu10) is < %lf m s^-1)" % (function, runParams.bias_k_value, runParams.bias_k_biology_value, runParams.bias_k_wind_value)
           else:
              print "\n%s Adding bias to chosen k (k_fdata) parameterisation data (bias value of - %lf percent being used, where biology (biology fdata) is > %lf mg m^-3 and wind speed (windu10) is < %lf m s^-1)" % (function, runParams.bias_k_value, runParams.bias_k_biology_value, runParams.bias_k_wind_value)


         ###############################
         # actual flux calculation 
         ###############################
        
        
         # determine flux based on kHO6_fdata derivation
         # CO2 flux = k * s * Delta_PCO2
         # s = salinity, Delta_PCO2 calculated above
          # solubility conversion between mol kg^-1 atm^-1 to g-C m^-3 uatm^-1
          # mol -> grams for CO2 x 12
          # kg=liter -> m^-3 = x 1000
          # atm -> uatm = /1000000
          # result = x12.0108/1000
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
        conc_factor = (12.0108/1000.0)
        
        self.add_empty_data_layer("FH06");
        self.add_empty_data_layer("concw");
        self.add_empty_data_layer("conca");
        self.add_empty_data_layer("FKo07");

        #If using rain wet deposition, calculate the gas solubility in distilled water
        if runParams.rain_wet_deposition_switch:
            self.add_empty_data_layer("solubility_distilled");
            calculate_solubility_distilled(self.data["solubility_distilled"].fdata, self.data["salinity"].fdata,
                                           runParams.rain_wet_deposition_switch, self.data["sstskin"], DeltaT_fdata, self.nx, self.ny);
        
        if ((runParams.kb_asymmetry != 1.0) and (runParams.k_parameterisation == 3)):
           print "%s kb asymetry has been enabled (runParams.kb_asymmetry:%lf and runParams.k_parameterisation:%d)" % (function, runParams.kb_asymmetry, runParams.k_parameterisation)


        #Main flux calculation loop
        for i in arange(nx * ny):
           self.data["FH06"].fdata[i] = missing_value
          
           if ( (self.data["solubility_fnd"].fdata[i] != missing_value) and (self.data["k"].fdata[i] != missing_value) and (self.data["pco2_sw_cor"].fdata[i] != missing_value) and (self.data["pco2_air_cor"].fdata[i] != missing_value) ):
               # original
              #flux_H06[i] = (k_H06[i] * 36.0 * 24.0 * pow(10.0,-2)) * (solubility_data[i]* (12.0/1000.0)) * DpCO2_cor_data[i]
               # expanded
            
              #mass boundary layer concentration (ie concentration in the water)
              self.data["concw"].fdata[i] = ((self.data["solubility_fnd"].fdata[i] * conc_factor) * self.data["pco2_sw_cor"].fdata[i])
              
              #interfacial concentration (ie at the interface between the ocean and the atmosphere)
              self.data["conca"].fdata[i] = ((self.data["solubility_skin"].fdata[i] * conc_factor) * self.data["pco2_air_cor"].fdata[i])
              
              #flux calculation
              if ((runParams.kb_asymmetry != 1.0) and (runParams.k_parameterisation == 3)):
                 kd_component = self.data["kd"].fdata[i] * k_factor
                 kb_component = self.data["kb"].fdata[i] * k_factor
                 self.data["FH06"].fdata[i] = (kd_component * (self.data["concw"].fdata[i] - self.data["conca"].fdata[i])) + (kb_component * (self.data["concw"].fdata[i] - (runParams.kb_asymmetry *self.data["conca"].fdata[i]) ) )
              else:
                 self.data["FH06"].fdata[i] = (self.data["k"].fdata[i] * k_factor) * (self.data["concw"].fdata[i] - self.data["conca"].fdata[i])
                    
              # using simplified flux calculation with no separation of temperature between airside and waterside CO2
              # assumes that the skin temperature dataset is the only temperature dataset
              if runParams.flux_calc == 3:         
                 self.data["FH06"].fdata[i] = (self.data["k"].fdata[i] * k_factor) * conc_factor * self.data["solubility_skin"].fdata[i] * (self.data["pco2_sw_cor"].fdata[i] - self.data["pco2_air_cor"].fdata[i])#IGA added conc_factor
                 # using Rik W's equation 6 from his new paper
                 #self.data["FH06"].fdata[i] = (7.7e-4*(12.0108/365.0)) * self.data["windu10_moment2"].fdata[i] * (self.data["pco2_sw_cor"].fdata[i] - self.data["pco2_air_cor"].fdata[i])
              
              
              # calculating and adding in flux component for wet deposition due to rain
              if runParams.rain_wet_deposition_switch:
                  # relationship from Komori et al.m 2007
                  # need solubility of CO2 in fresh water
                  # solubility calculated using sstkin
                  # 24.0/1000.0 = conversion from mm hr^-1 to m day^-1
                  # flux is then in g C m^-2 day^-1
                  # flux is always negative, ie going into the ocean
                 self.data["FKo07"].fdata[i] = -(self.data["rain"].fdata[i] * (24.0/1000.0)) * (conc_factor * self.data["solubility_distilled"].fdata[i]) * self.data["pco2_air_cor"].fdata[i]
                 if self.data["FH06"].fdata[i] != missing_value:
                    self.data["FH06"].fdata[i] += self.data["FKo07"].fdata[i]
                 else:
                    self.data["FH06"].fdata[i] = self.data["FKo07"].fdata[i]
              else:
                 self.data["FKo07"].fdata[i] = missing_value
           else:
              self.data["FH06"].fdata[i] = missing_value
              self.data["concw"].fdata[i] = missing_value
              self.data["conca"].fdata[i] = missing_value
              self.data["FKo07"].fdata[i] = missing_value
        
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
        self.add_empty_data_layer("dpco2_cor");
        self.add_empty_data_layer("dpconc_cor");
        for i in arange(self.nx * self.ny):
           if ( (self.data["pco2_sw_cor"].fdata[i] != missing_value) and (self.data["pco2_air_cor"].fdata[i] != missing_value) ):
              self.data["dpco2_cor"].fdata[i] = (self.data["pco2_sw_cor"].fdata[i] - self.data["pco2_air_cor"].fdata[i]) 
           else:
              self.data["dpco2_cor"].fdata[i] = missing_value
        
        for i in arange(self.nx * self.ny):
           if ( (self.data["concw"].fdata[i] != missing_value) and (self.data["conca"].fdata[i] != missing_value) ):
              self.data["dpconc_cor"].fdata[i] = (self.data["concw"].fdata[i] - self.data["conca"].fdata[i]) 
           else:
              self.data["dpconc_cor"].fdata[i] = missing_value

        #Calculate the total number of quality violations per grid location
        self.add_empty_data_layer("failed_quality")
        #self.d = {};
        for outputVar in ["FH06", "kt", "k", "kd", "kb", "salinity", "sstskinC", "concw", "conca", "dpco2_cor", "sstfndC", "pco2_sw_cor"]:
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
              self.data["pco2_sw"].fdata = array([missing_value] * self.nx*self.ny)
              self.data["pco2_sw_stddev"].fdata = array([missing_value] * self.nx*self.ny)
        
        #
        #procesing indictor attribute layers
        #
        #temp if clause using the -l flag for now #TODO: remove this, config should just specify required process indicator layers
        if self.runParams.use_sstfnd_switch == True:
            for piFunctor in self.processIndicatorFunctors:
                try:
                    for inputDataLayer in piFunctor.input_names(): #Check all required inputs exist.
                        if inputDataLayer not in self.data:
                            raise ValueError("%s: Missing input DataLayer (%s) for process indicator functor %s." % (function, inputDataLayer, piFunctor.name));
    
                    for outputDataLayer in piFunctor.output_names(): #Add any output DataLayers that don't already exist
                        if outputDataLayer not in self.data:
                            self.add_empty_data_layer(outputDataLayer);
                    #Execute the process indicator layer functor.
                    piFunctor();
                except ValueError as e:
                    print e.args;
                    print "Exiting...";
                    return 1;
        
        #
        # write out results
        #        
        #Temp refactoring:
        for name in ["krain", "kt", "kb", "kd"]:
            if name not in self.data:
                self.add_empty_data_layer(name);
        
        #write out the final ouput to netcdf
        write_netcdf(self);
        
        
        print "%s SUCCESS writing file %s" % (function, runParams.output_path)
        sys.stdout.flush();
        return 0;

