#! /usr/bin/env python

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
import os, sys
#from Numeric import *
import re # regular expression matching
from math import log, exp, sqrt, pow, isnan
import cerform.flux.coare3
 # numpy import
#from numpy import *
from numpy import size, flipud, mean, zeros, nonzero, array, resize, ma, ravel, arange, float64, dtype, transpose, reshape, ones, meshgrid, median, savetxt, pad
from random import seed, normalvariate

 # debug mode switches
DEBUG = False
DEBUG_PRODUCTS = False
TAKAHASHI_DRIVER = False # enables Takahashi data to drive code - used for verifying calculations
VERIFICATION_RUNS = False # forces flux calculations to use Takahashi SST_t as the SST dataset for the pco2 data
GAS = 'CO2'#IGA when using other gases, this can be used to turn off corrections and changes - vCO2 air is assumed to be the atmospheric concentrations and is not adjusted. No icrements are added as they are with CO2. - Short term solution - long term solution requires complete treatement of other gases
ATMGAS = 'V'#IGA - added for datasets that provide pco2 in air rather than vco2 in air. If vco2, change to ATMGAS = 'V'
print 'Calculations proceeding for ',GAS,' gas and using values in air as ',ATMGAS.lower(),GAS

 # missing value set for intermediate data sets and output dataset
missing_value = -999.0
fill_value = -999.0
missing_value_int = -999
fill_value_int = -999

 # valid data ranges for testing output products and flagging violations
 # orignally set in TS, some have been modified based on the actual input datasets
 # these are used in the NetCDF metadata and to check the valid data ranges to create the
 # OFA11 process indicator layer
OF_min = -0.75
OF_max = 0.6
OK1_min = 0.0
OK1_max = 100.0
OK3_min = 0.0
OK3_max = 100.0
OKD_min = 0.0
OKD_max = 100.0
OKR_min = 0.0
OKR_max = 100.0
OKB1_min = 0.0
OKB1_max = 100.0
OKS1_min = 15.0
OKS1_max = 40.0
OKT1_min = -1.8
OKT1_max = 30.5
OFWR_min = -0.75
OFWR_max = 0.6

 # these thresholds are in g-C m^-3 (rather than ppm as in the TS)
OSC1_min = 0.095
OSC1_max = 0.27
OIC1_min = 0.11
OIC1_max = 0.22

 # these thresholds are in uatm
OBPC_min = 200
OBPC_max = 800


 #
 # method definitions
 #
def GM12_kd_wind(windu10_fdata, windu10_moment2_fdata, windu10_moment3_fdata, scskin_fdata, nx, ny):
 # kd gas transfer from Goddijn-Murphy et al., JGR 2012
   
   kdwind_fdata = array([missing_value] * nx * ny)
   
   for i in arange(nx * ny):   
      kdwind_fdata[i] = missing_value
      if ( (windu10_fdata[i] != missing_value) and (windu10_moment2_fdata[i] != missing_value) and (windu10_moment3_fdata[i] != missing_value) and (scskin_fdata[i] != missing_value) and (scskin_fdata[i] > 0.0) ):
         kdwind_fdata[i] = (2.2 * windu10_fdata[i]) - 3.4
         kdwind_fdata[i] = kdwind_fdata[i] * sqrt(660.0/scskin_fdata[i])              
      else:
         kdwind_fdata[i] = missing_value
     
   return kdwind_fdata

def OceanFluxGHG_kt(kd_fdata, kb_fdata, nx, ny, kb_weighting, kd_weighting):

   # combining the Oceanflux kd and kb components
   ktotal_fdata = array([missing_value] * nx * ny)
   for i in arange(nx * ny):  
       # summing the results
       # units are in 10^-4 m/s
      if ( (kd_fdata[i] != missing_value) and (kb_fdata[i] != missing_value) ):
          # optional weighting in here
         ktotal_fdata[i] = (kb_weighting*kb_fdata[i]) + (kd_weighting*kd_fdata[i])
      elif ( (kd_fdata[i] != missing_value) and (kb_fdata[i] == missing_value) ): 
         ktotal_fdata[i] = kd_fdata[i]
      elif ( (kd_fdata[i] == missing_value) and (kb_fdata[i] != missing_value) ): 
         ktotal_fdata[i] = missing_value
      elif ( (kd_fdata[i] == missing_value) and (kb_fdata[i] == missing_value) ):
         ktotal_fdata[i] = missing_value
      else: # not possible to get into this state, but included for completeness
         ktotal_fdata[i] = missing_value

   return ktotal_fdata

def OceanFluxGHG_k(sigma0_fdata, sig_wv_ht_fdata, windu10_fdata, windu10_moment2_fdata, sstskinC_fdata, nx, ny):
    
    # determine the combined Goddijn-Murphy 2012 and Fangohr and Woolf parameterisation
   kinematic_fdata = array([missing_value] * nx*ny)
   CD_fdata = array([missing_value] * nx*ny)
   friction_fdata = array([missing_value] * nx*ny)
   
   kt_fdata = array([missing_value] * nx*ny)
   kd_fdata = array([missing_value] * nx*ny)
   kb_fdata = array([missing_value] * nx*ny)
   
   for i in arange(nx * ny):   

      # kinematic viscosity
     if ( (sstskinC_fdata[i] != missing_value) ):        
         # kinematic viscosity
        kinematic_fdata[i] = 0.00000183 * exp( (-(sstskinC_fdata[i])) / 36.0)
     else:
        kinematic_fdata[i] = missing_value
     pco2_sw_fdata.shape = (nx, ny)

      # wind drag coefficient
      # algorithm is only value for a wind speed of up to 26 ms^-1
     if ( (windu10_fdata[i] != missing_value) and (windu10_moment2_fdata[i] != missing_value) and (windu10_fdata[i] < 26.0) and (windu10_fdata[i] > 0.0)):
        if (windu10_fdata[i] >= 6.0):
           CD_fdata[i] = 0.60 + (0.070 * windu10_fdata[i])
        elif (windu10_fdata[i] < 6.0):
           CD_fdata[i] = 0.29 + (3.1 / windu10_fdata[i]) + ( 7.7 /(windu10_moment2_fdata[i]))
        else:
           CD_fdata[i] = missing_value
     else:
        CD_fdata[i] = missing_value
     
      # friction velocity in units of ?
     if ( (windu10_fdata[i] != missing_value) and (CD_fdata[i] != missing_value) ):
        friction_fdata[i] = sqrt(CD_fdata[i] * 0.001) * windu10_fdata[i]
     else:
        friction_fdata[i] = missing_value
     
      # kb component     
     if ( (sig_wv_ht_fdata[i] != missing_value) and (friction_fdata[i] != missing_value) and (kinematic_fdata[i] != missing_value) and  (friction_fdata[i] > 0.0) and (sig_wv_ht_fdata[i] > 0.0) ):  
         
        kb_fdata[i] = 0.00002 * ( ((sig_wv_ht_fdata[i] * 100.0) * (friction_fdata[i]) * 100 * 3600) / (kinematic_fdata[i] * 36000000.0) )
        kb_fdata[i] = kb_fdata[i]#/36.0 conversion from cm/h to 10^-4 m/s (100/3600) = 1/36
     else:
        kb_fdata[i] = missing_value
      
      # kd component                 
     if ( (scskin_fdata[i] != missing_value) and (sigma0_fdata[i] != missing_value) and (scskin_fdata[i] > 0.0) ):
        
         # conversion from sigma0 dB units to linear units
        sigma0_fdata[i] = sigma0_fdata[i]/10.0
        sigma0_fdata[i] = pow(10.0,sigma0_fdata[i])
        
        kd_fdata[i] = (( (2100/(sigma0_fdata[i]*sigma0_fdata[i]) ) + 0.1 ) * (pow( (scskin_fdata[i] / 600.0),-0.5))) 
         
        kd_fdata[i] = kd_fdata[i]# /36.0 # unit conversion from cm/h to 10^-4 m/s (100/3600) = 1/36
     else:
        kd_fdata[i] = missing_value
                
   return kd_fdata, kb_fdata

 # writing the final netcdf output
def write_netcdf(FH06_fdata, FKo07_fdata, k_fdata, kt_fdata, kd_fdata, kb_fdata, windu10_fdata, windu10_stddev_fdata, windu10_count_fdata, sig_wv_ht_fdata, sig_wv_ht_stddev_fdata, sig_wv_ht_count_fdata, sstskinC_fdata, sstskinK_stddev_fdata, sstskinK_count_fdata, sstfndC_fdata, sstfndK_stddev_fdata, sstfndK_count_fdata, pco2_air_cor_fdata, pco2_sw_cor_fdata, vco2_air_fdata, pco2_sw_fdata, sal_fdata, pres_fdata, conca_fdata, concw_fdata, rain_fdata, scskin_fdata, krain_fdata, whitecap_fdata, dpco2_cor_fdata, dpconc_cor_fdata,  solskin_fdata, solfnd_fdata, solskinDistilWater_fdata, dpCO2_diff_fdata, pCO2a_diff_fdata, pH2O_fdata, pH20_takahashi_fdata, humidity_fdata, pH20_diff_fdata, solskin_takadata, FH06_takadata, ice_fdata, lowwind_fdata, diurnalw_fdata, bioclass_fdata, bio_fdata, atlantic_ocean_fdata, pacific_ocean_fdata, southern_ocean_fdata, indian_ocean_fdata, longhurst_fdata, susp_particles_fdata, sstgrad_fdata_avg, failed_quality_fdata, k_standard_name, k_long_name, kb_standard_name, kb_long_name, kd_standard_name, kd_long_name, pco2_air_fdata):

    # open a new netCDF file for writing.
    #    need to set format type, defaults to NetCDF4
   ncfile = Dataset(outfile,'w',format='NETCDF3_CLASSIC') 


    # Assign units attributes to coordinate var data. This attaches a
    # text attribute to each of the coordinate variables, containing the
    # units.

    # create the lat and lon dimensions.
   if len(latitude_data.shape)<2:#IGA - If the initial latitude data was a vector, make latitude and longitude as dimensions
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
   secs[:] = time_data
   secs.valid_min = 0.0 
   secs.valid_max = 1.79769313486232e+308

    # Define the coordinate variables. They will hold the coordinate
    # information, that is, the latitudes and longitudes.
   if len(latitude_data.shape)<2:#IGA - If the initial latitude data was a vector, write data as vectors

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
       
       lats2[:] = latitude_data
       lons2[:] = longitude_data

       lats2.valid_min = -90.0
       lats2.valid_max = 90.0

       lons2.valid_min = -180.0
       lons2.valid_max = 180.0
   else:#IGA if the input lat/long was a grid, write output only as a grid.
       lats = ncfile.createVariable('latitude',dtype('float64').char,dims)
       lons = ncfile.createVariable('longitude',dtype('float64').char,dims)
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
       
       lats[:] = latitude_grid
       lons[:] = longitude_grid

       lats.valid_min = -90.0
       lats.valid_max = 90.0

       lons.valid_min = -180.0
       lons.valid_max = 180.0

    #
    # data layers
    #

    # fluxes from chosen k parameterisation
   OF_data = ncfile.createVariable('OF',dtype('float64').char,dims,fill_value=fill_value)
   OF_data[:] = FH06_fdata
   OF_data.units = 'g C m-2 day-1'
   OF_data.missing_value = missing_value
   OF_data.valid_min = -1000
   OF_data.valid_max = +1000
   #OF_data.scale_factor = 1.0
   #OF_data.add_offset = 0.0
   OF_data.standard_name = "air_to_sea_carbon_dioxide_flux"
   OF_data.long_name = "Air-sea CO2 flux using the %s %s gas transfer velocity (k)" % (k_standard_name, k_long_name)

   OW1_data = ncfile.createVariable('OW1',dtype('float64').char,dims,fill_value=fill_value)
   OW1_data[:] = whitecap_fdata
   OW1_data.units = '%'
   OW1_data.missing_value = missing_value
   OW1_data.valid_min = 0.0
   OW1_data.valid_max = 100
   #OW1_data.scale_factor = 1.0
   #OW1_data.add_offset = 0.0
   OW1_data.standard_name = "percentage_whitecap"
   OW1_data.long_name = "Percentage whitecapping from simple U10 wind parameterisation"

    # FKo07 flux due to wet deposition from rain
   if rain_wet_deposition:
      OFWR_data = ncfile.createVariable('OFWR',dtype('float64').char,dims,fill_value=fill_value)
      OFWR_data[:] = FKo07_fdata
      OFWR_data.units = 'g C m-2 day-1'
      OFWR_data.missing_value = missing_value
      OFWR_data.valid_min = OFWR_min
      OFWR_data.valid_max = OFWR_max
      #OFWR_data.scale_factor = 1.0
      #OFWR_data.add_offset = 0.0
      OFWR_data.standard_name = "air_to_sea_carbon_dioxide_flux_due_rain"
      OFWR_data.long_name = "Air-sea CO2 flux using the Komori et al., 2007 wet deposition term"


    # create the oceanflux GHG kt dataset
   OK1_data = ncfile.createVariable('OK1',dtype('float64').char,dims,fill_value=fill_value)
   OK1_data[:] = kt_fdata
   OK1_data.units = 'cm h-1'
   OK1_data.missing_value = missing_value
   OK1_data.valid_min = OK1_min
   OK1_data.valid_max = OK1_max
   #OK1_data.scale_factor = 1.0
   #OK1_data.add_offset = 0.0
   OK1_data.standard_name = kt_standard_name
   OK1_data.long_name = "OceanFluxGHG total (kd+kb) gas transfer velocity"

    # create the chosen k dataset
   OK3_data = ncfile.createVariable('OK3',dtype('float64').char,dims,fill_value=fill_value)
   OK3_data[:] = k_fdata
   OK3_data.units = 'cm h-1'
   OK3_data.missing_value = missing_value
   OK3_data.valid_min = OK3_min
   OK3_data.valid_max = OK3_max
   #OK3_data.scale_factor = 1.0
   #OK3_data.add_offset = 0.0
   OK3_data.standard_name = k_standard_name
   OK3_data.long_name = "Chosen Gas transfer velocity (%s)" % (k_long_name)

    # create the kb dataset
   OKB1_data = ncfile.createVariable('OKB1',dtype('float64').char,dims,fill_value=fill_value)
   OKB1_data[:] = kb_fdata
   OKB1_data.units = 'cm h-1'
   OKB1_data.missing_value = missing_value
   OKB1_data.valid_min = OKB1_min
   OKB1_data.valid_max = OKB1_max
   #OKB1_data.scale_factor = 1.0
   #OKB1_data.add_offset = 0.0
   OKB1_data.standard_name = kb_standard_name
   OKB1_data.long_name = "Bubble mediated gas transfer velocity (%s)" % (kb_long_name)

    # create the kd dataset
   OKD_data = ncfile.createVariable('OKD',dtype('float64').char,dims,fill_value=fill_value)
   OKD_data[:] = kd_fdata
   OKD_data.units = 'cm h-1'
   OKD_data.missing_value = missing_value
   OKD_data.valid_min = OKD_min
   OKD_data.valid_max = OKD_max
   #OKD_data.scale_factor = 1.0
   #OKD_data.add_offset = 0.0
   OKD_data.standard_name = kd_standard_name
   OKD_data.long_name = "Direct gas transfer velocity (%s)" % (kd_long_name)

    # create the shosen k dataset
   OKR_data = ncfile.createVariable('OKR',dtype('float64').char,dims,fill_value=fill_value)
   OKR_data[:] = krain_fdata
   OKR_data.units = 'cm h-1'
   OKR_data.missing_value = missing_value
   OKR_data.valid_min = OKR_min
   OKR_data.valid_max = OKR_max
   #OKR_data.scale_factor = 1.0
   #OKR_data.add_offset = 0.0
   OKR_data.standard_name = "gas_transfer_at_sea_surface_due_rain"
   OKR_data.long_name = "Gas transfer velocity due to rain from Ho et al, Tellus, 1997"

    # writing out the wind data
   WS1_mean_data = ncfile.createVariable('WS1_mean',dtype('float64').char,dims,fill_value=missing_value)
   WS1_mean_data[:] = windu10_fdata
   WS1_mean_data.units = 'm s-1'
   WS1_mean_data.missing_value = missing_value
   WS1_mean_data.valid_min = 0.0
   WS1_mean_data.valid_max = 100.0
   #WS1_mean_data.scale_factor = 1.0
   #WS1_mean_data.add_offset = 0.0
   WS1_mean_data.standard_name = "wind_speed"
   WS1_mean_data.long_name = "Wind speed at 10m used for determining k"

   WS1_stddev_data = ncfile.createVariable('WS1_stddev',dtype('float64').char,dims,fill_value=missing_value)
   WS1_stddev_data[:] = windu10_stddev_fdata
   WS1_stddev_data.units = 'm s-1'
   WS1_stddev_data.missing_value = missing_value
   WS1_stddev_data.valid_min = 0.0
   WS1_stddev_data.valid_max = 100.0
   #WS1_stddev_data.scale_factor = 1.0
   #WS1_stddev_data.add_offset = 0.0
   WS1_stddev_data.standard_name = "wind_speed_stddev"
   WS1_stddev_data.long_name = "Standard deviation of wind speed at 10m used for calculating the mean wind speed"

   WS1_N_data = ncfile.createVariable('WS1_N',dtype('float64').char,dims,fill_value=missing_value)
   WS1_N_data[:] = windu10_count_fdata
   WS1_N_data.units = 'm s-1'
   WS1_N_data.missing_value = missing_value
   WS1_N_data.valid_min = 0.0
   WS1_N_data.valid_max = 100.0
   #WS1_N_data.scale_factor = 1.0
   #WS1_N_data.add_offset = 0.0
   WS1_N_data.standard_name = "wind_speed_count"
   WS1_N_data.long_name = "Count of data points of wind speed at 10m used for calculating the mean wind speed"

    # writing out the significant wave height data
   SW1_mean_data = ncfile.createVariable('SW1_mean',dtype('float64').char,dims,fill_value=missing_value)
   SW1_mean_data[:] = sig_wv_ht_fdata
   SW1_mean_data.units = 'm'
   SW1_mean_data.missing_value = missing_value
   SW1_mean_data.valid_min = 0.0
   SW1_mean_data.valid_max = 100.0
   #SW1_mean_data.scale_factor = 1.0
   #SW1_mean_data.add_offset = 0.0
   SW1_mean_data.standard_name = "sea_surface_wave_significant_height"
   SW1_mean_data.long_name = "Significant wave height data used for determining kb"

    # writing out the significant wave height data
   SW1_stddev_data = ncfile.createVariable('SW1_stddev',dtype('float64').char,dims,fill_value=missing_value)
   SW1_stddev_data[:] = sig_wv_ht_stddev_fdata
   SW1_stddev_data.units = 'm'
   SW1_stddev_data.missing_value = missing_value
   SW1_stddev_data.valid_min = 0.0
   SW1_stddev_data.valid_max = 100.0
   #SW1_stddev_data.scale_factor = 1.0
   #SW1_stddev_data.add_offset = 0.0
   SW1_stddev_data.standard_name = "sea_surface_wave_significant_height_stddev"
   SW1_stddev_data.long_name = "Standard deviation of significant wave height data used to calculate mean significant wave height"

    # writing out the significant wave height data
   SW1_N_data = ncfile.createVariable('SW1_N',dtype('float64').char,dims,fill_value=missing_value)
   SW1_N_data[:] = sig_wv_ht_count_fdata
   SW1_N_data.units = 'm'
   SW1_N_data.missing_value = missing_value
   SW1_N_data.valid_min = 0.0
   SW1_N_data.valid_max = 100.0
   #SW1_N_data.scale_factor = 1.0
   #SW1_N_data.add_offset = 0.0
   SW1_N_data.standard_name = "sea_surface_wave_significant_height_count"
   SW1_N_data.long_name = "Count of data points of significant wave height data used to calculate mean significant wave height"


   # writing out the sst_skin data
   ST1_mean_data = ncfile.createVariable('ST1_mean',dtype('float64').char,dims,fill_value=fill_value)
   ST1_mean_data[:] = sstskinC_fdata
   ST1_mean_data.units = 'degreeC'
   ST1_mean_data.missing_value = missing_value
   ST1_mean_data.valid_min = -1.8
   ST1_mean_data.valid_max = 100.0
   #ST1_mean_data.scale_factor = 1.0
   #ST1_mean_data.add_offset = 0.0
   ST1_mean_data.standard_name = "sea_surface_skin_temperature"
   ST1_mean_data.long_name = "Sea surface skin temperature used for the flux calculations"

   ST1_stddev_data = ncfile.createVariable('ST1_stddev',dtype('float64').char,dims,fill_value=fill_value)
   ST1_stddev_data[:] = sstskinK_stddev_fdata
   ST1_stddev_data.units = 'degreeC'
   ST1_stddev_data.missing_value = missing_value
   ST1_stddev_data.valid_min = 0.0
   ST1_stddev_data.valid_max = 100.0
   #ST1_stddev_data.scale_factor = 1.0
   #ST1_stddev_data.add_offset = 0.0
   ST1_stddev_data.standard_name = "sea_surface_skin_temperature_stddev"
   ST1_stddev_data.long_name = "Standard deviation of Sea surface skin temperature used to create the mean sea surface skin data"

   ST1_N_data = ncfile.createVariable('ST1_N',dtype('float64').char,dims,fill_value=fill_value)
   ST1_N_data[:] = sstskinK_count_fdata
   ST1_N_data.units = 'degreeC'
   ST1_N_data.missing_value = missing_value
   ST1_N_data.valid_min = 0.0
   ST1_N_data.valid_max = 100.0
   #ST1_N_data.scale_factor = 1.0
   #ST1_N_data.add_offset = 0.0
   ST1_N_data.standard_name = "sea_surface_skin_temperature_count"
   ST1_N_data.long_name = "Count of data points of Sea surface skin temperature used to create the mean sea surface skin data"

   FT1_mean_data = ncfile.createVariable('FT1_mean',dtype('float64').char,dims,fill_value=fill_value)
   FT1_mean_data[:] = sstfndC_fdata
   FT1_mean_data.units ='degreeC'
   FT1_mean_data.missing_value = missing_value
   FT1_mean_data.valid_min = -1.8
   FT1_mean_data.valid_max = 100.0
   #FT1_mean_data.scale_factor = 1.0
   #FT1_mean_data.add_offset = 0.0
   FT1_mean_data.standard_name = "sea_surface_foundation_temperature"
   FT1_mean_data.long_name = "foundation SST used for the flux calculations"

   FT1_stddev_data = ncfile.createVariable('FT1_stddev',dtype('float64').char,dims,fill_value=fill_value)
   FT1_stddev_data[:] = sstfndK_stddev_fdata
   FT1_stddev_data.units ='degreeC'
   FT1_stddev_data.missing_value = missing_value
   FT1_stddev_data.valid_min = 0.0
   FT1_stddev_data.valid_max = 100.0
   #FT1_stddev_data.scale_factor = 1.0
   #FT1_stddev_data.add_offset = 0.0
   FT1_stddev_data.standard_name = "sea_surface_foundation_temperature_stddev"
   FT1_stddev_data.long_name = "standard deviation of data used to calculate the mean foundation SST"

   FT1_N_data = ncfile.createVariable('FT1_N',dtype('float64').char,dims,fill_value=fill_value)
   FT1_N_data[:] = sstfndK_count_fdata
   FT1_N_data.units ='degreeC'
   FT1_N_data.missing_value = missing_value
   FT1_N_data.valid_min = 0.0
   FT1_N_data.valid_max = 100.0
   #FT1_N_data.scale_factor = 1.0
   #FT1_N_data.add_offset = 0.0
   FT1_N_data.standard_name = "sea_surface_foundation_temperature_stddev"
   FT1_N_data.long_name = "Count of data points used to create the mean foundation sea surface temperature data"

   OAPC1_data = ncfile.createVariable('OAPC1',dtype('float64').char,dims,fill_value=fill_value)
   OAPC1_data[:] = pco2_air_cor_fdata
   OAPC1_data.units ='microatm'
   OAPC1_data.missing_value = missing_value
   OAPC1_data.valid_min = 0.0
   OAPC1_data.valid_max = 500.0
   #OAPC1_data.scale_factor = 1.0
   #OAPC1_data.add_offset = 0.0
   OAPC1_data.standard_name = "water_surface_partial_pressure_of_carbon_dioxide_in_air"
   OAPC1_data.long_name = "Water surface partial pressure (pCO2) in air from climtology corrected using modelled sea level pressure"

   OBPC_data = ncfile.createVariable('OBPC',dtype('float64').char,dims,fill_value=fill_value)
   OBPC_data[:] = pco2_sw_cor_fdata
   OBPC_data.units ='microatm'
   OBPC_data.missing_value = missing_value
   OBPC_data.valid_min = 0.0
   OBPC_data.valid_max = 1000.0
   #OBPC_data.scale_factor = 1.0
   #OBPC_data.add_offset = 0.0
   OBPC_data.standard_name = "sub_skin_partial_pressure_of_carbon_dioxide_in_sea_water"
   OBPC_data.long_name = "Sub skin partial pressure (pCO2) of carbon dioxide"

   SFUG_data = ncfile.createVariable('SFUG_krig',dtype('float64').char,dims,fill_value=fill_value)
   SFUG_data[:] = pco2_sw_fdata
   SFUG_data.units ='microatm'
   SFUG_data.missing_value = missing_value
   SFUG_data.valid_min = 0.0
   SFUG_data.valid_max = 1000.0
   #SFUG_data.scale_factor = 1.0
   #SFUG_data.add_offset = 0.0
   SFUG_data.standard_name = "surface_partial_pressure_of_carbon_dioxide_in_sea_water"
   SFUG_data.long_name = "pCO2 in sea water from climtology"

   SFUG_stddev_data = ncfile.createVariable('SFUG_stddev',dtype('float64').char,dims,fill_value=fill_value)
   SFUG_stddev_data[:] = pco2_sw_stddev_fdata
   SFUG_stddev_data.units ='microatm'
   SFUG_stddev_data.missing_value = missing_value
   SFUG_stddev_data.valid_min = 0.0
   SFUG_stddev_data.valid_max = 1000.0
   #SFUG_stddev_data.scale_factor = 1.0
   #SFUG_stddev_data.add_offset = 0.0
   SFUG_stddev_data.standard_name = "surface_partial_pressure_of_carbon_dioxide_in_sea_water"
   SFUG_stddev_data.long_name = "pCO2 in sea water from climtology"

   OKT1_data = ncfile.createVariable('OKT1',dtype('float64').char,dims,fill_value=fill_value)
   OKT1_data[:] = sstskinC_fdata
   OKT1_data.units ='degreeC' # count means dimensionless
   OKT1_data.missing_value = missing_value
   OKT1_data.valid_min = OKT1_min
   OKT1_data.valid_max = OKT1_max
   #OKT1_data.scale_factor = 1.0
   #OKT1_data.add_offset = 0.0
   OKT1_data.standard_name = "sea_surface_temperature_used_for_gas_transfer_velocity_and_schmidt"
   OKT1_data.long_name = "sea surface temperature used to calculate the gas transfer velocity and schmidt number"

   OKS1_data = ncfile.createVariable('OKS1',dtype('float64').char,dims,fill_value=fill_value)
   OKS1_data[:] = sal_fdata
   OKS1_data.units ='count' # count means dimensionless
   OKS1_data.missing_value = missing_value
   OKS1_data.valid_min = OKS1_min
   OKS1_data.valid_max = OKS1_max
   #OKS1_data.scale_factor = 1.0
   #OKS1_data.add_offset = 0.0
   OKS1_data.standard_name = "sea_water_salinity"
   OKS1_data.long_name = "salinity in sea water"

   R1_data = ncfile.createVariable('R1',dtype('float64').char,dims,fill_value=fill_value)
   R1_data[:] = rain_fdata
   R1_data.units ='mm hr-1' 
   R1_data.missing_value = missing_value
   R1_data.valid_min = 0.0
   R1_data.valid_max = 1000.0
   #R1_data.scale_factor = 1.0
   #R1_data.add_offset = 0.0
   R1_data.standard_name = "rainfall_rate"
   R1_data.long_name = "Precipitation Estimate from EO and/or satellite/gauge combined data set"

   SC_data = ncfile.createVariable('SC',dtype('float64').char,dims,fill_value=fill_value)
   SC_data[:] = scskin_fdata
   SC_data.units ='count' # count means dimensionless
   SC_data.missing_value = missing_value
   SC_data.valid_min = 0.0
   SC_data.valid_max = 4000.0
   #SC_data.scale_factor = 1.0
   #SC_data.add_offset = 0.0
   SC_data.standard_name = "schmidt_number_at_sea_skin"
   SC_data.long_name = "Schmidt number calculated using the sea surface skin temperature"

   OIC1_data = ncfile.createVariable('OIC1',dtype('float64').char,dims,fill_value=fill_value)
   OIC1_data[:] = conca_fdata
   OIC1_data.units ='g m-3'
   OIC1_data.missing_value = missing_value
   OIC1_data.valid_min = 0.0
   OIC1_data.valid_max = 400.0
   #OIC1_data.scale_factor = 1.0
   #OIC1_data.add_offset = 0.0
   OIC1_data.standard_name = "concentration_of_carbon_dioxide_at_the_sea_air_interface"
   OIC1_data.long_name = "Concentration of carbon dioxide at the sea water and air interface in g-C m^-3"

   OSFC_data = ncfile.createVariable('OSFC',dtype('float64').char,dims,fill_value=fill_value)
   OSFC_data[:] = concw_fdata
   OSFC_data.units ='g m-3'
   OSFC_data.missing_value = missing_value
   OSFC_data.valid_min = 0.0
   OSFC_data.valid_max = 500.0
   #OSFC_data.scale_factor = 1.0
   #OSFC_data.add_offset = 0.0
   OSFC_data.standard_name = "sub_skin_concentration_of_carbon_dioxide_in_sea_water"
   OSFC_data.long_name = "Sub skin concentration of carbon dioxide in sea water g-C m^-3"

   P1_data = ncfile.createVariable('P1',dtype('float64').char,dims,fill_value=fill_value)
   P1_data[:] = ice_fdata
   P1_data.units ='%'
   P1_data.missing_value = missing_value
   P1_data.valid_min = 0.0
   P1_data.valid_max = 100.0
   #P1_data.scale_factor = 1.0
   #P1_data.add_offset = 0.0
   P1_data.standard_name = "sea_ice_area_fraction"
   P1_data.long_name = "Percentage of sea area covered by ice (fraction if Takahashi driven)"


   if TAKAHASHI_DRIVER:
      taka_data = ncfile.createVariable('humidity',dtype('float64').char,dims,fill_value=fill_value)
      taka_data[:] = humidity_fdata
      taka_data.units ='%'
      taka_data.missing_value = missing_value
      taka_data.valid_min = -1500.0
      taka_data.valid_max = 1500.0
      taka_data.standard_name = "relative_humidity"
      taka_data.long_name = "relative humidity calculated using ratio of PH2O climatology and PH20 sst/salinity relationship"

      taka_data = ncfile.createVariable('pH2O_diff',dtype('float64').char,dims,fill_value=fill_value)
      taka_data[:] = pH2O_diff_fdata
      taka_data.units ='%'
      taka_data.missing_value = missing_value
      taka_data.valid_min = -1500.0
      taka_data.valid_max = 1500.0
      taka_data.standard_name = "water_vapour_pressure_difference"
      taka_data.long_name = "pH2O_takahashi minus PH2O calculated"

      taka_data = ncfile.createVariable('pH2O_takahashi',dtype('float64').char,dims,fill_value=fill_value)
      taka_data[:] = pH2O_takahashi_fdata
      taka_data.units ='mb'
      taka_data.missing_value = missing_value
      taka_data.valid_min = -1500.0
      taka_data.valid_max = 1500.0
      taka_data.standard_name = "water_vapour_pressure"
      taka_data.long_name = "pH2O_takahashi from climatology"
      
      taka_data = ncfile.createVariable('pCO2a_diff',dtype('float64').char,dims,fill_value=fill_value)
      taka_data[:] = pCO2a_diff_fdata
      taka_data.units ='microatm'
      taka_data.missing_value = missing_value
      taka_data.valid_min = -1500.0
      taka_data.valid_max = 1500.0
      taka_data.standard_name = "difference_pco2a"
      taka_data.long_name = "pco2a calculated from Wiess and Price 1980 minus Takahashi pCO2a"
      
      taka_data = ncfile.createVariable('dpCO2_diff',dtype('float64').char,dims,fill_value=fill_value)
      taka_data[:] = dpCO2_diff_fdata
      taka_data.units ='microatm'
      taka_data.missing_value = missing_value
      taka_data.valid_min = -1500.0
      taka_data.valid_max = 1500.0
      taka_data.standard_name = "difference_dpco2"
      taka_data.long_name = "dpco2 calculated minus Takahashi dpCO2"

      taka_data = ncfile.createVariable('taka_skin_solubility',dtype('float64').char,('time','latitude','longitude'),fill_value=fill_value)
      taka_data[:] = solskin_takadata
      taka_data.units ='mmol kg-1 atm-1'
      taka_data.missing_value = missing_value
      taka_data.valid_min = 0
      taka_data.valid_max = 1
      taka_data.standard_name = "solubility_of_carbon_dioxide_in_seawater_at_the_sea_air_interface"
      taka_data.long_name = "solubility as calculated using the sstskin (ST1) data"

      taka_data = ncfile.createVariable('taka_flux',dtype('float64').char,('time','latitude','longitude'),fill_value=fill_value)
      taka_data[:] = FH06_takadata
      taka_data.units ='g C m-2 mon-1'
      taka_data.missing_value = missing_value
      taka_data.valid_min = 0
      taka_data.valid_max = 1
      taka_data.standard_name = "air_to_sea_carbon_dioxide_flux"
      taka_data.long_name = "Air-sea CO2 flux using the %s %s gas transfer velocity (k)" % (k_standard_name, k_long_name)
      

   if DEBUG_PRODUCTS:
      PH2O_data = ncfile.createVariable('PH2O',dtype('float64').char,dims,fill_value=fill_value)
      PH2O_data[:] = pH2O_fdata
      PH2O_data.units ='mb'
      PH2O_data.missing_value = missing_value
      PH2O_data.valid_min = 0.0
      PH2O_data.valid_max = 1500.0
      PH2O_data.standard_name = "water_vapour_pressure"
      PH2O_data.long_name = "water vapour pressure from salinity and SST relationship of Weiss and Price 1980, Marine Chemistry"      
      
      pres_data = ncfile.createVariable('air_pressure',dtype('float64').char,dims,fill_value=fill_value)
      pres_data[:] = pres_fdata
      pres_data.units ='millibar'
      pres_data.missing_value = missing_value
      pres_data.valid_min = 0.0
      pres_data.valid_max = 1500.0
      #pres_data.scale_factor = 1.0
      #pres_data.add_offset = 0.0
      pres_data.standard_name = "air_pressure"
      pres_data.long_name = "Daily mean sea level air pressure derived from modelled air pressure data"

      vco2_air_data = ncfile.createVariable('VCO2',dtype('float64').char,dims,fill_value=fill_value)
      vco2_air_data[:] = vco2_air_fdata
      vco2_air_data.units ='micromol mol-1'
      vco2_air_data.missing_value = missing_value
      vco2_air_data.valid_min = 0.0
      vco2_air_data.valid_max = 500.0
      #vco2_air_data.scale_factor = 1.0
      #vco2_air_data.add_offset = 0.0
      vco2_air_data.standard_name = "dry_molecular_fraction_of_carbon_dioxide_in_atmosphere"
      vco2_air_data.long_name = "Concentration of CO2 in dry air in 2010 from NOAA ESLR (or dry molecular fraction of carbon dioxide in the atmosphere)"

      dpco2_cor_data = ncfile.createVariable('dpCO2',dtype('float64').char,dims,fill_value=fill_value)
      dpco2_cor_data[:] = dpco2_cor_fdata
      dpco2_cor_data.units ='microatm'
      dpco2_cor_data.missing_value = missing_value
      dpco2_cor_data.valid_min = -500.0
      dpco2_cor_data.valid_max = 500.0
      #dpco2_cor_data.scale_factor = 1.0
      #dpco2_cor_data.add_offset = 0.0
      dpco2_cor_data.standard_name = "difference_in_water_and_air_partial_pressure_of_carbon_dioxide"
      dpco2_cor_data.long_name = "Delta pCO2"

      dpconc_cor_data = ncfile.createVariable('Dconc',dtype('float64').char,dims,fill_value=fill_value)
      dpconc_cor_data[:] = dpconc_cor_fdata
      dpconc_cor_data.units ='g m-3'
      dpconc_cor_data.missing_value = missing_value
      dpconc_cor_data.valid_min = -500.0
      dpconc_cor_data.valid_max = 500.0
      #dpconc_cor_data.scale_factor = 1.0
      #dpconc_cor_data.add_offset = 0.0
      dpconc_cor_data.standard_name = "difference_in_mass_boundary_layer_and_interface_concentrations_carbon_dioxide"
      dpconc_cor_data.long_name = "difference between the interface and mass boundary layer concentrations"
      
      solskin_data = ncfile.createVariable('skin_solubility',dtype('float64').char,dims,fill_value=fill_value)
      solskin_data[:] = solskin_fdata
      solskin_data.units ='mol kg-1 atm-1'
      solskin_data.missing_value = missing_value
      solskin_data.valid_min = 0.0
      solskin_data.valid_max = 1.0
      #solskin_data.scale_factor = 1.0
      #solskin_data.add_offset = 0.0
      solskin_data.standard_name = "solubility_of_carbon_dioxide_in_seawater_at_the_sea_air_interface"
      solskin_data.long_name = "solubility as calculated using the sstskin (ST1) data"
      
      solfnd_data = ncfile.createVariable('fnd_solubility',dtype('float64').char,dims,fill_value=fill_value)
      solfnd_data[:] = solfnd_fdata
      solfnd_data.units ='mol kg-1 atm-1'
      solfnd_data.missing_value = missing_value
      solfnd_data.valid_min = 0.0
      solfnd_data.valid_max = 1.0
      #solfnd_data.scale_factor = 1.0
      #solfnd_data.add_offset = 0.0
      solfnd_data.standard_name = "solubility_of_carbon_dioxide_in_seawater_at_depth"
      solfnd_data.long_name = "solubility as calculated using the sstfnd (FT1) data"

      solskinDistil_data = ncfile.createVariable('solskinDistilWater',dtype('float64').char,dims,fill_value=fill_value)
      solskinDistil_data[:] = solskinDistilWater_fdata
      solskinDistil_data.units ='mol kg-1 atm-1'
      solskinDistil_data.missing_value = missing_value
      solskinDistil_data.valid_min = 0.0
      solskinDistil_data.valid_max = 1.0
      solskinDistil_data.standard_name = "solubility_of_carbon_dioxide_in_distilled_water"
      solskinDistil_data.long_name = "solubility as calculated using the sstskin (ST1) data in distilled water"
      
      pco2_air_d = ncfile.createVariable('pco2_air',dtype('float64').char,dims,fill_value=fill_value)
      pco2_air_d[:] = pco2_air_fdata
      pco2_air_d.units ='microatm?'
      pco2_air_d.missing_value = missing_value
      pco2_air_d.valid_min = 0.0
      pco2_air_d.valid_max = 500
      pco2_air_d.standard_name = "partial pressure of CO2 in air"
      pco2_air_d.long_name = "partial pressure of CO2 in air"
    #
    # process indicator attribute layers
    #
   if process_layers_off == 0:
       # OFA1 - low wind
      OFA1_data = ncfile.createVariable('OFA01',dtype('int32').char,dims,fill_value=fill_value_int)
      OFA1_data[:] = lowwind_fdata
      OFA1_data.units ='count'
      OFA1_data.missing_value = missing_value_int
      OFA1_data.valid_min = 0
      OFA1_data.valid_max = 1
      OFA1_data.scale_factor = 1
      OFA1_data.add_offset = 0
      #OFA1_data.standard_name = "Regions of low wind"
      OFA1_data.long_name = "Regions of low wind determined from the U10 input data"

       # OFA3 diurmal warming
      OFA3_data = ncfile.createVariable('OFA03',dtype('int32').char,dims,fill_value=fill_value_int)
      OFA3_data[:] = diurnalw_fdata
      OFA3_data.units ='count'
      OFA3_data.missing_value = missing_value_int
      OFA3_data.valid_min = 0
      OFA3_data.valid_max = 1
      OFA3_data.scale_factor = 1
      OFA3_data.add_offset = 0
      #OFA3_data.standard_name = "Regions of diurnal warming"
      OFA3_data.long_name = "Regions of diurnal warming where sstskin > sstfnd (as determined from the sstskin and sstfnd input data)"

       # OFA4 high bilogical activity
      OFA4_data = ncfile.createVariable('OFA04',dtype('int32').char,dims,fill_value=fill_value_int)
      OFA4_data[:] = bioclass_fdata
      OFA4_data.units ='count'
      OFA4_data.missing_value = missing_value_int
      OFA4_data.valid_min = 0
      OFA4_data.valid_max = 3
      OFA4_data.scale_factor = 1
      OFA4_data.add_offset = 0
      #OFA4_data.standard_name = "mass_concentration_of_chlorophyll_a_in_sea_water"
      OFA4_data.long_name = "Regions of low, medium and high biological activity"

       # OFA5 sst gradient layer
      OFA5_data = ncfile.createVariable('OFA05',dtype('float64').char,dims,fill_value=fill_value)
      OFA5_data[:] = sstgrad_fdata_avg
      OFA5_data.units ='count'
      OFA5_data.missing_value = missing_value
      OFA5_data.valid_min = 0.0
      OFA5_data.valid_max = 100.0
      #OFA5_data.scale_factor = 1.0
      #OFA5_data.add_offset = 0.0
      OFA5_data.standard_name = "sea_surface_temperature_gradients"
      OFA5_data.long_name = "climatology of regions of of strong sst gradients as determined from GHRSST OSTIA data"
   
       # OFA7 longhurst masks
      OFA6_data = ncfile.createVariable('OFA06',dtype('int32').char,dims,fill_value=fill_value_int)
      OFA6_data[:] = longhurst_fdata
      OFA6_data.units ='count'
      OFA6_data.missing_value = missing_value_int
      OFA6_data.valid_min = 0
      OFA6_data.valid_max = 51 # 51 provinces in total
      OFA6_data.scale_factor = 1
      OFA6_data.add_offset = 0
      #OFA6_data.standard_name = "provinces"
      OFA6_data.long_name = "Longhurst provinces"

       # OFA7 oceans masks
      OFA7_data = ncfile.createVariable('OFA07_atlantic',dtype('float64').char,dims,fill_value=fill_value)
      OFA7_data[:] = atlantic_ocean_fdata
      OFA7_data.units ='count'
      OFA7_data.missing_value = missing_value
      OFA7_data.valid_min = 0.0
      OFA7_data.valid_max = 1.0
      #OFA7_data.scale_factor = 1.0
      #OFA7_data.add_offset = 0.0
      #OFA7_data.standard_name = "oceanic basins"
      OFA7_data.long_name = "Atlantic Ocean mask"

       # OFA7 oceans masks
      OFA7_data = ncfile.createVariable('OFA07_pacific',dtype('float64').char,dims,fill_value=fill_value)
      OFA7_data[:] = pacific_ocean_fdata
      OFA7_data.units ='count'
      OFA7_data.missing_value = missing_value
      OFA7_data.valid_min = 0.0
      OFA7_data.valid_max = 1.0
      #OFA7_data.scale_factor = 1.0
      #OFA7_data.add_offset = 0.0
      #OFA7_data.standard_name = "oceanic basins"
      OFA7_data.long_name = "Pacific Ocean mask"

       # OFA7 oceans masks
      OFA7_data = ncfile.createVariable('OFA07_southern',dtype('float64').char,dims,fill_value=fill_value)
      OFA7_data[:] = southern_ocean_fdata
      OFA7_data.units ='count'
      OFA7_data.missing_value = missing_value
      OFA7_data.valid_min = 0.0
      OFA7_data.valid_max = 1.0
      #OFA7_data.scale_factor = 1.0
      #OFA7_data.add_offset = 0.0
      #OFA7_data.standard_name = "oceanic basins"
      OFA7_data.long_name = "Southern Ocean mask"

       # OFA7 oceans masks
      OFA7_data = ncfile.createVariable('OFA07_indian',dtype('float64').char,dims,fill_value=fill_value)
      OFA7_data[:] = indian_ocean_fdata
      OFA7_data.units ='count'
      OFA7_data.missing_value = missing_value
      OFA7_data.valid_min = 0.0
      OFA7_data.valid_max = 1.0
      #OFA7_data.scale_factor = 1.0
      #OFA7_data.add_offset = 0.0
      #OFA7_data.standard_name = "oceanic basins"
      OFA7_data.long_name = "Indian Ocean mask"

      # OFA10 remote sensing reflectance
      OFA10_data = ncfile.createVariable('OFA10',dtype('float64').char,dims,fill_value=fill_value)
      OFA10_data[:] = susp_particles_fdata
      OFA10_data.units ='sr-1'
      OFA10_data.missing_value = missing_value
      OFA10_data.valid_min = 0.0
      OFA10_data.valid_max = 1.0
      #OFA10_data.scale_factor = 1.0
      #OFA10_data.add_offset = 0.0
      OFA10_data.standard_name = "surface_ratio_of_upwelling_radiance_emerging_from_sea_water_to_downwelling_radiative_flux_in_air"
      OFA10_data.long_name = "Ratio of water leaving radiance to downwelling irradiance"

      # OFA12 chlorophyll-a
      OFA12_data = ncfile.createVariable('OFA12',dtype('float64').char,dims,fill_value=fill_value)
      OFA12_data[:] = bio_fdata
      OFA12_data.units ='mg m-3'
      OFA12_data.missing_value = missing_value
      OFA12_data.valid_min = 0
      OFA12_data.valid_max = 100
      OFA12_data.standard_name = "mass_concentration_of_chlorophyll_a_in_sea_water"
      OFA12_data.long_name = "mass concentration of chlorphyll-a at the surface from ESA Ocean Colour CCI data"


    # this process indicator layer is always written out as its part of the internal quality checks
   # OFA11 remote sensing reflectance
   OFA11_data = ncfile.createVariable('OFA11',dtype('float64').char,dims,fill_value=fill_value)
   OFA11_data[:] = failed_quality_fdata
   OFA11_data.units ='count'
   OFA11_data.missing_value = missing_value
   OFA11_data.valid_min = 0.0
   OFA11_data.valid_max = 1.0
   OFA11_data.scale_factor = 1.0
   OFA11_data.add_offset = 0.0
   #OFA11_data.standard_name = "data_values_failed_TS_ranges"
   OFA11_data.long_name = "Layer of data values that have failed the OceanFluxGHG TS data ranges"
 
 
    # set some global attributes
   setattr(ncfile, 'Conventions', 'CF-1.6') 
   setattr(ncfile, 'Institution', 'European Space Agency (ESA) OceanFlux Greenhouse Gases - Heriot Watt University; Plymouth Marine Laboratory (PML); Environmental Research Institute, Thurso; Institut Francais de Recherches pour l\'Exploitation de la Mer (IFREMER) and The UK National Oceanographic Centre (NOC)') 
   setattr(ncfile, 'contact', 'email: oceanfluxghg@pml.ac.uk') 

    # put in the contents of the configuration file into the netcdf
   setattr(ncfile, 'processing_hostname',hostname)
   setattr(ncfile, 'processing_time',processing_time)   
   setattr(ncfile, 'config_file',config_file)
   setattr(ncfile, 'year_of_run',year)
   setattr(ncfile, 'process_layers_off',process_layers_off)
   setattr(ncfile, 'flux_model',flux_model)
   setattr(ncfile, 'sst_gradients',sst_gradients)
   setattr(ncfile, 'use_sstskin',use_sstskin)
   setattr(ncfile, 'use_sstfnd',use_sstfnd)
   setattr(ncfile, 'saline_skin_value',saline_skin_value)
   setattr(ncfile, 'windu10_infile',windu10_infile)
   setattr(ncfile, 'windu10_prod',windu10_prod)
   setattr(ncfile, 'sstskin_infile', sstskin_infile)
   setattr(ncfile, 'sstskin_prod',sstskin_prod)
   setattr(ncfile, 'sstfnd_infile',sstfnd_infile)
   setattr(ncfile, 'sstfnd_prod',sstfnd_prod)
   setattr(ncfile, 'sigma0_infile',sigma0_infile)
   setattr(ncfile, 'sigma0_prod',sigma0_prod)
   setattr(ncfile, 'sig_wv_ht_infile',sig_wv_ht_infile)
   setattr(ncfile, 'sig_wv_ht_prod',sig_wv_ht_prod)
   setattr(ncfile, 'ice_infile',ice_infile)
   setattr(ncfile, 'ice_prod',ice_prod)
   setattr(ncfile, 'rain_infile',rain_infile)
   setattr(ncfile, 'rain_data_selection',rain_data_selection)
   setattr(ncfile, 'rain_prod',rain_prod)
   setattr(ncfile, 'pressure_infile',pressure_infile)
   setattr(ncfile, 'pressure_prod',pressure_prod)
   setattr(ncfile, 'salinity_data_selection',salinity_data_selection)
   setattr(ncfile, 'salinity_infile',salinity_infile)
   setattr(ncfile, 'salinity_prod',salinity_prod)
   setattr(ncfile, 'pco2_data_selection',pco2_data_selection)
   setattr(ncfile, 'pco2_infile',pco2_infile)
   setattr(ncfile, 'pco2_prod',pco2_prod)
   setattr(ncfile, 'pco2_sst_prod',pco2_sst_prod)
   setattr(ncfile, 'vco2_prod',vco2_prod)
   setattr(ncfile, 'bio_infile',bio_infile)
   setattr(ncfile, 'bio_prod',bio_prod)
   setattr(ncfile, 'susp_particles_infile',susp_particles_infile)
   setattr(ncfile, 'atlantic_ocean_maskfile',atlantic_ocean_maskfile)
   setattr(ncfile, 'pacific_ocean_maskfile',pacific_ocean_maskfile)
   setattr(ncfile, 'southern_ocean_maskfile',southern_ocean_maskfile)
   setattr(ncfile, 'indian_ocean_maskfile',indian_ocean_maskfile)
   setattr(ncfile, 'longhurst_maskfile',longhurst_maskfile)
   setattr(ncfile, 'sstgrad_infile',sstgrad_infile)
   setattr(ncfile, 'sstgrad_prod',sstgrad_prod)
   setattr(ncfile, 'random_noise_windu10',random_noise_windu10)
   setattr(ncfile, 'random_noise_sstskin',random_noise_sstskin)
   setattr(ncfile, 'random_noise_sstfnd',random_noise_sstfnd)
   setattr(ncfile, 'random_noise_pco2',random_noise_pco2)
   setattr(ncfile, 'bias_windu10',bias_windu10)
   setattr(ncfile, 'bias_windu10_value',bias_windu10_value)
   setattr(ncfile, 'bias_sstskin',bias_sstskin)
   setattr(ncfile, 'bias_sstskin_value',bias_sstskin_value)
   setattr(ncfile, 'bias_sstfnd',bias_sstfnd)
   setattr(ncfile, 'bias_sstfnd_value',bias_sstfnd_value)
   setattr(ncfile, 'bias_pco2',bias_pco2)
   setattr(ncfile, 'bias_pco2_value',bias_pco2_value)
   setattr(ncfile, 'bias_k',bias_k)
   setattr(ncfile, 'bias_k_value',bias_k_value)
   setattr(ncfile, 'bias_k_percent',bias_k_percent)
   setattr(ncfile, 'bias_k_biology_value',bias_k_biology_value)
   setattr(ncfile, 'bias_k_wind_value',bias_k_wind_value)
   setattr(ncfile, 'k_parameterisation',k_parameterisation)
   setattr(ncfile, 'k_generic_sc',k_generic_sc) 
   setattr(ncfile, 'k_generic_a0',k_generic_a0)
   setattr(ncfile, 'k_generic_a1',k_generic_a1)
   setattr(ncfile, 'k_generic_a2',k_generic_a2)
   setattr(ncfile, 'k_generic_a3',k_generic_a3)
   setattr(ncfile, 'kd_weighting',kd_weighting)
   setattr(ncfile, 'kb_weighting',kb_weighting)
   setattr(ncfile, 'kb_asymmetry',kb_asymmetry)
   setattr(ncfile, 'bias_sstskin_due_rain',bias_sstskin_due_rain)
   setattr(ncfile, 'bias_sstskin_due_rain_value',bias_sstskin_due_rain_value)
   setattr(ncfile, 'bias_sstskin_due_rain_intensity',bias_sstskin_due_rain_intensity)
   setattr(ncfile, 'bias_sstskin_due_rain_wind',bias_sstskin_due_rain_wind)
   setattr(ncfile, 'rain_wet_deposition',rain_wet_deposition)
   setattr(ncfile, 'k_rain_linear_ho1997',k_rain_linear_ho1997)
   setattr(ncfile, 'outfile',outfile)

    # close the file.
   ncfile.close()

   return 0


def add_noise(data, rmse_value, mean_value, no_elements):
# randomly adding noise to data array, based log normal distribution
# rmse value is used as the standard deviation of the noise function

   # intialise the random number generator
  for i in arange(no_elements):
     if ( (data[i] != missing_value) and (data[i] != 0.0) ):
      orig = data[i]
      value = log(data[i])
      stddev = float(rmse_value/orig)
      noise = normalvariate(0,stddev) # determines the random noise value based on the rain input data and the standard deviation of the uncertainty
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


def solubility(sstK, sal, DeltaT_fdata, nx, ny, flux_model):
 # solubility calculation
 # equation from Table A2 of Wanninkkhof, JGR, 1992
   
   sol = array([missing_value] * nx*ny)
   for i in arange(nx * ny):
      if ( (sstK[i] != missing_value) and (sal[i] != missing_value) and (sstK[i] > 0.0) ):
         #print "1sstskinK_fdata: (%d,%d) %d %f log:%f\n" %(nx, ny,i,sstK_fdata[i]/100.0,log(sstK_fdata[i]/100.0))
         sol[i] = -60.2409 + ( 93.4517*(100.0 / sstK[i]) ) + (23.3585 * (log(sstK[i]/100.0))) + (sal[i] * (0.023517 + ( (-0.023656)*(sstK[i]/100.0)) + (0.0047036*( (sstK[i]/100.0)*(sstK[i]/100.0) ) ) ) )
         sol[i] = exp(sol[i])
             # flux_model is a switch to remove Delta_T-Sb component - ie selects use of RAPID or EQUILIBRIUM flux models from Woolf et al., 2012
         if flux_model != 2:
            DeltaT_fdata[i] = 0.0
         sol[i] = sol[i] * (1 - (0.015*DeltaT_fdata[i]))
      else:
         sol[i] = missing_value
   return sol

def check_input(filename, dataset):
   function = "(check_input, main)"   
   if DEBUG:
      print "%s checking %s for %s data" % (function, filename, dataset)
   file = Dataset(filename,'r',clobber=False)
   for dimname, diminst in sorted(list(file.dimensions.iteritems())):
      if diminst.isunlimited():
         if DEBUG:
            print "%s %s dimension\t%s\t%s\tunlimited" % (function, dataset, dimname, len(diminst))
      else:
         if DEBUG:
            print "%s %s dimension\t%s\t%s" % (function, dataset, dimname, len(diminst))
              
   for varname, varinst in sorted(list(file.variables.iteritems())):
      if DEBUG:
         print "%s %s variable\t%s\t\t%s\t%s\t%s" % (function, dataset, varname, varinst.shape, varinst.dimensions, varinst.dtype)

   file.close()
   return 0

def check_output_dataset(data, name, failed_quality_fdata, nx, ny, min_range, max_range, sstskinC_fdata, sstfndC_fdata, windu10_fdata, sig_wv_ht_fdata, solfnd_fdata, solskin_fdata, k_fdata, concw_fdata, conca_fdata, sal_fdata):
   function = "(check_output, main)"
 # checking the contents of a dataset
   for i in arange(nx * ny):
      if ((data[i] != missing_value) and ((data[i] >max_range) or (data[i] <min_range))):
         #print "\n%s Dataset %s fails OceanFluxGHG TS table 7 valid limits, (defined min/max are: %lf/%lf, found %lf at grid point %d (%lf) exiting" % (function, name, min_range, max_range, data[i], i, i/nx)
         #print "\n%s Coincident data values sstskinC_fdata:%lf sstfndC_fdata:%lf windu10_fdata:%lf sig_wv_ht_fdata:%lf solfnd_fdata:%lf solskin_fdata:%lf k_fdata:%lf concw_fdata:%lf conca_fdata:%lf sal_fdata:%lf" % (function, sstskinC_fdata[i], sstfndC_fdata[i], windu10_fdata[i], sig_wv_ht_fdata[i], solfnd_fdata[i], solskin_fdata[i], k_fdata[i], concw_fdata[i], conca_fdata[i], sal_fdata[i])
         if failed_quality_fdata[i] == missing_value:
            failed_quality_fdata[i] = 1 # first entry so need to initiase
         else:
            failed_quality_fdata[i] += 1
         #sys.exit(1) # commented out to identify quantity of rogue values
   return 0

def average_sstgrad(data, n, missing_value):
    '''averages arr into superpixels each consisting of the mean of a n x n
    window in arr. Seems to be intended to go from a 0.5 x 0.5 degree grid to a
    1 x 1 degree grid, in which case n must be set to 2. Checks for and ignores
    values equal to missing_value.'''
    function = "(average_sstgrad, main)"
    print "%s Averaging sstgrad_fdata into 1x1 degree grid (N=%d)" % (function, n)
    nj0, ni0 = data.shape
    nj, ni = [nj0 / n, ni0 / n]
    if nj * n != nj0 or ni * n != ni0:
       print "Dimensions ", nj0, ni0, " indivisible by ", n
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
    return out

def check_dimensions(test_nx, test_ny, name, ref_nx, ref_ny):
   function = "(check_dimensions, main)"
   if test_nx == ref_nx and test_ny == ref_ny:
      if DEBUG:
         print "\n%s Input data (%s) have identical dimensions to reference values (%s, %s) "% (function, name, test_nx, test_ny)
   else :
      print "\n%s Input data have non-identical to reference, exiting (%s: %s, %s and reference dataset: %s, %s)" % (function, name, test_nx, test_ny, ref_nx, ref_ny)
      sys.exit(1)

def flip_data(data, this_variable):
    '''#IGA - for a netcdf data set, determine whether latitude orientation matches 'taka' and if not, flip the variable provided using flipud'''
    data_latitude_prod = [str(x) for x in data.variables.keys() if 'lat' in str(x).lower()] #finds the correct latitude name for data
    data_lat = data.variables[data_latitude_prod[0]]

    if len(data_lat.shape)<2 and data_lat[0]<0:#IGA - if true, it is a vector that is in opposite orientation to 'taka'
       flipped = True
       this_variable_out = flipud(this_variable)
    else: 
       this_variable_out = this_variable
       flipped = False

    return this_variable_out,flipped

 # function name definition for stderr and stdout messages
function = "(ofluxghg-flux-calc, main)"

 # checking arguments
if len(sys.argv)<81:
   print "%s Incorrect number of input parameters, exiting.\n" % (function)
   print "%s Usage: ofluxghg-flux-calc.py <sstskin-netcdf> <sstfnd-netcdf> <windu10-netcdf> <pressure-netcdf> <pco2-netcdf> <sal-netcdf> <sigma0-netcdf> <sig_wv_ht-netcdf> <rain-netcdf> <bio-netcdf> <susp_particles-netfcdf> <ice-netcdf> <atlantic-ocean-netcdf> <pacific-ocean-netcdf> <southern-ocean-netcdf> <indian-ocean-netcdf> <longhurst-netcdf> <sstgrad-netcdf> <output-netcdf> <sstskin-prod> <sstfnd-prod> <windu10-prod> <sigm0-prod> <sig_wv_ht-prod> <ice-prod> <pco2-prod> <pco2_sst-prod> <vco2-prod> <pressure-prod> <salinity-prod> <rain-prod> <bio-prod> <susp_particles-prod> <oceans_mask-prod> <longhurst_mask-prod> <sstgrad-prod> <year> <random_noise_windu10> <random_noise_sstskin> <random_noise_sstfnd> <random_noise_pco2> <bias_windu10> <bias_sstskin> <bias_sstfnd> <bias_pco2> <bias_k> <bias_windu10_value> <bias_sstskin_value> <bias_sstfnd_value> <bias_pco2_value> <bias_k_value> <bias_k_percent> <bias_k_biology_value> <bias_k_wind_value> <pco2_data_selection> <salinity_data_selection> <rain_data_selection> <k_sc_number> <k_parameterisation> <k_generic_a0> <k_generic_a1> <k_generic_a2> <k_generic_a3> <kd_weighting> <kb_weighting> <kb_asymmetry> <bias_rain_sstskin> <bias_rain_sstskin_value> <bias_rain_sstskin_intensity> <bias_rain_sstskin_wind> <rain_wet_deposition> <k_rain_linear_ho1997> <k_rain_nonlinear_h2012> <flux_model> <sst_gradients> <use_sstskin> <use_sstfnd> <saline_skin_value> <process_layers_off> <config_file> <hostname> <processing_time>\n" % (function)
   sys.exit(1)

 # passing arguments
sstskin_infile = sys.argv[1]
sstfnd_infile = sys.argv[2]
windu10_infile = sys.argv[3]
pressure_infile = sys.argv[4]
pco2_infile = sys.argv[5]
salinity_infile = sys.argv[6]
sigma0_infile = sys.argv[7]
sig_wv_ht_infile = sys.argv[8]
rain_infile = sys.argv[9]
bio_infile = sys.argv[10]
susp_particles_infile = sys.argv[11]
ice_infile = sys.argv[12]
atlantic_ocean_maskfile = sys.argv[13]
pacific_ocean_maskfile = sys.argv[14]
southern_ocean_maskfile = sys.argv[15]
indian_ocean_maskfile = sys.argv[16]
longhurst_maskfile = sys.argv[17]
sstgrad_infile = sys.argv[18]
outfile = sys.argv[19]
sstskin_prod = sys.argv[20]
sstfnd_prod = sys.argv[21]
windu10_prod = sys.argv[22]
sigma0_prod = sys.argv[23]
sig_wv_ht_prod = sys.argv[24]
ice_prod = sys.argv[25]
pco2_prod = sys.argv[26]
pco2_sst_prod = sys.argv[27]
vco2_prod = sys.argv[28]
pressure_prod = sys.argv[29]
salinity_prod = sys.argv[30]
rain_prod = sys.argv[31]
bio_prod = sys.argv[32]
susp_particles_prod = sys.argv[33]
oceans_mask_prod = sys.argv[34]
longhurst_mask_prod = sys.argv[35]
sstgrad_prod = sys.argv[36]
year = int(sys.argv[37])
random_noise_windu10 = int(sys.argv[38])
random_noise_sstskin = int(sys.argv[39])
random_noise_sstfnd = int(sys.argv[40])
random_noise_pco2 = int(sys.argv[41])
bias_windu10 = int(sys.argv[42])
bias_sstskin = int(sys.argv[43])
bias_sstfnd = int(sys.argv[44])
bias_pco2 = int(sys.argv[45])
bias_k = int(sys.argv[46])
bias_windu10_value = float(sys.argv[47])
bias_sstskin_value = float(sys.argv[48])
bias_sstfnd_value = float(sys.argv[49])
bias_pco2_value = float(sys.argv[50])
bias_k_value = float(sys.argv[51])
bias_k_percent = int(sys.argv[52])
bias_k_biology_value = float(sys.argv[53])
bias_k_wind_value = float(sys.argv[54])
pco2_data_selection = int(sys.argv[55])
salinity_data_selection = int(sys.argv[56])
rain_data_selection = int(sys.argv[57])
k_parameterisation = int(sys.argv[58])
k_generic_sc = float(sys.argv[59])
k_generic_a0 = float(sys.argv[60])
k_generic_a1 = float(sys.argv[61])
k_generic_a2 = float(sys.argv[62])
k_generic_a3 = float(sys.argv[63])
kd_weighting = float(sys.argv[64])
kb_weighting = float(sys.argv[65])
kb_asymmetry = float(sys.argv[66])
bias_sstskin_due_rain = int(sys.argv[67])
bias_sstskin_due_rain_value = float(sys.argv[68])
bias_sstskin_due_rain_intensity = float(sys.argv[69])
bias_sstskin_due_rain_wind = float(sys.argv[70])
rain_wet_deposition = int(sys.argv[71])
k_rain_linear_ho1997 = int(sys.argv[72])
k_rain_nonlinear_h2012 = int(sys.argv[73])
flux_model = int(sys.argv[74])
sst_gradients = int(sys.argv[75])
use_sstskin = int(sys.argv[76])
use_sstfnd = int(sys.argv[77])
saline_skin_value = float(sys.argv[78])
process_layers_off = int(sys.argv[79])
config_file = sys.argv[80]
hostname = sys.argv[81]
processing_time = sys.argv[82]


if DEBUG:
   print "\n%s Flux model: %d, k_parmaterisation: %d" % (function, flux_model, k_parameterisation)
   print "\n%s Using the following files %s (sstskin), %s (sstfnd) %s (windu10) %s (sal) %s (rain) %s (bio) %s (sstgrad)" % (function, sstskin_infile, sstfnd_infile, windu10_infile, salinity_infile, rain_infile, bio_infile, sstgrad_infile)
   print "\n%s Using the following products %s (sstskin), %s (sstfnd) %s (windu10) %s (msl) %s (sal) %s (rain) %s (bio) %s (sstgrad)" % (function, sstskin_prod, sstfnd_prod, windu10_prod, pressure_prod, salinity_prod, rain_prod, bio_prod, sstgrad_prod)


if TAKAHASHI_DRIVER:
   print "\n%s Using TAKAHASHI_DRIVER option: all input data are from the Takahashi et al., 2009 climatology" % (function)


 # checking that the input files exist

 # sstskin
try:
   os.stat(sstskin_infile)
except:
   print "\n%s sstskin-netcdf inputfile %s does not exist" % (function, sstskin_infile)
   sys.exit(1)

 # sstsfnd
try:
   os.stat(sstfnd_infile)
except:
   print "\n%s sstfnd-netcdf inputfile %s does not exist" % (function, sstfnd_infile)
   sys.exit(1)

 # windu10
try:
   os.stat(windu10_infile)
except:
   print "\n%s windu10-netcdf inputfile %s does not exist" % (function, windu10_infile)
   sys.exit(1)

try:
   os.stat(pco2_infile)
except:
   print "\n%s pco2-netcdf inputfile %s does not exist" % (function, pco2_infile)
   sys.exit(1)

 # salinity
try:
   os.stat(salinity_infile)
except:
   print "\n%s sal-netcdf inputfile %s does not exist" % (function, salinity_infile)
   sys.exit(1)

 # sigma0
try:
   os.stat(sigma0_infile)
except:
   print "\n%s sigma0-netcdf inputfile %s does not exist" % (function, sigma0_infile)
   sys.exit(1)

 # sig_wv_ht
try:
   os.stat(sig_wv_ht_infile)
except:
   print "\n%s sig_wv_ht-netcdf inputfile %s does not exist" % (function, sig_wv_ht_infile)
   sys.exit(1)


if process_layers_off == 0:
   # biological activity
   try:
      os.stat(bio_infile)
   except:
      print "\n%s bio-netcdf inputfile %s does not exist" % (function, bio_infile)
      sys.exit(1)

    # suspended particulates
   try:
      os.stat(susp_particles_infile)
   except:
      print "\n%s susp_particles-netcdf inputfile %s does not exist" % (function, susp_particles_infile)
      sys.exit(1)


 # atlantic oceans mask file
try:
   os.stat(atlantic_ocean_maskfile)
except:
   print "\n%s atlantic-ocean-netcdf inputfile %s does not exist" % (function, atlantic_ocean_maskfile)
   sys.exit(1)

 #pacific ocean mask
try:
   os.stat(pacific_ocean_maskfile)
except:
   print "\n%s pacific-ocean-netcdf inputfile %s does not exist" % (function, pacific_ocean_maskfile)
   sys.exit(1)

 #southern ocean mask
try:
   os.stat(southern_ocean_maskfile)
except:
   print "\n%s southern-ocean-netcdf inputfile %s does not exist" % (function, southern_ocean_maskfile)
   sys.exit(1)

 #indian ocean mask
try:
   os.stat(indian_ocean_maskfile)
except:
   print "\n%s indian-ocean-netcdf inputfile %s does not exist" % (function, indian_ocean_maskfile)
   sys.exit(1)

   
 # longhurst mask file
try:
   os.stat(longhurst_maskfile)
except:
   print "\n%s longhurst-netcdf inputfile %s does not exist" % (function, longhurst_maskfile)
   sys.exit(1)

 # sst gradient data
try:
   os.stat(sstgrad_infile)
except:
   print "\n%s sstgrad-netcdf inputfile %s does not exist" % (function, sstgrad_infile)
   sys.exit(1)

 # pressure data
try:
   os.stat(pressure_infile)
except:
   print "\n%s pressure-netcdf inputfile %s does not exist" % (function, pressure_infile)
   sys.exit(1)

 # rain data
try:
   os.stat(rain_infile)
except:
   print "\n%s rain-netcdf inputfile %s does not exist" % (function, rain_infile)
   sys.exit(1)

 # ice data
try:
   os.stat(ice_infile)
except:
   print "\n%s ice-netcdf inputfile %s does not exist" % (function, ice_infile)
   sys.exit(1)

print "\n%s Input files all exist, now checking contents" % (function)


 # display details of sstskin data input
check_input( sstskin_infile, sstskin_prod)
if use_sstfnd == 1:
   check_input(sstfnd_infile, sstfnd_prod)
check_input(ice_infile, ice_prod)
check_input(windu10_infile, windu10_prod)
check_input(pressure_infile, pressure_prod)
check_input(salinity_infile, salinity_prod)
check_input(sigma0_infile, sigma0_prod)
check_input(sig_wv_ht_infile, sig_wv_ht_prod)
check_input(rain_infile, rain_prod)
check_input(atlantic_ocean_maskfile, oceans_mask_prod)
check_input(pacific_ocean_maskfile, oceans_mask_prod)
check_input(southern_ocean_maskfile, oceans_mask_prod)
check_input(indian_ocean_maskfile, oceans_mask_prod)
check_input(longhurst_maskfile, longhurst_mask_prod)

if process_layers_off == 0:
   check_input(bio_infile, bio_prod)
   check_input(susp_particles_infile, susp_particles_prod)

 #
 # open datasets
 #
 # note some of the datasets have latitude information which is +90 to -90
 # and others have -90 to +90, hence use of numpy flipud method

 # if driving using TAKAHASHI_DRIVER, then changing skin input data
if TAKAHASHI_DRIVER:
   sstskin_prod = 'SST_t'
   sstskin_infile = pco2_infile
   latitude_prod = 'latitude'
   longitude_prod = 'longitude'
   windu10_prod = 'wind_t'
   windu10_infile = pco2_infile
   pressure_infile = pco2_infile
   pressure_prod ='pressure'
   ice_prod = 'sea_ice_coverage'
   ice_infile = pco2_infile   
else:
   latitude_prod = 'lat'
   longitude_prod = 'lon'



 # open sstskin dataset   
with Dataset(sstskin_infile,'r') as data:
   try:
    print sstskin_prod
     # passing sstskin_prod from command line
    sstskinK = data.variables[sstskin_prod]
    sstskinK_n, sstskinK_ny, sstskinK_nx = sstskinK.shape
    if DEBUG:
       print "\nsstskin data dimentions: %d %d %d" % (sstskinK_n, sstskinK_nx, sstskinK_ny)
    sstskinK_data = sstskinK[0,:,:]
    
    sstskinK_data,flipped = flip_data(data,sstskinK_data)#If necessary - data will be flipped using flipud
    
    try:
      sstskinK_fill_value = sstskinK._FillValue  
    except:
      sstskinK_fill_value = missing_value
    
     # assumes sstskin_prod ends in something like '_mean'
     # so need to load '_count and '_stddev' fields as well which are assumed to exist
    if TAKAHASHI_DRIVER != True:
       m = re.match(r"(\w+\_)\w+", sstskin_prod)
       sstskinK_stddev_prod = m.group(1) + 'stddev'
       sstskinK_stddev_data = data.variables[sstskinK_stddev_prod][0,:,:]
       sstskinK_count_prod = m.group(1) + 'count'
       sstskinK_count_data = data.variables[sstskinK_count_prod][0,:,:]
       if flipped:
          sstskinK_stddev_data = flipud(sstskinK_stddev_data)
          sstskinK_count_data = flipud(sstskinK_count_data)
    else:
       sstskinK_stddev_data = array([missing_value] * sstskinK_nx*sstskinK_ny)
       sstskinK_count_data = array([missing_value] * sstskinK_nx*sstskinK_ny)
   except:
    print "\n%s data variable '%s' or '%s' or '%s' is/are missing from sstskin-netcdf input (%s)" % (function, sstskin_prod, sstskinK_stddev_prod, sstskinK_count_prod, sstskin_infile)
    sys.exit(1)



   try:
      latitude_data = data.variables[latitude_prod][:]
      longitude_data = data.variables[longitude_prod][:]
   except:
      print "\n%s data variable 'lat' or 'lon' missing from sstskin-netcdf input (%s)" % (function, sstskin_infile)
      sys.exit(1)
   if len(latitude_data.shape)<2:#IGA If long and lat are vectors
      if latitude_data[0]<0:#IGA - it is a vector that is in opposite orientation to 'taka'
         latitude_data = flipud(latitude_data) # IFREMER latitude data goes from +ve to -ve, hence flip

      latitude_grid = meshgrid(ones((len(longitude_data))),latitude_data)[1]#IGA - gridding vectored geo data so that they can be treated the sam way as non-regular grids
      longitude_grid = meshgrid(longitude_data,ones((len(latitude_data))))[1]#IGA
   else:
      latitude_grid = latitude_data#IGA - geo data are already a grid
      longitude_grid = longitude_data#IGA

   try:
      time_data = data.variables['time'][:]
   except:
      print "\n%s data variable 'time' missing from sstskin-netcdf input (%s)" % (function, sstskin_infile)
      sys.exit(1)

 # open salinity dataset
with Dataset(salinity_infile,'r') as data:
   try:
      sal = data.variables[salinity_prod]
      sal_data = sal[0,:,:]#IGA
      sal_data,flipped = flip_data(data,sal_data)#If necessary - data will be flipped using flipud #This will error is in-situ data does not have associated latitude/longitude
      sal_ny, sal_nx = sal_data.shape
   except:
      print "\n%s data variable '%s' missing from sal-netcdf input (%s)" % (function, salinity_prod, salinity_infile)
      sys.exit(1)

 # open pco2 dataset
with Dataset(pco2_infile,'r') as data:
   try:
      vco2_air = data.variables[vco2_prod]
      pco2_n, pco2_ny, pco2_nx = vco2_air.shape
      vco2_air_data = vco2_air[0,:,:]

      vco2_air_data,flipped = flip_data(data,vco2_air_data)#If necessary - data will be flipped using flipud

   except:
      try:#IGA added to retrieve vco2 data from salinity infile if it does not exist in suggested location
        with Dataset(salinity_infile,'r') as data2:
            vco2_air = data2.variables[vco2_prod]
            pco2_n, pco2_ny, pco2_nx = vco2_air.shape
            vco2_air_data = vco2_air[0,:,:]
            if salinity_data_selection != 0: # Not takahashi data array so it will need flipping 
                vco2_air_data = flipud(vco2_air_data)
                flipped = True
      except:
          print "\n%s data variable %s missing from pco2-netcdf (%s) and salinity-netcdf input (%s)" % (function, vco2_prod, pco2_infile, salinity_infile)
          sys.exit(1)

   try:
      pco2_sw = data.variables[pco2_prod]
      pco2_sw_data = pco2_sw[0,:,:]
       # stddev information only exists for SOCAT data
      #print "\n\npco2_data_selection %s (%s)\n\n" % (pco2_data_selection, pco2_prod) 
      if pco2_data_selection == 2:
         m = re.match(r"(\w+\d?\_)\w+$", pco2_prod)
         try: 
             pco2_sw_stddev_prod = m.group(0) + '_std'
             pco2_sw_stddev_data = data.variables[pco2_sw_stddev_prod][0,:,:]
         except(KeyError): 
             pco2_sw_stddev_prod = m.group(1) + 'std'#IGA changed m.group(1) to m.group(0) as _mean not always present in these in-situ data
             pco2_sw_stddev_data = data.variables[pco2_sw_stddev_prod][0,:,:]

         if flipped:       
            pco2_sw_stddev_data = flipud(pco2_sw_stddev_data)
            pco2_sw_data = flipud(pco2_sw_data)
      if TAKAHASHI_DRIVER:
         pco2_air = data.variables['pCO2_air']
         pco2_air_data = pco2_air[0,:,:]
   except:
      print "\n%s data variable %s or _std variant, %s missing from pco2-netcdf input (%s) (pco2_data_selection %s)" % (function, pco2_prod, pco2_sw_stddev_prod, pco2_infile, pco2_data_selection)
      sys.exit(1)

   try:
      sstpco2 = data.variables[pco2_sst_prod]
      sstpco2_data = sstpco2[0,:,:]
      if flipped:
         sstpco2_data = flipud(sstpco2_data)
   except:
      print "\n%s data variable %s missing from pco2-netcdf input (%s)" % (function, pco2_sst_prod, pco2_infile)
      sys.exit(1)

 # loading XCO2 from 2000 for use with SOCAT data
 # assuming salinity file is Takahashi salinity
 #IGA - problem here understanding why this vco2_air is defined differently to that in the configuration file
# if pco2_data_selection != 0 and pco2_data_selection != 3:
#    vco2_prod = 'vCO2_air'
#    with Dataset(salinity_infile,'r') as data: 
#       try:
#          vco2_air = data.variables[vco2_prod]
#          pco2_n, pco2_ny, pco2_nx = vco2_air.shape
#          vco2_air_data = vco2_air[0,:,:]         
#          if VERIFICATION_RUNS:
#             # added for testing/verification, not for normal operation
#             # forces SST for pCO2 to be Takahashi SST_t
#             pco2_sst_prod = 'SST_t'
#             sstpco2 = data.variables[pco2_sst_prod]
#             pco2_n, pco2_ny, pco2_nx = sstpco2.shape
#             sstpco2_data = sstpco2[0,:,:]      
#       except:
#          print "\n%s data variable %s missing from pco2-netcdf input (%s)" % (function, vco2_prod, salinity_infile)
#          sys.exit(1)

 # open biological dataset
 # GlobColour and CCI datasets don't have a time dimension ?!
if process_layers_off == 0:
   with Dataset(bio_infile,'r') as data:
      try:
         bio = data.variables[bio_prod]
         bio_ny, bio_nx = bio.shape
         bio_data = bio[:,:]
         try:
             bio_fill_value = bio._FillValue
         except:
             bio_fill_value = bio[:].fill_value#IGA_fill_value
         bio_data,flipped = flip_data(data,bio_data)#If necessary - data will be flipped using flipud
      except:
         print "\n%s data variable '%s' missing from bio-netcdf input (%s)" % (function, bio_prod, bio_infile)
         sys.exit(1)

 # opening the atlantic ocean maskfile
with Dataset(atlantic_ocean_maskfile,'r') as data:
   try:
      atlantic = data.variables[oceans_mask_prod]
      atlantic_ocean_data = atlantic[0,:,:]
   except:
      print "\n%s data variable '%s' missing from atlantic-oceans-netcdf input (%s)" % (function, oceans_mask_prod, atlantic_ocean_maskfile)
      sys.exit(1)

 # opening the pacific ocean maskfile
with Dataset(pacific_ocean_maskfile,'r') as data:
   try:
      pacific = data.variables[oceans_mask_prod]
      pacific_ocean_data = pacific[0,:,:]
   except:
      print "\n%s data variable '%s' missing from pacific-oceans-netcdf input (%s)" % (function, oceans_mask_prod, pacific_ocean_maskfile)
      sys.exit(1)

 # opening the southern ocean maskfile
with Dataset(southern_ocean_maskfile,'r') as data:
   try:
      southern = data.variables[oceans_mask_prod]
      southern_ocean_data = southern[0,:,:]
   except:
      print "\n%s data variable '%s' missing from southern-oceans-netcdf input (%s)" % (function, oceans_mask_prod, southern_ocean_maskfile)
      sys.exit(1)

 # opening the indian ocean maskfile
with Dataset(indian_ocean_maskfile,'r') as data:
   try:
      indian = data.variables[oceans_mask_prod]
      indian_ocean_data = indian[0,:,:]
   except:
      print "\n%s data variable '%s' missing from indian-oceans-netcdf input (%s)" % (function, oceans_mask_prod, indian_ocean_maskfile)
      sys.exit(1)


 # opening the longhurst maskfile
with Dataset(longhurst_maskfile,'r') as data:
   try:
      longhurst = data.variables[longhurst_mask_prod]
      longhurst_data = longhurst[0,:,:]
   except:
      print "\n%s data variable '%s' missing from longhurst-netcdf input (%s)" % (function, longhurst_mask_prod, longhurst_maskfile)
      sys.exit(1)

# GlobColour and CCI datasets don't have a time dimension ?!
if process_layers_off == 0:
   with Dataset(susp_particles_infile,'r') as data:
      try:
         susp_particles = data.variables[susp_particles_prod]
         susp_particles_ny, susp_particles_nx = susp_particles.shape
         susp_particles_data = susp_particles[:,:]
         susp_particles_data,flipped = flip_data(data,susp_particles_data)#If necessary - data will be flipped using flipud
         try:
             susp_particles_fill_value = susp_particles._FillValue
         except:
             susp_particles_fill_value = susp_particles[:].fill_value#IGA_fill_value
      except:
         print "\n%s data variable '%s' missing from susp_particles-netcdf input (%s)" % (function, susp_particles_prod, susp_particles_infile)
         sys.exit(1)

 # open sst gradient dataset
with Dataset(sstgrad_infile,'r') as data:
   try:
       # passing sstgrad_prod from command line
      sstgrad = data.variables[sstgrad_prod]
      sstgrad_n, sstgrad_ny, sstgrad_nx = sstgrad.shape
      sstgrad_data = sstgrad[0,:,:]
      sstgrad_data,flipped = flip_data(data,sstgrad_data)#If necessary - data will be flipped using flipud
   except:
      print "\n%s data variable %s missing from sstgrad-netcdf input (%s)" % (function, sstgrad_prod, sstgrad_infile)
      sys.exit(1)

 # open pressure dataset
with Dataset(pressure_infile,'r') as data:
   try:
       # passing pressure_prod from command line
      pressure = data.variables[pressure_prod]
      pressure_n, pressure_ny, pressure_nx = pressure.shape
      pres_data = pressure[0,:,:]
      pres_data,flipped = flip_data(data,pres_data)#If necessary - data will be flipped using flipud
      try:
          pres_fill_value = pressure._FillValue
      except:
          pres_fill_value = pressure[:].fill_value#IGA_fill_value
   except:
      print "\n%s data variable 'pressure' missing from pressure-netcdf input (%s)" % (function, pressure_infile)
      sys.exit(1)


 # open sstfnd dataset, only if use_sstfnd ==1 (if it ==0 then we are using sstskin +0.14 for sstfnd
sstfndK_fill_value = -1 # intialised as -1 incase use_sstfnd == 0
if use_sstfnd == 1:
   with Dataset(sstfnd_infile,'r') as data:
      try:
          # passing sstfnd_prod from command line
         sstfndK = data.variables[sstfnd_prod]
         sstfndK_n, sstfndK_ny, sstfndK_nx = sstfndK.shape
         sstfndK_data =  sstfndK[0,:,:]
         if mean(sstfndK_data)<200:#IGA-added if to check whether data are celsius or Kelvin
           sstfndK_data =  sstfndK_data+273.15

         sstfndK_data,flipped = flip_data(data,sstfndK_data)#If necessary - data will be flipped using flipud
         try:
             sstfndK_fill_value = sstfndK._FillValue 
         except:
             sstfndK_fill_value = sstfndK_data[:].fill_value#IGA_fill_value

          # assumes windu10_prod ends in something like '_mean'
          # so need to load '_coun't and '_stddev' fields as well which are assumed to exist
         m = re.match(r"(\w+\_)\w+", sstfnd_prod)
         sstfndK_stddev_prod = m.group(1) + 'stddev'
         sstfndK_stddev_data = data.variables[sstfndK_stddev_prod][0,:,:]
         sstfndK_count_prod = m.group(1) + 'count'
         sstfndK_count_data = data.variables[sstfndK_count_prod][0,:,:]
         if flipped:
           sstfndK_stddev_data = flipud(sstfndK_stddev_data)
           sstfndK_count_data = flipud(sstfndK_count_data)      
      except:
         print "\n%s data variable '%s' or '%s' or '%s' is/are missing from sstfnd-netcdf input (%s)" % (function, sstfnd_prod, sstfndK_stddev_prod, sstfndK_count_prod, sstfnd_infile)
         sys.exit(1)


 # open windu10 dataset 
with Dataset(windu10_infile,'r') as data:
   try:
      windu10 = data.variables[windu10_prod]   
      windu10_n, windu10_ny, windu10_nx = windu10.shape
      windu10_data = windu10[0,:,:]
      windu10_data,flipped = flip_data(data,windu10_data)#If necessary - data will be flipped using flipud
      try:
          windu10_fill_value = windu10._FillValue
      except:
          windu10_fill_value = windu10[:].fill_value#iGA_fill_value

       # assumes windu10_prod ends in something like '_mean'
       # so need to load '_coun't and '_stddev' fields as well which are assumed to exist
      if TAKAHASHI_DRIVER != True:
         m = re.match(r"(\w+\_)\w+", windu10_prod)
         windu10_stddev_prod = m.group(1) + 'stddev'
         windu10_stddev_data = data.variables[windu10_stddev_prod][0,:,:]
         windu10_count_prod = m.group(1) + 'count'
         windu10_count_data = data.variables[windu10_count_prod][0,:,:]

          # loading the second and third order wind moments
         windu10_moment2_prod = m.group(1) + 'moment_2'
         windu10_moment2_data = data.variables[windu10_moment2_prod][0,:,:]
         windu10_moment2_fill_value = data.variables[windu10_moment2_prod]._FillValue
      
         windu10_moment3_prod = m.group(1) + 'moment_3'
         windu10_moment3_data = data.variables[windu10_moment3_prod][0,:,:]
         windu10_moment3_fill_value = data.variables[windu10_moment3_prod]._FillValue
         if flipped: 
           windu10_stddev_data = flipud(windu10_stddev_data)
           windu10_count_data = flipud(windu10_count_data)
           windu10_moment2_data = flipud(windu10_moment2_data)
           windu10_moment3_data = flipud(windu10_moment3_data)
      else:
         windu10_count_data = array([missing_value] * windu10_nx*windu10_ny)
         windu10_stddev_data = array([missing_value] * windu10_nx*windu10_ny)
         windu10_moment2_data = array([missing_value] * windu10_nx*windu10_ny)
         windu10_moment3_data = array([missing_value] * windu10_nx*windu10_ny)
         windu10_moment3_fill_value = missing_value
         windu10_moment2_fill_value = missing_value
     
     
   except:
      print "\n%s data variable '%s' or '%s' or '%s' or '%s' or '%s' is/are missing from windu10-netcdf input (%s)" % (function, windu10_prod, windu10_stddev_prod, windu10_count_prod, windu10_moment2_prod, windu10_moment3_prod, windu10_infile)
      sys.exit(1)
 

 # open sigma0 dataset 
with Dataset(sigma0_infile,'r') as data:
   try:
      sigma0 = data.variables[sigma0_prod]   
      sigma0_n, sigma0_ny, sigma0_nx = sigma0.shape
      sigma0_data = sigma0[0,:,:]
      sigma0_data,flipped = flip_data(data,sigma0_data)#If necessary - data will be flipped using flipud
      sigma0_fill_value = sigma0._FillValue

       # assumes sigma0_prod ends in something like '_mean'
       # so need to load '_coun't and '_stddev' fields as well which are assumed to exist
      #m = re.match(r"(\w+\_)\w+", sigma0_prod)
      #sigma0_stddev_prod = m.group(1) + 'stddev'
      #sigma0_stddev_data = data.variables[sigma0_stddev_prod][0,:,:]
      #sigma0_stddev_data = flipud(sigma0_stddev_data)
      #sigma0_count_prod = m.group(1) + 'count'
      #sigma0_count_data = data.variables[sigma0_count_prod][0,:,:]
      #sigma0_count_data = flipud(sigma0_count_data)      
   except:
      print "\n%s data variable '%s' is missing from sigma0-netcdf input (%s)" % (function, sigma0_prod, sigma0_infile)
      sys.exit(1)


 # open sig_wv_ht dataset 
with Dataset(sig_wv_ht_infile,'r') as data:
   try:
      sig_wv_ht = data.variables[sig_wv_ht_prod]   
      sig_wv_ht_n, sig_wv_ht_ny, sig_wv_ht_nx = sig_wv_ht.shape
      sig_wv_ht_data = sig_wv_ht[0,:,:]
      sig_wv_ht_data,flipped = flip_data(data,sig_wv_ht_data)#If necessary - data will be flipped using flipud
      sig_wv_ht_fill_value = sig_wv_ht._FillValue
      
       # assumes sig_wv_ht_prod ends in something like '_mean'
       # so need to load '_coun't and '_stddev' fields as well which are assumed to exist
      m = re.match(r"(\w+\_)\w+", sig_wv_ht_prod)
      sig_wv_ht_stddev_prod = m.group(1) + 'stddev'
      sig_wv_ht_stddev_data = data.variables[sig_wv_ht_stddev_prod][0,:,:]
      sig_wv_ht_count_prod = m.group(1) + 'count'
      sig_wv_ht_count_data = data.variables[sig_wv_ht_count_prod][0,:,:]
      if flipped:
        sig_wv_ht_stddev_data = flipud(sig_wv_ht_stddev_data)
        sig_wv_ht_count_data = flipud(sig_wv_ht_count_data)

   except:
      print "\n%s data variable '%s' or '%s' or '%s' is/are missing from sig_wv_ht-netcdf input (%s)" % (function, sig_wv_ht_prod, sig_wv_ht_stddev_prod, sig_wv_ht_count_prod, sig_wv_ht_infile)
      sys.exit(1)

 # open rain dataset
rain_missing_value = -999.0
with Dataset(rain_infile,'r') as data:
   try:
      if rain_data_selection == 0: # TRMM data are in mm hr^-1
         rain = data.variables[rain_prod]   
         rain_n, rain_ny, rain_nx = rain.shape
         rain_data = rain[0,:,:]
         rain_data,flipped = flip_data(data,rain_data)#If necessary - data will be flipped using flipud
         rain_fill_value = rain._FillValue
         rain_missing_value = rain.missing_value
      elif rain_data_selection == 1: # GPCP data are in mm day^-1
         rain = data.variables[rain_prod]
         rain_nx, rain_ny, rain_n = rain.shape
         rain_data = rain[:,:,0]
         rain_data = transpose(rain_data)
         rain_data,flipped = flip_data(data,rain_data)#If necessary - data will be flipped using flipud
         try:
             rain_fill_value = rain._FillValue
         except:
             rain_fill_value = sig_wv_ht_fill_value
     ##Ians additional code to bring in the error data###########
         # print '\n%sIGA adding rain error data'
         # rain_err = data.variables['Error']
         # rain_err_data = rain_err[:,:,0]
         # rain_err_data = transpose(rain_err_data)
         # rain_err_data = flipud(rain_err_data)
         # rain_err_fill_value = rain_err._FillValue
         # print '\n%sIGA adding rain error data END'
     ##Ians additional code to bring in the err or data END##############
     
      else:
         print "\n%s rain data selection is unrecognised, exiting " % (function)
         sys.exit(1)
   except:
      print "\n%s data variable '%s' missing from rain-netcdf input (%s)" % (function, rain_prod, rain_infile)
      sys.exit(1)

 # open ice dataset 
with Dataset(ice_infile,'r') as data:
   try:
      ice = data.variables[ice_prod]
      try: 
        ice_n, ice_ny, ice_nx = ice.shape
        ice_data = ice[0,:,:]
      except: #IGA Update - supports CERSAT ice data which does not have a time dimension
        ice_ny, ice_nx = ice.shape
        ice_data = ice[:].view(ma.MaskedArray)
        ice_data[ice_data<0]=ma.masked#CERSAT uses -1 for missing value. As it is a fraction/percentage, negative numbers will always be missing values
      ice_data,flipped = flip_data(data,ice_data)#If necessary - data will be flipped using flipud
      try:
          ice_fill_value = ice._FillValue
      except:
          ice_fill_value = ice.fill_value
   except:
      print "\n%s data variable '%s' is missing from ice-netcdf input (%s)" % (function, ice_prod, ice_infile)
      sys.exit(1)


 # year correction for pCO2/fCO2 data
 # if pco2_data_selection == 1, then using SOCAT data so no increment is added as SOCAT are normalised to 2010
 # if pco2_data_selection == 0, then using Takahashi et al data which are normalised to 2000
if pco2_data_selection == 0:
 # increment to be added to pCO2w and pCO2a in Takahashi climatology
 # 1.5uatm per year since 2000 as set by Takahashi et al 2009.
 # Takahashi data are normalised to the year 2000
   pco2_increment = (year - 2000.0) * 1.5
   print "%s year: %d Takahashi pCO2_increment: %lf (uatm) (pco2_data_selection: %d)" % (function, year, pco2_increment,pco2_data_selection)
elif (pco2_data_selection == 1): # signifies that SOCAT pco2 data are being used
    # need handling for SOCAT data for years before 2010
    # 1.5uatm per year since 2000 as set by Takahashi et al 2009.
    # so assuming -1.5uatm for each year prior to 2010
   pco2_increment = (year - 2010.0) * 1.5
   print "%s year: %d SOCAT pCO2_increment: %lf (uatm) (pco2_data_selection: %d)" % (function, year, pco2_increment,pco2_data_selection)
elif (pco2_data_selection == 2): # signifies that SOCAT fco2 data are being used
    # need handling for SOCAT data for years before 2010
    # 1.5uatm per year since 2000 as set by Takahashi et al 2009.
    # so assuming -1.5uatm for each year prior to 2010
   pco2_increment = (year - 2010.0) * 1.5
   print "%s year: %d SOCAT fCO2_increment: %lf (uatm) (pco2_data_selection: %d)" % (function, year, pco2_increment,pco2_data_selection)
elif (pco2_data_selection == 3): # signifies that in situ data are being used
   pco2_increment = 0.0
   print "%s year: %d insitu fCO2_increment: %lf (uatm) (pco2_data_selection: %d)" % (function, year, pco2_increment,pco2_data_selection)
else:
   pco2_increment = 0.0

 # all data dimensions need to be the same
nx = sstskinK_nx
ny = sstskinK_ny

 # checking dimensions of all datasets are identical
check_dimensions(windu10_nx, windu10_ny, "windu10", nx, ny)
check_dimensions(pco2_nx, pco2_ny, "pco2", nx, ny)
check_dimensions(sigma0_nx, sigma0_ny, "sigma0", nx, ny)
check_dimensions(sig_wv_ht_nx, sig_wv_ht_ny, "sig_wv_ht", nx, ny)
check_dimensions(rain_nx, rain_ny, "rain", nx, ny)

if use_sstfnd == 1:
    check_dimensions(sstfndK_nx, sstfndK_ny, "sstfndK", nx, ny)
check_dimensions(pressure_nx, pressure_ny, "pres", nx, ny)
check_dimensions(sal_nx, sal_ny, "sal", nx, ny)

 # process layers checks
if process_layers_off == 0:
   check_dimensions(bio_nx, bio_ny, "bio", nx, ny)
   check_dimensions(susp_particles_nx, susp_particles_ny, "susp_particles", nx, ny)
    # average gradient data into the same spatial resolution as all other datasets
   sstgrad_data_avg = average_sstgrad(sstgrad_data, 4, missing_value)
    # these aren't the same dimensions, as are spatially degraded on the next line
   sstgrd_data_avg_ny, sstgrd_data_avg_nx = sstgrad_data_avg.shape
   check_dimensions(sstgrd_data_avg_nx, sstgrd_data_avg_ny, "sstgrad", nx, ny)
 
 
 # changing data into 1d arrays ready for analysis

sstskinK_fdata = ravel(sstskinK_data)
sstskinK_count_fdata = ravel(sstskinK_count_data)
sstskinK_stddev_fdata = ravel(sstskinK_stddev_data)

if ( (sst_gradients == 1 and use_sstfnd == 1) or (sst_gradients == 0 and use_sstfnd == 1) ):
    # using actual sstfnd
   sstfndK_fdata = ravel(sstfndK_data)
   sstfndK_count_fdata = ravel(sstfndK_count_data)
   sstfndK_stddev_fdata = ravel(sstfndK_stddev_data)
elif ( (sst_gradients == 1 and use_sstfnd == 0) or (sst_gradients == 0 and use_sstfnd == 0) ): # initialising as needed later to be filled with sstskin +0.14
   sstfndK_fdata = array([missing_value] * nx*ny)
   sstfndK_count_fdata = array([missing_value] * nx*ny)
   sstfndK_stddev_fdata = array([missing_value] * nx*ny)
else:
   print "\n%s sst_gradients (%d) and use_sstfnd (%d) combination in configuration not recognised, exitiing." % (function, sst_gradients, use_sstfnd)
   sys.exit(1)
   
windu10_fdata = ravel(windu10_data)
windu10_count_fdata = ravel(windu10_count_data)
windu10_stddev_fdata = ravel(windu10_stddev_data)
windu10_moment2_fdata = ravel(windu10_moment2_data)
windu10_moment3_fdata = ravel(windu10_moment3_data)

if TAKAHASHI_DRIVER == True:
 # need to generate the moment2 and moment3 data
   for i in arange(nx * ny):
      if windu10_fdata[i] != missing_value:
         windu10_moment2_fdata[i] = windu10_fdata[i]*windu10_fdata[i]
         windu10_moment3_fdata[i] = windu10_fdata[i]*windu10_fdata[i]*windu10_fdata[i]
      else:
         windu10_moment2_fdata[i] = missing_value
         windu10_moment3_fdata[i] = missing_value

if TAKAHASHI_DRIVER == True:
   pco2_air_fdata = ravel(pco2_air_data)
else: 
   pco2_air_fdata = array([missing_value] * nx*ny)

sigma0_fdata = ravel(sigma0_data)

sig_wv_ht_fdata = ravel(sig_wv_ht_data)
sig_wv_ht_count_fdata = ravel(sig_wv_ht_count_data)
sig_wv_ht_stddev_fdata = ravel(sig_wv_ht_stddev_data)

rain_fdata = ravel(rain_data)
sal_fdata = ravel(sal_data)
ice_fdata = ravel(ice_data)

vco2_air_fdata = ravel(vco2_air_data)
pres_fdata = ravel(pres_data)
sstpco2_fdata = ravel(sstpco2_data)

 # process layers
 # bio data is always needed
#bio_fdata = ravel(bio_data)
if process_layers_off == 0:
   bio_fdata = ravel(bio_data)
   susp_particles_fdata = ravel(susp_particles_data)
   atlantic_ocean_fdata = ravel(atlantic_ocean_data)
   pacific_ocean_fdata = ravel(pacific_ocean_data)
   southern_ocean_fdata = ravel(southern_ocean_data)
   indian_ocean_fdata = ravel(indian_ocean_data)
   longhurst_fdata = ravel(longhurst_data)
   sstgrad_fdata_avg = ravel(sstgrad_data_avg)

 # some specific pco2 conditions
pco2_sw_fdata = ravel(pco2_sw_data)
if (pco2_data_selection == 1):
    # signifies that we're using SOCAT data 
    pco2_sw_stddev_fdata = ravel(pco2_sw_stddev_data)
else:
    # this data field doesn't exist if we aren't using SOCAT data
    pco2_sw_stddev_fdata = array([missing_value] * nx*ny)

  
 # initialising some data structures for the calculations
scskin_fdata = array([missing_value] * nx*ny) # 1d array of floats
scfnd_fdata = array([missing_value] * nx*ny)
k_fdata = array([missing_value] * nx*ny) # 1d array of floats
solskin_fdata = array([missing_value] * nx*ny)
solfnd_fdata = array([missing_value] * nx*ny)
pH2O_fdata = array([missing_value] * nx*ny)
pco2_air_cor_fdata = array([missing_value] * nx*ny)
pco2_sw_cor_fdata = array([missing_value] * nx*ny)
sstskinC_fdata = array([missing_value] * nx*ny)
sstfndC_fdata = array([missing_value] * nx*ny)
dpco2_cor_fdata = array([missing_value] * nx*ny)
dpconc_cor_fdata = array([missing_value] * nx*ny)


salskin_fdata = array([missing_value] * nx*ny)

  # ensuring salinity data have a missing_value entry.
  # as smos salinity data doesn't have a missing_value or _FillValue entry !?

for i in arange(nx * ny):
   if (sal_fdata[i] >= 0.0) and (sal_fdata[i] <= 50.0):
      sal_fdata[i] = sal_fdata[i]
      if ( (sal_fdata[i] + saline_skin_value) <= 50.0):
          salskin_fdata[i] = sal_fdata[i] + saline_skin_value
   else:
      sal_fdata[i] = missing_value
      salskin_fdata[i] = missing_value

 # conversion of missing values to standard value, rather than variations that seem to exist in some of these data
for i in arange(nx * ny):
    # windu10 fdata   
   if (windu10_fdata[i] == windu10_fill_value): 
      windu10_fdata[i] = missing_value
   if (windu10_moment2_fdata[i] == windu10_moment2_fill_value): 
      windu10_moment2_fdata[i] = missing_value
   if (windu10_moment3_fdata[i] == windu10_moment3_fill_value): 
      windu10_moment3_fdata[i] = missing_value
   
    # sig_wv_ht_fdata
   if (sig_wv_ht_fdata[i] == sig_wv_ht_fill_value): 
      sig_wv_ht_fdata[i] = missing_value
    # sigma0_fdata
   if (sigma0_fdata[i] == sigma0_fill_value): 
      sigma0_fdata[i] = missing_value
    # pres_fdata
   if (pres_fdata[i] == pres_fill_value):
      pres_fdata[i] = missing_value
   if (ice_fdata[i] == ice_fill_value):
      ice_fdata[i] = missing_value
    # biological data
    # working with preliminary ESA OC CCI data so missing values aren't always correctly assigned
    # so also sanity checking data with <0.0 test
   if process_layers_off == 0:
      if (bio_fdata[i] == bio_fill_value):
         bio_fdata[i] == missing_value
      if (susp_particles_fdata[i] == susp_particles_fill_value):
         susp_particles_fdata[i] = missing_value         
    # rain_fdata 

   if ((rain_fdata[i] == rain_fill_value) or (rain_fdata[i] == rain_missing_value)):
      rain_fdata[i] = missing_value
   else: # conversion of data to mm hr^-1 (from mm day^-1)
      rain_fdata[i] /= 24.0
   
if TAKAHASHI_DRIVER:
 # sea ice units in the Takahashi are in %age so need to be scaled
   for i in arange(nx * ny):
      if ice_fdata[i] != missing_value:
         ice_fdata[i] /= 100.0


  # interpreting fnd_data option
cool_skin_difference = 0.14 #Relationship and value from Donlon et al., 2002; page 358, first paragraph
if TAKAHASHI_DRIVER:
   cool_skin_difference = 0.0

if sst_gradients == 0 and use_sstskin == 0 and use_sstfnd == 1:
   print "%s SST gradient handling is off, using SSTfnd data selection in configuration file for all components of the flux calculation (ignoring SSTskin data in configuration file)." % (function)
    #  actually copy sstfnd data into the sstskin dataset to make sure
   for i in arange(nx * ny):
      if sstfndK_fdata[i] != missing_value:
         sstskinK_fdata[i] = sstfndK_fdata[i]
         sstskinK_stddev_fdata[i] =  sstfndK_stddev_fdata[i]
         sstskinK_count_fdata[i] = sstfndK_count_fdata[i]
      else:
         sstskinK_fdata[i] = missing_value

#IGA added for the case where only foundation is provided and gradients are on------------------------------
elif sst_gradients == 1 and use_sstskin == 0 and use_sstfnd == 1:
   print "%s SUsing SSTfnd data selection with correction for skin temperature (SSTskin = SSTfnd - 0.14)(ignoring SSTskin data in configuration file)." % (function)
    #  actually copy sstfnd data into the sstskin dataset to make sure
   for i in arange(nx * ny):
      if sstfndK_fdata[i] != missing_value:
         sstskinK_fdata[i] = sstfndK_fdata[i]-cool_skin_difference
         #sstskinC_fdata[i] = sstfndC_fdata[i]-cool_skin_difference #IGA_temp changed as sstskinC has not yet been defined
         sstskinK_stddev_fdata[i] =  sstfndK_stddev_fdata[i]
         sstskinK_count_fdata[i] = sstfndK_count_fdata[i]
      else:
         sstskinK_fdata[i] = missing_value
#IGA added for the case where only foundation is provided and gradients are on------------------------------

elif sst_gradients == 0 and use_sstskin == 1 and use_sstfnd == 0:
   print "%s SST gradient handling is off, using SSTskin to derive SSTfnd (SSTfnd = SSTskin + 0.14) for flux calculation (ignoring SSTfnd data in configuration file)." % (function)
    #setting sstfnd_ data fields to skin values   
   for i in arange(nx * ny):
      if sstskinK_fdata[i] != missing_value:
         sstfndK_fdata[i] = sstskinK_fdata[i] + cool_skin_difference
         sstskinK_fdata[i] = sstskinK_fdata[i] + cool_skin_difference
         #sstfndC_fdata[i] = sstskinC_fdata[i] + cool_skin_difference     #IGA_temp changed as sstskinC has not yet been defined
         sstfndK_stddev_fdata[i] =  sstskinK_stddev_fdata[i]
         sstfndK_count_fdata[i] = sstskinK_count_fdata[i]
      else:
         sstfndK_fdata[i] = missing_value
         sstskinK_fdata[i] = missing_value
         sstskinC_fdata[i] = missing_value
elif sst_gradients == 1 and use_sstskin == 1 and use_sstfnd == 0:
   print "%s SST gradient handling is on, using SSTskin and SSTfnd = SSTskin + 0.14 for flux calculation (ignoring SSTfnd data in configuration file)." % (function)
    #setting sstfnd_ data fields to skin values  
   for i in arange(nx * ny):
      if sstskinK_fdata[i] != missing_value:
         sstfndK_fdata[i] = sstskinK_fdata[i] + cool_skin_difference
         sstfndK_stddev_fdata[i] =  sstskinK_stddev_fdata[i]
         sstfndK_count_fdata[i] = sstskinK_count_fdata[i]
      else:
         sstfndK_fdata[i] = missing_value
elif sst_gradients == 1 and use_sstskin == 1 and use_sstfnd == 1:
   print "%s SST gradient handling is on, using SSTfnd and SSTskin from the configuration file." % (function)
else:
   print "\n%s sst_gradients (%d), use_sstskin (%d) and use_sstfnd (%d) combination in configuration not recognised, exitiing." % (function, sst_gradients, use_sstskin, use_sstfnd)
   sys.exit(1)


 # quality filtering and conversion of SST datasets
if TAKAHASHI_DRIVER != True:
   for i in arange(nx * ny):
       # convert SSTfnd to centrigrade and store the result
      if (sstfndK_fdata[i] != missing_value):
         sstfndC_fdata[i] = sstfndK_fdata[i] - 273.15
      else:
         sstfndC_fdata[i] = missing_value
       # convert SSTskin to centrigrade and store the result
      if (sstskinK_fdata[i] != missing_value):
         sstskinC_fdata[i] = sstskinK_fdata[i] - 273.15
      else:
         sstskinC_fdata[i] = missing_value

 # quality filtering if using TAKAHASHI_DRIVER
 # TAKAHASHI data are in oC (and not Kelvin), so we have to apply the reverse
if TAKAHASHI_DRIVER == True:
   for i in arange(nx * ny):
      if (sstfndK_fdata[i] != missing_value):
         sstfndC_fdata[i] = sstfndK_fdata[i]
         sstfndK_fdata[i] = sstfndC_fdata[i] + 273.15
      else:
         sstfndC_fdata[i] = missing_value

      if (sstskinK_fdata[i] != missing_value):
         sstskinC_fdata[i] = sstskinK_fdata[i]
         sstskinK_fdata[i] = sstskinC_fdata[i] + 273.15
      else:
         sstskinC_fdata[i] = missing_value

 # ability to randomly perturb the input datasets
 # needed for the ensemble analyses
 # stddev of noise is using published RMSE for each dataset
 # all datasets are considered to have bias=0, hence mean of noise=0
if (random_noise_windu10 == 1):
   add_noise(windu10_fdata, 0.44, 0.0, nx*ny)
   add_noise(windu10_moment2_fdata, 0.44, 0.0, nx*ny)
   add_noise(windu10_moment3_fdata, 0.44, 0.0, nx*ny)
   print "%s Adding random noise to windu10_mean, windu10_moment2 and windu10_moment3 (mean 0.0, stddev 0.44 ms^-1 - assuming using ESA GlobWave data)" % (function)

if (random_noise_sstskin == 1):
   add_noise(sstskinK_fdata, 0.14, 0.0, nx*ny)
   print "%s Adding random noise to sstskin (mean 0.0, stddev 0.14 ^oC - assuming using ESA CCI ARC data)" % (function)

if (random_noise_sstfnd == 1):
   add_noise(sstfndK_fdata, 0.6, 0.0, nx*ny)
   print "%s Adding random noise to sstfnd (mean 0.0, stddev 0.6 ^oC - assuming using OSTIA data)" % (function)

if (random_noise_pco2 == 1):
   print "/n%s Shape of pco2 data",pco2_sw_fdata.shape
   add_noise(pco2_sw_fdata, 2.0, 0.0, nx*ny)
   print "%s Adding random noise to pco2/fco2 data (mean 0.0, stddev 2.0 uatm - assuming using SOCAT flag A and flag B data)" % (function)

#Ians rain noise test
#if (random_noise_rain == 1):
# rain_data = reshape(rain_data,rain_data.shape[0]*rain_data.shape[1])
# rain_err_data = reshape(rain_err_data,rain_err_data.shape[0]*rain_err_data.shape[1])
# add_noise(rain_data, rain_err_data, 0.0, rain_nx*rain_ny)
# print "%s Adding random noise to rain data using variance field" % (function)

#Ians rain noise test end

 # bias values to be added here
if (bias_windu10 == 1):
   add_bias(windu10_fdata, bias_windu10_value, nx*ny)
    # makes no sense to add bias to second and third order moments, as any bias in the system 
    # would only impact on the mean (which is the 1st order moment)
   print "%s Adding bias to windu10_mean (not to second and third order moments) (value %lf ms^-1)" % (function, bias_windu10_value)

if (bias_sstskin == 1):
   add_bias(sstskinK_fdata, bias_sstskin_value, nx*ny)
   print "%s Adding bias noise to sstskin (value %lf ^oC)" % (function, bias_sstskin_value)

if (bias_sstfnd == 1):
   add_bias(sstfndK_fdata, bias_sstfnd_value, nx*ny)
   print "%s Adding bias noise to sstfnd (value %lf ^oC)" % (function, bias_sstfnd_value)

if (bias_pco2 == 1):
   add_bias(pco2_sw_fdata, bias_pco2_value, nx*ny)
   print "%s Adding bias noise to pco2w (value %lf uatm)" % (function, bias_pco2_value)


 # bias based on change in sstskin due to rain
if (bias_sstskin_due_rain == 1):
   add_sst_rain_bias(sstskinK_fdata, bias_sstskin_due_rain_value, bias_sstskin_due_rain_intensity, rain_fdata, bias_sstskin_due_rain_wind, windu10_fdata, nx*ny)


 # quality filtering of wind and Hs data
for i in arange(nx * ny):
   if (windu10_fdata[i] != missing_value):
       # valid range taken from GlobWave Product User Guide Phase 3 (PUG3) doucment
      if (windu10_fdata[i] > 50.0) or (windu10_fdata[i] < 0.0):
         windu10_fdata[i] = missing_value
   if (sig_wv_ht_fdata[i] != missing_value):
       # valid range taken from GlobWave Product User Guide Phase 3 (PUG3) doucment
      if (sig_wv_ht_fdata[i] > 20.0) or (sig_wv_ht_fdata[i] < 0.0):
         sig_wv_ht_fdata[i] = missing_value
   #if (sigma0_fdata[i] != missing_value):
   #   if (sigma0_fdata[i] > ??.0):
    #     sigma0_fdata[i] = missing_value



if (pco2_data_selection != 0 and VERIFICATION_RUNS != True):
   # signifies that we're using SOCAT data or in-situ data
   #print "%s Using the SOCAT data " % (function)
   #if sstpco2_fdata[i] != missing_value:
   for i in arange(nx * ny):
      if (isnan(sstpco2_fdata[i]) != True and sstpco2_fdata[i] > 0.0 ): 
         if sstpco2_fdata[i] > 260:
           sstpco2_fdata[i] = sstpco2_fdata[i] - 273.15#IGA - If statement added because in-situ SST data may not be in K!
      
         if sstpco2_fdata[i] > 30.5 or sstpco2_fdata[i] < -1.8:
            sstpco2_fdata[i] != missing_value
      else:
         sstpco2_fdata[i] = missing_value  

 # quality control/contrain all SST data
 # check all SST data are within -1.8 - 30.5^oC (or 271.35 - 303.65K)
for i in arange(nx * ny):
   if ((sstskinK_fdata[i] != missing_value) and (sstskinC_fdata[i] != missing_value) and (sstfndK_fdata[i] != missing_value) and (sstfndC_fdata[i] != missing_value)):
      if ( (sstskinC_fdata[i] > 30.5) or (sstfndC_fdata[i] > 30.5) or (sstskinC_fdata[i] < -1.8) or (sstfndC_fdata[i] < -1.8)):
         sstfndK_fdata[i] = missing_value
         sstskinK_fdata[i] = missing_value
         sstfndC_fdata[i] = missing_value
         sstskinC_fdata[i] = missing_value

    # ensure that the different SST data all cover the same spatial regions
    # convert all missing values into standard value, rather than variations that seem to exist in some of these data
for i in arange(nx * ny):
   if ((sstfndK_fdata[i] == sstfndK_fill_value) or (sstskinK_fdata[i] == sstskinK_fill_value)):
      sstfndK_fdata[i] = missing_value
      sstskinK_fdata[i] = missing_value  

 # converting pressure data from Pascals to millibar
if TAKAHASHI_DRIVER != True:
   for i in arange(nx * ny):
      if (pres_fdata[i] != missing_value):
         pres_fdata[i] = pres_fdata[i] * 0.01

 #
 # calculations for components of the flux calculation
 #
 
DeltaT_fdata = array([missing_value] * nx*ny)

if flux_model == 1:
   print "%s Using the RAPID model (from Woolf et al., 2012)" % (function)
elif flux_model == 2:
   print "%s Using the EQUILIBRIUM model (from Woolf et al., 2012)" % (function)
   for i in arange(nx * ny):
      if ( (sstfndC_fdata[i] != missing_value) & (sstskinC_fdata[i] != missing_value) ):
         DeltaT_fdata[i] = sstfndC_fdata[i] - sstskinC_fdata[i]
      else:
         DeltaT_fdata[i] = missing_value
elif flux_model == 3:
    print "%s Using the BULK model (ie Flux = k*S*delta_pCO2)" % (function)
    print "%s Note: this assumes that the SSTskin dataset is the only temperature dataset, and so this is what will be used to calculate k and the solubility" % (function)
else:
   print "\n%s flux_model from configuration not recognised, exitiing." % (function)
   sys.exit(1)

if (saline_skin_value != 0.0):
   print "%s Using the saline skin model (%lf psu added to skin salinities)" % (function, saline_skin_value)


# calculating the schmidt number at the skin and fnd
scskin_fdata = schmidt(sstskinC_fdata, nx, ny, GAS)
scfnd_fdata = schmidt(sstfndC_fdata, nx, ny, GAS)

 # calculating the skin solubility, using skin sst and salinity
solskin_fdata = solubility(sstskinK_fdata, salskin_fdata, DeltaT_fdata, nx, ny, True)

 # calculating the interfacial solubility
solfnd_fdata = solubility(sstfndK_fdata, sal_fdata, DeltaT_fdata, nx, ny, flux_model)


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
dpCO2_diff_fdata = array([missing_value] * nx*ny)
b11_fdata = array([missing_value] * nx*ny)
d12_fdata = array([missing_value] * nx*ny)

 # always using XCO2 from 2000 for SOCAT data
if pco2_data_selection != 0 and pco2_data_selection !=3:#IGA added and to account for in-situ data
   pco2_increment_air = (year - 2000.0) * 1.5
   print "%s year: %d fCO2/pCO2_increment_air: %lf (uatm)" % (function, year, pco2_increment_air)
else:
   pco2_increment_air = pco2_increment
   print "%s year: %d fCO2/pCO2_increment_air: %lf (uatm)" % (function, year, pco2_increment_air)

for i in arange(nx * ny):
   # if pco2_sw_fdata[i] != missing_value:
   #    print 'salskin_fdata[i],sstskinK_fdata[i],pres_fdata[i],vco2_air_fdata[i],sstfndK_fdata[i],sstpco2_fdata[i],pco2_sw_fdata[i],sstskinK_fdata[i]',salskin_fdata[i],sstskinK_fdata[i],pres_fdata[i],vco2_air_fdata[i],sstfndK_fdata[i],sstpco2_fdata[i],pco2_sw_fdata[i],sstskinK_fdata[i]
   if ( (salskin_fdata[i] != missing_value) and (sstskinK_fdata[i] != missing_value) and (pres_fdata[i] != missing_value) and (vco2_air_fdata[i] != missing_value) and (sstfndK_fdata[i] != missing_value) and (sstpco2_fdata[i] != missing_value) and (pco2_sw_fdata[i] != missing_value) and (sstskinK_fdata[i] !=0.0) ):
      
      #print "1sstskinK_fdata: (%d,%d) %d %f log:%f\n" %(nx, ny,i,sstskinK_fdata[i]/100.0,log(sstskinK_fdata[i]/100.0))
      pH2O_fdata[i] = 1013.25 * exp(24.4543 - (67.4509 * (100.0/sstskinK_fdata[i])) - (4.8489 * log(sstskinK_fdata[i]/100.0)) - 0.000544 * salskin_fdata[i])
      #print "2sstskinK_fdata: %d %f log:%f\n" %(i,sstskinK_fdata[i]/100.0,log(sstskinK_fdata[i]/100.0))
      
       # this may be needed when using SMOS salinity data
       # To-DO: awaiting info from Lonneke and David before implementing fully
      pCO2_salinity_term = 0
      #dSalinity = salinity_rmse
      #if (salinity_option == 1):                
      #  # need to determine dS/S using original salinity (its been modified above to add the salinity_rmse)
      #  # so (sal_fdata[i] - salinity_rmse) ensures that we are dealing with the original value of salinity
        #  # using Ys=1 as a global correction following Sarmiento and Gruber, 2006)
     #pCO2_salinity_term = 1.0*(dSalinity/(sal_fdata[i] - salinity_rmse) )
      #else:
      # pCO2_salinity_term = 0.0

       # correction to different years, correction is data and year specific.
       # note for 2010, correction for SOCAT isn't strictly required. However the contents of the exponential will collapse
       # to 1 (with some rounding error expected), so effectively no correction will be applied
      if GAS == 'CO2' and pco2_data_selection != 3:
          pco2_sw_cor_fdata[i] = pco2_increment + (pco2_sw_fdata[i] *exp( (0.0423*(sstfndC_fdata[i] - sstpco2_fdata[i])) - (0.0000435*((sstfndC_fdata[i]*sstfndC_fdata[i]) - (sstpco2_fdata[i]*sstpco2_fdata[i]) )) + pCO2_salinity_term) )
      else:
          pco2_sw_cor_fdata[i] = pco2_sw_fdata[i]
       
      # following disables temp correction when using SOCAT data
      #if pco2_data_selection != 0:
          # check added for testing; only valid for 2010 data, otherwise test is meaningless
          # e.g. (corrected value:15.115202, uncorrected value:350.351990); exponential should collapse to zero

          # if year == 2010 and pco2_sw_cor_fdata[i] != pco2_sw_fdata[i]:
      #         print "%s Using SOCAT data, seawater correction miscalculation, corrected version is not identical to uncorrected version (corrected value:%lf, uncorrected value:%lf, SSTfnd %lf, SSTco2 %lf); exponential should collapse to zero (pco2_increment %lf)" % (function, pco2_sw_cor_fdata[i],pco2_sw_fdata[i], sstfndC_fdata[i], sstpco2_fdata[i], pco2_increment)
      #         pco2_sw_cor_fdata[i] = pco2_sw_fdata[i]
      
       # vco2 in ppm * 1000000 = atm
       # result /1000 to then convert from atm to uatm
       # hence * 0.001 factor
      if GAS == 'CO2' and ATMGAS == 'V':
          pco2_air_cor_fdata[i] = (vco2_air_fdata[i] * 1e-6 * (pres_fdata[i] - pH2O_fdata[i]) / (1e-6 * 1013.25)) + (pco2_increment_air)
      else:
          pco2_air_cor_fdata[i] = vco2_air_fdata[i]

       #pco2_data_selection ==2 signifies SOCAT fCO2 data, so converting pCO2_air_cor_fdata to fCO2_air_cor_fdata      
      if pco2_data_selection == 2:
         # conversion of pCO2 to fCO2 from McGillis and Wanninkhof 2006, Marine chemistry with correction from Weiss 1974 (as the equation in 2006 paper has a set of brackets missing)
         b11_fdata[i] = -1636.75 + (12.0408*sstskinK_fdata[i]) - (0.0327957*sstskinK_fdata[i]*sstskinK_fdata[i]) + (3.16528e-5 * sstskinK_fdata[i]*sstskinK_fdata[i]*sstskinK_fdata[i])
         d12_fdata[i] = 57.7 - (0.118*sstskinK_fdata[i])
          # gas constant
         R = 82.0578 # in [cm^3 atm/(mol K)]
          # 1/0.987 = 1.0131712 - conversion between bar and atm, so *1013.25 is the conversion from millibar to atm.
          # the combination of the B11 and d12 terms are in cm^3/mol and so these cancel with the P/RT term (in mol/cm^3) so the whole of the exp term is dimensionless
         pco2_air_cor_fdata[i] = pco2_air_cor_fdata[i] * exp((b11_fdata[i] + (2*d12_fdata[i]) ) * 1e-6 * ((pres_fdata[i] * 1013.25)/(R*sstskinK_fdata[i]) ))
      
      #((vco2_air_fdata[i] * (pres_fdata[i] - pH2O_fdata[i]))*0.001) + (pco2_increment)
      
       #>>> V=373.769989/1e6
       #>>> pco2a = 374.220001*1e-6*1013.25

      if TAKAHASHI_DRIVER:           
         if pco2_air_fdata[i] != missing_value:
         
            pCO2a_diff_fdata[i] = pco2_air_cor_fdata[i] - pco2_air_fdata[i]
            
            #dpCO2_diff_fdata[i] = (pco2_sw_cor_fdata[i] - pco2_air_cor_fdata[i]) - (pco2_sw_cor_fdata[i] - pco2_air_fdata[i])
            
            pH2O_takahashi_fdata[i] = pres_fdata[i] -  (pco2_air_fdata[i] *1e-6 * 1013.25) / (vco2_air_fdata[i] * 1e-6)
            #-( 1013.25*((pco2_air_fdata[i]*1000000.0)/(vco2_air_fdata[i]/1000000.0)) - pres_fdata[i])
                
            humidity_fdata[i] = pH2O_takahashi_fdata[i]/pH2O_fdata[i]       
            
            pH2O_diff_fdata[i] = ((humidity_fdata[i])-1.0) * 100.0
            
         else:
            pH2O_takahashi_fdata[i] = missing_value
            humidity_fdata[i] = missing_value
            pH2O_diff_fdata[i] = missing_value
              # using Takahashi pCO2a to actual flux calculation
             #pco2_air_cor_fdata[i] = pco2_air_fdata[i]
          
   else:
      pH2O_fdata[i] = missing_value
      pco2_air_cor_fdata[i] = missing_value
      pco2_sw_cor_fdata[i] = missing_value

 # gas transfer velocity k parameterisation choice
k_fdata = array([missing_value] * nx * ny)


 # calculating OceanFlux GHG k parameterisation
 # default values for kb and kd
kb_standard_name="bubble_gas_transfer_velocity_of_carbon_dioxide"
kb_long_name="bubble mediated component of gas transfer velocity of carbon dioxide"
kd_standard_name="direct_gas_transfer_velocity_of_carbon_dioxide_from_backscatter"
kd_long_name="direct component of gas transfer velocity of carbon dioxide derived from radar backscatter"
kt_standard_name = "total_gas_transfer_velocity_of_carbon_dioxide"
kt_long_name="total (direct_from_backscatter plus bubble mediated) component of gas transfer velocity of carbon dioxide"

 # actual calculation
 # run by default as the field always appears in the output netcdf
(kd_fdata, kb_fdata) = OceanFluxGHG_k(sigma0_fdata, sig_wv_ht_fdata, windu10_fdata, windu10_moment2_fdata, sstskinC_fdata, nx, ny)

 # clculate the total kt
kt_fdata = OceanFluxGHG_kt(kd_fdata, kb_fdata, nx, ny, kb_weighting, kd_weighting)

 # intialisation of krain array
krain_fdata = array([missing_value] * nx*ny)

 # user selected k parameterisation for the flux calculation
k_standard_name="unset" # for netcdf metadata
k_long_name="unset" # or netcdf metadata

if (k_parameterisation == 1):
   k_standard_name="gas_transfer_velocity_of_carbon_dioxide" 
   k_long_name="Ho et al., 2006 (H06) gas transfer velocity"
   print "%s Using the Ho et al., 2006 (H06) k parameterisation (default option)" % (function)
 # determine the Ho et al 2006 k relationship
   for i in arange(nx * ny):   
      k_fdata[i] = missing_value
      if ((windu10_fdata[i] != missing_value) and (windu10_moment2_fdata[i] != missing_value) and (windu10_moment3_fdata[i] != missing_value) and (scskin_fdata[i] != missing_value) and (scskin_fdata[i] > 0.0) ):
         k_fdata[i] = 0.266 * windu10_moment2_fdata[i]
         k_fdata[i] = k_fdata[i] * sqrt(600.0/scskin_fdata[i])
                
         k_fdata[i] = k_fdata[i]#/36.0 # conversion from cm/h to 10^-4 m/s (100/3600) = 1/36
      else:
         k_fdata[i] = missing_value
    # resetting OceanFlux bubble and direct components to missing_values as not needed
   kb_fdata[:] = missing_value
   kd_fdata[:] = missing_value   
elif (k_parameterisation == 2):
   k_standard_name="gas_transfer_velocity_of_carbon_dioxide" 
   k_long_name="Nightingale et al., 2000 (N00) gas transfer velocity"
   print "%s Using the Nightingale et al., 2000 (N00) k parameterisation" % (function)
    # determine the Nightingale et al 2000 k relationship
   for i in arange(nx * ny):   
      k_fdata[i] = missing_value
      if ( (windu10_fdata[i] != missing_value) and (windu10_moment2_fdata[i] != missing_value) and (windu10_moment3_fdata[i] != missing_value) and (scskin_fdata[i] != missing_value) and (scskin_fdata[i] > 0.0) ):
         k_fdata[i] = 0.222 * windu10_moment2_fdata[i] + (0.333 * windu10_fdata[i])  
         k_fdata[i] = k_fdata[i] * sqrt(600.0/scskin_fdata[i])
          
         k_fdata[i] = k_fdata[i]# /36.0 # conversion from cm/h to 10^-4 m/s (100/3600) = 1/36
      else:
         k_fdata[i] = missing_value
    # resetting OceanFlux bubble and direct components to missing_values as not needed
   kb_fdata[:] = missing_value
   kd_fdata[:] = missing_value   
elif (k_parameterisation == 3):
    # using OceanFlux GHG kt approach
   k_fdata = kt_fdata
   k_standard_name = kt_standard_name
   k_long_name = kt_long_name
   print "%s Using the OceanFluxGHG kt parameterisation (kt = kd_backscatter + kb)" % (function)
elif (k_parameterisation == 4):
   k_standard_name="gas_transfer_velocity_of_carbon_dioxide" 
   k_long_name="Wanninkhof, 1992 (W92) gas transfer velocity"
   print "%s Using the Wanninkhof 1992 (W92) k parameterisation" % (function)
    # determine the Wanninkhof 1992 k relationship
   for i in arange(nx * ny):   
      k_fdata[i] = missing_value
      if ( (windu10_fdata[i] != missing_value) and (windu10_moment2_fdata[i] != missing_value) and (windu10_moment3_fdata[i] != missing_value) and (scskin_fdata[i] != missing_value) and (scskin_fdata[i] > 0.0) ):
         k_fdata[i] = 0.31 * windu10_moment2_fdata[i]    
         k_fdata[i] = k_fdata[i] * sqrt(660.0/scskin_fdata[i])
          
         k_fdata[i] = k_fdata[i]# /36.0 # conversion from cm/h to 10^-4 m/s (100/3600) = 1/36
      else:
         k_fdata[i] = missing_value
    # resetting OceanFlux bubble and direct components to missing_values as not needed
   kb_fdata[:] = missing_value
   kd_fdata[:] = missing_value   
elif (k_parameterisation == 5):
   k_standard_name="gas_transfer_velocity_of_carbon_dioxide" 
   k_long_name="McGillis et al., 2001 (M01) gas transfer velocity"
   print "%s Using the McGillis et al., 2001 (M01) k parameterisation" % (function)
    # determine the Wanninkhof and McGillis 1999 k relationship
   for i in arange(nx * ny):   
      k_fdata[i] = missing_value
      if ( (windu10_fdata[i] != missing_value) and (windu10_moment2_fdata[i] != missing_value) and (windu10_moment3_fdata[i] != missing_value) and (scskin_fdata[i] != missing_value) and (scskin_fdata[i] > 0.0) ):
         k_fdata[i] = 3.3 + (0.026 * windu10_moment3_fdata[i])
         k_fdata[i] = k_fdata[i] * sqrt((660.0/scskin_fdata[i]))
              
         k_fdata[i] = k_fdata[i]# /36.0 # conversion from cm/h to 10^-4 m/s (100/3600) = 1/36.0
      else:
         k_fdata[i] = missing_value
    # resetting OceanFlux bubble and direct components to missing_values as not needed
   kb_fdata[:] = missing_value
   kd_fdata[:] = missing_value   
elif (k_parameterisation == 6):
   k_standard_name="gas_transfer_velocity_of_carbon_dioxide_due_to_rain"
   k_long_name="Ho et al., 1997 (H97) gas transfer velocity"
   print "%s Using the Ho et al., 1997 (H97) k parameterisation (rain driven k)" % (function)
    # Ho et al, Tellus, 1997 rain component of k
    # note based on rain rate (Rn), but filtered based on wind component so still looking at the same datapoints
   for i in arange(nx * ny):   
      k_fdata[i] = missing_value
      if ( (windu10_fdata[i] != missing_value) and (windu10_moment2_fdata[i] != missing_value) and (windu10_moment3_fdata[i] != missing_value) and (scskin_fdata[i] != missing_value) and (scskin_fdata[i] > 0.0) and (rain_fdata[i] != missing_value) ):
         k_fdata[i] = 0.929 + (0.679 * rain_fdata[i]) - (0.0015*pow(rain_fdata[i], 2.0))
         k_fdata[i] = k_fdata[i] * sqrt(600.0/scskin_fdata[i])
              
         #k_fdata[i] = k_fdata[i]# /36.0 # conversion from cm/h to 10^-4 m/s (100/3600) = 1/36.0
          #copying the result into the krain array
         krain_fdata[i] = k_fdata[i]
      else:
         k_fdata[i] = missing_value
         krain_fdata[i] = missing_value
    # resetting OceanFlux bubble and direct components to missing_values as not needed
   kb_fdata[:] = missing_value
   kd_fdata[:] = missing_value   
elif (k_parameterisation == 7):
    # direct component of k using OceanFluxGHG radar backscatter parameterisation
   k_fdata = kd_fdata
   k_standard_name = kd_standard_name
   k_long_name = kd_long_name   
    # resetting bubble component to missing_values as not needed   
   kb_fdata[:] = missing_value
   print "%s Using the OceanFluxGHG kd-backscatter (direct component) parameterisation" % (function) 
elif (k_parameterisation == 8):
    # bubble mediated component of k using OceanFluxGHG parameterisation   
   k_fdata = kb_fdata
   k_standard_name = kb_standard_name
   k_long_name = kb_long_name 
    # resetting direct component to missing_values as not needed
   kd_fdata[:] = missing_value
   print "%s Using the OceanFluxGHG kb (bubble mediated) parameterisation" % (function)
elif (k_parameterisation == 9):
   # general form
   # kw = (600/Sc)0.5 [a0 + a1*U + a2*U2 + a3*U3]
   k_standard_name="generic_gas_transfer_velocity_formulation"
   k_long_name="User defined generic formulation of k (a0:%lf, a1:%lf, a2:%lf, a3:%lf)" % (k_generic_a0, k_generic_a1, k_generic_a2, k_generic_a3)
   
   if TAKAHASHI_DRIVER == True:
      k_generic_sc = 660.0   
   
   for i in arange(nx * ny):   
      k_fdata[i] = missing_value
      if ( (windu10_fdata[i] != missing_value) and (windu10_moment2_fdata[i] != missing_value) and (windu10_moment3_fdata[i] != missing_value) and (scskin_fdata[i] != missing_value) and (scskin_fdata[i] > 0.0) ):
          # general form kw = (600/Sc)0.5 [a0 + a1*U + a2*U2 + a3*U3]
         k_fdata[i] = k_generic_a0 + (k_generic_a1 * windu10_fdata[i]) +  (k_generic_a2 * windu10_moment2_fdata[i]) + (k_generic_a3 * windu10_moment3_fdata[i])
         k_fdata[i] = k_fdata[i] * sqrt((k_generic_sc/scskin_fdata[i]))
              
         k_fdata[i] = k_fdata[i]# /36.0 # conversion from cm/h to 10^-4 m/s (100/3600) = 1/36.0
      else:
         k_fdata[i] = missing_value  
   print "%s Using the general cubic form of k with user specified values (k_generic_a0:%lf k_generic_a1:%lf k_generic_a2:%lf k_generic_a3:%lf k_generic_sc:%lf)" % (function, k_generic_a0, k_generic_a1, k_generic_a2, k_generic_a3, k_generic_sc)
   kb_fdata[:] = missing_value
   kd_fdata[:] = missing_value
elif (k_parameterisation == 10):
    # direct component using the Goddijn-Murphy et al., 2012 wind parameterisation for kd
   k_standard_name="direct_gas_transfer_velocity_of_carbon_dioxide_from_wind_speed"
   k_long_name="direct component of gas transfer velocity of carbon dioxide derived from wind speed (Goddijn-Murphy et al., 2012)"

   print "%s Using the Goddijn-Murphy et al., 2012 kd-wind (direct-component) parameterisation" % (function)

   kd_fdata = GM12_kd_wind(windu10_fdata, windu10_moment2_fdata, windu10_moment3_fdata, scskin_fdata, nx, ny)
   
    # direct component of k using OceanFluxGHG kd-wind parameterisation
   k_fdata = kd_fdata
   
    # resetting bubble component to missing_values as not needed   
   kb_fdata[:] = missing_value
elif (k_parameterisation == 11):
    # using OceanFlux GHG kt approach with kd based on the wind parameterisation of Goddijn-Murphy et al., 2012
   k_standard_name = "total_gas_transfer_velocity_of_carbon_dioxide"
   k_long_name="total (direct_from_wind_speed plus bubble mediated) component of gas transfer velocity of carbon dioxide"
   print "%s Using the OceanFluxGHG kt parameterisation (kt = kd_wind_speed + kb)" % (function)
   
   kd_fdata = GM12_kd_wind(windu10_fdata, windu10_moment2_fdata, windu10_moment3_fdata, scskin_fdata, nx, ny)
   kt_fdata = OceanFluxGHG_kt(kd_fdata, kb_fdata, nx, ny, kb_weighting, kd_weighting)
   k_fdata[:] = kt_fdata
elif (k_parameterisation == 0):
    # rain wet deposition
   k_standard_name = "wet_deposition so no k parameterisation"
   k_long_name="wet deposition of DIC by rain so no k parameterisation"
    # setting all k values to 0.0
   k_fdata[:] = 0.0
elif (k_parameterisation == 12):
   k_standard_name="W14" 
   k_long_name="Wanninkhof 2014,  Limnol. Oceanogr.: Methods 12, 2014, 351-362"#IGA
   print "%s Using the Wanninkhof 2014 k parameterisation" % (function)
    # determine the k relationship
   for i in arange(nx * ny):   
      k_fdata[i] = missing_value
      if ( (windu10_fdata[i] != missing_value) and (windu10_moment2_fdata[i] != missing_value) and (windu10_moment3_fdata[i] != missing_value) and (scskin_fdata[i] != missing_value) and (scskin_fdata[i] > 0.0) ):
         k_fdata[i] = 0.251 * windu10_moment2_fdata[i]
         k_fdata[i] = k_fdata[i] * sqrt(660.0/scskin_fdata[i])
          
         k_fdata[i] = k_fdata[i]# /36.0 # conversion from cm/h to 10^-4 m/s (100/3600) = 1/36
      else:
         k_fdata[i] = missing_value
    # resetting OceanFlux bubble and direct components to missing_values as not needed
   kb_fdata[:] = missing_value
   kd_fdata[:] = missing_value

else:
   print "%s Chosen k parameterisation option is not recognised or available, exiting." % (function)
   sys.exit(1)




 # linear additive k for rain case
 # adding rain component to the results from the existing choice of parameterisation
if (k_rain_linear_ho1997 == 1):

    # add ho1997 value to data from chosen k parameterisation
    # note: data are filtered based on windU10 to ensure that they cover the same data as the wind based paramterisations  
   for i in arange(nx * ny):
      if ( (k_fdata[i] != missing_value) and (scskin_fdata[i] != missing_value) and (rain_fdata[i] != missing_value) and (scskin_fdata[i] > 0.0) and (windu10_fdata[i] != missing_value) and (windu10_moment2_fdata[i] != missing_value) and (windu10_moment3_fdata[i] != missing_value)):
         krain_fdata[i] = 0.929 + (0.679 * rain_fdata[i]) - (0.0015*pow(rain_fdata[i], 2.0))
         krain_fdata[i] = krain_fdata[i] * sqrt(600.0/scskin_fdata[i])
          #krain_fdata[i] = krain_fdata[i]# /36.0 # conversion from cm/h to 10^-4 m/s (100/3600) = 1/36.0
         k_fdata[i] = k_fdata[i] + krain_fdata[i]
      else:
         krain_fdata[i] = missing_value 

 # nonlinear rain and wind parameterisation from Harrison et al., JGR 2012, equations 11, 12, 13, 14
 # and air_density (in kg/m3) = air_pressure (in pascals) / R(gas constant for dry air) * T(in K)
 # R gas constant for dry air = 287.058 J/(kg K)
 # typical air pressure is 100 000 pascals
 # typical SST in kelvin in 300 K
 # typical density at sea level = 1.25 kg m^-3
 # need to select non-linear parameterisation and the generic k (k_parameterisation == 9)
 # untested 06/01/2014
if (k_rain_nonlinear_h2012 == 1 and k_parameterisation == 9):
   alpha_r = 0.3677 # from Harrison et al., 2012, equation 12
   beta_r = 0.0
   rho_a = 0.0
   # gas constant
   R = 287.058 # in J/(kg K)
           
   windstress_fdata = array([missing_value] * nx*ny)
   winddrag_fdata = array([missing_value] * nx*ny)   
    # need to calculate wind drag and then the wind stress
   for i in arange(nx * ny):
      if ( (windu10_fdata[i] != missing_value) and (windu10_moment2_fdata[i] != missing_value) ):
       # wind drag
       # Yelland and Taylor 1996
         if ( (windu10_fdata[i] >= 6.0) and (windu10_fdata[i] < 26.0) ):
            winddrag_fdata[i] = 0.60 + (0.070 * windu10_fdata[i])
         elif ( (windu10_fdata[i] < 6.0) and (windu10_fdata[i] > 0.0) ):
            winddrag_fdata[i] = 0.29 + (3.1 / windu10_fdata[i]) + ( 7.7 /(windu10_moment2_fdata[i]) )
         else:
            winddrag_fdata[i] = -999.0
         
          # wind stress
         if (winddrag_fdata[i] == -999.0 or windu10_moment2_fdata[i] == -999.0):
            windstress_fdata[i] = -999.0
         elif (winddrag_fdata[i] != 0.0 and windu10_moment2_fdata[i] != 0.0):
            windstress_fdata[i] = sqrt( (winddrag_fdata[i]/1000.0) * windu10_moment2_fdata[i] )
            #print "CD:%lf U^2:%lf U10:%lf stress: %lf" % (winddrag_fdata[i], windu10_moment2_fdata[i], windu10_fdata[i], windstress_fdata[i])
         else:
            windstress_fdata[i] = 0.0
      else:
         windstress_fdata[i] = missing_value
         winddrag_fdata[i] = missing_value
   
   for i in arange(nx * ny):
      if ( (k_fdata[i] != missing_value) and (scskin_fdata[i] != missing_value) and (rain_fdata[i] != missing_value) and (scskin_fdata[i] > 0.0) and (windu10_fdata[i] != missing_value) and (windu10_moment2_fdata[i] != missing_value) and (windu10_moment3_fdata[i] != missing_value) and (pres_fdata[i] != missing_value) and (sstskinK_fdata[i] != missing_value) and (windstress_fdata[i] != missing_value)):
          #guts in here
         # pres_fdata are in mb, need it in Pascals, so *100 to get Pascals (ie convert from mb to P)
         rho_a = (pres_fdata[i]*100.0) /(R*sstskinK_fdata[i])    
         beta_r = (0.0112*rain_fdata[i]) / (rho_a*pow(windstress_fdata[i],3)) 
         #print "beta_r:%lf stress^3:%lf rho_a:%lf pres:%lf R:%lf sstskinK:%lf" % (beta_r, pow(windstress_fdata[i],3), rho_a, pres_fdata[i], R, sstskinK_fdata[i])
         krain_fdata[i] = (1-exp(-(alpha_r*beta_r))) * (63.02 * pow(0.0112*rain_fdata[i], 0.6242))
         krain_fdata[i] = krain_fdata[i] * sqrt(600.0/scskin_fdata[i])
         k_fdata[i] = k_fdata[i] + krain_fdata[i]
      else:
         krain_fdata[i] = missing_value
         k_fdata[i] = missing_value


 # ability to investigate bias on k due to surface biology/slicks
 # assumes that bias values are realistic and that they won't cause the k_fdata to become unrealistic
 # also assumes that bio_fdata exist (ie it won't if they indicator layers are off)
if (bias_k == 1):
   add_bias_k_biology_wind(k_fdata, bias_k_value, bio_fdata, bias_k_biology_value, windu10_fdata, bias_k_wind_value, nx*ny, bias_k_percent)
   if bias_k_percent==0:
      print "\n%s Adding bias to chosen k (k_fdata) parameterisation data (bias value of %lf ms^-1 added, where biology (bio_fdata) is > %lf mg m^-3 and wind speed (windu10_fdata) is < %lf m s^-1)" % (function, bias_k_value, bias_k_biology_value, bias_k_wind_value)
   else:
      print "\n%s Adding bias to chosen k (k_fdata) parameterisation data (bias value of - %lf percent being used, where biology (bio_fdata) is > %lf mg m^-3 and wind speed (windu10_fdata) is < %lf m s^-1)" % (function, bias_k_value, bias_k_biology_value, bias_k_wind_value)
      

 #
 # actual flux calculation
 #


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
FH06_fdata = array([missing_value] * nx*ny)
concw_fdata = array([missing_value] * nx*ny)
conca_fdata = array([missing_value] * nx*ny)
FKo07_fdata = array([missing_value] * nx*ny)

 # creating a solubility dataset using zero salinity (ie distilled water)
solskinDistilWater_fdata = array([missing_value] * nx*ny)
if rain_wet_deposition:
   salDistil_fdata = array([missing_value] * nx*ny)
   for i in arange(nx * ny):
      if (sal_fdata[i] != missing_value):
         salDistil_fdata[i] = 0.0
      else:
         salDistil_fdata[i] = missing_value   
    # calculate the solubility of distilled water, using the zero salinity dataset
   solskinDistilWater_fdata = solubility(sstskinK_fdata, salDistil_fdata, DeltaT_fdata, nx, ny, True)

if ((kb_asymmetry != 1.0) and (k_parameterisation == 3)):
   print "%s kb asymetry has been enabled (kb_asymmetry:%lf and k_parameterisation:%d)" % (function, kb_asymmetry, k_parameterisation)
   
for i in arange(nx * ny):   
   FH06_fdata[i] = missing_value
  
   if ( (solfnd_fdata[i] != missing_value) and (k_fdata[i] != missing_value) and (pco2_sw_cor_fdata[i] != missing_value) and (pco2_air_cor_fdata[i] != missing_value) ):
       # original
      #flux_H06[i] = (k_H06[i] * 36.0 * 24.0 * pow(10.0,-2)) * (solubility_data[i]* (12.0/1000.0)) * DpCO2_cor_data[i]
       # expanded
      
      # mass boundary layer concentration (ie concentration in the water)
      concw_fdata[i] = ((solfnd_fdata[i] * conc_factor) * pco2_sw_cor_fdata[i])
      
       # interfacial concentration (ie at the interface between the ocean and the atmosphere)
      conca_fdata[i] = ((solskin_fdata[i] * conc_factor) * pco2_air_cor_fdata[i])
      
       #flux calculation
      if ((kb_asymmetry != 1.0) and (k_parameterisation == 3)):
         kd_component = kd_fdata[i] * k_factor
         kb_component = kb_fdata[i] * k_factor
         FH06_fdata[i] = (kd_component * (concw_fdata[i] - conca_fdata[i])) + (kb_component * (concw_fdata[i] - (kb_asymmetry *conca_fdata[i]) ) )
      else:
         FH06_fdata[i] = (k_fdata[i] * k_factor) * (concw_fdata[i] - conca_fdata[i])
            
       # using simplified flux calculation with no separation of temperature between airside and waterside CO2
       # assumes that the skin temperature dataset is the only temperature dataset
      if flux_model == 3:         
         FH06_fdata[i] = (k_fdata[i] * k_factor) * conc_factor * solskin_fdata[i] * (pco2_sw_cor_fdata[i] - pco2_air_cor_fdata[i])#IGA added conc_factor
          # using Rik W's equation 6 from his new paper
         #FH06_fdata[i] = (7.7e-4*(12.0108/365.0)) * windu10_moment2_fdata[i] * (pco2_sw_cor_fdata[i] - pco2_air_cor_fdata[i])
      
      
       # calculating and adding in flux component for wet deposition due to rain
      if rain_wet_deposition:
          # relationship from Komori et al.m 2007
          # need solubility of CO2 in fresh water
          # solubility calculated using sstkin
          # 24.0/1000.0 = conversion from mm hr^-1 to m day^-1
          # flux is then in g C m^-2 day^-1
          # flux is always negative, ie going into the ocean
         FKo07_fdata[i] = -(rain_fdata[i] * (24.0/1000.0)) * (conc_factor * solskinDistilWater_fdata[i]) * pco2_air_cor_fdata[i]
         if FH06_fdata[i] != missing_value:
            FH06_fdata[i] += FKo07_fdata[i]
         else:
            FH06_fdata[i] = FKo07_fdata[i]
      else:
         FKo07_fdata[i] = missing_value
   else:
      FH06_fdata[i] = missing_value
      concw_fdata[i] = missing_value
      conca_fdata[i] = missing_value
      FKo07_fdata[i] = missing_value
      
#Adding verification data outputs at same units as T09
if TAKAHASHI_DRIVER == True:
  solskin_takadata = array([missing_value] * nx * ny)
  FH06_takadata = array([missing_value] * nx * ny)
  for i in arange(nx * ny):
   if solskin_fdata[i] != missing_value:
      solskin_takadata[i] = solskin_fdata[i]*1000 #from 'mol kg-1 atm-1' to 'mmol kg-1 atm-1'
   if FH06_fdata[i] != missing_value:
      FH06_takadata[i] = FH06_fdata[i]/30.5 #from flux per day to flux per month
else:
  solskin_takadata = array([missing_value] * nx * ny)
  FH06_takadata = array([missing_value] * nx * ny)
  #
  # quality control of datasets, following TS (OceanFluxGHG_TS_D2-9_v1.8-signed.pdf) table 7
  #

 # array to log the total number of quality violations per pixel location  
failed_quality_fdata = array([missing_value] * nx*ny)

  # checking OF flux dataset
check_output_dataset(FH06_fdata, 'OF', failed_quality_fdata, nx, ny, OF_min, OF_max, sstskinC_fdata, sstfndC_fdata, windu10_fdata, sig_wv_ht_fdata, solfnd_fdata, solskin_fdata, k_fdata, concw_fdata, conca_fdata, sal_fdata)

  # checking oceanflux k total
check_output_dataset(kt_fdata, 'OK1', failed_quality_fdata, nx, ny, OK1_min, OK1_max, sstskinC_fdata, sstfndC_fdata, windu10_fdata, sig_wv_ht_fdata, solfnd_fdata, solskin_fdata, k_fdata, concw_fdata, conca_fdata, sal_fdata)

  # check user chosen k data
check_output_dataset(k_fdata, 'OK3', failed_quality_fdata, nx, ny, OK3_min, OK3_max, sstskinC_fdata, sstfndC_fdata, windu10_fdata, sig_wv_ht_fdata, solfnd_fdata, solskin_fdata, k_fdata, concw_fdata, conca_fdata, sal_fdata)

 # check limits of the direct kd
check_output_dataset(kd_fdata, 'OKD', failed_quality_fdata, nx, ny, OKD_min, OKD_max, sstskinC_fdata, sstfndC_fdata, windu10_fdata, sig_wv_ht_fdata, solfnd_fdata, solskin_fdata, k_fdata, concw_fdata, conca_fdata, sal_fdata)

 # check limits of the bubble mediated kb
check_output_dataset(kb_fdata, 'OKB1', failed_quality_fdata, nx, ny, OKB1_min, OKB1_max, sstskinC_fdata, sstfndC_fdata, windu10_fdata, sig_wv_ht_fdata, solfnd_fdata, solskin_fdata, k_fdata, concw_fdata, conca_fdata, sal_fdata)

 # checking salinity data
check_output_dataset(sal_fdata, 'OKS1', failed_quality_fdata, nx, ny, OKS1_min, OKS1_max, sstskinC_fdata, sstfndC_fdata, windu10_fdata, sig_wv_ht_fdata, solfnd_fdata, solskin_fdata, k_fdata, concw_fdata, conca_fdata, sal_fdata)

 # checking skin temperature
check_output_dataset(sstskinC_fdata, 'OKT1', failed_quality_fdata, nx, ny, OKT1_min, OKT1_max, sstskinC_fdata, sstfndC_fdata, windu10_fdata, sig_wv_ht_fdata, solfnd_fdata, solskin_fdata, k_fdata, concw_fdata, conca_fdata, sal_fdata)

 # checking interfacial concentration (ie water side)
check_output_dataset(concw_fdata, 'OSC1', failed_quality_fdata, nx, ny, OSC1_min, OSC1_max, sstskinC_fdata, sstfndC_fdata, windu10_fdata, sig_wv_ht_fdata, solfnd_fdata, solskin_fdata, k_fdata, concw_fdata, conca_fdata, sal_fdata)

 # checking mass boundary layer concentration (ie atmospheric side)
check_output_dataset(conca_fdata, 'OIC1', failed_quality_fdata, nx, ny, OIC1_min, OIC1_max, sstskinC_fdata, sstfndC_fdata, windu10_fdata, sig_wv_ht_fdata, solfnd_fdata, solskin_fdata, k_fdata, concw_fdata, conca_fdata, sal_fdata)


 # final checks - these are not stated within TS but adding for completeness
 # calculate takahashi style DpCO2 and check range
for i in arange(nx * ny):
   if ( (pco2_sw_cor_fdata[i] != missing_value) and (pco2_air_cor_fdata[i] != missing_value) ):
      dpco2_cor_fdata[i] = (pco2_sw_cor_fdata[i] - pco2_air_cor_fdata[i]) 
   else:
      dpco2_cor_fdata[i] = missing_value

for i in arange(nx * ny):
   if ( (concw_fdata[i] != missing_value) and (conca_fdata[i] != missing_value) ):
      dpconc_cor_fdata[i] = (concw_fdata[i] - conca_fdata[i]) 
   else:
      dpconc_cor_fdata[i] = missing_value



 # now check the data range      
check_output_dataset(dpco2_cor_fdata, 'DpCO2', failed_quality_fdata, nx, ny, -200, 200, sstskinC_fdata, sstfndC_fdata, windu10_fdata, sig_wv_ht_fdata, solfnd_fdata, solskin_fdata, k_fdata, concw_fdata, conca_fdata, sal_fdata)

 # checking foundation temperature range
check_output_dataset(sstfndC_fdata, 'FT1', failed_quality_fdata, nx, ny, OKT1_min, OKT1_max, sstskinC_fdata, sstfndC_fdata, windu10_fdata, sig_wv_ht_fdata, solfnd_fdata, solskin_fdata, k_fdata, concw_fdata, conca_fdata, sal_fdata)

 # checking the sub skin partial pressure or fugacity of CO2 in seawater
check_output_dataset(conca_fdata, 'OBPC', failed_quality_fdata, nx, ny, OBPC_min, OBPC_max, sstskinC_fdata, sstfndC_fdata, windu10_fdata, sig_wv_ht_fdata, solfnd_fdata, solskin_fdata, k_fdata, concw_fdata, conca_fdata, sal_fdata)



 # whitecapping data using relationship from TS and parameters from table 1 of Goddijn-Murphy et al., 2010, equation r1
whitecap_fdata = array([missing_value] * nx*ny)
for i in arange(nx * ny):
   if (windu10_fdata[i] != missing_value):
      whitecap_fdata[i] = 0.00159 * pow(windu10_fdata[i], 2.7)
   else:
      whitecap_fdata[i] = missing_value


  # pco2_data_selection == 1 signifies that we're using SOCAT data 
  # this dataset is filled with missing_values prior to output, otherwise non-socat data will appear in the SFUG field in the netcdf
  # which will confuse the user
  # enabled for TAKAHASHI_DRIVER to enable checking
if TAKAHASHI_DRIVER != True:
   if (pco2_data_selection == 0):
      pco2_sw_fdata = array([missing_value] * nx*ny)
      pco2_sw_stddev_fdata = array([missing_value] * nx*ny)


 # fill in the missing_values in the _N and _stddev datasets
for i in arange(nx * ny):
   if (windu10_fdata[i] == missing_value):
       windu10_count_fdata[i] = missing_value
       windu10_stddev_fdata[i] = missing_value
   if (sstskinK_fdata[i] == missing_value):
       sstskinK_count_fdata[i] = missing_value
       sstskinK_stddev_fdata[i] = missing_value
   if (sstfndK_fdata[i] == missing_value):
       sstfndK_count_fdata[i] = missing_value
       sstfndK_stddev_fdata[i] = missing_value
   if (sig_wv_ht_fdata[i] == missing_value):
       sig_wv_ht_count_fdata[i] = missing_value
       sig_wv_ht_stddev_fdata[i] = missing_value


  #
  # procesin indictor attribute layers
  #
if process_layers_off == 0:
   lowwind_fdata = array([missing_value] * nx*ny)
   bioclass_fdata = array([missing_value] * nx*ny)
   diurnalw_fdata = array([missing_value] * nx*ny)
  
    # low wind attribute layer
   for i in arange(nx * ny):
      lowwind_fdata[i] = missing_value_int
      if (windu10_fdata[i] != missing_value):
         if (windu10_fdata[i] <= 5.0):
           lowwind_fdata[i] = 1
         elif (windu10_fdata[i] > 5.0):
            lowwind_fdata[i] = 0
      else:
         lowwind_fdata[i] = missing_value_int

    # biological activity layer
    # adding in the > 0.0 condition as preliminary ESA CCI data are missing a load of data attributes
   for i in arange(nx * ny):
      bioclass_fdata[i] = missing_value_int
      if ((bio_fdata[i] != missing_value) and (bio_fdata[i] > 0.0) ):
         if  ( (bio_fdata[i] <= 0.5)):
            bioclass_fdata[i] = 1
         elif (bio_fdata[i] <= 1.0):
            bioclass_fdata[i] = 2
         elif( bio_fdata[i] > 1.0 ):
            bioclass_fdata[i] = 3
      else:
         bioclass_fdata[i] = missing_value_int

     # diurnal warming (sstskin > sstfnd)
   for i in arange(nx * ny):
      diurnalw_fdata[i] = missing_value_int
      if ( (sstskinK_fdata[i] != missing_value) and (sstfndK_fdata[i] != missing_value) ):
         sst_diff = sstskinK_fdata[i] - sstfndK_fdata[i]
         if ( (sstskinK_fdata[i] > sstfndK_fdata[i]) and (sst_diff > 0.5) ):
             # >0.5 difference condition focusses on strong gradients
            diurnalw_fdata[i] = 1
         else:
            diurnalw_fdata[i] = 0
      else:
         diurnalw_fdata[i] = missing_value_int


     # oceanic basins
     # re-assingin values and adding in missing_value_int entries
   for i in arange(nx * ny):
      if (atlantic_ocean_fdata[i] != missing_value):
         if  (atlantic_ocean_fdata[i] == 30.0 ):     
            atlantic_ocean_fdata[i] = 1.0
         else:
            atlantic_ocean_fdata[i] = missing_value
      if (pacific_ocean_fdata[i] != missing_value):
         if  (pacific_ocean_fdata[i] == 70.0 ):  
            pacific_ocean_fdata[i] = 1.0
         else:
            pacific_ocean_fdata[i] = missing_value
      if (southern_ocean_fdata[i] != missing_value):
         if  (southern_ocean_fdata[i] == 90.0 ):     
            southern_ocean_fdata[i] = 1.0
         else:
            southern_ocean_fdata[i] = missing_value
      if (indian_ocean_fdata[i] != missing_value):
         if  (indian_ocean_fdata[i] == 50.0 ):   
            indian_ocean_fdata[i] = 1.0
         else:
            indian_ocean_fdata[i] = missing_value
 
     # longhurst provinces
     # TO-DO: need the full longhurst mask file, currently only 3 provinces
   for i in arange(nx * ny):
      if (longhurst_fdata[i] != missing_value):
         if  ( longhurst_fdata[i] == 0.0 ):
            longhurst_fdata[i] = missing_value_int
         elif (longhurst_fdata[i] == 30.0): #
            longhurst_fdata[i] = 1
      else:       
         longhurst_fdata[i] = missing_value_int





 #   
 # write out results
 #

 # re-shaping datasets to 2d array ready for output
windu10_fdata.shape = (nx, ny)
windu10_stddev_fdata.shape = (nx, ny)
windu10_count_fdata.shape = (nx, ny)

sstskinK_stddev_fdata.shape = (nx, ny)
sstskinK_count_fdata.shape = (nx, ny)

sstfndK_stddev_fdata.shape = (nx, ny)
sstfndK_count_fdata.shape = (nx, ny)

sig_wv_ht_fdata.shape = (nx, ny)
sig_wv_ht_stddev_fdata.shape = (nx, ny)
sig_wv_ht_count_fdata.shape = (nx, ny)

sstskinC_fdata.shape = (nx, ny)
sstfndC_fdata.shape = (nx, ny)
FH06_fdata.shape = (nx, ny)
k_fdata.shape = (nx, ny)
concw_fdata.shape = (nx, ny)
conca_fdata.shape = (nx, ny)
FKo07_fdata.shape = (nx, ny)
scskin_fdata.shape = (nx * ny)
krain_fdata.shape = (nx * ny)
kt_fdata.shape = (nx *ny)
kd_fdata.shape = (nx *ny)
kb_fdata.shape = (nx *ny)
whitecap_fdata.shape = (nx *ny)

 # additional datasets for testing
pH2O_takahashi_fdata.shape = (nx *ny)
humidity_fdata.shape = (nx *ny)
pH2O_diff_fdata.shape = (nx *ny)
pH2O_fdata.shape = (nx *ny)
pCO2a_diff_fdata.shape = (nx *ny)
dpCO2_diff_fdata.shape = (nx *ny)

 # other parameters commented out - keep
solskin_fdata.shape = (nx, ny)
solfnd_fdata.shape = (nx, ny)
solskinDistilWater_fdata.shape = (nx, ny)
dpco2_cor_fdata.shape = (nx, ny)
dpconc_cor_fdata.shape = (nx, ny)
pco2_air_cor_fdata.shape = (nx, ny)
pco2_sw_cor_fdata.shape = (nx, ny)
vco2_air_fdata.shape = (nx, ny)
pco2_air_fdata.shape = (nx, ny)
solskin_takadata.shape = (nx, ny)
FH06_takadata.shape = (nx, ny)
pco2_sw_fdata.shape = (nx, ny)
pco2_sw_stddev_fdata.shape = (nx, ny)

pres_fdata.shape = (nx, ny)
sal_fdata.shape = (nx, ny)
rain_fdata.shape = (nx, ny)

 # process indicator attribute layers
if process_layers_off == 0:
   bio_fdata.shape = (nx, ny)
   lowwind_fdata.shape = (nx, ny)
   diurnalw_fdata.shape = (nx, ny)
   bioclass_fdata.shape = (nx, ny)
   atlantic_ocean_fdata.shape = (nx, ny)
   pacific_ocean_fdata.shape = (nx, ny)
   southern_ocean_fdata.shape = (nx, ny)
   indian_ocean_fdata.shape = (nx, ny)
   susp_particles_fdata.shape = (nx, ny)
   sstgrad_fdata_avg.shape = (nx, ny)
   longhurst_data.shape = (nx, ny)
else:
    # intialising all data arrays to missing_values so that they exist but are empty
   bio_fdata = array([missing_value] * nx*ny)
   lowwind_fdata = array([missing_value] * nx*ny)
   diurnalw_fdata = array([missing_value] * nx*ny)
   bioclass_fdata = array([missing_value] * nx*ny)
   atlantic_ocean_fdata = array([missing_value] * nx*ny)
   pacific_ocean_fdata = array([missing_value] * nx*ny)
   southern_ocean_fdata = array([missing_value] * nx*ny)
   indian_ocean_fdata = array([missing_value] * nx*ny)
   susp_particles_fdata = array([missing_value] * nx*ny)
   sstgrad_fdata_avg = array([missing_value] * nx*ny)
   longhurst_fdata = array([missing_value] * nx*ny)

 # need to always have this process attribute layer as its part of the internal quality check
failed_quality_fdata.shape = (nx, ny)

 # ice data
ice_fdata.shape = (nx, ny)
 # write out the final ouput to netcdf
write_netcdf(FH06_fdata, FKo07_fdata, k_fdata, kt_fdata, kd_fdata, kb_fdata, windu10_fdata, windu10_stddev_fdata, windu10_count_fdata, sig_wv_ht_fdata, sig_wv_ht_stddev_fdata, sig_wv_ht_count_fdata, sstskinC_fdata, sstskinK_stddev_fdata, sstskinK_count_fdata, sstfndC_fdata, sstfndK_stddev_fdata, sstfndK_count_fdata, pco2_air_cor_fdata, pco2_sw_cor_fdata, vco2_air_fdata, pco2_sw_fdata, sal_fdata, pres_fdata, conca_fdata, concw_fdata, rain_fdata, scskin_fdata, krain_fdata, whitecap_fdata, dpco2_cor_fdata, dpconc_cor_fdata, solskin_fdata, solfnd_fdata, solskinDistilWater_fdata, dpCO2_diff_fdata, pCO2a_diff_fdata, pH2O_fdata, pH2O_takahashi_fdata, humidity_fdata, pH2O_diff_fdata, solskin_takadata, FH06_takadata, ice_fdata, lowwind_fdata, diurnalw_fdata, bioclass_fdata, bio_fdata, atlantic_ocean_fdata, pacific_ocean_fdata, southern_ocean_fdata, indian_ocean_fdata, longhurst_fdata, susp_particles_fdata, sstgrad_fdata_avg, failed_quality_fdata, k_standard_name, k_long_name, kb_standard_name, kb_long_name, kd_standard_name, kd_long_name, pco2_air_fdata)


print "%s SUCCESS writing file %s" % (function, outfile)

