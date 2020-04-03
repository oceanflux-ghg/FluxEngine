#!/home/cercache/tools/environments/run_python_scibox_precise.sh
#-*- coding:utf-8 -*-

import sys
import copy
import numpy
import netCDF4
from . import netcdf_helper


def FromFilelist(filelist,output,weighting=None,outputtime=1e9,combiningregions=False):
   """
   Take a list of netCDF files and combine (mean average) all the same variables
   into a new netCDF file given by 'output'.

   Requires a variable named "count_nobs" - this should be updated and changed.

   This may not work as expected if some files have differently named (i.e. different) variables.
   """
   if weighting is not None:
      weighting_string=" %s."%(weighting)
   else:
      weighting_string=""

   dataitem_dims={}
   itemdata={}
   itemvar={}
   sumN=0

   if not isinstance(filelist,list):
      raise Exception("filelist should be a list: in FromFilelist.")

   #if filelist is empty then return as nothing to do
   if len(filelist)==0:
      return False

   #if output is not a string then raise exception
   if not isinstance(output,str):
      raise Exception("output should be a string filename to write output to: in FromFilelist.")

   for filename in filelist:
      #open the file
      fin=netCDF4.Dataset(filename)
      #get a list of all variables
      filevars=list(fin.variables.keys())
      filedims=list(fin.dimensions.keys())
      #read all the dimensions and
      #read all of the data variables
      for variable in filevars:
         if variable in filedims:
            #variable is a dimension - check if it is already in the dimension dictionary
            if variable not in list(dataitem_dims.keys()):
               dataitem_dims[variable]=copy.copy(fin.variables[variable][:])
            else:
               #test if the numpy arrays are the same (i.e. same dimensions)
               if not numpy.array_equal(dataitem_dims[variable],fin.variables[variable][:]) and variable!="time":
                  print("Dimensions are not the same for %s  what does this mean, they should be. "%variable)
                  print(dataitem_dims[variable],fin.variables[variable][:])
                  sys.exit(1)
               elif variable == "time":
                  #time could be a different value so this is OK
                  dataitem_dims['time']=numpy.vstack([dataitem_dims['time'],copy.copy(fin.variables[variable][:])])
         else: # not a dimension
            if variable in list(itemdata.keys()):
               #We set the filled value to 0 so that we can add without worrying about no-data-values
               itemdata[variable]=numpy.ma.vstack((itemdata[variable],fin.variables[variable][:]))
            else:
               #We set the filled value to 0 so that we can add without worrying about no-data-values
               itemdata[variable]=fin.variables[variable][:].astype(numpy.float64)
               itemvar[variable]=[fin.variables[variable].standard_name,fin.variables[variable]]
         #we only want to get the cells with data in once so we choose the count_nobs variable
         #and add onto the previous number to keep a running total
         if variable=="count_nobs":
            sumN+=numpy.where(fin.variables[variable][:].filled(0)!=0,1,0)

   #if combing regions to create the global product then we don't need to do much
   #other than sum up the products
   #as we assume that there are no overlap between regions
   if combiningregions == True:
      for variable in list(itemdata.keys()):
         itemdata[variable]=numpy.ma.sum(itemdata[variable],axis=0)
   #if we are not combining regions - i.e. we are combining cruises, then 
   #there could be overlap so we need to update the parameters appropriately
   else:
      #need to do std first as it depends on other variables
      for variable in list(itemdata.keys()):
         if variable.startswith("std_"):
            #we don't want to average these - we want to update them
            #but we want to update them by taking the stdev of the new data
            postfix=variable.replace("std_","")
            #"weighted" standard deviation (i.e. standard deviation based on these values only)
            itemdata[variable]=numpy.ma.std(itemdata[postfix],axis=0,ddof=1)

      for variable in list(itemdata.keys()):
         if variable.startswith("std_"):
            continue
         elif variable in ["count_nobs","count_ncruise"]: 
            # we don't want to average count_nobs or count_ncruise as these are number of observations
            itemdata[variable]=itemdata[variable].sum(axis=0)
         elif variable.startswith("max_"):
            #we don't want to average these - we want to update them
            itemdata[variable]=numpy.ma.max(itemdata[variable],axis=0)
         elif variable.startswith("min_"):
            #we don't want to average these - we want to update them
            itemdata[variable]=numpy.ma.min(itemdata[variable],axis=0)
         else:
            #we will take the mean of these - this is the weighted average where weighted means
            #we take into account any previous averaging of the data
            itemdata[variable]=numpy.ma.mean(itemdata[variable],axis=0)
            if weighting_string != "":
               itemvar[variable][0]="weighted_"+itemvar[variable][0]

   #Write out the concatenated variables to a new file
   with netCDF4.Dataset(output, 'w', format = 'NETCDF4') as ncfile:
      netcdf_helper.standard_setup_SOCAT(ncfile,timedata=outputtime,londata=dataitem_dims['longitude'],latdata=dataitem_dims['latitude'])
      for variable in list(itemdata.keys()):
         if variable not in filedims:
            tmpvar = ncfile.createVariable(itemvar[variable][0],'f4',('time','latitude','longitude'),
                                           fill_value=netcdf_helper.MISSINGDATAVALUE,zlib=True)
            tmpvar[:] = itemdata[variable][:]
            tmpvar.units = itemvar[variable][1].units
            tmpvar.missing_value = itemvar[variable][1].missing_value
            tmpvar.valid_min = itemvar[variable][1].valid_min
            tmpvar.valid_max = itemvar[variable][1].valid_max
            tmpvar.scale_factor = itemvar[variable][1].scale_factor
            tmpvar.add_offset = itemvar[variable][1].add_offset
            tmpvar.standard_name = itemvar[variable][0]
            if variable not in ["count_nobs","count_ncruise"] \
                        and not variable.startswith("min_") \
                        and not variable.startswith("max_"):
               tmpvar.long_name = itemvar[variable][1].long_name+weighting_string
            else:
               tmpvar.long_name = itemvar[variable][1].long_name
      if "count_ncruise" not in list(itemdata.keys()):
         #Also add a new variable for the number of cruises per cell
         ncruise = ncfile.createVariable("count_ncruise",'f4',('time','latitude','longitude'),fill_value=0,zlib=True)
         ncruise[:] = sumN
         ncruise.units = "count"
         ncruise.missing_value = netcdf_helper.MISSINGDATAVALUE
         ncruise.valid_min = 0
         ncruise.valid_max = 1e10
         ncruise.scale_factor = 1
         ncruise.add_offset = 0
         ncruise.standard_name = "count_ncruise"
         ncruise.long_name = "Number of cruises"
   return True

#This function is for wrapping the per-month global files into a single file
#using the time as a 3rd dimension for the arrays
def IntoTimeDimension(filelist,output):
   #We need to read in all the netcdf variables and create new one(s) which are concatenated
   if not isinstance(filelist,list):
      raise Exception("filelist should be a list: in FromFilelist.")

   #if filelist is empty then return as nothing to do
   if len(filelist)==0:
      return False

   #if output is not a string then raise exception
   if not isinstance(output,str):
      raise Exception("output should be a string filename to write output to: in FromFilelist.")

   dataitem_dims={}
   itemdata={}
   itemvar={}

   for filename in filelist:
      #open the file
      fin=netCDF4.Dataset(filename)
      #get a list of all variables
      filevars=list(fin.variables.keys())
      filedims=list(fin.dimensions.keys())
      #read all the dimensions and
      #read all of the data variables
      for variable in filevars:
         if variable in filedims:
            #variable is a dimension - check if it is already in the dimension dictionary
            if variable not in list(dataitem_dims.keys()):
               dataitem_dims[variable]=copy.copy(fin.variables[variable][:])
            else:
               #test if the numpy arrays are the same (i.e. same dimensions)
               if not numpy.array_equal(dataitem_dims[variable],fin.variables[variable][:]) and variable != "time":
                  print("Dimensions are not the same for %s  what does this mean, they should be. "%variable)
                  print(dataitem_dims[variable],fin.variables[variable][:])
                  sys.exit(1)
               if variable == "time":
                  dataitem_dims['time']=numpy.vstack([dataitem_dims['time'],copy.copy(fin.variables[variable][:])])
         else: # not a dimension
            if variable in list(itemdata.keys()):
               #itemdata[variable]=itemdata[variable]+fin.variables[variable][:].filled(0)
               itemdata[variable]=numpy.ma.vstack((itemdata[variable],fin.variables[variable][:]))
            else:
               itemdata[variable]=fin.variables[variable][:].astype(numpy.float64)
               itemvar[variable]=fin.variables[variable]

   #Write out the concatenated variables to a new file
   with netCDF4.Dataset(output, 'w', format = 'NETCDF4') as ncfile:
      netcdf_helper.standard_setup_SOCAT(ncfile,timedata=dataitem_dims['time'],londata=dataitem_dims['longitude'],latdata=dataitem_dims['latitude'])
      for variable in list(itemdata.keys()):
         if variable not in filedims:
            tmpvar = ncfile.createVariable(variable,'f4',('time','latitude','longitude'),
                                           fill_value=netcdf_helper.MISSINGDATAVALUE,zlib=True)
            tmpvar[:] = itemdata[variable][:]
            tmpvar.units = itemvar[variable].units
            tmpvar.missing_value = itemvar[variable].missing_value
            tmpvar.valid_min = itemvar[variable].valid_min
            tmpvar.valid_max = itemvar[variable].valid_max
            tmpvar.scale_factor = itemvar[variable].scale_factor
            tmpvar.add_offset = itemvar[variable].add_offset
            tmpvar.standard_name = itemvar[variable].standard_name
            tmpvar.long_name = itemvar[variable].long_name

def AddNewVariables(filename,newvars):
   print("AddNewVariables")
   with netCDF4.Dataset(filename, 'a', format = 'NETCDF4') as ncfile:
      fCO2_SST_data = ncfile.createVariable('unweighted_fCO2_SST','f4',('time','latitude','longitude'),
                                            fill_value=netcdf_helper.MISSINGDATAVALUE,zlib=True)
      fCO2_SST_data[:] = newvars['fCO2_SST'][:]
      fCO2_SST_data.units = 'uatm'
      fCO2_SST_data.missing_value = netcdf_helper.MISSINGDATAVALUE
      fCO2_SST_data.valid_min = 0.
      fCO2_SST_data.valid_max = 1e6
      fCO2_SST_data.scale_factor = 1.
      fCO2_SST_data.add_offset = 0.
      fCO2_SST_data.standard_name = "unweighted_fCO2_SST"
      fCO2_SST_data.long_name = "CO2 fugacity using SOCAT methodology (unweighted)"

      fCO2_Tym_data = ncfile.createVariable('unweighted_fCO2_Tym','f4',('time','latitude','longitude'),
                                            fill_value=netcdf_helper.MISSINGDATAVALUE,zlib=True)
      fCO2_Tym_data[:] = newvars['fCO2_Tym'][:]
      fCO2_Tym_data.units = 'uatm'
      fCO2_Tym_data.missing_value = netcdf_helper.MISSINGDATAVALUE
      fCO2_Tym_data.valid_min = 0.
      fCO2_Tym_data.valid_max = 1e6
      fCO2_Tym_data.scale_factor = 1.
      fCO2_Tym_data.add_offset = 0.
      fCO2_Tym_data.standard_name = "unweighted_fCO2_Tym"
      fCO2_Tym_data.long_name = "CO2 fugacity using OC-FLUX methodology (unweighted)"

      pCO2_SST_data = ncfile.createVariable('unweighted_pCO2_SST','f4',('time','latitude','longitude'),
                                            fill_value=netcdf_helper.MISSINGDATAVALUE,zlib=True)
      pCO2_SST_data[:] = newvars['pCO2_SST'][:]
      pCO2_SST_data.units = 'uatm'
      pCO2_SST_data.missing_value = netcdf_helper.MISSINGDATAVALUE
      pCO2_SST_data.valid_min = 0.
      pCO2_SST_data.valid_max = 1e6
      pCO2_SST_data.scale_factor = 1.
      pCO2_SST_data.add_offset = 0.
      pCO2_SST_data.standard_name = "unweighted_pCO2_SST"
      pCO2_SST_data.long_name = "CO2 partial pressure using SOCAT methodology (unweighted)"

      pCO2_Tym_data = ncfile.createVariable('unweighted_pCO2_Tym','f4',('time','latitude','longitude'),
                                            fill_value=netcdf_helper.MISSINGDATAVALUE,zlib=True)
      pCO2_Tym_data[:] = newvars['pCO2_Tym'][:]
      pCO2_Tym_data.units = 'uatm'
      pCO2_Tym_data.missing_value = netcdf_helper.MISSINGDATAVALUE
      pCO2_Tym_data.valid_min = 0.
      pCO2_Tym_data.valid_max = 1e6
      pCO2_Tym_data.scale_factor = 1.
      pCO2_Tym_data.add_offset = 0.
      pCO2_Tym_data.standard_name = "unweighted_pCO2_Tym"
      pCO2_Tym_data.long_name = "CO2 partial pressure using OC-FLUX methodology (unweighted)"

      fCO2_SST_data = ncfile.createVariable('unweighted_std_fCO2_SST','f4',('time','latitude','longitude'),
                                            fill_value=netcdf_helper.MISSINGDATAVALUE,zlib=True)
      fCO2_SST_data[:] = newvars['stds']['fCO2_SST'][:]
      fCO2_SST_data.units = 'uatm'
      fCO2_SST_data.missing_value = netcdf_helper.MISSINGDATAVALUE
      fCO2_SST_data.valid_min = 0.
      fCO2_SST_data.valid_max = 1e6
      fCO2_SST_data.scale_factor = 1.
      fCO2_SST_data.add_offset = 0.
      fCO2_SST_data.standard_name = "unweighted_std_fCO2_SST"
      fCO2_SST_data.long_name = "Standard deviation of CO2 fugacity using SOCAT methodology (unweighted)"

      fCO2_Tym_data = ncfile.createVariable('unweighted_std_fCO2_Tym','f4',('time','latitude','longitude'),
                                            fill_value=netcdf_helper.MISSINGDATAVALUE,zlib=True)
      fCO2_Tym_data[:] = newvars['stds']['fCO2_Tym'][:]
      fCO2_Tym_data.units = 'uatm'
      fCO2_Tym_data.missing_value = netcdf_helper.MISSINGDATAVALUE
      fCO2_Tym_data.valid_min = 0.
      fCO2_Tym_data.valid_max = 1e6
      fCO2_Tym_data.scale_factor = 1.
      fCO2_Tym_data.add_offset = 0.
      fCO2_Tym_data.standard_name = "unweighted_std_fCO2_Tym"
      fCO2_Tym_data.long_name = "Standard deviation of CO2 fugacity using OC-FLUX methodology (unweighted)"

      pCO2_SST_data = ncfile.createVariable('unweighted_std_pCO2_SST','f4',('time','latitude','longitude'),
                                            fill_value=netcdf_helper.MISSINGDATAVALUE,zlib=True)
      pCO2_SST_data[:] = newvars['stds']['pCO2_SST'][:]
      pCO2_SST_data.units = 'uatm'
      pCO2_SST_data.missing_value = netcdf_helper.MISSINGDATAVALUE
      pCO2_SST_data.valid_min = 0.
      pCO2_SST_data.valid_max = 1e6
      pCO2_SST_data.scale_factor = 1.
      pCO2_SST_data.add_offset = 0.
      pCO2_SST_data.standard_name = "unweighted_std_pCO2_SST"
      pCO2_SST_data.long_name = "Standard deviation of CO2 partial pressure using SOCAT methodology (unweighted)"

      pCO2_Tym_data = ncfile.createVariable('unweighted_std_pCO2_Tym','f4',('time','latitude','longitude'),
                                            fill_value=netcdf_helper.MISSINGDATAVALUE,zlib=True)
      pCO2_Tym_data[:] = newvars['stds']['pCO2_Tym'][:]
      pCO2_Tym_data.units = 'uatm'
      pCO2_Tym_data.missing_value = netcdf_helper.MISSINGDATAVALUE
      pCO2_Tym_data.valid_min = 0.
      pCO2_Tym_data.valid_max = 1e6
      pCO2_Tym_data.scale_factor = 1.
      pCO2_Tym_data.add_offset = 0.
      pCO2_Tym_data.standard_name = "unweighted_std_pCO2_Tym"
      pCO2_Tym_data.long_name = "Standard deviation of CO2 partial pressure using OC-FLUX methodology (unweighted)"

#Function for combining the global year/month files into an average for the whole stack of dates
def Climatology(filelist,output,varswewant):

   dataitem_dims={}
   itemdata={}
   itemvar={}
   sumN=0

   if not isinstance(filelist,list):
      raise Exception("filelist should be a list: in FromFilelist.")

   #if filelist is empty then return as nothing to do
   if len(filelist)==0:
      return False

   #if output is not a string then raise exception
   if not isinstance(output,str):
      raise Exception("output should be a string filename to write output to: in FromFilelist.")

   for filename in filelist:
      #open the file
      fin=netCDF4.Dataset(filename)
      #We only want the fCO2 and pCO2 from ocflux
      filevars=list(fin.variables.keys())
      filedims=list(fin.dimensions.keys())
      for variable in filevars:
         if variable in filedims:
            #variable is a dimension - check if it is already in the dimension dictionary
            if variable not in list(dataitem_dims.keys()):
               dataitem_dims[variable]=copy.copy(fin.variables[variable][:])
            else:
               #test if the numpy arrays are the same (i.e. same dimensions)
               if not numpy.array_equal(dataitem_dims[variable],fin.variables[variable][:]) and variable!="time":
                  print("Dimensions are not the same for %s  what does this mean, they should be. "%variable)
                  print(dataitem_dims[variable],fin.variables[variable][:])
                  sys.exit(1)
               elif variable == "time": 
                  #time could be a different value so this is OK
                  dataitem_dims['time']=numpy.vstack([dataitem_dims['time'],copy.copy(fin.variables[variable][:])])
         if variable in varswewant:
            if variable in list(itemdata.keys()):
               itemdata[variable]=numpy.ma.vstack((itemdata[variable],fin.variables[variable][:]))
            else:
               itemdata[variable]=fin.variables[variable][:].astype(numpy.float64)
               itemvar[variable]=fin.variables[variable]

   for variable in varswewant:
      itemdata[variable]=numpy.ma.mean(itemdata[variable],axis=0)

   with netCDF4.Dataset(output, 'w', format = 'NETCDF4') as ncfile:
      netcdf_helper.standard_setup_SOCAT(ncfile,timedata=1e9,londata=dataitem_dims['longitude'],latdata=dataitem_dims['latitude'])
      for variable in varswewant:
         if variable not in filedims:
            tmpvar = ncfile.createVariable(variable,'f4',('time','latitude','longitude'),
                                           fill_value=netcdf_helper.MISSINGDATAVALUE,zlib=True)
            tmpvar[:] = itemdata[variable][:]
            tmpvar.units = itemvar[variable].units
            tmpvar.missing_value = itemvar[variable].missing_value
            tmpvar.valid_min = itemvar[variable].valid_min
            tmpvar.valid_max = itemvar[variable].valid_max
            tmpvar.scale_factor = itemvar[variable].scale_factor
            tmpvar.add_offset = itemvar[variable].add_offset
            tmpvar.standard_name = itemvar[variable].standard_name
            tmpvar.long_name = itemvar[variable].long_name
   return True

