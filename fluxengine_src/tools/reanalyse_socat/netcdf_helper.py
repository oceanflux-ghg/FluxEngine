#!/home/cercache/tools/environments/run_python_scibox_precise.sh
#-*- coding:utf-8 -*-

import netCDF4
import numpy
import datetime

#Missing data value to be entered into all netcdf output
MISSINGDATAVALUE=-999.0

#As the code was repeated in places it has been pulled into this separate function.
#This function sets up the netCDF file 'ncfile' with the 3 standard dimensions
#of time,latitude,longitude that are used in the SOCAT conversion scripts
#Timedata should be either a datetime object or a float-able object
def standard_setup_SOCAT(ncfile,timedata,londata,latdata):
   #Create and add the time dimension
   if isinstance(timedata,datetime.datetime):
      #need to +1 as, apparently, 1st Jan is 1 day since 1st Jan ...
      #or to put it another way the, 1st Jan 1970 should have value = 1
      timeobj=timedata-datetime.datetime(1970,0o1,0o1)
      timeobj=timeobj.days+(timedata.hour/24.0) +1
      ntimes=1
   elif isinstance(timedata,numpy.ndarray):
      timeobj=timedata.astype(float)
      ntimes=timeobj.size
   elif isinstance(timedata,list):
      timeobj=[float(t) for t in timedata]
      ntimes=len(timedata)
   else:
      timeobj=float(timedata)
      ntimes=1
   ncfile.createDimension('time', ntimes)
   days = ncfile.createVariable('time','f8',('time',),zlib=True)
   days.units = 'Days since 1970-01-01 00:00:00'
   days.axis = "T"
   days.long_name = "Time - days since 1970-01-01 00:00:00"
   days.standard_name = "time"
   days[:] = timeobj
   days.valid_min = 0.0 
   days.valid_max = 1.79769313486232e+308
   #Create and add the latitude dimension
   ncfile.createDimension('latitude', 180)
   lats = ncfile.createVariable('latitude','f4',('latitude',),zlib=True)
   lats.units = 'degrees_north'
   lats.axis = "Y"
   lats.long_name = "latitude"
   lats.standard_name = "latitude"
   lats[:] = latdata
   lats.valid_min = -90.
   lats.valid_max = 90.
   #Create and add the longitude dimension
   ncfile.createDimension('longitude', 360)
   lons = ncfile.createVariable('longitude','f4',('longitude',),zlib=True)
   lons.units = 'degrees_east'
   lons.axis = "X"
   lons.long_name = "longitude"
   lons.standard_name = "longitude"
   lons[:] = londata
   lons.valid_min = -180.
   lons.valid_max = 180.


def write_netcdf_vars_to_ascii(ncfilein,asciifile,varlist,delimiter="\t"):
   """
   Function that will write out gridded netCDF variables to an ascii (tab) separated list.
   Longitude and latitude variable names are required to exist in the netCDF.
   Missing data are not output. Each variable will be a separate column, e.g.:
   #longitude  latitude    fCO2  pCO2
   -9.5     81.5    336.282         337.768
   -6.5     81.5    336.318         337.805
   -12.5    80.5    274.271         275.469
   ........................................
   """
   #Check varlist is a list of strings
   if not isinstance(varlist,list):
      raise Exception("varlist should be a list of the string variable names in: write_netcdf_vars_to_ascii")
   else:
      for item in varlist:
         if not isinstance(item,str):
            raise Exception("varlist should be a list of the string variable names in: write_netcdf_vars_to_ascii")

   #Open the netcdf file
   try:
      ncfile=netCDF4.Dataset(ncfilein)
      longitude=ncfile.variables['longitude'][:]
      latitude=ncfile.variables['latitude'][:]
      initialvar=ncfile.variables[varlist[0]][:][0].filled()
   except Exception as e:
      print("Problem getting latitude and longitude from the netCDF file.")
      raise
   #Get the inidices where there are data
   goodindices=numpy.where(initialvar!=MISSINGDATAVALUE)
   #Create an output array of correct size
   data=numpy.zeros((len(varlist)+2,len(goodindices[0])))
   #Fill the longitudes
   data[0,:]=longitude[goodindices[1]]
   #Fill the latitudes
   data[1,:]=latitude[goodindices[0]]
   #Create a header for the output file
   header="longitude%slatitude"%(delimiter)
   for i,var in enumerate(varlist):
      data[i+2,:]=ncfile.variables[var][:][0][goodindices]
      header=header+"%s%s"%(delimiter,var)
   #Transpose for output
   data=numpy.transpose(data)
   #Write out
   numpy.savetxt(asciifile,data,fmt='%.6g',delimiter=delimiter,newline='\n',header=header)
