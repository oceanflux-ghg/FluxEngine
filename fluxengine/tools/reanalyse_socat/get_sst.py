#!/home/cercache/tools/environments/run_python_scibox_precise.sh
#-*- coding:utf-8 -*-

import netCDF4
import numpy
import os
from scipy.ndimage import map_coordinates

def GetReynoldsSST(years, months, lons, lats, SSTdir, SSTtail):
   return GetSST(years, months, lons, lats, SSTdir, SSTtail,dataname='sst_mean',lonname='lon',latname='lat')

def GetAATSRSST(years, months, lons, lats, SSTdir, SSTtail):
   return GetSST(years, months, lons, lats, SSTdir, SSTtail,dataname='sst_skin_mean',lonname='lon',latname='lat')


def ReadSSTFile(filename,dataname='sst_skin_mean',lonname='lon',latname='lat'):
   #Data is 362 to deal with interpolation across dateline
   #As well as lon - to extend by 1 either side
   data=numpy.zeros([180,362])
   lons=numpy.zeros([362])
   with netCDF4.Dataset(filename,'r') as SST_file:
      data[:,1:361] = SST_file.variables[dataname][:]
      #The last line is a copy of the first (dateline)
      data[:,361] = data[:,1]
      data[:,0] = data[:,360]
      lons[1:361]=SST_file.variables[lonname][:]
      lons[0] = lons[360]
      lons[361]=lons[1]
      lats=SST_file.variables[latname][:]
      if SST_file.variables[dataname].units in ['degrees C']:
         #we want them in Kelvin (to be consistent with how the scripts were originally written)
         data=data+273.15
      else:
         print("SST units are assummed in Kelvin. If this is incorrect then convert the data to K in getsst.py (lines 30-34). ")
   return data,lons,lats

def GetSST(years, months, lons, lats, SSTdir, SSTtail,dataname,lonname,latname):
   """reads AATSR monthly climatology files and extracts data closest to the
      ship position:
      Arguments (all          #this is not good - suspect we are interpolating in an impossible area1D numpy arrays):
        years (INT)
        months (INT)
        lons - longitude in range (-180->180)
        lats - latitude in range (-90 -> 90)
        SSTdir - directory containing year subdirectories, each containing SST netcdfs
        SSTtail - remainder of the SST netcdf name after the year and month
      Return value: list of numpy arrays (equivalent inputs are overwritten):
        Tcl - climatological SST for each month and position
   """

   Tcl=numpy.ma.array([-999.] * years.size, fill_value = -999.)

   #Get a list of all year and month combinations from the data
   yrmon=list(zip(years,months))
   yrmon=sorted(list(set(yrmon)))

   #loop through each date set in turn
   for thisdate in yrmon:
      #Get the indices of points from this year and month
      indices=numpy.where((years==thisdate[0])&(months==thisdate[1]))[0]
      #Read in the data from the SST data file
      sstfilename=os.path.join(SSTdir,"%d"%thisdate[0],"%d%02d"%(thisdate[0],thisdate[1])+SSTtail)
      if os.path.isfile(sstfilename) == False: 
         print('%s: no SST file'%sstfilename)
         Tcl[indices] = -999
         continue
      else:
         sstdata,sstlons,sstlats=ReadSSTFile(sstfilename,dataname=dataname,lonname=lonname,latname=latname)
      #Now get the grid cell positions that relate to these
      #Start by rounding to the integer+0.5 
      X=numpy.floor(lons[indices])+0.5
      Y=numpy.floor(lats[indices])+0.5
      XY=list(zip(X,Y))
      lon_i=[]
      lat_i=[]
      #For each of the lon,lat pairs get the cell from the SST grid
      for i,xy in enumerate(XY):
         lonmatches=numpy.where(sstlons==xy[0])[0]
         if lonmatches.size==1:
            lon_i.append(lonmatches[0])
         else:
            #more than 1 match due to date line - where is the point we want to process
            #this decides which one to use
            if xy[0]-lons[indices[i]] >0:
               lon_i.append(lonmatches[1])
            else:
               lon_i.append(lonmatches[0])

         lat_i.append(numpy.where(sstlats==xy[1])[0][0])
      #Note we are zipping in order lat,lon
      grid_indices=numpy.array(list(zip(lat_i,lon_i)))
      #Now we need to get the fractional part of the grid index so that we can
      #interpolate to the position of the observation. Need to reshape due to 1 dim arrays.
      lonoffset=lons[indices]-sstlons[grid_indices[:,1]].reshape([-1])
      latoffset=lats[indices]-sstlats[grid_indices[:,0]].reshape([-1])
      #Add onto the grid cell indices
      grid_indices_offset=grid_indices+numpy.array(list(zip(latoffset,lonoffset)))
      #Now we can get the interpolated SST data for these points
      Tcl[indices]=map_coordinates(sstdata, grid_indices_offset.transpose(), order = 1)
      #We need to test the integrity of the interpolated data
      for index in range(indices.size):
         #Get the points used in the interpolation into a 1d array called window
         indexlat = int(numpy.floor(grid_indices_offset[index][0])); #TMH: converted to int
         indexlon = int(numpy.floor(grid_indices_offset[index][1])); #TMH: converted to int
         window=sstdata[indexlat:indexlat+2,indexlon:indexlon+2].reshape([-1])
         gooddata=numpy.where(window<9e9)[0]
         if gooddata.size == 4 or gooddata.size ==0:
            #all data are good so interpolation should be valid
            #or all data are bad and no interpolation can be done
            continue
         #Not all data were good - can we trust interpolation - need to test further
         window=window[gooddata]
         #Get all the other points that used exactly these data
         points=numpy.where((numpy.floor(grid_indices_offset[:,0])==indexlat)&
                            (numpy.floor(grid_indices_offset[:,1])==indexlon))
         #Get the weights as if doing a bi-linear interpolation
         a=grid_indices_offset[points,0][0]-indexlat
         b=grid_indices_offset[points,1][0]-indexlon
         a2=1-a
         b2=1-b
         weights=numpy.array([a2*b2,a2*b,a*b2,a*b])[gooddata].transpose()
         #Calculate the weighted value based only on the weights we DO have
         Tcl[indices[points]]=(window * weights).sum(axis=1)/weights.sum(axis=1)
         #Some points may be bad if there are no surrounding data
         bad=numpy.where(weights.max(axis=1) <= 0)
         Tcl[indices[bad]]=-999

   return Tcl

