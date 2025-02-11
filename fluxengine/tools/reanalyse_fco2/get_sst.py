#!/home/cercache/tools/environments/run_python_scibox_precise.sh
#-*- coding:utf-8 -*-

import netCDF4
import numpy
import os
from scipy.ndimage import map_coordinates

#DJF 06/02/2025 No longer required.
# def GetESACCISST(years, months, lons, lats, SSTdir, SSTtail,useaatsr=False,usereynolds=False,useESACCI=False,days = 0,daily = False,bias = 0.05 ):
#    #return GetSST(years, months, lons, lats, SSTdir, SSTtail,dataname='analysed_sst',lonname='longitude',latname='latitude',useaatsr=useaatsr,usereynolds=usereynolds,useESACCI=useESACCI,days=days,daily=daily,bias = bias )
#    return GetSST(years, months, lons, lats, SSTdir, SSTtail,dataname='analysed_sst',lonname='lon',latname='lat',useaatsr=useaatsr,usereynolds=usereynolds,useESACCI=useESACCI,days=days,daily=daily,bias = bias )
#
# def GetReynoldsSST(years, months, lons, lats, SSTdir, SSTtail,useaatsr=False,usereynolds=False,useESACCI=False,days = 0,daily = False):
#    return GetSST(years, months, lons, lats, SSTdir, SSTtail,dataname='sst_mean',lonname='lon',latname='lat',useaatsr=useaatsr,usereynolds=usereynolds,useESACCI=useESACCI,days=days,daily=daily,bias = bias )
#
# def GetAATSRSST(years, months, lons, lats, SSTdir, SSTtail,useaatsr=False,usereynolds=False,useESACCI=False,days = 0,daily = False):
#    return GetSST(years, months, lons, lats, SSTdir, SSTtail,dataname='sst_skin_mean',lonname='lon',latname='lat',useaatsr=useaatsr,usereynolds=usereynolds,useESACCI=useESACCI,days=days,daily=daily,bias = bias )


def ReadSSTFile(filename,dataname='sst_skin_mean',lonname='lon',latname='lat'):
   """
   Function to read the SST data from the netCDF file. We want it in (lat,lon) format.
   We assume the data is a global grid (-180 to 180), so we check this as well.

   DJF edited 21/01/2025 to remove the hardcoded 1 degree resolution, and added checks for dimensions,
   and lons being -180 to 180.
   """
   #DJF 21/01/2025: Old code commented out.
   #Data is 362 to deal with interpolation across dateline
   #As well as lon - to extend by 1 either side
   #data=numpy.zeros([180,362])
   #lons=numpy.zeros([362])
   with netCDF4.Dataset(filename,'r') as SST_file:
      lons_t = SST_file.variables[lonname][:] # Load longitude variable
      res = numpy.abs(lons_t[0]-lons_t[1])
      lons = numpy.zeros([len(lons_t)+2]) # Now we know the length we can setup a new longitude array

      lons[1:-1] = lons_t # Put the data in the middle
      # if daily: # For the daily we use a slightly different approach that allows us to use any resolution daily data
      lons[0] = lons_t[0]-res
      lons[-1] = lons_t[-1]+res
      # else:# So if its monthly data we use this approach.
      #     lons[0] = lons_t[-1] # Pad the dateline to deal with dateline interpolation
      #     lons[-1] = lons_t[0]
      del lons_t
      lats = SST_file.variables[latname][:]
      data = numpy.zeros([len(lats),len(lons)]) # Added 2 to the lons so we can do the date line addition as previously.
      sst_data_t = numpy.squeeze(SST_file.variables[dataname][:]) # Added a squeeze here as some daily diles have a time dimension of 1. This removes it, and does nothing if no 1 length dimensions exist.

      if sst_data_t.shape[0] != len(lats): # If the first dimension is not latitude the data must be (lon, lat) so we transpose.
          sst_data_t = sst_data_t.transpose()
      data[:,1:-1] = sst_data_t
      data[:,0] = sst_data_t[:,-1]
      data[:,-1] = sst_data_t[:,0]
      del sst_data_t

      #DJF 21/01/2025: Old code commented out
      #data[:,1:361] = SST_file.variables[dataname][:]
      #The last line is a copy of the first (dateline)
      # data[:,361] = data[:,1]
      # data[:,0] = data[:,360]
      # lons[1:361]=
      # lons[0] = lons[360]
      # lons[361]=lons[1]
      # lats=SST_file.variables[latname][:]
      if SST_file.variables[dataname].units in ['degrees C']:
         #we want them in Kelvin (to be consistent with how the scripts were originally written)
         data=data+273.15
      else:
         print("SST units are assummed in Kelvin. If this is incorrect then convert the data to K in getsst.py (lines 30-34). ")
   return data,lons,lats

def GetSST(years, months, lons, lats, SSTdir, SSTtail,dataname,lonname,latname,days=0,sst_bias=0,unc_extract = False,uncname=''): # DJF 06/02/2025: Removed the usaESACCI etc variables as no longer needed
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
   Tcl_unc=numpy.ma.array([-999.] * years.size, fill_value = -999.)

   # if daily:
   # Here we setup if we are going to run it in daily mode
   yrmon = list(zip(years,months,days)) # Getting  a list of all year, month, day combinations
   yrmon = sorted(list(set(yrmon)))
   print(yrmon)
   for thisdate in yrmon:
       print(thisdate)
       indices = numpy.where((years==thisdate[0]) & (months==thisdate[1]) & (days == thisdate[2]))[0]


       yearstring = thisdate[0]
       monthstring = thisdate[1]
       daystring = thisdate[2]
       # DJF 30/01/2025: Updated the sst file name generation to add flexibilty in the file names. Here we can now specify the filenaming with SST_Tail and allows
       # any combination of file paths to be used. For example, files can be in a '%Y/%m/%Y%m%d.nc' format for daily files or could be '%Y/%Y%m.nc' for monthly files.
       sstfilename=os.path.join(SSTdir,SSTtail.replace('%Y',str(yearstring)).replace('%m',format(monthstring, "02d")).replace('%d',format(daystring, "02d")))

       if os.path.isfile(sstfilename) == False:
           print('%s: no SST file'%sstfilename)
           Tcl[indices] = -999
           continue
       else:
           print('%s: File found, start loading...'%sstfilename)
           sstdata,sstlons,sstlats=ReadSSTFile(sstfilename,dataname=dataname,lonname=lonname,latname=latname)
           sstdata = sstdata+sst_bias # This allows a global bias to be applied to the data (i.e remove a cool bias in the data.)

           #DJF 09/02/2025: Adding ability to add unc information to the ASCII files.
           if unc_extract:
               print('Loading SST uncertainty data...')
               sstuncdata,sstlons,sstlats=ReadSSTFile(sstfilename,dataname=uncname,lonname=lonname,latname=latname)
           print('SST bias of '+ str(sst_bias) + ' applied')

       sst_data_res = numpy.abs(sstlons[0]-sstlons[1]) # Find the resolution of the sst data, we assume its the same in the latitude and longitudes

       # Now we work out the grid cells that each SOCAT observation corresponds to on a any dimension grid.
       X = lons[indices]
       Y = lats[indices]
       XY = list(zip(X,Y))
       lon_i = []
       lat_i = []
       #For each of the lon,lat pairs get the cell from the SST grid
       for i,xy in enumerate(XY):
           # Need to find the index in the longitude grid that is to the west of the point (to match the monthly code below...)
           t_lon = xy[0] - sstlons # Find the difference between all the sst_lons and the longitude of the SOCAT data

           lonmatches = numpy.where((t_lon >= 0) & (t_lon < sst_data_res))[0] # We find where this difference is greater or equal to 0 but up to a maximum of the resolution of the sst data

           if lonmatches.size == 1:
               lon_i.append(lonmatches[0])
           else:
               print('Multiple indexed values... suggests we fall exactly on the SST data longitude grid')
               print(t_lon)
               print(t_lon[lonmatches])
               print(sstlons[lonmatches])
               print(lonmatches)
               print(xy[0])
               tt = numpy.where(numpy.abs(t_lon) == numpy.min(numpy.abs(t_lon)))[0]
               print('TT = ' + str(tt))
               lon_i.append(tt[0])
               # raise Exception("Stopping for checks")

           t_lat = xy[1] - sstlats
           latmatches = numpy.where((t_lat >= 0) & (t_lat < sst_data_res))[0]
           if latmatches.size == 1:
               lat_i.append(latmatches[0])
           else:
               print('Multiple indexed values... suggests we fall exactly on the SST data latitude grid')
               print(t_lat)
               print(t_lat[latmatches])
               print(sstlats[latmatches])
               print(latmatches)
               print(xy[1])
               tt = numpy.where(numpy.abs(t_lat) == numpy.min(numpy.abs(t_lat)))[0]
               # print('TT = ' + str(tt))
               lat_i.append(tt[0])

       #Note we are zipping in order lat,lon
       grid_indices=numpy.array(list(zip(lat_i,lon_i)))
       #Now we need to get the fractional part of the grid index so that we can
       #interpolate to the position of the observation. Need to reshape due to 1 dim arrays.
       lonoffset=((lons[indices]-sstlons[grid_indices[:,1]])/sst_data_res).reshape([-1])
       latoffset=((lats[indices]-sstlats[grid_indices[:,0]])/sst_data_res).reshape([-1])
       #Add onto the grid cell indices
       grid_indices_offset=grid_indices+numpy.array(list(zip(latoffset,lonoffset)))

       Tcl[indices]=map_coordinates(sstdata, grid_indices_offset.transpose(), order = 1)
       if unc_extract:
           Tcl_unc[indices]=map_coordinates(sstuncdata, grid_indices_offset.transpose(), order = 1)
       #We need to test the integrity of the interpolated data
       for index in range(indices.size):
          #Get the points used in the interpolation into a 1d array called window
          indexlat = grid_indices[index,0]#int(numpy.floor(grid_indices_offset[index][0])); #TMH: converted to int
          indexlon = grid_indices[index,1]#int(numpy.floor(grid_indices_offset[index][1])); #TMH: converted to int
          window=sstdata[indexlat:indexlat+2,indexlon:indexlon+2].reshape([-1])
          if unc_extract:
              window_unc = sstuncdata[indexlat:indexlat+2,indexlon:indexlon+2].reshape([-1])
          gooddata=numpy.where((window<9e9) & (window > 0) & (numpy.isnan(window) == 0))[0] #DJF: Added second condition where the fill value is less than 0.
          if gooddata.size == 4 or gooddata.size ==0:
             #all data are good so interpolation should be valid
             #or all data are bad and no interpolation can be done
             continue
          #Not all data were good - can we trust interpolation - need to test further
          window=window[gooddata]
          if unc_extract:
              window_unc=window_unc[gooddata]
          # #Get all the other points that used exactly these data
          # points=numpy.where((numpy.floor(grid_indices_offset[:,0])==indexlat)&
          #                    (numpy.floor(grid_indices_offset[:,1])==indexlon))
          #Get the weights as if doing a bi-linear interpolation
          a=grid_indices_offset[index,0]-indexlat
          b=grid_indices_offset[index,1]-indexlon
          a2=1-a
          b2=1-b
          weights=numpy.array([a2*b2,a2*b,a*b2,a*b])[gooddata].transpose()
          #Calculate the weighted value based only on the weights we DO have
          Tcl[indices[index]]=(window * weights).sum()/weights.sum()
          if unc_extract:
              Tcl_unc[indices[index]] = (window_unc * weights).sum()/weights.sum()
          #Some points may be bad if there are no surrounding data
          if weights.max() <= 0:
            Tcl[indices[index]]=-999
            if unc_extract:
                Tcl_unc[indices[index]] = -999

   # else:
   #     #Get a list of all year and month combinations from the data
   #     yrmon=list(zip(years,months))
   #     yrmon=sorted(list(set(yrmon)))
   #
   #     #loop through each date set in turn
   #     for thisdate in yrmon:
   #        #Get the indices of points from this year and month
   #        indices=numpy.where((years==thisdate[0])&(months==thisdate[1]))[0]
   #        #Read in the data from the SST data file
   #        if useESACCI:
   #            yearstring=thisdate[0]
   #            monthstring=thisdate[1]
   #            sstfilename=os.path.join(SSTdir,"{0}/{0}{1}{2}".format(yearstring,format(monthstring, "02d"),SSTtail))
   #        else:
   #            sstfilename=os.path.join(SSTdir,"%d"%thisdate[0],"%d%02d"%(thisdate[0],thisdate[1])+SSTtail)
   #
   #        if os.path.isfile(sstfilename) == False:
   #           print('%s: no SST file'%sstfilename)
   #           Tcl[indices] = -999
   #           continue
   #        else:
   #           sstdata,sstlons,sstlats=ReadSSTFile(sstfilename,dataname=dataname,lonname=lonname,latname=latname)
   #        #Now get the grid cell positions that relate to these
   #        #Start by rounding to the integer+0.5
   #        X=numpy.floor(lons[indices])+0.5
   #        Y=numpy.floor(lats[indices])+0.5
   #        XY=list(zip(X,Y))
   #        lon_i=[]
   #        lat_i=[]
   #        #For each of the lon,lat pairs get the cell from the SST grid
   #        for i,xy in enumerate(XY):
   #           lonmatches=numpy.where(sstlons==xy[0])[0]
   #           if lonmatches.size==1:
   #              lon_i.append(lonmatches[0])
   #           else:
   #              #more than 1 match due to date line - where is the point we want to process
   #              #this decides which one to use
   #              if xy[0]-lons[indices[i]] >0:
   #                 lon_i.append(lonmatches[1])
   #              else:
   #                 lon_i.append(lonmatches[0])
   #
   #           lat_i.append(numpy.where(sstlats==xy[1])[0][0])
   #        #Note we are zipping in order lat,lon
   #        grid_indices=numpy.array(list(zip(lat_i,lon_i)))
   #        #Now we need to get the fractional part of the grid index so that we can
   #        #interpolate to the position of the observation. Need to reshape due to 1 dim arrays.
   #        lonoffset=lons[indices]-sstlons[grid_indices[:,1]].reshape([-1])
   #        latoffset=lats[indices]-sstlats[grid_indices[:,0]].reshape([-1])
   #        #Add onto the grid cell indices
   #        grid_indices_offset=grid_indices+numpy.array(list(zip(latoffset,lonoffset)))
   #        #Now we can get the interpolated SST data for these points
   #        Tcl[indices]=map_coordinates(sstdata, grid_indices_offset.transpose(), order = 1)
   #        #We need to test the integrity of the interpolated data
   #        for index in range(indices.size):
   #           #Get the points used in the interpolation into a 1d array called window
   #           indexlat = int(numpy.floor(grid_indices_offset[index][0])); #TMH: converted to int
   #           indexlon = int(numpy.floor(grid_indices_offset[index][1])); #TMH: converted to int
   #           window=sstdata[indexlat:indexlat+2,indexlon:indexlon+2].reshape([-1])
   #           gooddata=numpy.where((window<9e9) & (window > 0) & (numpy.isnan(window) == 0))[0] #DJF: Added second condition where the fill value is less than 0.
   #           if gooddata.size == 4 or gooddata.size ==0:
   #              #all data are good so interpolation should be valid
   #              #or all data are bad and no interpolation can be done
   #              continue
   #           #Not all data were good - can we trust interpolation - need to test further
   #           window=window[gooddata]
   #           #Get all the other points that used exactly these data
   #           points=numpy.where((numpy.floor(grid_indices_offset[:,0])==indexlat)&
   #                              (numpy.floor(grid_indices_offset[:,1])==indexlon))
   #           #Get the weights as if doing a bi-linear interpolation
   #           a=grid_indices_offset[points,0][0]-indexlat
   #           b=grid_indices_offset[points,1][0]-indexlon
   #           a2=1-a
   #           b2=1-b
   #           weights=numpy.array([a2*b2,a2*b,a*b2,a*b])[gooddata].transpose()
   #           #Calculate the weighted value based only on the weights we DO have
   #           Tcl[indices[points]]=(window * weights).sum(axis=1)/weights.sum(axis=1)
   #           #Some points may be bad if there are no surrounding data
   #           bad=numpy.where(weights.max(axis=1) <= 0)
   #           Tcl[indices[bad]]=-999

   return Tcl,Tcl_unc
