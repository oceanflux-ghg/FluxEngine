#IGA Jun 2016
#Created by Ian Ashton, University of Exeter
#iga202@ex.ac.uk

#Script to combine total alkalinity and dissolved inorganic carbon data into one netCDF file in preparation for the SeaCarb code
import numpy as np,csv
from netCDF4 import Dataset
from math import radians, cos, sin, asin, sqrt
import sys
sys.path.insert(0, '/Users/iga202/')

def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a)) 
    r = 6371 # Radius of earth in kilometers. Use 3956 for miles
    return c * r
      #IGA Jun 2016
      #Created by Ian Ashton, University of Exeter
      #iga202@ex.ac.uk


infilename = 'Ifremer/OCEANFLUX-SHARED/workspace/i.g.c.ashton_exeter.ac.uk/arctic_grid/arctic_mask.nc'
#outfilename = 'Ifremer/OCEANFLUX-SHARED/workspace/i.g.c.ashton_exeter.ac.uk/arctic_grid/arctic_grid_area_calc.nc'
outfilename = 'Documents/1RF/Arctic/arctic_grid_area_calc.nc'

with Dataset(outfilename,'w') as ncout:
      with Dataset(infilename,'r') as ncin:
            for dim in ncin.dimensions:
                ncout.createDimension(dim,ncin.dimensions[dim].size)#Assuming ncin and ncin have the same dimensions

            # #----------------------test start-------------------------------#
            # for key in ncin.variables:
            #       print key
            # #----------------------test end---------------------------------#

            #Create latitude and longitude variables
            latname = [nm for nm in ncin.variables.keys() if 'lat' in nm]
            latname = latname[0]
            lonname = [nm for nm in ncin.variables.keys() if 'lon' in nm]
            lonname = lonname[0]
            latvar = ncout.createVariable('latitude',ncin.variables[latname].datatype,(ncin.variables[latname].dimensions))
            latvar.units = ncin.variables[latname].units
            lat = ncin.variables[latname][:]
            latvar[:] = lat
            lonvar = ncout.createVariable('longitude',ncin.variables[lonname].datatype,(ncin.variables[lonname].dimensions))
            lonvar.units = ncin.variables[lonname].units
            lon = ncin.variables[lonname][:]
            lonvar[:] = lon

            #Calculate area
            gridArea = np.empty((448,304))

            for r in range(447):
                  for c in range(303):
                        horSide = haversine(lon[r,c], lat[r,c], lon[r,c+1], lat[r,c+1])
                        vertSide = haversine(lon[r,c], lat[r,c], lon[r+1,c], lat[r+1,c])
                        gridArea[r,c] = vertSide*horSide
            r = 447
            for c in range(303):
                  horSide = haversine(lon[r,c], lat[r,c], lon[r,c+1], lat[r,c+1])
                  vertSide = haversine(lon[r,c], lat[r,c], lon[r-1,c], lat[r-1,c])
                  gridArea[r,c] = vertSide*horSide
            c = 303
            for r in range(447):
                  horSide = haversine(lon[r,c], lat[r,c], lon[r,c-1], lat[r,c-1])
                  vertSide = haversine(lon[r,c], lat[r,c], lon[r+1,c], lat[r+1,c])
                  gridArea[r,c] = vertSide*horSide
            r = 447
            c = 303
            horSide = haversine(lon[r,c], lat[r,c], lon[r,c-1], lat[r,c-1])
            vertSide = haversine(lon[r,c], lat[r,c], lon[r-1,c], lat[r-1,c])
            gridArea[r,c] = vertSide*horSide

            areaVar = ncout.createVariable('area',ncin.variables['land_proportion'].datatype,(ncin.variables['land_proportion'].dimensions))
            areaVar.units = 'km2'
            areaVar[:] = gridArea

