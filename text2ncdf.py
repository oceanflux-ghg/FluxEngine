#! /usr/bin/env python
#Created by Peter Land, Plymouth Marine Laboratory.
#IGA branch to create files with different geographical grids (i.e. not necessarily regular)

'''Program to convert text files to netcdf files with the correct dimensions to be used as input for the Flux Engine software. 
Input inFiles should be csv, or tab delimited text file(s). Headers are interrogated for lat/long, fCO2, N2O, CH4, SST and dates.
IGA ??/??/2016 Updated variable names to include CH4 and N20
IGA 5/5/2016 Updated filenames to allow for Arctic (or other) grid
IGA 11/5/2016 Updated to include vCO2_air data as req'd in FE '''

#TO DO - Make generic for gas names so any gas can be used. 

import numpy as np, argparse, csv, sys
from netCDF4 import Dataset
from datetime import datetime, timedelta
from glob import iglob, glob

MISSING = -999.
defaultNcFile = "/Users/iga202/Ifremer/OCEANFLUX-SHARED/workspace/i.g.c.ashton_exeter.ac.uk/text2ncdfout.nc"
defaultRefNcFile = "/Users/iga202/Ifremer/OCEANFLUX-DATA/composites/air_pressure_at_sea_level/arctic/ecmwf/2008/200801_OCF-PRE-ARC-1M-25km-ECMWF.nc"#choose local file that has the correct dimensions for output.

delimiters = [',', '\t', 'None']
latNames = ["latitude", "lat (n)", "lat"]
lonNames = ["longitude", "lon (e)", "lon"]
fCO2Names = ["fco2_rec","fco2_recommended","pco2"]
pN2ONames = ["pn2o"]
pCH4Names = ["pch4"]
SSTNames = ["temp","sst"]
vco2Names = ['vco2', 'vco2_air', 'conc_air']
dateTimeNames = ["date/time", "date"]
yearNames = ["yr"]
monthNames = ["mon"]
dayNames = ["day"]
hourNames = ["hour"]
minuteNames = ["min"]
epoch = datetime(1970,1,1)
invalidlatlon_count = 0

def datetimeFormat(datestring, myFormat = None):
  '''try reading date/time string in various formats'''
  datestring = datestring.rstrip()
  formats = ["%d/%m/%Y", "%d/%m/%Y %H:%M", "%d/%m/%Y %H:%M:%S", "%d/%m/%Y%H:%M:%S", "%d/%m/%y", "%d/%m/%y %H:%M", "%d/%m/%y %H:%M:%S"]
  if myFormat is not None: # try current format first
    try:
      return datetime.strptime(datestring, myFormat), myFormat
    except:
      pass
  for myFormat in formats:
    try:
      return datetime.strptime(datestring, myFormat), myFormat
    except:
      pass
      print 'No format for date...',datestring
  return None

def stringIndex(inList, strings):
  '''return index of inList containing one of strings'''
  inList = [x.strip(' ') for x in inList]
  for string in strings:
    try:
      return inList.index(string)
    except:
      pass
  return None

parser = argparse.ArgumentParser()
parser.add_argument("inFiles", nargs="+", help="input csv file(s)")
parser.add_argument("-n", "--ncFile", help="output netCDF file")
parser.add_argument("-r", "--refNcFile", help="reference netCDF file")
parser.add_argument("-s", "--startTime",
  help="start date (and time) DD/MM/(YY)YY (HH:MM(:SS))")
parser.add_argument("-e", "--endTime",
  help="end date (and time) DD/MM/(YY)YY (HH:MM(:SS))")
parser.add_argument("-m", "--month", type = int,
  help="month to include (1-12)")
parser.add_argument("-l", "--limits", type = float, nargs = 4,
  help="[N latitude, S latitude, W longitude, E longitude]")
args = parser.parse_args()

if args.ncFile is None:
  ncFileName = defaultNcFile
else:
  ncFileName = args.ncFile
if args.refNcFile is None:
  refNcFileName = defaultRefNcFile
else:
  refNcFileName = args.refNcFile
justMonth = args.month is not None
if justMonth:
  month = args.month

#print 'ncFileName',ncFileName,'\nrefFileName',refNcFileName

with Dataset(ncFileName, 'w', format='NETCDF3_CLASSIC') as ncFile, \
    Dataset(refNcFileName) as refNcFile:
  
  # Create dimensions based on reference file
  dimlist = []#Will be list of dimensions
  for dim in refNcFile.dimensions.keys():
      dimlist.append(str(dim))
      ldim = refNcFile.dimensions[str(dim)].size
      ncFile.createDimension(str(dim), ldim)#IGA How to get the dimension

  xdim = [str(x) for x in refNcFile.dimensions.keys() if 'lon' in str(x).lower() or 'x' in str(x).lower()]
  ydim = [str(x) for x in refNcFile.dimensions.keys() if 'lat' in str(x).lower() or 'y' in str(x).lower()]
  xdim = xdim[0]
  ydim = ydim[0]
#Assuming time dimension is called time
  #print 'Dimensions read from input file -> ',dimlist

  #Check to see if RGB is dimension in reference file. If not, make it
  if 'RGB' not in refNcFile.dimensions.keys():
      ncFile.createDimension('RGB', 3)
      dimlist.append('RGB')

  #Find the long and lat variables
  latvarname = [str(x) for x in refNcFile.variables.keys() if 'lat' in str(x).lower()]
  lonvarname = [str(x) for x in refNcFile.variables.keys() if 'lon' in str(x).lower()]

  if not lonvarname or not latvarname:
    print "data variable 'latitude' missing from reference netcdf input"
    sys.exit(1)
  
  if len(lonvarname)>1 or len(latvarname)>1:
    print "Unable to determine Latitude and/or longitude variables from reference netcdf input (other variables with lat or lon in their names?"
    sys.exit(1)
  
  latvarname = latvarname[0]
  lonvarname = lonvarname[0]

  latitudeVar = refNcFile.variables[latvarname]
  latitudeData = latitudeVar[:]
  nlat = latitudeData.shape[0]
  
  longitudeVar = refNcFile.variables[lonvarname]
  longitudeData = longitudeVar[:]
  nlon = longitudeData.shape[0]

  if len(latitudeData.shape)<2:
     #Geographical data (lat-lon) are expressed as vectors
     lat0 = latitudeData[0]
     dlat = latitudeData[1] - lat0
     lon0 = longitudeData[0]
     dlon = longitudeData[1] - lon0
  else:
     nlon = longitudeData.shape[1]

  try:
    timeVar = refNcFile.variables['time']
  except:
    print "data variable 'time' missing from reference netcdf input"
    sys.exit(1)
  
  myFormat = "%d/%m/%Y %H:%M:%S"
  noTimes = True
  if args.startTime is None:
    startTime = datetime(1900,1,1)
  else:
    startTime, myFormat = datetimeFormat(args.startTime, myFormat)
    if startTime is None:
      print "Unable to parse start datetime ", args.startTime
      sys.exit(1)
  maxTime = startTime
  if args.endTime is None:
    endTime = datetime(2100,12,31)
  else:
    endTime, myFormat = datetimeFormat(args.endTime, myFormat)
    if len(myFormat) == 8: # d/m/y
      endTime += timedelta(1) # add a day, so if endTime=31/1, we include all times before 1/2 00:00
    if endTime is None:
      print "Unable to parse end datetime ", args.endTime
      sys.exit(1)
  minTime = endTime

#Create variables
  secs = ncFile.createVariable('time', 'd', ('time',))
  secs.units = 'seconds since 1970-01-01 00:00:00'
  secs.axis = "T"
  secs.long_name = "Time - seconds since 1970-01-01 00:00:00"
  secs.standard_name = "time"
  secs.valid_min = 0. 
  secs.valid_max = 1.79769313486232e+308

  print '\nDimensions in ncFile',ncFile.dimensions.keys()
  print '\nx and y dims',xdim,ydim
  if len(latitudeData.shape)<2:
  #     #Geographical data (lat-lon) are expressed as a vector not a grid
    lats = ncFile.createVariable('latitude','f',(ydim,))
    lons = ncFile.createVariable('longitude','f',(xdim,))
  else:
    lats = ncFile.createVariable('latitude','f',(ydim,xdim))
    lons = ncFile.createVariable('longitude','f',(ydim,xdim))
  lats.units = 'degrees_north'
  lats.axis = "Y"
  lats.long_name = "Latitude North"
  lats.standard_name = "latitude"
  lats[:] = latitudeData
  lats.valid_min = -90.
  lats.valid_max = 90.

  lons.units = 'degrees_east'
  lons.axis = "X"
  lons.long_name = "Longitude East"
  lons.standard_name = "longitude"
  lons[:] = longitudeData
  lons.valid_min = -180.
  lons.valid_max = 180.

  countVar = ncFile.createVariable("fCO2_count", 'i4', ('time',ydim,xdim),
    fill_value=0)
  countVar.units = "count"
  countVar.missing_value = 0
  countVar.valid_min = 0
  countVar.valid_max = 32767
  countVar.standard_name = "observation count"
  countVar.long_name = "Number of fCO2 observations per cell"
  count = np.zeros([nlat, nlon], dtype = 'int')
  print 'initial count size',count.shape

  seasonVar = ncFile.createVariable("seasonal", 'i1', (ydim,xdim,'RGB'))
  seasonVar.units = "seasonal"
  seasonVar.missing_value = 0
  seasonVar.valid_min = 0
  seasonVar.valid_max = 32767
  seasonVar.standard_name = "seasonal coverage"
  seasonVar.long_name = "cumulative seasonal coverage (R=Nov-Feb,G=Mar-Jun,B=Jul-Oct)"
  season = np.zeros([nlat, nlon, 3], dtype = 'byte')

  fCO2Used = False
  SSTUsed = False
  N2OUsed = False
  CH4Used = False
  vco2Used = False
  
  for inFiles in args.inFiles:
    for inFile in iglob(inFiles):
      fCO2min = float('inf')
      N2Omin = float('inf')
      CH4min = float('inf')
      vco2min = float('inf')
      with open(inFile, 'U') as inUnit:
        delimiterFound = False
        for delimiter in delimiters:
          if delimiter == 'None':
            print "Unrecognised file delimiter"
            exit(0)
          reader=csv.reader(inUnit, delimiter=delimiter)
          for header in reader: # keep going til we find a header line
            lowHeader = [item.lower() for item in header]
            latIndex = stringIndex(lowHeader, latNames)
            lonIndex = stringIndex(lowHeader, lonNames)
            if latIndex is None or lonIndex is None:
              print 'lat or long not recognised'
              continue # pre-header info
            delimiterFound = True
            fCO2Index = stringIndex(lowHeader, fCO2Names)
            if fCO2Index is not None and not fCO2Used:
              print 'fCO2 data read from ',fCO2Index#IGA
              fCO2Used = True
              fCO2s = ncFile.createVariable('fCO2_bin','f',('time',ydim,xdim))
              fCO2s.units = 'uatm'
              fCO2s.standard_name = "fCO2"
              fCO2s.long_name = "Carbon dioxide fugacity"
              fCO2s.valid_min = 0.
              fCO2s.valid_max = 1e6
              fCO2s.fill_value = -999.
              fCO2 = np.zeros([nlat, nlon]) # sum fCO2
              fCO2stds = ncFile.createVariable('fCO2_std','f',('time',ydim,xdim))
              fCO2stds.units = 'uatm'
              fCO2stds.standard_name = "fCO2"
              fCO2stds.long_name = "Carbon dioxide fugacity"
              fCO2stds.valid_min = 0.
              fCO2stds.valid_max = 1e6
              fCO2stds.fill_value = -999.
              fCO2std = np.zeros([nlat, nlon]) # sum fCO2^2
            SSTIndex = stringIndex(lowHeader, SSTNames)
            if SSTIndex is not None and not SSTUsed:
              print 'SST data read from ',SSTIndex#IGA
              SSTUsed = True
              SSTs = ncFile.createVariable('SST','f',('time',ydim,xdim))
              SSTs.units = 'Kelvins'
              SSTs.standard_name = "SST"
              SSTs.long_name = "Sea surface temperature"
              SSTs.valid_min = 273.15 - 1.8
              SSTs.valid_max = 273.15 + 30.5
              SSTs.fill_value = -999.
              SST = np.zeros([nlat, nlon]) # sum SST
            N2OIndex = stringIndex(lowHeader, pN2ONames)
            if N2OIndex is not None and not N2OUsed:
              print 'N2O data read from ',N2OIndex#IGA
              N2OUsed = True
              N2Os = ncFile.createVariable('N2O_bin','f',('time',ydim,xdim))
              N2Os.units = 'uatm'
              N2Os.standard_name = "N2O"
              N2Os.long_name = "Partial presure Nitrous Oxide"
              N2Os.valid_min = 0.
              N2Os.valid_max = 1e6
              N2Os.fill_value = -999.
              N2O = np.zeros([nlat, nlon]) # sum N2O
              N2Ostds = ncFile.createVariable('N2O_std','f',('time',ydim,xdim))
              N2Ostds.units = 'uatm'
              N2Ostds.standard_name = "N2O"
              N2Ostds.long_name = "Partial presure Nitrous Oxide"
              N2Ostds.valid_min = 0.
              N2Ostds.valid_max = 1e6
              N2Ostds.fill_value = -999.
              N2Ostd = np.zeros([nlat, nlon]) # sum N2O^2
            CH4Index = stringIndex(lowHeader, pCH4Names)
            if CH4Index is not None and not CH4Used:
              print 'CH4 data read from ',CH4Index#IGA
              CH4Used = True
              CH4s = ncFile.createVariable('CH4_bin','f',('time',ydim,xdim))
              CH4s.units = 'uatm'
              CH4s.standard_name = "CH4"
              CH4s.long_name = "Partial pressure CH4"
              CH4s.valid_min = 0.
              CH4s.valid_max = 1e6
              CH4s.fill_value = -999.
              CH4 = np.zeros([nlat, nlon]) # sum CH4
              CH4stds = ncFile.createVariable('CH4_std','f',('time',ydim,xdim))
              CH4stds.units = 'uatm'
              CH4stds.standard_name = "CH4"
              CH4stds.long_name = "Partial pressure CH4"
              CH4stds.valid_min = 0.
              CH4stds.valid_max = 1e6
              CH4stds.fill_value = -999.
              CH4std = np.zeros([nlat, nlon]) # sum CH4^2
            vco2Index = stringIndex(lowHeader, vco2Names)
            if vco2Index is not None and not vco2Used:
              print 'vCO2 data read from ',vco2Index#IGA
              vco2Used = True
              vco2s = ncFile.createVariable('vco2_bin','f',('time',ydim,xdim))
              vco2s.units = 'uatm'
              vco2s.standard_name = "vco2"
              vco2s.long_name = "Partial pressure vco2"
              vco2s.valid_min = 0.
              vco2s.valid_max = 1e6
              vco2s.fill_value = -999.
              vco2 = np.zeros([nlat, nlon]) # sum vco2
              vco2stds = ncFile.createVariable('vco2_std','f',('time',ydim,xdim))
              vco2stds.units = 'uatm'
              vco2stds.standard_name = "vco2"
              vco2stds.long_name = "Partial pressure vco2"
              vco2stds.valid_min = 0.
              vco2stds.valid_max = 1e6
              vco2stds.fill_value = -999.
              vco2std = np.zeros([nlat, nlon]) # sum vco2^2
            dateTimeIndex = stringIndex(lowHeader, dateTimeNames)
            if dateTimeIndex is not None:
              noTimes = False
            else:
              yearIndex = stringIndex(lowHeader, yearNames)
              if yearIndex is not None:
                noTimes = False
              monthIndex = stringIndex(lowHeader, monthNames)
              dayIndex = stringIndex(lowHeader, dayNames)
              hourIndex = stringIndex(lowHeader, hourNames)
              minuteIndex = stringIndex(lowHeader, minuteNames)
            break
          if delimiterFound:
            break
          else:
            inUnit.seek(0)
        for line in reader:
          try:
            lat, lon = float(line[latIndex]), float(line[lonIndex])
          except:
            invalidlatlon_count +=1
            continue#IGA changed from system exit to allow for missing or invalid rows
          if lat < -90. or lat > 90.:
            invalidlatlon_count +=1
            continue#IGA - changed from system exit to allow for missing or invalid rows
          while lon <= -180.:
            lon += 360.
          while lon > 180.:
            lon -= 360.
          try:
            if len(latitudeData.shape)>1:
                #Geographical data (lat-lon) are expressed as a grid
                deltapos = abs(lon-longitudeData)+abs(lat-latitudeData)#Identifying correct position in geographical grid
                indx,indy = np.where(deltapos==deltapos.min())
            else: #Geographical data (lat-lon) are expressed as vectors
               indx = int((lat - lat0)/dlat+.5)
               indy = int((lon - lon0)/dlon+.5)
               # indx = indx[0]
               # indy = indy[0]
          except:
              invalidlatlon_count +=1
              continue
          # if indx < 0 or indx >= nlat or indy < 0 or indy >= nlon:#IGA removed as not compatible with non-regular grids
          #   print 'indx,indy,nlat,nlon',indx,indy,nlat,nlon
          #   continue
          if not noTimes:
            if dateTimeIndex is not None:
              dateTime, myFormat = datetimeFormat(line[dateTimeIndex], myFormat)
            else:
              dateTime = datetime(int(line[yearIndex]),int(line[monthIndex]),
                int(line[dayIndex]),int(line[hourIndex]),int(line[minuteIndex]))
            if dateTime < startTime or dateTime >= endTime:
              continue
            if justMonth and dateTime.month != month:
              continue
            if dateTime < minTime:
              minTime = dateTime
            if dateTime > maxTime:
              maxTime = dateTime
            if dateTime.month < 3 or dateTime.month > 10:
              season[indx, indy, 0] = 255
            elif dateTime.month < 7:
              season[indx, indy, 1] = 255
            else:
              season[indx, indy, 2] = 255
          else:
            season[indx, indy, :] = 255
          if count[indx, indy] == 32767:
            continue
          if fCO2Index is not None:
            fCO2i = float(line[fCO2Index])
            if fCO2i == MISSING:
              continue
            if fCO2i < fCO2min:
              fCO2min = fCO2i
            fCO2[indx, indy] += fCO2i
            fCO2std[indx, indy] += fCO2i**2
          #N2O
          if N2OIndex is not None:
            N2Oi = float(line[N2OIndex])
            if N2Oi == MISSING:
              continue
            if N2Oi < N2Omin:
              N2Omin = N2Oi
            N2O[indx, indy] += N2Oi
            N2Ostd[indx, indy] += N2Oi**2
          #CH4
          if CH4Index is not None:
            try:
              CH4i = float(line[CH4Index])
            except ValueError:
              CH4i = float('NaN')
            if CH4i == MISSING:
              continue
            if CH4i < CH4min:
              CH4min = CH4i
            CH4[indx, indy] += CH4i
            CH4std[indx, indy] += CH4i**2
          #vco2
          if vco2Index is not None:
            try:
              vco2i = float(line[vco2Index])
            except ValueError:
              vco2i = float('NaN')
            if vco2i == MISSING:
              continue
            if vco2i < vco2min:
              vco2min = vco2i
            vco2[indx, indy] += vco2i
            vco2std[indx, indy] += vco2i**2

          if SSTIndex is not None:
            try:
                SST[indx, indy] += 273.15 + float(line[SSTIndex])
            except ValueError:
                SST[indx, indy] += float('NaN')
          count[indx, indy] += 1
          # if count[indx, indy] == 32767: # average the first 32767 values#IGA removed as not sure why arbitrary 32767 limit
          #   if fCO2Index is not None:
          #     fCO2std[indx, indy] = np.sqrt((fCO2std[indx, indy] - fCO2[indx, indy] / count[indx, indy]) / count[indx, indy])
          #     fCO2[indx, indy] /= count[indx, indy]
          #   #N2O
          #   if N2OIndex is not None:
          #     N2Ostd[indx, indy] = np.sqrt((N2Ostd[indx, indy] - N2O[indx, indy] / count[indx, indy]) / count[indx, indy])
          #     N2O[indx, indy] /= count[indx, indy]
          #   #CH4
          #   if CH4Index is not None:
          #     CH4std[indx, indy] = np.sqrt((CH4std[indx, indy] - CH4[indx, indy] / count[indx, indy]) / count[indx, indy])
          #     CH4[indx, indy] /= count[indx, indy]
          #   #vco2
          #   if vco2Index is not None:
          #     vco2std[indx, indy] = np.sqrt((vco2std[indx, indy] - vco2[indx, indy] / count[indx, indy]) / count[indx, indy])
          #     vco2[indx, indy] /= count[indx, indy]
          #   if SSTIndex is not None:
          #     SST[indx, indy] /= count[indx, indy]
      if fCO2min < 250.:
         print inFile, fCO2min

  filled = np.nonzero(count == 0)
  valid = np.nonzero(np.logical_and(count > 0,count < 9999999999))#IGA altered from < 32767
  if fCO2Used:
    fCO2[filled] = -999.
    fCO2std[filled] = -999.
    fCO2std[valid] = np.sqrt((fCO2std[valid] - fCO2[valid] / count[valid]) / count[valid])
    fCO2[valid] /= count[valid]
    fCO2s[:] = fCO2
    fCO2stds[:] = fCO2std
  if SSTUsed:
    SST[filled] = -999.
    SST[valid] /= count[valid]
    SSTs[:] = SST
  if CH4Used:
    CH4[filled] = -999.
    CH4std[filled] = -999.
    CH4std[valid] = np.sqrt((CH4std[valid] - CH4[valid] / count[valid]) / count[valid])
    CH4[valid] /= count[valid]
    CH4s[:] = CH4
    CH4stds[:] = CH4std
  if vco2Used:
    vco2[filled] = -999.
    vco2std[filled] = -999.
    vco2std[valid] = np.sqrt((vco2std[valid] - vco2[valid] / count[valid]) / count[valid])
    vco2[valid] /= count[valid]
    vco2s[:] = vco2
    vco2stds[:] = vco2std
  if N2OUsed:
    N2O[filled] = -999.
    N2Ostd[filled] = -999.
    N2Ostd[valid] = np.sqrt((N2Ostd[valid] - N2O[valid] / count[valid]) / count[valid])
    N2O[valid] /= count[valid]
    N2Os[:] = N2O
    N2Ostds[:] = N2Ostd

  countVar[:] = count
  seasonVar[:] = season

  if noTimes:
    if args.startTime is None or args.endTime is None:
      secs[:] = timeVar[:]
    else: # inner and outer brackets are timedeltas
      secs[:] = (startTime + (endTime - startTime) / 2 -epoch).total_seconds()
  else: # inner and outer brackets are timedeltas
    secs[:] = (minTime + (maxTime - minTime) / 2 - epoch).total_seconds()

   # set some global attributes
  setattr(ncFile, 'Conventions', 'CF-1.0') 
  setattr(ncFile, 'institutions', 'Plymouth Marine Laboratory Remote Sensing Group/University of Exeter Centre for Geog, Environment and Society') 
  setattr(ncFile, 'contact1', 'email: rsghelp@pml.ac.uk') 
  setattr(ncFile, 'contact2', 'email: iga202@ex.ac.uk') 
  setattr(ncFile, 'RSG_areacode', 'a7')
  setattr(ncFile, 'RSG_sensor', "NETCDF") 

print invalidlatlon_count," data had invalid lat/lon co-ordinates"
print "(text2ncdf, main) SUCCESS writing file %s" % ncFileName
