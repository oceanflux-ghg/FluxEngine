#! /usr/bin/env python
#Created by Peter Land, Plymouth Marine Laboratory.

'''Program to convert text files to netcdf files with the correct dimensions to be used as input for the Flux Engine software. 
Input inFiles should be csv, or tab delimited text file(s). Headers are interrogated for lat/long, fCO2, SST and dates.'''

import numpy as np, argparse, csv, sys
from netCDF4 import Dataset
from datetime import datetime, timedelta
from glob import iglob

MISSING = -999.
defaultNcFile = "default.nc"
defaultRefNcFile = "ref.nc"#choose local file that has the correct dimensions for output.
delimiters = [',', '\t', 'None']
latNames = ["latitude", "lat (n)", "lat"]
lonNames = ["longitude", "lon (e)", "lon"]
fCO2Names = ["fco2_rec","fco2_recommended"]
SSTNames = ["temp"]
dateTimeNames = ["date/time", "date"]
yearNames = ["yr"]
monthNames = ["mon"]
dayNames = ["day"]
hourNames = ["hour"]
minuteNames = ["min"]
epoch = datetime(1970,1,1)

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

with Dataset(ncFileName, 'w', format='NETCDF3_CLASSIC') as ncFile, \
    Dataset(refNcFileName) as refNcFile:

  try:
    latitudeVar = refNcFile.variables['latitude']
    latitudeData = latitudeVar[:]
    nlat = latitudeData.size
  except:
    print "data variable 'latitude' missing from reference netcdf input"
    sys.exit(1)
  lat0 = latitudeData[0]
  dlat = latitudeData[1] - lat0
  error = max(abs(latitudeData[1:] - latitudeData[:-1] - dlat))
  if error > 1e-5:
    print "data variable 'latitude' nonlinear", error
    sys.exit(1)
  try:
    longitudeVar = refNcFile.variables['longitude']
    longitudeData = longitudeVar[:]
    nlon = longitudeData.size
  except:
    print "data variable 'longitude' missing from reference netcdf input"
    sys.exit(1)
  lon0 = longitudeData[0]
  dlon = longitudeData[1] - lon0
  error = max(np.minimum(abs(longitudeData[1:] - longitudeData[:-1] - dlon),
    abs(longitudeData[1:] - longitudeData[:-1] + 360. - dlon)))
  if error > 1e-5:
    print "data variable 'longitude' nonlinear", error
    sys.exit(1)
  nlatlon = nlat * nlon
  print 'lat0,lon0,dlat,dlon ',lat0,lon0,dlat,dlon
  try:
    timeVar = refNcFile.variables['time']
  except:
    print "data variable 'time' missing from reference netcdf input"
    sys.exit(1)
  
  myFormat = None
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

  ncFile.createDimension('latitude', nlat)
  ncFile.createDimension('longitude', nlon)
  ncFile.createDimension('time', 1)
  ncFile.createDimension('RGB', 3)

  secs = ncFile.createVariable('time', 'd', ('time',))
  secs.units = 'seconds since 1970-01-01 00:00:00'
  secs.axis = "T"
  secs.long_name = "Time - seconds since 1970-01-01 00:00:00"
  secs.standard_name = "time"
  secs.valid_min = 0. 
  secs.valid_max = 1.79769313486232e+308

  lats = ncFile.createVariable('latitude','f',('latitude',))
  lons = ncFile.createVariable('longitude','f',('longitude',))
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

  countVar = ncFile.createVariable("fCO2_count", 'i4', ('time','latitude','longitude'),
    fill_value=0)
  countVar.units = "count"
  countVar.missing_value = 0
  countVar.valid_min = 0
  countVar.valid_max = 32767
  countVar.standard_name = "observation count"
  countVar.long_name = "Number of fCO2 observations per cell"
  count = np.zeros([nlat, nlon], dtype = 'int')

  seasonVar = ncFile.createVariable("seasonal", 'i1', ('latitude','longitude','RGB'))
  seasonVar.units = "seasonal"
  seasonVar.missing_value = 0
  seasonVar.valid_min = 0
  seasonVar.valid_max = 32767
  seasonVar.standard_name = "seasonal coverage"
  seasonVar.long_name = "cumulative seasonal coverage (R=Nov-Feb,G=Mar-Jun,B=Jul-Oct)"
  season = np.zeros([nlat, nlon, 3], dtype = 'byte')

  fCO2Used = False
  SSTUsed = False

  for inFiles in args.inFiles:
    for inFile in iglob(inFiles):
      fCO2min = float('inf')
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
              continue # pre-header info
            delimiterFound = True
            fCO2Index = stringIndex(lowHeader, fCO2Names)
            if fCO2Index is not None and not fCO2Used:
              fCO2Used = True
              fCO2s = ncFile.createVariable('fCO2_bin','f',('time','latitude','longitude'))
              fCO2s.units = 'uatm'
              fCO2s.standard_name = "fCO2"
              fCO2s.long_name = "Carbon dioxide fugacity"
              fCO2s.valid_min = 0.
              fCO2s.valid_max = 1e6
              fCO2s.fill_value = -999.
              fCO2 = np.zeros([nlat, nlon]) # sum fCO2
              fCO2stds = ncFile.createVariable('fCO2_std','f',('time','latitude','longitude'))
              fCO2stds.units = 'uatm'
              fCO2stds.standard_name = "fCO2"
              fCO2stds.long_name = "Carbon dioxide fugacity"
              fCO2stds.valid_min = 0.
              fCO2stds.valid_max = 1e6
              fCO2stds.fill_value = -999.
              fCO2std = np.zeros([nlat, nlon]) # sum fCO2^2
            SSTIndex = stringIndex(lowHeader, SSTNames)
            if SSTIndex is not None and not SSTUsed:
              SSTUsed = True
              SSTs = ncFile.createVariable('SST','f',('time','latitude','longitude'))
              SSTs.units = 'Kelvins'
              SSTs.standard_name = "SST"
              SSTs.long_name = "Sea surface temperature"
              SSTs.valid_min = 273.15 - 1.8
              SSTs.valid_max = 273.15 + 30.5
              SSTs.fill_value = -999.
              SST = np.zeros([nlat, nlon]) # sum SST
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
            sys.exit(1)
          if lat < -90. or lat > 90.:
            print "Invalid latitude ", line[latIndex]
            sys.exit(1)
          while lon <= -180.:
            lon += 360.
          while lon > 180.:
            lon -= 360.
          latPos = int((lat - lat0)/dlat+.5)
          lonPos = int((lon - lon0)/dlon+.5)
          if latPos < 0 or latPos >= nlat or lonPos < 0 or lonPos >= nlon:
            print 'latPos,lonPos,nlat,nlon',latPos,lonPos,nlat,nlon
            continue
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
              season[latPos, lonPos, 0] = 255
            elif dateTime.month < 7:
              season[latPos, lonPos, 1] = 255
            else:
              season[latPos, lonPos, 2] = 255
          else:
            season[latPos, lonPos, :] = 255
          if count[latPos, lonPos] == 32767:
            continue
          if fCO2Index is not None:
            fCO2i = float(line[fCO2Index])
            if fCO2i == MISSING:
              continue
            if fCO2i < fCO2min:
              fCO2min = fCO2i
            fCO2[latPos, lonPos] += fCO2i
            fCO2std[latPos, lonPos] += fCO2i**2
          if SSTIndex is not None:
            SST[latPos, lonPos] += 273.15 + float(line[SSTIndex])
          count[latPos, lonPos] += 1
          if count[latPos, lonPos] == 32767: # average the first 32767 values
            if fCO2Index is not None:
              fCO2std[latPos, lonPos] = np.sqrt((fCO2std[latPos, lonPos] - fCO2[latPos, lonPos] / count[latPos, lonPos]) / count[latPos, lonPos])
              fCO2[latPos, lonPos] /= count[latPos, lonPos]
            if SSTIndex is not None:
              SST[latPos, lonPos] /= count[latPos, lonPos]
      if fCO2min < 250.:
         print inFile, fCO2min

  filled = np.nonzero(count == 0)
  valid = np.nonzero(np.logical_and(count > 0,count < 32767))
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
  setattr(ncFile, 'institution', 'Plymouth Marine Laboratory Remote Sensing Group') 
  setattr(ncFile, 'contact1', 'email: rsghelp@pml.ac.uk') 
  setattr(ncFile, 'RSG_areacode', 'a7')
  setattr(ncFile, 'RSG_sensor', "NETCDF") 

print "(text2ncdf, main) SUCCESS writing file %s" % ncFileName
