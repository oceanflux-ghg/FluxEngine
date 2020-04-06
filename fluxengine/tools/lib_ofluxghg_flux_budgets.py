#!/usr/bin/env python2

'''utility to calculate yearly and monthly net fluxes in regions'''

# ofluxghg-flux-budgets.py
# utility to load various input netcdf datasets and determine air-sea flux of
# CO2 based on parameterisations in the ESA STSE OceanFlux Greenhouse Gases
# Technical Specification (TS) 


# History
# date, description, author, email, institution
# v0 08/2013  Basic structure and output, Peter Land, peland@pml.ac.uk, Plymouth Marine Laboratory.

# TO DO
# account for months in which there are no valid pixels in a region


from netCDF4 import Dataset
from glob import glob
import numpy as np, sys, csv, calendar
from numpy import zeros
from math import sqrt, cos, sin, radians
from scipy.integrate import quad

def run_flux_budgets(args):
    np.seterr(all = 'raise', invalid = 'ignore', under='ignore')#SOCATv4 - added under='ignore'
    NaN = float('nan')
    
    baseTitles = ['YEAR', 'MONTH', 'NET FLUX TgC', 'MISSING FLUX TgC',
       '30.5 DAY FLUX TgC', 'GROSS DOWNWARD FLUX TgC',
       'MISSING GROSS DOWNWARD FLUX TgC', '30.5 DAY GROSS DOWNWARD FLUX TgC',
       'GROSS UPWARD FLUX TgC', 'MISSING GROSS UPWARD FLUX TgC',
       '30.5 DAY GROSS UPWARD FLUX TgC']
    
     # unit conversion g->Tg, km^2->m^2 (day->month needs doing separately)
    unitConversion = 1e-12 * 1e6
     # kw conversion cm->m, /hr->/day
    kwConversion = 24. / 100.
    daysInMeanMonth = 30.5
    
     # Earth equatorial and polar radii in km
    re, rp = 6378.137, 6356.7523
    dtor = radians(1.)
    
    def checkShape(item):
       name, shape = item[0], item[1].shape
       if shape != inshape:
          print(name, ' data have a different shape: ', shape, inshape)
          sys.exit(1)
    
    def oneDegreeArea(latDegrees):
       # area of a 1x1 degree box at a given latitude in radians
       latRadians = latDegrees * dtor
       cosLat, sinLat = cos(latRadians), sin(latRadians)
       rc, rs = re * cosLat, rp * sinLat
       r2c, r2s = re * re * cosLat, rp * rp * sinLat
       earth_radius = sqrt((r2c * r2c + r2s * r2s) / (rc * rc + rs * rs))
       erd = earth_radius * dtor
       return erd * erd * cosLat
    
     # function name definition for stderr and stdout messages
    function = "(ofluxghg-flux-budgets, main)"
    
    
    
    indir, refdir, outroot = args.dir, args.ref, args.outroot
    fluxdataset, kwdataset = args.fluxdataset, args.kwdataset
    cidataset, cwdataset = args.cidataset, args.cwdataset
    icedataset = args.icedataset
    maskfile, maskdatasets = args.maskfile, args.maskdatasets
    landfile, landdataset = args.landfile, args.landdataset
    regions, LooseIce, icePercent = args.regions, args.LooseIce, args.icePercent
    window, places, verbosity = args.window, args.places, args.verbosity
    gridArea, gridAreafile, gridAreaName = args.gridarea, args.gridareafile, args.gridareadataset
    
    referencing = refdir is not None
    windowing = window is not None
    if windowing and len(window) != 4:
       print('window argument has length ', len(window), ' (should be 4)')
       sys.exit(1)
    mytype = '|S' + str(places)
    
    if indir[-1:] != '/':
       indir += '/'
    if outroot[-1:] != '_':
       outroot += '_'
    
    nRegions = len(regions) # sort out regions and masks
    nMask = len(maskdatasets)
    if nMask == 0:
       if nRegions == 0:
          print("Using default region 'global'")
       nMask, maskdatasets = nRegions, regions
    elif nRegions == 0:
       nRegions, regions = nMask, maskdatasets
    if nMask != nRegions:
       print('Number of regions does not match number of mask datasets')
       sys.exit(1)
    addedGlobal = not (('global' in regions) or (nRegions != 0 and windowing))
    if addedGlobal:
       if windowing:
          print("WARNING: 'global' data restricted to window")
       print('global output file ', outroot + 'global.txt')
       printGlobalArea = True
       globalLun = open(outroot + 'global.txt', 'w')
       globalWriter = csv.writer(globalLun)
       if referencing:
          print('global reference output file ', outroot + 'globalRef.txt')
          refGlobalLun = open(outroot + 'globalRef.txt', 'w')
          refGlobalWriter = csv.writer(refGlobalLun)
    luns = []
    writers = []
    if referencing:
       refLuns = []
       refWriters = []
    for region in regions:
       print(region + ' output file ', outroot + region + '.txt')
       lun = open(outroot + region + '.txt', 'w')
       luns.append(lun)
       writers.append(csv.writer(lun))
       if referencing:
          print(region + ' reference output file ', outroot + region + 'Ref.txt')
          refLun = open(outroot + region + 'Ref.txt', 'w')
          refLuns.append(lun)
          refWriters.append(csv.writer(lun))
    
    yeardirs = sorted(glob(indir + '[0-9][0-9][0-9][0-9]/'))
    nYears = len(yeardirs)
    if nYears == 0:
       print('No year directories found!')
    if referencing:
       refyeardirs = sorted(glob(refdir + '[0-9][0-9][0-9][0-9]/'))
       if len(refyeardirs) != nYears:
          print("Year directories don't match: ", nYears, len(refyeardirs))
          sys.exit(1)
    refFMask = False
    for yearNo in range(nYears):
       yeardir = yeardirs[yearNo]
       syear = yeardir[-5:-1]
       year = int(syear)
       monthdirs = sorted(glob(yeardir + '[0-9][0-9]/'))
       nMonths = len(monthdirs)
       if referencing:
          refyeardir = refyeardirs[yearNo]
          if refyeardir[-5:-1] != syear:
             print("Years don't match: ", syear, refyeardir[-5:-1])
             sys.exit(1)
          refmonthdirs = sorted(glob(refyeardir + '[0-9][0-9]/'))
          if len(refmonthdirs) != nMonths:
             print("Month directories don't match: ", nMonths, len(refmonthdirs))
             sys.exit(1)
       if nMonths == 0:
          print('No month directories found in ', year, ' - continuing')
          continue
       for monthNo in range(nMonths):
          monthdir = monthdirs[monthNo]
          smonth = monthdir[-3:-1]
          month = int(smonth)
          if referencing:
             refmonthdir = refmonthdirs[monthNo]
             if refmonthdir[-3:-1] != smonth:
                print("Months don't match: ", syear, smonth, refmonthdir[-3:-1])
                sys.exit(1)
          daysInMonth = calendar.monthrange(year, month)[1]
          ncFile = glob(monthdir + '*.nc')
          nFiles = len(ncFile)
          if referencing:
             refNcFile = glob(refmonthdir + '*.nc')
             if len(refNcFile) != nFiles:
                print("Files don't match: ", syear, smonth, nFiles, len(refNcFile))
                sys.exit(1)
          if nFiles == 0:
             print('No .nc files found in ', year, ',', month, ' - continuing')
             continue
          if nFiles != 1:
             print(nFiles, ' .nc files found in ', year, ',', month)
             sys.exit(1)
          try:
             ncdata = Dataset(ncFile[0], 'r')
          except:
             print('Unable to open ', ncFile[0])
             sys.exit(1)
          if referencing:
             try:
                refncdata = Dataset(refNcFile[0], 'r')
             except:
                print('Unable to open ', refNcFile[0])
                sys.exit(1)
          variables = ncdata.variables
          nVariables = len(variables)
          try:
             variable = variables[fluxdataset]
          except:
             print('Unable to open ', fluxdataset, ' dataset in ', ncFile[0])
             sys.exit(1)
          if referencing:
             refVariables = refncdata.variables
             if len(refVariables) != nVariables:
                print("Variables don't match: ", nVariables, len(refVariables))
                sys.exit(1)
             try:
                refVariable = refVariables[fluxdataset]
             except:
                print('Unable to open ', fluxdataset, ' dataset in ', refNcFile[0])
                sys.exit(1)
          if 'inshape' in locals():
             checkShape((fluxdataset, variable))
             if referencing:
                checkShape((fluxdataset, refVariable))
          else:
             inshape = variable.shape
             if inshape[0] != 1:
                print(fluxdataset, ' data have an odd shape: ', inshape)
                sys.exit(1)
             if referencing and refVariable.shape != inshape:
                print(fluxdataset, ' inconsistent data: ', inshape,
                   refVariable.shape)
                sys.exit(1)
             if windowing:
                y0, y1, x0, x1 = window
             else:
                y0 = x0 = 0
                y1, x1 = inshape[1], inshape[2]
             ny, nx = y1 - y0, x1 - x0
             if verbosity > 0:
                print('dimensions ', ny, nx, y0, y1, x0, x1)
             extraNames = []
             for item in variables.items(): # find extra variable names
                shape = item[1].shape
                if shape == inshape:
                   extraNames.append(item[0])
             nExtras = len(extraNames)
             extraData = np.ma.MaskedArray(np.empty((nExtras, ny, nx)), mask=False)
             extraData.fill(float('nan'))
             if referencing:
                refExtraData = np.ma.MaskedArray(np.empty((nExtras, ny, nx)),
                   mask=False)
                refExtraData.fill(float('nan'))
             titles = np.append(baseTitles, extraNames)
             #titles = np.array([title.encode() for title in titles]); #TMH: convert from str to byte array
             if addedGlobal:
                globalWriter.writerow(titles)
                if referencing:
                   refGlobalWriter.writerow(titles)
             for regionNo in range(nRegions):
                writers[regionNo].writerow(titles)
                if referencing:
                   refWriters[regionNo].writerow(titles)
             try:
                landDataset = Dataset(landfile, 'r')
             except:
                print('Unable to open land file ', landfile)
                sys.exit(1)
             try:
                land_variable = landDataset.variables[landdataset]
             except:
                print('Unable to open ', landdataset, ' dataset in ', landfile)
                sys.exit(1)
             checkShape((landdataset, land_variable))
             ocean = 1. - land_variable[0, y0:y1, x0:x1]
             landDataset.close()
             try:
                mask_dataset = Dataset(maskfile, 'r')
             except:
                print('Unable to open ', maskfile)
                sys.exit(1)
             for maskdataset in maskdatasets:
                try:
                   mask_variable = mask_dataset.variables[maskdataset]
                except:
                   try:
                      maskdataset = maskdataset.replace("_"," ")
                      mask_variable = mask_dataset.variables[maskdataset]
                   except:
                      print('Unable to open ', maskdataset, ' dataset in ', maskfile)
                      sys.exit(1)
                checkShape((maskdataset, mask_variable))
                if 'masks' in locals():
                   masks = np.ma.vstack((masks,
                      mask_variable[0, y0:y1, x0:x1].reshape(1, ny, nx)))
                else:
                   masks = mask_variable[0, y0:y1, x0:x1].reshape(1, ny, nx)
                # masks = np.fliplr(masks)Uncomment if using Europe200.nc mask or similar
             mask_dataset.close()
    
             #  # calculate area of a 1x1 degree box at a given latitude
             # cellAreas = np.empty(ny)
             # cellAreas.fill(float('nan'))
             # for j in xrange(y0,y1):
             #    cellAreas[j-y0] = quad(oneDegreeArea, 89. - j, 90. - j)[0]
    
    #IGA----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    #gridArea, gridAreafile, gridAreaName
             if gridArea == 0 and gridAreafile == 'no_file':
                print('No grid areas given - assuming a 1deg x 1deg grid and calculating net flux accordingly.')
                cellAreas = np.empty(ny)
                cellAreas.fill(float('nan'))
                for j in range(y0,y1):
                   cellAreas[j-y0] = quad(oneDegreeArea, 89. - j, 90. - j)[0]
             elif gridAreafile != 'no_file':
                print('Reading grid area from file', gridAreafile)
                with Dataset(gridAreafile, 'r') as ncArea:
                   cellAreas = ncArea.variables[gridAreaName][:]
                cellAreas = cellAreas[0,y0:y1, x0:x1]
             elif gridArea != 0 and gridAreafile == 'no_file':
                print('Reading grid area from integer input = ', gridArea)
                cellAreas = np.empty(ny)
                cellAreas.fill(float(gridArea))
    
                # if len(gridAreadataset) == 1:#Wrong if argument
                #    cellAreas = np.empty(ny)
                #    cellAreas.fill(float(gridAreadataset))
                # else:
                #    with Dataset(gridAreadataset,'r') as ga_in:
                #     cellAreas = ga_in['cell_area']
    #IGA-----------------------------------------------------------------------------
    
             if addedGlobal:
                yearlyGlobalFlux = np.array(0.) # to divide by zero without errors
                yearlyGlobalDownFlux = np.array(0.)
                yearlyGlobalUpFlux = np.array(0.)
                yearlyGlobalArea = np.array(0.)
                yearlyGlobalMissingArea = np.array(0.)
                yearlyGlobalExtras = np.zeros(nExtras)
                yearlyGlobalExtraAreas = np.zeros(nExtras)
                if referencing:
                   refYearlyGlobalFlux = np.array(0.)
                   refYearlyGlobalDownFlux = np.array(0.)
                   refYearlyGlobalUpFlux = np.array(0.)
                   refYearlyGlobalArea = np.array(0.)
                   refYearlyGlobalExtras = np.zeros(nExtras)
                   refYearlyGlobalExtraAreas = np.zeros(nExtras)
             yearlyFluxes = np.zeros(nRegions)
             yearlyDownFluxes = np.zeros(nRegions)
             yearlyUpFluxes = np.zeros(nRegions)
             yearlyAreas = np.zeros(nRegions)
             yearlyMissingAreas = np.zeros(nRegions)
             yearlyExtras = np.zeros((nExtras,nRegions))
             yearlyExtraAreas = np.zeros((nExtras,nRegions))
             if referencing:
                refYearlyFluxes = np.zeros(nRegions)
                refYearlyDownFluxes = np.zeros(nRegions)
                refYearlyUpFluxes = np.zeros(nRegions)
                refYearlyAreas = np.zeros(nRegions)
                refYearlyExtras = np.zeros((nExtras,nRegions))
                refYearlyExtraAreas = np.zeros((nExtras,nRegions))
          checkShape((fluxdataset, variable))
          flux = variable[0, :, :]
          fluxMask = flux.mask[y0:y1, x0:x1]
          flux = flux[y0:y1, x0:x1]
          if referencing:
             checkShape((fluxdataset, refVariable))
             refFlux = refVariable[0, :, :]
             refFluxMask = refFlux.mask[y0:y1, x0:x1]
             refFlux = refFlux[y0:y1, x0:x1]
          try:
             variable = variables[kwdataset]
          except:
             print('Unable to open ', kwdataset, ' dataset in ', ncFile[0])
             sys.exit(1)
          checkShape((kwdataset, variable))
          kw = variable[0, y0:y1, x0:x1]
          if referencing:
             try:
                refVariable = refVariables[kwdataset]
             except:
                print('Unable to open ', kwdataset, ' dataset in ', refNcFile[0])
                sys.exit(1)
             checkShape((kwdataset, refVariable))
             refKw = refVariable[0, y0:y1, x0:x1]
          try:
             variable = variables[cidataset]
          except:
             print('Unable to open ', cidataset, ' dataset in ', ncFile[0])
             sys.exit(1)
          checkShape((cidataset, variable))
          ci = variable[0, y0:y1, x0:x1]
          if referencing:
             try:
                refVariable = refVariables[cidataset]
             except:
                print('Unable to open ', cidataset, ' dataset in ', refNcFile[0])
                sys.exit(1)
             checkShape((cidataset, refVariable))
             refCi = refVariable[0, y0:y1, x0:x1]
          try:
             variable = variables[cwdataset]
          except:
             print('Unable to open ', cwdataset, ' dataset in ', ncFile[0])
             sys.exit(1)
          checkShape((cwdataset, variable))
          cw = variable[0, y0:y1, x0:x1]
          if referencing:
             try:
                refVariable = refVariables[cwdataset]
             except:
                print('Unable to open ', cwdataset, ' dataset in ', refNcFile[0])
                sys.exit(1)
             checkShape((cwdataset, refVariable))
             refCw = refVariable[0, y0:y1, x0:x1]
          try:
             variable = variables[icedataset]
          except:
             print('Unable to open ', icedataset, ' dataset in ', ncFile[0])
             sys.exit(1)
          checkShape((fluxdataset, variable))
          ice = variable[0, :, :]
          
          try:
              if isinstance(ice.mask, np.bool_): #If the mask is a single boolean value, convert it to a filled matrix of that value (otherwise window subsetting does not work)
                  ice.mask = np.full(ice.data.shape, np.bool_(ice.mask));
              iceMask = ice.mask[y0:y1, x0:x1]
          except AttributeError:
             ice = ice[:].view(np.ma.MaskedArray)
             ice.mask = zeros(ice.shape)
             iceMask = ice.mask[y0:y1, x0:x1]
    
          if icePercent:
             ice /= 100. # NB -999 now becomes -99.9, shouldn't matter though
          ice = ice[y0:y1, x0:x1]
          iceMask[ice<0.]=True
          if referencing:
             try:
                refVariable = refVariables[icedataset]
             except:
                print('Unable to open ', icedataset, ' dataset in ', refNcFile[0])
                sys.exit(1)
             checkShape((fluxdataset, refVariable))
             refIce = refVariable[0, y0:y1, x0:x1]
          for extraNo in range(nExtras):
             extraName = extraNames[extraNo]
             try:
                variable = variables[extraName]
             except:
                print('Unable to open ', extraName, ' dataset in ', ncFile[0])
                sys.exit(1)
             checkShape((fluxdataset, variable))
             extraData[extraNo, :, :] = variable[0, y0:y1, x0:x1]
             if referencing:
                try:
                   refVariable = refVariables[extraName]
                except:
                   print('Unable to open ', extraName, ' dataset in ', refNcFile[0])
                   sys.exit(1)
                checkShape((fluxdataset, refVariable))
                refExtraData[extraNo, :, :] = refVariable[0, y0:y1, x0:x1]
          ncdata.close()
          if referencing:
             refncdata.close()
          if addedGlobal:
             globalFlux = np.array(0.) # so you can divide by zero without errors
             globalDownFlux = np.array(0.)
             globalUpFlux = np.array(0.)
             globalArea = np.array(0.)
             globalMissingArea = np.array(0.)
             globalExtras = np.zeros(nExtras)
             globalExtraAreas = np.zeros(nExtras)
             if referencing:
                refGlobalFlux = np.array(0.)
                refGlobalDownFlux = np.array(0.)
                refGlobalUpFlux = np.array(0.)
                refGlobalExtras = np.zeros(nExtras)
                refGlobalExtraAreas = np.zeros(nExtras)
          fluxes = np.zeros(nRegions)
          downFluxes = np.zeros(nRegions)
          upFluxes = np.zeros(nRegions)
          areas = np.zeros(nRegions)
          missingAreas = np.zeros(nRegions)
          extras = np.zeros((nExtras,nRegions))
          extraAreas = np.zeros((nExtras,nRegions))
          if referencing:
             refFluxes = np.zeros(nRegions)
             refDownFluxes = np.zeros(nRegions)
             refUpFluxes = np.zeros(nRegions)
             refExtras = np.zeros((nExtras,nRegions))
             refExtraAreas = np.zeros((nExtras,nRegions))
          if verbosity > 0:
             print('Processing ', month, ' ', year)
          for j in range(ny):
             for i in range(nx):
                if ocean[j, i] == 0.:
                   continue
                f = flux[j, i]
                if abs(f) < 1e-300:
                   f = 0.
                if verbosity > 1 and not (fluxMask[j, i] or np.isfinite(f)):
                   print(month, year, i, j, 'flux', f)
                if len(cellAreas.shape)>1:
                   tCellArea = cellAreas[j,i]
                else:
                   tCellArea = cellAreas[j]#If the cellAreas are regular, they have only one dimension. Otherwise, they are x and y varying and picked accordingly
                if iceMask[j, i]:
                   iceCover = 0.
                else:
                   iceCover = ice[j, i]
                if verbosity > 1 and not (fluxMask[j, i] or np.isfinite(iceCover)):
                   print(month, year, i, j, 'ice', iceCover)
                if not (fluxMask[j, i]) and (iceCover < 0. or iceCover > 1.):
                   print(month, year, i, j, 'ice', iceCover)
                if referencing and refIce[j,i] != iceCover:
                   print("ice cover doesn't match:", j, i, iceCover, refIce[j,i])
                   sys.exit(1)
                if iceCover > .1: # Takahashi neglects ice<10%
                   icefreeArea = max([1. - iceCover, .1]) # Takahashi icefreeArea>=.1
                   if LooseIce:
                      icefreeArea = pow(icefreeArea, .4)
                   icefreeArea *= tCellArea
                else: # this includes missing values, which are -999
                   icefreeArea = tCellArea
                try:
                   ag = icefreeArea * ocean[j, i]
                except ArithmeticError:
                   print(month, year, i, j, 'ag overflow', iceCover, icefreeArea, ocean[j, i])
                if nRegions != 0:
                   mask = masks[:, j, i]
                   try:
                      regionNos = np.arange(nRegions)[
                          np.logical_and(~mask.mask, mask > 0.)]
                   except:#IGA_netfluxBugTestJun2017
                      mask = np.ma.MaskedArray(mask)
                      regionNos = np.arange(nRegions)[
                          np.logical_and(~mask.mask, mask > 0.)]
                else:
                   regionNos = []
                if not (fluxMask[j, i]):
                   try:
                      fAg = f * ag
                   except ArithmeticError:
                      print(month, year, i, j, 'fAg overflow', f, iceCover, icefreeArea, ag)
                kAg = kw[j, i] * kwConversion * ag
                downFAg = kAg * ci[j, i]
                upFAg = kAg * cw[j,i]
                es = extraData[:, j, i]
                if not (fluxMask[j, i]):
                   for extraNo in np.arange(nExtras)[~es.mask]:
                      if abs(es[extraNo]) < 1e-300:
                         es[extraNo] = 0.
                      try:
                         eAg = float(es[extraNo]) * ag
                      except ArithmeticError:
                         print(month, year, i, j, 'eAg overflow', f, iceCover, icefreeArea, ag)
                      if np.isnan(eAg) == False: #TMH: Only sum if it's not nan.
                          if addedGlobal:
                             globalExtras[extraNo] += eAg
                             globalExtraAreas[extraNo] += ag
                          for regionNo in regionNos:
                             extras[extraNo, regionNo] += eAg
                             extraAreas[extraNo, regionNo] += ag
                if referencing:
                   refFAg = refFlux[j, i] * ag
                   refDownFAg = refKw[j, i] * refCi[j, i] * ag
                   refUpFAg = refKw[j, i] * refCw[j, i] * ag
                   refFMask = refFluxMask[j, i]
                   refEs = refExtraData[:, j, i]
                   for extraNo in np.arange(nExtras)[~refEs.mask]:
                      e = float(refEs[extraNo])
                      if addedGlobal:
                         refGlobalExtras[extraNo] += eAg
                         refGlobalExtraAreas[extraNo] += ag
                      for regionNo in regionNos:
                         refExtras[extraNo, regionNo] += eAg
                         refExtraAreas[extraNo, regionNo] += ag
                if fluxMask[j,i] or refFMask:
                   if addedGlobal:
                      globalMissingArea += ag
                   for regionNo in regionNos:
                      missingAreas[regionNo] = min(icefreeArea * mask[regionNo], ag)
                   continue
                #TMH: added if clause to prevent adding of nans
                if np.isnan(fAg)==False:
                    if addedGlobal:
                       globalFlux += fAg
                       globalDownFlux += downFAg
                       globalUpFlux += upFAg
                       globalArea += ag
                       if referencing:
                          refGlobalFlux += refFAg
                          refGlobalDownFlux += refDownFAg
                          refGlobalUpFlux += refUpFAg
                    for regionNo in regionNos:
                       a = min(icefreeArea * mask[regionNo], ag)
                       fluxes[regionNo] += fAg
                       downFluxes[regionNo] += downFAg
                       upFluxes[regionNo] += upFAg
                       areas[regionNo] += ag
                       if referencing:
                          refFluxes[regionNo] += refF * ag
                          refDownFluxes[regionNo] += refDownFAg
                          refUpFluxes[regionNo] += refUpFAg
          if addedGlobal:
             if printGlobalArea:
                print('Global area: ', globalArea, globalMissingArea,
                   globalArea + globalMissingArea)
                printGlobalArea = False
             yearlyGlobalFlux += globalFlux
             yearlyGlobalDownFlux += globalDownFlux
             yearlyGlobalUpFlux += globalUpFlux
             yearlyGlobalArea += globalArea
             yearlyGlobalMissingArea += globalMissingArea
             yearlyGlobalExtras += globalExtras
             yearlyGlobalExtraAreas += globalExtraAreas
             if referencing:
                refYearlyGlobalFlux += refGlobalFlux
                refYearlyGlobalDownFlux += refGlobalDownFlux
                refYearlyGlobalUpFlux += refGlobalUpFlux
                refYearlyGlobalExtras += refGlobalExtras
                refYearlyGlobalExtraAreas += refGlobalExtraAreas
          yearlyFluxes += fluxes
          yearlyDownFluxes += downFluxes
          yearlyUpFluxes += upFluxes
          yearlyAreas += areas
          yearlyMissingAreas += missingAreas
          yearlyExtras += extras
          yearlyExtraAreas += extraAreas
          if referencing:
             refYearlyFluxes += refFluxes
             refYearlyDownFluxes += refDownFluxes
             refYearlyUpFluxes += refUpFluxes
             refYearlyExtras += refExtras
             refYearlyExtraAreas += refExtraAreas
          if addedGlobal:
             globalFlux *= unitConversion # convert to TgC/day
             globalDownFlux *= unitConversion
             globalUpFlux *= unitConversion
             if referencing:
                refGlobalFlux *= unitConversion
                refGlobalDownFlux *= unitConversion
                refGlobalUpFlux *= unitConversion
             if globalMissingArea != 0.:
                meanFlux = globalFlux / globalArea # in TgC/km2/day
                meanDownFlux = globalDownFlux / globalArea
                meanUpFlux = globalUpFlux / globalArea
                print('Mean flux: ',
                  np.array([meanFlux, meanDownFlux, meanUpFlux]) / unitConversion)
                missingFlux = meanFlux * globalMissingArea
                missingDownFlux = meanDownFlux * globalMissingArea
                missingUpFlux = meanUpFlux * globalMissingArea
                globalFlux += missingFlux
                globalDownFlux += missingDownFlux
                globalUpFlux += missingUpFlux
                if referencing:
                   meanFlux = refGlobalFlux / globalArea
                   meanDownFlux = refGlobalDownFlux / globalArea
                   meanUpFlux = refGlobalUpFlux / globalArea
                   refMissingFlux = meanFlux * globalMissingArea
                   refMissingDownFlux = meanDownFlux * globalMissingArea
                   refMissingUpFlux = meanUpFlux * globalMissingArea
                   refGlobalFlux += refMissingFlux
                   refGlobalDownFlux += refMissingDownFlux
                   refGlobalUpFlux += refMissingUpFlux
             else:
                missingFlux = missingDownFlux = missingUpFlux = 0.
                if referencing:
                   refMissingFlux = refMissingDownFlux = refMissingUpFlux = 0.
             standardFlux = globalFlux * daysInMeanMonth # in TgC/month
             standardDownFlux = globalDownFlux * daysInMeanMonth
             standardUpFlux = globalUpFlux * daysInMeanMonth
             globalFlux *= daysInMonth
             globalDownFlux *= daysInMonth
             globalUpFlux *= daysInMonth
             missingFlux *= daysInMonth
             missingDownFlux *= daysInMonth
             missingUpFlux *= daysInMonth
             globalWriter.writerow(np.append([syear, smonth], np.append([globalFlux,
                missingFlux, standardFlux, globalDownFlux, missingDownFlux,
                standardDownFlux, globalUpFlux, missingUpFlux, standardUpFlux],
                np.divide(globalExtras, globalExtraAreas)).astype(mytype)))
             if referencing:
                standardFlux = refGlobalFlux * daysInMeanMonth
                standardDownFlux = refGlobalDownFlux * daysInMeanMonth
                standardUpFlux = refGlobalUpFlux * daysInMeanMonth
                refGlobalFlux *= daysInMonth
                refGlobalDownFlux *= daysInMonth
                refGlobalUpFlux *= daysInMonth
                refMissingFlux *= daysInMonth
                refMissingDownFlux *= daysInMonth
                refMissingUpFlux *= daysInMonth
                refGlobalWriter.writerow(np.append([syear, smonth],
                   np.append([refGlobalFlux, refMissingFlux, standardFlux,
                   refGlobalDownFlux, refMissingDownFlux, standardDownFlux,
                   refGlobalUpFlux, refMissingUpFlux, standardUpFlux],
                   np.divide(refGlobalExtras, refGlobalExtraAreas)).astype(mytype)))
          if nRegions == 0:
             continue
          fluxes *= unitConversion
          downFluxes *= unitConversion
          upFluxes *= unitConversion
          meanFluxes = fluxes / areas
          meanDownFluxes = downFluxes / areas
          meanUpFluxes = upFluxes / areas
          if np.any(missingAreas != 0.):
             missingFluxes = meanFluxes * missingAreas
             missingDownFluxes = meanDownFluxes * missingAreas
             missingUpFluxes = meanUpFluxes * missingAreas
             fluxes += missingFluxes
             downFluxes += missingDownFluxes
             upFluxes += missingUpFluxes
          else:
             print("year, month, missing areas ", year, month, missingAreas)
             missingFluxes = np.zeros(nRegions)
             missingDownFluxes = np.zeros(nRegions)
             missingUpFluxes = np.zeros(nRegions)
          standardFluxes = fluxes * daysInMeanMonth # flux from a standard month
          standardDownFluxes = downFluxes * daysInMeanMonth
          standardUpFluxes = upFluxes * daysInMeanMonth
          fluxes *= daysInMonth # flux from the actual month (less in February)
          downFluxes *= daysInMonth
          upFluxes *= daysInMonth
          missingFluxes *= daysInMonth
          missingDownFluxes *= daysInMonth
          missingUpFluxes *= daysInMonth
          extras = np.divide(extras, extraAreas)
          if referencing:
             refFluxes *= unitConversion
             refDownFluxes *= unitConversion
             refUpFluxes *= unitConversion
             meanFluxes = refFluxes / areas
             meanDownFluxes = refDownFluxes / areas
             meanUpFluxes = refUpFluxes / areas
             if np.any(missingAreas != 0.):
                refMissingFluxes = meanFluxes * missingAreas
                refMissingDownFluxes = meanDownFluxes * missingAreas
                refMissingUpFluxes = meanUpFluxes * missingAreas
                refFluxes += refMissingFluxes
                refDownFluxes += refMissingDownFluxes
                refUpFluxes += refMissingUpFluxes
             else:
                refMissingFluxes = np.zeros(nRegions)
                refMissingDownFluxes = np.zeros(nRegions)
                refMissingUpFluxes = np.zeros(nRegions)
             refStandardFluxes = refFluxes * daysInMeanMonth
             refStandardDownFluxes = refDownFluxes * daysInMeanMonth
             refStandardUpFluxes = refUpFluxes * daysInMeanMonth
             refFluxes *= daysInMonth
             refDownFluxes *= daysInMonth
             refUpFluxes *= daysInMonth
             refMissingFluxes *= daysInMonth
             refMissingDownFluxes *= daysInMonth
             refMissingUpFluxes *= daysInMonth
             refExtras = np.divide(refExtras, refExtraAreas)
          for regionNo in range(nRegions):
             writers[regionNo].writerow(np.append([syear, smonth],
                np.append([fluxes[regionNo], missingFluxes[regionNo],
                standardFluxes[regionNo], downFluxes[regionNo],
                missingDownFluxes[regionNo], standardDownFluxes[regionNo],
                upFluxes[regionNo], missingUpFluxes[regionNo],
                standardUpFluxes[regionNo]], extras[:, regionNo]).astype(mytype)))
             if referencing:
                refWriters[regionNo].writerow(np.append([syear, smonth],
                np.append([refFluxes[regionNo], refMissingFluxes[regionNo],
                refStandardFluxes[regionNo], refDownFluxes[regionNo],
                refMissingDownFluxes[regionNo], refStandardDownFluxes[regionNo],
                refUpFluxes[regionNo], refMissingUpFluxes[regionNo],
                refStandardUpFluxes[regionNo]],
                refExtras[:, regionNo]).astype(mytype)))
             print('Region number and area',regionNo, areas[regionNo])
       if addedGlobal:
          yearlyGlobalFlux *= unitConversion
          yearlyGlobalDownFlux *= unitConversion
          yearlyGlobalUpFlux *= unitConversion
          if referencing:
             refYearlyGlobalFlux *= unitConversion
             refYearlyGlobalDownFlux *= unitConversion
             refYearlyGlobalUpFlux *= unitConversion
          if yearlyGlobalMissingArea != 0.:
             meanFlux = yearlyGlobalFlux / yearlyGlobalArea
             meanDownFlux = yearlyGlobalDownFlux / yearlyGlobalArea
             meanUpFlux = yearlyGlobalUpFlux / yearlyGlobalArea
             print('Yearly mean flux ('+str(year)+"):", np.array([meanFlux, meanDownFlux, meanUpFlux]) / unitConversion)
             missingFlux = meanFlux * yearlyGlobalMissingArea
             missingDownFlux = meanDownFlux * yearlyGlobalMissingArea
             missingUpFlux = meanUpFlux * yearlyGlobalMissingArea
             yearlyGlobalFlux += missingFlux
             yearlyGlobalDownFlux += missingDownFlux
             yearlyGlobalUpFlux += missingUpFlux
             if referencing:
                meanFlux = refYearlyGlobalFlux / yearlyGlobalArea
                meanDownFlux = refYearlyGlobalDownFlux / yearlyGlobalArea
                meanUpFlux = refYearlyGlobalUpFlux / yearlyGlobalArea
                refMissingFlux = meanFlux * yearlyGlobalMissingArea
                refMissingDownFlux = meanDownFlux * yearlyGlobalMissingArea
                refMissingUpFlux = meanUpFlux * yearlyGlobalMissingArea
                refYearlyGlobalFlux += refMissingFlux
                refYearlyGlobalDownFlux += refMissingDownFlux
                refYearlyGlobalUpFlux += refMissingUpFlux
          else:
             missingFlux = missingDownFlux = missingUpFlux = 0.
             if referencing:
                refMissingFlux = refMissingDownFlux = refMissingUpFlux = 0.
          standardGlobalFlux = yearlyGlobalFlux * daysInMeanMonth
          standardGlobalDownFlux = yearlyGlobalDownFlux * daysInMeanMonth
          standardGlobalUpFlux = yearlyGlobalUpFlux * daysInMeanMonth
          yearlyGlobalFlux *= daysInMeanMonth
          yearlyGlobalDownFlux *= daysInMeanMonth
          yearlyGlobalUpFlux *= daysInMeanMonth
          missingFlux *= daysInMeanMonth
          missingDownFlux *= daysInMeanMonth
          missingUpFlux *= daysInMeanMonth
          globalWriter.writerow(np.append([syear, 'ALL'],
             np.append([yearlyGlobalFlux, missingFlux, standardGlobalFlux,
             yearlyGlobalDownFlux, missingDownFlux, standardGlobalDownFlux,
             yearlyGlobalUpFlux, missingUpFlux, standardGlobalUpFlux],
             yearlyGlobalExtras / yearlyGlobalExtraAreas).astype(mytype)))
          yearlyGlobalFlux = np.array(0.) # to divide by zero without errors
          yearlyGlobalDownFlux = np.array(0.)
          yearlyGlobalUpFlux = np.array(0.)
          yearlyGlobalArea = np.array(0.)
          yearlyGlobalMissingArea = np.array(0.)
          yearlyGlobalExtras = np.zeros(nExtras)
          yearlyGlobalExtraAreas = np.zeros(nExtras)
          if referencing:
             refStandardGlobalFlux = refYearlyGlobalFlux * daysInMeanMonth
             refStandardGlobalDownFlux = refYearlyGlobalDownFlux * daysInMeanMonth
             refStandardGlobalUpFlux = refYearlyGlobalUpFlux * daysInMeanMonth
             refYearlyGlobalFlux *= daysInMeanMonth
             refYearlyGlobalDownFlux *= daysInMeanMonth
             refYearlyGlobalUpFlux *= daysInMeanMonth
             refMissingFlux *= daysInMeanMonth
             refMissingDownFlux *= daysInMeanMonth
             refMissingUpFlux *= daysInMeanMonth
             refGlobalWriter.writerow(np.append([syear, 'ALL'],
                np.append([refYearlyGlobalFlux, refMissingFlux,
                refStandardGlobalFlux, refYearlyGlobalDownFlux, refMissingDownFlux,
                refStandardGlobalDownFlux, refYearlyGlobalUpFlux, refMissingUpFlux,
                refStandardGlobalUpFlux],
                refYearlyGlobalExtras / refYearlyGlobalExtraAreas).astype(mytype)))
             refYearlyGlobalFlux = np.array(0.) # to divide by zero without errors
             refYearlyGlobalDownFlux = np.array(0.)
             refYearlyGlobalUpFlux = np.array(0.)
             refYearlyGlobalArea = np.array(0.)
             refYearlyGlobalMissingArea = np.array(0.)
             refYearlyGlobalExtras = np.zeros(nExtras)
             refYearlyGlobalExtraAreas = np.zeros(nExtras)
       if nRegions == 0:
          continue
       yearlyFluxes *= unitConversion
       yearlyDownFluxes *= unitConversion
       yearlyUpFluxes *= unitConversion
       meanFluxes = yearlyFluxes / yearlyAreas
       meanDownFluxes = yearlyDownFluxes / yearlyAreas
       meanUpFluxes = yearlyUpFluxes / yearlyAreas
       if np.any(yearlyMissingAreas != 0.):
          missingFluxes = meanFluxes * yearlyMissingAreas
          missingDownFluxes = meanDownFluxes * yearlyMissingAreas
          missingUpFluxes = meanUpFluxes * yearlyMissingAreas
          yearlyFluxes += missingFluxes
          yearlyDownFluxes += missingDownFluxes
          yearlyUpFluxes += missingUpFluxes
       else:
          missingFluxes = missingDownFluxes = missingUpFluxes = np.zeros(nRegions)
       standardYearlyFluxes = yearlyFluxes * daysInMeanMonth
       standardYearlyDownFluxes = yearlyDownFluxes * daysInMeanMonth
       standardYearlyUpFluxes = yearlyUpFluxes * daysInMeanMonth
       yearlyFluxes *= daysInMeanMonth
       yearlyDownFluxes *= daysInMeanMonth
       yearlyUpFluxes *= daysInMeanMonth
       yearlyExtras /= yearlyExtraAreas
       if referencing:
          refYearlyFluxes *= unitConversion
          refYearlyDownFluxes *= unitConversion
          refYearlyUpFluxes *= unitConversion
          meanFluxes = refYearlyFluxes / yearlyAreas
          meanDownFluxes = refYearlyDownFluxes / yearlyAreas
          meanUpFluxes = refYearlyUpFluxes / yearlyAreas
          if np.any(yearlyMissingAreas != 0.):
             refMissingFluxes = meanFluxes * yearlyMissingAreas
             refMissingDownFluxes = meanDownFluxes * yearlyMissingAreas
             refMissingUpFluxes = meanUpFluxes * yearlyMissingAreas
             refYearlyFluxes += refMissingFluxes
             refYearlyDownFluxes += refMissingDownFluxes
             refYearlyUpFluxes += refMissingUpFluxes
          else:
             refMissingFluxes = np.zeros(nRegions)
             refMissingDownFluxes = refMissingUpFluxes = np.zeros(nRegions)
          refStandardYearlyFluxes = refYearlyFluxes * daysInMeanMonth
          refStandardYearlyDownFluxes = refYearlyDownFluxes * daysInMeanMonth
          refStandardYearlyUpFluxes = refYearlyUpFluxes * daysInMeanMonth
          refYearlyFluxes *= daysInMeanMonth
          refYearlyDownFluxes *= daysInMeanMonth
          refYearlyUpFluxes *= daysInMeanMonth
          refYearlyExtras /= refYearlyExtraAreas
       for regionNo in range(nRegions):
          writers[regionNo].writerow(np.append([syear, 'ALL'],
             np.append([yearlyFluxes[regionNo], missingFluxes[regionNo],
             standardYearlyFluxes[regionNo], yearlyDownFluxes[regionNo],
             missingDownFluxes[regionNo], standardYearlyDownFluxes[regionNo],
             yearlyUpFluxes[regionNo], missingUpFluxes[regionNo],
             standardYearlyUpFluxes[regionNo]],
             yearlyExtras[:, regionNo]).astype(mytype)))
          if referencing:
             writers[regionNo].writerow(np.append([syear, 'ALL'],
                np.append([refYearlyFluxes[regionNo], refMissingFluxes[regionNo],
                refStandardYearlyFluxes[regionNo], refYearlyDownFluxes[regionNo],
                refMissingDownFluxes[regionNo],
                refStandardYearlyDownFluxes[regionNo], refYearlyUpFluxes[regionNo],
                refMissingUpFluxes[regionNo], refStandardYearlyUpFluxes[regionNo]],
                refYearlyExtras[:, regionNo]).astype(mytype)))
       yearlyFluxes = np.zeros(nRegions)
       yearlyDownFluxes = np.zeros(nRegions)
       yearlyUpFluxes = np.zeros(nRegions)
       yearlyAreas = np.zeros(nRegions)
       yearlyMissingAreas = np.zeros(nRegions)
       yearlyExtras = np.zeros((nExtras,nRegions))
       yearlyExtraAreas = np.zeros((nExtras,nRegions))
       if referencing:
          refYearlyFluxes = np.zeros(nRegions)
          refYearlyDownFluxes = np.zeros(nRegions)
          refYearlyUpFluxes = np.zeros(nRegions)
          refYearlyExtras = np.zeros((nExtras,nRegions))
          refYearlyExtraAreas = np.zeros((nExtras,nRegions))
    if addedGlobal:
       globalLun.close()
    for lun in luns:
       lun.close()
    if referencing:
       if addedGlobal:
          refGlobalLun.close()
       for lun in refLuns:
          lun.close()


