#!/home/cercache/tools/environments/run_python_scibox_precise.sh
#-*- coding:utf-8 -*-

# Script to convert the SOCAT(v2 and v3) files into netCDF monthly grids

#Standard + 3rd Party libs
import os
import numpy
import numpy.lib.recfunctions
import datetime
import argparse
import netCDF4
import glob
import multiprocessing
import pandas as pd;

#Internal tool libs
from . import get_sst
from . import datenum
from . import v2_f_conversion
from . import netcdf_helper
from . import combine_nc_files

#These are command line defaults. Set here (rather than in command line) as they can
#then be used as defaults when called from another script instead of command line
cldefaults={'start':1991,
            'end':2010,
            'prefix':None,
            'outputdir':'.',
            'sstdir':"/home/cerdata/provider/oceanflux/composites/sea_surface_temperature/global/arc/atsr/data/",
            'ssttail':"01_OCF-SST-GLO-1M-100-ATS-ARC.nc",
            'notperyear':False,
            'extrapolatetoyear':None,
            'socatversion':2,
            'asciioutput':False,
            'percruisedir':None,
            'coastalfile':None,
            'useaatsr': False,
            'usereynolds': False,
            'keepduplicates' : False,
            }

#list of variables for which we want std,min,max in outputfiles
statvariables=['fCO2_Tym','pCO2_Tym','fCO2_SST','pCO2_SST']
coastaldata=None


def GetCommandline():
   """
   Deal with command line parameters (if script is run from the command line)
   """
   #Deal with the command line arguments with some extra information and error checking
   parser = argparse.ArgumentParser(description="Convert SOCAT v2 or v3 files into netCDF files.")
   parser.add_argument('--start',type=int, metavar='<int>',help ='The year to start from for conversion calculation',default=cldefaults['start'])
   parser.add_argument('--end', type=int,metavar='<int>',help ='The year to end for conversion calculation',default=cldefaults['end'])
   parser.add_argument('--prefix', type=str,metavar='<string>',help ='A prefix to use for the output netCDF files.',default=cldefaults['prefix'])
   parser.add_argument('--outputdir',type=str, metavar='<path>',help ='A directory to write the output files to.',default=cldefaults['outputdir'])
   parser.add_argument('--percruisedir',type=str, metavar='<path>',help ='A directory to write the per-cruise output files to. If nothing given will not do cruise weighted averaging.',default=cldefaults['outputdir'])
   parser.add_argument('--sstdir', type=str,metavar='<path>',help ='A directory to read the SST data from',default=cldefaults['sstdir'])
   parser.add_argument('--ssttail',type=str, metavar='<path>',help ='The SST file prefix (after yyyymm).',default=cldefaults['ssttail'])
   parser.add_argument('--notperyear',dest='notperyear',action='store_true',help="Use this if you don't want separate files per year",default=cldefaults['notperyear'])
   parser.add_argument('--extrapolatetoyear',type=int,dest='extrapolatetoyear',help="Extrapolate to given year using Takahashi trend.",default=cldefaults['extrapolatetoyear'])
   parser.add_argument('--socatversion',type=int,dest='socatversion',help="The version of SOCAT data files to read.",default=cldefaults['socatversion'])
   parser.add_argument('--useaatsr',type=int,dest='useaatsr',help="To use the AATSR SST data.",default=cldefaults['useaatsr'])
   parser.add_argument('--usereynolds',type=int,dest='usereynolds',help="To use the Reynolds SST data.",default=cldefaults['usereynolds'])
   parser.add_argument('inputfile',metavar='<filename>',help ='The input SOCAT file to use')
   parser.add_argument('--no-grid-output',action='store_true',dest='asciioutput',help="To output as an ascii list rather than gridded netcdf.",default=cldefaults['asciioutput'])
   parser.add_argument('--coastalfile', type=str,metavar='<path>',help ='The file with the SOCAT coastal data in',default=cldefaults['coastalfile'])
   parser.add_argument('--keepduplicates',action='store_true',metavar='boolean',help ='Whether to keep duplicate data',default=cldefaults['keepduplicates'])
   commandline=parser.parse_args()

   #expand user incase path has been entered in style '~/'
   commandline.outputdir=os.path.expanduser(commandline.outputdir)

   #Set the output prefix to be the input filename
   if commandline.prefix is None:
      commandline.prefix=os.path.basename(commandline.inputfile)

   return commandline

#This wrapper function essentially does everything. If this script is used as a library
#rather than run on the command line then this is the function to call
#TMH: noConversion, if set to True, bypasses any reanalysis and instead writes raw SOCAT data to the output files.
def DoConversion(inputfile, columnInfo, startyr=cldefaults['start'],endyr=cldefaults['end'],prefix=cldefaults['prefix'],
                  outputdir=cldefaults['outputdir'],notperyear=cldefaults['notperyear'],sstdir=cldefaults['sstdir'],
                  ssttail=cldefaults['ssttail'],extrapolatetoyear=cldefaults['extrapolatetoyear'],
                  ASCIIOUT=cldefaults['asciioutput'], socatversion=6,
                  percruisedir=cldefaults['asciioutput'],coastalfile=cldefaults['coastalfile'],
                  useaatsr=cldefaults['useaatsr'],usereynolds=cldefaults['usereynolds'],
                  removeduplicates=not cldefaults['keepduplicates']):
   """
   Does all the hard work - if run as a library then call this function.
      inputfile - filename of the SOCAT ascii csv file
      startyr - start year of the range to process
      endyr - final year of the range to process
      prefix - prefix of the output file name to write
      outputdir - directory to write the output files to
      notperyear - if true then convert the range together, else if false convert years in the range individually
      sstdir - directory that contains the SST monthly climatology files
      ssttail - remainder of the SST netcdf name after the year and month
      extrapolatetoyear - year to extrapolate data to (or None)
      version - version of SOCAT data to use (2 or 3)
      ASCIIOUT - to write out as ASCII rather than netCDF
      percruisedir - directory to write per cruise files to
      coastalfile - file that contains the SOCAT coastal data
   """
   year_ranges_to_process=GetYearsToProcess(startyr,endyr,notperyear)
   #Read in the data (including coastal if requested)
   data=ReadInData(inputfile=inputfile,columnInfo=columnInfo, socatversion=socatversion)
   #read in the coastal data but use a global variable to
   #store the data so that this is only done once and not 
   #per run of function
   global coastaldata
   if coastalfile is not None and coastaldata is None:
      coastaldata=ReadInData(inputfile=coastalfile,columnInfo=columnInfo, socatversion=socatversion)

   #The coastal data needs to be handled differently to other data
   if coastaldata is not None:
      print("Updating coastal data for this region.")
      #If there are coastaldata then combine that to the data array
      #But only for cruises that exist in the data array and only for data that fall in a cell already occupied
      #These cruise data then need to be removed from the coastal data so that they are not added twice.
      #Get the lats and lons from the coastaldata and data and floor to ints
      clons=numpy.floor(coastaldata['longitude'])
      clats=numpy.floor(coastaldata['latitude'])
      dlons=numpy.floor(data['longitude'])
      dlats=numpy.floor(data['latitude'])
      data_to_add=[]
      #For each cruise in the data array
      for cr in list(set(data['expocode'])):
         #get the indices in the coastaldata for this cruise
         cindices=numpy.where(coastaldata['expocode']==cr)
         #Get the lon,lat pairs for this cruise
         cllpairs=list(set(zip(clons[cindices],clats[cindices])))
         #get the indices in the data array for this cruise
         dindices=numpy.where(data['expocode']==cr)
         #get the lon,lat pairs for the data array for this cruise
         dllpairs=list(set(zip(dlons[dindices],dlats[dindices])))
         #for each pair in the coastaldata if that pair is in the data array also
         #get the indices of those data and append to a list
         for pair in cllpairs:
            if pair in dllpairs:
               dllpairs.remove(pair)
               indices_to_move=numpy.where((numpy.floor(coastaldata[cindices]['longitude'])==pair[0])&
                                   (numpy.floor(coastaldata[cindices]['latitude'])==pair[1]))
               if len(indices_to_move)!=0:
                  data_to_add.append(cindices[0][indices_to_move[0]])
      #If data have been found in the coastaldata then add these onto the data array
      if len(data_to_add)!=0:
         all_indices_to_add=numpy.hstack(data_to_add)
         data=numpy.lib.recfunctions.stack_arrays((data,coastaldata[all_indices_to_add]),asrecarray=True,usemask=False)
         #now remove these data from the coastaldata array
         coastaldata=numpy.delete(coastaldata,all_indices_to_add)
   
   #Data has now had relevant coastal points appended to it.
   #Run the following for each pair of years in the range to process
   total_number_of_data_points=0
   total_duplicates=[]
   for year_range in year_ranges_to_process:
      number_of_data_points,duplicates=ConvertYears(data,year_range,sstdir,ssttail,prefix,outputdir,
                                                extrapolatetoyear,version=socatversion,ASCIIOUT=ASCIIOUT,
                                                percruisedir=percruisedir,removeduplicates=removeduplicates,
                                                useaatsr=useaatsr,usereynolds=usereynolds)
      total_number_of_data_points+=number_of_data_points
      total_duplicates.extend(duplicates)
   if len(total_duplicates)>0:
      print("Duplicates were found in %s and removed:"%inputfile)
      for dupe in total_duplicates:
         print(dupe)

   print("Total number of data points imported and used from file: %s : %d"%(inputfile,total_number_of_data_points))


def FinalCoastalConversion(startyr=cldefaults['start'],endyr=cldefaults['end'],prefix=cldefaults['prefix'],
                  outputdir=cldefaults['outputdir'],notperyear=cldefaults['notperyear'],sstdir=cldefaults['sstdir'],
                  ssttail=cldefaults['ssttail'],extrapolatetoyear=cldefaults['extrapolatetoyear'],
                  version=cldefaults['socatversion'],ASCIIOUT=cldefaults['asciioutput'],
                  percruisedir=cldefaults['asciioutput'],useaatsr=False,usereynolds=False,removeduplicates=True):
   #Run the following for each pair of years in the range to process
   print("Working on remaining coastal data.")
   total_number_of_data_points=0
   total_duplicates=[]
   year_ranges_to_process=GetYearsToProcess(startyr,endyr,notperyear)
   for year_range in year_ranges_to_process:
      number_of_data_points,duplicates=ConvertYears(coastaldata,year_range,sstdir,ssttail,prefix,outputdir,
                                                extrapolatetoyear,version=version,ASCIIOUT=ASCIIOUT,
                                                percruisedir=percruisedir,removeduplicates=removeduplicates,
                                                useaatsr=useaatsr,usereynolds=usereynolds)
      total_number_of_data_points+=number_of_data_points
      total_duplicates.extend(duplicates)
   if len(total_duplicates)>0:
      print("Duplicates were found in remaining coastal data and removed:")
      for dupe in total_duplicates:
         print(dupe)

   print("Total number of data points imported and used from remaining coastal data: %d"%(total_number_of_data_points))

def GetYearsToProcess(start,end,notperyear):
   """
   Get the range of years to process.
      start - the year at the start of the range (inclusive)
      end - the year at the end of the range (inclusive)
      notperyear - True if want the range to be converted together.
                 - False if want each year within this range to be converted individually
   """
   #Get the range of years to process so that we can process all years in one run of script
   #So, we create a list of [start,end] years dependent on notperyear: 
   # notperyear = True - then the start,end year range are converted together
   # notperyear = False - then for each year in the range (start,end) there is one run
   if notperyear is False:
      year_ranges=[[Y,Y] for Y in range(start,end+1)]
   else:
      year_ranges=[[start,end]]

   return year_ranges

def convert_column_id_to_index(header, columnIdentifier):
    if columnIdentifier == None:
        return None;
    else:
        try:
            if int(columnIdentifier) < len(header):
                index = int(columnIdentifier);
            else:
                raise IndexError("Index (%d) is out of range for header length of %d." % (int(columnIdentifier), len(header)));
        except ValueError:
            try:
                index = header.index(columnIdentifier);
            except ValueError:
                raise ValueError("column index '%s' could not be determined in header %s" % (columnIdentifier, str(header)));
    return index;

#columnInfo countains tuples for each column (standardName, dtype, identifier), where identifier is a string column name or int index corresponding
#   to the input files's header
#   set socatversion to None if using insitu
def ReadInData(inputfile, columnInfo, socatversion, delimiter='\t'):
    """
    Reads in the data from the SOCAT v2 or v3 files
       inputfile - the input SOCAT ascii csv file
       delimiter - the delimiter in the file (default to tabs)
    """
    with open(inputfile) as FILE:
       linestoskip=0;
       #Find how many lines we want to skip in the input SOCAT file by searching for 2 column headings
       if socatversion != None:
           for preline in FILE:
               if 'Expocode' not in preline or 'yr' not in preline:
                   linestoskip+=1
               else:
                   break;
           #Extract the header
           header = [s.strip() for s in preline.strip().split(delimiter)];
       else: #Using insitu, assume header is first line.
           header = [s.strip() for s in FILE.readline().split(delimiter)];
    

    columnInfoToExtract = [info for info in columnInfo if info[2] != None]; #These will be extracted from the datafile
    
    #Convert columns into indices
    indicesToExtract = [convert_column_id_to_index(header, info[2]) for info in columnInfoToExtract];
    namesOfExtracted = [info[0] for info in columnInfoToExtract];
    order = numpy.argsort(indicesToExtract); #pandas ignores the column order so we need to rearrange the column names accordingly
    namesOfExtracted = [namesOfExtracted[i] for i in order];
    
    print(namesOfExtracted);
    
    #dtypesOfExtracted = [info[1] for info in columnInfoToExtract];

    #Read in the columns we want into a data array
    print("Reading in data from SOCAT file: %s"%inputfile)
    print("This can take a while for large datasets.\n");

    data = pd.read_table(inputfile, skiprows=linestoskip+1, sep=delimiter, engine='c', usecols=indicesToExtract, names=namesOfExtracted, low_memory=False);#, dtype=dtypesOfExtracted);
    
    #Now insert additional columns filled with nan
    colNamesToInsert = [info[0] for info in columnInfo if info[2] == None];
    toInsert = numpy.full((len(data), len(colNamesToInsert)), numpy.nan);
    toInsert = pd.DataFrame(toInsert, columns=colNamesToInsert);
    data = data.join(toInsert);
    
    #Reorder columns
    orderedColNames = [info[0] for info in columnInfo];
    data = data[orderedColNames];

    #attach metadata to the data
    data = data.to_records(index=False);
    
    
    
    
    if "fCO2_qc_flag" in colNamesToInsert: #nan values must be float type, fCO2_qc_flag is usually int, so have to split here.
        data = data.astype([('expocode', 'S24'), ('year','<i8'), ('month','<i8'), ('day','<i8'), ('hour','<i8'), ('minute','<i8'), ('second','<i8'),
                     ('longitude','<f8'), ('latitude','<f8'), ('salinity', '<f8'),
                    ('SST', '<f8'), ('T_equ', '<f8'), ('air_pressure', '<f8'), ('air_pressure_equ', '<f8'),
                    ('salinity_sub', '<f8'), ('air_pressure_sub', '<f8'), ('fCO2', '<f8'), ('fCO2_qc_flag', '<f8')]);
    else: #fCO2_qc_flag was provided by the data file, so set type as int.
        data = data.astype([('expocode', 'S24'), ('year','<i8'), ('month','<i8'), ('day','<i8'), ('hour','<i8'), ('minute','<i8'), ('second','<i8'),
                         ('longitude','<f8'), ('latitude','<f8'), ('salinity', '<f8'),
                        ('SST', '<f8'), ('T_equ', '<f8'), ('air_pressure', '<f8'), ('air_pressure_equ', '<f8'),
                        ('salinity_sub', '<f8'), ('air_pressure_sub', '<f8'), ('fCO2', '<f8'), ('fCO2_qc_flag', '<i8')]);
    
    return data;

#def ReadInData(inputfile, columns, delimiter='\t'):
#   """
#   Reads in the data from the SOCAT v2 or v3 files
#      inputfile - the input SOCAT ascii csv file
#      delimiter - the delimiter in the file (default to tabs)
#   """
#   #Find how many lines we want to skip in the input SOCAT file by searching for 2 column headings
#   with open(inputfile) as SOCAT:
#      linestoskip=0
#      for preline in SOCAT:
#         if 'Expocode' not in preline or 'yr' not in preline:
#            linestoskip+=1
#         else:
#            break
#
#   #Read in the columns we want into a data array
#   print "Reading in data from SOCAT file: %s"%inputfile
#   print "This can take a while for large datasets.\n";
#   data=iter_loadtxt(filename=inputfile, columns=columns, delimiter=delimiter, skiprows=linestoskip)
#
#   return data

#def iter_loadtxt(filename, columns, delimiter='\t', skiprows=0):
#   """
#   Function to read in the data from the csv into numpy array
#
#   """
#   #This function is to replace the numpy.genfromtxt since that uses too much RAM
#   #for the machines we need to run on. This method uses less but is more complex
#
##   #Which columns do we want from the file - note 0 is the first column index
##   #Use these for SOCAT v2
##   SOCATv2_cols=[1,2,3,4,5,6,7,8,10,11,12,13,14,15,16,20,22]
##   #Use these for SOCAT v3
##   SOCATv3_cols=[4,5,6,7,8,9,10,11,13,14,15,16,17,18,19,23,25]
##      
##   if version==2:
##      usecols=SOCATv2_cols
##   elif version==3:
##      usecols=SOCATv3_cols
##   elif version==4: #SOCATv4, 5 and 6 use the same as v3.
##      usecols=SOCATv3_cols
##   elif version==5:
##      usecols=SOCATv3_cols
##   elif version==6:
##      usecols=SOCATv3_cols
##   else:
##      raise Exception("Specified version %d cannot be loaded. Currently supports versions 2, 3, 4, 5 and 6."%version)
#
#   #Function to use in numpy.fromiter()
#   def iter_func(cols,dtype,filename,skiprows,delimiter):
#      with open(filename) as SOCAT:
#         for _ in range(skiprows):#+1 to skip the column names also
#            next(SOCAT)
#         
#         #extract header and interpret columns as either indices or colnames.
#         iter_loadtxt.header = next(SOCAT).split("\t");
#         iter_loadtxt.colIndices = [];
#         for col in cols:
#             try:
#                 iter_loadtxt.colIndices.append(int(col));
#             except ValueError:
#                 try:
#                     iter_loadtxt.colIndices.append(iter_loadtxt.header.index(col));
#                 except ValueError:
#                     raise ValueError("column index '%s' could not be determined in file %s" % (col, filename));
#         
#         #start reading data from here
#         for line in SOCAT:
#            newline=[]
#            #strip off whitespace and split the line into delimited columns
#            line=line.rstrip().split(delimiter)
#            #for each column in the usecols list, extract that value into the newline
#            for colIndex in iter_loadtxt.colIndices:
#               newline.append(line[colIndex])
#            #for each item in the newline convert to a float
#            for item in newline:
#               try:
#                  yield dtype(item)
#               except Exception, e:
#                  raise Exception("Data value found in column of SOCAT data cannot be converted to float: %s"%(item)+str(e))
#      iter_loadtxt.rowlength=len(newline)
#   
#   #Create a numpy array using the above function
#   data=numpy.fromiter(iter_func(cols=columns,dtype=float,filename=filename,skiprows=skiprows,delimiter=delimiter),dtype=float)
#   #Reshape it
#   data=data.reshape((-1,iter_loadtxt.rowlength))
#   #Convert to a record array to be consistent with the genfromtxt result
#   data=data.view(numpy.recarray)
#   #Also get the expocode in case user wants to write to ascii list rather than grid
#   expocode=numpy.fromiter(iter_func([0],dtype=str,filename=filename,skiprows=skiprows,delimiter=delimiter),dtype='S24') #TMH: was S16 but some Expocodes are now longer than this
#   #now we need to add the field names - read from file
#   with open(filename) as SOCAT:
#      for _ in range(skiprows):
#         next(SOCAT)
#      #Now read in the column names removing the first one 
#      #as we don't keep that column in the data array below
#      names=next(SOCAT).rstrip().split(delimiter)
#
#   #Tidy up the names to match those generated from genfromtxt
#   #remove spaces,commas,brackets etc
#   names=[x.replace(" ","_").replace('[','').replace(']','').replace('.','').replace('/','') for x in names]
#   new_names=[]
#   for column in iter_loadtxt.colIndices:
#      new_names.append(names[column])
#   names=new_names
#
#   #Create a new dtype for the recarray - use all as float here for simplicity
#   dtype_list=[(x,float) for x in names]
#   dtype_list = numpy.dtype(dtype_list); #TMH: newer versions of numpy do not automatically do this conversion.
#   data.dtype=dtype_list
#
#   #Had to hard code here - this is not ideal. Could instead update the
#   #rest of the code to expect floats for all columns of data.
#   #Convert the datatypes to what genfromtxt would give 
#   data=data.astype([('yr','<i8'), ('mon','<i8'), ('day','<i8'), ('hh','<i8'), ('mm','<i8'), ('ss','<i8'),
#                     ('longitude_decdegE','<f8'), ('latitude_decdegN','<f8'), ('sal', '<f8'),
#                    ('SST_degC', '<f8'), ('Tequ_degC', '<f8'), ('PPPP_hPa', '<f8'), ('Pequ_hPa', '<f8'),
#                    ('WOA_SSS', '<f8'), ('NCEP_SLP_hPa', '<f8'), ('fCO2rec_uatm', '<f8'), ('fCO2rec_flag', '<i8')])
#   data=numpy.lib.recfunctions.append_fields(data, 'expocode', expocode,
#                                             dtypes=expocode.dtype, usemask=False, asrecarray=True)
#
#   return data

def ConvertYears(data,year_range,sstdir,ssttail,prefix,outputdir,extrapolatetoyear,version,
                 percruisedir=None,ASCIIOUT=False,removeduplicates=True,useaatsr=False,usereynolds=False):
   """
   Convert the data from the year range into netcdf files
      data - structured numpy array containing the SOCAT data 
      year_range - the range list(start,end) defining the range of years to convert
      sstdir - directory that contains the SST monthly climatology files
      ssttail - remainder of the SST netcdf name after the year and month
      prefix - prefix to use for the output file
      outputdir - directory to write output files to
      percruisedir - directory to write per-cruise files to
      ASCIIOUT - write out as ascii data
      withcoastal - True if using data from coastal file also
   """
   
   data_subset=[]
   #subset the year(s) we want
   print("Subsetting data for year range: %d %d"%(year_range[0],year_range[1]))
   #print len(data[numpy.where((data['year'] >= year_range[0]) & (data['year'] <= year_range[1]))]);
   data_subset=data[numpy.where((data['year'] >= year_range[0]) & (data['year'] <= year_range[1]))]
   

   #Test if there are any data - if not then return
   if data_subset.size==0:
      print('No data available for these years: %d %d'%(year_range[0],year_range[1]))
      return 0,[]

   #remove rows that fail quality checks
   #fco2 quality flag (if a SOCAT fCO2 quality control flag is present, use it to remove data that failed the check)
   if numpy.all(numpy.isnan(data_subset['fCO2_qc_flag'])) == False:
       data_subset=data_subset[numpy.where(data_subset['fCO2_qc_flag'] == 2)]
   
   #check if fco2 is not nan
   data_subset=data_subset[numpy.where(numpy.isfinite(data_subset['fCO2']))]

   #Now do some data conversion
   #we want lon in range [-180,180) so -360 from longitudes greater or equal to 180
   data_subset['longitude']=numpy.where(data_subset['longitude']>=180,
                                                data_subset['longitude']-360,data_subset['longitude'])
   #Make an array of all the dates julian days
   jds=datenum.datenum_array(data_subset['year'], data_subset['month'], data_subset['day'], data_subset['hour'], data_subset['minute'], 0)

   #Sort the data to aid in removing duplicates
   sortedindices=numpy.argsort(data_subset,order=['expocode','year','month','day','hour','minute','second','longitude','latitude'])
   data_subset=data_subset[sortedindices]
   #Also sort the jds array too
   jds=jds[sortedindices]
   duplicate_data=[]

   if removeduplicates:
      duplicates=[]
      #need to loop through and remove items which are the same in expocode, time and location
      for i in range(1,data_subset.size):
         duplicate=True
         #because the data are ordered we only need to check if consequetive records are the same
         #we will check if every column is identical - ignoring nans as they are not identical (undefined)
         for name in list(data_subset.dtype.names):
            #if data_subset[i][name]!=data_subset[i-1][name] and not numpy.isnan(data_subset[i][name]):
            if data_subset[i][name].dtype == float or data_subset[i][name].dtype == int: #Only float values can be nan (and things that can be cast to float)
               if data_subset[i][name]!=data_subset[i-1][name] and not numpy.isnan(data_subset[i][name]):
                    duplicate=False;
                    break;
            else: #not float
               if data_subset[i][name]!=data_subset[i-1][name]:
                   duplicate=False;
                   break;
                  
               #there is a difference so these cannot be duplicates - exit the loop
               duplicate=False
               break
         if duplicate:
            #duplicate has been found - append to a list for removal at end of loop
            duplicates.append(i)
            duplicate_data.append(data_subset[i])
      #now remove the duplicates
      l = len(duplicates);
      print("num duplicates is: ", l);
      if len(duplicates)>0:
         #get an array of all indices and then remove the ones marked as duplicates
         allindices=list(range(data_subset.size))
         for dupe in duplicates:
            allindices.remove(dupe)
         data_subset=data_subset[allindices]
         jds=jds[allindices]
   
   #Finally we can remove columns which were only required for quality checks
   #these are days,hours,mins,fCO2rec_flag
   names=list(data_subset.dtype.names)
   [names.remove(x) for x in ['fCO2_qc_flag']] #'day','mm','hh','ss',
   names = [str(x) for x in names]; #Unicode strings can't be used to index, but are returned as column names. Weird.
   data_subset=data_subset[names]

   #Optionally update names to be more sensible - could do this at data import instead
   #Note this ASSUMES the order of the column naming - FIXME: should do this more robustly #TMH: Fixed: Order is now determined by the columnInfo object after parsing command line args.
   #newnames=['yr','mon','day','hh','mm','ss','lon','lat','sal','SST_C','Teq_C','P','Peq','sal_woa','P_ncep','fCO2_rec','expocode']
   newnames=["expocode", "year", "month", "day", "hour", "minute", "second", "longitude", "latitude", "salinity", "sst", "T_equ", "air_pressure", "air_pressure_equ", "salinity_sub", "air_pressure_sub", "fCO2"];
   data_subset.dtype.names=newnames
   #keep track of the number of data points that are used (for output info only)
   number_of_data_points=data_subset.shape
   #Now do some actual conversion of the data
   #Get temperature from SST climatology
   if useaatsr and not usereynolds:
      #Use the AATSR data to get the SST
      Tcls = get_sst.GetAATSRSST(data_subset['year'], data_subset['month'], data_subset['longitude'], 
                              data_subset['latitude'],sstdir, ssttail)
      if numpy.all(Tcls==-999):
         print("All Temperature data are no-data-values - skipping for this year.")
         return 0,[]
      #Convert the temperature from skin to subskin
      Tcls += 0.17
   elif usereynolds and not useaatsr:
      #Use the Reynolds data to get the SST
      Tcls = get_sst.GetReynoldsSST(data_subset['year'], data_subset['month'], data_subset['longitude'], 
                              data_subset['latitude'],sstdir, ssttail)
      if numpy.all(Tcls==-999):
         print("All Temperature data are no-data-values - skipping for this year/month combination.")
         return 0,[]
      #Temperature is already (kind-of) subskin so no need to convert it
   else:
      raise Exception("No SST data specified. Currently must be one (and only one) of either AATSR or Reynolds.")

   #Update the pressure?
   Peq_cls = data_subset['air_pressure_sub'] + 3.
   #Extract the expocodes here as they get removed in the conversion
   expocodes=data_subset['expocode']
   #Recalculate the fugacity and partial pressure
   conversion = v2_f_conversion.v2_f_conversion_wrap(jds,data_subset,Tcls,Peq_cls,extrapolatetoyear)
   if conversion is None:
      #There were no good data to use
      return 0,[]
   
   #Write out the data into gridded monthly netCDF files - this will be easier if we append all arrays and use numpy
   #First convert jd_y into month - write a lambda function to do this so we can use numpy arrays
   ##convert_to_month=numpy.vectorize(lambda x: datetime.datetime.fromordinal(x).month)
   #apply function to jd_y (as integer)
   ##months=convert_to_month((conversion['jd']).astype(numpy.int))
   #months = conversion['mon'];
   #now concatenate arrays into one single array for easier/efficient splicing
   #conversion=numpy.lib.recfunctions.append_fields(conversion, 'month', months,
   #                                                      dtypes=months.dtype, usemask=False, asrecarray=True)

   #Update the file prefix to contain year info
   if extrapolatetoyear is not None:
      prefix=prefix+'_'+str(extrapolatetoyear)

   #For every month in turn
   for m in range(1,13): #i.e. for m 1-12 inclusive
      #create an output file for this months data
      if version != None:
          outputfile=prefix+'_from_%s_to_%s_%02d_v%d.nc'%(year_range[0],year_range[1],m,version)
      else:
          outputfile=prefix+'_from_%s_to_%s_%02d_custom_insitu.nc'%(year_range[0],year_range[1],m)
      outputfilepath=os.path.join(outputdir,"%02d"%m,outputfile)
      #get the data for this month
      month_data=conversion[numpy.where(conversion['mon']==m)]

      #Also extract the expocodes for the same data points and append on month_data
      expocodes_month=expocodes[numpy.where(conversion['mon']==m)]
#      month_data=numpy.lib.recfunctions.append_fields(month_data, 'expocode', expocodes_month,
#                                                   dtypes=expocodes_month.dtype, usemask=False, asrecarray=True)
      month_data=numpy.lib.recfunctions.append_fields(month_data, 'expocode', expocodes_month, dtypes=expocodes_month.dtype, usemask=False, asrecarray=True)
      
      
      #Get this month into a datetime object - use the average of year to get a centre point
      #Note in most usual cases the year range is a single year so averaging does nothing strange
      if m!=12:
         half_days_in_month=(datetime.datetime((year_range[0]+year_range[1])//2,m+1,1)- #TMH: / to // on switch to Python3
                          datetime.datetime((year_range[0]+year_range[1])//2,m,1)).days/2.0
         fraction=half_days_in_month-int(half_days_in_month)
         hours=int(fraction*24)
         minutes=0 # we wont bother with minutes
         seconds=0 #we wont bother with seconds
      else:
         #December needs to use the next year (as m+1==13) - so simpler to hard code as 15.5 days
         half_days_in_month=15.5
      datadate=datetime.datetime((year_range[0]+year_range[1])//2,m,int(half_days_in_month),hours,minutes,seconds) #TMH: / to // on switch to Python3
      if ASCIIOUT is False:
         #Now we want to average by cruise first and then bin all the cruises
         #this in effect takes a cruise-weighted average of the data
         if percruisedir is not None:
            percruiseoutput=os.path.join(outputdir,"%02d"%m,"per_cruise")
            #Loop through each expocode and average the data for that cruise
            for expo in set(month_data['expocode']):
               output_cruise_file=outputfile.replace('.nc','-%s.nc'%expo)
               output_cruise_file=os.path.join(percruiseoutput,output_cruise_file)
               expo_indices=numpy.where(month_data['expocode']==expo)
               expo_data=month_data[expo_indices]
               #If we write out as a grid we call this function
               #multiprocessing is required due to a suspected bug in netcdf library
               #which makes writing a netcdf fail after X variables have been written (in total)
               #running in a separate thread seems to reset the netcdf library for each run
               variabledictionary=CreateBinnedData(expo_data)
               myprocess=multiprocessing.Process(target=WriteOutToNCAsGrid,args=(variabledictionary,output_cruise_file,extrapolatetoyear,datadate))
               myprocess.start()
               myprocess.join()
            #Now we bin all the per-cruise data files
            #and write out as a netcdf file
            #(if there are any)
            common_prefix=outputfile.replace('.nc','')
            cruise_files=glob.glob("%s/%s*"%(percruiseoutput,common_prefix))
            if len(cruise_files) !=0:
               print()
               print("Combining cruises from region, year and month: ",cruise_files)
               print()
               combine_nc_files.FromFilelist(filelist=cruise_files,output=outputfilepath,
                                             weighting="cruise-weighted",outputtime=datadate)

               #get the binned data for the whole month to add to the nc file as other variables
               allnewvars=CreateBinnedData(month_data)
               newvars={v : allnewvars[v] for v in statvariables+['stds']}
               combine_nc_files.AddNewVariables(filename=outputfilepath,newvars=newvars)
         elif percruisedir is None:
            #give every observation an equal weighting - i.e. bin all data together
            variabledictionary=CreateBinnedData(month_data)
            WriteOutToNCAsGrid(variabledictionary,outputfilepath,extrapolatetoyear,outputtime=datadate)
      else:
         #If we write out as a list we call this function
         WriteOutToAsciiList(month_data,outputfilepath,extrapolatetoyear)
   return number_of_data_points[0],duplicate_data

def WriteOutToAsciiList(month_data,outputfile,extrapolatetoyear):

    #Calculate differences of data for this month
    #Difference in temperature
    dT=month_data['Tcl_C'] - month_data['SST_C']
    #Difference in fugacity
    dF=month_data['fCO2_Tym'] - month_data['fCO2_SST']
    #Difference in partial pressure
    dP=month_data['pCO2_Tym'] - month_data['pCO2_SST']
    
    outputfile=outputfile.replace('.nc','.txt')
    #Write out the data into a netCDF file
    #Test directory exists
    if not os.path.exists(os.path.dirname(outputfile)):
        raise Exception("Directory to write file to does not exist: %s"%(os.path.dirname(outputfile)))
    
    output_data=month_data

    if output_data.size > 0:
        print("Writing to: %s"%outputfile)
        numpy.savetxt(outputfile,output_data,fmt="%.7f,%d,%d,%d,%d,%d,%d,%.6f,%.6f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%d,%s",
                      header=",".join(output_data.dtype.names),delimiter=',')

def CreateBinnedData(month_data):
   #import pandas as pd;
   #allData = pd.DataFrame(month_data);
   
   #grid information
   nlon = 360 # number of longitude pixels
   lon0 = -180. # start longitude
   lon1 = 180. # end longitude
   nlat = 180 # number of latitude pixels
   lat0 = 90. # start latitude
   lat1 = -90. # end latitude
   dlon = (lon1 - lon0) / nlon
   dlat = (lat1 - lat0) / nlat

   #Set up arrays with default values
   dTs = numpy.zeros((nlat,nlon))+netcdf_helper.MISSINGDATAVALUE
   fCO2_Tyms = numpy.zeros((nlat,nlon))+netcdf_helper.MISSINGDATAVALUE
   fCO2_SSTs = numpy.zeros((nlat,nlon))+netcdf_helper.MISSINGDATAVALUE
   dFs = numpy.zeros((nlat,nlon))+netcdf_helper.MISSINGDATAVALUE
   pCO2_Tyms = numpy.zeros((nlat,nlon))+netcdf_helper.MISSINGDATAVALUE
   pCO2_SSTs = numpy.zeros((nlat,nlon))+netcdf_helper.MISSINGDATAVALUE
   dPs = numpy.zeros((nlat,nlon))+netcdf_helper.MISSINGDATAVALUE
   ndata = numpy.zeros((nlat,nlon))# keep track of multiple entries

   maximums={}
   minimums={}
   stds={}
   #Other data we want to keep track of
   for var in statvariables:
      maximums[var]=numpy.zeros((nlat,nlon))+netcdf_helper.MISSINGDATAVALUE
      minimums[var]=numpy.zeros((nlat,nlon))+netcdf_helper.MISSINGDATAVALUE
      stds[var]=numpy.zeros((nlat,nlon))+netcdf_helper.MISSINGDATAVALUE

   #Calculate differences of data for this month
   #Difference in temperature
   dT=month_data['Tcl_C'] - month_data['SST_C']
   #Difference in fugacity
   dF=month_data['fCO2_Tym'] - month_data['fCO2_SST']
   #Difference in partial pressure
   dP=month_data['pCO2_Tym'] - month_data['pCO2_SST']

   #Have changed method here as with previous one if lat data was an integer
   #then it fell into the wrong grid cell.
   #i.e. grid cells were 0<x<=1, 1<x<=2, ...
   #and we wanted 0<=x<1, 1<=x<2, ...
   #so round down to integer first, then subtract 1 [for lats only]
   ilons=numpy.int_((numpy.floor(month_data['lon'])-lon0) / dlon)
   ilats=numpy.int_((numpy.floor(month_data['lat'])-lat0) / dlat)-1

   #indicies where the ilons,ilats are within the grid bounds
   w=numpy.where((ilons >= 0)&(ilons < nlon)&(ilats >= 0)&(ilats < nlat))
   w=w[0]
   if w.size==0: return
   #update ilons,ilats to ones which fall in grid (all of them?)
   if len(w)>1:
      ilons, ilats = ilons[w], ilats[w]
   else:
      #edge case where w only has 1 element - in which case ilons/ilats are scalars not arrays
      #to fix it we convert them to a list of 1 element
      ilons=[ilons]
      ilats=[ilats]
   #update these data points to 0 (from -999) for all output arrays
   dTs[ilats,ilons]=0
   fCO2_Tyms[ilats, ilons] = 0.
   fCO2_SSTs[ilats,ilons]=0
   pCO2_Tyms[ilats, ilons] = 0.
   pCO2_SSTs[ilats,ilons]=0
   for var in statvariables:
      maximums[var][ilats,ilons]=0
      minimums[var][ilats,ilons]=0
      stds[var][ilats,ilons]=0

   #get a list of all the lat,lon pairs we have
   indices=list(set(zip(ilats,ilons)))
   for index in indices:
      #for each index in this list find all the points that fall into that grid cell index
      points=numpy.where(((ilats==index[0])&(ilons==index[1])))
      #Now we can bin these data into the grid cell
      fCO2_Tyms[index] += numpy.mean(month_data['fCO2_Tym'][points])
      pCO2_Tyms[index] += numpy.mean(month_data['pCO2_Tym'][points])
      fCO2_SSTs[index] += numpy.mean(month_data['fCO2_SST'][points])
      pCO2_SSTs[index] += numpy.mean(month_data['pCO2_SST'][points])
      dTs[index] += numpy.mean(dT[points])
      ndata[index] += points[0].size

      dFs[index]=fCO2_Tyms[index] - fCO2_SSTs[index]
      dPs[index]=pCO2_Tyms[index] - pCO2_SSTs[index]
      for var in statvariables:
         maximums[var][index]=numpy.nanmax(month_data[var][points])
         minimums[var][index]=numpy.nanmin(month_data[var][points])
         #we use ddof=1 to get the unbiased estimator (i.e. divide by N-1 in stdev formula)
         #unless only 1 element then set std = NAN (so that std are consistent)
         if len(points[0])>1:
            stds[var][index]=numpy.std(month_data[var][points],ddof=1)
         else:
            stds[var][index]=numpy.nan

   vardict={}
   vardict['fCO2_Tym']=fCO2_Tyms
   vardict['pCO2_Tym']=pCO2_Tyms
   vardict['fCO2_SST']=fCO2_SSTs
   vardict['pCO2_SST']=pCO2_SSTs
   vardict['dT']=dTs
   vardict['dF']=dFs
   vardict['dP']=dPs
   vardict['ndata']=ndata
   vardict['maximums']=maximums
   vardict['minimums']=minimums
   vardict['stds']=stds

   return vardict

def WriteOutToNCAsGrid(vardict,outputfile,extrapolatetoyear,outputtime=1e9):
   #grid information
   nlon = 360 # number of longitude pixels
   lon0 = -180. # start longitude
   lon1 = 180. # end longitude
   nlat = 180 # number of latitude pixels
   lat0 = 90. # start latitude
   lat1 = -90. # end latitude
   dlon = (lon1 - lon0) / nlon
   dlat = (lat1 - lat0) / nlat

   #Write out the data into a netCDF file
   print("Writing to: %s"%outputfile)
   #Test directory exists
   if not os.path.exists(os.path.dirname(outputfile)):
      raise Exception("Directory to write file to does not exist: %s"%(os.path.dirname(outputfile)))


   #update ndata in order to put as a variable in output file
   vardict['ndata']=numpy.where(vardict['ndata']==0,netcdf_helper.MISSINGDATAVALUE,vardict['ndata'])

   with netCDF4.Dataset(outputfile, 'w', format = 'NETCDF4') as ncfile:
      #Add the standard dims and variables
      netcdf_helper.standard_setup_SOCAT(ncfile,timedata=outputtime,londata=numpy.arange(lon0 + 0.5 * dlon, lon1, dlon),
                                          latdata=numpy.arange(lat0 + 0.5 * dlat, lat1, dlat))
      #set output names to allow extrapolation to any year or no extrapolation
      #so that the variable names make sense in the netcdf file
      if extrapolatetoyear is None:
         varext=""
         nameext=""
      else:
         varext="_%d"%extrapolatetoyear
         nameext=" extrapolated to %d"%extrapolatetoyear

      #Add the newly calculated fugacity, partial pressures and differences
      fCO2_SST_data = ncfile.createVariable('fCO2_SST','f4',('time','latitude','longitude'),fill_value=netcdf_helper.MISSINGDATAVALUE,zlib=True)
      fCO2_SST_data[:] = vardict['fCO2_SST']
      fCO2_SST_data.units = 'uatm'
      fCO2_SST_data.missing_value = netcdf_helper.MISSINGDATAVALUE
      fCO2_SST_data.valid_min = 0.
      fCO2_SST_data.valid_max = 1e6
      fCO2_SST_data.scale_factor = 1.
      fCO2_SST_data.add_offset = 0.
      fCO2_SST_data.standard_name = "fCO2_SST"
      fCO2_SST_data.long_name = "CO2 fugacity using SOCAT methodology"

      fCO2_Tym_data = ncfile.createVariable('fCO2_Tym'+varext,'f4',('time','latitude','longitude'),fill_value=netcdf_helper.MISSINGDATAVALUE,zlib=True)
      fCO2_Tym_data[:] = vardict['fCO2_Tym']
      fCO2_Tym_data.units = 'uatm'
      fCO2_Tym_data.missing_value = netcdf_helper.MISSINGDATAVALUE
      fCO2_Tym_data.valid_min = 0.
      fCO2_Tym_data.valid_max = 1e6
      fCO2_Tym_data.scale_factor = 1.
      fCO2_Tym_data.add_offset = 0.
      fCO2_Tym_data.standard_name = "fCO2_Tym"+varext
      fCO2_Tym_data.long_name = "CO2 fugacity using OC-FLUX methodology"+nameext

      pCO2_SST_data = ncfile.createVariable('pCO2_SST','f4',('time','latitude','longitude'),fill_value=netcdf_helper.MISSINGDATAVALUE,zlib=True)
      pCO2_SST_data[:] = vardict['pCO2_SST']
      pCO2_SST_data.units = 'uatm'
      pCO2_SST_data.missing_value = netcdf_helper.MISSINGDATAVALUE
      pCO2_SST_data.valid_min = 0.
      pCO2_SST_data.valid_max = 1e6
      pCO2_SST_data.scale_factor = 1.
      pCO2_SST_data.add_offset = 0.
      pCO2_SST_data.standard_name = "pCO2_SST"
      pCO2_SST_data.long_name = "CO2 partial pressure using SOCAT methodology"

      pCO2_Tym_data = ncfile.createVariable('pCO2_Tym'+varext,'f4',('time','latitude','longitude'),fill_value=netcdf_helper.MISSINGDATAVALUE,zlib=True)
      pCO2_Tym_data[:] = vardict['pCO2_Tym']
      pCO2_Tym_data.units = 'uatm'
      pCO2_Tym_data.missing_value = netcdf_helper.MISSINGDATAVALUE
      pCO2_Tym_data.valid_min = 0.
      pCO2_Tym_data.valid_max = 1e6
      pCO2_Tym_data.scale_factor = 1.
      pCO2_Tym_data.add_offset = 0.
      pCO2_Tym_data.standard_name = "pCO2_Tym"+varext
      pCO2_Tym_data.long_name = "CO2 partial pressure using OC-FLUX methodology"+nameext

      dT_data = ncfile.createVariable('dT','f4',('time','latitude','longitude'),fill_value=netcdf_helper.MISSINGDATAVALUE,zlib=True)
      dT_data[:] = vardict['dT']
      dT_data.units = 'Degree C'
      dT_data.missing_value = netcdf_helper.MISSINGDATAVALUE
      dT_data.valid_min = -999.
      dT_data.valid_max = 999.
      dT_data.scale_factor = 1.
      dT_data.add_offset = 0.
      dT_data.standard_name = "dT"
      dT_data.long_name = "difference Tym - SST"

      dfCO2_data = ncfile.createVariable('dfCO2','f4',('time','latitude','longitude'),fill_value=netcdf_helper.MISSINGDATAVALUE,zlib=True)
      dfCO2_data[:] = vardict['dF']
      dfCO2_data.units = 'uatm'
      dfCO2_data.missing_value = netcdf_helper.MISSINGDATAVALUE
      dfCO2_data.valid_min = -999.
      dfCO2_data.valid_max = 999.
      dfCO2_data.scale_factor = 1.
      dfCO2_data.add_offset = 0.
      dfCO2_data.standard_name = "dfCO2"
      dfCO2_data.long_name = "difference fCO2,Tym - fCO2,SST"

      dpCO2_data = ncfile.createVariable('dpCO2','f4',('time','latitude','longitude'),fill_value=netcdf_helper.MISSINGDATAVALUE,zlib=True)
      dpCO2_data[:] = vardict['dP']
      dpCO2_data.units = 'uatm'
      dpCO2_data.missing_value = netcdf_helper.MISSINGDATAVALUE
      dpCO2_data.valid_min = -999.
      dpCO2_data.valid_max = 999.
      dpCO2_data.scale_factor = 1.
      dpCO2_data.add_offset = 0.
      dpCO2_data.standard_name = "dpCO2"
      dpCO2_data.long_name = "difference pCO2,Tym - pCO2,SST"

      N_data = ncfile.createVariable('count_nobs','f4',('time','latitude','longitude'),fill_value=netcdf_helper.MISSINGDATAVALUE,zlib=True)
      N_data[:] = vardict['ndata']
      N_data.units = 'count'
      N_data.missing_value = netcdf_helper.MISSINGDATAVALUE
      N_data.valid_min = 0.
      N_data.valid_max = 10000000.
      N_data.scale_factor = 1.
      N_data.add_offset = 0.
      N_data.standard_name = "count_nobs"
      N_data.long_name = "Number of observations mean-averaged in cell"

      for var in statvariables:
         if var in ['fCO2_Tym','pCO2_Tym']:
            stats_varext=varext
         else:
            stats_varext=""
         minf_data = ncfile.createVariable('min_'+var+stats_varext,'f4',('time','latitude','longitude'),fill_value=netcdf_helper.MISSINGDATAVALUE,zlib=True)
         minf_data[:] = vardict['minimums'][var]
         minf_data.units = 'uatm'
         minf_data.missing_value = netcdf_helper.MISSINGDATAVALUE
         minf_data.valid_min = 0.
         minf_data.valid_max = 1000000.
         minf_data.scale_factor = 1.
         minf_data.add_offset = 0.
         minf_data.standard_name = "min_"+var+stats_varext
         minf_data.long_name = "Minimum "+var+stats_varext+" occupying binned cell"

         maxf_data = ncfile.createVariable('max_'+var+stats_varext,'f4',('time','latitude','longitude'),fill_value=netcdf_helper.MISSINGDATAVALUE,zlib=True)
         maxf_data[:] = vardict['maximums'][var]
         maxf_data.units = 'uatm'
         maxf_data.missing_value = netcdf_helper.MISSINGDATAVALUE
         maxf_data.valid_min = 0.
         maxf_data.valid_max = 1000000.
         maxf_data.scale_factor = 1.
         maxf_data.add_offset = 0.
         maxf_data.standard_name = "max_"+var+stats_varext
         maxf_data.long_name = "Maximum "+var+stats_varext+" occupying binned cell"

         stdf_data = ncfile.createVariable('std_'+var+stats_varext,'f4',('time','latitude','longitude'),fill_value=netcdf_helper.MISSINGDATAVALUE,zlib=True)
         stdf_data[:] = vardict['stds'][var]
         stdf_data.units = 'uatm'
         stdf_data.missing_value = netcdf_helper.MISSINGDATAVALUE
         stdf_data.valid_min = 0.
         stdf_data.valid_max = 1000000.
         stdf_data.scale_factor = 1.
         stdf_data.add_offset = 0.
         stdf_data.standard_name = "stdev_"+var+stats_varext
         stdf_data.long_name = "Standard deviation of "+var+stats_varext+" occupying binned cell"

#if the script is run stand alone (i.e. from command line)
def Main():
   print("%s started at: %s "%(os.path.basename(__file__),str(datetime.datetime.now())))
   cl=GetCommandline()
   DoConversion(inputfile=cl.inputfile,startyr=cl.start,endyr=cl.end,notperyear=cl.notperyear,
                sstdir=cl.sstdir,ssttail=cl.ssttail,prefix=cl.prefix,outputdir=cl.outputdir,
                extrapolatetoyear=cl.extrapolatetoyear,version=cl.socatversion,ASCIIOUT=cl.asciioutput,
                percruisedir=cl.percruisedir,coastalfile=cl.coastalfile,useaatsr=cl.useaatsr,
                usereynolds=cl.usereynolds,removeduplicates=not cl.keepduplicates)
   print("%s ended at: %s "%(os.path.basename(__file__),str(datetime.datetime.now())))

if __name__=="__main__":
   Main()

