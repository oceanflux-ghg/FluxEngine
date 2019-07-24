#!/home/cercache/tools/environments/run_python_scibox_precise.sh
#-*- coding:utf-8 -*-

#NOTE. At some point between Numpy v1.8.2 and v1.11.0 the numpy library interface for
#rec arrays changed and the code no longer works for version 1.11.0. It worked for version 1.8.2.
#   TMH: Fixed issue with rec arrays (2018-06-02).
#works well on Fedora 23

#common python libs
import subprocess
import argparse
import string
import os
import glob
import datetime
import multiprocessing

#internal libs
from . import v2_convert_f_SSH
from . import combine_xy
from . import combine_nc_files
from . import netcdf_helper

#Returns a dictionary mapping region codes (keys) to input file names
def GenerateRegionFileMap(socatfiles, regions):
    if len(socatfiles) == 1 and regions == None:
        regionFileMap = {"GL":socatfiles[0]};
        return regionFileMap;
    
    #Multiple regions defined
    elif len(socatfiles) == len(regions):
        regionFileMap = {regions[i]:socatfiles[i] for i in range(0, len(regions))};
        if len(list(regionFileMap.keys())) != len(regions): #Checks we're not overwriting files by using the same region key more than once.
            raise ValueError("More than one file per region is not currently supported. To get around this you can define additional regions or concatenate input files manually.");
        if len(regions) > 1:
            print("**********\nWarning using multiple regions: Region boundaries must be defined such that more than one region does not overlap a grid cell. If data from more than one region fall into one grid cell the values corresponding to this grid cell will be inaccurate. This constraint does not apply to regions which are specified as coastal using 'usecoastal'.");
            input("Press a key to continue...");
        return regionFileMap;
    #Not global and there are more files than regions or more regions than files.
    else:
        raise ValueError("GenerateRegionFileMap: number of region codes suppled do not match the number of input files.");
        

#Output directories for the various data: $dirname will be replaced with the output directory 
#at a later point in the processing
outputdirectory={'socat' : string.Template('$dirname/reanalysed_data'),
                 'socatmonth' : string.Template('$dirname/reanalysed_data/%02d'),
                 'socatmonthpercruise' : string.Template('$dirname/reanalysed_data/%02d/per_cruise'),
                 'joined' : string.Template('$dirname/joined'),
                 'combined' : string.Template('$dirname/combined'),
                 'fpred' : string.Template('$dirname/interpolated/fCO2/predictions'),
                 'ppred' : string.Template('$dirname/interpolated/pCO2/predictions'),
                 'fvar' : string.Template('$dirname/interpolated/fCO2/variances'),
                 'pvar' : string.Template('$dirname/interpolated/pCO2/variances'),
                 'temp' : string.Template('$dirname/temp_workspaces'),
                 'reanalysedglobal' : string.Template('$dirname/reanalysed_global'),
                 'reanalysedglobalmonth' : string.Template('$dirname/reanalysed_global/%02d'),
                }

#Use of global means we need to explicitly reset this between calls. TODO: Refactor code to that outputdirectory is not global.
def reset_outputdirectory():
    global outputdirectory;
    outputdirectory={'socat' : string.Template('$dirname/reanalysed_data'),
                 'socatmonth' : string.Template('$dirname/reanalysed_data/%02d'),
                 'socatmonthpercruise' : string.Template('$dirname/reanalysed_data/%02d/per_cruise'),
                 'joined' : string.Template('$dirname/joined'),
                 'combined' : string.Template('$dirname/combined'),
                 'fpred' : string.Template('$dirname/interpolated/fCO2/predictions'),
                 'ppred' : string.Template('$dirname/interpolated/pCO2/predictions'),
                 'fvar' : string.Template('$dirname/interpolated/fCO2/variances'),
                 'pvar' : string.Template('$dirname/interpolated/pCO2/variances'),
                 'temp' : string.Template('$dirname/temp_workspaces'),
                 'reanalysedglobal' : string.Template('$dirname/reanalysed_global'),
                 'reanalysedglobalmonth' : string.Template('$dirname/reanalysed_global/%02d'),
                };

#Template for the SOCAT string filename
#SOCATinfile='SOCATv%d_%s.tsv'
#output_global_monthly_file='OCF-CO2-GLO-M1-100-SOCAT-CONV.nc'

#Location of GSTAT executable - would be better in the path?
GSTATEXE='/home/cercache/project/oceanflux-shared/software/gstat-2.4.0/gstat'
#Location of the DIVA3D executable wrapper
DIVA3D="%s/run_diva.sh"%os.path.dirname(os.path.realpath(__file__))
#Location of the default diva parameter file
divaparamdefault="/home/cercache/project/oceanflux-shared/workspace/mwarren/default_parameters.param"

#These are the columns of the 'joined' ascii files that contain the fCO2 and pCO2 values
#that are to be interpolated (indexed from 1)
FCO2=3
PCO2=4

#This function handles command line parameters
def GetCommandline():
   """
   Function that handles the command line
   """
   parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
   arg_datalocationgroup=parser.add_argument_group(title='Data location')
   arg_datalocationgroup.add_argument('--sstdir', metavar='<path>',help ='Location of the directory containing the SST data files',
                        type=str,default=None)
   arg_datalocationgroup.add_argument('--ssttail',type=str, metavar='<path>',help ='The SST file prefix (after yyyymm).',
                        default=v2_convert_f_SSH.cldefaults['ssttail'])
   arg_datalocationgroup.add_argument('--vco2dir', metavar='<path>',help ='Location of the directory containing the vCO2 data files',
                        type=str,default=None)
   arg_datalocationgroup.add_argument('--socatdir', metavar='<path>',help ='Location of the directory containing the SOCAT data files',
                        type=str,default=None)
   arg_datalocationgroup.add_argument('--socatfiles', nargs='*', help='List of files containing data',
                        type=str,default=None);
   parser.add_argument('--output', metavar='<path>',help ='Location of the top level directory to write output - this will be created so should not already exist.',
                        type=str,default="./output")
   parser.add_argument('--startyr', metavar='<year>',type=int,help ='The year to start importing SOCAT data for.',default="1991")
   parser.add_argument('--endyr', metavar='<year>',type=int,help ='The year to end importing SOCAT data for.',default="2010")
   parser.add_argument('--regions', metavar='<keyword_list>',type=str,nargs='*',help ='A list of region prefixes to process.',
                        default=None)
   arg_interpolationgroup=parser.add_argument_group(title='Interpolation methods')
   #arg_interpolationgroup.add_argument('--none',action='store_true',help='Only process up to the interpolation stage and then exit.',default=False)
   arg_interpolationgroup.add_argument('--methodused',help='Assume the interpolation stage has been completed and carry on from that stage. Give the interpolation method used: DIVA or GSTAT',default=None)
   arg_interpolationgroup.add_argument('--diva',action='store_true',help='Use the DIVA interpolation routine',default=False)
   arg_interpolationgroup.add_argument('--gstatcmds', metavar='<path>',help ='The directory containing gstat command files for kriging with gstat.',default=None)
   arg_convertgroup=parser.add_argument_group(title='Subscript parameters')
   arg_convertgroup.add_argument('--notperyear',dest='notperyear',action='store_true',help="Use this if you don't want separate files per year",default=v2_convert_f_SSH.cldefaults['notperyear'])
   arg_convertgroup.add_argument('--extrapolatetoyear',type=int,dest='extrapolatetoyear',help="Use this flag to extrapolate to a given year (using the Takahashi trend).",default=None)
   arg_convertgroup.add_argument('--socatversion',type=int,dest='socatversion',help="The version of the SOCAT data files to read",default=2)
   arg_convertgroup.add_argument('--asciioutput',dest='asciioutput',action='store_true',help="To output data as ascii lists rather than gridded netcdf.",default=v2_convert_f_SSH.cldefaults['asciioutput'])
   arg_convertgroup.add_argument('--withcoastal',dest='withcoastal',type=str,help="The region code (defined using --regions) which corresponds to coastal data. Coastal data is appended to other regions where gridcells overlap and any remaining data is analysed seperately. Do not specify if no coastal data is used. Default is None (no coastal data used).",default=None)
   arg_convertgroup.add_argument('--useaatsr',dest='useaatsr',action='store_true',help="To use the AATSR SST data.",default=False)
   arg_convertgroup.add_argument('--usereynolds',dest='usereynolds',action='store_true',help="To use the Reynolds SST data.",default=False)
   arg_convertgroup.add_argument('--keepduplicates',dest='keepduplicates',action='store_true',help="To use duplicate data.",default=False)

   cl=parser.parse_args()

   return cl;

def CreateOutputTree(topdir,continuing=False):
   """
   Function that creates the output directory tree within the given top level directory.
   The top level directory is also created here, and so should not already exist.
      topdir - the name to give the output directory
      continuing - if continuing (i.e. a second run) after interpolation then True else False
   """
   #Do some path munging to get a full path name (deal with '~')
   topdir=os.path.abspath(os.path.expanduser(topdir))

   #For ease in later functions update the output directory variables to include topdir
   for key in list(outputdirectory.keys()):
      #if os.path.isabs(outputdirectory[key]) == False: #TMH: when using insitu data this isn't already appended.
          outputdirectory[key]=outputdirectory[key].safe_substitute(dirname=topdir)

   #If not continuing then we want to create a new output directory
   #else if continuing is True then we already have output directory tree existing
   if not continuing:
      #Check if topdir exists - if so then exit (to avoid accidental overwrites)
      if os.path.exists(topdir): 
         raise Exception("Output directory: %s already exists - please specify a new one that does not exist."%topdir)

      #Make the topdir directory
      os.makedirs(topdir)

      #Make all the other sub directories
      os.mkdir(outputdirectory['socat']) #SOCAT outputs
      os.mkdir(outputdirectory['joined']) #Joined files
      os.mkdir(outputdirectory['combined']) #Combined files
      os.makedirs(outputdirectory['fpred']) #interpolated fCO2 predictions
      os.makedirs(outputdirectory['ppred']) #interpolated pCO2 predictions
      os.makedirs(outputdirectory['fvar']) #interpolated fCO2 variances
      os.makedirs(outputdirectory['pvar']) #interpolated pCO2 variances
      os.mkdir(outputdirectory['temp']) #temporary workspaces
      os.mkdir(outputdirectory['reanalysedglobal'])#global monthly reanalysed product
      #Make month subdirs in SOCAT output directory
      for i in range(1,13):
         os.mkdir(outputdirectory['socatmonth']%(i))
         os.mkdir(outputdirectory['socatmonthpercruise']%(i))
         os.mkdir(outputdirectory['reanalysedglobalmonth']%(i))

def KrigeWithGstat(gstatcommdir):
   """
   Experimental function to allow kriging using the gstat application 
      gstatcommdir: directory containing all the gstat cmd files for each month
   """
   #Get all the files that match the 'glob' in the gstatcommdir
   gstatcommfiles=glob.glob('%s/*.cmd'%gstatcommdir)
   print("Found %d gstat command files."%len(gstatcommfiles))
   #For each command file - run it with gstat and log the output
   for filename in sorted(gstatcommfiles):
      command=[GSTATEXE,filename]
      #Run the command
      process=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
      #wait until process is finished and get the output messages
      stdout,stderr=process.communicate()
      #write output messages to a log file
      logfile="%s-log.txt"%filename
      fout=open(logfile,'w')
      fout.writelines(stdout+stderr)
      fout.close()

def InterpolateWithDiva(inputfilename,destinationpath,variable,divaparamdefault):
   """
   Function to call the third party GPLv3 DIVA interplation routines: 
      http://modb.oce.ulg.ac.be/mediawiki/index.php/DIVA
   The routines are written in Fortran and bash so will be called via subprocess.
   """
   #Build up the command that we need to run
   command=[]
   command.extend([DIVA3D])
   #add the command line arguments
   command.extend(['-d',destinationpath])
   command.extend(['-f',inputfilename])
   command.extend(['-c',str(variable)])
   command.extend(['-p',divaparamdefault])
   command.extend(['-t',outputdirectory['temp']])

   #use the python subprocess module to run the scripts
   process=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
   stdout,stderr=process.communicate()
   #write output messages to a log file
   logfile="%s-log.txt"%os.path.join(destinationpath,os.path.basename(inputfilename))
   fout=open(logfile,'w')
   fout.writelines(stdout+stderr)
   fout.close()


#Returns standard column names, dtypes and columns as string names or string column numbers   
def construct_column_info(year_col, month_col, day_col, hour_col, minute_col, second_col, longitude_col, latitude_col, \
                                   salinity_col, salinity_sub_col, SST_C_col, Tequ_col, air_pressure_col, air_pressure_sub_col, air_pressure_equ_col, \
                                   fCO2_col, expocode_col, socatversion, notsocatformat):
    stndColNames = ["expocode", "year", "month", "day", "hour", "minute", "second", "longitude", "latitude", "salinity", "sst", "T_equ", \
                    "air_pressure", "air_pressure_equ", "salinity_sub", "air_pressure_sub", "fCO2", "fCO2_qc_flag"];
    colDTypes = ['U24', '<i8', '<i8', '<i8', '<i8', '<i8', '<i8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<i8'];
    #colDTypes = [str, int, int, int, int, int, int, float, float, float, float, float, float, float, float, float, float, int];
    
    if notsocatformat: #Specify the columns as given by the command line parameters.
        colIdentifiers = [expocode_col, year_col, month_col, day_col, hour_col, minute_col, second_col, longitude_col, latitude_col, \
                                   salinity_col, SST_C_col, Tequ_col, air_pressure_col, air_pressure_equ_col, salinity_sub_col, air_pressure_sub_col, \
                                   fCO2_col, None];
    else: #using socat so use socatversion to determine correct columns
        if socatversion == 2:
            colIdentifiers = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '10', '11', '12', '13', '14', '15', '16', '20', '22'];
#            return ['SOCAT_DOI',
#             'QC_ID',
#             'yr',
#             'mon',
#             'day',
#             'hh',
#             'mm',
#             'ss',
#             'latitude [dec.deg.N]',
#             'sample_depth [m]',
#             'sal',
#             'SST [deg.C]',
#             'Tequ [deg.C]',
#             'PPPP [hPa]',
#             'Pequ [hPa]',
#             'd2l [km]',
#             'fCO2rec [uatm]'];
        elif socatversion in [3, 4, 5, 6]:
            colIdentifiers = ['0', '4', '5', '6', '7', '8', '9', '10', '11', '13', '14', '15', '16', '17', '18', '19', '23', '25'];
            #colIdentifiers = ['0', '4', '5', '6', '7', '8', '9', '10', '11', '13', '18', '14', '15', '16', '19', '17', '23', '25'];
#            return ['yr',
#             'mon',
#             'day',
#             'hh',
#             'mm',
#             'ss',
#             'longitude [dec.deg.E]',
#             'latitude [dec.deg.N]',
#             'sal',
#             'SST [deg.C]',
#             'Tequ [deg.C]',
#             'PPPP [hPa]',
#             'Pequ [hPa]',
#             'WOA_SSS',
#             'NCEP_SLP [hPa]',
#             'fCO2rec [uatm]',
#             'fCO2rec_flag']
        else:
            raise ValueError("No value columns could be generated. Only socat version 2, 3, 4, 5, and 6 are supported.");
        
    columnInfo = [];
    for i in range(len(colIdentifiers)):
        columnInfo.append( (stndColNames[i], colDTypes[i], colIdentifiers[i]) );
    return columnInfo;

#For argument descriptions, see GetCommandLine function in this file.
def RunReanalyseSocat(socatdir=None, socatfiles=None, sstdir=None, ssttail=None, vco2dir=None, output="./output", startyr=2010,
                      endyr=2010, regions=None,  methodused=None, diva=False, gstatcmds=None,
                      notperyear=v2_convert_f_SSH.cldefaults['notperyear'], extrapolatetoyear=None, keepduplicates=False,
                      asciioutput=v2_convert_f_SSH.cldefaults['asciioutput'], withcoastal=False, useaatsr=False, usereynolds=False,
                      socatversion=6,
                      notsocatformat=False,
                      year_col=None,
                      month_col=None,
                      day_col=None,
                      hour_col=None,
                      minute_col=None,
                      second_col=None,
                      longitude_col=None,
                      latitude_col=None,
                      salinity_col=None,
                      salinity_sub_col=None,
                      SST_C_col=None,
                      Tequ_col=None,
                      air_pressure_col=None,
                      air_pressure_sub_col=None,
                      air_pressure_equ_col=None,
                      fCO2_col=None,
                      expocode_col=None
                      ):
   #Dodgy use of global means we need to explicitly reset 'outputdirectory'this between calls.
   reset_outputdirectory();
   
   
   columnInfo = construct_column_info(year_col, month_col, day_col, hour_col, minute_col, second_col, longitude_col, latitude_col, \
                                   salinity_col, salinity_sub_col, SST_C_col, Tequ_col, air_pressure_col, air_pressure_sub_col, air_pressure_equ_col, \
                                   fCO2_col, expocode_col, socatversion, notsocatformat);
   
   
   #Dictionary mapping region codes with names.
   regionFileMap=GenerateRegionFileMap(socatfiles, regions);
   
   ###Check for valid parameters
   if useaatsr and usereynolds:
       raise Exception("Can only use one of AATSR or Reynolds SST products not both.")
   if not useaatsr and not usereynolds:
       raise Exception("Must specify one of AATSR or Reynolds SST products.")

   if usereynolds:
      print()
      print("ATTENTION: selected to use Reynolds SST. This assumes that the sstdir points to the Reynolds data (if not you need to change it).")
      print()
   if useaatsr:
      print()
      print("ATTENTION: selected to use AATSR SST. This assumes that the sstdir points to the AATSR data (if not you need to change it).")
      print()
   

   if methodused is not None and (diva==True or gstatcmds is not None):
      print("Can only select one of --methodused, --diva or --gstatcmds.")
      return 1;

   elif diva==True and gstatcmds is not None:
      print("Can only select one of --methodused, --diva or --gstatcmds.")
      return 1;

   if methodused is not None:
      allowedmethods=["DIVA","GSTAT"]
      if methodused not in allowedmethods:
         print("Interpolation method should be one of: ", allowedmethods)
         return 1;

   #Expand and convert to absolute file paths
   if sstdir != None: sstdir = os.path.abspath(os.path.expanduser(sstdir));
   if vco2dir != None: vco2dir = os.path.abspath(os.path.expanduser(vco2dir));
   if socatdir != None: socatdir = os.path.abspath(os.path.expanduser(socatdir));
   output = os.path.abspath(os.path.expanduser(output));
   
   #Check if paths exist
   for path in [sstdir, vco2dir, socatdir]:
       if path != "" and path != None:
           if not os.path.exists(os.path.expanduser(path)):
               print("File path does not exist: %s" % path)
               return 2;
   
   
   #######
   ###Start main logic
   #######

   #Create a new variable as part to specify whether we exit at the interpolation stage or not
   exitatinterpolation=False;
   if methodused is None and diva==False and gstatcmds is None:
      #No interpolation has been specified. Exit when reach interpolation stage
      exitatinterpolation=True;
      print("No interpolation option specified so will exit at interpolation stage.");
   
   #Get command line variables
   interpmethod=methodused
   #Create the output directory tree (if methodused is false)
   if methodused is None:
      cont=False
   else:
      cont=True
   CreateOutputTree(output,cont)

   removeduplicates = not keepduplicates

   #This bit is a little experimental ...
   #only do this if methodused has not been specified
   if methodused is None:
      if withcoastal != None:
         coastalfile=os.path.join(socatdir, regionFileMap[withcoastal]);
      else:
         coastalfile=None
      #Run the conversion of the SOCAT ascii files to netCDF
      for region in regions:
         #we now skip coastal regions as that data is read in at the same
         #time as the other data to make cruise weighting possible
         #unless --asciioutput is asked for.
         if region == "CO":
            if not asciioutput:
               continue

         v2_convert_f_SSH.DoConversion(os.path.join(socatdir, regionFileMap[region]),
                                                 columnInfo,
                                                 startyr,
                                                 endyr,
                                                 region,
                                                 outputdirectory['socat'],
                                                 notperyear,
                                                 sstdir,
                                                 ssttail,
                                                 extrapolatetoyear,
                                                 asciioutput,
                                                 socatversion,
                                                 outputdirectory['socatmonthpercruise'],
                                                 coastalfile,
                                                 useaatsr,
                                                 usereynolds,
                                                 removeduplicates)

      if withcoastal != None:
         v2_convert_f_SSH.FinalCoastalConversion(startyr,
                                                 endyr,
                                                 "CO",#region, #TMH: otherwise overwrites last region used.
                                                 outputdirectory['socat'],
                                                 notperyear,
                                                 sstdir,
                                                 ssttail,
                                                 extrapolatetoyear,
                                                 asciioutput,
                                                 outputdirectory['socatmonthpercruise'],
                                                 useaatsr,
                                                 usereynolds,
                                                 removeduplicates)

      #Get the variable names based on whether we have extrapolated or not
      #These are (only) needed when we write out to ASCII
      if extrapolatetoyear is not None:
         fvar='weighted_fCO2_Tym_%d'%extrapolatetoyear
         pvar='weighted_pCO2_Tym_%d'%extrapolatetoyear
      else:
         fvar='weighted_fCO2_Tym'
         pvar='weighted_pCO2_Tym'

      #Here we can combine the regional netCDFs into a global product
      #we will do it per month per year
      #Only do if netCDFs were produced in previous step 
      #TODO FIXME AND if more than 1 regions selected
      if not asciioutput:
         print("Creating global netCDF products ...")
         for month in range(1,13):
            for year in range(startyr,endyr+1):
               global_output_filename="%s/%d%02d01-OCF-CO2-GLO-1M-100-SOCAT-CONV.nc"%(
                                                      outputdirectory['reanalysedglobalmonth']%(month),year,month)
               files_for_global=glob.glob("%s/%02d/*%s_to_%s*.nc"%(outputdirectory['socat'],month,str(year),str(year)))
               if month!=12:
                  half_days_in_month=(datetime.datetime(year,month+1,1)-
                                 datetime.datetime(year,month,1)).days/2.0
                  fraction=half_days_in_month-int(half_days_in_month)
                  hours=int(fraction*24)
               else:
                  #December needs to use the next year (as m+1==13) - so simpler to hard code as 15.5 days
                  half_days_in_month=15.5
               datadate=datetime.datetime(year,month,int(half_days_in_month),hours,0,0)
               #Using multiprocess as netcdf4 lib throws HDF error if there are too many
               #files/variables written out in loops
               #NOTE Have changed the last parameter to False - this is because in SOCATv4 it didn't appear
               #to be the case that the regions where mutually exclusive.
               myprocess=multiprocessing.Process(target=combine_nc_files.FromFilelist,
                                                 args=(files_for_global,global_output_filename,None,
                                                       datadate,False))
               myprocess.start()
               myprocess.join()

            #Now we can "join" the files from all the years to get an "average" over the period
            #This replaces the old join_xy functionality
            joined_ncfile=os.path.join(outputdirectory['joined'],'%02d_joined.nc'%month)
            dofilesexist=combine_nc_files.Climatology(glob.glob("%s/*.nc"%(outputdirectory['reanalysedglobalmonth']%(month))),
                                          joined_ncfile,[fvar,pvar])
            #convert the file to ASCII also (for interpolation) - we only output the fCO2 and pCO2 variables
            if dofilesexist:
               netcdf_helper.write_netcdf_vars_to_ascii(joined_ncfile,joined_ncfile.replace(".nc",".txt"),[fvar,pvar])

         ##Now we will wrap all the global monthly files into a single file
         ##containing data for each month of each year
         ##We need a filelist ordered by time (e.g. 01 1991, 02 1991, 03 1991 ... 10 2010, 11 2010, 12 2010)
         #globalfilelist=[]
         #for month in range(1,13):
            #globpath="%s/*OCF-CO2-GLO-1M-100-SOCAT-CONV.nc"%outputdirectory['reanalysedglobalmonth']%(month)
            #globalfilelist.extend(glob.glob(globpath))
         ##this is ordered yearly by month. We need to change it to monthly by year
         #pathfilelist=[os.path.split(f) for f in globalfilelist]
         #sortedlist=[os.path.join(x,y) for x,y in sorted(pathfilelist,key=lambda pair: pair[1])]
         #combine_nc_files.IntoTimeDimension(sortedlist,
                                            #"%s/%s"%(outputdirectory['reanalysedglobal'],output_global_monthly_file))

      if exitatinterpolation:
         print("Have finished the pre-interpolation (reanalysis) stage.")
         print("If you are going to run with post-interpolation, please ensure that interpolated data files are named as expected and in the expected directories.")
         return 0;

      #Test which interpolation methods to use
      if gstatcmds is not None:
         #Now do the kriging
         interpmethod="GSTAT"
         KrigeWithGstat(gstatcmds)
      elif diva:
         interpmethod="DIVA"
         for filename in sorted(glob.glob(outputdirectory['joined']+"/*.txt")):
            destdir=outputdirectory['fpred']
            InterpolateWithDiva(filename,destdir,FCO2,divaparamdefault)
            destdir=outputdirectory['ppred']
            InterpolateWithDiva(filename,destdir,PCO2,divaparamdefault)
   
   #We do the following after interpolation
   #Now do the combining of files 
   if extrapolatetoyear is not None:
      output_prefix="%s/%4d"%(outputdirectory['combined'],extrapolatetoyear)
      tcl_filename_base="%s/%s/%4d%%02d01_OCF-SST-GLO-1M-100-ATS-ARC.nc"%(sstdir,extrapolatetoyear,extrapolatetoyear)
   else:
      output_prefix="%s/"%(outputdirectory['combined'])
      tcl_filename_base=None

   for month in range(1,13):
      if tcl_filename_base is not None:
         tcl_filename=tcl_filename_base%(month)
      else:
         tcl_filename=None
      if useaatsr:
         tcltype="aatsr"
      elif usereynolds:
         tcltype="reynolds"

      #Due to the file name of the vco2 data it is non-trivial to auto get for each month
      #but we can convert month to string representation and get the file with that in the name
      vco2_pattern="%s/*01%s011200*nc"%(vco2dir,datetime.datetime(2010,month,1).strftime("%b").lower())
      vco2_filename=glob.glob(vco2_pattern)
      if len(vco2_filename) == 0:
         # no file found
         print("Cannot find vco2 filename based on filename %s"%(vco2_pattern))
         vco2_filename=None
      elif len(vco2_filename) > 1:
         # more than one file found
         print("Multiple vco2 files found based on filename %s"%(vco2_pattern))
         print("Selecting first one: %s"%(vco2_filename[0]))
         vco2_filename=vco2_filename[0]
      else:
         # one file found
         vco2_filename=vco2_filename[0]

      combine_xy.DoTheCombining(fpred="%s/%02d_joined.txt.anl"%(outputdirectory['fpred'],month),
                        fvari="%s/%02d_joined.txt-error.anl"%(outputdirectory['fpred'],month),
                        ppred="%s/%02d_joined.txt.anl"%(outputdirectory['ppred'],month),
                        pvari="%s/%02d_joined.txt-error.anl"%(outputdirectory['ppred'],month),
                        tcl=tcl_filename,
                        vco2=vco2_filename,
                        output=output_prefix+"%02d_01_OCF-CO2-GLO-1M-100-INTRP-CLIM.nc"%(month),
                        plot=True,month=month,tclyear=str(extrapolatetoyear),
                        extrapolatedyear=extrapolatetoyear,
                        method=interpmethod,tcltype=tcltype)
      
      return 0;


if __name__ == "__main__":
    #Get command line variables
    cl = GetCommandline();
    
    exitCode = RunReanalyseSocat(**vars(cl));
    
    if exitCode == 0:
        print("Reanalysis of socat data completed successfully.");
    else:
        print("An error occured while performing reanalysis. Exited with exit-code:", exitCode);
    


