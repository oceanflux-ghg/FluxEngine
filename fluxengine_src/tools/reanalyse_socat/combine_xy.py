#!/home/cercache/tools/environments/run_python_scibox_precise.sh
#-*- coding:utf-8 -*-

import csv
import numpy as np
import os
import sys
import netCDF4
from . import netcdf_helper
import argparse
import matplotlib.pyplot as plt

#NOTE THAT THERE ARE HARD CODED VALUES AND ATTRIBUTES THROUGHOUT THIS AND SHOULD BE MOVED TO BE ENTERED
#VIA COMMAND LINE TO BE MORE ROBUST. E.G. VCO2 DATA ATTRIBUTES AND SST VARIABLES

#Function to make a graph and write to file
def MakeGraph(filename,x,y,z,title="",xlabel="Longtitude",ylabel="Latitude",axis=[-180., 180., -90., 90.],scale=None):
   """
   Inputs: 
      filename - output filename to write to
      x,y,z - the X,Y and Data points that make the mesh
      title,xlabel,ylabel - the graph and axes labels
      axis - list of 4 floats describing [xmin,xmax,ymin,ymax] axis limits
      scale - the limits to use for the colour scale, tuple of two e.g. (5,50)
   """
   plt.figure(1)
   plt.pcolormesh(x, y, z)
   plt.axis(axis)
   plt.xlabel(xlabel)
   plt.ylabel(ylabel)
   plt.suptitle(title)
   plt.colorbar()
   if scale is not None:
      plt.clim(*scale)
   plt.axes().set_aspect('equal')
   plt.savefig(filename)
   plt.clf()

def GetCommandline():
   """
   Deal with command line parameters (if script is run from the command line)
   """
   parser = argparse.ArgumentParser(description="Combine multiple data sets into grids in same netCDF file.")
   #TODO/FIXME 2010 is used for Tcl and output prefix filename
   parser.add_argument('--output', metavar='<filename>',help ='The output file',default="./2010_%s_01_OCF-CO2-GLO-1M-100-KRG-CLIM.nc")
   parser.add_argument('--fpred', metavar='<filename>',help ='The input interpolated prediction of fCO2 file',default=None)
   parser.add_argument('--fvari', metavar='<filename>',help ='The input interpolated variation of fCO2 file',default=None)
   parser.add_argument('--ppred', metavar='<filename>',help ='The input interpolated prediction of pCO2 file',default=None)
   parser.add_argument('--pvari', metavar='<filename>',help ='The input interpolated variation of pCO2 file',default=None)
   parser.add_argument('--tcl', metavar='<filename>',help ='The input Tcl file',default=None)
   parser.add_argument('--tcltype', metavar='<keyname>',help ='The SST data type (e.g. reynolds, aatsr)',default="aatsr")
   parser.add_argument('--tclyear', type=str,metavar='<year>',help ='The year of the Tcl data file',default=None)
   parser.add_argument('--vco2', metavar='<filename>',help ='The input vCO2 file',default=None)
   parser.add_argument('--plot',action='store_true',help ='Create plots of the interpolated fCO2 data',default=False)
   parser.add_argument('--month', metavar='<string>',type=str,help ='The month to insert into plot title - defaults to filename',default=None)
   parser.add_argument('--extrapolatedyear',metavar='year',type=int,help="The year (if any) that data has been extrapolated to",default=None)
   parser.add_argument('--interpmethod', metavar='<string>',type=str,help ='The interpolation method used: DIVA or GSTAT',default=None,required=True)
   commandline=parser.parse_args()

   #expand linux user if present
   commandline.output = os.path.expanduser(commandline.output)
   #expand users on input files
   paths=[commandline.fpred,commandline.ppred,commandline.fvari,commandline.pvari,commandline.tcl,commandline.vco2]
   for path in paths:
      if path is not None:
         path=os.path.expanduser(path)

   if commandline.tclyear is None:
      print("Error! Must supply the year the Tcl data is relevant for.")
      sys.exit(1)

   if commandline.tcltype not in ["aatsr","reynolds"]:
      print("Error! tcltype must be either aatsr or reynolds")
      sys.exit(1)

   return commandline

def CheckFilesExist(fpred,ppred,tcl,vco2,fvari,pvari):
   """
   Check if given files exist and throw exception if they don't
   """
   vardict={'fpred': fpred,'ppred': ppred,'vco2':vco2}
   for fnamekey in list(vardict.keys()):
      if vardict[fnamekey] is None:
         raise Exception("Must specify filename for: %s"%fnamekey)
      if not os.path.exists(vardict[fnamekey]):
         raise Exception("Cannot open file as it does not exist: %s"%(vardict[fnamekey]))

   #Now check optional files
   inputlist=[fvari,pvari,tcl]
   for fname in inputlist:
      if fname is not None:
         if not os.path.exists(fname):
            raise Exception("Cannot open file as it does not exist: %s"%(fname))

def LoadKrigedAscii(filename):
   """
   Function to load in the data from a kriged ASCII file as output from GSTAT
   First 6 lines are header, e.g.:
      NCOLS              360
      NROWS              180
      XLLCORNER -180.008333333
      YLLCORNER          -90
      CELLSIZE  1.0000462963
      NODATA_VALUE     -9999
   """
   with open(filename,'r') as filein:
      for line in filein.readlines():
         if line.find("NODATA_VALUE")!=-1:
            nodatavalue=line.split()[1]
            break
   #Read into a numpy array skipping the 6 line header
   data = np.loadtxt(filename, skiprows=6)
   data = data[::-1,:] #rows in reverse order
   data[np.where(data==float(nodatavalue))] = netcdf_helper.MISSINGDATAVALUE
   return data

def LoadDIVAAscii(filename,nodatavalue=-999):
   """
   Function to load in the ASCII output from the DIVA interpolated files
   First row is row index
   First column is column index
   Then repeats in blocks defined by the grid size
   """
   data=np.zeros([180,360])
   with open(filename,'r') as csvfile:
      reader=csv.reader(csvfile,delimiter=' ',skipinitialspace=True)
      try:
         #repeat for all lines in the file
         #get the row indices
         rindex=next(reader)
         #convert to floats
         rindex=list(map(float,rindex))
         cindex=[]
         while True:
            #update rindex to be the last cindex read in from loop below
            if cindex != []:
               rindex=cindex
            #set the previous column index = 0
            prev_cindex=0
            #read in the next row (column index + data)
            cindex=next(reader)
            #convert to floats
            cindex=list(map(float,cindex))
            #for each block of data - i.e. until col index is not increasing
            while(float(cindex[0])==prev_cindex+1):
               prev_cindex=cindex[0]
               #insert the data in to the data array
               for i in range(1,len(rindex)):
                  data[rindex[i]-1,cindex[0]-1]=float(cindex[i])
               #get the next row from the file
               cindex=next(reader)
               #convert to floats
               cindex=list(map(float,cindex))
      except StopIteration:
         #this occurs if try to read past end of file
         pass
      except Exception as err:
         #catch any other errors
         raise Exception("An error occurred whilst reading %s: Due to error: %s"%(filename,err.message))

   #data = np.flipud(data) #rows in reverse order
   data[np.where(data==nodatavalue)] = netcdf_helper.MISSINGDATAVALUE
   return data

def DoTheCombining(fpred,ppred,tcl,vco2,output,fvari=None,pvari=None,plot=False,month=None,
                   extrapolatedyear=None,tclyear=None,method="DIVA",tcltype=None):
   """
   Main function that does all the work - call this if importing as a library.
      fpred - interpolated fCO2 prediction filename
      ppred - interpolated pCO2 prediction filename
      tcl - atsr climatology filename
      vco2 - vCO2 (from Takahashi) filename
      output - the combined netcdf output filename
      fvari - fCO2 variance filename (output from interpolation)
      pvari - pCO2 variance filename (output from interpolation)
      plot - boolean, whether to create plots of the data (True or False)
      month - the month of which the data refers to - only required for title of plots
      extrapolatedyear - the year data has been extrapolated to 
      tclyear - the year of the ATSR data (should be the same as extrapolatedyear?)
      method - the interpolation method that has been used (DIVA, GSTAT)
   """
   #First check if all input files exist
   CheckFilesExist(fpred,ppred,tcl,vco2,fvari,pvari)

   if method == "DIVA":
      DataLoader=LoadDIVAAscii
   elif method == "GSTAT":
      DataLoader=LoadKrigedAscii
   else:
      raise Exception("Unrecognised interpolation method in combine function: %s"%method)

   f_pred=DataLoader(fpred)
   p_pred=DataLoader(ppred)
   
   if fvari is None:
      f_std = np.array([[netcdf_helper.MISSINGDATAVALUE] * 360] * 180)
   else:
      f_vari = DataLoader(fvari)
      if method == "GSTAT":
         f_std = np.where(f_vari>=0,f_vari**0.5,netcdf_helper.MISSINGDATAVALUE)
         #f_std[np.where(f_vari==-9999.)] = netcdf_helper.MISSINGDATAVALUE
      elif method == "DIVA":
         f_std=f_vari

   if pvari is None:
      p_std = np.array([[netcdf_helper.MISSINGDATAVALUE] * 360] * 180)
   else:
      p_vari = DataLoader(pvari)
      if method == "GSTAT":
         p_std = np.where(p_vari>=0,p_vari**0.5,netcdf_helper.MISSINGDATAVALUE)
         #p_std[np.where(p_vari==-9999.)] = netcdf_helper.MISSINGDATAVALUE
      elif method == "DIVA":
         p_std=p_vari

   # Get Tcl in 'tclyear' for each cell
   if tcl is not None:
      if tcltype == "aatsr":
         sstkeyname="sst_skin_mean"
         sstdataname="ARC_ATSR"
      elif tcltype == "reynolds":
         sstkeyname="sst_mean"
         sstdataname="Reynolds"
      else:
         raise "Unknown SST data type: %s"
      with netCDF4.Dataset(tcl) as sst_file:
         Tcl = sst_file.variables[sstkeyname]
         Tcl = Tcl[0,:,:]
         if sst_file.variables[sstkeyname].units in ['degrees C']:
            #we want them in Kelvin (to be consistent with how the scripts were originally written)
            Tcl=Tcl+273.15
         elif sst_file.variables[sstkeyname].units not in ['Kelvin','kelvin','K']:
            print("Unsure of data units for temperature - aassuming kelvin.")
   else:
      Tcl = None

   # Get vCO2  for each cell
   with netCDF4.Dataset(vco2) as vCO2_file:
      vCO2_data = vCO2_file.variables['vCO2_2010_grid_up']
      vCO2_data = vCO2_data[0,:,:]
   #Lon/lat values for dimension arrays in netCDF
   xi = np.linspace(-179.5,179.5,360)   
   yi = np.linspace(-89.5, 89.5, 180)

   #Write out the data to the netCDF file 
   #Test directory exists
   if not os.path.exists(os.path.dirname(output)) and os.path.dirname(output) != "":
      raise Exception("Directory to write file to does not exist: %s"%(os.path.dirname(output)))

   if extrapolatedyear is not None:
      extrapyear=str(extrapolatedyear)+'_'
      nameext=" extrapolated to %d"%extrapolatedyear
   else:
      extrapyear=""
      nameext=""

   with netCDF4.Dataset(output, 'w', format = 'NETCDF4') as ncfile:
      #Add standard dimensions and arrays
      netcdf_helper.standard_setup_SOCAT(ncfile,timedata=1e9,londata=xi,latdata=yi)
      #Now add other grid arrays
      fCO2_interp_pred = ncfile.createVariable('fCO2_'+extrapyear+'interpolated_pred',
                                               'f4',('time','latitude','longitude'),
                                               fill_value=netcdf_helper.MISSINGDATAVALUE,zlib=True)
      fCO2_interp_pred[:] = f_pred
      fCO2_interp_pred.units = 'uatm'
      fCO2_interp_pred.missing_value = netcdf_helper.MISSINGDATAVALUE
      fCO2_interp_pred.valid_min = 0.
      fCO2_interp_pred.valid_max = 1e6
      fCO2_interp_pred.scale_factor = 1.
      fCO2_interp_pred.add_offset = 0.
      fCO2_interp_pred.standard_name = "fCO2_"+extrapyear+"interpolated_pred"
      fCO2_interp_pred.long_name = "Fugacity of CO2 using SOCAT methodology"+nameext+" and interpolated"
      fCO2_interp_error = ncfile.createVariable('fCO2_'+extrapyear+'interpolated_error',
                                                'f4',('time','latitude','longitude'),
                                                fill_value=netcdf_helper.MISSINGDATAVALUE,zlib=True)
      fCO2_interp_error[:] = f_std
      fCO2_interp_error.units = 'uatm'
      fCO2_interp_error.missing_value = netcdf_helper.MISSINGDATAVALUE
      fCO2_interp_error.valid_min = 0.
      fCO2_interp_error.valid_max = 1e6
      fCO2_interp_error.scale_factor = 1.
      fCO2_interp_error.add_offset = 0.
      fCO2_interp_error.standard_name = "fCO2_"+extrapyear+"interpolated_std"
      fCO2_interp_error.long_name = "Fugacity of CO2 using SOCAT methodology"+nameext+" error of interpolation"
      pCO2_interp_pred = ncfile.createVariable('pCO2_'+extrapyear+'interpolated_pred',
                                               'f4',('time','latitude','longitude'),
                                               fill_value=netcdf_helper.MISSINGDATAVALUE,zlib=True)
      pCO2_interp_pred[:] = p_pred
      pCO2_interp_pred.units = 'uatm'
      pCO2_interp_pred.missing_value = netcdf_helper.MISSINGDATAVALUE
      pCO2_interp_pred.valid_min = 0.
      pCO2_interp_pred.valid_max = 1e6
      pCO2_interp_pred.scale_factor = 1.
      pCO2_interp_pred.add_offset = 0.
      pCO2_interp_pred.standard_name = "pCO2_"+extrapyear+"interpolated_pred"
      pCO2_interp_pred.long_name = "Partial pressure of CO2 using SOCAT methodology"+nameext+" and interpolated"
      pCO2_interp_error = ncfile.createVariable('pCO2_'+extrapyear+'interpolated_error',
                                                'f4',('time','latitude','longitude'),
                                                fill_value=netcdf_helper.MISSINGDATAVALUE,zlib=True)
      pCO2_interp_error[:] = p_std
      pCO2_interp_error.units = 'uatm'
      pCO2_interp_error.missing_value = netcdf_helper.MISSINGDATAVALUE
      pCO2_interp_error.valid_min = 0.
      pCO2_interp_error.valid_max = 1e6
      pCO2_interp_error.scale_factor = 1.
      pCO2_interp_error.add_offset = 0.
      pCO2_interp_error.standard_name = "fCO2_"+extrapyear+"interpolated_error"
      pCO2_interp_error.long_name = "Partial pressure of CO2 using SOCAT methodology"+nameext+" error of interpolation"
      if tcl is not None:
         Tcl_year = ncfile.createVariable('Tcl_'+tclyear,'f4',('time','latitude','longitude'),
                                          fill_value=netcdf_helper.MISSINGDATAVALUE,zlib=True)
         Tcl_year[:] = Tcl
         Tcl_year.units = 'Kelvin'
         Tcl_year.missing_value = netcdf_helper.MISSINGDATAVALUE
         Tcl_year.valid_min = 0.
         Tcl_year.valid_max = 1e6
         Tcl_year.scale_factor = 1.
         Tcl_year.add_offset = 0.
         Tcl_year.standard_name = "Tcl_"+tclyear
         Tcl_year.long_name = "Climatologial temperature from "+sstdataname+" in "+tclyear 
      vCO2 = ncfile.createVariable('vCO2','f4',('time','latitude','longitude'),
                                   fill_value=netcdf_helper.MISSINGDATAVALUE,zlib=True)
      vCO2[:] = vCO2_data
      vCO2.units = 'ppm'
      vCO2.missing_value = netcdf_helper.MISSINGDATAVALUE
      vCO2.valid_min = 0.
      vCO2.valid_max = 1e6
      vCO2.scale_factor = 1.
      vCO2.add_offset = 0.
      vCO2.standard_name = "vCO2"
      vCO2.long_name = "Concentration of CO2 in dry air in 2000 from NOAA ESRL"
      vCO2.data_citation="Dlugokencky, E.J., K.A. Masarie, P.M. Lang, and P.P. Tans (2014), NOAA Greenhouse Gas Reference from Atmospheric Carbon Dioxide Dry Air Mole Fractions from the NOAA ESRL Carbon Cycle Cooperative Global Air Sampling Network. Data Path: ftp://aftp.cmdl.noaa.gov/data/trace_gases/co2/flask/surface/."
      vCO2.data_contacts="Ed.Dlugokencky_AT_noaa.gov, Pieter.Tans_AT_noaa.gov"

   if plot is True:
      #if month has not been passed then use the output netcdf filename in the graph title
      if month is None:
         month=output
      #create a plot of the predicted fCO2
      fpredfilename=output+'_pred_fCO2_'+extrapyear+'.png'
      fpredtitle = 'pred fCO2 (uatm), %s, year %s'%(month,extrapyear)
      goodindices=np.where(~np.isnan(f_pred)&(f_pred!=netcdf_helper.MISSINGDATAVALUE))
      scale=(f_pred[goodindices].min(),f_pred[goodindices].max())
      #scale=(200,600)
      MakeGraph(filename=fpredfilename,x=xi,y=yi,z=f_pred,title=fpredtitle,scale=scale)

      #create a plot of the standard deviation of fCO2 prediction if it has been included in the combination netcdf file
      if fvari is not None:
         fvarfilename=output+'_std_fCO2_'+extrapyear+'.png'
         fvartitle = 'std fCO2 (uatm), %s, year %s'%(month,extrapyear)
         MakeGraph(filename=fvarfilename,x=xi,y=yi,z=f_std,title=fvartitle,scale=(5,50))


def Main():
   commandline=GetCommandline()
   DoTheCombining(commandline.fpred,commandline.ppred,commandline.tcl,commandline.vco2,commandline.output,
                  commandline.fvari,commandline.pvari,commandline.plot,commandline.month,commandline.extrapolatedyear,
                  tclyear=commandline.tclyear,tcltype=commandline.tcltype,method=commandline.interpmethod)

if __name__=='__main__':
   Main()
