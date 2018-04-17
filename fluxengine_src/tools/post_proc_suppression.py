#Script to read FluxEngine output and apply a variable suppression factor to the OF

import numpy as np,csv, argparse
from netCDF4 import Dataset
import sys
import os
import calendar
import shutil
sys.path.insert(0, '/Users/iga202/')

missing_value = -999.
# read arguments
parser = argparse.ArgumentParser(description = 'Apply variable suppression factor to FE results')
parser.add_argument('-p','--path', default = 'Documents/1RF/surfactants/Mar2017/Results/', action='store',
   nargs='*', help="Path for FE results")
parser.add_argument('-y','--year', default = 2014, action='store',
   nargs='*', help="Year for calculation - default 2014")
# parser.add_argument('-o','--output', default = None, action='store',
#    nargs='*', help="Output fileName. Default - the same as the input name but with _supp at the end")
parser.add_argument('-eq','--equation',default = 'power', action='store', nargs='*',
   help="Equation choice either 'linear' or 'power'")
args = parser.parse_args()
path = args.path[0]
year = args.year
equation = args.equation[0]

#path = 'Ifremer/OCEANFLUX-SHARED/workspace/i.g.c.ashton_exeter.ac.uk/'
a = [i for i,name in enumerate(os.listdir(path)) if name==str(year)]
if a:
   inpath = path+'/'+str(year)+'/'
else:
   print('No year directories found in ', inpath)
   sys.exit(1)
for rn in range(2,21):
  for month in range(1,13):
     mnth = '%02d' % month
     outpath = path+str(rn)+'/'+str(year)+'/'+mnth
     if not os.path.exists(outpath):
       os.makedirs(outpath)
     fileName = inpath+mnth+'/OceanFluxGHG-month'+mnth+'-'+calendar.month_abbr[month].lower()+'-'+str(year)+'-v0.nc'#Define the filename according to year, month and naming convention
     outFile = outpath+'/OceanFluxGHG-month'+mnth+'-'+calendar.month_abbr[month].lower()+'-'+str(year)+'-v'+str(rn)+'.nc'#Define the filename according to year, month and version 
     shutil.copy2(fileName,outFile)#Create output file as a copy of input file
     print fileName
     ncin = Dataset(outFile,'a')#open the results file for appending
     ncin.setncattr("_Fill_Value", missing_value)
     OFin = ncin['OF'][:]
     OFin[OFin.mask] = missing_value
     SST = ncin['ST1_mean'][:]
     SST[SST.mask] = missing_value
     t,nx,ny = OFin.shape
     #print nx,ny
     
     OFin_f = np.ravel(OFin)
     SST_f = np.ravel(SST)
    #    IGA - Adding surfactant treatment ----------------------------------------------------------------------------START
     supp = np.array([missing_value] * nx*ny)
     supp_fac = np.ravel(supp)
     supp_fac = supp_fac.astype(float)
     #print 'shapes supp, OF, SST',supp_fac.shape,OFin_f.shape,SST_f.shape
     for i in np.arange(nx * ny):
       if not any([OFin_f[i] == missing_value,SST_f[i] == missing_value]):
        if equation.lower() == 'linear':
         #Linear
         supp_fac[i] = 1.-(((1.18*SST_f[i])-11.04)/100.)#IGA - model
         #k_surf[i] = (2.108*sstskinC_fdata[i])-9.94#IGA - +95%
         #k_surf[i] = (0.2522*sstskinC_fdata[i])-32.03#IGA - +95%
        elif equation.lower() == 'power':
         #Power model
         try: 
           supp_fac[i] = 1.-(0.0046*pow(SST_f[i],2.5673)/100.)#IGA - power model
         except:
           supp_fac[i] = 0.#IGA - avoids complex numbers when SST<0 (This is also outside of the range of measurements taken at this stage)
           #print supp_fac
    #IGA - Adding errors in surfactant ------------------------------------------------------------------------------START
        if rn>1.:
            if equation.lower() == 'power':
              supp_fac[i] = supp_fac[i]+np.random.normal(0.,0.066486,1)
            elif equation.lower() == 'linear':
              supp_fac[i] = supp_fac[i]+np.random.normal(0.,0.070684,1)

    #IGA - Adding errors in surfactant ------------------------------------------------------------------------------END
        if (supp_fac[i] >1):
           supp_fac[i] = 1.#No error values
        if (supp_fac[i] <0):
           supp_fac[i] = 0.#No negative values
        OFin_f[i] = OFin_f[i]*supp_fac[i]
     else:
        OFin_f[i] = missing_value
     #IGA - Adding surfactant treatment ------------------------------------------------------------------------------END

     #Write output to file
     OFin_f.shape = (nx, ny)

     OFout = ncin.createVariable('OF_suppression',ncin.variables['OF'].datatype,(ncin.variables['OF'].dimensions))
     OFout.units = ncin.variables['OF'].units
     OFout.long_name = ncin.variables['OF'].long_name
     OFout.missing_value = missing_value
     OFin_f[np.isnan(OFin_f)]=missing_value
     OFin_f[abs(OFin_f)>990]=missing_value
     OFout[:] = OFin_f[:]

     supp_fac.shape = (nx, ny)
     supp_out = ncin.createVariable('suppression_fac',ncin.variables['OF'].datatype,(ncin.variables['OF'].dimensions))
     supp_out.units = 'None'
     supp_out.long_name = 'Suppression factor aplied to flux'
     supp_out.missing_value = missing_value
     supp_fac[np.isnan(supp_fac)]=missing_value
     supp_fac[abs(supp_fac)>990]=missing_value
     supp_out[:] = supp_fac[:]

