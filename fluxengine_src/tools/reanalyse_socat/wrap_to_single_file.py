#!/home/cercache/tools/environments/run_python_scibox_precise.sh
#-*- coding:utf-8 -*-

from . import combine_nc_files
import glob
import os
import argparse

fileglob="OCF-CO2-GLO-1M-100-SOCAT-CONV.nc"

def WrapToSingleFile(inputdir,outputfile):
   #Now we will wrap all the global monthly files into a single file
   #containing data for each month of each year
   #We need a filelist ordered by time (e.g. 01 1991, 02 1991, 03 1991 ... 10 2010, 11 2010, 12 2010)
   globalfilelist=[]
   for month in range(1,13):
      globpath="%s/%02d/*%s"%(inputdir,month,fileglob)
      globalfilelist.extend(glob.glob(globpath))
   #this is ordered yearly by month. We need to change it to monthly by year
   pathfilelist=[os.path.split(f) for f in globalfilelist]
   sortedlist=[os.path.join(x,y) for x,y in sorted(pathfilelist,key=lambda pair: pair[1])]
   combine_nc_files.IntoTimeDimension(sortedlist,outputfile)

if __name__ == "__main__":
   parser = argparse.ArgumentParser(description="Wrap per month netCDF files into a single file.")
   parser.add_argument('--output', metavar='<filename>',type=str,help ='The output file',required=True)
   parser.add_argument('--inputdir',metavar='<path>',type=str,help="The input directory for the month dirs",default=None)
   commandline=parser.parse_args()

   WrapToSingleFile(commandline.inputdir,commandline.output)