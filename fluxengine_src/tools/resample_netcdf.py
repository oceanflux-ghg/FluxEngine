#! /usr/bin/env python
'''resamples 1x1 degree netcdfs to 4 lat x 5 lon (Takahashi)'''

import sys, glob, os, argparse, csv, numpy as np, netCDF4 as nc

 # checking arguments
parser = argparse.ArgumentParser(description='Resample 1x1 degree netCDFs to 4 lat x 5 lon (Takahashi)')
parser.add_argument('infiles', nargs = '?',
    help='file or files (with wildcards) to resample')
group = parser.add_mutually_exclusive_group()
group.add_argument('-o','--outfile', nargs = '?', help='output file name')
group.add_argument('-e','--extension', nargs = '?', default = '_Takahashi',
   help='string to append to file name')
parser.add_argument('-c','--csvfile', action = 'store_true',
   help="create csv file with extension '.csv' instead of netCDF")
args = parser.parse_args()
infiles = args.infiles
outfilename = args.outfile
extension = args.extension
csvfile = args.csvfile

for infilename in glob.glob(sys.argv[1]):
   print("INPUT FILE:",infilename)
   infile = nc.Dataset(infilename, 'r')
   if outfilename == None:
      filename, extension2 = os.path.splitext(infilename)
      if csvfile:
         extension2 = '.csv'
      outfilename = filename + extension + extension2
   print("OUTPUT FILE:",outfilename)
   if csvfile:
      outfile = open(outfilename, 'wb')
      writer = csv.writer(outfile)
      headings = []
      outdata = np.ma.empty([0,3240])
   else:
      outfile = nc.Dataset(outfilename, 'w', format='NETCDF3_CLASSIC',
         clobber = True)
      for attname in infile.ncattrs():
         setattr(outfile, attname, getattr(infile,attname))
      for dimname, diminst in list(infile.dimensions.items()):
         print("copying dimension ", dimname)
         if dimname == 'latitude':
            if len(diminst) != 180:
               print("(resample_netcdf, main) latitude length", len(diminst))
               infile.close()
               outfile.close()
               os.remove(outfilename)
               sys.exit(1)
            outfile.createDimension('latitude', 45)
         elif dimname == 'longitude':
            if len(diminst) != 360:
               print("(resample_netcdf, main) longitude length", len(diminst))
               infile.close()
               outfile.close()
               os.remove(outfilename)
               sys.exit(1)
            outfile.createDimension('longitude', 72)
         else:
            outfile.createDimension(dimname, len(diminst))
   for varname, varinst in list(infile.variables.items()):
      print("copying variable ", varname)
      if csvfile:
         if varname == 'time' or varname == 'latitude_longitude':
            continue
         headings.append(varname)
         writer.writerow([varname, varinst.units])
         if varname == 'latitude':
            var = np.repeat(np.reshape(np.arange(88., -89., -4.), [45, 1]), 72,
               axis = 1)
         elif varname == 'longitude':
            var = np.repeat(np.reshape(np.arange(-177.5, 180., 5.), [1, 72]),
               45, axis = 0)
         else:
            if varinst.shape != (1, 180, 360):
               print('(resample_netcdf, main) odd shape ', varname, varinst.shape)
               infile.close()
               outfile.close()
               os.remove(outfilename)
               sys.exit(1)
            indata = varinst[:]
            var = np.ma.resize(np.ma.array([0.], mask = True), (45, 72))
            for j in range(45):
               j0=j*4
               j1=j0+4
               for i in range(72):
                  var[j, i] = np.ma.mean(indata[0, j0:j1, i*5:i*5+5])
         outdata = np.ma.vstack((outdata, np.ma.ravel(var)))
      else:
         atts = varinst.ncattrs()
         fill_found = '_FillValue' in atts
         if fill_found:
            var = outfile.createVariable(varname, varinst.dtype,
               varinst.dimensions, fill_value = varinst.getncattr('_FillValue'))
            atts.remove('_FillValue')
         else:
            var = outfile.createVariable(varname, varinst.dtype,
               varinst.dimensions)
         for attname in atts:
            # python converts to integer if these are present
            if not (attname in ['scale_factor', 'add_offset']):
               setattr(var, attname, getattr(varinst,attname))
         if varname == 'latitude':
            var[:] = np.arange(88., -89., -4.)
         elif varname == 'longitude':
            var[:] = np.arange(-177.5, 180., 5.)
         elif varname == 'time' or varname == 'latitude_longitude':
            var[:] = varinst[:]
         else:
            if varinst.shape != (1, 180, 360):
               print('(resample_netcdf, main) odd shape ', varname, varinst.shape)
               infile.close()
               outfile.close()
               os.remove(outfilename)
               sys.exit(1)
            indata = varinst[:]
            outdata = np.ma.resize(np.ma.array([0.], mask = True), (1, 45, 72))
            for j in range(45):
               j0=j*4
               j1=j0+4
               for i in range(72):
                  outdata[0, j, i] = np.mean(indata[0, j0:j1, i*5:i*5+5])
            var[:] = outdata
   infile.close()
   if csvfile:
      writer.writerow(headings)
      writer.writerows(outdata.T)
   outfile.close()
