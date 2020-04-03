# README #

Code to reanalyse SOCAT partial pressure using OceanFlux GHG methodology 

# Contents #

**reanalyse_socat_v2.py**: Executable script to control the reanalysis [run this]

**get_sst**: code to read monthly climatology files and extract data closest to 
   the ship position

**combine_nc_files.py**: code to combine arrays from netcdf files into new netcdf files

**combine_xy.py**: code to combine data files after spatial interpolation into a 
   single netCDF file per month

**wrap_to_single_file.py**: script to run at end that will wrap all global files
   into a single file with 3D arrays.

**datenum.py**: code to replicate the matlab datenum function

**netcdf_helper.py**: common functions used when dealing with the netcdf files

**v2_convert_f_2010_SSH.py**: code to read in from the SOCAT files and write to 
   netCDF

**v2_f_conversion_2010.py**: code to recalculate the CO2 flux


