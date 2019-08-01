#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Tue Jul 30 10:14:16 2019
Driver for the lib_ofluxghg_flux_budgets tools (originally written by Peter Land, PML).
Calculates the net integrated flux budgets from monthly global or regional FluxEngine outputs.
"""
import argparse;
from fluxengine.tools.lib_ofluxghg_flux_budgets import run_flux_budgets;

if __name__ == "__main__":
     # default directories
    defaultDir = '/.../OceanFlux_example-climatology'
    defaultRef = '/.../OceanFlux_example-climatology-ref'
    defaultOutRoot = '/OceanFlux_output'
    defaultMaskfile = '/.../World_Seas-final-complete.nc'
    defaultLandfile = '/.../onedeg_land.nc'
    defaultFluxName = 'OF'
    defaultKwName = 'OK3'
    defaultCiName = 'OIC1'
    defaultCwName = 'OSFC' # NB mislabelled in the netCDF as atmospheric conc
    defaultIceName = 'P1'
    defaultLandName = 'land_proportion'
    defaultGridAreafile = 'no_file'
    defaultGridAreaName = 'area'
    
     # checking arguments
    parser = argparse.ArgumentParser(description='Calculate yearly and monthly '+
       'net fluxes in regions')
    parser.add_argument('-d','--dir', nargs='?', default=defaultDir, help=
       'directory to process, containing YYYY/MM/*.nc '+
       '(will object if any month directory contains more than one .nc file)')
    parser.add_argument('-R','--ref', nargs='?',
       help='optional reference directory, to be compared with DIR. '+
       'WARNING: setting this may change fluxes')
    parser.add_argument('-fd', '--fluxdataset', nargs='?', default=defaultFluxName,
       help='flux dataset name')
    parser.add_argument('-kwd', '--kwdataset', nargs='?', default=defaultKwName,
       help='gas transfer velocity dataset name')
    parser.add_argument('-cid', '--cidataset', nargs='?', default=defaultCiName,
       help='interfacial CO2 concentration dataset name')
    parser.add_argument('-cwd', '--cwdataset', nargs='?', default=defaultCwName,
       help='interfacial CO2 concentration dataset name')
    parser.add_argument('-id', '--icedataset', nargs='?', default=defaultIceName,
       help='ice dataset name')
    parser.add_argument('-o', '--outroot', nargs='?', default=defaultOutRoot,
       help="root of output .txt files containing processed data, one per region."+
       "Output files will be <outroot>_<region>.txt. If 'global' is not included, "+
       "<outroot>_global.txt will also be created")
    parser.add_argument('-mf', '--maskfile', nargs='?', default=defaultMaskfile,
       help='netCDF file containing mask information')
    parser.add_argument('-md', '--maskdatasets', nargs='+', default=[],
        help='netCDF mask dataset name, one per region (default REGIONS if '+
        'specified, otherwise global)')
    parser.add_argument('-r', '--regions', nargs='+', default=[],
       help='region names for output files (default MASKDATASETS if specified, '+
       'otherwise global)')
    parser.add_argument('-lf', '--landfile', nargs='?', default=defaultLandfile,
       help='netCDF file containing global land proportion')
    parser.add_argument('-ld', '--landdataset', nargs='?', default=defaultLandName,
       help='netCDF global land dataset name')
    parser.add_argument('-ga', '--gridarea', type=int, nargs='?', default=0,
       help='Grid area filename or value')
    parser.add_argument('-gaf', '--gridareafile', nargs='?', default=defaultGridAreafile,
       help='netCDF file containing area of grid cells')
    parser.add_argument('-gad', '--gridareadataset', nargs='?', default=defaultGridAreaName,
       help='netCDF grid area dataset name')
    parser.add_argument('-L', '--LooseIce', default=False,
       help='use Loose et al ice parameterization instead of Takahashi')
    parser.add_argument('-ip', '--icePercent', default=False,
       help='Ice coverage values are in % instead of 0-1')
    parser.add_argument('-w', '--window', type=int, nargs=4,
       help='image window to process (zero relative) '+
       '[startLine,endLine,startSample,endSample]')
    parser.add_argument('-p', '--places', type=int, nargs='?', default=10,
       help='maximum number of characters in string formatting of floats')
    parser.add_argument("-v", "--verbosity", action="count",
       help="increase output verbosity")
    args = parser.parse_args()
    
    run_flux_budgets(args);