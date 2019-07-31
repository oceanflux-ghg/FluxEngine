#! /usr/bin/env python
# -*- coding: utf-8 -*-
#Created 2019/07/30
#Simple driving script for text2ncdflib

import argparse;
from fluxengine.tools.lib_text2ncdf import convert_text_to_netcdf;

#Parse and set variables from command line arguments. Returns an object containing each command line argument as an attribute.
def parse_cl_arguments():
    description = u"""Converts text encoded data (e.g. csv, tsv) to netCDF3 data which is compatible for use with FluxEngine.
    The text input file must have a header containing the column names.
    Column names must start with a letter (a-z or A-Z) and can only contain letters, numbers, spaces, underscores and certain symbols.
    For information on allowed symbols consult the netCDF3 documentation.
    Symbol requirements do not apply to units (e.g. "Temp [°C]" or "fCO2 [µatm]" are valid column name if unit parsing is on.)""";
    
    parser = argparse.ArgumentParser(description=description);
    parser.add_argument("inFiles", nargs="+", help="list of paths to input text file(s) containing the data you want to use. These must all have the same header configuration and column names. Standard Unix glob patterns are supported.");
    parser.add_argument("-n", "--ncOutPath", help="Path to the output netCDF file(s). Note that if more than one netCDF file will be created a date/time stamp will be appended to this filename indicating the start of this temporal 'bin'.");
    parser.add_argument("-s", "--startTime",
        help="start date (and time). Format: YYYY[[[-MM-DD] hh:mm]:ss]. If only the year is supplied, the first second of the year is used.");
    parser.add_argument("-e", "--endTime",
        help="end date (and time). Format: YYYY[[[-MM-DD] hh:mm]:ss]. If only the year is supplied, the first second of the year is used.");
    parser.add_argument("-t", "--temporalResolution", default="monthly",
        help="Temporal resolution to use in the output file. Defaults to monthly. Format required (days hours:minutes): D hh:mm");
    parser.add_argument("-k", "--temporalChunking", type=int, default=1,
        help="Temporal chunking: How many time points to store in each file. Defaults to one time point per file.");
    
    parser.add_argument("--latProd", default="Latitude",
        help="latitude product name. This should match the column name in the input text file. Default is 'Latitude'");
    parser.add_argument("--lonProd", default="Longitude",
        help="longitude product name. This should match the column name in the input text file. Default is 'Longitude'");
    parser.add_argument("--latResolution", default=1.0, type=float, help="spatial resolution of the grid (latitude)");
    parser.add_argument("--lonResolution", default=1.0, type=float, help="spatial resolution of the grid (longitude)");
    parser.add_argument("--delim", default="\t",
        help="delimiter token used to seperate data in the input text file. Default is a tab.")
    parser.add_argument("-d", "--dateIndex", type=int, default=0,
        help="The column number (starting from 0) which contains the date/time field in your input text file. Default is a 0.");
    parser.add_argument("-c", "--numCommentLines", nargs="+", type=int, default=[0],
        help="Number of comment lines before the header to be ignored). Default is a 0.");
    parser.add_argument("--cols", nargs="+", default=None,
        help="List of column names or numbers which will be added to the output netCDF file. This should not include longitude, latitude or data/time columns. Default will add all suitable columns.");
    parser.add_argument("-m", "--missing_value", default="nan", help="Indicates a missing value in the input text file. Default is 'nan'.");
    parser.add_argument("-u", "--parse_units", action="store_true",
                          help="will automatically parse units in the header column names if follow the variable name and are formatted between square brackets (e.g. Depth [m])");
    parser.add_argument("--encoding", default="utf-8", help="Encoding of the input text file. Default it 'utf-8'");
    parser.add_argument("-l", "--limits", type=float, nargs = 4, default = [-90.0, 90.0, -180.0, 180],
                        help="Coordinates which define the grid limits given in latitude and longitude: South North West East, e.g. -45 -30 -45 10 Defaults to -90 90 -180 180")
    
    parser.add_argument("--dateFormatDayFirst", help="Specifies that the date/time column in the input data is formatted day-first (international/European format).", action="store_true", default=False);

    clArgs = parser.parse_args();
    
    return clArgs;

if __name__ == "__main__":
    #parse command line arguments
    print("Parsing command line arguments.");
    args = parse_cl_arguments();
    
    convert_text_to_netcdf(args.inFiles, args.startTime, args.endTime, args.ncOutPath,
                           limits=args.limits, latResolution=args.latResolution, lonResolution=args.lonResolution,
                           temporalResolution=args.temporalResolution, temporalChunking=args.temporalChunking,
                           delim=args.delim, numCommentLines=args.numCommentLines, encoding=args.encoding,
                           textFileMissingValue=args.missing_value, parseUnits=args.parse_units,
                           dateIndex=args.dateIndex, colNames=args.cols, latProd=args.latProd, lonProd=args.lonProd, dateFormatDayFirst=args.dateFormatDayFirst); 

