#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Commandline driver script for ncdf2text.
Converts a netCDF output file produced by FluxEngine into a text file with one
column for each data layer and one row for each non-missing value grid cell.

Created on 2019-07-31

@author: tom holding
"""

import argparse;
from fluxengine.tools.lib_ncdf2text import convert_ncdf_to_text;

def parse_cl_arguments():
    description = str("""Converts netCDF3 data produced as output from FluxEngine into text data (e.g. csv, tsv).
    The text output will contain one column for each data layer, and one row for each cell of the grid (if there is no missing data).
    One text file will be created for each netCDF3 file.
    """, 'utf-8');
    
    parser = argparse.ArgumentParser(description=description);
    parser.add_argument("inFiles", nargs="+",
                        help="list of paths to netCDF file(s) containing output from FluxEngine. Standard Unix glob patterns are supported.");
    parser.add_argument("-n", "--textOutPath", default="",
                        help="Path to store the produced text file(s). This should be a directory. Text files will be named after the netCDF file used to generate them. Defaults to the current working directory.");
    parser.add_argument("--delim", default="\t",
                        help="delimiter token used to seperate data in the input text file. Default is a tab.")
    parser.add_argument("-c", "--commentChar", default="#",
                        help="Character prefix which indicates a comment. Default is a '#'. Set to an empty string to turn comments off.");
    parser.add_argument("-m", "--missingValue", default="nan",
                        help="Value used to indicate a missing value in the input text file. Default is 'nan'.");
    parser.add_argument("--encoding", default="utf-8",
                        help="Encoding of the input text file. Default it 'utf-8'");
    clArgs = parser.parse_args();
    
    return clArgs.inFiles, clArgs.textOutPath, clArgs.delim, clArgs.missingValue, clArgs.encoding;


if __name__ == "__main__":
    #Parse commandline arguments
    inFiles, textOutPath, delim, commentChar, missingValue, encoding = parse_cl_arguments();
    
    #Run
    convert_ncdf_to_text(inFiles=inFiles, textOutPath=textOutPath, delim=delim, commentChar=commentChar, missingValue=missingValue, encoding=encoding);
    










