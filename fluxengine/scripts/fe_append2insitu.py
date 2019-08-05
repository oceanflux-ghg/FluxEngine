#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Commandline driver for the append2insitu tool.

Converts a netCDF output file produced by FluxEngine into a text file with one
column for each data layer and one row for each non-missing value grid cell.

Created on 2019-07-30

@author: tom holding
"""

import argparse;
from fluxengine.tools.lib_append2insitu import append_to_in_situ;

#Parsed commandline arguments and returns an argparse.Namespace object.
def parse_cl_arguments():
    description = """Converts netCDF3 data produced as output from FluxEngine into text data (e.g. csv, tsv).
    The text output will contain one column for each data layer, and one row for each cell of the grid (if there is no missing data).
    One text file will be created for each netCDF3 file.
    """;
    
    parser = argparse.ArgumentParser(description=description);
    parser.add_argument("feOutput",
                        help="path to a FluxEngine output file");
    parser.add_argument("insituData",
                        help="path to a plain text in situ data file");
    parser.add_argument("outPath",
                        help="Path to store the produced text file.");
    parser.add_argument("--varsToAppend", nargs="+", default=["OF", "OK1"],
                        help="List of netCDF variable names which will be appended to the in situ data file. Default is the air-sea gas flux (OF) and the gas transfer velocity (K)");
    parser.add_argument("--delim", default="\t",
                        help="delimiter token used to seperate data in the in situ data file. Default is a tab.");
    parser.add_argument("--latCol", default="Latitude",
                        help="Name or column number (starting from 0) of the latitude column in the in situ input file. Default is a 'latitude'.");
    parser.add_argument("--lonCol", default="Longitude",
                        help="Name or column number (starting from 0) of the longitude column in the in situ input file. Default is 'longitude'.");
    #parser.add_argument("-c", "--commentChar", default="#",
    #                    help="Character prefix which indicates a comment. Default is a '#'. Set to an empty string to turn comments off.");
    parser.add_argument("-m", "--missingValue", default="nan",
                        help="Value used to indicate a missing value in the input text file. Default is 'nan'.");
    parser.add_argument("--encoding", default="utf-8",
                        help="Encoding of the input text file. Default it 'utf-8'");
    parser.add_argument("-d", "--dateIndex", type=int, default=0,
        help="The column number (starting from 0) which contains the date/time field in the in situ input data file. Default is a 0.");
    parser.add_argument("-n", "--numCommentLines", nargs="+", type=int, default=[],
        help="Number of comment lines before the header in your in situ input data file. These lines will be ignored). Default is a 0.");
    clArgs = parser.parse_args();
    
    return clArgs


if __name__ == "__main__":
    #Parse commandline arguments
    clArgs = parse_cl_arguments();
    feOutputPath = clArgs.feOutput;
    insituDataPath = clArgs.insituData;
    outputPath = clArgs.outPath;
    varsToAppend = clArgs.varsToAppend;
    delim = clArgs.delim;
    dateIndex = clArgs.dateIndex;
    rowsToSkip = clArgs.numCommentLines;
    missingValue = clArgs.missingValue;
    encoding = clArgs.encoding;
    latCol = clArgs.latCol;
    lonCol = clArgs.lonCol;
    
    #append selected variables as columns to the in situ data file
    append_to_in_situ(feOutputPath, insituDataPath, outputPath, varsToAppend, delim, latCol, lonCol, dateIndex, rowsToSkip, missingValue, encoding, verbose=True);
    
    print("Finished merging", varsToAppend, "from", feOutputPath, "with", insituDataPath);
    print("Output written to", outputPath);












