#! /usr/bin/env python
# -*- coding: utf-8 -*-
#Created by Peter Land, Plymouth Marine Laboratory.
#IGA branch to create files with different geographical grids (i.e. not necessarily regular)
#TMH option for specifying different grid resolutions

'''Program to convert text files to netcdf files with the correct dimensions to be used as input for the Flux Engine software. 
Input inFiles should be csv, or tab delimited text file(s). Headers are interrogated for lat/long, fCO2, N2O, CH4, SST and dates.
IGA ??/??/2016 Updated variable names to include CH4 and N20
IGA 5/5/2016 Updated filenames to allow for Arctic (or other) grid
IGA 11/5/2016 Updated to include vCO2_air data as req'd in FE
TMH 2018/05/17: Updated to increase flexibility and reflext new functionality of FluxEngine v3.0
TMH 2018/09/04: Updated to support use of temporal dimension in output files'''


import numpy as np, argparse;
from netCDF4 import Dataset;
from datetime import datetime, timedelta;
import pandas as pd;
from glob import glob;
#import pathlib;


#############
# Function definitions
#############

#Parse and set variables from command line arguments. Returns an object containing each command line argument as an attribute.
def parse_cl_arguments():
    description = unicode("""Converts text encoded data (e.g. csv, tsv) to netCDF3 data which is compatible for use with FluxEngine.
    The text input file must have a header containing the column names.
    Column names must start with a letter (a-z or A-Z) and can only contain letters, numbers, spaces, underscores and certain symbols.
    For information on allowed symbols consult the netCDF3 documentation.
    Symbol requirements do not apply to units (e.g. "Temp [°C]" or "fCO2 [µatm]" are valid column name if unit parsing is on.)""", 'utf-8');
    
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


#Try to parse a string to create a datatime object using a list of formats.
#Returns the datetime object and format used.
#If none of the formats match a ValueError exception is raised.
def parseDateTimeString(dateString, formats, verbose=False):
    dateString = dateString.strip();
    
    #formats = ["%d/%m/%Y", "%d/%m/%Y %H:%M", "%d/%m/%Y %H:%M:%S", "%d/%m/%Y%H:%M:%S", "%d/%m/%y", "%d/%m/%y %H:%M", "%d/%m/%y %H:%M:%S", "%Y-%m-%dT%H:%M:%S"]
    #Try to parse the string with each format in turn.
    for currentFormat in formats:
        try:
            return datetime.strptime(dateString, currentFormat), currentFormat;
        except:
            pass; #Try the next.
            
    #No format matched, so return None tuple.
    raise ValueError("datetimeFormat: No matching format for dateString '"+dateString+"'.");


# Returns a tuple containing the (x, y) grid coordinate that corresponds to a given latitude and longitude
#   latResolution: spatial resolution in latitude
#   lonResolution: spatial resolution in longitude
#   refLat: reference latitude at 0, 0 grid point
#   refLon: reference longitude at 0, 0 grid point
def get_grid_coordinates(latitude, longitude, latResolution, lonResolution, ref0Lat, ref0Lon):
    #Difference in lat/long
    #print "lat, lon:", latitude, longitude;
    #print "ref, ref:", ref0Lat, ref0Lon;
    dLatitude = latitude-ref0Lat;
    #print "dLat: ", dLatitude;
    dLongitude = longitude-ref0Lon;
    #print "dLon: ", dLongitude;
    if dLatitude < 0 or dLongitude < 0:
        raise ValueError("difference between latitude or longitude and reference point was negative.");
    
    #bin coordinates to the correct grid point
    x = int(dLatitude/latResolution);
    y = int(dLongitude/lonResolution);
    
    return (x, y);

#Calculates and returns the temporal coordinate (first dimension) diven a time, temporal resolution and reference data.
#   dateTime: date and time to calculate for
#   temporalResolution: temporal resolution
#   refDateTime: reference date time which corresponds to an index/coordinate of 0.
def get_temporal_coordinate(dateTime, temporalResolution, refDateTime):
    if temporalResolution != "monthly":
        diff = dateTime-refDateTime;
        index = int(np.floor(diff.total_seconds() / temporalResolution.total_seconds()));
        return index;
    else: #monthly
        yearDiff = dateTime.year - refDateTime.year;
        monthDiff = dateTime.month - refDateTime.month;
        index = (yearDiff*12) + monthDiff;
        return index;

#Calculates and sets standard deviations in outputArray. outputArray is written in place and no value is returned.
#   allValues: array containing a list of each value for each grid point
#   mask: boolean matrix indicating which grid points to calculate for
#   outputArray: Array to write output to
def calc_standard_deviation(allValues, mask, output):
    for coords, vals in np.ndenumerate(allValues):
        if mask[coords] == True:
            #print "sd: ", np.std(vals);
            output[coords] = np.nanstd(vals);

#Calculates and sets the mean in outputArray. outputArray is written in place and no value is returned.
#   allValues: array containing a list of each value for each grid point
#   mask: boolean matrix indicating which grid points to calculate for
#   outputArray: Array to write output to
def calc_mean(allValues, mask, output):
    for coords, vals in np.ndenumerate(allValues):
        if mask[coords] == True:
            output[coords] = np.nanmean(vals);

#Converts column indices to names, converts strings correct encoding. clArgs is modified in place.
#   dataFrame: data frame containing the full dataset
#   clArgs: parsed commandline arguments
def process_cols(dataFrame, columnNames, encoding):
    if columnNames == None: #If no columns are specified, then add them all.
        columnNames = dataFrame.keys();
    else: #Convert integers indexed to strings column names
        for i, colName in enumerate(columnNames):
            try:
                colIndex = int(colName);
                try:
                    columnNames[i] = dataFrame.keys()[colIndex];
                except IndexError:
                    print "WARNING: Invalid column number specified (%d) because there are not this many labels in the data's header. This will be ignored." % colIndex;
            except ValueError:
                pass;
    
    #Check each colName matches a column name in the header
    for colName in columnNames:
        if colName not in dataFrame.keys():
            print "WARNING: Unrecognised column name specified (%s). Column names must exactly match the header of your data file. Column names are case sensitive, and must be surrounded by quotes if they include spaces." %colName;
    
    #Finally convert them to the correct encoding
    for i, colName in enumerate(columnNames):
        if not isinstance(colName, unicode):
            columnNames[i] = unicode(colName, encoding);

#Create arrays to store and process data in.
#   colNames: list of column names to use (mean, count, stddev and value matrices will be created for each column)
#   temporalDimLength: length of the time dimension
#   gridDimLengthX: length of the grid (latitude)
#   gridDimLengthY: length of the grid (longitude)
def initialise_data_storage(colNames, temporalDimLength, gridDimLengthX, gridDimLengthY):
    allArrays = {}; #Store numpy arrays separately until we copy them over to the netCDF variable.
    print "The following columns will be extracted:";
    for colName in colNames:
        print "\t", colName;
        
        #somewhere to store the mean, count and standard deviation of each grid cell
        allArrays[colName+"_mean"] = np.full([temporalDimLength, gridDimLengthX, gridDimLengthY], 0.0, dtype='f');
        allArrays[colName+"_count"] = np.zeros([temporalDimLength, gridDimLengthX, gridDimLengthY], dtype='i');
        allArrays[colName+"_stddev"] = np.full([temporalDimLength, gridDimLengthX, gridDimLengthY], 0.0, dtype='f');
        
        #3D matrix of empty lists to store each value in.
        allArrays[colName+"_vals"] = np.empty([temporalDimLength, gridDimLengthX, gridDimLengthY], dtype=object);
        for coords, val in np.ndenumerate(allArrays[colName+"_vals"]):
            allArrays[colName+"_vals"][coords] = [];
    return allArrays;


#Creates the required netCDF dimensions, and adds latitude and longitude data
#   netCDFFile: the netCDF file to be written to
#   latData: latitude data for the grid
#   lonData: longitude data for the grid
def create_netCDF_dimensions(netCDFFile, latData, lonData, DIM_NAME_LAT, DIM_NAME_LON, DIM_NAME_TIME, startTime, timeDimLength, temporalResolution):
    if temporalResolution == "monthly":
        raise ValueError("'monthly' temporal resolution not yet supported in create_netCDF_dimensions.");
    
    #Create dimensions for the output netCDF file
    netCDFFile.createDimension(DIM_NAME_LAT, len(latData));
    netCDFFile.createDimension(DIM_NAME_LON, len(lonData));
    netCDFFile.createDimension(DIM_NAME_TIME, timeDimLength);
    
    #Create longitude and latitude variables
    lats = netCDFFile.createVariable('latitude','f',(DIM_NAME_LAT,));
    lats.units = 'degrees_north';
    lats.axis = "Y";
    lats.long_name = "Latitude North";
    lats.standard_name = "latitude";
    lats[:] = latData;
    lats.valid_min = -90.0;
    lats.valid_max = 90.0;
    
    lons = netCDFFile.createVariable('longitude','f',(DIM_NAME_LON,));
    lons.units = 'degrees_east';
    lons.axis = "X";
    lons.long_name = "Longitude East";
    lons.standard_name = "longitude";
    lons[:] = lonData;
    lons.valid_min = -180.0;
    lons.valid_max = 180.0;
    
    secs = netCDFFile.createVariable('time', 'f', ('time',));
    secs.units = 'seconds since 1970-01-01 00:00:00'
    secs.axis = "T"
    secs.long_name = "Time - seconds since 1970-01-01 00:00:00"
    secs.standard_name = "time"
    
    vals = [];
    curPoint = (startTime - datetime(1970, 1, 1)).total_seconds();
    for i in range(0, timeDimLength):
        vals.append(curPoint);
        curPoint += temporalResolution.total_seconds();
    secs[:] = vals;
        
    
    #firstPoint = (startTime - datetime(1970, 1, 1));
    #lastPoint = firstPoint + (timeDimLength * temporalResolution);
    #secs[:] = np.linspace(firstPoint.total_seconds(), lastPoint.total_seconds(), timeDimLength);
    secs.valid_min = 0.0 
    secs.valid_max = 1.79769313486232e+308
    


#Creates the required netCDF variables (mean, count, stddev) for each selected column and returns them
#   netCDFFile: the netCDF file to be written to
#   colNames: a list of column names (must match input data header)
#   parseUnits: if true, units will be split from the column name and set correctly (assumes format is "col name [units]")
def create_netCDF_variables(netCDFFile, colNames, parseUnits, DIM_NAMES, MISSING_VALUE):
    netCDFVariables = {};
    for colName in colNames:
        if parseUnits == True:
            unitStr = colName[colName.find("[")+1 : colName.find("]")];
            nameStr = colName.split("[")[0].strip();
        else:
            unitStr = "";
            nameStr = colName;
        
        newVar = netCDFFile.createVariable(nameStr+"_mean", 'f', DIM_NAMES);
        newVar.units = unitStr;
        newVar.long_name = colName;
        newVar.standard_name = colName;
        newVar.fill_value = MISSING_VALUE;
        netCDFVariables[colName+"_mean"] = newVar;
        #newVar.valid_min = -100;
        #newVar.valid_max = 1000;
        
        newVar = netCDFFile.createVariable(nameStr+"_count", 'i', DIM_NAMES);
        newVar.units = "count";
        newVar.long_name = colName;
        newVar.standard_name = colName;
        newVar.fill_value = MISSING_VALUE;
        netCDFVariables[colName+"_count"] = newVar;
        
        newVar = netCDFFile.createVariable(nameStr+"_stddev", 'f', DIM_NAMES);
        newVar.units = unitStr;
        newVar.long_name = colName;
        newVar.standard_name = colName;
        newVar.fill_value = MISSING_VALUE;
        netCDFVariables[colName+"_stddev"] = newVar;
    return netCDFVariables;

#Returns filePath with a timestamp appended to then end of the filename.
#The timestamp will correspond to the start of the given timestep.
#   filePath: path to append to
#   startTime: start of the first timestep
#   temporalResolution: length of each timestep (must be either deltatime object or the string "monthly")
#   temporalIndex: current timestep
def append_timestamp_to_filename(filePath, startTime, temporalResolution, temporalIndex):
    if filePath.endswith(".nc") == False:
        filePath += ".nc";
    
    if temporalResolution != "monthly":
        timestamp = startTime + (temporalResolution*temporalIndex);
        timestamp = datetime.strftime(timestamp, "%Y-%m-%d_%H-%M-%S");
    else: #monthly temporal resolution
        months = startTime.month-1 + temporalIndex;
        newYear = startTime.year + (months // 12);
        newMonth = (months % 12) + 1;
        timestamp = datetime.strftime(datetime(newYear, newMonth, 1), "%Y-%m-%d");
    
    parts = filePath.rsplit('.', 1);
    return parts[0]+"_"+timestamp+"."+parts[1];


def convert_text_to_netcdf(inFiles, startTime, endTime, ncOutPath,
                           limits=[-90.0, 90.0, -180.0, 180], latResolution=1.0, lonResolution=1.0,
                           temporalResolution="monthly", temporalChunking=1,
                           delim="\t", numCommentLines=0, encoding='utf-8', parseUnits=True, textFileMissingValue=np.nan,
                           dateIndex=0, colNames=None, latProd="Latitude", lonProd="Longitude", dateFormatDayFirst=False):
    ##############
    # Constant definitions
    ##############
    MISSING_VALUE = np.nan;#-999.0 #Missing value to be used in output netCDF file
    CL_DATE_FORMATS = ["%Y", "%Y-%m-%d", "%Y-%m-%d %H:%M", "%Y-%m-%d %H:%M:%S"]; #Formats which can be used to specify start and end dates in the commandline arguments
    
    #Dimension names for the output netCDF will not change, so define them here. Use to create variables in the output file.
    DIM_NAME_TIME = "time";
    DIM_NAME_LAT = "latitude";
    DIM_NAME_LON = "longitude";
    DIM_NAMES = (DIM_NAME_TIME, DIM_NAME_LAT, DIM_NAME_LON);
    
    ############
    # Process arguments to ensure correct formatting etc.
    ############
    
    #Convert YYYY start time to be last day of the year.
    try:
        intYear = int(startTime);
        startTime = str(intYear)+"-12-31 23:59:59";
    except ValueError:
        pass; #It's not a YYYY year format, so no need to convert.
        
    #Parse start and end times and create datetime objects
    try:
        startTime, formatUsed = parseDateTimeString(startTime, CL_DATE_FORMATS)
        startTime = startTime;
    except ValueError as e:
        print e.clArgs;
        raise SystemExit("Unable to parse start datetime "+startTime);
    try:
        endTime, formatUsed = parseDateTimeString(endTime, CL_DATE_FORMATS);
        endTime = endTime;
    except ValueError as e:
        print e.args;
        raise SystemExit("Unable to parse end datetime "+endTime);
    
    #Parse temporalResolution
    if temporalResolution != "monthly":
        days, time = temporalResolution.split(" ");
        hours, minutes = time.split(":");
        temporalResolution = timedelta(days=int(days), hours=int(hours), minutes=int(minutes));
    
    #Expand file globs
    expandedGlobs = [];
    for filepath in inFiles:
        expandedGlobs += glob(filepath);
    if len(expandedGlobs) != 0:
        inFiles = expandedGlobs;
#    expandedGlobs = [];
#    for inFile in inFiles:
#        expandedGlobs += pathlib.Path(".").glob(inFile);
#    expandedGlobs = [str(p) for p in expandedGlobs]; #Convert from pathlib.Path to string
    
    #numCommentLines should be a list
    if isinstance(numCommentLines, list) == False:
        numCommentLines = [numCommentLines];
    
    ###########
    # Define dimension data (lat, lon and time), grid size and resolution
    ###########
    print "Calculating dimensions."
    latitudeData = np.arange(limits[0]+(latResolution/2.0), limits[1]+(latResolution/2.0), latResolution); #+latResolution so it uses an inclusive range
    longitudeData = np.arange(limits[2]+(lonResolution/2.0), limits[3]+(lonResolution/2.0), lonResolution); #+lonResolution so it uses an inclusive range
    gridDimLengthX = len(latitudeData);
    gridDimLengthY = len(longitudeData);
    southLimit = limits[0]; #min lat
    northLimit = limits[1]; #max lat
    westLimit = limits[2]; #min lon
    eastLimit = limits[3]; #max lon
    temporalDimLength = get_temporal_coordinate(endTime, temporalResolution, startTime)+1; #Total number of time points (across all output files)
    
    ###########
    # Read and process data
    ###########
    #Read and process each input file
    allArrays = {};
    for iFile, inFile in enumerate(inFiles):
        #Parse data from text file into a pandas dataframe
        if len(numCommentLines) == 1: #Use the same numCommentLines value for all files
            rowsToSkip = numCommentLines[0];
        else: #Should have a different numCommentLines value for each file
            rowsToSkip = numCommentLines[iFile];
        df = pd.read_table(inFile, sep=delim, skiprows=rowsToSkip, parse_dates=[dateIndex], encoding=encoding, dayfirst=dateFormatDayFirst);
        
        
        #If this is the first input file to be read there are some additional things to setup
        if iFile==0:
            process_cols(df, colNames, encoding); #Process selected columns to ensure that they are valid column names/indices.
            allArrays = initialise_data_storage(colNames, temporalDimLength, gridDimLengthX, gridDimLengthY); #Create arrays to accumulate and process data in.
        
        #Loop through each row in the current file and process data
        print "Processing data in file", inFile;
        missingValueRows = []; #Missing values
        skippedRows = []; #Skipped due to being outside lon/lat limits
        for i, row in df.iterrows():
            #Check that the latitude and longitude are within limits
            if (row[latProd] > northLimit or row[latProd] < southLimit) or (row[lonProd] > eastLimit or row[lonProd] < westLimit):
                skippedRows.append( (inFile, i) );
                continue;
            #Check temporal limits
            if (row[dateIndex] < startTime or row[dateIndex] > endTime):
                skippedRows.append( (inFile, i) );
                continue;
            
            #Calculate coordinates
            xCoord, yCoord = get_grid_coordinates(row[latProd], row[lonProd], latResolution, lonResolution, southLimit, westLimit);
            temporalCoord = get_temporal_coordinate(row[dateIndex], temporalResolution, startTime);
            if temporalCoord < 0:
                raise IndexError("Trying to use time index of %d. Time index cannot be negative. Are you specifying the correct startTime?"%temporalCoord);
            
            #Process data for each variable from the current row
            for colName in colNames:
                curValue = row[colName];
                if str(curValue) == textFileMissingValue:
                    missingValueRows.append( (inFile, i) );
                    continue;
                
                allArrays[colName+"_count"][temporalCoord, xCoord, yCoord] += 1;
                allArrays[colName+"_vals"][temporalCoord, xCoord, yCoord].append(curValue);
    
    #Now that all the data files have been parsed, calculate the stuff we're interested in.
    #Final processing for each column: calculate mean, std. dev., and set missing values
    for colName in colNames:
        #Set missing values
        validMask = allArrays[colName+"_count"][:] != 0;
        allArrays[colName+"_mean"][validMask==False] = MISSING_VALUE;
        allArrays[colName+"_stddev"][validMask==False] = MISSING_VALUE;
        
        #Calculation standard deviation
        calc_standard_deviation(allArrays[colName+"_vals"], validMask, allArrays[colName+"_stddev"]);
        
        #Calculate mean
        calc_mean(allArrays[colName+"_vals"], validMask, allArrays[colName+"_mean"]);
    

    ##############
    # Create and write output netCDF files
    ##############
    #A seperate netCDF file will be created for each index in the time dimension (unless chunking is being used).
    #Loop through the temporal index and create the netCDF files
    print "Writing output netCDF file(s)...";
    currentChunk = 0;
    for temporalIndex in range(0, temporalDimLength):
        #Create netCDF file if required
        if currentChunk == 0:
            if temporalDimLength == 1 or temporalDimLength >= temporalChunking: #if only one file
                outFilePath = ncOutPath;
            else: #must append the start timestamp of the current temporal step to the filename
                outFilePath = append_timestamp_to_filename(ncOutPath, startTime, temporalResolution, temporalIndex);
            ncOutput = Dataset(outFilePath, 'w');
        
            #Create dimensions set lat/lon data
            create_netCDF_dimensions(ncOutput, latitudeData, longitudeData, DIM_NAME_LAT, DIM_NAME_LON, DIM_NAME_TIME, startTime, temporalChunking, temporalResolution);
        
            #Create the required netCDF variables (mean, count, stddev)
            allVariables = create_netCDF_variables(ncOutput, colNames, parseUnits, DIM_NAMES, MISSING_VALUE);
        
        #Copy relevant time slice of each array to the netCDF
        for colName in colNames:
            allVariables[colName+"_mean"][currentChunk, :, :] = allArrays[colName+"_mean"][temporalIndex,:,:];
            allVariables[colName+"_stddev"][currentChunk, :, :] = allArrays[colName+"_stddev"][temporalIndex,:,:];
            allVariables[colName+"_count"][currentChunk, :, :] = allArrays[colName+"_count"][temporalIndex,:,:];
        
        currentChunk = int((currentChunk+1)%temporalChunking);
        
        #Close netCDF file
        if currentChunk == 0:
            ncOutput.close();
    
    
    ############
    # Output summary
    ############
    
    #Output some info to the user
    print "Finished converting text file to netCDF3. There were %d values which fell outside the specified lat/lon or start/stop time boundaries." % len(skippedRows)
    #if len(invalidRowIndices) != 0:
    #    print "Rows with invalid values:", invalidRowIndices;
    print "Number of missing values found:", len(missingValueRows);
    #if len(missingValueRows) != 0:
    #    print "Rows with missing values:", missingValueRows;
    
    

if __name__ == "__main__":
    #parse command line arguments
    print "Parsing command line arguments.";
    args = parse_cl_arguments();
    
    convert_text_to_netcdf(args.inFiles, args.startTime, args.endTime, args.ncOutPath,
                           limits=args.limits, latResolution=args.latResolution, lonResolution=args.lonResolution,
                           temporalResolution=args.temporalResolution, temporalChunking=args.temporalChunking,
                           delim=args.delim, numCommentLines=args.numCommentLines, encoding=args.encoding,
                           textFileMissingValue=args.missing_value, parseUnits=args.parse_units,
                           dateIndex=args.dateIndex, colNames=args.cols, latProd=args.latProd, lonProd=args.lonProd, dateFormatDayFirst=args.dateFormatDayFirst); 

