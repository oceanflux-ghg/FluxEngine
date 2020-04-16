#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 12 10:39:52 2018

Contains functions which are used to setup and run FluxEngine simulations.
This includes:
    Parsing and validating config files
    Intrepretting config variables
    Copying run parameters for particular years/months

@author: tomholding
"""

from os import path, listdir, makedirs;
import inspect;
import time;
import socket; #for gethostname
import calendar;
from datetime import timedelta, datetime;
import fnmatch; #matching file globs
from numpy import cumsum;

from fluxengine.core import rate_parameterisation as k_params; #This is where k parameterisation logic is kept.
from fluxengine.core import data_preprocessing as data_preprocessing; #preprocessing functions
from fluxengine.core import process_indicator_layers as indicator_layers; #Process indicator layer functors
from fluxengine.core import fe_core as fluxengine;


#Gets the root fluxengine directory
def get_fluxengine_root():
    selfPath = inspect.stack()[0][1];
    rootPath = path.dirname(path.dirname(selfPath)); #Assumes this file is in the root/core/ directory.
    return rootPath;

#Parses .conf files and returns a dictionary containing variable values.
#Makes no checks of the validity of these variables, simply splits based on '=' sign:
#<name> = <value> #comment
#Warns if duplicate definitions found but takes the first value.
def read_config_file(configPath, verbose=False):
    function = inspect.stack()[0][1]+", "+inspect.stack()[0][3];
    #Check file exists
    if path.isfile(configPath) == False:
        raise ValueError("%s: Could not find config file at path: %s" % (function, configPath));
    configVariables = {}
    try:
        for line in open(configPath):
            line = line.split('#', 1)[0]; #Ignores comments
            if "=" in line: #The format we're looking for is "<name> = <value>" format
                try: #Extracting they name value pair
                    name, value = line.split('=', 1);
                    name = name.strip();
                    value = value.strip();
                    if (name in configVariables) == False: #Duplicate definition of same variable?
                        configVariables[name] = value;
                    else:
                        print("Warning: Duplicate definition in config file for '%s'. The first definition will be used." % name);
                except: #Will only happen if there isn't both left and righthand sides to the '=' assignment
                    print("%s: Error parsing config file line:\n"%function, line);
                    return None;
    except Exception as e:
        print("Error while parsing config file at path:", configPath);
        print(type(e), e.args);
        return None;
    
    if verbose:
        print("\nParsed configurable variables as follows:");
        for key in list(configVariables.keys()):
            print(key, configVariables[key]);
    return configVariables;


#Given a string, parse it for configuration file version information
def parse_config_version_tag(firstline):
    if "#?FluxEngineConfigVersion:" in firstline:
        return float(firstline.split(":")[1]);
    else:
        raise ValueError("No version tag found in configuration file. Instead found:\n\t"+firstline);

#Load and return metadata about the variables contained in config files.
#Returns a dictionary of dictionarys allowing indexing by variablename then attribute, e.g.:
#   varMetadata[varName][attribute]
def read_config_metadata(settingsPath, verbose=False):
    function = inspect.stack()[0][1]+", "+inspect.stack()[0][3];
    
    import xml.etree.ElementTree as ET; #Only needed here.
    
    if verbose:
        print("Parsing settings file at:", settingsPath);
    tree = ET.parse(settingsPath);
    root = tree.getroot();
    
    #Read config metadata for each variable
    varMetadata = {};
    configParametersElement = root.find("ConfigParameters");
    if configParametersElement != None:
        for element in configParametersElement:
            try:
                if element.tag == "Variable":
                    #Multioption: Multioptions must also contain a list of options
                    if element.attrib["type"] == "multioption":
                        varMetadata[element.attrib["name"]] = element.attrib;
                        options = {};
                        for child in element:
                            if child.tag == "Option":
                                try:
                                    options[child.attrib["str"]] = int(child.attrib["value"]);
                                except ValueError as e:
                                    print("%s: Invalid multioption child 'Option' element in %s. 'str' and 'value' attributes must be specified." % (function, element.attrib["name"]))
                                    print(e.args);
                        varMetadata[element.attrib["name"]]["options"] = options;
                    
                    #DataLayer: Datalayers require several variables to define them, so create these here.
                    if element.attrib["type"] == "DataLayer":
                        varMetadata[element.attrib["name"]+"_path"] = {"name": element.attrib["name"],
                                                                       "required": element.attrib["required"],
                                                                       "type": "DataLayerPath"};
                        varMetadata[element.attrib["name"]+"_prod"] = {"name": element.attrib["name"],
                                                                       "required": element.attrib["required"],
                                                                       "type": "string"};
                        varMetadata[element.attrib["name"]+"_netCDFName"] = {"name": element.attrib["name"],
                                                                       "required": "false",
                                                                       "type": "string"}
                        varMetadata[element.attrib["name"]+"_units"] = {"name": element.attrib["name"],
                                                                       "required": "false",
                                                                       "type": "string"}
                        varMetadata[element.attrib["name"]+"_minBound"] = {"name": element.attrib["name"],
                                                                       "required": "false",
                                                                       "type": "float"}
                        varMetadata[element.attrib["name"]+"_maxBound"] = {"name": "false",
                                                                       "required": element.attrib["required"],
                                                                       "type": "float"}
                        varMetadata[element.attrib["name"]+"_standardName"] = {"name": element.attrib["name"],
                                                                       "required": "false",
                                                                       "type": "string"}
                        varMetadata[element.attrib["name"]+"_longName"] = {"name": element.attrib["name"],
                                                                       "required": "false",
                                                                       "type": "string"}
                        
                        varMetadata[element.attrib["name"]+"_temporalChunking"] = {"name": element.attrib["name"],
                                                                       "required": "false",
                                                                       "type": "integer"}
                        varMetadata[element.attrib["name"]+"_temporalSkipInterval"] = {"name": element.attrib["name"],
                                                                       "required": "false",
                                                                       "type": "integer"}
                        varMetadata[element.attrib["name"]+"_timeDimensionName"] = {"name": element.attrib["name"],
                                                                       "required": "false",
                                                                       "type": "string"}
                    else: #just a straightforward variable, so add it.
                        varMetadata[element.attrib["name"]] = element.attrib;
            except ValueError as e:
                print("%s: 'name' or 'type' attribute missing from Variable element in settings xml file. Specified settings file is invalid: %s" % (function, settingsPath));
                print(e.args);
                
    else:
        raise IOError("%s: Couldn't find parent 'ConfigParameters' element in %s."%(function, settingsPath));
    
    return varMetadata;


#Verifies config variables contain valid values. Reports missing variables or invalid values.
#Converts variables to appropriate types expected by FluxEngine.
#Modified configVariables in place.
    #configVariables is a dictionary detailing name:value for each variable defined in the config file
    #metadata is a dictionary or dictionaries containing the metadata for each variable
def verify_config_variables(configVariables, metadata, verbose=False):
    function = inspect.stack()[0][1]+", "+inspect.stack()[0][3];
    
    #Add any parameters with default values if they're not already specified
    for varName in metadata:
        if "default" in metadata[varName] and varName not in configVariables:
            configVariables[varName] = metadata[varName]["default"];
    
    #Identify custom vs standard variables and make seperate lists of them.
    customVars = []; #Custom variables are used by extension components (e.g. k parameterisations).
    standardVars = [];
    for varName in configVariables:
        if varName in metadata:
            standardVars.append(varName);
        else:
            customVars.append(varName);
    
    #Process each standard variable to check for validity and consistency.
    #Convert to the data types fluxengine expects while we're at it.
    for varName in standardVars:
        #Different data types require processing differently and have different constraints.
        #Paths: Convert all paths to absolute paths.
#        if metadata[varName]["type"] == "path" or metadata[varName]["type"] == "DataLayerPath":
#            if path.isabs(configVariables[varName]) == False:
#                configVariables[varName] = path.abspath(path.expanduser(configVariables[varName]));
        
        #Multioption: Check that the config value matches at least one option.
        #             Store string value with _name suffix overwrite original with integer value.
        if metadata[varName]["type"] == "multioption":
            if configVariables[varName] not in metadata[varName]["options"]:
                raise ValueError("%s: Invalid option specified in config file for %s. Valid options are: %s" % (function, varName, str(list(metadata[varName]["options"].keys()))));
            configVariables[varName+"_name"] = configVariables[varName]; #name is useful for writing messages to user
            configVariables[varName] = metadata[varName]["options"][configVariables[varName]]; #Overwrite with integer option
        
        #Switch: convert 'yes' and 'no' strings to True and False.
        elif metadata[varName]["type"] == "switch":
            #Convert 'no' to False and 'yes' to True.
            if configVariables[varName].lower() == "yes":
                configVariables[varName] = True;
            elif configVariables[varName].lower() == "no":
                configVariables[varName] = False;
            else:
                raise ValueError("%s: Invalid value for %s in config file. Accepted values are 'yes' or 'no' but found '%s'." % (function, varName, configVariables[varName]));

        #Float: string to float.
        elif metadata[varName]["type"] == "float":
            try:
                configVariables[varName] = float(configVariables[varName]);
            except ValueError:
                print("%s: Config variable '%s' requires a decimal / floating point number. Got %s instead." % (function, varName, configVariables[varName]));
        
        #Int: string to int
        elif metadata[varName]["type"] == "integer":
            try:
                configVariables[varName] = float(configVariables[varName]);
            except ValueError:
                print("%s: Config variable '%s' requires a whole / integer number. Got %s instead." % (function, varName, configVariables[varName]));
        
        #Strings: Do nothing, strings are just strings.
        elif metadata[varName]["type"] == "string":
            pass;

    #Check that each datalayer path has at least a corresponding prod.
    for varName in standardVars:
        #Datalayers must have at least a datalayername_path and a datalayername_prod defined.
        #If we get any DataLayerPath types then check for the corresponding prod.
        if metadata[varName]["type"] == "DataLayerPath":
            if varName[0:-4]+"prod" not in configVariables:
                raise ValueError("%s: DataLayer path is specified (%s) but there is no corresponding prod. '%s' must be also be defined in the configuration file." % (function, varName, varName[0:-4]+"prod"));
    
    #TODO: add datetime data type to settings.xml?
    #deltatime is a special case
    #Parse delta time config attribute
    try:
        days, time = configVariables["temporal_resolution"].split(" ");
        hours, minutes = time.split(":");
        configVariables["temporal_resolution"] = timedelta(days=int(days), hours=int(hours), minutes=int(minutes));
        
        #tmpTime = datetime.strptime(configVariables["temporal_resolution"], "%d %H:%M");
        #configVariables["temporal_resolution"] = timedelta(days=tmpTime.day, hours=tmpTime.hour, minutes=tmpTime.minute);
        if verbose:
            print("Temporal resolution set to: ", configVariables["temporal_resolution"]);
    except ValueError:
        if verbose:
            print("Temporal resolution has been set to default (monthly)");
        configVariables["temporal_resolution"] = None;
    
    #Now process custom vars. Try to convert them to a float, but if they fail assume they're supposed to be a string.
    for varName in customVars:
        try:
            configVariables[varName] = float(configVariables[varName]);
        except ValueError:
            pass; #Leave it as a string.


#Substitutes various time tokens into the input string.
#Does not modify the original arguments.
#Accounts for leap years correctly.
#curDatetime should be a datetime object.
def substitute_tokens(inputStr, curDatetime):
    year = curDatetime.year;
    month = curDatetime.month;
    
    day = curDatetime.day;
    hour = curDatetime.hour;
    minute = curDatetime.minute;

    outputStr = inputStr;
    outputStr = outputStr.replace("<YYYY>", str(year));
    outputStr = outputStr.replace("<YY>", str(year)[-2:]);
    outputStr = outputStr.replace("<MM>", "%02d"%(month));
    outputStr = outputStr.replace("<M>", str(month)); #numeric month with no 0 padding
    outputStr = outputStr.replace("<MMM>", calendar.month_abbr[month].upper());
    outputStr = outputStr.replace("<Mmm>", calendar.month_abbr[month]);
    outputStr = outputStr.replace("<mmm>", calendar.month_abbr[month].lower());
    outputStr = outputStr.replace("<DD>", "%02d"%day); #<DD> day of the month
        
    #<DDD> day of the year
    cumulativeDaysByMonth = [0] + list(cumsum([calendar.monthrange(year, m)[1] for m in range(1,12)])); #Cumulative days in each month of the year (accounting for leap years)
    outputStr = outputStr.replace("<DDD>", "%03d"%(day+cumulativeDaysByMonth[month-1])); #Day number in year, starting at 1 (up to 365 or 366 on leap years)
    
    outputStr = outputStr.replace("<hh>", "%02d"%hour); #<hh> hours in 24 hour time.
    outputStr = outputStr.replace("<mm>", "%02d"%minute) #<mm> mintue past the hour
    
    #FluxEngine root
    outputStr = outputStr.replace("<FEROOT>", get_fluxengine_root());

    return outputStr;


#Returns a list of filepaths which match a glob
def match_filenames(givenPath, glob):
    #Make sure the glob doesn't include any folders here. If it does add them to the path.
    #(Otherwise fnmatch crashes the kernal for some reason...)
    splits = glob.split("/");
    if len(splits) > 1:
        givenPath = path.join(givenPath, *splits[:-1]);
        glob = splits[-1];
    matches = [];
    for filename in listdir(givenPath): #Search for a match
        if fnmatch.fnmatch(filename, glob):
            matches.append(path.join(givenPath, filename));    
    return matches;


#Creates a dictionary of parameters needed to run the flux engine for a particular month and year.
#Generates paths to input files using supplied globs
#Converts between inconsistencies in the config file naming convention and FluxEngine naming convention, e.g.:
#   path -> infile
#   appending _switch to bool variables
#Arguments:
#   configVariables:    Parsed .conf file
#   varMetadata:        Metadata for configuration files
#   curTimePoint:       Current datetime of the time point to run
#   executionCount:     Number of previous executions in this set (used to group output files temporally).
def create_run_parameters(configVariables, varMetadata, curTimePoint, executionCount, processTimeStr, configFile, processIndicatorLayersOff, previousRunParams):
    function = inspect.stack()[0][1]+", "+inspect.stack()[0][3];
    
    runParams = {};
    #Copy over misc. parameters
    #runParams["year"] = 2000 if clArgs.use_takahashi_validation==True else runParams["year"] = year;
    runParams["year"] = curTimePoint.year;
    runParams["month"] = curTimePoint.month;
    runParams["day"] = curTimePoint.day;
    runParams["hour"] = curTimePoint.hour;
    runParams["minute"] = curTimePoint.minute;
    runParams["second"] = curTimePoint.second;
    runParams["run_count"] = executionCount;
    runParams["hostname"] = socket.gethostname();
    runParams["config_file"] = configFile;
    runParams["processing_time"] = processTimeStr;
    runParams["monthyear_processing_time"] = time.strftime("%d/%m/%Y %H:%M:%S");
    runParams["process_layers_off"] = int(processIndicatorLayersOff);
    
    #infer whether each sst
    
    #Copy over all config parameters:
    missingFiles = []; #Much more useful to the user to have a complete list of missing files instead of one at a time.
    for varName in configVariables:
        #Different variable types need to be copied differently
        #DataLayerPaths: These can change based on month and year.
        if varName in varMetadata:
            if varMetadata[varName]["type"] == "DataLayerPath":
                pathGlob = substitute_tokens(configVariables[varName], curTimePoint); #substitute time tokens into the path/glob
                pathGlob = path.abspath(path.expanduser(pathGlob)); #Expand to absolute path only AFTER substituting tokens
                curDir = path.dirname(pathGlob); #Directory to search in.
                curGlob = path.basename(pathGlob); #File glob to match against.
                #print "path:", curPath
                #print "glob:", curGlob;
                try:
                    matches = match_filenames(curDir, curGlob);
                    matches = [match for match in matches if match.rsplit('/', 1)[-1][0] != '.']; #Remove hidden metadata files created by MACOS on some file systems. Required if sharing folders between MACOS and LINUX/WINDOWS.
                    
                    #If there's exactly one match for the file glob, add the data layer.
                    if len(matches) == 1: #Store the file path.
                        dataLayerName = varName[0:-5];
                        runParams[dataLayerName+"_infile"] = path.join(curDir, matches[0]);
                        runParams[dataLayerName+"_prod"] = configVariables[dataLayerName+"_prod"];
                        runParams[dataLayerName+"_stddev_prod"] = configVariables[dataLayerName+"_stddev_prod"] if dataLayerName+"_stddev_prod" in configVariables else None;
                        runParams[dataLayerName+"_count_prod"] = configVariables[dataLayerName+"_count_prod"] if dataLayerName+"_count_prod" in configVariables else None;
                        
                        if dataLayerName+"_preprocessing" in configVariables:
                            runParams[dataLayerName+"_preprocessing"] = configVariables[dataLayerName+"_preprocessing"];
                        else:
                            runParams[dataLayerName+"_preprocessing"] = "no_preprocessing";
                        
                    else: #There isn't exactly 1 matching file, so we have a problem, add it to the list.
                        missingFiles.append( (varName, path.join(curDir, curGlob) ) );
                except Exception as e:
                    print("Invalid filepath when checking for %s: %s" % (varName, curDir));
                    print(type(e), e.args);
                    return None;
        
            #Switches: _switch suffix added for clarity
            if varMetadata[varName]["type"] == "switch":
                runParams[varName+"_switch"] = configVariables[varName];
        
            #Every other type is simple - just copy over... (paths, strings, floats, multioption (value and name)...)
            else:
                runParams[varName] = configVariables[varName];
        
        else: #varName has no metadata associated with it so just copy as is
            runParams[varName] = configVariables[varName];
        
        #if there were missing files, raise them.
        if len(missingFiles) != 0:
            print("Missing files:");
            for missingFile in missingFiles:
                print(missingFile[0]+": "+missingFile[1]);
            raise ValueError("%s: Error finding unique file matches for one or more input data layers." % (function));
            return None;
            

    #Output file/path
    outputChunk = executionCount%runParams["output_temporal_chunking"];
    #runParams["output_chunk"] = int(outputChunk); #if 0 a new file will be created, otherwise it will be append to previous chunk
    
    if outputChunk == 0: #Will need to create a new file
        outputRoot = substitute_tokens(configVariables["output_dir"], curTimePoint);
        outputStructure = substitute_tokens(configVariables["output_structure"], curTimePoint);
        outputFile = substitute_tokens(configVariables["output_file"], curTimePoint);
    
        runParams["output_dir"] = path.join(outputRoot, outputStructure); #root output directory
        runParams["output_path"] = path.join(outputRoot, outputStructure, outputFile);
    else: #Output will be appended to previous file
        runParams["output_dir"] = previousRunParams["output_dir"];
        runParams["output_path"] = previousRunParams["output_path"];
    
    return runParams;
    

def list_input_datalayers(runParameters, metadata):
    dataLayers = [];
    for varName in runParameters:
        if varName in metadata and metadata[varName]["type"] == "DataLayerPath":
            dataLayers.append(varName[0:-5]); #Truncate the '_path' suffix
    return dataLayers;
   

#returns an object of the appropriate class, as specified by k_parameterisation in the config file.
#Classes are checked so that they derive from KCalculationBase.
#If no suitable class is found, or there is some other error None is returned.
def build_k_functor(runParameters):
    function = inspect.stack()[0][1]+", "+inspect.stack()[0][3];
    try:
        selectedParameterisation = runParameters["k_parameterisation"];
        
        #Search for the specified k parameterisation class
        for name, obj in inspect.getmembers(k_params):
            if (name == selectedParameterisation) and (inspect.isclass(obj) == True):
                #Get a handle to the class type
                try:
                    ClassHandle = getattr(k_params, name);
                except AttributeError as e:
                    print("%s: Could not find a k parameterisation functor which corresponds to the specified k_parameterisation (%s). If this specified correctly in the config file?"%(function, name));
                    print(e.args);
                    return None;
                
                #Check it derives from KCalculateBase - this is used as a simple way to filter out irrelevant classes
                if issubclass(ClassHandle, k_params.KCalculationBase):
                    #To initialise the class we need to know the arguments it's __init__ function uses.
                    initialiserArgNames = [arg for arg in inspect.getargspec(ClassHandle.__init__).args if arg != "self"];
                    
                    try:
                        #Create a dictionary of name:value pairs for the __init__ arguments. Assumes arguments match config file names.
                        argDict = {key : runParameters[key] for key in initialiserArgNames if key in runParameters}; #create a dictionary of arguments name:value pairs
                    except KeyError as e:
                        print("%s: Could not find all the required initialiser arguments. Are they specified correctly in the config file?\nExpected arguments are: "%function, initialiserArgNames);
                        print("KeyError.args: ", e.args);
                        return None;
                    
                    #Finally create and return the k functor instance
                    return ClassHandle(**argDict);
    
    except ImportError as e:
        print("%s: Cannot find k_parameterisation module. Check FluxEngine is installed correctly.", e.args);
        return None;
    return None;

#Returns a list of preprocessing functions (matched by name from fluxengine_src.data_preprocessing.py)
#funcNamesArg is a string of function names delimited by commas (i.e. parsed from the config file)
def get_preprocessing_funcs(funcNamesArg):
    function = inspect.stack()[0][1]+", "+inspect.stack()[0][3];
    functionList = [];

    #get a list of each funcName by splitting on commas
    funcNames = [s.strip() for s in funcNamesArg.split(',')];
    for funcName in funcNames:
        #Search for the specified function and return is
        for name, obj in inspect.getmembers(data_preprocessing):
            if name == funcName:
                functionList.append(obj);
                break;
        #if we get to here we found no matching function so raise an exception
        if name != funcName:
            raise ValueError("%s: Couldn't find preprocessing function matching '%s'." % (function, funcName));
    
    return functionList;


#Constructs a fluxengine object based on a set of run parameters.
#The resulting fe object is ready to run.
def fe_obj_from_run_parameters(runParameters, metadata, processLayersOff=True, verbose=False):
    fe = fluxengine.FluxEngine(runParameters);
            
    #Add the k parameterisation functor
    kFunctor = build_k_functor(runParameters);
    if kFunctor != None:
        fe.add_k_parameterisation_component(kFunctor);
    else:
        raise SystemExit("Could not create the k_parameterisation object. Does the config file specifiy a valid k_parameterisation (%s is used) and specify all the required initialisers?" % (runParameters["k_parameterisation_name"]));
    
    #TODO: these should be specifed by name in the config file so that they can be individually turned on/off and their required inputs checked etc.
    #Set process indicator layers
    if processLayersOff == False:
        fe.add_process_indicator_functor(indicator_layers.low_wind_indicator());
        fe.add_process_indicator_functor(indicator_layers.bioclass_indicator());
        fe.add_process_indicator_functor(indicator_layers.diurnal_warming_indicator());
        fe.add_process_indicator_functor(indicator_layers.oceanic_basins_indicator());
        fe.add_process_indicator_functor(indicator_layers.longhurst_provinces_indicator());
    
    
    #TODO: Do this in a similar way to the data preprocessing; a list separated by commas allowing sequential definition
    #      to be built up
    ###Add additional k parameterisation components:                      
    #linear additive k/krain for rain case - adding rain component to the results from the existing choice of parameterisation
    if (runParameters["k_rain_linear_ho1997_switch"] == 1):
        fe.add_k_parameterisation_component(k_params.AddKRainLinearHo1997());
    
    #nonlinear changes to k/krain for rain and wind - parameterisation from Harrison et al., JGR 2012, equations 11, 12, 13, 14
    if (runParameters["k_rain_nonlinear_h2012_switch"] == 1 and runParameters["k_parameterisation"] == 9): #9==generic_k
        fe.add_k_parameterisation_component(k_params.AddKRainNonlinearHarrison2012());
    
    
    ###Add each datalayer
    inputDataLayers = list_input_datalayers(runParameters, metadata);
    for dataLayerName in inputDataLayers:
        if dataLayerName+"_infile" in runParameters:
            preprocessingFuncs = get_preprocessing_funcs(runParameters[dataLayerName+"_preprocessing"]);
            status = fe.add_data_layer(dataLayerName, runParameters[dataLayerName+"_infile"],
                                       runParameters[dataLayerName+"_prod"],
                                       stddevProd=runParameters[dataLayerName+"_stddev_prod"], #Note: if no stddev or count data then these are set to None
                                       countProd=runParameters[dataLayerName+"_count_prod"],
                                       preprocessing=preprocessingFuncs);
            if status == False:
                break;
    
    if status == False:
        print("Error adding data layers.");
        fe = None;
    
    return fe;

#Takes a string start date and a string end date and returns a list of date objects between them (inclusive of start and end)
#deltaTime is a datetime.deltatime object which specifies the interval between timepoints. If set to None then intervals will be monthly.
#singleDate returns a list containing only the first date (used for testing).
#Valid date/time string formats: YYYY, YYYY-MM-DD, YYYY-MM-DD hh:mm
def generate_datetime_points(startStr, endStr, deltaTime=None, singleDate=False):
    #Parse start/end date strings
    try: #Parsing start date string
        startYear = int(startStr);
        startDate = datetime(startYear, 1, 1);
    except ValueError: #startStr is not a single year, so try it as a datetime string
        try:
            startDate = datetime.strptime(startStr, "%Y-%m-%d %H:%M");
        except ValueError: #startStr is not a single year, so try it as a date string (with no time component)
            startDate = datetime.strptime(startStr, "%Y-%m-%d");
    
    try: #Parsing end date string
        endYear = int(endStr);
        endDate = datetime(endYear, 12, 31);
    except ValueError: #endStr is not a single year, so try it as a datetime string
        try:
            endDate = datetime.strptime(endStr, "%Y-%m-%d %H:%M");
        except ValueError: #endStr is not a single year, so try it as a date string (with no time component)
            endDate = datetime.strptime(endStr, "%Y-%m-%d");
    
    
    if startDate > endDate:
        raise ValueError("Start date (%s) is after end data (%s)." % (str(startDate), str(endDate)));
    
    #Create datetime objects for each date between the start and end dates.
    timePoints = [];
    currentDate = startDate;
    while currentDate <= endDate:
        timePoints.append(currentDate);
        if deltaTime == None: #Increment by the number of days in the current month (taking into account leap years)
            currentDate += timedelta(days=calendar.monthrange(currentDate.year, currentDate.month)[1]);
        else: #Increment by deltatime
            currentDate += deltaTime;

    if singleDate: #Discard all other time points and just return the first one.
        return [timePoints[0]];
    else:
        return timePoints;


#Takes a config file, a list of years and months, and runs the flux engine for each month/year combination.
def run_fluxengine(configFilePath, startDate, endDate, singleRun=False, verbose=False, processLayersOff=True,
                   takahashiDriver=False, pco2DirOverride=None, outputDirOverride=None, dailyResolution=False):
    function = inspect.stack()[0][1]+", "+inspect.stack()[0][3];
    hostname = socket.gethostname();
    rootPath = path.abspath(path.expanduser(path.join(__file__, "../..")));
    if verbose:
        print("Hostname identified as: ", hostname);
        print("Working directory is: ", rootPath);
    
    #Parse config file
    configPath = path.join(configFilePath); #Don't make configFilePath absolute!
    configVariables = read_config_file(configPath, verbose=verbose);
    
    #Parse settings file for default metadata about the config variables
    settingsPath = path.join(rootPath, "core", "settings.xml");
    metadata = read_config_metadata(settingsPath, verbose=verbose);
    
    #Append custom datalayers to metadata file, this means they will be automatically added as a datalayer
    for varName in list(configVariables.keys()):
        if (varName not in metadata) and (varName[-5:] == "_path"):
            varBase = varName[:-5];
            metadata[varName] = {"required":"false", "type":"DataLayerPath", "name":varBase};
    
    #Substitute commandline override arguments (-pco2_dir_override, -output_dir_override)
    if (pco2DirOverride != None):
        print("Using optional override for pCO2w data directory. '%s' will be set to '%s' (ie overriding both directory of pco2 data and selection in the configuration file)." % (configVariables["pco2"], pco2DirOverride));
        configVariables["pco2_path"] = path.abspath(path.expanduser(pco2DirOverride));
    if (outputDirOverride != None):
        configVariables["output_dir"] = path.abspath(path.expanduser(outputDirOverride));
        print("Using optional override for output directory, output_dir will be set to %s (overriding the configuration file value)." % (outputDirOverride));
    
    #Printing some feedback...
    if processLayersOff == True and verbose:
        print("Switching off generation of processing indicator layers. This reduces the processing time by appx. 50% (switch: process_layers_off).");
        #configVariables["process_indicator_layers"] = None;
    if takahashiDriver == True and verbose:
        print("This is a takahashi validation run. Ensure config is the configs/takahashi09_validation.conf file supplied with FluxEngine.");
    
    #Using the metadata, process the config variables appropriately
    #Checks types are valid, converts strings to the data types required by fluxengine
    #Checks data layers contain at least a prod and path, etc.
    verify_config_variables(configVariables, metadata, verbose=verbose);
    
    ############
    ##TODO:
    #Alert user to any unused variables (can we infer all required variables yet?).    
    #Error if required variable hasn't been specified.
    ############
    
    
    #Begin main execution logic
    processTimeStr = time.strftime("%d/%m/%Y %H:%M:%S");
    print("Executing on '%s' at %s" % (hostname, processTimeStr));
    
    #Generate a list of datetime objects overwhich to run the simulations
    timePoints = generate_datetime_points(startDate, endDate, deltaTime=configVariables["temporal_resolution"], singleDate=singleRun);
    #print "num timepoints:", len(timePoints);
    #input("key to continue...");
    for i, timePoint in enumerate(timePoints):
        #Run parameters can vary depending on the current month and year (e.g. paths and filenames,
        #So these must be generated on a per-month/year basis.
        try:
            if i==0:
                runParameters = None; #Creates an empty previous runParameters... TODO: tidy this and the create_run_parameters function
            runParameters = create_run_parameters(configVariables, metadata, timePoint, i, processTimeStr, configFilePath, processLayersOff, runParameters);
            
            #TODO: Takahashi driver switch should be removed from future releases and moved to the configuration file.
            if takahashiDriver == True:
                runParameters["TAKAHASHI_DRIVER"] = True;
            else:
                runParameters["TAKAHASHI_DRIVER"] = False;
        except ValueError as e:
            print(e.args);
            raise e;
        except OSError as e:
            print(e.args);
            raise e;
        
        #Create output file path
        try:
            if path.exists(runParameters["output_dir"]) == False:
                makedirs(runParameters["output_dir"]);
        except OSError as e:
            print("Couldn't create output directory '%s'. Do you have write access?" % runParameters["output_dir"]);
            print(type(e), e.args);
        
        #Create fluxengine object to use runParameters
        fe = fe_obj_from_run_parameters(runParameters, metadata, processLayersOff, verbose=False);
        
        
        #Run fluxengine            
        if fe != None:
            returnCode = fe.run();
        
        #Check for successful run, if one fails don't run the rest.
        if returnCode != 0:
            print(("There was an error running flux engine: "+returnCode+"\n"%function));
            print("Exiting...");
            return (returnCode, fe);
        else:
            print("Flux engine exited with exit code:", returnCode);
            print("%02d"%timePoint.day, calendar.month_abbr[timePoint.month], timePoint.year, "%02d:%02d:%02d completed successfully.\n"%(timePoint.hour, timePoint.minute, timePoint.second));

    return (returnCode, fe); #return code, FluxEngine object.


##returns a dictionary the initialisation requirements for each rate parameterisation object (classes derived from KCalculationBase).
#def list_available_rate_parameterisations():
#    #function = inspect.stack()[0][1]+", "+inspect.stack()[0][3];
#    
#    rateParameterisations = {};
#    #try:
#    #iterate through each element of k_params and append their names
#    for name, obj in inspect.getmembers(k_params):
#        if inspect.isclass(obj) == True:
#            #Get a handle to the class type
#            #try:
#            #check it is a rate parameterisation class
#            ClassHandle = getattr(k_params, name);
#            if issubclass(ClassHandle, k_params.KCalculationBase):
#                #get initialiser arguments
#                initialiserArgNames = [arg for arg in inspect.getargspec(ClassHandle.__init__).args if arg != "self"];
#                rateParameterisations[name] = initialiserArgNames;
#                        
#    return rateParameterisations;



