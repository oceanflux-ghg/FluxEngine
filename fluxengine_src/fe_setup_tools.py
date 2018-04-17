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

from os import path, listdir;
import inspect;
import time;
import socket; #for gethostname
import calendar;
import fnmatch; #matching file globs

import fe_core as fluxengine
import rate_parameterisation as k_params; #This is where k parameterisation logic is kept.
import data_preprocessing as data_preprocessing; #preprocessing functions
import process_indicator_layers as indicator_layers; #Process indicator layer functors


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
                        print "Warning: Duplicate definition in config file for '%s'. The first definition will be used." % name;
                except: #Will only happen if there isn't both left and righthand sides to the '=' assignment
                    print "%s: Error parsing config file line:\n"%function, line;
                    return None;
    except Exception as e:
        print "Error while parsing config file at path:", configPath;
        print type(e), e.args;
        return None;
    
    if verbose:
        print "\nParsed configurable variables as follows:";
        for key in configVariables.keys():
            print key, configVariables[key];
    return configVariables;


#Load and return metadata about the variables contained in config files.
#Returns a dictionary of dictionarys allowing indexing by variablename then attribute, e.g.:
#   varMetadata[varName][attribute]
def read_config_metadata(settingsPath, verbose=False):
    function = inspect.stack()[0][1]+", "+inspect.stack()[0][3];
    
    import xml.etree.ElementTree as ET; #Only needed here.
    
    if verbose:
        print "Parsing settings file at:", settingsPath;
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
                                    print "%s: Invalid multioption child 'Option' element in %s. 'str' and 'value' attributes must be specified." % (function, element.attrib["name"])
                                    print e.args;
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
                    else: #just a straightforward variable, so add it.
                        varMetadata[element.attrib["name"]] = element.attrib;
            except ValueError as e:
                print "%s: 'name' or 'type' attribute missing from Variable element in settings xml file. Specified settings file is invalid: %s" % (function, settingsPath);
                print e.args;
                
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
        if metadata[varName]["type"] == "path" or metadata[varName]["type"] == "DataLayerPath":
            if path.isabs(configVariables[varName]) == False:
                configVariables[varName] = path.abspath(configVariables[varName]);
        
        #Multioption: Check that the config value matches at least one option.
        #             Store string value with _name suffix overwrite original with integer value.
        elif metadata[varName]["type"] == "multioption":
            if configVariables[varName] not in metadata[varName]["options"]:
                raise ValueError("%s: Invalid option specified in config file for %s. Valid options are: %s" % (function, varName, str(metadata[varName]["options"].keys())));
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
                print "%s: Config variable '%s' requires a decimal / floating point number. Got %s instead." % (function, varName, configVariables[varName]);
        
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
    
    #Now process custom vars. Try to convert them to a float, but if they fail assume they're supposed to be a string.
    for varName in customVars:
        try:
            configVariables[varName] = float(configVariables[varName]);
        except ValueError:
            pass; #Leave it as a string.
    
    #Various misc. checks
    if (configVariables["use_sstfnd"]==True) and ("sstfnd_path" not in configVariables):
        raise ValueError("%s: use_sstfnd is set but no sstfnd inputfile (sstfnd_path was specified in the config file." % function);
    
    if (configVariables["use_sstskin"]==True) and ("sstskin_path" not in configVariables):
        raise ValueError("%s: use_sstskin is set but no sstskin inputfile (sstskin_path) was specified in the config file." % function);
    
    if (configVariables["sst_gradients"]==True) and ("sstgrad_path" not in configVariables):
        raise ValueError("%s: use_gradients is set but no sstgrad inputfile (sstgrad_path) was specified in the config file: " % function);

#Substitutes month/year tokens to create a valid glob.
#Returns the a tuple containing (searchDirectory, fileGlob)
def generate_glob(globSig, year, monthNum):
    glob = globSig;
    glob = glob.replace("<YYYY>", str(year));
    glob = glob.replace("<YY>", str(year)[-2:]);
    glob = glob.replace("<MM>", "%02d"%(monthNum+1));
    glob = glob.replace("<MMM>", calendar.month_abbr[monthNum+1].upper());
    glob = glob.replace("<Mmm>", calendar.month_abbr[monthNum+1]);
    glob = glob.replace("<mmm>", calendar.month_abbr[monthNum+1].lower());
    
    glob = glob.rsplit("/", 1);
    #print "Split result:", glob;
    
    return (glob[0], glob[1]);


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
#Converts between inconsistencies in the config file naming convention and FluxEngine nameing convention, e.g.:
#   path -> infile
#   appending _switch to bool variables
def create_run_parameters(configVariables, varMetadata, year, monthNum, processTimeStr, configFile, processIndicatorLayersOff):
    function = inspect.stack()[0][1]+", "+inspect.stack()[0][3];
    
    runParams = {};
    #Copy over misc. parameters
    #runParams["year"] = 2000 if clArgs.use_takahashi_validation==True else runParams["year"] = year;
    runParams["year"] = year;
    runParams["hostname"] = socket.gethostname();
    runParams["config_file"] = configFile;
    runParams["src_home"] = configVariables["src_home"];
    runParams["processing_time"] = processTimeStr;
    runParams["monthyear_processing_time"] = time.strftime("%d/%m/%Y %H:%M:%S");
    runParams["process_layers_off"] = int(processIndicatorLayersOff);
    
    #Copy over all config parameters:
    missingFiles = []; #Much more useful to the user to have a complete list of missing files instead of one at a time.
    for varName in configVariables:
        #Different variable types need to be copied differently
        #DataLayerPaths: These can change based on month and year.
        if varName in varMetadata:
            if varMetadata[varName]["type"] == varMetadata[varName]["type"] == "DataLayerPath":
                curPath, curGlob = generate_glob(configVariables[varName], year, monthNum);
                #print "path:", curPath
                #print "glob:", curGlob;
                try:
                    matches = match_filenames(curPath, curGlob);
                    matches = [match for match in matches if match.rsplit('/', 1)[-1][0] != '.']; #Remove hidden metadata files created by MACOS on some file systems. Required if sharing folders between MACOS and LINUX/WINDOWS.
                    
                    #If there's exactly one match for the file glob, add the data layer.
                    if len(matches) == 1: #Store the file path.
                        dataLayerName = varName[0:-5];
                        runParams[dataLayerName+"_infile"] = path.join(curPath, matches[0]);
                        runParams[dataLayerName+"_prod"] = configVariables[dataLayerName+"_prod"];
                        runParams[dataLayerName+"_stddev_prod"] = configVariables[dataLayerName+"_stddev_prod"] if dataLayerName+"_stddev_prod" in configVariables else None;
                        runParams[dataLayerName+"_count_prod"] = configVariables[dataLayerName+"_count_prod"] if dataLayerName+"_count_prod" in configVariables else None;
                        
                        if dataLayerName+"_preprocessing" in configVariables:
                            runParams[dataLayerName+"_preprocessing"] = configVariables[dataLayerName+"_preprocessing"];
                        else:
                            runParams[dataLayerName+"_preprocessing"] = "no_preprocessing";
                        
                    else: #There isn't exactly 1 matching file, so we have a problem, add it to the list.
                        missingFiles.append( (varName, path.join(curPath, curGlob) ) );
                except Exception as e:
                    print "Invalid filepath when checking for %s: %s" % (varName, curPath);
                    print type(e), e.args;
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
            print "Missing files:";
            for missingFile in missingFiles:
                print missingFile[0]+": "+missingFile[1];
            raise ValueError("%s: Error finding unique file matches for one or more input data layers." % (function));
            return None;
            

    #Output file
    outFile = "OceanFluxGHG-month%02d-%s-%d-v0.nc" % (monthNum+1, calendar.month_abbr[monthNum+1].lower(), year);
#    if clArgs.use_takahashi_validation==True:
#        ##TODO: This needs to be corrected by passing "start_year 2000 end_year 2000" in the command line arguments once the correct ice and pressure data are being used.
#        outDir = path.join(configVariables["output_dir"], str(2000), "%02d"%(monthNum+1));
#    else:
    outDir = path.join(configVariables["output_dir"], str(year), "%02d"%(monthNum+1));

    runParams["output_dir"] = outDir;
    runParams["output_path"] = path.join(outDir, outFile);
    
    
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
                    print "%s: Could not find a k parameterisation functor which corresponds to the specified k_parameterisation (%s). If this specified correctly in the config file?"%(function, name);
                    print e.args;
                    return None;
                
                #Check it derives from KCalculateBase - this is used as a simple way to filter out irrelevant classes
                if issubclass(ClassHandle, k_params.KCalculationBase):
                    #To initialise the class we need to know the arguments it's __init__ function uses.
                    initialiserArgNames = [arg for arg in inspect.getargspec(ClassHandle.__init__).args if arg != "self"];
                    
                    try:
                        #Create a dictionary of name:value pairs for the __init__ arguments. Assumes arguments match config file names.
                        argDict = {key : runParameters[key] for key in initialiserArgNames}; #create a dictionary of arguments name:value pairs
                    except KeyError as e:
                        print "%s: Could not find all the required initialiser arguments. Are they specified correctly in the config file?\nExpected arguments are: "%function, initialiserArgNames;
                        print "KeyError.args: ", e.args;
                        return None;
                    
                    #Finally create and return the k functor instance
                    return ClassHandle(**argDict);
    
    except ImportError as e:
        print "%s: Cannot find k_parameterisation module. Check FluxEngine is installed correctly.", e.args;
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
        fe.add_k_parameterisation_component(k_params.AddKRainLinearHo19997());
    
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
        print "Error adding data layers.";
        fe = None;
    
    return fe;


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



