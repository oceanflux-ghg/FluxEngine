#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from netCDF4 import Dataset
from numpy import flipud, ma, ravel, transpose, full, array, squeeze;
from numpy import all as npall;
from .debug_tools import calc_mean;

#Simple class for storing metadata about a datalayer.
#These values are used as the default values for datalayers, unless they're overwritten by the config file.
class DataLayerMetaData:
#    def __init__(self, argDict):
#        for argName in argDict:
#            setattr(self, argName, argDict[argName]);
    def __init__(self, name, netCDFName=None, units=None, minBound=None, maxBound=None, standardName=None, longName=None, fillValue=None):
        function = "(DataLayerMetaData.__init__)";
        self.name = name;
        self.netCDFName = netCDFName if netCDFName != None else name;
        self.units = units if units != None else name;
        self.standardName = standardName if standardName != None else name;
        self.longName = longName if longName != None else name;
        self.fillValue = fillValue; #None is allowed
        
        try:
            self.minBound = float(minBound) if minBound != None else None; #None is allowed
        except ValueError:
            print("%s: Invalid minimum value supplied (%s) for DataLayer %s, defaulting to None." % {function, minBound, self.name});
            self.minBound = None;
        
        try:
            self.maxBound = float(maxBound) if maxBound != None else None; #None is allowed
        except ValueError:
            print("%s: Invalid maximum value supplied (%s) for DataLayer %s, defaulting to None." % (function, maxBound, self.name));
            self.maxBound = None;
        
        self.temporalChunking = 1;
        self.temporalSkipInterval = 0;
        self.timeDimensionName = "time";


#Encapsulates of an input data layer, and manages checking it for integrity etc.
class DataLayer:
    #Class fields
    missing_value = -999.0;
    fill_value = -999.0;
    missing_value_int = -999;
    fill_value_int = -999;
    DEBUG = False;
    
    @classmethod
    def create_empty_datalayer(cls, name, nx, ny, metadata, fillValue=None):
        if fillValue == None:
            fillValue = DataLayer.missing_value;
        data = full((ny, nx), fillValue);
        return cls(name, data, metadata, fillValue);
    
    #preprocessing is a list of functions to modify fdata
    #Throws IOError is netCDF file isn't found
    #Throws KeyError if variable (prod name) isn't found
    #Throws ValueError if an unexpected number of dimensions are found
    #TODO: transposeData should be handled as a preprocessing function
    @classmethod
    def create_from_file(cls, name, infile, prod, metadata, timeIndex, transposeData=False, preprocessing=None):
        function = "(DataLayer.create_from_file)"
        
        #Open netCDF file
        try:
            dataset = Dataset(infile);
        except IOError as e:
            print("\n%s: %s inputfile %s does not exist" % (function, name, infile))
            print(e.args);
        
        #Check netCDF file: Prints some info when in DEBUG mode.
        check_input(infile, prod, DataLayer.DEBUG);
        
        #Open netCDF file and check variable exists.
        dataset = Dataset(infile);
        ncVariable = dataset.variables[prod];
        
        #Find the right time dimension index and slice/copy the data appropriately
        dims = ncVariable.dimensions;
        
        #Two spatial dimensions and a time dimensions
        if len(dims) == 3:
            if metadata.timeDimensionName in dims:
                if dims.index(metadata.timeDimensionName) == 0:
                    data = ncVariable[timeIndex, :, :];
                elif dims.index(metadata.timeDimensionName) == 1:
                    data = ncVariable[:, timeIndex, :];
                elif dims.index(metadata.timeDimensionName) == 2:
                    data = ncVariable[:, :, timeIndex];
            else:
                raise RuntimeError("Time dimension name ('%s') for Datalayer '%s' was not found. Try setting this manually in the configuration file using (for example) datalayername_timeDimensionName = time"%(metadata.timeDimensionName, name));
        #No time dimension anyway
        elif len(dims) == 2:
            data = ncVariable[:];
        else: #
            raise RuntimeError("Invalid number of dimensions (%d) when reading datalayer '%s' from '%s'"%(len(dims), name, infile));
        
        #TODO: APPLY PREPROCESSING HERE instead of later.
        

        #Extract just the dimensions we want.
        #requiredDims = [None if v in ['latitude', 'lat', 'longitude', 'lon'] else 0 for v in ncVariable.dimensions]
        #data = squeeze(ncVariable[slice(*requiredDims)]);
        if data.shape == (1, 1, 1): #Remove temporal dimension
            data.shape = (1, 1);
        
        if data.shape != (1, 1): #Don't squeeze if there is a single lon/lat point or we'll end up with an empty array
            data = squeeze(data); #Remove any 1d dimensions.

        #check number of dimensions
        dataDims = data.shape;
        if len(dataDims) != 2:
            raise ValueError("\n%sError: Unexpected number of dimensions (%d) in %s when reading in %s variable from %s" % (function, len(dataDims), name, prod, infile));
        
        #Convert from a masked array (np.ma.array) to a plain np.array
        data = array(data);
        
        #Seems to be a bug in netCDF4 which sometimes causes unmasked variables to be completely masked: https://github.com/Unidata/netcdf4-python/issues/707
        #    Needs looking into, but for now this work-around fixes things:
        #    TODO: now irrelevent as we convert to standard arrays?
        if ma.isMaskedArray(data) and (npall(data.mask) == True):
              data.mask = False; #Remove the mask...
                
        #If necessary flip the data #TODO: remove this as this should be handled in the pre-processing functions by the user
        data, flipped = flip_data(dataset, data, name); #If different from takahashi orientation, flip data.

        #Extract fill value from netCDF if it exists. Note that this will overwrite fill default or config specified fill value.
        if hasattr(ncVariable, "_FillValue"):
            fillValue = ncVariable._FillValue;
        elif hasattr(ncVariable, "fill_value"):
            fillValue = ncVariable.fill_value;
        else:
            fillValue = DataLayer.missing_value;
        
        return cls(name, data, metadata, fillValue, preprocessing=preprocessing);
        
    
    def __init__(self, name, data, metadata, fillValue, preprocessing=None):
        #function = "(DataLayer.__init__)";
        self.name = name; #Human readable name for the data layer
        self.data = data;
        self.ny = data.shape[0];
        self.nx = data.shape[1];
        self.missing_value = DataLayer.missing_value;
        
        #Copy over all metadata attributes. Note: this also copies over name - not sure this is a good thing
        for attribute in vars(metadata):
            setattr(self, attribute, getattr(metadata, attribute));
        
        #TODO: this only applies to data from file, so move to create_from_file
        #replace fill value with standardised missing_value
        if fillValue != self.missing_value:
            self.data[self.data == fillValue] = self.missing_value;
        
        #Create a view of 'data' which is 1D. Note: MAY sometimes copy data. This should be checked...
        self.calculate_fdata();

        #Apply any preprocessing to the data (e.g. unit conversion)
        if preprocessing != None:
            for preprocessingFunction in preprocessing:
                preprocessingFunction(self); #Modifies in place

        #Should we recalculate_fdata after preprocessing to guard against badly written custom preprocessing functions. 
        #self.calculate_fdata(); WARNING: cannot do both this and propagate_fdata_to_data()
        
        #required in case of automatic resampling of large datalayers (which adjusts data)
        #self.propagate_fdata_to_data();
        
        
        #Replace anything outside of the valid range with missing_value
        validate_range(self.fdata, self.minBound, self.maxBound, self.missing_value);
    
    #Create a view of 'data' which is 1D. Note: Sometimes makes a copy of the data so this shouldn't re relied on.
    def calculate_fdata(self):
        self.fdata = ravel(self.data);
        if ma.is_masked(self.fdata):
            self.fdata.unshare_mask(); #TMH: masked array behaviour is changing in future versions of numpy. This avoids ambiguity between current and future behaviour.
    
    #def propagate_fdata_to_data():
    #   pass;


#Replaces elements in a matrix who's values are outside the specified range [minBound, maxBound] with missing value.
def validate_range(fdata, minBound, maxBound, missingValue):
    if minBound != None and maxBound != None: #Must handle each case seperately to support any value of missingValue.
        for i in range(0, len(fdata)):
            if fdata[i] < minBound or fdata[i] > maxBound:
                fdata[i] = missingValue;
    elif minBound != None and maxBound == None:
        for i in range(0, len(fdata)):
            if fdata[i] < minBound:
                fdata[i] = missingValue;
    elif minBound == None and maxBound != None:
        for i in range(0, len(fdata)):
            if fdata[i] > maxBound:
                fdata[i] = missingValue;
        
        


def flip_data(dataset, this_variable, name):
    '''#IGA - for a netcdf data set, determine whether latitude orientation matches 'taka' and if not, flip the variable provided using flipud'''
    try:
        data_latitude_prod = [str(x) for x in list(dataset.variables.keys()) if 'lat' in str(x).lower()] #finds the correct latitude name for data
        if len(data_latitude_prod) != 1:
            raise IndexError("Ambiguous or no latitude variable: "+str(list(dataset.variables.keys())));
        
        data_lat = dataset.variables[data_latitude_prod[0]]
        if len(data_lat.shape)<2 and data_lat[0]<0:#IGA - if true, it is a vector that is in opposite orientation to 'taka'
            flipped = True
            this_variable_out = flipud(this_variable)
        else:
            this_variable_out = this_variable
            flipped = False
    
        return this_variable_out, flipped
    except IndexError:
        print("Assuming correct orientation for %s. Variable has not been flipped." % name);
        return this_variable, False;
        

def check_input(filename, dataset, DEBUG=False):
   function = "(check_input, main)"   
   if DEBUG:
      print("%s checking %s for %s data" % (function, filename, dataset))
   file = Dataset(filename,'r',clobber=False)
   for dimname, diminst in sorted(list(file.dimensions.items())):
      if diminst.isunlimited():
         if DEBUG:
            print("%s %s dimension\t%s\t%s\tunlimited" % (function, dataset, dimname, len(diminst)))
      else:
         if DEBUG:
            print("%s %s dimension\t%s\t%s" % (function, dataset, dimname, len(diminst)))
              
   for varname, varinst in sorted(list(file.variables.items())):
      if DEBUG:
         print("%s %s variable\t%s\t\t%s\t%s\t%s" % (function, dataset, varname, varinst.shape, varinst.dimensions, varinst.dtype))

   file.close()
   return 0