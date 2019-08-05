#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 23 09:37:32 2018



@author: tomholding
"""
from .datalayer import DataLayerMetaData;

#A place to store default information about DataLayers.
class Settings:
    #filename
    def __init__(self, filename, verbose=False):
        import xml.etree.ElementTree as ET; #Only needed here.
        self.allDataLayers = {};
        
        print("Parsing settings file at:", filename);
        tree = ET.parse(filename);
        root = tree.getroot();
        
        #Read input data layers
        dataLayersElement = root.find("DataLayers");
        if dataLayersElement != None:
            for element in dataLayersElement:
                if element.tag == "DataLayer":
                    metaData = DataLayerMetaData(**element.attrib);
                    self.allDataLayers[metaData.name] = metaData;
        else:
            raise IOError("Couldn't find parent 'DataLayers' element in settings.xml because the file is improperly formatted.");