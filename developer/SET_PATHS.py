#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 20 10:00:52 2018
Generates a complete set of reference runs:
    Takahashi09 - all inputs to reproduce results from Takahashi 2009
    Takahashi09 pCO2: #Not implemented
        no_gradients_N00 #Not implemented
        sst_gradients_N00 #Not implemented
        sst_gradients_k0 #Not implemented
        sst_gradients_k1 #Not implemented
        sst_gradients_k2 #Not implemented
        sst_gradients_k3 #Not implemented
    SOCATv4 pCO2:
        no_gradients_N00
        sst_gradients_N00
        sst_gradients_k0
        sst_gradients_k1
        sst_gradients_k2
        sst_gradients_k3
    
@author: tomholding
"""

#Determines the root directory of the particular version of FE being ran.
#Extracts source directory and inserts it as a priority search path
#to ensure the current in-development version is imported and not
#any other installed versions.
def set_paths(setPath=True, verbose=False):
    print "__FILE__:", __file__;
    import sys;
    #import inspect;
    from os import path, chdir;
    
    #Developer scripts add fluxengine path to sys.path to ensure the in-development version of fluxengine is being used
    #rather than any other installed version.
    #rootPath = path.dirname(path.dirname(path.abspath(__file__))); #Must make abs first or when ran in interpretter it fails.
    rootPath = path.abspath(path.join(__file__, path.abspath(".."))); #Must make __file__ abs first or when ran in interpretter it fails.
    #rootPath = path.dirname(path.dirname((inspect.stack()[0][1])));
    srcPath = path.join(rootPath, "fluxengine_src");
    if setPath==True and srcPath not in sys.path:
        sys.path.insert(1, srcPath);
    if verbose:
        print "rootPath: ", rootPath;
        print "srcPath: ", srcPath;
        print "srcPath added to sys.path\n";
    
    #TODO: This is not great. We need a more systematic way to manage resources and paths.
    chdir(rootPath); #when running developer scripts, set the cwd to the FluxEngine root folder...
    
    return (rootPath, srcPath);

#import fe_setup_tools as setup;




