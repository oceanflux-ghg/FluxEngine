#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 15:15:53 2019

@author: Tom Holding
"""

import setuptools;
import os.path as path;

#List the standalone scripts
scripts = ["append2insitu.py", "compare_net_budgets.py",
           "fluxengine_tutorials.py", "ncdf2text.py",
           "ofluxghg_flux_budgets.py", "ofluxghg_run.py",
           "reanalyse_socat_driver.py", "resample_netcdf.py",
           "full_fluxengine_verification.py", "text2ncdf.py",
           "verify_socatv4.py", "verify_takahashi09.py"];
for i in range(len(scripts)): #Append the directory to each script name
    scripts[i] = "fluxengine/scripts/"+scripts[i];

#list the dependencies
dependencies = ["numpy>=1.16.4",
                "pandas>=0.24.2",
                "matplotlib>=3.1.0",
                "netCDF4>=1.4.2",
                "argparse>=1.1",
                "scipy>=1.3.0",
                "jupyter>=1.0.0",
                ]

#Create a list of package data files
dataPaths = setuptools.findall(dir="fluxengine/data"); #list all data files
configPaths = setuptools.findall(dir="fluxengine/configs"); #list all configuration files
tutorialPaths = setuptools.findall(dir="fluxengine/tutorials"); #list all tutorial files
allPackageDataPaths = dataPaths+configPaths+tutorialPaths; #Combine to get a list of all the data paths
for i in range(len(allPackageDataPaths)): #remove the "fluxengine" parent directory
    allPackageDataPaths[i] = path.join( *(allPackageDataPaths[i].split(path.sep)[1:]));

#Append special cases to allPackageDataPaths
allPackageDataPaths.append("core/settings.xml");


#Main setup object.
setuptools.setup(
    name="FluxEngine",
    version="4.0.dev0",
    author="Tom Holding, Jamie Shutler and others",
    author_email="t.m.holding@exeter.ac.uk, j.d.shutler@exeter.ac.uk",
    description="Open-source toolkit for calculating atmosphere-ocean gas transfer",
    url="https://github.com/oceanflux-ghg/FluxEngine",
    license='Creative Commons Attribution-Noncommercial-Share Alike license',
    project_urls={"Source Code":"https://github.com/oceanflux-ghg/FluxEngine",
                  "Ancillary Tools":"https://github.com/oceanflux-ghg/FluxEngineAncillaryTools",
                  "More Info":"http://www.oceanflux-ghg.org/Products/FluxEngine"},
    
    packages=setuptools.find_packages("."),  # include all packages in this dir "fluxengine"
    package_dir={"":"."},   # tell distutils packages in this dir "fluxengine"
    scripts=scripts,#["fluxengine/scripts/text2ncdf.py", "fluxengine/scripts/verify_takahashi09.py"],
    include_package_data=True,
    package_data={"fluxengine": allPackageDataPaths},
    
    install_requires=dependencies,
)



#from distutils.core import setup;
#
#scripts = ["append2insity.py", "compare_net_budgets.py",
#           "ofluxghg_flux_budgets.py", "ofluxghg_run.py",
#           "reanalyse_socat_driver.py","text2ncdf.py", "text2ncdf_examples.sh",
#           "verify_socat4_sst_salinity_gradients_N00.py",
#           "verify_takahashi09.py"];
#for i in range(len(scripts)): #Append the directory to each script name
#    scripts[i] = "fluxengine/scripts/"+scripts[i];
#
#setup(
#    name='FluxEngine',
#    version='4.0dev',
#    packages=['fluxengine',],
#    scripts=[scripts],
#    license='Creative Commons Attribution-Noncommercial-Share Alike license',
#    long_description=open('README.md').read(),
#    include_package_data=True,
#    install_requires=[],
#    
#)
