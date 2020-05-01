#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 15:15:53 2019

@author: Tom Holding
"""

import setuptools;
import os.path as path;

#List the standalone scripts
scripts = ["fe_append2insitu.py", "fe_compare_net_budgets.py",
           "fe_tutorials.py", "fe_ncdf2text.py",
           "fe_calc_budgets.py", "fe_run.py",
           "fe_reanalyse_fco2_driver.py", "fe_resample_netcdf.py",
           "fe_full_verification.py", "fe_text2ncdf.py",
           "fe_verify_socatv4.py", "fe_verify_takahashi09.py",
           "fe_update_config.py"];
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
                ];

#Create a list of package data files
dataPaths = setuptools.findall(dir="fluxengine/data"); #list all data files
configPaths = setuptools.findall(dir="fluxengine/configs"); #list all configuration files
tutorialPaths = setuptools.findall(dir="fluxengine/tutorials"); #list all tutorial files
allPackageDataPaths = dataPaths+configPaths+tutorialPaths; #Combine to get a list of all the data paths
for i in range(len(allPackageDataPaths)): #remove the "fluxengine" parent directory
    allPackageDataPaths[i] = path.join( *(allPackageDataPaths[i].split(path.sep)[1:]));

#Append special cases to allPackageDataPaths
allPackageDataPaths.append("core/settings.xml");


#Read description
with open(path.join('README.md'), encoding='utf-8') as f:
    longDescription = f.read();

#Main setup object.
setuptools.setup(
    name="fluxengine",
    version="4.0.0",
    author="Tom Holding and Jamie Shutler",
    author_email="t.m.holding@exeter.ac.uk, j.d.shutler@exeter.ac.uk",
    description="Open-source toolkit for calculating atmosphere-ocean gas transfer",
    long_description=longDescription,
    long_description_content_type="text/markdown",
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
    );