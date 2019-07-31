#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 15:15:53 2019

@author: Tom Holding
"""

from setuptools import setup, find_packages

#List the standalone scripts
scripts = ["append2insity.py", "compare_net_budgets.py", "fluxengine_tutorials.sh",
           "netcdf2text.py", "ofluxghg_flux_budgets.py", "ofluxghg_run.py",
           "reanalyse_socat_driver.py", "resample_netcdf.py",
           "full_fluxengine_verification.py", "text2ncdf.py",
           "verify_socatv4_sst_salinity_gradients_N00.py",
           "verify_takahashi09.py"];
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

setup(
    name="FluxEngine",
    version="4.0dev",
    author="Tom Holding, Jamie Shutler and others",
    author_email="t.m.holding@exeter.ac.uk, j.d.shutler@exeter.ac.uk",
    description="Open-source toolkit for calculating atmosphere-ocean gas transfer",
    url="https://github.com/oceanflux-ghg/FluxEngine",
    project_urls={"Source Code":"https://github.com/oceanflux-ghg/FluxEngine",
                  "Ancillary Tools":"https://github.com/oceanflux-ghg/FluxEngineAncillaryTools",
                  "More Info":"http://www.oceanflux-ghg.org/Products/FluxEngine"},
    
    packages=find_packages("fluxengine"),  # include all packages under src
    package_dir={"":"fluxengine"},   # tell distutils packages are under src
    scripts=scripts,#["fluxengine/scripts/text2ncdf.py", "fluxengine/scripts/verify_takahashi09.py"],
    
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
