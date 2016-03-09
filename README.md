FluxEngine
==========

open source atmosphere-ocean gas flux data processing tools.

FluxEngine source code and associated tools are contained within this repository

For more in formation, a journal paper describing the toolbox has now been accepted for publication and an early version of the manuscript is avaialable here: http://journals.ametsoc.org/doi/abs/10.1175/JTECH-D-14-00204.1

The FluxEngine and OceanFlux Greenhouse Gases project team (9th March 2016)

FE_license - license file
Oceans-combined-final-complete.nc - netCDF file containing the ocean regions. Serves as a mask file for the flux-budgets tool
ofluxghg-example_config.conf - an example config file used to run FE
ofluxghg-example_flux-budgets_runfile.sh - example runfile for calculating net values using the flux-budgets tool. Requires onedeg_land.nc and Oceans-combined-final-complete.nc as land file and mask file respectively.
ofluxghg-flux-budgets.py - flux-budgets tool used to calculate the net global and regional flux using the results from the flux-calc tool.  
ofluxghg-flux-calc.py - flux calc tool is the main FE code to calculate air-sea flux
ofluxghg-run-climatology.pl - PERL wrapper to run flux-calc tool from configuration file
ofluxghg-example-flux-calcs_runfile.sh - runfile to run flux-calc tool using config file and PERL wrapper
resample_netcdf.py - tool to resample in-situ data to the correct grid sizing for FE
text2ncdf.py - tool to read in-situ data and export in correct format for FE (netCDF)
