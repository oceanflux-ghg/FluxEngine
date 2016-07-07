FluxEngine
==========

v1.0
The FluxEngine open source atmosphere-ocean gas flux data processing tools. The license for this software toolbox can be found within this github repository.
Please reference the publication linked below when using this toolbox and presenting its output in any publications.
A journal paper describing the toolbox has been published here: http://journals.ametsoc.org/doi/abs/10.1175/JTECH-D-14-00204.1
Please send any feedback and comments to Jamie Shutler, email: j.d.shutler@exeter.ac.uk
The FluxEngine software was developed by The European Space Agency OceanFlux Greenhouse Gases Evolution project team
Jamie Shutler (09 March 2016)

Updates July 2016:
v2.0
These updates have been verified against Takahashi (2009) using the verification options within the code. All results were identical to those derived from v1.0
The updates included contribute to further publications in preparation and details will be posted here following publication.
Updates include: 
    •   Improved Handling for irregular grids,
    •   handling for different gases including O2, N2O and CH4, 
    •   handling for in-situ data.

Specifically, data on irregular grids can now be handled through the FE calculations. Note that the ofluxghg-flux-budgets.py tool is only valid for regular (1deg x 1deg) grids. 
In-situ data should be put in separate netCDF files whereby the last two digits of the filename represent the month as a two digit number. e.g. Jan -> ’01’. 
To operate the system with different gases, the appropriate switch should be changed in ofluxghg-flux-calc.py
