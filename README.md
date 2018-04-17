FluxEngine
==========

v1.0 (09 March 2016)
----
The FluxEngine open source atmosphere-ocean gas flux data processing tools. The license for this software toolbox can be found within this github repository.
Please reference the publication linked below when using this toolbox and presenting its output in any publications.
A journal paper describing the toolbox has been published here: http://journals.ametsoc.org/doi/abs/10.1175/JTECH-D-14-00204.1
Please send any feedback and comments to Jamie Shutler, email: j.d.shutler@exeter.ac.uk
The FluxEngine software was developed by The European Space Agency OceanFlux Greenhouse Gases Evolution project team.

v2.0 (July 2016)
----
These updates have been verified against Takahashi (2009) using the verification options within the code. All results were identical to those derived from v1.0.
The updates included contribute to further publications in preparation and further details will be posted here following publication.
The updates include improved:

    •   handling for irregular grids,
    •   handling for different gases including O2, N2O and CH4, 
    •   handling for in-situ data.

Specifically, data on irregular grids can now be handled through the main flux calculations. Note: the ofluxghg-flux-budgets.py tool is only valid for regular (1deg x 1deg) grids. 
In-situ data should be put in separate netCDF files and the last two digits of the filename needs to represent the month of interest as a two digit number. e.g. January -> ’01’. 
To operate the system with different gases, the appropriate switch should be changed in ofluxghg-flux-calc.py. Please use ofluxghg-flux-calc.py --help for further information.

v3.0 (April 2018)
----
These updates have been verified against reference runs using SOCATv4 pCO2 and all results were identical to those produced using FluxEngine v2.0. Additionally verification has been performed using references runs of the Takahashi et al. (2009) dataset as described in Shutler et al. (2015). Results were identical to those produced using FluxEngine v1 and FluxEngine v2.
Changes include:

    •   A more flexible way of specifying input data in the configuration files,
    •   Data pre-processing options (e.g. unit conversion),
    •   Python is used for all tools, allowing a more streamlined workflow,
    •   A move toward an API-like toolkit, beyond a simple set of commandline tools
    •   A more modularised structure to the code including modular k parameterisation and data pre-processing to facilitate easy extension
    •   Metadata and default options specified in an xml file (settings.xml)
    •   Automatic validation scripts for SOCATv4 and Takahashi09


