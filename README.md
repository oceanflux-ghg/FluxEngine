FluxEngine
==========


v3.0 (first release April 2018, updated September 2018, February 2019, April 2019, June 2019, July 2019)
Version 3.0 fixed as of 02 August 2019.
----
These updates have been verified against reference runs using SOCATv4 pCO2 and all results were identical to those produced using FluxEngine v2.0. A more comprehensive verification has been performed using references runs of the Takahashi et al. (2009) dataset as described in Shutler et al. (2016) http://journals.ametsoc.org/doi/abs/10.1175/JTECH-D-14-00204.1. All results were identical to those produced using v1.0 and v2.0. A journal paper describing the v3.0 updates is now available, Holding et al., (in-review) and can be found here https://www.ocean-sci-discuss.net/os-2019-45/ . The FluxEngine v3.0 updates and extensions were funded by the European Space Agency and two European Union research projects, Ringo and Integral, which are preparatory projects for the European Integrated Carbon Observing System (ICOS).  v3.0 extensions to the toolbox include:

    •   A more flexible way of specifying input data in the configuration files.
    •   Data pre-processing options (e.g. unit conversion).
    •   Python is used for all tools, allowing a more streamlined workflow.
    •   A move toward an API-like toolkit, beyond a simple set of commandline tools.
    •   A more modularised structure to the code including modular k parameterisation and data pre-processing options.
    •   Metadata and default options specified in an xml file (settings.xml).
    •   Automatic verification scripts that use SOCATv4 and Takahashi09 reference datasets.
    •   Tools for simplifying analysis of in situ data (e.g. SOCAT format data from research cruises and fixed stations).
    •   Improvements for calculating N2O and CH4 gas fluxes (now using MOMENTO data format).


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


v1.0 (09 March 2016)
----
The FluxEngine open source atmosphere-ocean gas flux data processing tools. The license for this software toolbox can be found within this github repository.
Please reference the publication linked below when using this toolbox and presenting its output in any publications.
A journal paper describing the toolbox has been published here: Shutler et al., (2016) http://journals.ametsoc.org/doi/abs/10.1175/JTECH-D-14-00204.1
Please send any feedback and comments to Jamie Shutler, email: j.d.shutler@exeter.ac.uk
The FluxEngine software was originally developed by The European Space Agency OceanFlux Greenhouse Gases and Evolution project teams.


Publications which use FluxEngine and/or FluxEngine outputs:
----

News Article: 'How Northern European waters soak up carbon dioxide' https://www.bbc.co.uk/news/science-environment-35654938

Animation explaining the concepts of atmosphere-ocean gas exchange, why this is important and what the FluxEngine enables 
http://due.esrin.esa.int/stse/videos/page_video013.php

1. Holding, T., Ashton, I. G., Shutler, J. D., Land, P. E., Nightingale, P. D., Rees, A. P., Brown, I., Piolle, J.-F., Kock, A., Bange, H. W., Woolf, D. K., Goddijn-Murphy, L., Pereira, R., Paul, F., Girand-Ardhuin, F., Chapron, B., Rehder, G., Ardhuin, F., Donlon, C. J. (in-review) The FluxEngine air-sea gas flux toolbox: simplified interface and extensions for in situ analyses and multiple sparingly soluble gases, Ocean Sci. Discuss., https://doi.org/10.5194/os-2019-45

2. Henson SA, Humphreys MP, Land PE, Shutler JD, Goddijn-Murphy L, Warren M (2018). Controls on open-ocean North Atlantic ΔpCO2 at seasonal and interannual timescales are different. Geophysical Research Letters, doi:10.1029/2018GL078797

3. Pereira R, Ashton, I, Sabbaghzadeh, B, Shutler, JD and Upstill-Goddard RC (2018). Reduced air–sea CO2 exchange in the Atlantic Ocean due to biological surfactants. Nature Geoscience, 1. doi: 10.1038/s41561-018-0136-2

4. Holding T, Ashton I, Woolf DK, Shutler JD (2018): FluxEngine v2.0 and v3.0 reference and verification data, PANGAEA, doi: 10.1594/PANGAEA.890118

5. Wrobel, I. (2017) Monthly dynamics of carbon dioxide exchange across the sea surface of the Arctic Ocean in response to changes in gas transfer velocity and partial pressure of CO2 in 2010. Oceanologia, 59(4), 445-459, doi: 10.1016/j.oceano.2017.05.001.

6. Ashton IG, Shutler JD, Land PE, Woolf DK, Quartly GD (2016), A Sensitivity Analysis of the Impact of Rain on Regional and Global Sea-Air Fluxes of CO2. PLoS ONE 11(9): e0161105. doi: 10.1371/journal.pone.0161105.

7. Wrobel I, Piskozub J (2016) Effect of gas-transfer velocity parameterization choice on air–sea CO2 fluxes in the North Atlantic Ocean and the European Arctic, Ocean Science, 12, 1091-1103, doi: 10.5194/os-12-1091-2016.

8. Shutler JD, Land PE, Piolle J-F, Woolf DK, Goddijn-Murphy L, Paul F, Girard-Ardhuin F, Chapron B, Donlon CJ (2016), FluxEngine: a flexible processing system for calculating atmosphere-ocean carbon dioxide gas fluxes and climatologies, Journal of Atmospheric and Oceanic Technology, doi: 10.1175/JTECH-D-14-00204.1

9. Rödenbeck C, Bakker DCE, Gruber N, Iida Y, Jacobson AR, Jones S, Landschützer P, Metzl N, Nakaoka S, Olsen A, Park G-H, Peylin P, Rodgers KB, Sasse TP, Schuster U, Shutler JD, Valsala V, Wanninkhof R, and Zeng J (2015) Data-based estimates of the ocean carbon sink variability – first results of the Surface Ocean pCO2 Mapping intercomparison (SOCOM), Biogeosciences, 12, 7251-7278, doi: 10.5194/bg-12-7251-2015.



