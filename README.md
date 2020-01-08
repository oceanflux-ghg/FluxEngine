FluxEngine
==========

The FluxEngine is an open source atmosphere-ocean gas flux data processing toolbox. The toolbox has so far contributed to 12 different journal publications, resulting in 4 press releases, contributed to 2 completed PhDs and 1 ongoing PhD, has been used within 5 UK and EU research projects and has been used in undergraduate and masters level teaching. It is now being used within the European Integrated Carbon Observing System (ICOS).  This work collectively identifies and quantifies the importance of the oceans in regulating and storing carbon. 

Version 3.0 (static as of 02 August 2019).
----
v3.0 (first release April 2018, updated September 2018, February 2019, April 2019, June 2019, July 2019, static 02 August 2019)

Version 3 (v3.0) has been verified against reference runs using SOCATv4 pCO2 and all results were identical to those produced using FluxEngine v2.0. A more comprehensive verification has been performed using references runs of the Takahashi et al. (2009) dataset as described in Shutler et al. (2016) http://journals.ametsoc.org/doi/abs/10.1175/JTECH-D-14-00204.1. All results were identical to those produced using v1.0 and v2.0. A journal paper describing the v3.0 updates is now available, Holding et al., (2019) and can be found here https://doi.org/10.5194/os-15-1707-2019. 

Please reference these journal publications when using this toolbox and presenting its output in any publications.

The FluxEngine v3.0 updates and extensions were funded by the European Space Agency (ESA) research projects (OceanFlux Evolution, SKIM SciSoc) and two European Union (EU) research projects (Ringo and Integral). The two EU studies are preparatory projects for the European Integrated Carbon Observing System (ICOS).  v3.0 additions to the toolbox include:

    •   A more flexible way of specifying input data in the configuration files.
    •   Data pre-processing options (e.g. unit conversion).
    •   Python is used for all tools, allowing a more streamlined workflow.
    •   A move toward an API-like toolkit, beyond a simple set of commandline tools.
    •   A more modularised structure to the code including modular k parameterisation and data pre-processing options.
    •   Metadata and default options specified in an xml file (settings.xml).
    •   Automatic verification scripts that use SOCATv4 and Takahashi09 reference datasets.
    •   Tools for simplifying analysis of in situ data (e.g. SOCAT format data from research cruises and fixed stations).
    •   Improvements for calculating N2O and CH4 gas fluxes (now using MOMENTO data format).

Example news articles (resulting from research performed using the FluxEngine)
----
1. Phys.org (2019) 'Satellites are key to monitoring ocean carbon' https://phys.org/news/2019-11-satellites-key-ocean-carbon.html

1. ESA news (2019), 'Can oceans turn the tide on the climate crisis' https://www.esa.int/Our_Activities/Observing_the_Earth/Can_oceans_turn_the_tide_on_the_climate_crisis

2. The Guardian (2018), 'Invisible scum on sea cuts CO2 exchange with air by up to 50%' https://www.theguardian.com/environment/2018/may/28/invisible-scum-on-sea-cuts-co2-exchange-with-air-by-up-to-50

3. BBC News (2016), 'How Northern European waters soak up carbon dioxide' https://www.bbc.co.uk/news/science-environment-35654938


Animation 
----
A short animation explaining the concepts of atmosphere-ocean gas exchange, why this is important, and what the FluxEngine enables 
http://due.esrin.esa.int/stse/videos/page_video013.php


Journal publications (which use FluxEngine and/or FluxEngine outputs)
----
1. Kitidis, V., Shutler, JD., Ashton, I. et al (2020) Winter weather controls net influx of atmospheric CO2 on the north-west European shelf, Scientific Reports, 9(20153), https://www.nature.com/articles/s41598-019-56363-5

2. Shutler, JD, Wanninkhof, R, Nightingale, PD, Woolf, DK, Bakker, DCE, Watson, A, Ashton, I, Holding, T, Chapron, B, Quilfen, Y, Fairall, C, Schuster, U, Nakajima, M, Donlon, D (2019) Satellites are critical for addressing critical science priorities for quantifying ocean carbon, Frontiers in Ecology and Environment,  https://doi.org/10.1002/fee.2129

3. Woolf DK, Shutler JD, Goddijn‐Murphy L, Watson AJ, Chapron B, Nightingale PD, Donlon CJ, Piskozub J, Yelland MJ, Ashton I, et al (2019). Key Uncertainties in the Recent Air‐Sea Flux of CO2. Global Biogeochemical Cycles, https://doi.org/10.1029/2018GB006041

4. Holding, T., Ashton, I. G., Shutler, J. D., Land, P. E., Nightingale, P. D., Rees, A. P., Brown, I., Piolle, J.-F., Kock, A., Bange, H. W., Woolf, D. K., Goddijn-Murphy, L., Pereira, R., Paul, F., Girand-Ardhuin, F., Chapron, B., Rehder, G., Ardhuin, F., Donlon, C. J. (2019) The FluxEngine air–sea gas flux toolbox: simplified interface and extensions for in situ analyses and multiple sparingly soluble gases, Ocean Sci., https://doi.org/10.5194/os-15-1707-2019

5. Henson SA, Humphreys MP, Land PE, Shutler JD, Goddijn-Murphy L, Warren M (2018). Controls on open-ocean North Atlantic ΔpCO2 at seasonal and interannual timescales are different. Geophysical Research Letters, doi:10.1029/2018GL078797

6. Pereira R, Ashton, I, Sabbaghzadeh, B, Shutler, JD and Upstill-Goddard RC (2018). Reduced air–sea CO2 exchange in the Atlantic Ocean due to biological surfactants. Nature Geoscience, 1. doi: 10.1038/s41561-018-0136-2

7. Holding T, Ashton I, Woolf DK, Shutler JD (2018): FluxEngine v2.0 and v3.0 reference and verification data, PANGAEA, doi: 10.1594/PANGAEA.890118

8. Wrobel, I. (2017) Monthly dynamics of carbon dioxide exchange across the sea surface of the Arctic Ocean in response to changes in gas transfer velocity and partial pressure of CO2 in 2010. Oceanologia, 59(4), 445-459, doi: 10.1016/j.oceano.2017.05.001.

9. Ashton IG, Shutler JD, Land PE, Woolf DK, Quartly GD (2016), A Sensitivity Analysis of the Impact of Rain on Regional and Global Sea-Air Fluxes of CO2. PLoS ONE 11(9): e0161105. doi: 10.1371/journal.pone.0161105.

10. Wrobel I, Piskozub J (2016) Effect of gas-transfer velocity parameterization choice on air–sea CO2 fluxes in the North Atlantic Ocean and the European Arctic, Ocean Science, 12, 1091-1103, doi: 10.5194/os-12-1091-2016.

11. Shutler JD, Land PE, Piolle J-F, Woolf DK, Goddijn-Murphy L, Paul F, Girard-Ardhuin F, Chapron B, Donlon CJ (2016), FluxEngine: a flexible processing system for calculating atmosphere-ocean carbon dioxide gas fluxes and climatologies, Journal of Atmospheric and Oceanic Technology, doi: 10.1175/JTECH-D-14-00204.1

12. Rödenbeck C, Bakker DCE, Gruber N, Iida Y, Jacobson AR, Jones S, Landschützer P, Metzl N, Nakaoka S, Olsen A, Park G-H, Peylin P, Rodgers KB, Sasse TP, Schuster U, Shutler JD, Valsala V, Wanninkhof R, and Zeng J (2015) Data-based estimates of the ocean carbon sink variability – first results of the Surface Ocean pCO2 Mapping intercomparison (SOCOM), Biogeosciences, 12, 7251-7278, doi: 10.5194/bg-12-7251-2015.



Information about older versions of the FluxEngine toolbox
----

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

