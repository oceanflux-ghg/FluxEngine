#! /bin/csh -f
# To run -> /home/cercache/project/oceanflux-shared/src/climatology/GitHub/validation_test/net-flux-validation.sh
 # OceanFlux GHG uncertainties
 # calculating net fluxes IGA, adapted from Jamie Shutler

set RESULTS=/Your filepath
set LOGS=$RESULTS/logs/
set LAND_FILE=/Your fliepath/onedeg_land.nc
set MASK_FILE=/Your Filepath/World_Seas-final-complete.nc

 # setting up python
 # access virtual env if required
#source /your filepath/scientific_toolbox_cloudphys_precise/bin/activate.csh

/your filepath/GitHub/ofluxghg-flux-budgets.py -d $RESULTS -v -lf $LAND_FILE -mf $MASK_FILE -o $RESULTS/ofluxghg-example