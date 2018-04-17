#! /bin/csh -f
 # OceanFlux GHG uncertainties
 # calculating net fluxes IGA, adapted from Jamie Shutler

set RESULTS=/Your filepath
set LOGS=$RESULTS/logs/
set LAND_FILE=/Your fliepath/onedeg_land.nc
set MASK_FILE=/Your Filepath/World_Seas-final-complete.nc


/your filepath/ofluxghg-flux-budgets.py -d $RESULTS -v -lf $LAND_FILE -mf $MASK_FILE -o $RESULTS/ofluxghg-example
