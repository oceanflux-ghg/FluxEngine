#! /bin/csh -f

 # setting up python
  #If necessary, set up a Python environment and the relevant packages available. 
#source <your setup file>.csh

#set location of configutration file
set CONFIG_SOURCE=/.../...

# Set location for results (year/month directory structure will be created)
set RESULTS=/.../...

 # actual reference climatology run
time perl $CONFIG_SOURCE/ofluxghg-run-climatology.pl --config $CONFIG_SOURCE/ofluxghg-example_config.conf -l -o $RESULTS/ofluxghg_test -s 2000 -e 2000 #>&! $RESULTS/log #Optional term to direct log to a file
