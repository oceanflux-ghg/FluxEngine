#Example calls using text2ncdf (uncomment as needed)

#Show help and options
python text2ncdf.py -h

#Single input file
#python text2ncdf.py raw_data/Kitidis2017/74JC20131009_CO2_underway_SOCATv5.tab.tsv -n testout.nc -s "2013-10-09 19:42:07" -e "2013-11-08 21:32:06" -c 32 --parse_units --cols 4 5 7 13 14 -t "monthly"

#Two input files:
#python text2ncdf.py raw_data/Kitidis2017/74JC20131009_CO2_underway_SOCATv5.tab.tsv raw_data/Kosugi2017/49HH20111201_CO2_underway_SOCATv5.tab.tsv -n testout.nc -s "2011-12-01 09:19:00" -e "2013-11-08 21:32:06" -c 32 33 --parse_units --cols 4 5 7 13 14 -t "monthly"

#Lon/Lat limits:
#python text2ncdf.py raw_data/Kitidis2017/74JC20131009_CO2_underway_SOCATv5.tab.tsv -n testout.nc -s "2013-10-09 19:42:07" -e "2013-11-08 21:32:06" -c 32 --parse_units --cols 4 5 7 13 14 -t "monthly" -l -35.5 47.5 -36.5 -9.5

#Specify temporal resolution
#python text2ncdf_6.py raw_data/Kitidis2017/74JC20131009_CO2_underway_SOCATv5.tab.tsv -n testout.nc -s "2013-10-09 19:42:07" -e "2013-11-08 21:32:06" -c 32 --parse_units --cols 4 5 7 13 14 -t "6 00:00" -l -35.5 47.5 -36.5 -9.5


