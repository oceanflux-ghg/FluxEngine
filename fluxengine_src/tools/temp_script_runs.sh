#Global SOCATv5 netCDF
python reanalyse_socat_driver.py -socat_dir ~/data/SOCAT_ascii/SOCATv5/ -socat_files SOCATv5.tsv -sst_dir ~/data/ocean_flux_ftp/SST_reynolds_avhrr/ -sst_tail 01_OCF-SST-GLO-1M-100-REYNOLDS.nc -output_dir ~/Files/Tasks/Reanalyse_SOCATv5_PANGAEA/output_SOCATv5_netCDF_2 -socatversion 5 -usereynolds -startyr 1957 -endyr 2017 -keepduplicates -keeptempfiles

#Global SOCATv5 ascii
python reanalyse_socat_driver.py -socat_dir ~/data/SOCAT_ascii/SOCATv5/ -socat_files SOCATv5.tsv -sst_dir ~/data/ocean_flux_ftp/SST_reynolds_avhrr/ -sst_tail 01_OCF-SST-GLO-1M-100-REYNOLDS.nc -output_dir ~/Files/Tasks/Reanalyse_SOCATv5_PANGAEA/output_SOCATv5_ascii_2 -socatversion 5 -usereynolds -startyr 1957 -endyr 2017 -keepduplicates -keeptempfiles -asciioutput



#Global SOCATv6 netCDF
python reanalyse_socat_driver.py -socat_dir ~/data/SOCAT_ascii/SOCATv6/ -socat_files SOCATv6.tsv -sst_dir ~/data/ocean_flux_ftp/SST_reynolds_avhrr/ -sst_tail 01_OCF-SST-GLO-1M-100-REYNOLDS.nc -output_dir ~/Files/Tasks/Reanalyse_SOCATv6_PANGAEA/output_SOCATv6_netCDF_2 -socatversion 6 -usereynolds -startyr 1957 -endyr 2017 -keepduplicates -keeptempfiles

#Global SOCATv6 ascii
python reanalyse_socat_driver.py -socat_dir ~/data/SOCAT_ascii/SOCATv6/ -socat_files SOCATv6.tsv -sst_dir ~/data/ocean_flux_ftp/SST_reynolds_avhrr/ -sst_tail 01_OCF-SST-GLO-1M-100-REYNOLDS.nc -output_dir ~/Files/Tasks/Reanalyse_SOCATv6_PANGAEA/output_SOCATv6_ascii_2 -socatversion 6 -usereynolds -startyr 1957 -endyr 2017 -keepduplicates -keeptempfiles -asciioutput


