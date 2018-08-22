#!/home/cercache/tools/environments/run_python_scibox_precise.sh
from datetime import datetime
import numpy

def datenum(y,m,d,h=0,mn=0,s=0):
   """duplicates matlab datenum function."""
   res=datetime(y,m,d,h,mn,s)-datetime(1,1,1)
   #Have changed the below from +367 to +1 - not sure what the +367 was for but had the
   #effect of adding a year and a day to most of the results (depending on leap year)
   return res.days+res.seconds/(60.*60.*24)+1. # +1 to get starting from 1 rather than 0 (consistent with datetime.fromordinal)

#Make the datenum function work with numpy arrays
datenum_array=numpy.frompyfunc(datenum,6,1)



#from datetime import timedelta;
##reverses the datenum function
#def undatenum(jds):
#    baseline = datetime(1,1,1);
#    dt = timedelta(days=jds-1);
#    
#    return baseline+dt;
#
#j = datenum(2010, 5, 11, 6, 0, 0);
#d = undatenum(j);
#
#print j
#print d