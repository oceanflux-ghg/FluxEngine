#!/home/cercache/tools/environments/run_python_scibox_precise.sh
#-*- coding:utf-8 -*-

import numpy as np
from . import datenum
import tempfile

# Constants
R = 82.0578 # cm^3 atm/(mol K)
hPa2atm = 100. * 9.867E-06
EPS = 2.2204460492503131e-16
trend = 1.5e-06 # pCO2 [atm/year] (Takahashi et al., 2009) 
tempdir=tempfile.mkdtemp(prefix="ignore_points_")

def v2_f_conversion_wrap(jds,data_array,Tcls,Peq_cls,extrapolatetoyear=None):
   """
    Wrapper function to run v2_f_conversion but from a structured array as input. 
    Also returns result as a structured array.
   """
   #Run the conversion function
   jd, yr, mon, day, hh, mm, ss, lon, lat, SST_C, Tcl_C, fCO2_SST, fCO2_Tym_final, pCO2_SST, pCO2_Tym_final, qf = v2_f_conversion(jds, data_array['year'],data_array['month'],data_array['day'],data_array['hour'],data_array['minute'],data_array['second'],
                                                      data_array['longitude'], data_array['latitude'], data_array['sst'],data_array['salinity'], data_array['T_equ'], 
                                                      data_array['air_pressure'], data_array['air_pressure_equ'], data_array['salinity_sub'],data_array['air_pressure_sub'], 
                                                      data_array['fCO2'], Tcls, Peq_cls,extrapolatetoyear);

   if jd is None:
      #this is only if there were no usable data after the validity checks
      return None
   else:
      #concatenate arrays into single, structured array for return
      result=np.recarray((jd.size,),dtype=[('jd',np.float),
                                           ('yr',np.int32),
                                           ('mon',np.int32),
                                           ('day', np.int32),
                                           ('hh', np.int32),
                                           ('mm', np.int32),
                                           ('ss', np.int32),
                                           ('lat',np.float),
                                           ('lon',np.float),
                                           ('SST_C',np.float),
                                           ('Tcl_C',np.float),
                                           ('fCO2_SST',np.float),
                                           ('fCO2_Tym',np.float),
                                           ('pCO2_SST',np.float),
                                           ('pCO2_Tym',np.float),
                                           ('qf',np.int)]);
      result['jd']=jd
      result['yr']=yr;
      result['mon']=mon;
      result['day']=day;
      result['hh']=hh;
      result['mm']=mm;
      result['ss']=ss;
      result['lat']=lat
      result['lon']=lon
      result['SST_C']=SST_C
      result['Tcl_C']=Tcl_C
      result['fCO2_SST']=fCO2_SST
      result['fCO2_Tym']=fCO2_Tym_final
      result['pCO2_SST']=pCO2_SST
      result['pCO2_Tym']=pCO2_Tym_final
      result['qf']=qf

   return result


def v2_f_conversion(jds, yrs, mons, days, hhs, mms, sss, lons, lats, SST_Cs, sals, Teq_Cs, Ps, Peqs, sal_woas1, \
   P_nceps, fCO2_recs, Tcls, Peq_cls,extrapolatetoyear):
   """Recalculates CO2 flux from the ocean:
      Arguments (all numpy arrays):
        jds - days since 0/1/0000 (from datenum.py)
        lons - longitudes
        lats - latitudes
        SST_Cs - sea surface temperatures (temp) [deg C]
        sals - salinities (salinity) [PSU]
        Teq_Cs - water temperature at equilibrator (temperature_eyui) [deg C]
        Ps - atmospheric pressure (Pressure_atm) [hPa]
        Peqs - pressure at the equilibrator (Pressure_equi) [hPa]
        sal_woas - salinity estimated from World Ocean Atlas 2005 (woa_sss) [PSU]
        P_nceps - sea level pressure(hPa) estimated from NCEP/NCAR (ncep_slp) [hPa]
        fCO2_recs - fCO2 "recomputed" (fCO2_rec) [uatm] (-999 is invalid data)
        Tcls - climatological sea surface sub skin temperature [K] to convert to
        Peq_cls - pressure at the equilibrator estimated from climatological atmospheric pressure [hPa] (use P_ncep + 3 hPa if not available)
        extrapolatetoyear - The year we want to extrapolate to using the Takahashi trend
      Inputs come from Surface Ocean CO2 Atlas (SOCAT) version 1.5.
      Return value: list of numpy arrays (equivalent inputs are overwritten):
        jd - jds for which valid data exist
        lon - lons for which valid data exist
        lat - lats for which valid data exist
        SST_C - SST_Cs for which valid data exist
        Tcl_C - Tcls for which valid data exist [deg C]
        fCO2_SST - fCO2 recomputed by SOCAT for SST_C (uatm)
        fCO2_Tym - fCO2 recomputed for Tcl_C (uatm)
        qf - quality flag?"""

   
   #Because this function changes the values of the sal_woas array we should copy it and change the copy instead
   sal_woas=sal_woas1.copy()
   # only use records where SST_Cs, fCO2_recs and Tcls are valid 
   goodpoints=np.where((np.isfinite(SST_Cs)) & (np.isfinite(fCO2_recs)) & (Tcls >= 0) & (Tcls < 1000) & (np.isfinite(Peq_cls)))# some Tcl data = 9.96921e+36 were found for ATS-ARC
   badpoints=np.where(~(np.isfinite(SST_Cs)) | ~(np.isfinite(fCO2_recs)) | (Tcls <0) | (Tcls >= 1000) | ~(np.isfinite(Peq_cls))) #TMH: updated so that good and bad points are mutually exclusive and the whole domain
#   aa = len(np.where(~(np.isfinite(SST_Cs)))[0]);
#   bb = len(np.where(~(np.isfinite(fCO2_recs)))[0]);
#   cc = len(np.where((Tcls <0) | (Tcls >= 1000) )[0]);
#   dd = len(np.where(~(np.isfinite(Peq_cls)) )[0]);
   
   
   if badpoints[0].size != 0:
      fileout=tempfile.mkstemp(prefix="%s/ignored_points_"%tempdir)[1]
      print("Writing ignored points to temp file: %s"%fileout)
      fout=open(fileout,'w')
      fout.write("jd,lat,lon,SST_C,fCO2_rec,Tcl,Peq_cl\n")
      for item in badpoints[0]:
         fout.write("%s,%s,%s,%s,%s,%s,%s\n"%(jds[item],lats[item],lons[item],SST_Cs[item],fCO2_recs[item],Tcls[item],Peq_cls[item]))
   if goodpoints[0].size == 0:
      #there are no records with valid  Tcls, fCO2_recs and SST_cs
      #raise Exception("No data records with valid Tcls, fCO2_recs and SST_Cs. Cannot reanalyse these data.")
      #return a list of Nones of the length that needs to be unpacked
      return [None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None];

   jd, lon, lat, SST_C, sal = jds[goodpoints], lons[goodpoints], lats[goodpoints], SST_Cs[goodpoints], sals[goodpoints]
   yr, mon, day, hh, mm, ss = yrs[goodpoints], mons[goodpoints], days[goodpoints], hhs[goodpoints], mms[goodpoints], sss[goodpoints];
   Teq_C, P, Peq, sal_woa = Teq_Cs[goodpoints], Ps[goodpoints], Peqs[goodpoints], sal_woas[goodpoints]
   P_ncep, fCO2_rec, Tcl, Peq_cl = P_nceps[goodpoints], fCO2_recs[goodpoints], Tcls[goodpoints], Peq_cls[goodpoints]
   n = np.size(jd)
   Tcl_C  = Tcl - 273.15
   # if sal_woas is invalid, use 35.
   w = np.where(np.isnan(sal_woas))
   sal_woas[w] = 35.
   # if sal is invalid, use WOA 2005
   w = np.where(np.isnan(sal))
   sal[w] = np.take(sal_woa,w)
   f1 = np.zeros(n)
   f1[w] = 1
   # if Teq is invalid, use SST
   w = np.where(np.isnan(Teq_C))
   Teq_C[w] = np.take(SST_C,w)
   f2 = np.zeros(n)
   f2[w] = 1
   # if P is invalid, use NCEP/NCAR
   w = np.where(np.isnan(P))
   P[w] = np.take(P_ncep,w)
   f3 = np.zeros(n)
   f3[w] = 1
   # if Peq is invalid, add 3 hPa to P
   w = np.where(np.isnan(Peq))
   Peq[w] = np.take(P,w) + 3
   f4 = np.zeros(n)
   f4[w] = 1
   qf = f1 + f2 + f3 + f4
   P *= hPa2atm # [atm]
   Peq *= hPa2atm # [atm]
   Peq_cl *= hPa2atm  # [atm]
   fCO2_SST = fCO2_rec * 1E-06 # [atm]
   # Recalculation of original measurements
   SST = SST_C + 273.15 # [deg C]
   delta = 57.7 - 0.118 * SST # cm^3/mol
   B = -1636.75 + SST * (12.0408 + SST * (-3.27957E-02 + 3.16528E-05*SST)) # cm^3/mol
   pCO2_SST = fCO2_SST.copy() # initial first guess of pCO2_SST
   dT = SST_C - Teq_C
   y = [0]
   while np.any(np.absolute(pCO2_SST - y) > EPS * pCO2_SST):
#      print "mean(y)", np.mean(y);
#      print "y:", y[0];
#      print "DT:", dT[0];
#      print "pCO2_SST", pCO2_SST[0];
#      print "Peq:", Peq[0];
#      print "fCO2_SST:", fCO2_SST[0];
#      print "SST:", SST[0];
#      print "B", B[0];
#      print "delta:", delta[0];
#      print "R:", R; # cm^3 atm/(mol K)
      
      
      y = pCO2_SST
      pCO2_Teq = pCO2_SST * np.exp(-0.0423 * dT)
      XCO2_Teq = pCO2_Teq / Peq # wet XCO2_Teq
      pCO2_SST = fCO2_SST * np.exp(-(B + 2 * delta * (1 - XCO2_Teq) ** 2) * Peq / (R * SST))
      
#      raw_input("modified...");

   pCO2_Teq = pCO2_SST * np.exp(-0.0423 * dT)
   # Recalculation to climatological values
   pCO2_Tym = pCO2_Teq * np.exp(0.0433 * (Tcl_C - Teq_C) - \
      4.35E-05 * (Tcl_C ** 2 - Teq_C ** 2))

   if extrapolatetoyear is not None:
      # extrapolate pCO2_Tym to given year using the Takahashi et al 2009 trend
      #TODO FIXME this should be updated to take into account leap years
      dt = (datenum.datenum(extrapolatetoyear,1,1,0,0,0)-jd)/365.0
      pCO2_Tym_final = pCO2_Tym + trend*dt
   else:
      # rename the variable so that the succeeding script will work and continue 
      # would be better not doing this (from intuitive point of view) but easiest way.
      pCO2_Tym_final = pCO2_Tym

   # Convert fro pCO2 to fCO2
   delta = 57.7 - .118 * Tcl # cm^3/mol
   B = -1636.75 + Tcl * (12.0408 + Tcl * (-3.27957E-02 + 3.16528E-05 * Tcl)) # cm^3/mol
   exponent=np.exp((B + 2. * delta * (1 - pCO2_Teq / Peq_cl) ** 2) * Peq_cl / (R * Tcl))
   fCO2_Tym_final = pCO2_Tym_final * exponent

   # conversion from atm to uatm
   fCO2_Tym_final *= 1E+06 # uatm 
   pCO2_Tym_final *= 1E+06 # uatm 
   fCO2_SST *= 1E+06 # uatm (this is the same as fCO2_rec)
   pCO2_SST *= 1E+06 # uatm

   return [jd, yr, mon, day, hh, mm, ss, lon, lat, SST_C, Tcl_C, fCO2_SST, fCO2_Tym_final, pCO2_SST, pCO2_Tym_final, qf];

