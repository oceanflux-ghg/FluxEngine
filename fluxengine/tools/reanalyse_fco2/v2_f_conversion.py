#!/home/cercache/tools/environments/run_python_scibox_precise.sh
#-*- coding:utf-8 -*-

import numpy as np
from . import datenum
import tempfile
import PyCO2SYS as pyco2 # DJF 20/12/2024: Import pyCO2sys as we want to use the temperature correction procedure described in Humpreys (2024; )
                         # If were importing pyCO2sys we may as well use it for all components, and it gives us flexibility to use the Takahashi et al. approaches too.

# Constants
R = 82.0578 # cm^3 atm/(mol K)
hPa2atm = 100. * 9.867E-06
EPS = 2.2204460492503131e-16
trend = 1.5e-06 # pCO2 [atm/year] (Takahashi et al., 2009)
tempdir=tempfile.mkdtemp(prefix="ignore_points_")

def v2_f_conversion_wrap(jds,data_array,Tcls,Tcls_unc,Peq_cls,extrapolatetoyear=None,temperature_handling=1):
   """
    Wrapper function to run v2_f_conversion but from a structured array as input.
    Also returns result as a structured array.
   """
   #Run the conversion function
   jd, yr, mon, day, hh, mm, ss, lon, lat, SST_C, Tcl_C, Tcl_C_unc, fCO2_SST, fCO2_Tym_final, pCO2_SST, pCO2_Tym_final, qf = v2_f_conversion(jds, data_array['year'],data_array['month'],data_array['day'],data_array['hour'],data_array['minute'],data_array['second'],
                                                      data_array['longitude'], data_array['latitude'], data_array['sst'],data_array['salinity'], data_array['T_equ'],
                                                      data_array['air_pressure'], data_array['air_pressure_equ'], data_array['salinity_sub'],data_array['air_pressure_sub'],
                                                      data_array['fCO2'], Tcls,Tcls_unc, Peq_cls,extrapolatetoyear,temperature_handling);

   if jd is None:
      #this is only if there were no usable data after the validity checks
      return None
   else:
      #concatenate arrays into single, structured array for return
      result=np.recarray((jd.size,),dtype=[('jd',float),
                                           ('yr',np.int32),
                                           ('mon',np.int32),
                                           ('day', np.int32),
                                           ('hh', np.int32),
                                           ('mm', np.int32),
                                           ('ss', np.int32),
                                           ('lat',float),
                                           ('lon',float),
                                           ('SST_C',float),
                                           ('Tcl_C',float),
                                           ('Tcl_C_unc',float),
                                           ('fCO2_SST',float),
                                           ('fCO2_Tym',float),
                                           ('pCO2_SST',float),
                                           ('pCO2_Tym',float),
                                           ('qf',int)]);
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
      result['Tcl_C_unc'] = Tcl_C_unc
      result['fCO2_SST']=fCO2_SST
      result['fCO2_Tym']=fCO2_Tym_final
      result['pCO2_SST']=pCO2_SST
      result['pCO2_Tym']=pCO2_Tym_final
      result['qf']=qf

   return result


def v2_f_conversion(jds, yrs, mons, days, hhs, mms, sss, lons, lats, SST_Cs, sals, Teq_Cs, Ps, Peqs, sal_woas1, \
   P_nceps, fCO2_recs, Tcls, Tcls_unc, Peq_cls,extrapolatetoyear,temperature_handling):
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

   print('FluxEngine using temperature handler number: ' + str(temperature_handling) + ' - refer to pyCO2sys for details')
   #Because this function changes the values of the sal_woas array we should copy it and change the copy instead
   sal_woas=sal_woas1.copy()
   # only use records where SST_Cs, fCO2_recs and Tcls are valid
   goodpoints=np.where((np.isfinite(SST_Cs)) & (np.isfinite(fCO2_recs)) & (Tcls >= 0) & (Tcls < 1000) & (np.isfinite(Tcls)) & (np.isfinite(Peq_cls)))# some Tcl data = 9.96921e+36 were found for ATS-ARC
   badpoints=np.where(~(np.isfinite(SST_Cs)) | ~(np.isfinite(fCO2_recs)) | (Tcls <0) | (Tcls >= 1000) | ~(np.isfinite(Tcls)) | ~(np.isfinite(Peq_cls))) #TMH: updated so that good and bad points are mutually exclusive and the whole domain
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
   Tcl_unc = Tcls_unc[goodpoints]
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
   # if Teq is invalid, use SST (provided in SOCAT files)
   w = np.where(np.isnan(Teq_C))
   Teq_C[w] = np.take(SST_C,w) # If Teq isn't avaiable we make this the same as SST_C (so we don't recorrect back to Teq as its not provided, and then go forwards from here).
   f2 = np.zeros(n)
   f2[w] = 1
   # if P is invalid, use NCEP/NCAR provided in SOCAT fiels
   w = np.where(np.isnan(P))
   P[w] = np.take(P_ncep,w)
   f3 = np.zeros(n)
   f3[w] = 1
   # if Peq is invalid, add 3 hPa to P (Due to overpressure maintained generally in ships)
   w = np.where(np.isnan(Peq))
   Peq[w] = np.take(P,w) + 3
   f4 = np.zeros(n)
   f4[w] = 1
   qf = f1 + f2 + f3 + f4
   P *= hPa2atm # [atm]
   Peq *= hPa2atm # [atm]
   Peq_cl *= hPa2atm  # [atm]

   # DJF 20/12/2024: Don't need to do these conversions as this will be handled by PyCO2sys
   # fCO2_SST = fCO2_rec * 1E-06 # [atm]
   # Recalculation of original measurements
   # SST = SST_C + 273.15 # [deg C]
   # delta = 57.7 - 0.118 * SST # cm^3/mol
   # B = -1636.75 + SST * (12.0408 + SST * (-3.27957E-02 + 3.16528E-05*SST)) # cm^3/mol

   # DJF 20/12/2024: This conversion of fCO2 back to Teq and then recalculated forward to fCO2sw will be done with CO2sys
#    pCO2_SST = fCO2_SST.copy() # initial first guess of pCO2_SST
#    dT = SST_C - Teq_C
#    y = [0]
#    while np.any(np.absolute(pCO2_SST - y) > EPS * pCO2_SST):
# #      print "mean(y)", np.mean(y);
# #      print "y:", y[0];
# #      print "DT:", dT[0];
# #      print "pCO2_SST", pCO2_SST[0];
# #      print "Peq:", Peq[0];
# #      print "fCO2_SST:", fCO2_SST[0];
# #      print "SST:", SST[0];
# #      print "B", B[0];
# #      print "delta:", delta[0];
# #      print "R:", R; # cm^3 atm/(mol K)
#
#
#       y = pCO2_SST
#       pCO2_Teq = pCO2_SST * np.exp(-0.0423 * dT)
#       XCO2_Teq = pCO2_Teq / Peq # wet XCO2_Teq
#       pCO2_SST = fCO2_SST * np.exp(-(B + 2 * delta * (1 - XCO2_Teq) ** 2) * Peq / (R * SST))
#
# #      raw_input("modified...");

   # DJF 20/12/2024: We assume here that SOCAT data providers will have converted their fCO2(sw) at SST from the Teq
   # using the Takahashi approach as the newer method is only recent (i.e 2024). So we correct the fCO2(sw) back to
   # Teq using Takahashi then go forwards to the subskin temp (from satellite etc) using the method defined.
   # Either Takahashi linear, Takahashi quadratic or the Humpreys 2024 apporahc (new default)

   pyco2_equil = pyco2.sys(
            par1 = fCO2_rec, # SOCAT recommended fCO2(sw)
            par1_type = 5, # Specifying that above input is fCO2(sw) in uatm (5) - could also be pCO2sw in uatm (4) or xCO2sw in ppm (9)
            par2 = None, # Specifying to CO2sys that we're not using another parameter (so not solving carbonate system) just doing conversions
            par2_type = None,
            salinity = sal, # Provides the salinity, or the WOA salinity or a fixed value (35.0) - salinity not needed for this correction
            temperature = SST_C, # The temperature that the fCO2sw above is at.
            temperature_out = Teq_C, # The temperature that we want the output to be at (i.e the Teq)
            pressure_atmosphere = P, # This will be the sea level pressure
            pressure_atmosphere_out = Peq, # And we want the output to be at the equlibrator pressure
            opt_adjust_temperature = 5 # This sets pyCO2sys to use the Takahashi et al. (1993) linear temperature correction -
            # THis is the approach used by SOCAT to get to fCO2rec (Bakker et al. 2016) - so we are coverting back to then move forwards with the method selected in the fucntion.
   )

   pyco2_subskin = pyco2.sys(
            par1 = pyco2_equil['pCO2_out'], #SOCAT fCO2swrec corrected back to equlibrator, and output as pCO2sw
            par1_type = 4, # Giving as pCO2sw
            par2 = None,# Specifying to CO2sys that we're not using another parameter (so not solving carbonate system) just doing conversions
            par2_type = None,
            salinity = sal, # Provides the salinity, or the WOA salinity or a fixed value (35.0) - salinity maybe needed
            temperature = Teq_C,# The temperature that the data is at in the equilibrator (i.e the Teq)
            temperature_out = Tcl_C,# The temperature that we want the output to be at (i.e the satellite subskin data)
            pressure_atmosphere = Peq, #Pressure at equlibrator
            pressure_atmosphere_out = P, # Pressure at sea level
            opt_adjust_temperature = temperature_handling #
   )

   pCO2_Tym = pyco2_subskin['pCO2_out']

   #DJF: This is no longer needed and covered by pyco2sys calls above
   # pCO2_Teq = pCO2_SST * np.exp(-0.0423 * dT)
   # # Recalculation to climatological values
   # pCO2_Tym = pCO2_Teq * np.exp(0.0433 * (Tcl_C - Teq_C) - \
   #    4.35E-05 * (Tcl_C ** 2 - Teq_C ** 2))

   if extrapolatetoyear is not None:
      # extrapolate pCO2_Tym to given year using the Takahashi et al 2009 trend
      #TODO FIXME this should be updated to take into account leap years
      dt = (datenum.datenum(extrapolatetoyear,1,1,0,0,0)-jd)/365.0
      pCO2_Tym_final = pCO2_Tym + trend*dt
   else:
      # rename the variable so that the succeeding script will work and continue
      # would be better not doing this (from intuitive point of view) but easiest way.
      pCO2_Tym_final = pCO2_Tym

   # Convert from subskin pCO2sw to subskin fCO2sw
   pyco2_subskin_fco2 = pyco2.sys(
            par1 = pCO2_Tym_final, #pCO2sw corrected to the subskin, with extrapolation to year above applied
            par1_type = 4, # Giving as pCO2sw
            par2 = None,# Specifying to CO2sys that we're not using another parameter (so not solving carbonate system) just doing conversions
            par2_type = None,
            salinity = sal, # Provides the salinity, or the WOA salinity or a fixed value (35.0) - salinity maybe needed
            temperature = Tcl_C,# The temperature that the data is at in the equilibrator (i.e the Teq)
            pressure_atmosphere = P, #Pressure at equlibrator
   )

   fCO2_Tym_final = pyco2_subskin_fco2['fCO2']
   fCO2_SST = fCO2_rec

   # Converting fCO2swrec to pCO2swrec
   pyco2_pCO2_SST = pyco2.sys(
            par1 = fCO2_SST, # SOCAT recommended fCO2(sw)
            par1_type = 5, # Specifying that above input is fCO2(sw) in uatm (5) - could also be pCO2sw in uatm (4) or xCO2sw in ppm (9)
            par2 = None, # Specifying to CO2sys that we're not using another parameter (so not solving carbonate system) just doing conversions
            par2_type = None,
            salinity = sal, # Provides the salinity, or the WOA salinity or a fixed value (35.0) - salinity not needed for this correction
            temperature = SST_C, # The temperature that the fCO2sw above is at.
            pressure_atmosphere = P, # This will be the sea level pressure
   )

   pCO2_SST = pyco2_pCO2_SST['pCO2']

   #DJF 20/12/2024: These are no longer required as all components are calculated with the pyCO2sys calls above - and everything is output as uatm instead of atm
   # delta = 57.7 - .118 * Tcl # cm^3/mol
   # B = -1636.75 + Tcl * (12.0408 + Tcl * (-3.27957E-02 + 3.16528E-05 * Tcl)) # cm^3/mol
   # exponent=np.exp((B + 2. * delta * (1 - pCO2_Teq / Peq_cl) ** 2) * Peq_cl / (R * Tcl))
   # fCO2_Tym_final = pCO2_Tym_final * exponent

   # # conversion from atm to uatm
   # fCO2_Tym_final *= 1E+06 # uatm
   # pCO2_Tym_final *= 1E+06 # uatm
   # fCO2_SST *= 1E+06 # uatm (this is the same as fCO2_rec)
   # pCO2_SST *= 1E+06 # uatm

   return [jd, yr, mon, day, hh, mm, ss, lon, lat, SST_C, Tcl_C, Tcl_unc, fCO2_SST, fCO2_Tym_final, pCO2_SST, pCO2_Tym_final, qf];
