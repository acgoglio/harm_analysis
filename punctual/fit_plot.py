#
# Script for HARMONIC ANALYSIS AND POST_PROC
#
# imports
import matplotlib.pyplot as plt
import matplotlib as mpl # Palettes
import numpy as np
import netCDF4 as NC
import os
import sys
import warnings
#
warnings.filterwarnings("ignore") # Avoid warnings
#
from scipy.optimize import curve_fit
from scipy import stats
import collections
import pandas as pd
import csv
import math
import datetime
from datetime import datetime
from operator import itemgetter 
import plotly
from plotly import graph_objects as go # for bar plot
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LogNorm
from operator import itemgetter # to order lists
from statsmodels.distributions.empirical_distribution import ECDF # empirical distribution functions
#
# Import ttide code for Forman harmonic analysis
import ttide
# Import TPXO and Literature values
from lit_tpxo import *
#
#
# by AC Goglio (CMCC)
# annachiara.goglio@cmcc.it
#
# Written: 06/2019
# Last Modified: 14/10/2021
#
# Script to fit and plot the harmonic analisi reults wrt tide gauges data, tpxo model and literature values
#
#################################################################
# The user should modify the following lines to set his run
#################################################################
# General run parameters:
#---------------------
# Work dir path:
# WARNING: the inputs must be here, the outputs will be moved to subdirs   
workdir= '/work/oda/ag15419/tmp/tides8_v31_newHA/' 
# input files:
emodnettg_coo_file = '/users_home/oda/ag15419/harm_analysis/punctual/emodnet_TGb_newTGb_all.coo'
model_bathy='/work/oda/ag15419/PHYSW24_DATA/TIDES/DATA0/bathy_meter.nc'
#

# Domain (Med or AtlBox)
where_box='AtlBox'

# Option for phase plots
cos_pha = 0 
# =1 to compare cosine values of phases instead of abs values (to avoid +-360*n corrections..)

########################################################
# DO NOT CHANGE THE CODE BELOW THIS LINE!!!
########################################################
# Parameter setting
#--------------------
# MODEL DATASET
# WARNING: this must be the same as in p_extr.ini file (var ANA_INTAG)
mod_file_template='Tides8_v31' #'Tides8_v31' 'eas6'

# Fields to be analized
grid = 'T' # Choose T, U, V or uv2t grid

if grid == 'T':
   # 2D fields:
   var_2d_mod='sossheig'
   field_2d_units='m'
   # WARNING: the sossheig field is supposed to be given in meters 
   # You can provide two different name for the var because of different name of time attributes 
   # for different output periods
   time_var_mod='time'
   time_var_mod2='time_counter'
   lat_var_mod='lat'

#--------------------
# OBS DATASET

if grid == 'T':
   var_2d_obs=var_2d_mod # Only for EMODnet dataset
   time_var_obs='time' # Only for new EMODnet dataset
   time_var_obs2='TIME' # Old EMODnet dataset (e.g. IOC_* files )
   fieldobs_2d_units=('cm','m')
   # WARNING: the sossheig field is supposed to be given 
   # in centimeters for ispra and in meters for emodnet  
   var_2d_qf_obs='sossheig_qf'
   pos_2d_qf_obs='pos_qf'

#--------------------
# OTHER PARAM

# Tidal components (WARNING: the script is set to work with the following 8 constituents in the following order!)
tidal_comp=['M2','S2','K1','O1','N2','P1','Q1','K2']

# Colors for each sub area

# MED subregions are: Gibraltar, Adriatic, Messina, East, Other
# WARNING: you cannot change the region order!!
subregions_color=['red','blue','orange','magenta','green','cyan','deeppink','tab:olive']
subregions_labels=['G','A','M','E','O','TA','TG','T']
# Atl Box subregions are: Biscay Bay, Gibraltar Atlantic side, Portugal  
# WARNING: you cannot change the region order!!
#subregions_color.append['cyan','deeppink','tab:olive']
#
# Function to distinghish relevant regions
def which_region(longitude,latitude):
       # Regions 
       # EAST MED:
       if longitude > 30.000:
           color=subregions_color[3]
       # ADRIACTIC SEA:
       elif (longitude < 20.000 and longitude > 12.000 and latitude < 46.000 and latitude > 42.000 ) or (longitude < 20.000 and longitude > 14.000 and latitude < 42.000 and latitude > 41.000 ) or (longitude < 20.000 and longitude > 16.000 and latitude < 42.000 and latitude > 40.000):
           color=subregions_color[1]
       # MESSINA STRAIT AREA:
       elif longitude < 16.500 and longitude > 14.500 and latitude < 38.200 and latitude > 37.800:
           color=subregions_color[2]
       # GIBRALTAR STRAIT AREA:
       elif longitude < -1.800 and longitude > -6.000 and latitude < 37.2630:
           color=subregions_color[0]
       # GIBRALTAR ATLANTIC BOX SIDE
       elif longitude < -6.000 and longitude > -10.000 and latitude < 37.2630:
           color=subregions_color[6]
       # PORUGAL ATLANTIC BOX
       elif (longitude < -6.000 and latitude > 37.2630 and latitude < 43.2100) or ( longitude < -10.000 ):
           color=subregions_color[7]
       # BISCAY BAY ATLANTIC BOX
       elif longitude < 0.000 and latitude > 43.2100:
           color=subregions_color[5]
       # OTHER MED REGIONS
       else:
           color=subregions_color[4]
       return color


# Function to divide Med basin from Atlantic Box
# +------+--+-+-----+-+
# |      |AB| |     | |
# |  AB  +--+-+ MED +-+
# |      |MED |       |
# +------+----+-------+
#
def which_domain(longitude,latitude):
    # Check if the TG is inside the whole domain
    if longitude < 36.29 and longitude > -18.12 and latitude < 45.98 and latitude > 30.18:
       if longitude > 0.000:
          domain='Med'
       elif longitude < -6.000:
          domain='AtlBox'
       else:
          if latitude < 42.000:
             domain='Med'
          else:
             if longitude < -2.000:
                domain='AtlBox'
             else:
                domain='OUT' 
    else:
       domain='OUT'

    # Remove Black Sea
    if longitude < 36.29 and longitude > 27.20 and latitude > 40.00 and latitude < 45.98:
       domain='OUT'

    return domain

# Signal-to-noise ration threshold
# WARNING: if any component shows a fit snr lower than this the fit is not perfomed
snr_thershold=0.1

# OPTIONS ON NAMES AND ERROR BARS
# Flag to avoid TG names/numbers in the lin reg plots (set linreg_name_flag = 0 to avoid the strings..)
linreg_name_flag = 1
# To avoid fit error bars in Amp, Pha and lin reg plots set errbar_flag = 0 
errbar_flag = 1

# Option to invert the order in the time-series plot
# =0 MOD/OBS (for Med); !=0 OBS/MOD (for AtlBox)
if where_box == 'Med':
   revts_flag = 1
else:
   revts_flag = 0 

# Option to use the full names of the TGs instead of numbers for lit analysis type only (the other use numbers) 
# = 0 to use numbers; = 1 to use full names  
lit_fullname = 0

# FLAG or loop on analysis type
# Options:
# lit       --> Compare the common datasets with respect to literature 
# anatpxo   --> Apply all the analysis and compare datasets with TPXO model results
# all       --> Linear regression concerning all avalilable tide-gauges 
for anatype_flag in ('all','anatpxo','lit'): #'all','lit','anatpxo'

   # Buil the dir and move in it
   workdir_path = workdir+'/'+anatype_flag+'_'+where_box+'/'
   try:
       os.stat(workdir_path)
   except:
       os.mkdir(workdir_path)

   print ('##############################')
   # for each type of analysis set the parameters and print the function
   if anatype_flag=='anatpxo':
      tpxo_flag = 1 # to compare also wrt TPXO data
      flag_15stats = 0 # to compare results with literature
      print ('Comparison wrt obs and TPXO model results.. Results in ',workdir_path)
      fontsize_tg=40
   elif anatype_flag=='lit':
      tpxo_flag = 0 # to compare also wrt TPXO data
      flag_15stats = 1 # to compare results with literature
      print ('Comparison wrt obs and literature results.. Results in ',workdir_path)
      fontsize_tg=40
   elif anatype_flag=='all':
      tpxo_flag = 0 # to compare also wrt TPXO data
      flag_15stats = 0 # to compare results with literature
      print ('Comparison wrt all available obs.. Results in ',workdir_path)
      if where_box=='Med':
         fontsize_tg=20
      elif where_box=='AtlBox':
         fontsize_tg=40

   # Check on Domain and fix the bdy
   if where_box=='AtlBox':
      if anatype_flag == 'all':
         tpxo_flag = 0
         flag_15stats = 0
      elif anatype_flag=='anatpxo':
         tpxo_flag = 1
         flag_15stats = 0
      else:
         exit(0)

   
   # ======== EMODnet TG ========== 
   
   # STZ INFO
   
   # Read tide-gauges infos from .coo file
   
   tg_name=[]  # Name of the TG
   tg_col=[]   # Color corresponding to the subregion belonging of the TG
   tg_lon=[]   # Longitude of the TG
   tg_lat=[]   # Latitude of the TG
   tg_sdate=[] # Start date of the obs time series

   # Open file and get values
   fT1_coo = pd.read_csv(emodnettg_coo_file,sep=';',comment='#',header=None)
   #
   tg_inlat = fT1_coo[0][:] 
   tg_inlon = fT1_coo[1][:] 
   tg_inname = fT1_coo[2][:] 
   anatpxo_inflag = fT1_coo[5][:] 
   lit_inflag = fT1_coo[6][:] 
 
   if where_box == 'Med':
  
      if anatype_flag=='lit':
         for idx_tg in range (0,len(lit_inflag)):
             if lit_inflag[idx_tg] == 1 :
                if which_domain(tg_inlon[idx_tg],tg_inlat[idx_tg]) == 'Med':
                   tg_name.append(tg_inname[idx_tg])
                   tg_lon.append(tg_inlon[idx_tg])
                   tg_lat.append(tg_inlat[idx_tg])
                   tg_col.append(which_region(tg_inlon[idx_tg],tg_inlat[idx_tg]))

      elif anatype_flag=='anatpxo':
         for idx_tg in range (0,len(anatpxo_inflag)):
             if anatpxo_inflag[idx_tg] == 1 :
                if which_domain(tg_inlon[idx_tg],tg_inlat[idx_tg]) == 'Med':
                   tg_name.append(tg_inname[idx_tg])
                   tg_lon.append(tg_inlon[idx_tg])
                   tg_lat.append(tg_inlat[idx_tg])
                   tg_col.append(which_region(tg_inlon[idx_tg],tg_inlat[idx_tg]))
   
      elif anatype_flag=='all':
         
         for idx_tg in range (0,len(tg_inname)):
            if which_domain(tg_inlon[idx_tg],tg_inlat[idx_tg]) == 'Med':
               tg_name.append(tg_inname[idx_tg])
               tg_lon.append(tg_inlon[idx_tg])
               tg_lat.append(tg_inlat[idx_tg])
               tg_col.append(which_region(tg_inlon[idx_tg],tg_inlat[idx_tg]))
   
   elif where_box == 'AtlBox':

       if anatype_flag=='all':
         for idx_tg in range (0,len(tg_inname)):
            if which_domain(tg_inlon[idx_tg],tg_inlat[idx_tg]) == 'AtlBox':
               tg_name.append(tg_inname[idx_tg])
               tg_lon.append(tg_inlon[idx_tg])
               tg_lat.append(tg_inlat[idx_tg])
               tg_col.append(which_region(tg_inlon[idx_tg],tg_inlat[idx_tg]))
   
       elif anatype_flag=='anatpxo':
         for idx_tg in range (0,len(anatpxo_inflag)):
             if anatpxo_inflag[idx_tg] == 1 :
                if which_domain(tg_inlon[idx_tg],tg_inlat[idx_tg]) == 'AtlBox':
                   tg_name.append(tg_inname[idx_tg])
                   tg_lon.append(tg_inlon[idx_tg])
                   tg_lat.append(tg_inlat[idx_tg])
                   tg_col.append(which_region(tg_inlon[idx_tg],tg_inlat[idx_tg]))
   
   tg_num=len(tg_name)
   print("EmodNET Stations num:",tg_num)
   
   # EMODNET FIT by means of ttide code
   
   # Allocate space for arrays
   # Amplitudes and Phases for model 
   M2_amp_mod=['Nan']*len(tg_name); M2_pha_mod=['Nan']*len(tg_name); K1_amp_mod=['Nan']*len(tg_name); K1_pha_mod=['Nan']*len(tg_name)
   O1_amp_mod=['Nan']*len(tg_name); O1_pha_mod=['Nan']*len(tg_name); S2_amp_mod=['Nan']*len(tg_name); S2_pha_mod=['Nan']*len(tg_name)
   P1_amp_mod=['Nan']*len(tg_name); P1_pha_mod=['Nan']*len(tg_name); N2_amp_mod=['Nan']*len(tg_name); N2_pha_mod=['Nan']*len(tg_name)
   Q1_amp_mod=['Nan']*len(tg_name); Q1_pha_mod=['Nan']*len(tg_name); K2_amp_mod=['Nan']*len(tg_name); K2_pha_mod=['Nan']*len(tg_name)
   # Amplitudes and Phases for obs
   M2_amp_obs=['Nan']*len(tg_name); M2_pha_obs=['Nan']*len(tg_name); K1_amp_obs=['Nan']*len(tg_name); K1_pha_obs=['Nan']*len(tg_name)
   O1_amp_obs=['Nan']*len(tg_name); O1_pha_obs=['Nan']*len(tg_name); S2_amp_obs=['Nan']*len(tg_name); S2_pha_obs=['Nan']*len(tg_name)
   P1_amp_obs=['Nan']*len(tg_name); P1_pha_obs=['Nan']*len(tg_name); N2_amp_obs=['Nan']*len(tg_name); N2_pha_obs=['Nan']*len(tg_name)
   Q1_amp_obs=['Nan']*len(tg_name); Q1_pha_obs=['Nan']*len(tg_name); K2_amp_obs=['Nan']*len(tg_name); K2_pha_obs=['Nan']*len(tg_name)
   # Signal to noise values for the fit of mod
   M2_snr_mod=['Nan']*len(tg_name); S2_snr_mod=['Nan']*len(tg_name); K1_snr_mod=['Nan']*len(tg_name); O1_snr_mod=['Nan']*len(tg_name);
   N2_snr_mod=['Nan']*len(tg_name); P1_snr_mod=['Nan']*len(tg_name); Q1_snr_mod=['Nan']*len(tg_name); K2_snr_mod=['Nan']*len(tg_name);
   # Signal to noise values for the fit of obs
   M2_snr_obs=['Nan']*len(tg_name); S2_snr_obs=['Nan']*len(tg_name); K1_snr_obs=['Nan']*len(tg_name); O1_snr_obs=['Nan']*len(tg_name);
   N2_snr_obs=['Nan']*len(tg_name); P1_snr_obs=['Nan']*len(tg_name); Q1_snr_obs=['Nan']*len(tg_name); K2_snr_obs=['Nan']*len(tg_name);
   # Amplitude and phase errors for mod 
   M2_aerr_mod=['Nan']*len(tg_name); M2_perr_mod=['Nan']*len(tg_name); K1_aerr_mod=['Nan']*len(tg_name); K1_perr_mod=['Nan']*len(tg_name)
   O1_aerr_mod=['Nan']*len(tg_name); O1_perr_mod=['Nan']*len(tg_name); S2_aerr_mod=['Nan']*len(tg_name); S2_perr_mod=['Nan']*len(tg_name)
   P1_aerr_mod=['Nan']*len(tg_name); P1_perr_mod=['Nan']*len(tg_name); N2_aerr_mod=['Nan']*len(tg_name); N2_perr_mod=['Nan']*len(tg_name)
   Q1_aerr_mod=['Nan']*len(tg_name); Q1_perr_mod=['Nan']*len(tg_name); K2_aerr_mod=['Nan']*len(tg_name); K2_perr_mod=['Nan']*len(tg_name)
   # Amplitude and phase errors for obs
   M2_aerr_obs=['Nan']*len(tg_name); M2_perr_obs=['Nan']*len(tg_name); K1_aerr_obs=['Nan']*len(tg_name); K1_perr_obs=['Nan']*len(tg_name)
   O1_aerr_obs=['Nan']*len(tg_name); O1_perr_obs=['Nan']*len(tg_name); S2_aerr_obs=['Nan']*len(tg_name); S2_perr_obs=['Nan']*len(tg_name)
   P1_aerr_obs=['Nan']*len(tg_name); P1_perr_obs=['Nan']*len(tg_name); N2_aerr_obs=['Nan']*len(tg_name); N2_perr_obs=['Nan']*len(tg_name)
   Q1_aerr_obs=['Nan']*len(tg_name); Q1_perr_obs=['Nan']*len(tg_name); K2_aerr_obs=['Nan']*len(tg_name); K2_perr_obs=['Nan']*len(tg_name)
   
   # For salish sea method
   M2phazero=['Nan']*len(tg_name); S2phazero=['Nan']*len(tg_name); K1phazero=['Nan']*len(tg_name); O1phazero=['Nan']*len(tg_name);
   N2phazero=['Nan']*len(tg_name); P1phazero=['Nan']*len(tg_name); Q1phazero=['Nan']*len(tg_name); K2phazero=['Nan']*len(tg_name);
   M2ampnc=['Nan']*len(tg_name); S2ampnc=['Nan']*len(tg_name); K1ampnc=['Nan']*len(tg_name); O1ampnc=['Nan']*len(tg_name);
   N2ampnc=['Nan']*len(tg_name); P1ampnc=['Nan']*len(tg_name); Q1ampnc=['Nan']*len(tg_name); K2ampnc=['Nan']*len(tg_name);
   M2phanc=['Nan']*len(tg_name); S2phanc=['Nan']*len(tg_name); K1phanc=['Nan']*len(tg_name); O1phanc=['Nan']*len(tg_name);
   N2phanc=['Nan']*len(tg_name); P1phanc=['Nan']*len(tg_name); Q1phanc=['Nan']*len(tg_name); K2phanc=['Nan']*len(tg_name);
   ss_amp=[[ 0 for i in range(0,len(tidal_comp))] for j in range(0,len(tg_name))]
   ss_pha=[[ 0 for i in range(0,len(tidal_comp))] for j in range(0,len(tg_name))]
   ssmod_amp=[[ 0 for i in range(0,len(tidal_comp))] for j in range(0,len(tg_name))]
   ssmod_pha=[[ 0 for i in range(0,len(tidal_comp))] for j in range(0,len(tg_name))]
  
   # Loop on EMODnet tide-gauges
   for stn in range (0,len(tg_name)): 
       print(tg_name[stn])
   
       # Read mod and obs 
       fT1_mod = NC.Dataset(workdir+tg_name[stn]+'_mod_'+mod_file_template+'.nc','r') 
       fT1_obs = NC.Dataset(workdir+tg_name[stn]+'_obs.nc','r')
   
       # EMODnet field extraction
       try:
          xin = fT1_obs.variables[var_2d_obs][:]*100.0 # want cm not meters
       except:
          xin = fT1_obs.variables[var_2d_obs][:,0]*100.0 # want cm not meters
       xin_mean=np.nanmean(xin)
       xin_obs_sub=xin-xin_mean
       ## Read quality flag
       #xin_qf_obs=fT1_obs.variables[var_2d_qf_obs][:,0]
       # Read latitude from model file
       latitudes=fT1_mod.variables[lat_var_mod][0]
       ## Read position quality flag
       #pos_qf_obs=fT1_obs.variables[pos_2d_qf_obs][:,0]
       # Set the time
       try:
          starttime=fT1_obs.variables[time_var_obs][0]
          time_var_units = fT1_obs.variables[time_var_obs].getncattr('units')
       except:
          starttime=fT1_obs.variables[time_var_obs2][0]
          time_var_units = fT1_obs.variables[time_var_obs2].getncattr('units')

       time=datetime(NC.num2date(starttime,time_var_units).year,NC.num2date(starttime,time_var_units).month,NC.num2date(starttime,time_var_units).day,NC.num2date(starttime,time_var_units).hour,NC.num2date(starttime,time_var_units).minute,NC.num2date(starttime,time_var_units).second)
       tg_sdate.append(datetime(NC.num2date(starttime,time_var_units).year,NC.num2date(starttime,time_var_units).month,NC.num2date(starttime,time_var_units).day))
       print ('Extracting EMODnet obs.. TG LAT TIME', tg_name[stn] , latitudes, time )
   
       # run ttide script for the emodnet obs
       harmanOUT = ttide.t_tide(np.squeeze(np.array(xin_obs_sub)), dt=1, stime=time, lat=latitudes, constitnames=tidal_comp, out_style=None, outfile=None,synth=snr_thershold)
       tideconout=harmanOUT['tidecon']
       snr=harmanOUT['snr']

       # Write values in the arrays
       # Amplitudes:
       M2_amp_obs[stn]=tideconout[5][0]
       S2_amp_obs[stn]=tideconout[6][0]
       K1_amp_obs[stn]=tideconout[3][0]
       O1_amp_obs[stn]=tideconout[1][0]
       N2_amp_obs[stn]=tideconout[4][0]
       P1_amp_obs[stn]=tideconout[2][0]
       Q1_amp_obs[stn]=tideconout[0][0]
       K2_amp_obs[stn]=tideconout[7][0]
       # Phases:
       M2_pha_obs[stn]=tideconout[5][2]
       S2_pha_obs[stn]=tideconout[6][2]
       K1_pha_obs[stn]=tideconout[3][2]
       O1_pha_obs[stn]=tideconout[1][2]
       N2_pha_obs[stn]=tideconout[4][2]
       P1_pha_obs[stn]=tideconout[2][2]
       Q1_pha_obs[stn]=tideconout[0][2]
       K2_pha_obs[stn]=tideconout[7][2]
       # Signal-to-noise ratio:
       M2_snr_obs[stn]=snr[5]
       S2_snr_obs[stn]=snr[6]
       K1_snr_obs[stn]=snr[3]
       O1_snr_obs[stn]=snr[1]
       N2_snr_obs[stn]=snr[4]
       P1_snr_obs[stn]=snr[2]
       Q1_snr_obs[stn]=snr[0]
       K2_snr_obs[stn]=snr[7]

       if errbar_flag != 0:
          # Amplitude errors
          M2_aerr_obs[stn]=tideconout[5][1]
          S2_aerr_obs[stn]=tideconout[6][1]
          K1_aerr_obs[stn]=tideconout[3][1]
          O1_aerr_obs[stn]=tideconout[1][1]
          N2_aerr_obs[stn]=tideconout[4][1]
          P1_aerr_obs[stn]=tideconout[2][1]
          Q1_aerr_obs[stn]=tideconout[0][1]
          K2_aerr_obs[stn]=tideconout[7][1]
          # Phase errors
          M2_perr_obs[stn]=tideconout[5][3]
          S2_perr_obs[stn]=tideconout[6][3]
          K1_perr_obs[stn]=tideconout[3][3]
          O1_perr_obs[stn]=tideconout[1][3]
          N2_perr_obs[stn]=tideconout[4][3]
          P1_perr_obs[stn]=tideconout[2][3]
          Q1_perr_obs[stn]=tideconout[0][3]
          K2_perr_obs[stn]=tideconout[7][3]

       else:
          # Amplitude errors
          M2_aerr_obs[stn]=0
          S2_aerr_obs[stn]=0
          K1_aerr_obs[stn]=0
          O1_aerr_obs[stn]=0
          N2_aerr_obs[stn]=0
          P1_aerr_obs[stn]=0
          Q1_aerr_obs[stn]=0
          K2_aerr_obs[stn]=0
          # Phase errors
          M2_perr_obs[stn]=0
          S2_perr_obs[stn]=0
          K1_perr_obs[stn]=0
          O1_perr_obs[stn]=0
          N2_perr_obs[stn]=0
          P1_perr_obs[stn]=0
          Q1_perr_obs[stn]=0
          K2_perr_obs[stn]=0
     
#       ####### Salish Sea #######
#
#       freq=harmanOUT['fu']
#       phazero=harmanOUT['pha']
#       ampnc=harmanOUT['f']
#       phanc=harmanOUT['vu']
#
#       # Read tidal components frequencies from ttide
#       M2freq=(freq[5]*360)*np.pi/180
#       S2freq=(freq[6]*360)*np.pi/180
#       K1freq=(freq[3]*360)*np.pi/180
#       O1freq=(freq[1]*360)*np.pi/180
#       N2freq=(freq[4]*360)*np.pi/180
#       P1freq=(freq[2]*360)*np.pi/180
#       Q1freq=(freq[0]*360)*np.pi/180
#       K2freq=(freq[7]*360)*np.pi/180
#
#       # Read phase0 from ttide
#       M2phazero[stn]=phazero[5][0]
#       S2phazero[stn]=phazero[6][0]
#       K1phazero[stn]=phazero[3][0]
#       O1phazero[stn]=phazero[1][0]
#       N2phazero[stn]=phazero[4][0]
#       P1phazero[stn]=phazero[2][0]
#       Q1phazero[stn]=phazero[0][0]
#       K2phazero[stn]=phazero[7][0]
#
#       # Read amplitude nodal correction
#       M2ampnc[stn]=ampnc[5]
#       S2ampnc[stn]=ampnc[6]
#       K1ampnc[stn]=ampnc[3]
#       O1ampnc[stn]=ampnc[1]
#       N2ampnc[stn]=ampnc[4]
#       P1ampnc[stn]=ampnc[2]
#       Q1ampnc[stn]=ampnc[0]
#       K2ampnc[stn]=ampnc[7]
#
#       # Read pha nodal correction
#       M2phanc[stn]=phanc[5]
#       S2phanc[stn]=phanc[6]
#       K1phanc[stn]=phanc[3]
#       O1phanc[stn]=phanc[1]
#       N2phanc[stn]=phanc[4]
#       P1phanc[stn]=phanc[2]
#       Q1phanc[stn]=phanc[0]
#       K2phanc[stn]=phanc[7]
#
#       def octuple_obs(x, M2_ssa_obs, M2_ssp_obs, K1_ssa_obs, K1_ssp_obs, O1_ssa_obs, O1_ssp_obs, S2_ssa_obs, S2_ssp_obs, P1_ssa_obs, P1_ssp_obs, N2_ssa_obs, N2_ssp_obs, Q1_ssa_obs, Q1_ssp_obs, K2_ssa_obs, K2_ssp_obs):
#           return (M2ampnc[stn]*M2_ssa_obs*np.cos(M2freq*x-M2_ssp_obs*np.pi/180.+M2phanc[stn]+M2phazero[stn])+
#            K1ampnc[stn]*K1_ssa_obs*np.cos(K1freq*x-K1_ssp_obs*np.pi/180.+K1phanc[stn]+K1phazero[stn])+
#            O1ampnc[stn]*O1_ssa_obs*np.cos(O1freq*x-O1_ssp_obs*np.pi/180.+O1phanc[stn]+O1phazero[stn])+
#            S2ampnc[stn]*S2_ssa_obs*np.cos(S2freq*x-S2_ssp_obs*np.pi/180.+S2phanc[stn]+S2phazero[stn])+
#            P1ampnc[stn]*P1_ssa_obs*np.cos(P1freq*x-P1_ssp_obs*np.pi/180.+P1phanc[stn]+P1phazero[stn])+
#            N2ampnc[stn]*N2_ssa_obs*np.cos(N2freq*x-N2_ssp_obs*np.pi/180.+N2phanc[stn]+N2phazero[stn])+
#            Q1ampnc[stn]*Q1_ssa_obs*np.cos(Q1freq*x-Q1_ssp_obs*np.pi/180.+Q1phanc[stn]+Q1phazero[stn])+
#            K2ampnc[stn]*K2_ssa_obs*np.cos(K2freq*x-K2_ssp_obs*np.pi/180.+K2phanc[stn]+K2phazero[stn]))
#
#       sstime_obs=(fT1_obs.variables[time_var_obs]-starttime)*24
#       xinnotnan=[]
#       sstimenotnan_obs=[]
#       for idx_rmnan in range (0,len(xin_obs_sub)):
#           if xin_obs_sub[idx_rmnan] != 'nan':
#              xinnotnan.append(xin_obs_sub[idx_rmnan])
#              sstimenotnan_obs.append(sstime_obs[idx_rmnan])
#       try:
#         fitted_obs, cov_obs = curve_fit(octuple_obs,sstimenotnan_obs,np.array(xinnotnan))
#         for idx_atobecorr in (0,2,4,6,8,10,12,14):
#           if fitted_obs[idx_atobecorr]<0:
#              fitted_obs[idx_atobecorr]=-fitted_obs[idx_atobecorr]
#              fitted_obs[idx_atobecorr+1]=fitted_obs[idx_atobecorr+1]+180
#         ss_amp[stn]=[fitted_obs[0],fitted_obs[6],fitted_obs[2],fitted_obs[4],fitted_obs[10],fitted_obs[8],fitted_obs[12],fitted_obs[14]]
#
#         for idx_ptobecorr in (1,3,5,7,9,11,13,15):
#           while fitted_obs[idx_ptobecorr]<0:  
#                 fitted_obs[idx_ptobecorr]=fitted_obs[idx_ptobecorr]+360
#           while fitted_obs[idx_ptobecorr]>360:
#                 fitted_obs[idx_ptobecorr]=fitted_obs[idx_ptobecorr]-360 
#         ss_pha[stn]=[fitted_obs[1],fitted_obs[7],fitted_obs[3],fitted_obs[5],fitted_obs[11],fitted_obs[9],fitted_obs[13],fitted_obs[15]]
#       except:
#         ss_amp[stn]=[0 for i in range(0,8)]
#         ss_pha[stn]=[0 for i in range(0,8)]
      
       ##################### END Salish Sea Method #############       

       # MODEL field extraction (for EMODnet TG)
       xin = fT1_mod.variables[var_2d_mod][:,0,0] *100.0 # want cm not meters
       xin=np.array(xin)
       xin_mean=np.nanmean(xin)
       xin_mod_sub=xin-xin_mean
       latitudes=fT1_mod.variables[lat_var_mod][0]
       # Set the time
       try:
         starttime=fT1_mod.variables[time_var_mod][0]
         time_var_units = fT1_mod.variables[time_var_mod].getncattr('units')
       except:
         starttime=fT1_mod.variables[time_var_mod2][0]
         time_var_units = fT1_mod.variables[time_var_mod2].getncattr('units')

       time=datetime(NC.num2date(starttime,time_var_units).year,NC.num2date(starttime,time_var_units).month,NC.num2date(starttime,time_var_units).day,NC.num2date(starttime,time_var_units).hour,NC.num2date(starttime,time_var_units).minute,NC.num2date(starttime,time_var_units).second)
       print ('Extracting.. TG LAT TIME', tg_name[stn] , latitudes, time)
   
       # run ttide script for the model
       harmanOUT = ttide.t_tide(np.array(xin_mod_sub), dt=1, stime=time, lat=latitudes, constitnames=tidal_comp, out_style=None, outfile=None,synth=snr_thershold)
       tideconout=harmanOUT['tidecon']
       snr=harmanOUT['snr']

       # Write values in the arrays
       # Amplitudes
       M2_amp_mod[stn]=tideconout[5][0]
       S2_amp_mod[stn]=tideconout[6][0]
       K1_amp_mod[stn]=tideconout[3][0]
       O1_amp_mod[stn]=tideconout[1][0]
       N2_amp_mod[stn]=tideconout[4][0]
       P1_amp_mod[stn]=tideconout[2][0]
       Q1_amp_mod[stn]=tideconout[0][0]
       K2_amp_mod[stn]=tideconout[7][0]
       # Phases:
       M2_pha_mod[stn]=tideconout[5][2]
       S2_pha_mod[stn]=tideconout[6][2]
       K1_pha_mod[stn]=tideconout[3][2]
       O1_pha_mod[stn]=tideconout[1][2]
       N2_pha_mod[stn]=tideconout[4][2]
       P1_pha_mod[stn]=tideconout[2][2]
       Q1_pha_mod[stn]=tideconout[0][2]
       K2_pha_mod[stn]=tideconout[7][2]
       # Signal-to-noise ratio:
       M2_snr_mod[stn]=snr[5]
       S2_snr_mod[stn]=snr[6]
       K1_snr_mod[stn]=snr[3]
       O1_snr_mod[stn]=snr[1]
       N2_snr_mod[stn]=snr[4]
       P1_snr_mod[stn]=snr[2]
       Q1_snr_mod[stn]=snr[0]
       K2_snr_mod[stn]=snr[7]
       if errbar_flag != 0:
          # Amplitude errors
          M2_aerr_mod[stn]=tideconout[5][1]
          S2_aerr_mod[stn]=tideconout[6][1]
          K1_aerr_mod[stn]=tideconout[3][1]
          O1_aerr_mod[stn]=tideconout[1][1]
          N2_aerr_mod[stn]=tideconout[4][1]
          P1_aerr_mod[stn]=tideconout[2][1]
          Q1_aerr_mod[stn]=tideconout[0][1]
          K2_aerr_mod[stn]=tideconout[7][1]
          # Phase errors
          M2_perr_mod[stn]=tideconout[5][3]
          S2_perr_mod[stn]=tideconout[6][3]
          K1_perr_mod[stn]=tideconout[3][3]
          O1_perr_mod[stn]=tideconout[1][3]
          N2_perr_mod[stn]=tideconout[4][3]
          P1_perr_mod[stn]=tideconout[2][3]
          Q1_perr_mod[stn]=tideconout[0][3]
          K2_perr_mod[stn]=tideconout[7][3]

       else:
          # Amplitude errors
          M2_aerr_mod[stn]=0
          S2_aerr_mod[stn]=0
          K1_aerr_mod[stn]=0
          O1_aerr_mod[stn]=0
          N2_aerr_mod[stn]=0
          P1_aerr_mod[stn]=0
          Q1_aerr_mod[stn]=0
          K2_aerr_mod[stn]=0
          # Phase errors
          M2_perr_mod[stn]=0
          S2_perr_mod[stn]=0
          K1_perr_mod[stn]=0
          O1_perr_mod[stn]=0
          N2_perr_mod[stn]=0
          P1_perr_mod[stn]=0
          Q1_perr_mod[stn]=0
          K2_perr_mod[stn]=0

#       ####### Salish Sea #######
#
#       freq=harmanOUT['fu']
#       phazero=harmanOUT['pha']
#       ampnc=harmanOUT['f']
#       phanc=harmanOUT['vu']
#
#       # Read tidal components frequencies from ttide
#       M2freq=(freq[5]*360)*np.pi/180
#       S2freq=(freq[6]*360)*np.pi/180
#       K1freq=(freq[3]*360)*np.pi/180
#       O1freq=(freq[1]*360)*np.pi/180
#       N2freq=(freq[4]*360)*np.pi/180
#       P1freq=(freq[2]*360)*np.pi/180
#       Q1freq=(freq[0]*360)*np.pi/180
#       K2freq=(freq[7]*360)*np.pi/180
#
#       # Read phase0 from ttide
#       M2phazero[stn]=phazero[5][0]
#       S2phazero[stn]=phazero[6][0]
#       K1phazero[stn]=phazero[3][0]
#       O1phazero[stn]=phazero[1][0]
#       N2phazero[stn]=phazero[4][0]
#       P1phazero[stn]=phazero[2][0]
#       Q1phazero[stn]=phazero[0][0]
#       K2phazero[stn]=phazero[7][0]
#
#       # Read amplitude nodal correction
#       M2ampnc[stn]=ampnc[5]
#       S2ampnc[stn]=ampnc[6]
#       K1ampnc[stn]=ampnc[3]
#       O1ampnc[stn]=ampnc[1]
#       N2ampnc[stn]=ampnc[4]
#       P1ampnc[stn]=ampnc[2]
#       Q1ampnc[stn]=ampnc[0]
#       K2ampnc[stn]=ampnc[7]
#
#       # Read pha nodal correction
#       M2phanc[stn]=phanc[5]
#       S2phanc[stn]=phanc[6]
#       K1phanc[stn]=phanc[3]
#       O1phanc[stn]=phanc[1]
#       N2phanc[stn]=phanc[4]
#       P1phanc[stn]=phanc[2]
#       Q1phanc[stn]=phanc[0]
#       K2phanc[stn]=phanc[7]
#
#       def octuple_mod(x, M2_ssa_mod, M2_ssp_mod, K1_ssa_mod, K1_ssp_mod, O1_ssa_mod, O1_ssp_mod, S2_ssa_mod, S2_ssp_mod, P1_ssa_mod, P1_ssp_mod, N2_ssa_mod, N2_ssp_mod, Q1_ssa_mod, Q1_ssp_mod, K2_ssa_mod, K2_ssp_mod):
#           return (M2ampnc[stn]*M2_ssa_mod*np.cos(M2freq*x-M2_ssp_mod*np.pi/180.+M2phanc[stn]+M2phazero[stn])+
#            K1ampnc[stn]*K1_ssa_mod*np.cos(K1freq*x-K1_ssp_mod*np.pi/180.+K1phanc[stn]+K1phazero[stn])+
#            O1ampnc[stn]*O1_ssa_mod*np.cos(O1freq*x-O1_ssp_mod*np.pi/180.+O1phanc[stn]+O1phazero[stn])+
#            S2ampnc[stn]*S2_ssa_mod*np.cos(S2freq*x-S2_ssp_mod*np.pi/180.+S2phanc[stn]+S2phazero[stn])+
#            P1ampnc[stn]*P1_ssa_mod*np.cos(P1freq*x-P1_ssp_mod*np.pi/180.+P1phanc[stn]+P1phazero[stn])+
#            N2ampnc[stn]*N2_ssa_mod*np.cos(N2freq*x-N2_ssp_mod*np.pi/180.+N2phanc[stn]+N2phazero[stn])+
#            Q1ampnc[stn]*Q1_ssa_mod*np.cos(Q1freq*x-Q1_ssp_mod*np.pi/180.+Q1phanc[stn]+Q1phazero[stn])+
#            K2ampnc[stn]*K2_ssa_mod*np.cos(K2freq*x-K2_ssp_mod*np.pi/180.+K2phanc[stn]+K2phazero[stn]))
#
#       try:
#          sstime_mod=(fT1_mod.variables[time_var_mod]-starttime)/60.0
#       except:
#          sstime_mod=(fT1_mod.variables[time_var_mod2]-starttime)/60.0
#       xinnotnan=[]
#       sstimenotnan_mod=[]
#       for idx_rmnan in range (0,len(xin_mod_sub)):
#           if xin_mod_sub[idx_rmnan] != 'nan':
#              xinnotnan.append(xin_mod_sub[idx_rmnan])
#              sstimenotnan_mod.append(sstime_mod[idx_rmnan])
#       try:
#         fitted_mod, cov_mod = curve_fit(octuple_mod,sstimenotnan_mod,np.array(xinnotnan))
#         for idx_atobecorr in (0,2,4,6,8,10,12,14):
#           if fitted_mod[idx_atobecorr]<0:
#              fitted_mod[idx_atobecorr]=-fitted_mod[idx_atobecorr]
#              fitted_mod[idx_atobecorr+1]=fitted_mod[idx_atobecorr+1]+180
#         ssmod_amp[stn]=[fitted_mod[0],fitted_mod[6],fitted_mod[2],fitted_mod[4],fitted_mod[10],fitted_mod[8],fitted_mod[12],fitted_mod[14]]
#
#         for idx_ptobecorr in (1,3,5,7,9,11,13,15):
#           while fitted_mod[idx_ptobecorr]<0:
#                 fitted_mod[idx_ptobecorr]=fitted_mod[idx_ptobecorr]+360
#           while fitted_mod[idx_ptobecorr]>360:
#                 fitted_mod[idx_ptobecorr]=fitted_mod[idx_ptobecorr]-360
#         ssmod_pha[stn]=[fitted_mod[1],fitted_mod[7],fitted_mod[3],fitted_mod[5],fitted_mod[11],fitted_mod[9],fitted_mod[13],fitted_mod[15]]
#       except:
#         ssmod_amp[stn]=[0 for i in range(0,8)]
#         ssmod_pha[stn]=[0 for i in range(0,8)]

       ######################
       # Plot time series and quality flag
       if revts_flag == 0:
          plotname=workdir_path+'ts_'+tg_name[stn]+'.jpg'
       else:
          plotname=workdir_path+'ts_'+tg_name[stn]+'_r.jpg'
       # Fig
       plt.figure(figsize=(20,5)) # (20,10) or (20,5) to compare with AtBox add plt.ylim(-250,250)
       plt.rc('font', size=16) # size=9 or size=16
       #plt.subplot(2,1,1)
       # Plot Title
       plt.title ('Time-series TG: '+tg_name[stn]+' Period: '+str(tg_sdate[stn])+' Lon/Lat: '+str(tg_lon[stn])+'/'+str(tg_lat[stn]))
       if revts_flag == 0:
          plt.plot(xin_mod_sub, '-', label = tg_name[stn]+' MOD')
          plt.plot(xin_obs_sub, '-', label = tg_name[stn]+' OBS')
       else:
          plt.plot(xin_obs_sub, '-', color='#ff7f0e',label = tg_name[stn]+' OBS')
          plt.plot(xin_mod_sub, '-', color='#1f77b4',label = tg_name[stn]+' MOD')
       plt.legend( loc='upper left',fontsize = 'large' )
       plt.grid ()
       plt.ylabel ('Sea Level [cm]')
       plt.xlabel ('Time [hours]')
       # 
       #plt.subplot(2,1,2)
       ## Plot Title
       #plt.title ('OBS Quality Flag')
       #plt.plot(xin_qf_obs, '--', label = comp+'Sea Leavel quality flag')
       #plt.plot(pos_qf_obs, '--', label = comp+'Position quality flag')
       #plt.legend( loc='upper left',fontsize = 'large' )
       #plt.grid ()
       #plt.ylabel ('Quality flag values')
       #plt.xlabel ('Time [hours]')
       #plt.ylim(-1, 11)
       # Save and close 
       plt.savefig(plotname)
       plt.clf()

       ########## Plot periodgrams (spectra)
       plotname=workdir_path+'spt_'+tg_name[stn]+'.jpg'

       spt_nc_obs=[]
       freq_nc_obs=[]
       freq_mod=[]
       spt_mod=[]

       rate_s=3600 # Hourly data
       spect_log=1 # log scale
       smooth_flag=1 # to smooth the spectrum
       spt_max=100000 # Max Spt

       # Compute the periodgrams
       spt_nc_obs=abs(np.fft.fft(xin_obs_sub[np.logical_not(np.isnan(xin_obs_sub))]))
       freq_nc_obs = abs(np.fft.fftfreq(xin_obs_sub[np.logical_not(np.isnan(xin_obs_sub))].shape[-1],rate_s))
       spt_mod=abs(np.fft.fft(xin_mod_sub[np.logical_not(np.isnan(xin_mod_sub))]))
       freq_mod = abs(np.fft.fftfreq(xin_mod_sub[np.logical_not(np.isnan(xin_mod_sub))].shape[-1],rate_s))

       # Create the plot
       plt.figure(figsize=(16,8))
       plt.rc('font', size=10)
       plt.grid ()

       plt.xlim(0.0000001,0.00015)
       plt.xscale('log')
       plt.xlabel ('Log(Frequency [Hz])')
       if spect_log == 1 :
          plt.yscale('log')
          plt.ylabel ('Log (Spectrum Amplitude)')
          plt.ylim(1,spt_max)
       elif spect_log == 0 :
          plt.ylabel ('Spectrum Amplitude')
          plt.ylim(1,20000)
       plt.title ('SSH'+' Spectrum - TG: '+tg_name[stn]+' Period: '+str(tg_sdate[stn])+' Lon/Lat: '+str(tg_lon[stn])+'/'+str(tg_lat[stn])+'\n')
       # Add tidal constituents freqs
       if spect_log == 0 :
          text_vertical_position=[50000,48000,50000,48000,44000,42000,42000,44000]
       elif spect_log == 1 :
          text_vertical_position=[10,8,10,8,4,2,2,4]
       vlines_colors=['black','black','green','green','green','green','black','black']
       # Add hour axes 
       freqs2hours_labels=['50','40','30','25','20','15','17','12','10','9','8','7','Period [Hours]']
       freqs2hours_val=[50.0,40.0,30.0,25.0,20.0,15.0,17.0,12.0,10.0,9.0,8.0,7.0,6.9]
       freqs2sec=np.asarray(freqs2hours_val)*3600
       freqs2hours_xposition=1.0/(freqs2sec)
       freqs2hours_yposition=[spt_max,spt_max,spt_max,spt_max,spt_max,spt_max,spt_max,spt_max,spt_max,spt_max,spt_max,spt_max,spt_max]
       # Add day axes 
       freqs2hours_labels=['100d','28d','15d','10d','8d','6d','4d','2d','24h','12h','10h','8h','6h','4h','2h']
       freqs2hours_xposition=[0.0000001157,0.0000004134,0.0000007716,0.000001157,0.000001447,0.000001929,0.000002894,0.000005787,0.00001157,0.00002315,0.00002778,0.00003472,0.00004630,0.00006944,0.0001389]
       freqs2hours_yposition=[spt_max,spt_max,spt_max,spt_max,spt_max,spt_max,spt_max,spt_max,spt_max,spt_max,spt_max,spt_max,spt_max,spt_max,spt_max,spt_max]
       for hidx in range(0,len(freqs2hours_labels)):
           plt.text(freqs2hours_xposition[hidx],freqs2hours_yposition[hidx],freqs2hours_labels[hidx],size=12)
           plt.axvline(x=freqs2hours_xposition[hidx], color='black') #, linestyle="dashed"

       if smooth_flag == 0:
          # Plot mod and obs spectrum
          plt.plot(freq_mod,spt_mod,color=spt_color,label=tg_name[stn]+' MOD')
          plt.plot(freq_nc_obs,spt_nc_obs,label=tg_name[stn]+' OBS')
       # Spectra Smoothing
       else:
          def moving_average(x, w):
              return np.convolve(x, np.ones(w), 'valid') / w
          num_avgpoints=11
          num_avgint=(num_avgpoints-1)/2
          mod_avg=moving_average(spt_mod, num_avgpoints)
          nc_obs_avg=moving_average(spt_nc_obs, num_avgpoints)
          #
          plt.plot(freq_mod[int(num_avgint):-int(num_avgint)],mod_avg,label=tg_name[stn]+' MOD')
          plt.plot(freq_nc_obs[int(num_avgint):-int(num_avgint)],nc_obs_avg,label=tg_name[stn]+' OBS')

       plt.legend( loc='upper right' )
       plt.savefig(plotname)
       plt.clf()

       ########## Plot empirical distrib functions
       plotname=workdir_path+'edf_'+tg_name[stn]+'.jpg'

       ecdf_nc_obs=[]
       ecdf_mod=[]
       # EDFs computation
       ecdf_nc_obs = ECDF(xin_obs_sub[np.logical_not(np.isnan(xin_obs_sub))])
       ecdf_mod = ECDF(xin_mod_sub[np.logical_not(np.isnan(xin_mod_sub))])
       # Plot
       plt.figure(figsize=(20,10))
       plt.rc('font', size=12)
       plt.grid ()
       #
       plt.title ('Empirical Distribution Function - TG: '+tg_name[stn]+' Period: '+str(tg_sdate[stn])+' Lon/Lat: '+str(tg_lon[stn])+'/'+str(tg_lat[stn]))
       plt.xlabel (var_2d_mod+field_2d_units)
       plt.ylabel ('EDF')
       plt.axhline(y=0.5, color='black', linestyle="dashed")
       plt.plot(ecdf_mod.x,ecdf_mod.y,label=tg_name[stn]+' MOD')
       plt.plot(ecdf_nc_obs.x,ecdf_nc_obs.y,label=tg_name[stn]+' OBS')
       plt.legend( loc='upper left' )
       #
       plt.savefig(plotname)
       plt.clf()

   # Store values in proper arrays
   for comp in tidal_comp:
       
       nameA_obs=comp+'_amp_obs'
       nameA_mod=comp+'_amp_mod'
   
       nameP_obs=comp+'_pha_obs'
       nameP_mod=comp+'_pha_mod'

       nameR_obs=comp+'_snr_obs'
       nameR_mod=comp+'_snr_mod'
   
   ##################################
   
   # Initialize the tables for Amp and Pha statistics
       # Table for TEX Amp
   if where_box=='Med':
      Amp_file = open(workdir_path+"amp_stats.txt","w")
   elif where_box=='AtlBox':
      Amp_file = open(workdir_path+"amp_stats_AB.txt","w")
   print('\\begin{table}'+'\n',file=Amp_file)
   print('\\footnotesize'+'\n',file=Amp_file)
   print('\hspace{-2.5cm}'+'\n',file=Amp_file)
   print('\\begin{tabular}{||c||c|c|c||c|c|c||}'+'\n',file=Amp_file)
   print('     \hline'+'\n',file=Amp_file)
   print('     \hline'+'\n',file=Amp_file)
   print('     & \multicolumn{3}{c||}{Absolute bias $\|Mod-Obs\|$} & \multicolumn{3}{c||}{Relative bias $\|\\frac{Mod-Obs}{Obs}\|\cdot 100$ } \\\\'+'\n',file=Amp_file)
   print('      & \multicolumn{3}{c||}{}& \multicolumn{3}{c||}{} \\\\'+'\n',file=Amp_file)
   print('     \hline'+'\n',file=Amp_file)
   print('     Tidal & Max & 95th percentile& Mean & Max & 95th percentile& Mean \\\\'+'\n',file=Amp_file)
   print('     components & & & & & & \\\\'+'\n',file=Amp_file)
   print('     \hline'+'\n',file=Amp_file)
   print('     \hline'+'\n',file=Amp_file)
   
       # Table for TEX Pha
   if where_box=='Med':
      Pha_file = open(workdir_path+"pha_stats.txt","w")
   elif where_box=='AtlBox':
      Pha_file = open(workdir_path+"pha_stats_AB.txt","w")
   print('\\begin{table}'+'\n',file=Pha_file)
   print('\\begin{tabular}{||c||c|c|c||}'+'\n',file=Pha_file)
   print('     \hline'+'\n',file=Pha_file)
   print('     \hline'+'\n',file=Pha_file)
   print('     & \multicolumn{3}{c||}{Absolute bias $\|Mod-Obs\|$} \\\\'+'\n',file=Pha_file)
   print('     \hline'+'\n',file=Pha_file)
   print('     Tidal & Max & 95th percentile & Mean \\\\'+'\n',file=Pha_file)
   print('     components & & & \\\\'+'\n',file=Pha_file)
   print('     \hline'+'\n',file=Pha_file)
   print('     \hline'+'\n',file=Pha_file)
   
   # Initialize the tables for lin reg statistics
   if where_box=='Med':
      LinReg_file = open(workdir_path+"linreg_stats.txt","w")
   elif where_box=='AtlBox':
      LinReg_file = open(workdir_path+"linreg_stats_AB.txt","w")
   print('\\begin{table}'+'\n',file=LinReg_file)
   print('\\begin{tabular}{||c||c|c||c|c||}'+'\n',file=LinReg_file)
   print('     \hline'+'\n',file=LinReg_file)
   print('     \hline'+'\n',file=LinReg_file)
   print('     & \multicolumn{2}{c||}{Amplitudes} & \multicolumn{2}{c||}{Phases} \\\\'+'\n',file=LinReg_file)
   print('     Tidal components & Slope & $R^{2}$ & Slope & $R^{2}$ \\\\'+'\n',file=LinReg_file)
   print('     \hline'+'\n',file=LinReg_file)
   print('     \hline'+'\n',file=LinReg_file)
   
   # Initialize global lists of arrays (with all the constituents) 
   N_comp=len(tidal_comp)
   N_stz=len(tg_name)
   GLOB_A_mod=[[ 0 for i in range(N_comp+1) ] for j in range(N_stz+1) ]
   GLOB_P_mod=[[ 0 for i in range(N_comp+1) ] for j in range(N_stz+1) ] 
   GLOB_A_obs=[[ 0 for i in range(N_comp+1) ] for j in range(N_stz+1) ] 
   GLOB_P_obs=[[ 0 for i in range(N_comp+1) ] for j in range(N_stz+1) ] 
   
   # Initialize global matrix for Foreman distances and Root Mean Square misfits
   d_foreman=[[ 0 for i in range(N_comp+1) ] for j in range(N_stz+1) ]
   RMSm=[0 for i in range(N_comp)]
   
   ################### SORT the TG, PLOT MAP and TABLES #####################

   # Get the info and sort the tg on longitude order
   ALL_tg_name=np.append(tg_name,[])
   ALL_tg_col=np.append(tg_col,[])
   ALL_tg_lon=np.append(tg_lon,[])
   ALL_tg_lat=np.append(tg_lat,[])
   ALL_tg_sdate=np.append(tg_sdate,[])

   lonsort_idx, tg_lon_sorted = zip(*sorted(enumerate(ALL_tg_lon), key=itemgetter(1)))
   tg_lat_sorted = np.take(ALL_tg_lat,lonsort_idx)
   tg_name_sorted = np.take(ALL_tg_name,lonsort_idx)
   tg_col_sorted = np.take(ALL_tg_col,lonsort_idx)
   tg_sdate_sorted = np.take(ALL_tg_sdate,lonsort_idx)

   # Define and sort the Amp and Pha arrays
   comp_idx=0
   for comp in tidal_comp:
      nameA_obs=comp+'_amp_obs'
      nameA_mod=comp+'_amp_mod'
      nameP_obs=comp+'_pha_obs'
      nameP_mod=comp+'_pha_mod'
      nameR_obs=comp+'_snr_obs'
      nameR_mod=comp+'_snr_mod'
      nameEA_obs=comp+'_aerr_obs'
      nameEA_mod=comp+'_aerr_mod'
      nameEP_obs=comp+'_perr_obs'
      nameEP_mod=comp+'_perr_mod'
      ALL_A_obs=np.append(np.array(globals()[nameA_obs]),[])
      ALL_A_mod=np.append(np.array(globals()[nameA_mod]),[])
      ALL_P_obs=np.append(np.array(globals()[nameP_obs]),[])
      ALL_P_mod=np.append(np.array(globals()[nameP_mod]),[])
      ALL_R_obs=np.append(np.array(globals()[nameR_obs]),[])
      ALL_R_mod=np.append(np.array(globals()[nameR_mod]),[])
      ALL_EA_obs=np.append(np.array(globals()[nameEA_obs]),[])
      ALL_EA_mod=np.append(np.array(globals()[nameEA_mod]),[])
      ALL_EP_obs=np.append(np.array(globals()[nameEP_obs]),[])
      ALL_EP_mod=np.append(np.array(globals()[nameEP_mod]),[])
      sort_A_obs=np.take(ALL_A_obs,lonsort_idx)
      sort_A_mod=np.take(ALL_A_mod,lonsort_idx)
      sort_P_obs=np.take(ALL_P_obs,lonsort_idx)
      sort_P_mod=np.take(ALL_P_mod,lonsort_idx)
      sort_R_obs=np.take(ALL_R_obs,lonsort_idx)
      sort_R_mod=np.take(ALL_R_mod,lonsort_idx)
      sort_EA_obs=np.take(ALL_EA_obs,lonsort_idx)
      sort_EA_mod=np.take(ALL_EA_mod,lonsort_idx)
      sort_EP_obs=np.take(ALL_EP_obs,lonsort_idx)
      sort_EP_mod=np.take(ALL_EP_mod,lonsort_idx)
      globals()['ALL_A_obs'+comp]=sort_A_obs
      globals()['ALL_A_mod'+comp]=sort_A_mod
      globals()['ALL_P_obs'+comp]=sort_P_obs
      globals()['ALL_P_mod'+comp]=sort_P_mod
      globals()['ALL_R_obs'+comp]=sort_R_obs
      globals()['ALL_R_mod'+comp]=sort_R_mod
      globals()['ALL_EA_obs'+comp]=sort_EA_obs
      globals()['ALL_EA_mod'+comp]=sort_EA_mod
      globals()['ALL_EP_obs'+comp]=sort_EP_obs
      globals()['ALL_EP_mod'+comp]=sort_EP_mod

#      # Salish Sea Method
#      globals()['SS_A_OBS'+comp]=[ss_amp[i][comp_idx] for i in range(0,N_stz)]
#      sort_ssa_obs=np.take(globals()['SS_A_OBS'+comp],lonsort_idx)
#      globals()['SS_A_OBS'+comp]=sort_ssa_obs
#
#      globals()['SS_A_MOD'+comp]=[ssmod_amp[i][comp_idx] for i in range(0,N_stz)]
#      sort_ssa_mod=np.take(globals()['SS_A_MOD'+comp],lonsort_idx)
#      globals()['SS_A_MOD'+comp]=sort_ssa_mod

      comp_idx=comp_idx+1

   # Define TG LABELS to be used in plots
   tg_lab_sorted=[]
   if anatype_flag=='lit':
      if lit_fullname == 1: # If you have enough space for the names!
         tg_lab_sorted = tg_name_sorted 
      else: # use the numbers
         for tg_idx in range(0,len(tg_name_sorted)):
             shift_idx=tg_idx+1
             tg_lab_sorted.append(str(shift_idx))
   elif anatype_flag=='anatpxo':
      for tg_idx in range(0,len(tg_name_sorted)):
          shift_idx=tg_idx+1
          tg_lab_sorted.append(str(shift_idx))
   elif anatype_flag=='all':
      for tg_idx in range(0,len(tg_name_sorted)):
          shift_idx=tg_idx+1
          tg_lab_sorted.append(str(shift_idx)) 

   ALL_tg_lat=tg_lat_sorted
   ALL_tg_lon=tg_lon_sorted
   ALL_tg_lab=tg_lab_sorted
   ALL_tg_name=tg_name_sorted
   ALL_tg_col=tg_col_sorted
   ALL_tg_sdate=tg_sdate_sorted

   # PLOT THE MAP of tide-gauge location
   nc2open3=model_bathy # tidal bathimetry
   model3 = NC.Dataset(nc2open3,'r')
   vals_bathy=model3.variables['Bathymetry'][:]
   lons = model3.variables['nav_lon'][:]
   lats = model3.variables['nav_lat'][:]
   if anatype_flag == 'lit':
      plotname=workdir_path+anatype_flag+'_tg.jpg'
   elif anatype_flag == 'anatpxo':
      plotname=workdir_path+anatype_flag+'_tg.jpg'
   elif anatype_flag == 'all':
      plotname=workdir_path+anatype_flag+'_tg.jpg'
   # Fig
   plt.figure(figsize=(20,10))
   plt.rc('font', size=12)
   # Plot Title
   plt.title ('Bathymetry [m] and Tide-Gauges location')
   lon_0 = lons.mean()
   llcrnrlon = lons.min()
   urcrnrlon = lons.max()
   if where_box=='Med':
     llcrnrlon = -10.0
     urcrnrlon = lons.max()
   elif where_box=='AtlBox':
     llcrnrlon = -18.0
     urcrnrlon = 0.0
   lat_0 = lats.mean()
   llcrnrlat = lats.min()
   urcrnrlat = lats.max()
   # Create the map
   m = Basemap(llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat,urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,resolution='c',projection='merc',lat_0=lat_0,lon_0=lon_0)
   xi, yi = m(lons, lats)
   # Plot the frame to the map
   plt.rcParams["axes.linewidth"]  = 1.25
   m.drawparallels(np.arange(30., 46., 5.), labels=[1,0,0,0], fontsize=10)
   m.drawmeridians(np.arange(-20., 40., 10.), labels=[0,0,0,1], fontsize=10)
   contourf = plt.contour(xi,yi,np.squeeze(vals_bathy),0.0,colors='black')
   # Plot the bathy
   cmap = mpl.cm.Blues(np.linspace(0,1,20))
   cmap = mpl.colors.ListedColormap(cmap[5:,:-1])
   cmap =  cmap.reversed()
   cs = m.pcolor(xi,yi,-np.squeeze(vals_bathy),cmap=cmap,vmax=-5000,vmin=0)
   contourf = plt.contourf(xi,yi,np.squeeze(vals_bathy),[-1000,0.0],colors='gray')
   # Plot the legend and its label
   cbar = m.colorbar(cs, location='right', pad="10%")
   bar_label_string='Bathymetry [m]'
   cbar.set_label(bar_label_string)
   # Add tide-gauges
   for tg2plot_idx in range(0,len(ALL_tg_name)) :
     xp, yp = m(ALL_tg_lon[tg2plot_idx],ALL_tg_lat[tg2plot_idx])
     plt.text(xp,yp,ALL_tg_lab[tg2plot_idx], fontsize=12,backgroundcolor=ALL_tg_col[tg2plot_idx],alpha=1,color='black')
   # Save and close 
   plt.savefig(plotname)
   plt.clf()

   # Write TAB with names, labels and coordinates
   TGTab_filename=workdir_path+anatype_flag+'_TGtab.txt'
   TGTab_file = open(TGTab_filename,"w")
   print('\\begin{table}'+'\n',file=TGTab_file)
   print('\\# & Name & Label & Longitude & Latitude & TS Start Date \\\\'+'\n',file=TGTab_file)
   print('\\hline'+'\n',file=TGTab_file)
   for tgtab_idx in range(0,len(ALL_tg_name)):
       shift_idx=tgtab_idx+1
       print(str(shift_idx)+' & '+str(ALL_tg_name[tgtab_idx])+' & '+'{\color{'+str(ALL_tg_col[tgtab_idx])+'}{'+str(ALL_tg_lab[tgtab_idx])+'}} & '+str(round(ALL_tg_lon[tgtab_idx],4))+' & '+str(round(ALL_tg_lat[tgtab_idx],4))+' & '+str(ALL_tg_sdate[tgtab_idx])+'\\\\'+'\n',file=TGTab_file)
   print('\\hline'+'\n',file=TGTab_file)
   TGTab_file.close()

   # Write TAB with signal-to-noise ratios
   SNRTab_filename=workdir_path+anatype_flag+'_SNRtab.txt'
   SNRTab_file = open(SNRTab_filename,"w")
   print('\\begin{table}'+'\n',file=SNRTab_file)
   print('\\# & M2 SNR obs & M2 SNR mod & S2 SNR obs & S2 SNR mod & K1 SNR obs & K1 SNR mod & O1 SNR obs & O1 SNR mod & N2 SNR obs & N2 SNR mod & P1 SNR obs & P1 SNR mod & Q1 SNR obs & Q1 SNR mod & K2 SNR obs & K2 SNR mod \\\\'+'\n',file=SNRTab_file)
   print('\\hline'+'\n',file=SNRTab_file)
   for tgtab_idx in range(0,len(ALL_tg_name)):
       shift_idx=tgtab_idx+1
       print(str(shift_idx)+'&'+str(ALL_R_obsM2[tgtab_idx])+'&'+str(ALL_R_modM2[tgtab_idx])+'&'+str(ALL_R_obsS2[tgtab_idx])+'&'+str(ALL_R_modS2[tgtab_idx])+'&'+str(ALL_R_obsK1[tgtab_idx])+'&'+str(ALL_R_modK1[tgtab_idx])+'&'+str(ALL_R_obsO1[tgtab_idx])+'&'+str(ALL_R_modO1[tgtab_idx])+'&'+str(ALL_R_obsN2[tgtab_idx])+'&'+str(ALL_R_modN2[tgtab_idx])+'&'+str(ALL_R_obsP1[tgtab_idx])+'&'+str(ALL_R_modP1[tgtab_idx])+'&'+str(ALL_R_obsQ1[tgtab_idx])+'&'+str(ALL_R_modQ1[tgtab_idx])+'&'+str(ALL_R_obsK2[tgtab_idx])+'&'+str(ALL_R_modK2[tgtab_idx])+'\\\\'+'\n',file=SNRTab_file)
   print('\\hline'+'\n',file=SNRTab_file)
   TGTab_file.close()

   # Signal-to-noise ratio plots
   comp_idx=0
   for comp in tidal_comp:
      #plotname=workdir_path+comp+'_snr.jpg'
      plotname=workdir_path+comp+'_fitdiag.jpg'
      # Fig
      plt.figure(figsize=(20,10))
      plt.rc('font', size=9)
      plt.subplot(2,1,1)
      # Plot Title
      plt.title ('Signal-to-noise ratio harmonic analysis fit --- tidal component '+comp)
      plt.plot(ALL_tg_lab,globals()['ALL_R_mod'+comp], '-o', label = comp+' SNR mod')
      plt.plot(ALL_tg_lab,globals()['ALL_R_obs'+comp], '-o', label = comp+' SNR obs')
      plt.legend( loc='upper left',fontsize = 'large' )
      plt.grid ()
      plt.yscale('log')
      plt.ylabel ('Signal-to-noise ratio')
      plt.xlabel ('Tide-gauges')
      # Save and close 
      plt.savefig(plotname)
      plt.clf()

#    # Methods comparison plot
#   #comp_idx=0
#   #for comp in tidal_comp:
#   #   plotname=workdir_path+comp+'_methods.jpg'
#   #   # Fig
#   #   plt.figure(figsize=(20,8))
#   #   plt.rc('font', size=9)
#      plt.subplot(2,1,2)
#      # Plot Title
#      plt.title ('Harmonic analysis Foreman Vs Salish Sea method results --- tidal component '+comp)
#      plt.plot(ALL_tg_lab,np.abs(globals()['SS_A_MOD'+comp]-globals()['ALL_A_mod'+comp]), '--o', label = comp+'||SalishSea-Foreman|| Amp Diff MOD')
#      plt.plot(ALL_tg_lab,np.abs(globals()['SS_A_OBS'+comp]-globals()['ALL_A_obs'+comp]), '--o', label = comp+'||SalishSea-Foreman|| Amp Diff OBS')
#      #plt.plot(ALL_tg_lab,globals()['ALL_A_obs'+comp], 'r--o', label = comp+' Foreman Amp OBS')
#      #plt.plot(ALL_tg_lab,globals()['ALL_A_mod'+comp], 'm-o', label = comp+' Foreman Amp MOD')
#      plt.legend( loc='upper left',fontsize = 'large' )
#      plt.grid ()
#      plt.ylabel ('Tidal Amplitudes Diff [cm]')
#      plt.xlabel ('Tide-gauges')
#      # Save and close 
#      plt.savefig(plotname)
#      plt.clf() 

      comp_idx=comp_idx+1

   ##############################################################
   # Loop on the 8 components 
   print ('Working on single tidal components..')
   comp_idx=0
   for comp in tidal_comp:
       print (comp)
   
       for regi_lab in subregions_labels:
           globals()['TOT_tg_lab_'+regi_lab+'area']=[]; globals()['TOT_tg_name_'+regi_lab+'area']=[]
           globals()['TOT_A_obs_'+regi_lab+'area']=[]; globals()['TOT_P_obs_'+regi_lab+'area']=[]
           globals()['TOT_A_mod_'+regi_lab+'area']=[]; globals()['TOT_P_mod_'+regi_lab+'area']=[]
           globals()['TOT_EA_obs_'+regi_lab+'area']=[]; globals()['TOT_EP_obs_'+regi_lab+'area']=[]
           globals()['TOT_EA_mod_'+regi_lab+'area']=[]; globals()['TOT_EP_mod_'+regi_lab+'area']=[]

       TOT_tg_lab_ord=[]
       TOT_tg_name_ord=[]
       TOT_A_obs_ord=[]
       TOT_P_obs_ord=[]
       TOT_A_mod_ord=[]
       TOT_P_mod_ord=[]
   
       # Pass array to the analysis looop
       # WARNING: the values of this arrays change in each loop, thus here they should be defined!
       TOT_tg_lat=ALL_tg_lat
       TOT_tg_lon=ALL_tg_lon
       TOT_tg_col=ALL_tg_col
       TOT_tg_name=ALL_tg_name
       TOT_tg_lab=ALL_tg_lab
       #
       TOT_A_obs=globals()['ALL_A_obs'+comp]
       TOT_A_mod=globals()['ALL_A_mod'+comp]
       TOT_P_obs=globals()['ALL_P_obs'+comp]
       TOT_P_mod=globals()['ALL_P_mod'+comp]
       TOT_EA_obs=globals()['ALL_EA_obs'+comp]
       TOT_EA_mod=globals()['ALL_EA_mod'+comp]
       TOT_EP_obs=globals()['ALL_EP_obs'+comp]
       TOT_EP_mod=globals()['ALL_EP_mod'+comp]


       # OBS Adjustment for linear regressions and phase plots 
       for idx_diff_mo in range(0,len(TOT_P_mod)):

           if TOT_P_mod[idx_diff_mo]-TOT_P_obs[idx_diff_mo] > 200:
              TOT_P_obs[idx_diff_mo]=TOT_P_obs[idx_diff_mo]+360
              if TOT_P_obs[idx_diff_mo]>400 or TOT_P_mod[idx_diff_mo]>400:
                 TOT_P_obs[idx_diff_mo]=TOT_P_obs[idx_diff_mo]-360
                 TOT_P_mod[idx_diff_mo]=TOT_P_mod[idx_diff_mo]-360

           elif TOT_P_mod[idx_diff_mo]-TOT_P_obs[idx_diff_mo] < -200 :
              TOT_P_obs[idx_diff_mo]=TOT_P_obs[idx_diff_mo]-360
              if TOT_P_obs[idx_diff_mo]<-50 or TOT_P_mod[idx_diff_mo]<-50:
                 TOT_P_obs[idx_diff_mo]=TOT_P_obs[idx_diff_mo]+360
                 TOT_P_mod[idx_diff_mo]=TOT_P_mod[idx_diff_mo]+360

       #############################
       # Split in different datasets: Gibraltar area, Adriatic area and others
       # Count stz per area (color based)
       regi_idx=0
       for regi_lab in subregions_labels:
         globals()['howmany_'+regi_lab+'area']=0
         for idx_area in range(0,len(TOT_tg_col)):
           if TOT_tg_col[idx_area] == subregions_color[regi_idx]: 
              globals()['TOT_tg_lab_'+regi_lab+'area'].append(TOT_tg_lab[idx_area])
              globals()['TOT_tg_name_'+regi_lab+'area'].append(TOT_tg_name[idx_area])
              globals()['TOT_A_obs_'+regi_lab+'area'].append(TOT_A_obs[idx_area])
              globals()['TOT_P_obs_'+regi_lab+'area'].append(TOT_P_obs[idx_area])   
              globals()['TOT_A_mod_'+regi_lab+'area'].append(TOT_A_mod[idx_area])
              globals()['TOT_P_mod_'+regi_lab+'area'].append(TOT_P_mod[idx_area])
              globals()['TOT_EA_obs_'+regi_lab+'area'].append(TOT_EA_obs[idx_area])
              globals()['TOT_EP_obs_'+regi_lab+'area'].append(TOT_EP_obs[idx_area])
              globals()['TOT_EA_mod_'+regi_lab+'area'].append(TOT_EA_mod[idx_area])
              globals()['TOT_EP_mod_'+regi_lab+'area'].append(TOT_EP_mod[idx_area])
              globals()['howmany_'+regi_lab+'area']=globals()['howmany_'+regi_lab+'area']+1
         regi_idx=regi_idx+1

       # Reorder the array 
       # WARNING: this cannot be changed because the plots etc follows this order..
       outregion_order=['Oarea','Marea','Garea','Aarea','Tarea','TAarea','TGarea','Earea']

       TOT_A_obs_ord=globals()['TOT_A_obs_'+outregion_order[0]]
       TOT_P_obs_ord=globals()['TOT_P_obs_'+outregion_order[0]]
       TOT_A_mod_ord=globals()['TOT_A_mod_'+outregion_order[0]]
       TOT_P_mod_ord=globals()['TOT_P_mod_'+outregion_order[0]]
       TOT_EA_obs_ord=globals()['TOT_EA_obs_'+outregion_order[0]]
       TOT_EP_obs_ord=globals()['TOT_EP_obs_'+outregion_order[0]]
       TOT_EA_mod_ord=globals()['TOT_EA_mod_'+outregion_order[0]]
       TOT_EP_mod_ord=globals()['TOT_EP_mod_'+outregion_order[0]]
       TOT_tg_lab_ord=globals()['TOT_tg_lab_'+outregion_order[0]]
       TOT_tg_name_ord=globals()['TOT_tg_name_'+outregion_order[0]]
       for outregi_idx in outregion_order[1:]:
           TOT_A_obs_ord.extend(globals()['TOT_A_obs_'+outregi_idx])
           TOT_P_obs_ord.extend(globals()['TOT_P_obs_'+outregi_idx])
           TOT_A_mod_ord.extend(globals()['TOT_A_mod_'+outregi_idx])
           TOT_P_mod_ord.extend(globals()['TOT_P_mod_'+outregi_idx])
           TOT_EA_obs_ord.extend(globals()['TOT_EA_obs_'+outregi_idx])
           TOT_EP_obs_ord.extend(globals()['TOT_EP_obs_'+outregi_idx])
           TOT_EA_mod_ord.extend(globals()['TOT_EA_mod_'+outregi_idx])
           TOT_EP_mod_ord.extend(globals()['TOT_EP_mod_'+outregi_idx])
           TOT_tg_lab_ord.extend(globals()['TOT_tg_lab_'+outregi_idx])
           TOT_tg_name_ord.extend(globals()['TOT_tg_name_'+outregi_idx])

   
       TOT_A_obs=[]
       TOT_P_obs=[]
       TOT_A_mod=[]
       TOT_P_mod=[]
       TOT_EA_obs=[]
       TOT_EP_obs=[]
       TOT_EA_mod=[]
       TOT_EP_mod=[]
       TOT_tg_lab=[]
       TOT_tg_name=[]
       TOT_tg_orcol=[]
   
       TOT_A_obs=np.array(TOT_A_obs_ord)
       TOT_P_obs=np.array(TOT_P_obs_ord)
       TOT_A_mod=np.array(TOT_A_mod_ord)
       TOT_P_mod=np.array(TOT_P_mod_ord)
       TOT_EA_obs=np.array(TOT_EA_obs_ord)
       TOT_EP_obs=np.array(TOT_EP_obs_ord)
       TOT_EA_mod=np.array(TOT_EA_mod_ord)
       TOT_EP_mod=np.array(TOT_EP_mod_ord)
       TOT_tg_lab=TOT_tg_lab_ord
       TOT_tg_name=TOT_tg_name_ord
       
       # Shift the color order along with the region order
       oc_idx=(4,2,0,1,7,5,6,3)
       an_idx=0
       for a_idx in outregion_order:
           howmany_idx=globals()['howmany_'+a_idx]
           for c_idx in range(0,howmany_idx):
               new_col=oc_idx[an_idx]
               TOT_tg_orcol.append(subregions_color[new_col])
           an_idx=an_idx+1

       ###### TPXO and LIT extraction ######
       if tpxo_flag == 1:
          globals()['TPXO_'+comp]=[]
          globals()['TPXO_P_'+comp]=[]
          for tg2beadded in TOT_tg_name:
              A_where=[x for x in globals()['INTPXO_'+comp] if tg2beadded in x][0][1]
              P_where=[x for x in globals()['INTPXO_P_'+comp] if tg2beadded in x][0][1]
              globals()['TPXO_'+comp].append(A_where)
              globals()['TPXO_P_'+comp].append(P_where)
   
          # TPXO phase Adjustment for linear regressions and phase plots 
          for idx_tpxoadj in range(0,len(TOT_tg_name)):
              if TOT_P_mod[idx_tpxoadj]-globals()['TPXO_P_'+comp][idx_tpxoadj] > 200:
                 globals()['TPXO_P_'+comp][idx_tpxoadj]=globals()['TPXO_P_'+comp][idx_tpxoadj]+360
              elif TOT_P_mod[idx_tpxoadj]-globals()['TPXO_P_'+comp][idx_tpxoadj] < -200 :
                 globals()['TPXO_P_'+comp][idx_tpxoadj]=globals()['TPXO_P_'+comp][idx_tpxoadj]-360

       if comp_idx < 4 and (flag_15stats == 1 or tpxo_flag == 1):
   
          globals()['PALMA_'+comp]=[]
          globals()['PALMA_P_'+comp]=[]
          globals()['PALMA_d_'+comp]=[]
          globals()['TSIMPLIS_'+comp]=[]
          globals()['TSIMPLIS_P_'+comp]=[]
          globals()['TSIMPLIS_d_'+comp]=[]
   
          for idx_tglit,tg2beadded in enumerate(TOT_tg_name):
            try:
              A_where=[x for x in globals()['INPALMA_'+comp] if tg2beadded in x][0][1] 
              P_where=[x for x in globals()['INPALMA_P_'+comp] if tg2beadded in x][0][1] 
              d_where=[x for x in globals()['INPALMA_d_'+comp] if tg2beadded in x][0][1] 
              # LITERATURE phase Adjustment for phase plots
              if TOT_P_mod[idx_tglit]-P_where > 200:
                 P_where=P_where+360
              elif TOT_P_mod[idx_tglit]-P_where < -200:
                 P_where=P_where-360

              globals()['PALMA_'+comp].append(A_where)
              globals()['PALMA_P_'+comp].append(P_where)
              globals()['PALMA_d_'+comp].append(d_where)
              #
              A_where=[x for x in globals()['INTSIMPLIS_'+comp] if tg2beadded in x][0][1]
              P_where=[x for x in globals()['INTSIMPLIS_P_'+comp] if tg2beadded in x][0][1]
              d_where=[x for x in globals()['INTSIMPLIS_d_'+comp] if tg2beadded in x][0][1]

              # LITERATURE phase Adjustment for phase plots
              if TOT_P_mod[idx_tglit]-P_where > 200:
                 P_where=P_where+360
              elif TOT_P_mod[idx_tglit]-P_where < -200:
                 P_where=P_where-360

              globals()['TSIMPLIS_'+comp].append(A_where)
              globals()['TSIMPLIS_P_'+comp].append(P_where)
              globals()['TSIMPLIS_d_'+comp].append(d_where)
 
              #
            except:
              globals()['PALMA_'+comp].append('')
              globals()['PALMA_P_'+comp].append('')
              globals()['PALMA_d_'+comp].append('')
              globals()['TSIMPLIS_'+comp].append('')
              globals()['TSIMPLIS_P_'+comp].append('')
              globals()['TSIMPLIS_d_'+comp].append('')

       #####################################
   
       # Diff and diff_mean computation
       # diff Amp
       x_textA=[]
       y_textA=[]
       x_textA=TOT_A_obs
       y_textA=TOT_A_mod
       diffA_mo=y_textA-x_textA
   
   
       diffA_mo_Oarea=diffA_mo[0:howmany_Oarea]
       diffA_mo_Marea=diffA_mo[howmany_Oarea:howmany_Oarea+howmany_Marea]
       diffA_mo_Garea=diffA_mo[howmany_Oarea+howmany_Marea:howmany_Oarea+howmany_Marea+howmany_Garea]
       diffA_mo_Aarea=diffA_mo[howmany_Oarea+howmany_Marea+howmany_Garea:howmany_Oarea+howmany_Marea+howmany_Garea+howmany_Aarea]
       diffA_mo_Tarea=diffA_mo[howmany_Oarea+howmany_Marea+howmany_Garea+howmany_Aarea:howmany_Oarea+howmany_Marea+howmany_Garea+howmany_Aarea+howmany_Tarea]
       diffA_mo_TAarea=diffA_mo[howmany_Oarea+howmany_Marea+howmany_Garea+howmany_Aarea+howmany_Tarea:howmany_Oarea+howmany_Marea+howmany_Garea+howmany_Aarea+howmany_Tarea+howmany_TAarea]
       diffA_mo_TGarea=diffA_mo[howmany_Oarea+howmany_Marea+howmany_Garea+howmany_Aarea+howmany_Tarea+howmany_TAarea:howmany_Oarea+howmany_Marea+howmany_Garea+howmany_Aarea+howmany_Tarea+howmany_TAarea+howmany_TGarea]
       diffA_mo_Earea=diffA_mo[howmany_Oarea+howmany_Marea+howmany_Garea+howmany_Aarea+howmany_Tarea+howmany_TAarea+howmany_TGarea:howmany_Oarea+howmany_Marea+howmany_Garea+howmany_Aarea+howmany_Tarea+howmany_TAarea+howmany_TGarea+howmany_Earea]
   
   
       pdiffA_mo=(y_textA-x_textA)*100.0/(0.5+x_textA) # do not consider values below 0.5 cm
       mean_diffA=np.mean(diffA_mo)
   
       # STATISTICS x TAB AMP
       meanabs_diffA=np.mean(abs(diffA_mo))
       max_diffA=np.max(abs(diffA_mo))
       wheremax_AmpAbs=TOT_tg_name[np.argmax(abs(diffA_mo))]
       wheremax_AmpAbs_col=TOT_tg_orcol[np.argmax(abs(diffA_mo))]
       perc95_diffA=np.percentile(abs(diffA_mo),95)
       mean_pdiffA=np.mean(abs(pdiffA_mo))
       meanabs_pdiffA=np.mean(abs(pdiffA_mo))
       max_pdiffA=np.max(abs(pdiffA_mo))
       wheremax_AmpPerc=TOT_tg_name[np.argmax(abs(pdiffA_mo))]
       wheremax_AmpPerc_col=TOT_tg_orcol[np.argmax(abs(pdiffA_mo))]
       perc95_pdiffA=np.percentile(abs(pdiffA_mo),95)
   
       # Table for TEX Amp
       print(comp,' &',str(round(max_diffA,2)),' cm ({\color{',wheremax_AmpAbs_col,'}{',wheremax_AmpAbs,'}})&',str(round(perc95_diffA,2)),' cm &',str(round(meanabs_diffA,1)),' cm &',str(round(max_pdiffA,2)),' \% ({\color{',wheremax_AmpPerc_col,'}{',wheremax_AmpPerc,'}})&',str(round(perc95_pdiffA,2)),' \% &',str(round(meanabs_pdiffA,1)),'\% \\\\ '+'\n', file=Amp_file)
       print('\hline'+'\n', file=Amp_file)
   
   
       # diff Pha
       x_textP=[]
       y_textP=[]
       x_textP=TOT_P_obs
       y_textP=TOT_P_mod
       if cos_pha == 0:
          diffP_mo=y_textP-x_textP
          pdiffP_mo=(y_textP-x_textP)*100.0/(2+x_textP) # do not consider values below 2 deg
          # Constrain the pha diff and pdiff angles in the interval [-180,180):
          for Pdiff_idx in range(0,len(diffP_mo)):
              while diffP_mo[Pdiff_idx] >= 180:
                    diffP_mo[Pdiff_idx]=diffP_mo[Pdiff_idx]-360
              while diffP_mo[Pdiff_idx] <-180:
                    diffP_mo[Pdiff_idx]=diffP_mo[Pdiff_idx]+360
          for Ppdiff_idx in range(0,len(pdiffP_mo)):
              while pdiffP_mo[Ppdiff_idx] >= 180:
                    pdiffP_mo[Ppdiff_idx]=pdiffP_mo[Ppdiff_idx]-360
              while pdiffP_mo[Ppdiff_idx] <-180:
                    pdiffP_mo[Ppdiff_idx]=pdiffP_mo[Ppdiff_idx]+360
       elif cos_pha == 1:
          diffP_mo=np.cos(np.array(y_textP)*np.pi/180)-np.cos(np.array(x_textP)*np.pi/180)
          pdiffP_mo=(np.cos(np.array(y_textP)*np.pi/180)-np.cos(np.array(x_textP))*np.pi/180)*100.0/(np.cos(np.array(x_textP)*np.pi/180)+np.cos(10*np.pi/180))
   
       diffP_mo_Oarea=diffP_mo[0:howmany_Oarea]
       diffP_mo_Marea=diffP_mo[howmany_Oarea:howmany_Oarea+howmany_Marea]
       diffP_mo_Garea=diffP_mo[howmany_Oarea+howmany_Marea:howmany_Oarea+howmany_Marea+howmany_Garea]
       diffP_mo_Aarea=diffP_mo[howmany_Oarea+howmany_Marea+howmany_Garea:howmany_Oarea+howmany_Marea+howmany_Garea+howmany_Aarea]
       diffP_mo_Tarea=diffP_mo[howmany_Oarea+howmany_Marea+howmany_Garea+howmany_Aarea:howmany_Oarea+howmany_Marea+howmany_Garea+howmany_Aarea+howmany_Tarea]
       diffP_mo_TAarea=diffP_mo[howmany_Oarea+howmany_Marea+howmany_Garea+howmany_Aarea+howmany_Tarea:howmany_Oarea+howmany_Marea+howmany_Garea+howmany_Aarea+howmany_Tarea+howmany_TAarea]
       diffP_mo_TGarea=diffP_mo[howmany_Oarea+howmany_Marea+howmany_Garea+howmany_Aarea+howmany_Tarea+howmany_TAarea:howmany_Oarea+howmany_Marea+howmany_Garea+howmany_Aarea+howmany_Tarea+howmany_TAarea+howmany_TGarea]
       diffP_mo_Earea=diffP_mo[howmany_Oarea+howmany_Marea+howmany_Garea+howmany_Aarea+howmany_Tarea+howmany_TAarea+howmany_TGarea:howmany_Oarea+howmany_Marea+howmany_Garea+howmany_Aarea+howmany_Tarea+howmany_TAarea+howmany_TGarea+howmany_Earea]
   
          
       mean_diffP=np.mean(diffP_mo)

       # STATISTICS x TAB PHA
       max_diffP=np.max(abs(diffP_mo))
       meanabs_diffP=np.mean(abs(diffP_mo))
       wheremax_PhaAbs=TOT_tg_name[np.argmax(abs(diffP_mo))]
       wheremax_PhaAbs_col=TOT_tg_orcol[np.argmax(abs(diffP_mo))]
       perc95_diffP=np.percentile(abs(diffP_mo),95)
       mean_pdiffP=np.mean(abs(pdiffP_mo))
       meanabs_pdiffP=np.mean(abs(pdiffP_mo))
       max_pdiffP=np.max(abs(pdiffP_mo))
       perc95_pdiffP=np.percentile(abs(pdiffP_mo),95)

       # Table for TEX Pha
       print(comp,' & ',round(max_diffP,2),'$^{\circ}$ ({\color{',wheremax_PhaAbs_col,'}{',wheremax_PhaAbs,'}})& ',round(perc95_diffP,2),'$^{\circ}$ &',round(meanabs_diffP,1),'$^{\circ}$ \\\\ '+'\n', file=Pha_file)
       print('\hline'+'\n', file=Pha_file)
       print(' '+'\n', file=Pha_file)
   
       # Plot 1 ( A and P x stz )
       plt.figure(figsize=(24,10))
       plt.rc('font', size=12)
       # Amp 
       plt.subplot(2,1,1)
       plt.xticks(fontsize=12)
       plt.errorbar(TOT_tg_lab,TOT_A_obs,yerr=np.array(TOT_EA_obs),fmt='-s', color = 'black' ,label = 'Obs')
       
       if tpxo_flag == 1:
          TPXO_AMP=globals()['TPXO_'+comp]
          TPXO_PHA=globals()['TPXO_P_'+comp]
          TPXO_AMP=np.multiply(TPXO_AMP,100) # Want cm not m!
          plt.plot(TPXO_AMP, '--v', color = 'black' ,label = 'TPXO')
   
       if howmany_Oarea != 0:
          plt.errorbar(TOT_tg_lab_Oarea[0:howmany_Oarea], np.array(TOT_A_mod_Oarea[0:howmany_Oarea]),yerr=np.array(TOT_EA_mod_Oarea[0:howmany_Oarea]),fmt='-o', color=subregions_color[4],label = 'Mod (Other areas)')
       if howmany_Marea != 0:
          plt.errorbar(TOT_tg_lab_Marea, np.array(TOT_A_mod_Marea),yerr=np.array(TOT_EA_mod_Marea),fmt='-o', color=subregions_color[2], label = 'Mod (Messina area)')
       if howmany_Garea != 0:
          plt.errorbar(TOT_tg_lab_Garea, np.array(TOT_A_mod_Garea),yerr=np.array(TOT_EA_mod_Garea),fmt='-o', color=subregions_color[0], label = 'Mod (Gibraltar area)')
       if howmany_Aarea != 0:
          plt.errorbar(TOT_tg_lab_Aarea, np.array(TOT_A_mod_Aarea),yerr=np.array(TOT_EA_mod_Aarea),fmt='-o', color=subregions_color[1], label = 'Mod (Adriatic area)')
       if howmany_Tarea != 0:
          plt.errorbar(TOT_tg_lab_Tarea, np.array(TOT_A_mod_Tarea),yerr=np.array(TOT_EA_mod_Tarea),fmt='-o', color=subregions_color[7], label = 'Mod (Atlantic Box - Portugal)')
       if howmany_TAarea != 0:
          plt.errorbar(TOT_tg_lab_TAarea, np.array(TOT_A_mod_TAarea),yerr=np.array(TOT_EA_mod_TAarea),fmt='-o', color=subregions_color[5], label = 'Mod (Atlantic Box - Biscay Bay)')
       if howmany_TGarea != 0:
          plt.errorbar(TOT_tg_lab_TGarea, np.array(TOT_A_mod_TGarea),yerr=np.array(TOT_EA_mod_TGarea),fmt='-o', color=subregions_color[6], label = 'Mod (Atlantic Box - Gibraltar strait)')
       if howmany_Earea != 0:
          plt.errorbar(TOT_tg_lab_Earea, np.array(TOT_A_mod_Earea),yerr=np.array(TOT_EA_mod_Earea),fmt='-o', color=subregions_color[3], label = 'Mod (East Med area)')
       plt.title(comp+' Amplitude [cm] ')
       plt.legend( loc='upper left',fontsize = 'large' )
       plt.grid ()
       plt.ylabel ('Amplitude [cm]')
       plt.ylim(bottom=0.0)
       # Amp diff
       plt.subplot(2,1,2)
       plt.xticks(fontsize=12)
       if howmany_Oarea != 0:
          plt.errorbar(TOT_tg_lab_Oarea[0:howmany_Oarea], diffA_mo_Oarea,yerr=np.sqrt(np.array(TOT_EA_mod_Oarea[0:howmany_Oarea])**2+np.array(TOT_EA_obs_Oarea[0:howmany_Oarea])**2),fmt='-o', color=subregions_color[4],label = 'Mod - Obs (Other areas)')
       if howmany_Marea != 0:
          plt.errorbar(TOT_tg_lab_Marea, diffA_mo_Marea,yerr=np.sqrt(np.array(TOT_EA_mod_Marea)**2+np.array(TOT_EA_obs_Marea)**2),fmt='-o', color=subregions_color[2], label = 'Mod - Obs (Messina area)')
       if howmany_Garea != 0:
          plt.errorbar(TOT_tg_lab_Garea, diffA_mo_Garea,yerr=np.sqrt(np.array(TOT_EA_mod_Garea)**2+np.array(TOT_EA_obs_Garea)**2),fmt='-o', color=subregions_color[0] ,label = 'Mod - Obs (Gibraltar area)')
       if howmany_Aarea != 0:
          plt.errorbar(TOT_tg_lab_Aarea, diffA_mo_Aarea,yerr=np.sqrt(np.array(TOT_EA_mod_Aarea)**2+np.array(TOT_EA_obs_Aarea)**2),fmt='-o', color=subregions_color[1],label = 'Mod - Obs (Adriatic area)')
       if howmany_Tarea != 0:
          plt.errorbar(TOT_tg_lab_Tarea, diffA_mo_Tarea,yerr=np.sqrt(np.array(TOT_EA_mod_Tarea)**2+np.array(TOT_EA_obs_Tarea)**2),fmt='-o', color=subregions_color[7], label = 'Mod - Obs (Atlantic Box - Portugal)')
       if howmany_TAarea != 0:
          plt.errorbar(TOT_tg_lab_TAarea, diffA_mo_TAarea,yerr=np.sqrt(np.array(TOT_EA_mod_TAarea)**2+np.array(TOT_EA_obs_TAarea)**2),fmt='-o', color=subregions_color[5], label = 'Mod - Obs (Atlantic Box - Biscay Bay)')
       if howmany_TGarea != 0:
          plt.errorbar(TOT_tg_lab_TGarea, diffA_mo_TGarea,yerr=np.sqrt(np.array(TOT_EA_mod_TGarea)**2+np.array(TOT_EA_obs_TGarea)**2),fmt='-o', color=subregions_color[6], label = 'Mod - Obs (Atlantic Box  - Gibraltar strait)')
       if howmany_Earea != 0:
          plt.errorbar(TOT_tg_lab_Earea, diffA_mo_Earea,yerr=np.sqrt(np.array(TOT_EA_mod_Earea)**2+np.array(TOT_EA_obs_Earea)**2),fmt='-o', color=subregions_color[3],label = 'Mod - Obs (East Med area)')
       if tpxo_flag == 1 and flag_15stats==0 :
          plt.plot(TPXO_AMP-TOT_A_obs, '--v', color = 'black' ,label = 'TPXO-OBS')
       plt.title(comp+' Amplitude diff (Mod - Obs) [cm] ')
       plt.grid ()
       plt.axhline(y=0, color='black')
       App_mean_diffA=round(mean_diffA,2)
       plt.axhline(y=mean_diffA, color='grey', linestyle='dashed' ,label = '<Amp Diff>='+str(App_mean_diffA)+'[cm]')
       plt.ylabel ('Mod - Obs [cm]')
       plt.legend( loc='upper left',fontsize = 'large' )
       # Add %
       i=0
       for word in pdiffA_mo:
            word_approx=round(word,1)
            string_to=str(word_approx)+'%'
            plt.text(TOT_tg_lab[i],diffA_mo[i]+.03,string_to,fontsize=12,color = 'black')
            i=i+1
       
       if where_box=='Med' and tpxo_flag == 1 :
          plt.savefig(workdir_path+comp+'_Atpxo.jpg')
       elif where_box=='Med' and tpxo_flag == 0 :
          plt.savefig(workdir_path+comp+'_A.jpg')
       elif where_box=='AtlBox':
          plt.savefig(workdir_path+comp+'_A_AB.jpg')
       plt.clf()
          ###
       # Pha
       if cos_pha == 0:
   
        plt.figure(figsize=(24,10))
        plt.rc('font', size=12)
        plt.subplot(2,1,1)
        plt.xticks(fontsize=12)
        plt.errorbar(TOT_tg_lab,TOT_P_obs,yerr=np.array(TOT_EP_obs),fmt='-s', color = 'black' , label = 'Obs')
        if tpxo_flag == 1:
           plt.plot(TPXO_PHA, '--v', color = 'black' ,label = 'TPXO')
    
        if howmany_Oarea != 0:
           plt.errorbar(TOT_tg_lab_Oarea[0:howmany_Oarea], np.array(TOT_P_mod_Oarea[0:howmany_Oarea]),yerr=np.array(TOT_EP_mod_Oarea[0:howmany_Oarea]),fmt='-o', color=subregions_color[4],label = 'Mod (Other areas)')
        if howmany_Marea != 0:
           plt.errorbar(TOT_tg_lab_Marea, np.array(TOT_P_mod_Marea),yerr=np.array(TOT_EP_mod_Marea),fmt='-o',color=subregions_color[2], label = 'Mod (Messina area)')
        if howmany_Garea != 0:
           plt.errorbar(TOT_tg_lab_Garea, np.array(TOT_P_mod_Garea),yerr=np.array(TOT_EP_mod_Garea),fmt='-o', color=subregions_color[0], label = 'Mod (Gibraltar area)')
        if howmany_Aarea != 0:
           plt.errorbar(TOT_tg_lab_Aarea, np.array(TOT_P_mod_Aarea),yerr=np.array(TOT_EP_mod_Aarea),fmt='-o', color=subregions_color[1],label = 'Mod (Adriatic area)')
        if howmany_Tarea != 0:
           plt.errorbar(TOT_tg_lab_Tarea, np.array(TOT_P_mod_Tarea),yerr=np.array(TOT_EP_mod_Tarea),fmt='-o', color=subregions_color[7], label = 'Mod (Atlantic Box - Portugal)')
        if howmany_TAarea != 0:
           plt.errorbar(TOT_tg_lab_TAarea, np.array(TOT_P_mod_TAarea),yerr=np.array(TOT_EP_mod_TAarea),fmt='-o', color=subregions_color[5],label = 'Mod (Atlantic Box - Biscay Bay)')
        if howmany_TGarea != 0:
           plt.errorbar(TOT_tg_lab_TGarea, np.array(TOT_P_mod_TGarea),yerr=np.array(TOT_EP_mod_TGarea),fmt='-o', color=subregions_color[6], label = 'Mod (Atlantic Box - Gibraltar strait)')
        if howmany_Earea != 0:
           plt.errorbar(TOT_tg_lab_Earea, np.array(TOT_P_mod_Earea),yerr=np.array(TOT_EP_mod_Earea),fmt='-o', color=subregions_color[3],label = 'Mod (East Med area)')
        plt.title(comp+' Phase [deg] ')
        plt.grid ()
        plt.ylabel ('Phase [deg]')
        plt.ylim(-50.0, 450.0)
        plt.legend( loc='upper left',fontsize = 'large')
        # Pha diff
        plt.subplot(2,1,2)
        plt.xticks(fontsize=12)
        if howmany_Oarea != 0:
           plt.errorbar(TOT_tg_lab_Oarea[0:howmany_Oarea], diffP_mo_Oarea,yerr=np.sqrt(np.array(TOT_EP_mod_Oarea[0:howmany_Oarea])**2+np.array(TOT_EP_obs_Oarea[0:howmany_Oarea])**2),fmt='-o', color=subregions_color[4],label = 'Mod - Obs (Other areas)')
        if howmany_Marea != 0:
           plt.errorbar(TOT_tg_lab_Marea, diffP_mo_Marea,yerr=np.sqrt(np.array(TOT_EP_mod_Marea)**2+np.array(TOT_EP_obs_Marea)**2),fmt='-o', color=subregions_color[2], label = 'Mod - Obs (Messina area)')
        if howmany_Garea != 0:
           plt.errorbar(TOT_tg_lab_Garea, diffP_mo_Garea,yerr=np.sqrt(np.array(TOT_EP_mod_Garea)**2+np.array(TOT_EP_obs_Garea)**2),fmt='-o', color=subregions_color[0] ,label = 'Mod - Obs (Gibraltar area)')
        if howmany_Aarea != 0:
           plt.errorbar(TOT_tg_lab_Aarea, diffP_mo_Aarea,yerr=np.sqrt(np.array(TOT_EP_mod_Aarea)**2+np.array(TOT_EP_obs_Aarea)**2),fmt='-o', color=subregions_color[1],label = 'Mod - Obs (Adriatic area)')
        if howmany_Tarea != 0:
           plt.errorbar(TOT_tg_lab_Tarea, diffP_mo_Tarea,yerr=np.sqrt(np.array(TOT_EP_mod_Tarea)**2+np.array(TOT_EP_obs_Tarea)**2),fmt='-o', color=subregions_color[7], label = 'Mod - Obs (Atlantic Box - Portugal)')
        if howmany_TAarea != 0:
           plt.errorbar(TOT_tg_lab_TAarea, diffP_mo_TAarea,yerr=np.sqrt(np.array(TOT_EP_mod_TAarea)**2+np.array(TOT_EP_obs_TAarea)**2),fmt='-o', color=subregions_color[5],label = 'Mod - Obs (Atlantic Box - Biscay Bay)')
        if howmany_TGarea != 0:
           plt.errorbar(TOT_tg_lab_TGarea, diffP_mo_TGarea,yerr=np.sqrt(np.array(TOT_EP_mod_TGarea)**2+np.array(TOT_EP_obs_TGarea)**2),fmt='-o', color=subregions_color[6], label = 'Mod - Obs (Atlantic Box - Gibraltar strait)')
        if howmany_Earea != 0:
           plt.errorbar(TOT_tg_lab_Earea, diffP_mo_Earea,yerr=np.sqrt(np.array(TOT_EP_mod_Earea)**2+np.array(TOT_EP_obs_Earea)**2),fmt='-o', color=subregions_color[3],label = 'Mod - Obs (East Med area)')
        if tpxo_flag == 1 and flag_15stats==0 :
           plt.plot(TPXO_PHA-TOT_P_obs, '--v', color = 'black' ,label = 'TPXO-OBS')
        plt.title(comp+' Phase diff (Mod - Obs) [deg] ')
        plt.grid ()
        plt.ylabel ('Mod - Obs [deg]')
        App_mean_diffP=round(mean_diffP,2)
        plt.axhline(y=mean_diffP, color='gray' , linestyle='dashed' ,label = '<Pha Diff>='+str(App_mean_diffP)+'[deg]')
        plt.axhline(y=0, color='black')
        plt.legend( loc='upper left',fontsize = 'large' )
        # Add %
        i=0
        for word in pdiffP_mo:
             word_approx=round(word,1)
             string_to=str(word_approx)+'%'
             plt.text(TOT_tg_lab[i],diffP_mo[i]+.03,string_to,fontsize=12,color = 'black')
             i=i+1
        #
        if where_box=='Med' and tpxo_flag == 1 :
          plt.savefig(workdir_path+comp+'_Ptpxo.jpg')
        elif where_box=='Med' and tpxo_flag == 0 :
          plt.savefig(workdir_path+comp+'_P.jpg')
        elif where_box=='AtlBox':
          plt.savefig(workdir_path+comp+'_P_AB.jpg')
        plt.clf()
   
       elif cos_pha == 1:
        # WARNING: in this case should be added the Pha error..
        plt.figure(figsize=(24,10))
        plt.rc('font', size=12)
        plt.subplot(2,1,1)
        plt.xticks(fontsize=12)
        plt.plot(np.cos(TOT_P_obs*np.pi/180), '-s', color = 'black' , label = 'Obs')
        if tpxo_flag == 1:
           plt.plot(np.cos(TPXO_PHA), '--v', color = 'black' ,label = 'TPXO')
    
        if howmany_Oarea != 0:
           plt.plot(TOT_tg_lab_Oarea[0:howmany_Oarea], np.cos(np.array(TOT_P_mod_Oarea[0:howmany_Oarea])*np.pi/180), '-o', color=subregions_color[4],label = 'Mod (Other areas)')
        if howmany_Marea != 0:
           plt.plot(TOT_tg_lab_Marea, np.cos(np.array(TOT_P_mod_Marea)*np.pi/180), '-o',color=subregions_color[2], label = 'Mod (Messina area)')
        if howmany_Garea != 0:
           plt.plot(TOT_tg_lab_Garea, np.cos(np.array(TOT_P_mod_Garea)*np.pi/180), '-o', color=subregions_color[0], label = 'Mod (Gibraltar area)')
        if howmany_Aarea != 0:
           plt.plot(TOT_tg_lab_Aarea, np.cos(np.array(TOT_P_mod_Aarea)*np.pi/180), '-o', color=subregions_color[1],label = 'Mod (Adriatic area)')
        if howmany_Tarea != 0:
           plt.plot(TOT_tg_lab_Tarea, np.cos(np.array(TOT_P_mod_Tarea)*np.pi/180), '-o', color=subregions_color[7], label = 'Mod (Atlantic Box - Portugal)')
        if howmany_TAarea != 0:
           plt.plot(TOT_tg_lab_TAarea, np.cos(np.array(TOT_P_mod_TAarea)*np.pi/180), '-o', color=subregions_color[5],label = 'Mod (Atlantic Box - Biscay Bay)')
        if howmany_TGarea != 0:
           plt.plot(TOT_tg_lab_TGarea, np.cos(np.array(TOT_P_mod_TGarea)*np.pi/180), '-o', color=subregions_color[6], label = 'Mod (Atlantic Box - Gibraltar strait)')
        if howmany_Earea != 0:
           plt.plot(TOT_tg_lab_Earea, np.cos(np.array(TOT_P_mod_Earea)*np.pi/180), '-o', color=subregions_color[3],label = 'Mod (East Med area)')
        plt.title(comp+' cos(Phase) ')
        plt.grid ()
        plt.ylabel ('cos(Phase)')
        plt.ylim(-1, 1)
        plt.legend( loc='upper left',fontsize = 'large')
        # Pha diff
        plt.subplot(2,1,2)
        plt.xticks(fontsize=12)
        if howmany_Oarea != 0:
           plt.plot(TOT_tg_lab_Oarea[0:howmany_Oarea], diffP_mo_Oarea, '-o', color=subregions_color[4],label = 'Mod - Obs (Other areas)')
        if howmany_Marea != 0:
           plt.plot(TOT_tg_lab_Marea, diffP_mo_Marea, '-o', color=subregions_color[2], label = 'Mod - Obs (Messina area)')
        if howmany_Garea != 0:
           plt.plot(TOT_tg_lab_Garea, diffP_mo_Garea, '-o', color=subregions_color[0], label = 'Mod - Obs (Gibraltar area)')
        if howmany_Aarea != 0:
           plt.plot(TOT_tg_lab_Aarea, diffP_mo_Aarea, '-o', color=subregions_color[1], label = 'Mod - Obs (Adriatic area)')
        if howmany_Tarea != 0:
           plt.plot(TOT_tg_lab_Tarea, diffP_mo_Tarea, '-o', color=subregions_color[7], label = 'Mod - Obs (Atlantic Box - Portugal)')
        if howmany_TAarea != 0:
           plt.plot(TOT_tg_lab_TAarea, diffP_mo_TAarea, '-o', color=subregions_color[5], label = 'Mod - Obs (Atlantic Box - Biscay Bay)')
        if howmany_TGarea != 0:
           plt.plot(TOT_tg_lab_TGarea, diffP_mo_TGarea, '-o', color=subregions_color[6], label = 'Mod - Obs (Atlantic Box - Gibraltar strait)')
        if howmany_Earea != 0:
           plt.plot(TOT_tg_lab_Earea, diffP_mo_Earea, '-o', color=subregions_color[3], label = 'Mod - Obs (East Med area)')
        if tpxo_flag == 1:
           plt.plot(np.cos(np.array(TPXO_PHA)*np.pi/180)-np.cos(np.array(TOT_P_obs)*np.pi/180), '--v', color = 'black' ,label = 'TPXO-OBS')
        plt.title(comp+' cos(Phase) diff (Mod - Obs) ')
        plt.grid ()
        plt.ylabel ('Mod - Obs')
        App_mean_diffP=round(mean_diffP,2)
        plt.axhline(y=mean_diffP, color='gray' , linestyle='dashed' ,label = '<Pha Diff>='+str(App_mean_diffP)+'[deg]')
        plt.axhline(y=0, color='black')
        plt.legend( loc='upper left',fontsize = 'large' )
        # Add %
        i=0
        for word in pdiffP_mo:
             word_approx=round(word,1)
             string_to=str(word_approx)+'%'
             i=i+1
        #
        if where_box=='Med' and tpxo_flag == 1 :
          plt.savefig(workdir_path+comp+'_cosPtpxo.jpg')
        elif where_box=='Med' and tpxo_flag == 0 :
          plt.savefig(workdir_path+comp+'_cosP.jpg')
        elif where_box=='AtlBox':
          plt.savefig(workdir_path+comp+'_cosP_AB.jpg')
        plt.clf()
   
   
   
   ######### LINEAR REGRESSIONS
   
      # Plot ( Lin reg mod Vs obs x A and P x stz )
       if cos_pha == 0:
          plt.figure(figsize=(6,12))
       elif cos_pha == 1:
          plt.figure(figsize=(7,12))
       plt.rc('font', size=12)
       #
       plt.subplot(2,1,1)
       plt.title(comp+' Amplitude [cm] ')
       plt.grid ()
       # Arrays defn
       x_text=[]
       y_text=[]
       y_text2=[]
       x_text=np.array(TOT_A_obs)
       y_text=np.array(TOT_A_mod)
       top=np.maximum(x_text,y_text)
       top=max(top[:])
       # Linear regression
       slope, intercept, r_value, p_value, std_err = stats.linregress(x_text,y_text)
       m_A=[]
       q_A=[]
       fitted_A=[]
       cov_A=[]
       def line_A(x, m_A, q_A):
           return (m_A*x+q_A)
       fitted_A, cov_A = curve_fit(line_A,x_text[:],y_text[:])
       m_A=fitted_A[0]
       q_A=fitted_A[1]
       rx=np.linspace(0.0,top)
       retta_A=m_A*rx+q_A
       m_A_approx=round(m_A,2)
       perr = np.abs(np.diag(cov_A))
       m_Ae_approx=round(perr[0],2)
       r_A=round(r_value,2)
       lr_leg_str='( Slope='+str(m_A_approx)+'; R2='+str(r_A)+')'
   
       # Arrays defn
       x_text=[]
       x_text=np.array(TOT_A_obs)
       # Plot points
       if howmany_Oarea != 0:
          plt.errorbar(np.array(TOT_A_obs_Oarea), np.array(TOT_A_mod_Oarea),xerr=np.array(TOT_EA_obs_Oarea),yerr=np.array(TOT_EA_mod_Oarea),fmt='go', label = 'Other areas')
       if howmany_Marea != 0:
          plt.errorbar(np.array(TOT_A_obs_Marea), np.array(TOT_A_mod_Marea),xerr=np.array(TOT_EA_obs_Marea),yerr=np.array(TOT_EA_mod_Marea),fmt='o',color=subregions_color[2], label = 'Messina area')
       if howmany_Garea != 0:
          plt.errorbar(np.array(TOT_A_obs_Garea), np.array(TOT_A_mod_Garea),xerr=np.array(TOT_EA_obs_Garea),yerr=np.array(TOT_EA_mod_Garea),fmt='ro', label = 'Gibraltar area')
       if howmany_Aarea != 0:
          plt.errorbar(np.array(TOT_A_obs_Aarea), np.array(TOT_A_mod_Aarea),xerr=np.array(TOT_EA_obs_Aarea),yerr=np.array(TOT_EA_mod_Aarea),fmt='bo', label = 'Adriatic area')
       if howmany_Tarea != 0:
          plt.errorbar(np.array(TOT_A_obs_Tarea), np.array(TOT_A_mod_Tarea),xerr=np.array(TOT_EA_obs_Tarea),yerr=np.array(TOT_EA_mod_Tarea),fmt='o',color=subregions_color[7], label = 'Atlantic Box - Portugal')
       if howmany_TAarea != 0:
          plt.errorbar(np.array(TOT_A_obs_TAarea), np.array(TOT_A_mod_TAarea),xerr=np.array(TOT_EA_obs_TAarea),yerr=np.array(TOT_EA_mod_TAarea),fmt='co', label = 'Atlantic Box - Biscay Bay')
       if howmany_TGarea != 0:
          plt.errorbar(np.array(TOT_A_obs_TGarea), np.array(TOT_A_mod_TGarea),xerr=np.array(TOT_EA_obs_TGarea),yerr=np.array(TOT_EA_mod_TGarea),fmt='o', color=subregions_color[6], label = 'Atlantic Box - Gibraltar strait')
       if howmany_Earea != 0:
          plt.errorbar(np.array(TOT_A_obs_Earea), np.array(TOT_A_mod_Earea),xerr=np.array(TOT_EA_obs_Earea),yerr=np.array(TOT_EA_mod_Earea),fmt='mo', label = 'East Med area')
       # Ax settings
       plt.xlabel ('OBS Amplitude [cm]')
       plt.ylabel ('MOD Amplitude [cm]')
       bottom_ax, top_ax = plt.xlim()
       bottom_ay, top_ay = plt.ylim()
       top_a=np.maximum(top_ax,top_ay)
       bottom_a=np.minimum(bottom_ax,bottom_ay)
       plt.ylim(bottom_a, top_a)
       plt.xlim(bottom_a, top_a)
       # Plot Lines
       plt.plot(rx,retta_A,color = 'red',label=lr_leg_str)
       plt.plot([bottom_a,top_a], [bottom_a,top_a], 'k-', color = 'black')
       # Legend
       plt.legend( loc='upper left' )
       # point label (stn names)
       if linreg_name_flag != 0: 
         i=0
         for word in TOT_tg_lab:
             plt.text(x_text[i]+.03,y_text[i]+.03,word,fontsize=12,color = 'black')
             i=i+1
   
       plt.subplot(2,1,2)
       ### Pha linear reg
       if cos_pha == 0:
        plt.title(comp+' Phase [deg] ')
        plt.grid ()
        if where_box == 'AtlBox':
           plt.ylim(0.0, 360.0)
        else:
           plt.ylim(-50.0, 400.0)
        plt.xlim(-50.0, 400.0)
        plt.plot([-50.0, 400.0], [-50.0, 400.0], 'k-', color = 'black')
        x_text=[]
        y_text=[]
        x_text=np.array(TOT_P_obs)
        y_text=np.array(TOT_P_mod)
        slopeP, interceptP, r_valueP, p_valueP, std_errP = stats.linregress(x_text,y_text)
        rx=np.linspace(-50.0,400.0)
        retta_P=slopeP*rx+interceptP
        lr_leg_str='( Slope='+str(round(slopeP,2))+'; R2='+str(round(r_valueP,2))+')'
        plt.plot(rx,retta_P,color = 'red',label=lr_leg_str)
        if howmany_Oarea != 0:
           plt.errorbar(np.array(TOT_P_obs_Oarea), np.array(TOT_P_mod_Oarea),xerr=np.array(TOT_EP_obs_Oarea),yerr=np.array(TOT_EP_mod_Oarea),fmt='go', label = 'Other areas')
        if howmany_Marea != 0:
           plt.errorbar(np.array(TOT_P_obs_Marea), np.array(TOT_P_mod_Marea),xerr=np.array(TOT_EP_obs_Marea),yerr=np.array(TOT_EP_mod_Marea),fmt='o',color=subregions_color[2], label = 'Messina area')
        if howmany_Garea != 0:
           plt.errorbar(np.array(TOT_P_obs_Garea), np.array(TOT_P_mod_Garea),xerr=np.array(TOT_EP_obs_Garea),yerr=np.array(TOT_EP_mod_Garea),fmt='ro', label = 'Gibraltar area')
        if howmany_Aarea != 0:
           plt.errorbar(np.array(TOT_P_obs_Aarea), np.array(TOT_P_mod_Aarea),xerr=np.array(TOT_EP_obs_Aarea),yerr=np.array(TOT_EP_mod_Aarea),fmt='bo', label = 'Adriatic area')
        if howmany_Tarea != 0:
           plt.errorbar(np.array(TOT_P_obs_Tarea), np.array(TOT_P_mod_Tarea),xerr=np.array(TOT_EP_obs_Tarea),yerr=np.array(TOT_EP_mod_Tarea),fmt='o',color=subregions_color[7], label = 'Atlantic Box - Portugal')
        if howmany_Tarea != 0:
           plt.errorbar(np.array(TOT_P_obs_TAarea), np.array(TOT_P_mod_TAarea),xerr=np.array(TOT_EP_obs_TAarea),yerr=np.array(TOT_EP_mod_TAarea),fmt='co', label = 'Atlantic Box - Biscay Bay')
        if howmany_TGarea != 0:
           plt.errorbar(np.array(TOT_P_obs_TGarea), np.array(TOT_P_mod_TGarea),xerr=np.array(TOT_EP_obs_TGarea),yerr=np.array(TOT_EP_mod_TGarea),fmt='o', color=subregions_color[6], label = 'Atlantic Box - Gibraltar strait')
        if howmany_Earea != 0:
           plt.errorbar(np.array(TOT_P_obs_Earea), np.array(TOT_P_mod_Earea),xerr=np.array(TOT_EP_obs_Earea),yerr=np.array(TOT_EP_mod_Earea),fmt='mo', label = 'East Med area') 
    
        # Axes
        plt.xlabel ('OBS Phase [deg]')
        plt.ylabel ('MOD Phase [deg]')
        # Legend
        plt.legend( loc='upper left' )
        if linreg_name_flag != 0:
           i=0
           for word in TOT_tg_lab:
               plt.text(x_text[i]+.03,y_text[i]+.03,word,fontsize=12,color = 'black')
               i=i+1
   
        if where_box=='Med':
          plt.savefig(workdir_path+comp+'_lr.jpg')
        elif where_box=='AtlBox':
          plt.savefig(workdir_path+comp+'_lr_AB.jpg')
        plt.clf()
   
       elif cos_pha == 1:
        plt.title(comp+' cos(Phase)')
        plt.grid ()
        plt.ylim(-1, 1)
        plt.xlim(-1, 1)
        plt.plot([-1, 1], [-1, 1], 'k-', color = 'black')
        x_text=[]
        y_text=[]
        x_text=np.cos(np.array(TOT_P_obs)*np.pi/180)
        y_text=np.cos(np.array(TOT_P_mod)*np.pi/180)
        slopeP, interceptP, r_valueP, p_valueP, std_errP = stats.linregress(x_text,y_text)
        rx=np.linspace(-1,1)
        retta_P=slopeP*rx+interceptP
        lr_leg_str='( Slope='+str(round(slopeP,2))+'; R2='+str(round(r_valueP,2))+')'
        plt.plot(rx,retta_P,color = 'red',label=lr_leg_str)
        if howmany_Oarea != 0:
           plt.plot(np.cos(np.array(TOT_P_obs_Oarea)*np.pi/180), np.cos(np.array(TOT_P_mod_Oarea)*np.pi/180), 'go', label = 'Other areas')
        if howmany_Marea != 0:
           plt.plot(np.cos(np.array(TOT_P_obs_Marea)*np.pi/180), np.cos(np.array(TOT_P_mod_Marea)*np.pi/180), 'o',color=subregions_color[2], label = 'Messina area')
        if howmany_Garea != 0:
           plt.plot(np.cos(np.array(TOT_P_obs_Garea)*np.pi/180), np.cos(np.array(TOT_P_mod_Garea)*np.pi/180), 'ro', label = 'Gibraltar area')
        if howmany_Aarea != 0:
           plt.plot(np.cos(np.array(TOT_P_obs_Aarea)*np.pi/180), np.cos(np.array(TOT_P_mod_Aarea)*np.pi/180), 'bo', label = 'Adriatic area')
        if howmany_Tarea != 0:
           plt.plot(np.cos(np.array(TOT_P_obs_Tarea)*np.pi/180), np.cos(np.array(TOT_P_mod_Tarea)*np.pi/180), 'co', label = 'Atlantic Box - Portugal')
        if howmany_TAarea != 0:
           plt.plot(np.cos(np.array(TOT_P_obs_TAarea)*np.pi/180), np.cos(np.array(TOT_P_mod_TAarea)*np.pi/180), 'o',color=subregions_color[7], label = 'Atlantic Box - Biscay Bay')
        if howmany_TGarea != 0:
           plt.plot(np.cos(np.array(TOT_P_obs_TGarea)*np.pi/180), np.cos(np.array(TOT_P_mod_TGarea)*np.pi/180), 'o', color=subregions_color[6], label = 'Atlantic Box - Gibraltar strait')
        if howmany_Earea != 0:
           plt.plot(np.cos(np.array(TOT_P_obs_Earea)*np.pi/180), np.cos(np.array(TOT_P_mod_Earea)*np.pi/180), 'mo', label = 'East Med area')
    
        # Axes
        plt.xlabel ('OBS cos(Phase)')
        plt.ylabel ('MOD cos(Phase)')
        # Legend
        plt.legend( loc='lower right' )
        if linreg_name_flag != 0:
           i=0
           for word in TOT_tg_lab:
               plt.text(x_text[i]+.03,y_text[i]+.03,word,fontsize=12,color = 'black')
               i=i+1
   
        if where_box=='Med':
          plt.savefig(workdir_path+comp+'_lr_cos.jpg')
        elif where_box=='AtlBox':
          plt.savefig(workdir_path+comp+'_lr_cos_AB.jpg')
        plt.clf()
   
   
       print(comp,' & ',str(m_A_approx),' & ',str(r_A),' & ',str(round(slopeP,2)),' & ',str(round(r_valueP,2)),' \\\\ '+'\n', file=LinReg_file)
   
   
       # Save val in GLOBAL arrays 
   
       # Components names
       GLOB_A_mod[0][comp_idx+1]=comp
       GLOB_P_mod[0][comp_idx+1]=comp
       GLOB_A_obs[0][comp_idx+1]=comp
       GLOB_P_obs[0][comp_idx+1]=comp
   
       d_foreman[0][comp_idx+1]=comp 
   
       # Stations labels
       for nnn_stz in range (0,len(TOT_tg_lab_ord)): 
         GLOB_A_mod[nnn_stz+1][0]=TOT_tg_lab_ord[nnn_stz]
         GLOB_P_mod[nnn_stz+1][0]=TOT_tg_lab_ord[nnn_stz]
         GLOB_A_obs[nnn_stz+1][0]=TOT_tg_lab_ord[nnn_stz]
         GLOB_P_obs[nnn_stz+1][0]=TOT_tg_lab_ord[nnn_stz]
   
         d_foreman[nnn_stz+1][0]=TOT_tg_lab_ord[nnn_stz]
   
       # Put right numbers in the matrix!
       for nnn_AP in range (0,len(TOT_tg_lab_ord)):
         GLOB_A_mod[nnn_AP+1][comp_idx+1]=TOT_A_mod_ord[nnn_AP]
         GLOB_P_mod[nnn_AP+1][comp_idx+1]=TOT_P_mod_ord[nnn_AP]
         GLOB_A_obs[nnn_AP+1][comp_idx+1]=TOT_A_obs_ord[nnn_AP]
         GLOB_P_obs[nnn_AP+1][comp_idx+1]=TOT_P_obs_ord[nnn_AP]
   
   
         # Compute Distances in the complex plane [Foreman et al. 93]
         d_foreman[nnn_AP+1][comp_idx+1]=np.sqrt((TOT_A_obs_ord[nnn_AP]*np.cos((np.pi/180.0)*TOT_P_obs_ord[nnn_AP])-(TOT_A_mod_ord[nnn_AP]*np.cos((np.pi/180.0)*TOT_P_mod_ord[nnn_AP])))**2+((TOT_A_obs_ord[nnn_AP]*np.sin((np.pi/180.0)*TOT_P_obs_ord[nnn_AP])-(TOT_A_mod_ord[nnn_AP]*np.sin((np.pi/180.0)*TOT_P_mod_ord[nnn_AP])))**2))
   
         # Root Mean Square misfits
         RMSm[comp_idx]=RMSm[comp_idx]+(TOT_A_obs_ord[nnn_AP]*np.cos((np.pi/180.0)*TOT_P_obs_ord[nnn_AP])-(TOT_A_mod_ord[nnn_AP]*np.cos((np.pi/180.0)*TOT_P_mod_ord[nnn_AP])))**2+((TOT_A_obs_ord[nnn_AP]*np.sin((np.pi/180.0)*TOT_P_obs_ord[nnn_AP])-(TOT_A_mod_ord[nnn_AP]*np.sin((np.pi/180.0)*TOT_P_mod_ord[nnn_AP])))**2)
   
       RMSm[comp_idx]=np.sqrt((1/(2*len(TOT_tg_lab_ord)))*RMSm[comp_idx]) 
   
       comp_idx=comp_idx+1
   
   
   ###################### END OF THE LOOP ON COMPONENTS
   # Close tables
   Amp_file.close()
   Pha_file.close() 
   LinReg_file.close()
   
   
   ##################### ANALYSIS ON THE WHOLE SET OF COMPONENTS
   print ('Working on the whole set of tidal components..')
   comp=tidal_comp
   num_comp_1=len(tidal_comp)+1
   
   # Write tables
   if tpxo_flag == 1:
    for j in range (1,num_comp_1): # Loop on components
     if tpxo_flag != 0:
      TPXO_AMP=globals()['TPXO_'+comp[j-1]]
      TPXO_PHA=globals()['TPXO_P_'+comp[j-1]]
      for i in range (0,len(TPXO_PHA)):
              while TPXO_PHA[i]<0:
                 TPXO_PHA[i]=TPXO_PHA[i]+360
     if where_box=='Med':
       Tab_filename=workdir_path+'tab_'+comp[j-1]+'.txt'
     elif where_box=='AtlBox':
       Tab_filename=workdir_path+'tab_'+comp[j-1]+'_AB.txt'
   
     # TABLE A_mod/A_obs P_mod/P_obs alpha(ratio Amp) dvec Per stz
     Tab_file = open(Tab_filename,"w")
     Tab_file.write('\\begin{table}'+'\n')
     Tab_file.write('\\footnotesize'+'\n')
     Tab_file.write('\hspace{-2.5cm}'+'\n')
     Tab_file.write('\\begin{tabular}{||c||c|c|c|c|c||c|c|c|c||c|c||}'+'\n')
     Tab_file.write('     \hline'+'\n')
     Tab_file.write('     \hline'+'\n')
     Tab_file.write('{\\bf '+str(comp)+'} & \multicolumn{5}{|c||}{Amplitudes} & \multicolumn{4}{|c||}{Phases} & \multicolumn{2}{|c||}{Vectorial distances} \\\\'+'\n')
     Tab_file.write('     \hline'+'\n')
     Tab_file.write(' TG Name & $A_{mod}$ [cm] & $A_{obs}$ [cm] & $A_{lit}$ [cm] & $A_{tpxo}$ [cm] & $\\alpha$ & $P_{mod}$ [deg] & $P_{obs}$ [deg] & $P_{lit}$ [deg]  & $P_{tpxo}$ [deg]& $d_{mod/obs}$ [cm]  & $d_{lit}$ [cm] \\\\'+'\n')
     Tab_file.write('     \hline'+'\n')
     Tab_file.write('     \hline'+'\n')
   
     for i in range(1,N_stz+1): # Loop on stz
       ALPHA=GLOB_A_mod[i][j]/GLOB_A_obs[i][j]
       print (TOT_tg_name_ord[i-1],'&',end =" ",file=Tab_file)
       if j > 4:
          print(np.round(np.array(GLOB_A_mod[i][j]),1),'&',np.round(np.array(GLOB_A_obs[i][j]),1),'&',np.round(np.array(TPXO_AMP[i-1])*100,1),'&',np.round(ALPHA,2),'&',np.round(np.array(GLOB_P_mod[i][j]),1),'&',np.round(np.array(GLOB_P_obs[i][j]),1),'&',np.round(np.array(TPXO_PHA[i-1]),1),'&',np.round(np.array(d_foreman[i][j]),1),'\\\\'+'\n',file=Tab_file)
       else:
          A_lit1=globals()['TSIMPLIS_'+str(comp[j-1])]
          P_lit1=globals()['TSIMPLIS_P_'+str(comp[j-1])]
          d_lit1=globals()['TSIMPLIS_d_'+str(comp[j-1])]
          A_lit2=globals()['PALMA_'+str(comp[j-1])]
          P_lit2=globals()['PALMA_P_'+str(comp[j-1])]
          d_lit2=globals()['PALMA_d_'+str(comp[j-1])]
          #
          #
          print(np.round(np.array(GLOB_A_mod[i][j]),1),'&',np.round(np.array(GLOB_A_obs[i][j]),1),'&',A_lit1[i-1],'/',A_lit2[i-1],'&',np.round(np.array(TPXO_AMP[i-1])*100,1),'&',np.round(ALPHA,2),'&',np.round(np.array(GLOB_P_mod[i][j]),1),'&',np.round(np.array(GLOB_P_obs[i][j]),1),'&',P_lit1[i-1],'/',P_lit2[i-1],'&',np.round(np.array(TPXO_PHA[i-1]),1),'&',np.round(np.array(d_foreman[i][j]),1),'&',d_lit1[i-1],'/',d_lit2[i-1],'\\\\'+'\n',file=Tab_file)
     Tab_file.write('\hline'+'\n')
     print ('RMSm [cm] &&&&&&&&&&&',np.round(RMSm[j-1],2),'\\\\'+'\n',file=Tab_file)
     Tab_file.close() 
   
   ################### GLOBAL PLOTS ############
   # Colors array
   TOT_color_stz=[]
   for idx_color in range (0,howmany_Oarea):
       TOT_color_stz.append(subregions_color[4])
   for idx_color in range (howmany_Oarea,howmany_Oarea+howmany_Marea):
       TOT_color_stz.append(subregions_color[2])
   for idx_color in range (howmany_Oarea+howmany_Marea,howmany_Oarea+howmany_Marea+howmany_Garea):
       TOT_color_stz.append(subregions_color[0])
   for idx_color in range (howmany_Oarea+howmany_Marea+howmany_Garea,howmany_Oarea+howmany_Marea+howmany_Garea+howmany_Aarea):
       TOT_color_stz.append(subregions_color[1])
   for idx_color in range (howmany_Oarea+howmany_Marea+howmany_Garea+howmany_Aarea,howmany_Oarea+howmany_Marea+howmany_Garea+howmany_Aarea+howmany_Tarea):
       TOT_color_stz.append(subregions_color[7])
   for idx_color in range (howmany_Oarea+howmany_Marea+howmany_Garea+howmany_Aarea+howmany_Tarea,howmany_Oarea+howmany_Marea+howmany_Garea+howmany_Aarea+howmany_Tarea+howmany_TAarea):
       TOT_color_stz.append(subregions_color[5])
   for idx_color in range (howmany_Oarea+howmany_Marea+howmany_Garea+howmany_Aarea+howmany_Tarea+howmany_TAarea,howmany_Oarea+howmany_Marea+howmany_Garea+howmany_Aarea+howmany_Tarea+howmany_TAarea+howmany_TGarea):
       TOT_color_stz.append(subregions_color[6])
   for idx_color in range (howmany_Oarea+howmany_Marea+howmany_Garea+howmany_Aarea+howmany_Tarea+howmany_TAarea+howmany_TGarea,howmany_Oarea+howmany_Marea+howmany_Garea+howmany_Aarea+howmany_Tarea+howmany_TAarea+howmany_TGarea+howmany_Earea):
       TOT_color_stz.append(subregions_color[3])
   
   
   # AMPLITUDES Global plots
   
   plt.figure(figsize=(80,24)) 
   plt.rc('font', size=40) 
   
   labels=[ GLOB_A_mod[i][0] for i in range(1,N_stz+1) ]
   mod_M2 = [ GLOB_A_mod[j][1] for j in range (1,N_stz+1)]
   obs_M2 = [ GLOB_A_obs[j][1] for j in range (1,N_stz+1)]
   mod_S2 = [ GLOB_A_mod[j][2] for j in range (1,N_stz+1)]
   obs_S2 = [ GLOB_A_obs[j][2] for j in range (1,N_stz+1)]
   mod_K1 = [ GLOB_A_mod[j][3] for j in range (1,N_stz+1)]
   obs_K1 = [ GLOB_A_obs[j][3] for j in range (1,N_stz+1)]
   mod_O1 = [ GLOB_A_mod[j][4] for j in range (1,N_stz+1)]
   obs_O1 = [ GLOB_A_obs[j][4] for j in range (1,N_stz+1)]
   #
   mod_N2 = [ GLOB_A_mod[j][5] for j in range (1,N_stz+1)]
   obs_N2 = [ GLOB_A_obs[j][5] for j in range (1,N_stz+1)]
   mod_P1 = [ GLOB_A_mod[j][6] for j in range (1,N_stz+1)]
   obs_P1 = [ GLOB_A_obs[j][6] for j in range (1,N_stz+1)]
   mod_Q1 = [ GLOB_A_mod[j][7] for j in range (1,N_stz+1)]
   obs_Q1 = [ GLOB_A_obs[j][7] for j in range (1,N_stz+1)]
   mod_K2 = [ GLOB_A_mod[j][8] for j in range (1,N_stz+1)]
   obs_K2 = [ GLOB_A_obs[j][8] for j in range (1,N_stz+1)]
   
   x = np.arange(len(labels))  # the label locations
   width = 0.45  # the width of the bars
   
   fig,ax=plt.subplots( figsize=(80,24)) 
   #
   rects1_M2 = ax.bar(x - width/2, mod_M2, width-0.05, color='#1f77b4', label='M2')
   rects2_M2 = ax.bar(x + width/2, obs_M2, width-0.05, color='#1f77b4')
   
   topbottom_M=mod_M2
   topbottom_O=obs_M2
   rects1_S2 = ax.bar(x - width/2, mod_S2, width-0.05, bottom=topbottom_M, color='#ff7f03',label='S2')
   rects2_S2 = ax.bar(x + width/2, obs_S2, width-0.05, bottom=topbottom_O, color='#ff7f03')
   
   zipped_lists_M = zip(topbottom_M,mod_S2)
   zipped_lists_O = zip(topbottom_O,obs_S2)
   topbottom_M=[x + y for (x, y) in zipped_lists_M] 
   topbottom_O=[x + y for (x, y) in zipped_lists_O] 
   rects1_K1 = ax.bar(x - width/2, mod_K1, width-0.05, bottom=topbottom_M, color='#2ca02c', label='K1')
   rects2_K1 = ax.bar(x + width/2, obs_K1, width-0.05, bottom=topbottom_O, color='#2ca02c')
   
   zipped_lists_M = zip(topbottom_M,mod_K1)
   zipped_lists_O = zip(topbottom_O,obs_K1)
   topbottom_M=[x + y for (x, y) in zipped_lists_M]
   topbottom_O=[x + y for (x, y) in zipped_lists_O]
   rects1_O1 = ax.bar(x - width/2, mod_O1, width-0.05, bottom=topbottom_M, color='#d62728', label='O1')
   rects2_O1 = ax.bar(x + width/2, obs_O1, width-0.05, bottom=topbottom_O, color='#d62728')
   #
   zipped_lists_M = zip(topbottom_M,mod_O1)
   zipped_lists_O = zip(topbottom_O,obs_O1)
   topbottom_M=[x + y for (x, y) in zipped_lists_M]
   topbottom_O=[x + y for (x, y) in zipped_lists_O]
   rects1_N2 = ax.bar(x - width/2, mod_N2, width-0.05, bottom=topbottom_M, color='#bcdb22', label='N2')
   rects2_N2 = ax.bar(x + width/2, obs_N2, width-0.05, bottom=topbottom_O, color='#bcdb22')
   
   zipped_lists_M = zip(topbottom_M,mod_N2)
   zipped_lists_O = zip(topbottom_O,obs_N2)
   topbottom_M=[x + y for (x, y) in zipped_lists_M]
   topbottom_O=[x + y for (x, y) in zipped_lists_O]
   rects1_P1 = ax.bar(x - width/2, mod_P1, width-0.05, bottom=topbottom_M, color='#17becf', label='P1')
   rects2_P1 = ax.bar(x + width/2, obs_P1, width-0.05, bottom=topbottom_O, color='#17becf')
   
   zipped_lists_M = zip(topbottom_M,mod_P1)
   zipped_lists_O = zip(topbottom_O,obs_P1)
   topbottom_M=[x + y for (x, y) in zipped_lists_M]
   topbottom_O=[x + y for (x, y) in zipped_lists_O]
   rects1_Q1 = ax.bar(x - width/2, mod_Q1, width-0.05, bottom=topbottom_M, color='#9467bd', label='Q1')
   rects2_Q1 = ax.bar(x + width/2, obs_Q1, width-0.05, bottom=topbottom_O, color='#9467bd')
   
   zipped_lists_M = zip(topbottom_M,mod_Q1)
   zipped_lists_O = zip(topbottom_O,obs_Q1)
   topbottom_M=[x + y for (x, y) in zipped_lists_M]
   topbottom_O=[x + y for (x, y) in zipped_lists_O]
   rects1_K2 = ax.bar(x - width/2, mod_K2, width-0.05, bottom=topbottom_M, color='#e377c2', label='K2')
   rects2_K2 = ax.bar(x + width/2, obs_K2, width-0.05, bottom=topbottom_O, color='#e377c2')
   
   zipped_lists_M = zip(topbottom_M,mod_K2)
   zipped_lists_O = zip(topbottom_O,obs_K2)
   topbottom_M=[x + y for (x, y) in zipped_lists_M]
   topbottom_O=[x + y for (x, y) in zipped_lists_O]
   
   # Add some text for labels, title and custom x-axis tick labels, etc.
   ax.set_ylabel('Tidal Amplitudes [cm]')
   ax.set_title('Tidal Amplitudes')
   ax.set_xticks(x)
   ax.set_xticklabels(labels,fontweight='bold')
   for xtick, color in zip(ax.get_xticklabels(), TOT_color_stz):
       xtick.set_color(color)
       xtick.set_fontsize(fontsize_tg)
   ax.legend()
   
   plt.grid ()
   
   def autolabel1(rects,g_height):
       bin_idx=0
       for rect in rects:
           height = g_height
           ax.annotate('M',
                       xy=(rect.get_x() + rect.get_width() / 2, height[bin_idx]),
                       xytext=(0, 3),  # 3 points vertical offset
                       textcoords="offset points",
                       ha='center', va='bottom')
           bin_idx=bin_idx+1
   
   def autolabel2(rects,g_height):
       bin_idx=0
       for rect in rects:
           height = g_height
           ax.annotate('O',
                       xy=(rect.get_x() + rect.get_width() / 2, height[bin_idx]),
                       xytext=(0, 3),  # 3 points vertical offset
                       textcoords="offset points",
                       ha='center', va='bottom')
           bin_idx=bin_idx+1
   
   autolabel1(rects1_K2,topbottom_M)
   autolabel2(rects2_K2,topbottom_O)
   
   if where_box=='Med':
          plt.savefig(workdir_path+'GLOBAL_A.jpg')
   elif where_box=='AtlBox':
          plt.savefig(workdir_path+'GLOBAL_A_AB.jpg')
   plt.clf()
   
   # AMPLITUDES Global plots ONLY OBS
   x = np.arange(len(labels))  # the label locations
   width = 0.8  # the width of the bars

   fig,ax=plt.subplots( figsize=(80,24))
   plt.rc('font', size=40)

   rects2_M2 = ax.bar(x, obs_M2, width-0.05, color='#1f77b4', label='M2')

   topbottom_O=obs_M2
   rects2_S2 = ax.bar(x, obs_S2, width-0.05, bottom=topbottom_O, color='#ff7f03',label='S2')

   zipped_lists_O = zip(topbottom_O,obs_S2)
   topbottom_O=[x + y for (x, y) in zipped_lists_O]
   rects2_K1 = ax.bar(x, obs_K1, width-0.05, bottom=topbottom_O, color='#2ca02c', label='K1')

   zipped_lists_O = zip(topbottom_O,obs_K1)
   topbottom_O=[x + y for (x, y) in zipped_lists_O]
   rects2_O1 = ax.bar(x, obs_O1, width-0.05, bottom=topbottom_O, color='#d62728', label='O1')
   #
   zipped_lists_O = zip(topbottom_O,obs_O1)
   topbottom_O=[x + y for (x, y) in zipped_lists_O]
   rects2_N2 = ax.bar(x, obs_N2, width-0.05, bottom=topbottom_O, color='#bcdb22', label='N2')

   zipped_lists_O = zip(topbottom_O,obs_N2)
   topbottom_O=[x + y for (x, y) in zipped_lists_O]
   rects2_P1 = ax.bar(x, obs_P1, width-0.05, bottom=topbottom_O, color='#17becf', label='P1')

   zipped_lists_O = zip(topbottom_O,obs_P1)
   topbottom_O=[x + y for (x, y) in zipped_lists_O]
   rects2_Q1 = ax.bar(x, obs_Q1, width-0.05, bottom=topbottom_O, color='#9467bd', label='Q1')

   zipped_lists_O = zip(topbottom_O,obs_Q1)
   topbottom_O=[x + y for (x, y) in zipped_lists_O]
   rects2_K2 = ax.bar(x, obs_K2, width-0.05, bottom=topbottom_O, color='#e377c2', label='K2')

   zipped_lists_O = zip(topbottom_O,obs_K2)
   topbottom_O=[x + y for (x, y) in zipped_lists_O]

   # Add some text for labels, title and custom x-axis tick labels, etc.
   ax.set_ylabel('Tidal Amplitudes [cm]')
   ax.set_title('Tidal Amplitudes')
   ax.set_xticks(x)
   ax.set_xticklabels(labels,fontweight='bold')
   for xtick, color in zip(ax.get_xticklabels(), TOT_color_stz):
       xtick.set_color(color)
       xtick.set_fontsize(fontsize_tg)
   ax.legend()

   plt.grid ()

   if where_box=='Med':
          plt.savefig(workdir_path+'GLOBAL_A_onlyobs.jpg')
   elif where_box=='AtlBox':
          plt.savefig(workdir_path+'GLOBAL_A_AB_onlyobs.jpg')
   plt.clf()


   # AMPLITUDES Global plots ONLY MOD
   x = np.arange(len(labels))  # the label locations
   width = 0.8  # the width of the bars

   fig,ax=plt.subplots( figsize=(80,24))
   plt.rc('font', size=40)
   #
   rects1_M2 = ax.bar(x, mod_M2, width-0.05, color='#1f77b4', label='M2')

   topbottom_M=mod_M2
   rects1_S2 = ax.bar(x, mod_S2, width-0.05, bottom=topbottom_M, color='#ff7f03',label='S2')

   zipped_lists_M = zip(topbottom_M,mod_S2)
   topbottom_M=[x + y for (x, y) in zipped_lists_M]
   rects1_K1 = ax.bar(x, mod_K1, width-0.05, bottom=topbottom_M, color='#2ca02c', label='K1')

   zipped_lists_M = zip(topbottom_M,mod_K1)
   topbottom_M=[x + y for (x, y) in zipped_lists_M]
   rects1_O1 = ax.bar(x, mod_O1, width-0.05, bottom=topbottom_M, color='#d62728', label='O1')
   #
   zipped_lists_M = zip(topbottom_M,mod_O1)
   topbottom_M=[x + y for (x, y) in zipped_lists_M]
   rects1_N2 = ax.bar(x, mod_N2, width-0.05, bottom=topbottom_M, color='#bcdb22', label='N2')

   zipped_lists_M = zip(topbottom_M,mod_N2)
   topbottom_M=[x + y for (x, y) in zipped_lists_M]
   rects1_P1 = ax.bar(x, mod_P1, width-0.05, bottom=topbottom_M, color='#17becf', label='P1')

   zipped_lists_M = zip(topbottom_M,mod_P1)
   topbottom_M=[x + y for (x, y) in zipped_lists_M]
   rects1_Q1 = ax.bar(x, mod_Q1, width-0.05, bottom=topbottom_M, color='#9467bd', label='Q1')

   zipped_lists_M = zip(topbottom_M,mod_Q1)
   topbottom_M=[x + y for (x, y) in zipped_lists_M]
   rects1_K2 = ax.bar(x, mod_K2, width-0.05, bottom=topbottom_M, color='#e377c2', label='K2')

   zipped_lists_M = zip(topbottom_M,mod_K2)
   topbottom_M=[x + y for (x, y) in zipped_lists_M]

   # Add some text for labels, title and custom x-axis tick labels, etc.
   ax.set_ylabel('Tidal Amplitudes [cm]')
   ax.set_title('Tidal Amplitudes')
   ax.set_xticks(x)
   ax.set_xticklabels(labels,fontweight='bold')
   for xtick, color in zip(ax.get_xticklabels(), TOT_color_stz):
       xtick.set_color(color)
       xtick.set_fontsize(fontsize_tg)
   ax.legend()

   plt.grid ()

   if where_box=='Med':
          plt.savefig(workdir_path+'GLOBAL_A_onlymod.jpg')
   elif where_box=='AtlBox':
          plt.savefig(workdir_path+'GLOBAL_A_AB_onlymod.jpg')
   plt.clf()
 
   # PHASES Global plots
   
   plt.figure(figsize=(80,24)) 
   plt.rc('font', size=40) 
   
   labels=[ GLOB_P_mod[i][0] for i in range(1,N_stz+1) ]
   mod_M2 = [ GLOB_P_mod[j][1] for j in range (1,N_stz+1)]
   obs_M2 = [ GLOB_P_obs[j][1] for j in range (1,N_stz+1)]
   mod_S2 = [ GLOB_P_mod[j][2] for j in range (1,N_stz+1)]
   obs_S2 = [ GLOB_P_obs[j][2] for j in range (1,N_stz+1)]
   mod_K1 = [ GLOB_P_mod[j][3] for j in range (1,N_stz+1)]
   obs_K1 = [ GLOB_P_obs[j][3] for j in range (1,N_stz+1)]
   mod_O1 = [ GLOB_P_mod[j][4] for j in range (1,N_stz+1)]
   obs_O1 = [ GLOB_P_obs[j][4] for j in range (1,N_stz+1)]
   #
   mod_N2 = [ GLOB_P_mod[j][5] for j in range (1,N_stz+1)]
   obs_N2 = [ GLOB_P_obs[j][5] for j in range (1,N_stz+1)]
   mod_P1 = [ GLOB_P_mod[j][6] for j in range (1,N_stz+1)]
   obs_P1 = [ GLOB_P_obs[j][6] for j in range (1,N_stz+1)]
   mod_Q1 = [ GLOB_P_mod[j][7] for j in range (1,N_stz+1)]
   obs_Q1 = [ GLOB_P_obs[j][7] for j in range (1,N_stz+1)]
   mod_K2 = [ GLOB_P_mod[j][8] for j in range (1,N_stz+1)]
   obs_K2 = [ GLOB_P_obs[j][8] for j in range (1,N_stz+1)]
   
   x = np.arange(len(labels))  # the label locations
   width = 0.45  # the width of the bars
   
   fig,ax=plt.subplots( figsize=(80,24)) 
   #
   rects1_M2 = ax.bar(x - width/2, mod_M2, width-0.05, color='#1f77b4', label='M2')
   rects2_M2 = ax.bar(x + width/2, obs_M2, width-0.05, color='#1f77b4')
   
   topbottom_M=mod_M2
   topbottom_O=obs_M2
   rects1_S2 = ax.bar(x - width/2, mod_S2, width-0.05, bottom=topbottom_M, color='#ff7f03',label='S2')
   rects2_S2 = ax.bar(x + width/2, obs_S2, width-0.05, bottom=topbottom_O, color='#ff7f03')
   
   zipped_lists_M = zip(topbottom_M,mod_S2)
   zipped_lists_O = zip(topbottom_O,obs_S2)
   topbottom_M=[x + y for (x, y) in zipped_lists_M]
   topbottom_O=[x + y for (x, y) in zipped_lists_O]
   rects1_K1 = ax.bar(x - width/2, mod_K1, width-0.05, bottom=topbottom_M, color='#2ca02c', label='K1')
   rects2_K1 = ax.bar(x + width/2, obs_K1, width-0.05, bottom=topbottom_O, color='#2ca02c')
   
   zipped_lists_M = zip(topbottom_M,mod_K1)
   zipped_lists_O = zip(topbottom_O,obs_K1)
   topbottom_M=[x + y for (x, y) in zipped_lists_M]
   topbottom_O=[x + y for (x, y) in zipped_lists_O]
   rects1_O1 = ax.bar(x - width/2, mod_O1, width-0.05, bottom=topbottom_M, color='#d62728', label='O1')
   rects2_O1 = ax.bar(x + width/2, obs_O1, width-0.05, bottom=topbottom_O, color='#d62728')
   #
   zipped_lists_M = zip(topbottom_M,mod_O1)
   zipped_lists_O = zip(topbottom_O,obs_O1)
   topbottom_M=[x + y for (x, y) in zipped_lists_M]
   topbottom_O=[x + y for (x, y) in zipped_lists_O]
   rects1_N2 = ax.bar(x - width/2, mod_N2, width-0.05, bottom=topbottom_M, color='#bcdb22', label='N2')
   rects2_N2 = ax.bar(x + width/2, obs_N2, width-0.05, bottom=topbottom_O, color='#bcdb22')
   
   zipped_lists_M = zip(topbottom_M,mod_N2)
   zipped_lists_O = zip(topbottom_O,obs_N2)
   topbottom_M=[x + y for (x, y) in zipped_lists_M]
   topbottom_O=[x + y for (x, y) in zipped_lists_O]
   rects1_P1 = ax.bar(x - width/2, mod_P1, width-0.05, bottom=topbottom_M, color='#17becf', label='P1')
   rects2_P1 = ax.bar(x + width/2, obs_P1, width-0.05, bottom=topbottom_O, color='#17becf')
   
   zipped_lists_M = zip(topbottom_M,mod_P1)
   zipped_lists_O = zip(topbottom_O,obs_P1)
   topbottom_M=[x + y for (x, y) in zipped_lists_M]
   topbottom_O=[x + y for (x, y) in zipped_lists_O]
   rects1_Q1 = ax.bar(x - width/2, mod_Q1, width-0.05, bottom=topbottom_M, color='#9467bd', label='Q1')
   rects2_Q1 = ax.bar(x + width/2, obs_Q1, width-0.05, bottom=topbottom_O, color='#9467bd')
   
   zipped_lists_M = zip(topbottom_M,mod_Q1)
   zipped_lists_O = zip(topbottom_O,obs_Q1)
   topbottom_M=[x + y for (x, y) in zipped_lists_M]
   topbottom_O=[x + y for (x, y) in zipped_lists_O]
   rects1_K2 = ax.bar(x - width/2, mod_K2, width-0.05, bottom=topbottom_M, color='#e377c2', label='K2')
   rects2_K2 = ax.bar(x + width/2, obs_K2, width-0.05, bottom=topbottom_O, color='#e377c2')
   
   zipped_lists_M = zip(topbottom_M,mod_K2)
   zipped_lists_O = zip(topbottom_O,obs_K2)
   topbottom_M=[x + y for (x, y) in zipped_lists_M]
   topbottom_O=[x + y for (x, y) in zipped_lists_O]
   
   # Add some text for labels, title and custom x-axis tick labels, etc.
   ax.set_ylabel('Tidal Phases [deg]')
   ax.set_title('Tidal Phases')
   ax.set_xticks(x)
   ax.set_xticklabels(labels,fontweight='bold')
   for xtick, color in zip(ax.get_xticklabels(), TOT_color_stz):
       xtick.set_color(color)
       xtick.set_fontsize(fontsize_tg)
   ax.legend()
   
   plt.grid ()
   
   def autolabel1(rects,g_height):
       bin_idx=0
       for rect in rects:
           height = g_height
           ax.annotate('M',
                       xy=(rect.get_x() + rect.get_width() / 2, height[bin_idx]),
                       xytext=(0, 3),  # 3 points vertical offset
                       textcoords="offset points",
                       ha='center', va='bottom')
           bin_idx=bin_idx+1
   
   def autolabel2(rects,g_height):
       bin_idx=0
       for rect in rects:
           height = g_height
           ax.annotate('O',
                       xy=(rect.get_x() + rect.get_width() / 2, height[bin_idx]),
                       xytext=(0, 3),  # 3 points vertical offset
                       textcoords="offset points",
                       ha='center', va='bottom')
           bin_idx=bin_idx+1
   
   autolabel1(rects1_K2,topbottom_M)
   autolabel2(rects2_K2,topbottom_O)
   if where_box=='Med':
          plt.savefig(workdir_path+'GLOBAL_P.jpg')
   elif where_box=='AtlBox':
          plt.savefig(workdir_path+'GLOBAL_P_AB.jpg')
   plt.clf()
   
   
   
   # GLOBAL VECTORIAL DISTANCES PLOT
   
   plt.figure(figsize=(80,24)) 
   plt.rc('font', size=40)  
   
   labels=[ d_foreman[i][0] for i in range(1,N_stz+1) ]
   mod_M2 = [ d_foreman[j][1] for j in range (1,N_stz+1)]
   mod_S2 = [ d_foreman[j][2] for j in range (1,N_stz+1)]
   mod_K1 = [ d_foreman[j][3] for j in range (1,N_stz+1)]
   mod_O1 = [ d_foreman[j][4] for j in range (1,N_stz+1)]
   #
   mod_N2 = [ d_foreman[j][5] for j in range (1,N_stz+1)]
   mod_P1 = [ d_foreman[j][6] for j in range (1,N_stz+1)]
   mod_Q1 = [ d_foreman[j][7] for j in range (1,N_stz+1)]
   mod_K2 = [ d_foreman[j][8] for j in range (1,N_stz+1)]
   
   # Compute and print statistics on vectorial distances
   if flag_15stats == 1:
      # Table for Comp with literature on 15 TGs
      Lit_file = open(workdir_path+"stats_lit.txt","w")
      print ('# TG dataset: ',labels,file=Lit_file)
      print ('# Mean vec dist [cm] Max vec dist [cm] (where max vec dist)',file=Lit_file)
      print ('M2 EAS: ', round(np.mean(mod_M2),2), round(np.max(mod_M2),2),'(',labels[np.argmax(mod_M2)],')',file=Lit_file)
      print ('M2 PALMA: ', round(np.mean(PALMA_d_M2),2), round(np.max(PALMA_d_M2),2),'(',labels[np.argmax(PALMA_d_M2)],')',file=Lit_file)
      print ('M2 TSIMPLIS: ', round(np.mean(TSIMPLIS_d_M2),2), round(np.max(TSIMPLIS_d_M2),2),'(',labels[np.argmax(TSIMPLIS_d_M2)],')',file=Lit_file)
      print ('---',file=Lit_file)
      print ('S2 EAS: ', round(np.mean(mod_S2),2), round(np.max(mod_S2),2),'(',labels[np.argmax(mod_S2)],')',file=Lit_file)
      print ('S2 PALMA: ', round(np.mean(PALMA_d_S2),2), round(np.max(PALMA_d_S2),2),'(',labels[np.argmax(PALMA_d_S2)],')',file=Lit_file)
      print ('S2 TSIMPLIS: ', round(np.mean(TSIMPLIS_d_S2),2), round(np.max(TSIMPLIS_d_S2),2),'(',labels[np.argmax(TSIMPLIS_d_S2)],')',file=Lit_file)
      print ('---',file=Lit_file)
      print ('K1 EAS: ', round(np.mean(mod_K1),2), round(np.max(mod_K1),2),'(',labels[np.argmax(mod_K1)],')',file=Lit_file)
      print ('K1 PALMA: ', round(np.mean(PALMA_d_K1),2), round(np.max(PALMA_d_K1),2),'(',labels[np.argmax(PALMA_d_K1)],')',file=Lit_file)
      print ('K1 TSIPLIS: ', round(np.mean(TSIMPLIS_d_K1),2), round(np.max(TSIMPLIS_d_K1),2),'(',labels[np.argmax(TSIMPLIS_d_K1)],')',file=Lit_file)
      print ('---',file=Lit_file)
      print ('O1 EAS: ', round(np.mean(mod_O1),2), round(np.max(mod_O1),2),'(',labels[np.argmax(mod_O1)],')',file=Lit_file)
      print ('O1 PALMA: ', round(np.mean(PALMA_d_O1),2), round(np.max(PALMA_d_O1),2),'(',labels[np.argmax(PALMA_d_O1)],')',file=Lit_file)
      print ('K1 TSIPLIS: ', round(np.mean(TSIMPLIS_d_O1),2), round(np.max(TSIMPLIS_d_O1),2),'(',labels[np.argmax(TSIMPLIS_d_O1)],')',file=Lit_file)
      #try:
      #   print ('O1 TSIMPLIS: ', round(np.mean(i for i in TSIMPLIS_d_O1 if i != ''),2), round(np.max(i for i in TSIMPLIS_d_O1 if i != ''),2),'(',labels[np.nanargmax(i for i in TSIMPLIS_d_O1 if i != '')],')',file=Lit_file)
      #except:
      #   print ('ERROR reading O1 TSIMPLIS values..')
      print ('---',file=Lit_file)
      Lit_file.close()

   # Plot vectorial distances
   x = np.arange(len(labels))  # the label locations
   width = 0.45*2  # the width of the bars
   
   fig,ax=plt.subplots( figsize=(80,24)) 
   #
   rects1_M2 = ax.bar(x, mod_M2, width-0.05, color='#1f77b4', label='M2')
   
   topbottom_M=mod_M2
   rects1_S2 = ax.bar(x, mod_S2, width-0.05, bottom=topbottom_M, color='#ff7f03',label='S2')
   
   zipped_lists_M = zip(topbottom_M,mod_S2)
   topbottom_M=[x + y for (x, y) in zipped_lists_M]
   rects1_K1 = ax.bar(x, mod_K1, width-0.05, bottom=topbottom_M, color='#2ca02c', label='K1')
   
   zipped_lists_M = zip(topbottom_M,mod_K1)
   topbottom_M=[x + y for (x, y) in zipped_lists_M]
   rects1_O1 = ax.bar(x, mod_O1, width-0.05, bottom=topbottom_M, color='#d62728', label='O1')
   #
   zipped_lists_M = zip(topbottom_M,mod_O1)
   topbottom_M=[x + y for (x, y) in zipped_lists_M]
   rects1_N2 = ax.bar(x, mod_N2, width-0.05, bottom=topbottom_M, color='#bcdb22', label='N2')
   
   zipped_lists_M = zip(topbottom_M,mod_N2)
   topbottom_M=[x + y for (x, y) in zipped_lists_M]
   rects1_P1 = ax.bar(x, mod_P1, width-0.05, bottom=topbottom_M, color='#17becf', label='P1')
   
   zipped_lists_M = zip(topbottom_M,mod_P1)
   topbottom_M=[x + y for (x, y) in zipped_lists_M]
   rects1_Q1 = ax.bar(x, mod_Q1, width-0.05, bottom=topbottom_M, color='#9467bd', label='Q1')
   
   zipped_lists_M = zip(topbottom_M,mod_Q1)
   topbottom_M=[x + y for (x, y) in zipped_lists_M]
   rects1_K2 = ax.bar(x, mod_K2, width-0.05, bottom=topbottom_M, color='#e377c2', label='K2')
   
   zipped_lists_M = zip(topbottom_M,mod_K2)
   topbottom_M=[x + y for (x, y) in zipped_lists_M]
   
   # Add some text for labels, title and custom x-axis tick labels, etc.
   ax.set_ylabel('Vectorial Distances [cm]')
   ax.set_title('Vectorial Distances per TG per tidal component') 
   ax.set_xticks(x)
   ax.set_xticklabels(labels,fontweight='bold')
   for xtick, color in zip(ax.get_xticklabels(), TOT_color_stz):
       xtick.set_color(color)
       xtick.set_fontsize(fontsize_tg) 
   ax.legend()
   
   plt.grid ()
   
   if where_box=='Med':
          plt.savefig(workdir_path+'GLOBAL_d_Foreman.jpg')
   elif where_box=='AtlBox':
          plt.savefig(workdir_path+'GLOBAL_d_Foreman_AB.jpg')
   plt.clf()
   
   
   # GLOBAL RMS DISTANCES PLOT
   
   plt.figure(figsize=(20,10))
   plt.rc('font', size=20)
   
   labels=[ d_foreman[0][j] for j in range(1,N_comp+1) ]
   
   x = np.arange(len(labels))  # the label locations
   width = 0.45*2  # the width of the bars
   
   fig,ax=plt.subplots( figsize=(20,10))
   #
   comp_color=['#1f77b4','#ff7f03','#2ca02c','#d62728','#bcdb22','#17becf','#9467bd','#e377c2']
   rects1 = ax.bar(x, RMSm, width-0.05,color=comp_color)
   
   
   # Add some text for labels, title and custom x-axis tick labels, etc.
   ax.set_ylabel('RMSm [cm]')
   ax.set_title('Root Mean Square misfits - '+where_box)
   ax.set_xticks(x)
   ax.set_xticklabels(labels,fontweight='bold')
   
   plt.grid ()
   
   
   def autolabel(rects):
       bin_idx=0
       for rect in rects:
           height = rect.get_height()
           ax.annotate(round(RMSm[bin_idx],1),
                       xy=(rect.get_x() + rect.get_width() / 2, height ),
                       xytext=(0, 3),  # 3 points vertical offset
                       textcoords="offset points",
                       ha='center', va='bottom')
           bin_idx=bin_idx+1
   
   autolabel(rects1)
   if where_box=='Med':
          plt.savefig(workdir_path+'GLOBAL_RMSm.jpg')
   elif where_box=='AtlBox':
          plt.savefig(workdir_path+'GLOBAL_RMSm_AB.jpg')
   plt.clf()

print ('Output path: ',workdir_path)
