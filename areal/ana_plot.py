#
# imports
from __future__ import print_function, division
import os 
import time 
import sys
import warnings
#
warnings.filterwarnings("ignore") # Avoid warnings
#
import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as NC
from scipy.optimize import curve_fit
from scipy import stats
import collections
import pandas as pd
import csv
import math
import datetime
from operator import itemgetter 
#
from numpy import ma
from scipy.stats import ks_2samp 
from statsmodels.distributions.empirical_distribution import ECDF
from mpl_toolkits.basemap import Basemap 
from matplotlib.colors import LogNorm 
from pyresample.geometry import SwathDefinition
from pyresample.kd_tree import resample_nearest
from scipy import interpolate
#
#
# by AC Goglio (CMCC)
# annachiara.goglio@cmcc.it
#
# Written: 01/11/2019
# Modified: 12/03/2021
#
# Script to analize and plot the harmonic analisi reults on area
#
#################################################################
# The user should modify the following lines to set his run
#################################################################
# General run parameters:
#---------------------
# work dir path and bathymetry path/name
workdir_path = '/work/oda/ag15419/tmp/atl_amph/'
model_bathy='/work/oda/ag15419/PHYSW24_DATA/TIDES/DATA0/bathy_meter.nc'
model_meshmask='/work/oda/ag15419/PHYSW24_DATA/TIDES/DATA0/mesh_mask.nc'
#
# Dates
# Choose start and end dates of the period (format dd/mm/yyyy)
inidate = '01/07/2017'
enddate = '31/12/2017'

# TPXO9 path
tpxo9_path='/work/oda/ag15419/OTPS/OTPSnc/DATA/NC/'

# FLAGS for different analysis (set 1 to activate the tasks)

# For Amplitude/Phase maps (1 map per tidal component)
amppha_flag=0
ampha_var='AtlBox' # NEEDED if ampha_flag=1; and if tpxo2eas_flag=1 or ampha_tpxo=1 (for this cases only AtlBox Vs other is implemented)
# Values can be 
# ) AmpPha to plot amplitude/phase maps
# ) AtlBox to plot amplitude/phase maps with AtlOcean proper palette
# ) Amp to plot only amplitude maps
# ) Pha to plot only phase maps
# ) Pha_Pal to compare phase maps wrt Palma et al
# ) Amp_Pal to compare amplitude maps wrt Palma et al 
# ) Amp_Ar to compare amplitude maps wrt Arabelos et al
# ) Pha_Ar to compare phase maps wrt Arabelos et al 
# ) AmpPha_Ag to compare amplitude/phase maps wrt Agresti 

ampha_tpxo=0
# For TPXO9 Amplitude/Phase maps on TPXO grid (1 map per tidal component)

pha_tpxo=1
# For TPXO9 Phase maps on TPXO grid (1 map per tidal component)

ampha_tpxo_atlamph=0
# For TPXO9 Amplitude/Phase maps on TPXO grid in North Atlantic (ONLY M2 component, to be extended..)

doseong_flag=0 # TO compute Do-Seong factor for EAS system
# For the following maps (1 map per factor):
# Tidal Form Factor [Do-Seong] F=(A_K1+A_O1)/(A_M2+A_S2)
# Tidal Envelope Factor [Do-Seong] E=(A_M2+A_N2)/(A_M2+A_S2)
# Tidal envelope asymmetric factor [Do-Seong] Ea=cos(P_K1+P_O1-P_M2)
# WARNING: the following fields are needed: A_K1 A_O1 A_M2 A_S2 A_N2 P_K1 P_O1 P_M2
doseong_tpxo=0
# For computing DoSeong factor from TPXO9 model on TPXO grid

tpxo2eas_flag=0
# To interpolate tpxo9 1/30 to MED24 grid and plot Amplitude/Phase maps

diff_tpxoeas_flag=0
# For amplitude diffs between eas and tpxo on MED24 grid

vectorial_dist_flag=0
# For vectorial distances between eas and tpxo9 on MED24 grid

bathy_diff_flag=0
# For diffs between bathymethries eas vs tpxo on MED24 grid

########################################################
# DO NOT CHANGE THE CODE BELOW THIS LINES
########################################################

# Parameter setting
#--------------------
# MODEL DATASETS
model_path=workdir_path
model_fileprename='amppha' # DO NOT change this
model_postname='mod_Tides8_v31' # WARNING: Use the same string as in fit_marea.py

bathylim4RMSE=0 # Bathymetry threshold for RMSE (only grid points with bathy>bathylim4RMSE are taken into account)

dates_label=inidate[6:10]+inidate[3:5]+inidate[0:2]+'_'+enddate[6:10]+enddate[3:5]+enddate[0:2]

# TPXO9 model file names and label
tpxo9_fileprename='h_'
tpxo9_postname='_tpxo9_atlas_30_v2'
tpxo9_label='TPXO9'

print ('Whole time interval: ',dates_label)
print ('Model file templates: ',model_path,model_fileprename,'?D_','%periods_of_analysis%','_%vlev%_','%field%',dates_label,'_',model_postname,'.nc')
if ampha_tpxo == 1 or doseong_tpxo == 1 or diff_tpxoeas_flag == 1 or vectorial_dist_flag == 1 or bathy_diff_flag == 1:
  print ('TPXO9 Model file templates: ',tpxo9_path,tpxo9_fileprename,'%tidal_component%',tpxo9_postname,'.nc')

# field(s) to be plotted
#-----------------------
# 2D fields:
field_2d_name=['M2_Amp','K1_Amp','O1_Amp','S2_Amp','P1_Amp','N2_Amp','Q1_Amp','K2_Amp']
P_field_2d_name=['M2_Pha','K1_Pha','O1_Pha','S2_Pha','P1_Pha','N2_Pha','Q1_Pha','K2_Pha']
field_2d_units=['cm','cm','cm','cm','cm','cm','cm','cm']
P_field_2d_units=['deg','deg','deg','deg','deg','deg','deg','deg']

# Lat/Lon from BATHYMETRY (WARNING: the name must be the one in the bathy_meter.nc, NOT in NEMO outputs)
latitude_name='nav_lat'
longitude_name='nav_lon'

# Move to the work dir 
if os.path.exists(workdir_path) : 
 os.chdir(workdir_path) # mv to work dir
 print('I am moving in the directory: ',os.getcwd()) # pwd
 os.listdir() # ls
else:
 os.mkdir(workdir_path)
 os.chdir(workdir_path)
 print('I am in a new work directory: ',os.getcwd()) 


# List input/output cases  
print ('Fields 2D: ',field_2d_name)

##################################
# 2D VARIABLES :
##################################

# ===============================
# Loop on 2D VARS :
# ===============================
# Amplitude/Phase maps (1 map per tidal component)
if amppha_flag == 1:
   for idx_2d in range(0,len(field_2d_name)):
            var_2d=field_2d_name[idx_2d]
            var_2d_udm=field_2d_units[idx_2d]

            P_var_2d=P_field_2d_name[idx_2d]
            P_var_2d_udm=P_field_2d_units[idx_2d]

            # Build the path/name of the nc file and open it 
            nc2open=model_path+model_fileprename+'2D_'+'0_'+'sossheig'+'_'+dates_label+'_'+model_postname+'.nc'
            print ('Input file = ',nc2open)
            model = NC.Dataset(nc2open,'r')

            # Read lat, lon and fild values 
            nc2open3=model_bathy # tidal bathimetry
            model3 = NC.Dataset(nc2open3,'r')
            nc2open4=model_meshmask # mesh mask
            model4 = NC.Dataset(nc2open4,'r')

            vals_bathy=model3.variables['Bathymetry'][:]
            lons = model3.variables[longitude_name][:]
            lats = model3.variables[latitude_name][:]
            vals = model.variables[var_2d][:]*100 # Want cm not meters!
            P_vals = model.variables[P_var_2d][:]
         
            vals_land=model4.variables['tmask'][0,0,:,:]
            vals_land=np.squeeze(vals_land)

            #thresh = 0.0000
            #mask = np.abs(vals) == thresh
            #vals_ma = np.ma.masked_where(mask, vals)
            vals_ma = np.ma.masked_where(vals_land, vals)

            # Plot the map and save in the path/name

            if ampha_var == 'Pha_Pal':
               plotname=workdir_path+model_fileprename+'2D_'+'0_'+var_2d+'_'+dates_label+'_'+model_postname+'_Pha_VsPalma.jpg' 
            elif ampha_var == 'Amp_Pal':
               plotname=workdir_path+model_fileprename+'2D_'+'0_'+var_2d+'_'+dates_label+'_'+model_postname+'_Amp_VsPalma.jpg'  
            elif ampha_var == 'AmpPha_Ag':
               plotname=workdir_path+model_fileprename+'2D_'+'0_'+var_2d+'_'+dates_label+'_'+model_postname+'_AmpPha_VsAgresti.jpg'
            elif ampha_var == 'AmpPha':
               plotname=workdir_path+model_fileprename+'2D_'+'0_'+var_2d+'_'+dates_label+'_'+model_postname+'_AmpPha.jpg'
            elif ampha_var == 'Amp':
               plotname=workdir_path+model_fileprename+'2D_'+'0_'+var_2d+'_'+dates_label+'_'+model_postname+'_Amp.jpg'
            elif ampha_var == 'Pha':
               plotname=workdir_path+model_fileprename+'2D_'+'0_'+var_2d+'_'+dates_label+'_'+model_postname+'_Pha.jpg'
            elif ampha_var == 'Amp_Ar':
               plotname=workdir_path+model_fileprename+'2D_'+'0_'+var_2d+'_'+dates_label+'_'+model_postname+'_Amp_VsArabelos.jpg'
            elif ampha_var == 'Pha_Ar':
               plotname=workdir_path+model_fileprename+'2D_'+'0_'+var_2d+'_'+dates_label+'_'+model_postname+'_Pha_VsArabelos.jpg'
            elif ampha_var == 'AtlBox':
               plotname=workdir_path+model_fileprename+'2D_'+'0_'+var_2d+'_'+dates_label+'_'+model_postname+'_AmpPha_AtlBox.jpg'

            print ('Plot path/name: ',plotname)

            plt.figure(figsize=(20,10))
            plt.rc('font', size=16)
            # Plot Title
            if ampha_var == 'Amp' or ampha_var == 'Amp_Ar' or ampha_var == 'Amp_Pal':
               plt.title (var_2d[:2]+' Amplitude - Harmonic Analysis '+dates_label)
            elif ampha_var == 'Pha' or ampha_var == 'Pha_Ar' or ampha_var == 'Pha_Pal':
               plt.title (var_2d[:2]+' Phase - Harmonic Analysis '+dates_label)
            elif ampha_var == 'AmpPha'or ampha_var == 'AmpPha_Ag' or ampha_var == 'AtlBox':
               plt.title (var_2d[:2]+' Amplitude and Phase - Harmonic Analysis '+dates_label)

            # Read the coordinates for the plot 
            lon_0 = lons.mean()
            llcrnrlon = lons.min()
            urcrnrlon = lons.max()
            lat_0 = lats.mean()
            llcrnrlat = lats.min()
            urcrnrlat = lats.max()

            # Create the map
            m = Basemap(llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat,urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,resolution='c',projection='merc',lat_0=lat_0,lon_0=lon_0)
            xi, yi = m(lons, lats)

            # Plot the frame to the map
            plt.rcParams["axes.linewidth"]  = 1.25

            # Plot the Amplitude or Amplitude/Phase map 
            if ampha_var == 'AmpPha' or ampha_var == 'Amp':
               if var_2d == 'M2_Amp':
                   med_range=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25]
               elif var_2d == 'K1_Amp' or var_2d == 'S2_Amp':
                   med_range=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
               elif var_2d == 'O1_Amp' or var_2d == 'P1_Amp' or var_2d == 'N2_Amp':
                   med_range=[0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5]
               elif var_2d == 'K2_Amp':
                   med_range=[0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3]
               elif var_2d == 'Q1_Amp':
                   med_range=[0,0.25,0.5,0.75,1,1.25,1.5,1.75,2]
               # Amp plot:
               cs = plt.contourf(xi,yi,np.abs(np.squeeze(vals)),med_range,cmap='jet',extend='max')

            # Pha:
            elif ampha_var == 'Pha':
                  P_vals=np.where(P_vals>=360.0,P_vals-360,P_vals)
                  P_vals=np.where(P_vals>=360.0,P_vals-360,P_vals)
                  P_vals=np.where(P_vals<0.0,P_vals+360,P_vals)
                  P_vals=np.where(P_vals<0.0,P_vals+360,P_vals)
                  cs = plt.contourf(xi,yi,np.abs(np.squeeze(P_vals)),[0,15,30,45,60,75,90,105,120,135,150,165,180,195,210,225,240,255,270,285,300,315,330,345,360],cmap='hsv')

            # Amplitude and Phase Comparison wrt [V.Agresti; 2018]
            elif ampha_var == 'AmpPha_Ag':
               if var_2d == 'M2_Amp':
                   med_range=[0,5,10,15,20,25,30,35,40,45,50]
               elif var_2d == 'K1_Amp' or var_2d == 'S2_Amp':
                   med_range=[0,2.5,5,7.5,10,12.5,15,17.5,20,22.5,25]
               elif var_2d == 'O1_Amp' or var_2d == 'P1_Amp' or var_2d == 'N2_Amp':
                   med_range=[0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5]
               elif var_2d == 'K2_Amp':
                   med_range=[0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3]
               elif var_2d == 'Q1_Amp':
                   med_range=[0,0.25,0.5,0.75,1,1.25,1.5,1.75,2]
               # Amp plot:
               cs = plt.contourf(xi,yi,np.abs(np.squeeze(vals)),med_range,cmap='jet',extend='max')

            # Amplitude Comparison with [Palma et al; 2020]
            elif ampha_var == 'Amp_Pal':
               if var_2d == 'M2_Amp':
                  med_range=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21]
               elif var_2d == 'K1_Amp':
                  med_range=[0,1,2,3,4,5,6,7,8,9,10,11,12,13]
               cs = plt.contourf(xi,yi,np.abs(np.squeeze(vals)),med_range,cmap='jet',extend='max')

            # Amplitude Comparison with [Arabelos et al; 2011]
            elif ampha_var == 'Amp_Ar':
               if var_2d == 'M2_Amp':
                  med_range=[0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58]
               elif var_2d == 'K1_Amp':
                  med_range=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18]
               elif var_2d == 'S2_Amp':
                  med_range=[0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40]
               elif var_2d == 'O1_Amp':
                  med_range=[0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7]
               cs = plt.contourf(xi,yi,np.abs(np.squeeze(vals)),med_range,cmap='jet',extend='max')

            # Pha Comparison with [Palma et al; 2020]:
            elif ampha_var == 'Pha_Pal':
                  P_vals=np.where(P_vals>360.0,P_vals-360,P_vals)
                  P_vals=np.where(P_vals>360.0,P_vals-360,P_vals)
                  P_vals=np.where(P_vals<0.0,P_vals+360,P_vals)
                  P_vals=np.where(P_vals<0.0,P_vals+360,P_vals)
                  cs = plt.contourf(xi,yi,np.abs(np.squeeze(P_vals)),[0,30,60,90,120,150,180,210,240,270,300,330,360,390],cmap='jet') 
            
            # Pha Comparison with [Arabelos et al; 2011]:
            elif ampha_var == 'Pha_Ar':
                  P_vals=np.where(P_vals>360.0,P_vals-360,P_vals)
                  P_vals=np.where(P_vals>360.0,P_vals-360,P_vals)
                  P_vals=np.where(P_vals<0.0,P_vals+360,P_vals)
                  P_vals=np.where(P_vals<0.0,P_vals+360,P_vals)
                  cs = plt.contourf(xi,yi,np.abs(np.squeeze(P_vals)),[0,30,60,90,120,150,180,210,240,270,300,330,360],cmap='jet')
            
            # Add the Phase to amplitude maps when required     
            if ampha_var == 'AmpPha' or ampha_var == 'AmpPha_Ag' or ampha_var == 'AtlBox':
               amphid_p = plt.contour(xi,yi,np.squeeze(P_vals),[0,20,40,60,80,100,120,140,160,180,200,220,240,260,280,300,320,340,360],colors='white')
               amphid_n = plt.contour(xi,yi,np.squeeze(P_vals),[-360,-340,-320,-300,-280,-260,-240,-220,-200,-180,-160,-140,-120,-100,-80,-60,-40,-20],colors='white',linestyles='dashed')

            # Adjust the palettes for AtlBox case
            if ampha_var == 'AtlBox':
               if var_2d == 'M2_Amp':
                   #med_range=[0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150]
                   med_range=[40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150]
               elif var_2d == 'S2_Amp':
                   #med_range=[0,2.5,5,7.5,10,12.5,15,17.5,20,22.5,25,27.5,30,32.5,35,37.5,40,42.5,45,47.5,50]
                   med_range=[15,17.5,20,22.5,25,27.5,30,32.5,35,37.5,40,42.5,45,47.5,50]
               elif var_2d == 'N2_Amp':
                   med_range=[0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34]
                   med_range=[8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34]
               elif var_2d == 'K2_Amp':
                   #med_range=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
                   med_range=[4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10,10.5,11,11.5,12,12.5,13,13.5,14,14.5,15,15.5,16]

               elif var_2d == 'K1_Amp' or var_2d == 'O1_Amp':
                   #med_range=[0,1,2,3,4,5,6,7,8,9,10]
                   med_range=[3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10]
               elif var_2d == 'P1_Amp' or var_2d == 'Q1_Amp':
                   #med_range=[0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5]
                   med_range=[0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.15,2.5,2.75,3,3.15,3.5]
               cs = plt.contourf(xi,yi,np.abs(np.squeeze(vals)),med_range,cmap='jet',extend='both')

            # Add the grid
            m.drawparallels(np.arange(30., 46., 5.), labels=[1,0,0,0], fontsize=10)
            m.drawmeridians(np.arange(-20., 40., 10.), labels=[0,0,0,1], fontsize=10)

            # Land from bathymetry or mesh mask file
            contourf1 = plt.contourf(xi,yi,np.abs(np.squeeze(vals_bathy)),[0.00000,0.00000001], colors='gray')
            contour1 = plt.contour(xi,yi,np.abs(np.squeeze(vals_bathy)),0.00000,colors='black')

            # Plot the legend 
            cbar = m.colorbar(cs, location='bottom', pad="10%",extend='max')
            if ampha_var == 'Amp' or ampha_var == 'Amp_Ar' or ampha_var == 'Amp_Pal':
               bar_label_string=var_2d[:2]+' Amplitude ['+var_2d_udm+']'
            elif ampha_var == 'Pha' or ampha_var == 'Pha_Ar' or ampha_var == 'Pha_Pal':
               bar_label_string=var_2d[:2]+' Phase [deg]'
            elif ampha_var == 'AmpPha' or ampha_var == 'AmpPha_Ag' or ampha_var == 'AtlBox':
               bar_label_string=var_2d[:2]+' Amplitude'+' ['+var_2d_udm+']'
            cbar.set_label(bar_label_string)

            # Add Max Amp value computed on the whole domain
            if ampha_var == 'Amp' or ampha_var == 'AmpPha' or ampha_var == 'AmpPha_Ag' or ampha_var == 'Amp_Pal' or ampha_var == 'AtlBox':
               thresh = 0.0000
               mask = np.abs(vals) == thresh
               vals_ma = np.ma.masked_where(mask, vals)
               vals_max=np.amax(abs(vals_ma))
               vals_min=0
               text_max_x,text_max_y= m(32,29.0)
               plt.text(text_max_x,text_max_y,'max='+str(round(vals_max,1))+var_2d_udm, fontsize=12)

            plt.savefig(plotname)
            plt.clf()




# TPXO Amp Pha
if ampha_tpxo == 1:
   print ('# 2D VAR: Amplitude and Phase from TPXO9 model')
   # Build the pathname of the nc file and open it 
   idx_comp = 0
   for tidal_component in ('m2','k1','o1','s2','p1','n2','q1','k2'):
       if tidal_component == 'm2':
          med_range=[0,5,10,15,20,25]
       elif tidal_component == 'k1' or tidal_component == 's2':
          med_range=[0,2.5,5,7.5,10,12.5,15]
       elif tidal_component == 'o1' or tidal_component == 'p1' or tidal_component == 'n2':
          med_range=[0,1,2,3,4,5]
       elif tidal_component == 'k2':
          med_range=[0,0.5,1,1.5,2,2.5,3]
       elif tidal_component == 'q1':
          med_range=[0,0.25,0.5,0.75,1,1.25,1.5,1.75,2]

       # Adjust the palettes for AtlBox case
       if ampha_var == 'AtlBox':
               if tidal_component == 'm2':
                   #med_range=[0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150]
                   med_range=[40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150]
               elif tidal_component == 's2':
                   #med_range=[0,2.5,5,7.5,10,12.5,15,17.5,20,22.5,25,27.5,30,32.5,35,37.5,40,42.5,45,47.5,50]
                   med_range=[15,17.5,20,22.5,25,27.5,30,32.5,35,37.5,40,42.5,45,47.5,50]
               elif tidal_component == 'n2':
                   #med_range=[0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34]
                   med_range=[8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34]
               elif tidal_component == 'k2':
                   #med_range=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
                   med_range=[4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10,10.5,11,11.5,12,12.5,13,13.5,14,14.5,15,15.5,16]
               elif tidal_component == 'k1' or tidal_component == 'o1':
                   #med_range=[0,1,2,3,4,5,6,7,8,9,10]
                   med_range=[3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10]
               elif tidal_component == 'p1' or tidal_component == 'q1':
                   #med_range=[0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5]
                   med_range=[0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.15,2.5,2.75,3,3.15,3.5]

       nc2open=tpxo9_path+tpxo9_fileprename+tidal_component+tpxo9_postname+'.nc'
       print ('Input file = ',nc2open)
       model = NC.Dataset(nc2open,'r')
       # Read lat, lon and fild values 
       lons_tpxo = model.variables['lon_z'][:]
       lats_tpxo = model.variables['lat_z'][:]
       lons_tpxo=np.squeeze(lons_tpxo)
       lats_tpxo=np.squeeze(lats_tpxo)
       name_var='vals_tpxo_A'+tidal_component
       globals()[name_var] = (np.sqrt(model.variables['hRe'][:]**2+model.variables['hIm'][:]**2))/10 # Want cm not millimiters!
       name_pha='vals_tpxo_P'+tidal_component
       globals()[name_pha] = np.arctan2(-model.variables['hIm'][:],model.variables['hRe'][:])*180/np.pi
       #
       field2plot=globals()[name_var]
       P_vals_tpxo=globals()[name_pha]
       field2plot=np.transpose(field2plot)
       P_vals_tpxo=np.transpose(P_vals_tpxo)

       land_value = 0.00000
       if idx_comp == 0:
          vals_tpxo_Am2=np.transpose(vals_tpxo_Am2)      
       thresh = 0.0000
       mask = np.abs(vals_tpxo_Am2) == thresh
       vals_tpxo_ma = np.ma.masked_where(mask, vals_tpxo_Am2)

       # Plot the map and save in the path/name
       var_2d=str(globals()[name_var])
       var_2d_udm=field_2d_units[idx_comp]
       plotname=workdir_path+tpxo9_fileprename+'2D_'+'0_'+dates_label+tpxo9_postname+'_'+tidal_component+'_Amp_nopha.jpg'
       # Adjust the plot name for AtlBox case
       if ampha_var == 'AtlBox':
          plotname=workdir_path+tpxo9_fileprename+'2D_'+'0_'+dates_label+tpxo9_postname+'_'+tidal_component+'_AmpPha_AtlBox.jpg'
       print ('Plot path/name: ',plotname)

       fig, ax = plt.subplots(1, 3, sharey=True,figsize=(14,4)) 
       fig.subplots_adjust(wspace=0)
       plt.rc('font', size=12)
       # Plot Title
       plt.suptitle (tidal_component+'Amplitude ['+var_2d_udm+'] --- TPXO9 ---'+dates_label)
       # Read the coordinates for the plot 
       lon_0 = lons_tpxo.mean()
       llcrnrlon = lons_tpxo.min()
       urcrnrlon = lons_tpxo.max()
       lat_0 = lats_tpxo.mean()
       llcrnrlat = lats_tpxo.min()
       urcrnrlat = lats_tpxo.max()
       # Plot the frame to the map
       plt.rcParams["axes.linewidth"]  = 1.25
       vals_tpxo_max=np.amax(field2plot)
       vals_tpxo_min=np.amin(field2plot)
       lowth=0

       plt.subplot(1,3,1)
       plt.ylim(30, 46)
       plt.xlim(340, 360)
       if med_range[0] == 0:
          factor_contour1 = plt.contourf(lons_tpxo,lats_tpxo,np.squeeze(field2plot),med_range,extend='max',cmap='jet') 
       else:
          factor_contour1 = plt.contourf(lons_tpxo,lats_tpxo,np.squeeze(field2plot),med_range,extend='both',cmap='jet')
       contour1p = plt.contour(lons_tpxo,lats_tpxo,np.squeeze(P_vals_tpxo),[0,30,60,90,120,150,180,210,240,270,300,330,360,390],colors='white')
       contour1n = plt.contour(lons_tpxo,lats_tpxo,np.squeeze(P_vals_tpxo),[-390,-360,-330,-300,-270,-240,-210,-180,-150,-120,-90,-60,-30],colors='white',linestyles='dashed')
       contourf1 = plt.contourf(lons_tpxo,lats_tpxo,np.abs(np.squeeze(vals_tpxo_Am2)),[land_value,0.00000001], colors='gray')
       contour1 = plt.contour(lons_tpxo,lats_tpxo,np.abs(np.squeeze(vals_tpxo_Am2)),land_value,colors='black')
       ax[0].set_xticks(np.arange(340,360,10))
       plt.tick_params(labelcolor='none', top=False, bottom=True, left=False, right=False)
       plt.grid(color='black',linestyle='--')
       #
       plt.subplot(1,3,(2,3))
       plt.ylim(30, 46)
       plt.xlim(0, 40)
       plt.tick_params(labelcolor='none', top=False, bottom=True, left=False, right=False)
       if med_range[0] == 0:
          factor_contour2 = plt.contourf(lons_tpxo,lats_tpxo,np.squeeze(field2plot),med_range,extend='max',cmap='jet')
       else:
          factor_contour2 = plt.contourf(lons_tpxo,lats_tpxo,np.squeeze(field2plot),med_range,extend='both',cmap='jet')
       contour2p = plt.contour(lons_tpxo,lats_tpxo,np.squeeze(P_vals_tpxo),[0,30,60,90,120,150,180,210,240,270,300,330,360,390],colors='white')
       contour2n = plt.contour(lons_tpxo,lats_tpxo,np.squeeze(P_vals_tpxo),[-390,-360,-330,-300,-270,-240,-210,-180,-150,-120,-90,-60,-30],colors='white',linestyles='dashed')
       contourf2 = plt.contourf(lons_tpxo,lats_tpxo,np.abs(np.squeeze(vals_tpxo_Am2)),[land_value,0.00000001], colors='gray')
       contour2 = plt.contour(lons_tpxo,lats_tpxo,np.abs(np.squeeze(vals_tpxo_Am2)),land_value,colors='black')
       ax[1].get_yaxis().set_visible(False)
       plt.grid(color='black',linestyle='--')
       cbar = plt.colorbar(factor_contour2, extend='max',label=tidal_component+' Amplitude [cm]') 
       # Save and close 
       plt.savefig(plotname)
       plt.clf()

       idx_comp=idx_comp+1

# TPXO Phase maps
if pha_tpxo == 1:
   print ('# 2D VAR: Phase maps from TPXO9 model')
   # Build the pathname of the nc file and open it 
   idx_comp = 0
   for tidal_component in ('m2','k1','o1','s2','p1','n2','q1','k2'):

       nc2open=tpxo9_path+tpxo9_fileprename+tidal_component+tpxo9_postname+'.nc'
       print ('Input file = ',nc2open)
       model = NC.Dataset(nc2open,'r')
       # Read lat, lon and fild values 
       lons_tpxo = model.variables['lon_z'][:]
       lats_tpxo = model.variables['lat_z'][:]
       lons_tpxo=np.squeeze(lons_tpxo)
       lats_tpxo=np.squeeze(lats_tpxo)
       name_var='vals_tpxo_A'+tidal_component
       globals()[name_var] = (np.sqrt(model.variables['hRe'][:]**2+model.variables['hIm'][:]**2))/10 # Want cm not millimiters!
       name_pha='vals_tpxo_P'+tidal_component
       globals()[name_pha] = np.arctan2(-model.variables['hIm'][:],model.variables['hRe'][:])*180/np.pi
       #
       field2plot=globals()[name_var]
       P_vals_tpxo=globals()[name_pha]
       field2plot=np.transpose(field2plot)
       P_vals_tpxo=np.transpose(P_vals_tpxo)
 
       # Pha range:
       P_vals_tpxo=np.where(P_vals_tpxo>=360.0,P_vals_tpxo-360,P_vals_tpxo)
       P_vals_tpxo=np.where(P_vals_tpxo>=360.0,P_vals_tpxo-360,P_vals_tpxo)
       P_vals_tpxo=np.where(P_vals_tpxo<0.0,P_vals_tpxo+360,P_vals_tpxo)
       P_vals_tpxo=np.where(P_vals_tpxo<0.0,P_vals_tpxo+360,P_vals_tpxo)

       # Set the land 
       land_value = 0.00000
       if idx_comp == 0:
          vals_tpxo_Am2=np.transpose(vals_tpxo_Am2)
       thresh = 0.0000
       mask = np.abs(vals_tpxo_Am2) == thresh
       vals_tpxo_ma = np.ma.masked_where(mask, vals_tpxo_Am2)

       # Set the correct var to be plottet in this case
       name_var=name_pha
       field2plot=P_vals_tpxo

       # Plot the map and save in the path/name
       var_2d=str(globals()[name_var])
       var_2d_udm='deg'
       plotname=workdir_path+tpxo9_fileprename+'2D_'+'0_'+dates_label+tpxo9_postname+'_'+tidal_component+'_Pha.jpg'
       print ('Plot path/name: ',plotname)

       fig, ax = plt.subplots(1, 3, sharey=True,figsize=(14,4))
       fig.subplots_adjust(wspace=0)
       plt.rc('font', size=12)
       # Plot Title
       plt.suptitle (tidal_component+' Amplitude --- TPXO9 ---'+dates_label)
       # Read the coordinates for the plot 
       lon_0 = lons_tpxo.mean()
       llcrnrlon = lons_tpxo.min()
       urcrnrlon = lons_tpxo.max()
       lat_0 = lats_tpxo.mean()
       llcrnrlat = lats_tpxo.min()
       urcrnrlat = lats_tpxo.max()
       # Plot the frame to the map
       plt.rcParams["axes.linewidth"]  = 1.25
       vals_tpxo_max=np.amax(field2plot)
       vals_tpxo_min=np.amin(field2plot)
       lowth=0

       plt.subplot(1,3,1)
       plt.ylim(30, 46)
       plt.xlim(340, 360)
       factor_contour1 = plt.contourf(lons_tpxo,lats_tpxo,np.abs(np.squeeze(P_vals_tpxo)),[0,15,30,45,60,75,90,105,120,135,150,165,180,195,210,225,240,255,270,285,300,315,330,345,360],cmap='hsv')
       contourf1 = plt.contourf(lons_tpxo,lats_tpxo,np.abs(np.squeeze(vals_tpxo_Am2)),[land_value,0.00000001], colors='gray')
       contour1 = plt.contour(lons_tpxo,lats_tpxo,np.abs(np.squeeze(vals_tpxo_Am2)),land_value,colors='black')
       ax[0].set_xticks(np.arange(340,360,10))
       plt.tick_params(labelcolor='none', top=False, bottom=True, left=False, right=False)
       plt.grid(color='black',linestyle='--')
       #
       plt.subplot(1,3,(2,3))
       plt.ylim(30, 46)
       plt.xlim(0, 40)
       plt.tick_params(labelcolor='none', top=False, bottom=True, left=False, right=False)
       factor_contour2 = plt.contourf(lons_tpxo,lats_tpxo,np.abs(np.squeeze(P_vals_tpxo)),[0,15,30,45,60,75,90,105,120,135,150,165,180,195,210,225,240,255,270,285,300,315,330,345,360],cmap='hsv')
       contourf2 = plt.contourf(lons_tpxo,lats_tpxo,np.abs(np.squeeze(vals_tpxo_Am2)),[land_value,0.00000001], colors='gray')
       contour2 = plt.contour(lons_tpxo,lats_tpxo,np.abs(np.squeeze(vals_tpxo_Am2)),land_value,colors='black')
       ax[1].get_yaxis().set_visible(False)
       plt.grid(color='black',linestyle='--')
       cbar = plt.colorbar(factor_contour2, extend='max',label=tidal_component+' Phase [deg]')
       # Save and close 
       plt.savefig(plotname)
       plt.clf()

       idx_comp=idx_comp+1



# TPXO Amp Pha in North Atlantic Ocean (only M2 component, to be extended..)
if ampha_tpxo_atlamph == 1:
   print ('# 2D VAR: Amplitude and Phase in the North Atlantic Ocean from TPXO9 model')
   # Build the pathname of the nc file and open it 
   idx_comp = 0
   for tidal_component in ('m2','m2'):
       if tidal_component == 'm2':
          med_range=[0,20,40,60,80,100,120,140,160,180,200]

       nc2open=tpxo9_path+tpxo9_fileprename+tidal_component+tpxo9_postname+'.nc'
       print ('Input file = ',nc2open)
       model = NC.Dataset(nc2open,'r')
       # Read lat, lon and fild values 
       lons_tpxo = model.variables['lon_z'][:]
       lats_tpxo = model.variables['lat_z'][:]
       lons_tpxo=np.squeeze(lons_tpxo)
       lats_tpxo=np.squeeze(lats_tpxo)
       name_var='vals_tpxo_A'+tidal_component
       globals()[name_var] = (np.sqrt(model.variables['hRe'][:]**2+model.variables['hIm'][:]**2))/10 # Want cm not millimiters!
       name_pha='vals_tpxo_P'+tidal_component
       globals()[name_pha] = np.arctan2(-model.variables['hIm'][:],model.variables['hRe'][:])*180/np.pi
       #
       field2plot=globals()[name_var]
       P_vals_tpxo=globals()[name_pha]
       field2plot=np.transpose(field2plot)
       P_vals_tpxo=np.transpose(P_vals_tpxo)

       land_value = 0.00000
       if idx_comp == 0:
          vals_tpxo_Am2=np.transpose(vals_tpxo_Am2)
       thresh = 0.0000
       mask = np.abs(vals_tpxo_Am2) == thresh
       vals_tpxo_ma = np.ma.masked_where(mask, vals_tpxo_Am2)

       # Plot the map and save in the path/name
       var_2d=str(globals()[name_var])
       var_2d_udm=field_2d_units[idx_comp]
       plotname=workdir_path+tpxo9_fileprename+'2D_'+'0_'+dates_label+tpxo9_postname+'_'+tidal_component+'_Amppha_AtlAmph.jpg'
       print ('Plot path/name: ',plotname)

       fig, ax = plt.subplots(1, 1,figsize=(8,4)) 
       plt.rc('font', size=12)
       # Plot Title
       plt.suptitle (tidal_component+' Amplitude ['+var_2d_udm+'] --- TPXO9 ---'+dates_label)
       # Read the coordinates for the plot 
       lon_0 = lons_tpxo.mean()
       llcrnrlon = lons_tpxo.min()
       urcrnrlon = lons_tpxo.max()
       lat_0 = lats_tpxo.mean()
       llcrnrlat = lats_tpxo.min()
       urcrnrlat = lats_tpxo.max()
       # Plot the frame to the map
       plt.rcParams["axes.linewidth"]  = 1.25
       vals_tpxo_max=np.amax(field2plot)
       vals_tpxo_min=np.amin(field2plot)
       lowth=0

       plt.subplot(1,1,1) 
       plt.ylim(30,60) 
       plt.xlim(300,360) # Or (0,60) for Med  
       factor_contour1 = plt.contourf(lons_tpxo,lats_tpxo,np.squeeze(field2plot),med_range,extend='max',cmap='jet')
       contour1p = plt.contour(lons_tpxo,lats_tpxo,np.squeeze(P_vals_tpxo),[0,30,60,90,120,150,180,210,240,270,300,330,360,390],colors='white')
       contour1n = plt.contour(lons_tpxo,lats_tpxo,np.squeeze(P_vals_tpxo),[-390,-360,-330,-300,-270,-240,-210,-180,-150,-120,-90,-60,-30],colors='white',linestyles='dashed')
       contour = plt.contour(lons_tpxo,lats_tpxo,np.abs(np.squeeze(field2plot)),[land_value,0.00000001], colors='black')
       contourf = plt.contourf(lons_tpxo,lats_tpxo,np.abs(np.squeeze(field2plot)),[land_value,0.00000001], colors='gray')
       plt.tick_params(labelcolor='black', top=False, bottom=True, left=True, right=False)
       plt.grid(color='black',linestyle='--')
       #
       cbar = plt.colorbar(factor_contour1, extend='max', label=tidal_component+' Amplitude [cm]')
       #plt.clabel(contour1p, fmt = '%2.1d', colors = 'k', fontsize=12)
       #plt.clabel(contour1n, fmt = '%2.1d', colors = 'k', fontsize=12)
       # Save and close 
       plt.savefig(plotname)
       plt.clf()

       idx_comp=idx_comp+1

########
# Intepolation of global TPXO9 fields on MED24 grid (1/30 --> 1/24)
if tpxo2eas_flag ==  1: 
        # Read new grid structure from MED24 mesh_mask files
        # Open the mesh mask file and read lat/lon 
        nc2open3=model_bathy
        model3 = NC.Dataset(nc2open3,'r')
        nc2open4=model_meshmask # mesh mask
        model4 = NC.Dataset(nc2open4,'r')

        lons = model3.variables[longitude_name][:]
        lats = model3.variables[latitude_name][:]
        vals_bathy=model3.variables['Bathymetry'][:]
        vals_land=model4.variables['tmask'][0,0,:,:]
        vals_land=np.squeeze(vals_land)

        x=lons[0] # Array 1D longitudes
        y=[ el[0] for el in lats] #  Array 1D latitudes
        y=np.asarray(y)

        # Loop on TPXO9 fields to be extracted and plot range setting  
        print ('Reading Amplitude and Phase fields from TPXO9 model...')
        # Build the pathname of the nc file and open it 
        idx_comp = 0
        # Loop on tidal components
        for tidal_component in ['m2','k1','o1','s2','p1','n2','q1','k2']: 
            if tidal_component == 'm2':
               med_range=[0,5,10,15,20,25]
            elif tidal_component == 'k1' or tidal_component == 's2':
               med_range=[0,2.5,5,7.5,10,12.5,15]
            elif tidal_component == 'o1' or tidal_component == 'p1' or tidal_component == 'n2':
               med_range=[0,1,2,3,4,5]
            elif tidal_component == 'k2':
               med_range=[0,0.5,1,1.5,2,2.5,3]
            elif tidal_component == 'q1':
               med_range=[0,0.25,0.5,0.75,1,1.25,1.5,1.75,2]

            # Adjust the palettes for AtlBox case
            if ampha_var == 'AtlBox':
               if tidal_component == 'm2':
                   #med_range=[0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150]
                   med_range=[40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150]
               elif tidal_component == 's2':
                   #med_range=[0,2.5,5,7.5,10,12.5,15,17.5,20,22.5,25,27.5,30,32.5,35,37.5,40,42.5,45,47.5,50]
                   med_range=[15,17.5,20,22.5,25,27.5,30,32.5,35,37.5,40,42.5,45,47.5,50]
               elif tidal_component == 'n2':
                   #med_range=[0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34]
                   med_range=[8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34]
               elif tidal_component == 'k2':
                   #med_range=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
                   med_range=[4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10,10.5,11,11.5,12,12.5,13,13.5,14,14.5,15,15.5,16]
               elif tidal_component == 'k1' or tidal_component == 'o1':
                   #med_range=[0,1,2,3,4,5,6,7,8,9,10]
                   med_range=[3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10]
               elif tidal_component == 'p1' or tidal_component == 'q1':
                   #med_range=[0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5]
                   med_range=[0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.15,2.5,2.75,3,3.15,3.5]

            nc2open=tpxo9_path+tpxo9_fileprename+tidal_component+tpxo9_postname+'.nc'
            print ('Input file = ',nc2open)
            model = NC.Dataset(nc2open,'r')

            # Read lat, lon fields only from the first file 
            if idx_comp == 0:
                lons_tpxo = model.variables['lon_z'][:]
                lats_tpxo = model.variables['lat_z'][:]
                lons_tpxo=np.squeeze(lons_tpxo)
                lats_tpxo=np.squeeze(lats_tpxo)

            # Compute Amplitudes and Phases from Re and Im parts             
            name_var='vals_tpxo_A'+tidal_component
            globals()[name_var] = (np.sqrt(model.variables['hRe'][:]**2+model.variables['hIm'][:]**2))/10 # Want cm not millimiters!
            name_pha='vals_tpxo_P'+tidal_component
            globals()[name_pha] = np.arctan2(-model.variables['hIm'][:],model.variables['hRe'][:])*180/np.pi

            field2interp=globals()[name_var]
            P_vals_tpxo=globals()[name_pha]
            field2interp=np.reshape(field2interp,len(lons_tpxo)*len(lats_tpxo))
            P_vals_tpxo=np.reshape(P_vals_tpxo,len(lons_tpxo)*len(lats_tpxo))

            # Interpolation from the old to the new grid 
            f_tmp = interpolate.interp2d(lons_tpxo,lats_tpxo,field2interp)
            p_tmp = interpolate.interp2d(lons_tpxo,lats_tpxo,P_vals_tpxo)
            var_new = f_tmp(x+360,y)
            pha_new = p_tmp(x+360,y)

            # Plot the map and save in the path/name
            plotname=workdir_path+'interp_tpxo2med24_'+tidal_component+'.jpg'
            # Adjust the plot name for AtlBox case
            if ampha_var == 'AtlBox':
               plotname=workdir_path+'interp_tpxo2med24_'+tidal_component+'AtlBox.jpg'

            print ('Plot path/name: ',plotname)

            fig, ax = plt.subplots(1, 3, sharey=True,figsize=(14,4)) 
            fig.subplots_adjust(wspace=0)
            plt.rc('font', size=12)
            # Plot Title
            plt.suptitle ('TPXO9 Amplitude and Phase interp. on MED24 grid --- '+tidal_component)
            # Plot the frame to the map
            plt.rcParams["axes.linewidth"]  = 1.25
            # Plot the map and add the amphidromes
            plt.subplot(1,3,1)
            plt.ylim(30, 46)
            plt.xlim(-20, 0)
            if med_range[0] == 0:
               cs = plt.contourf(x,y,np.abs(np.squeeze(var_new)),med_range,cmap='jet',extend='max')
            else:
               cs = plt.contourf(x,y,np.abs(np.squeeze(var_new)),med_range,cmap='jet',extend='both')
            contour1p = plt.contour(x,y,np.squeeze(pha_new),[0,30,60,90,120,150,180,210,240,270,300,330,360,390],colors='white')
            contour1n = plt.contour(x,y,np.squeeze(pha_new),[-390,-360,-330,-300,-270,-240,-210,-180,-150,-120,-90,-60,-30],colors='white',linestyles='dashed')
            contourf1 = plt.contourf(x,y,np.abs(np.squeeze(vals_bathy)),[0.00000,0.00000001], colors='gray') 
            contourf = plt.contour(x,y,np.abs(np.squeeze(vals_bathy)),[0.00000,0.00000001], colors='black')
            ax[0].set_xticks(np.arange(340,360,10))
            plt.tick_params(labelcolor='none', top=False, bottom=True, left=False, right=False)
            plt.grid(color='black',linestyle='--')
            #
            var_new2 = f_tmp(x,y)
            pha_new2 = p_tmp(x,y)
            #
            plt.subplot(1,3,(2,3))
            plt.ylim(30, 46)
            plt.xlim(0, 40)
            plt.tick_params(labelcolor='none', top=False, bottom=True, left=False, right=False)
            if med_range[0] == 0:
               cs = plt.contourf(x,y,np.abs(np.squeeze(var_new2)),med_range,cmap='jet',extend='max')
            else:
               cs = plt.contourf(x,y,np.abs(np.squeeze(var_new2)),med_range,cmap='jet',extend='both')
            contour1p = plt.contour(x,y,np.squeeze(pha_new2),[0,30,60,90,120,150,180,210,240,270,300,330,360,390],colors='white')
            contour1n = plt.contour(x,y,np.squeeze(pha_new2),[-390,-360,-330,-300,-270,-240,-210,-180,-150,-120,-90,-60,-30],colors='white',linestyles='dashed')
            contourf1 = plt.contourf(x,y,np.abs(np.squeeze(vals_bathy)),[0.00000,0.00000001], colors='gray') 
            contourf = plt.contour(x,y,np.abs(np.squeeze(vals_bathy)),[0.00000,0.00000001], colors='black')
            ax[1].get_yaxis().set_visible(False)
            plt.grid(color='black',linestyle='--')
            bar = plt.colorbar(cs,extend='max', label=tidal_component+' Amplitude [cm]')

            # Save and close 
            plt.savefig(plotname)
            plt.clf()

            idx_comp=idx_comp+1


# Intepolation of global TPXO9 fields on MED24 grid and differences  (1/30 --> 1/24)
if diff_tpxoeas_flag ==  1:
            # Read new grid structure from MED24 mesh_mask files
            # Open the mesh mask file and read lat/lon 
            nc2open3=model_bathy
            model3 = NC.Dataset(nc2open3,'r')
            nc2open4=model_meshmask # mesh mask
            model4 = NC.Dataset(nc2open4,'r')

            lons = model3.variables[longitude_name][:]
            lats = model3.variables[latitude_name][:]
            vals_bathy=model3.variables['Bathymetry'][:]
            vals_land=model4.variables['tmask'][0,0,:,:]
            vals_land=np.squeeze(vals_land)

            x=lons[0] # Array 1D longitudes
            y=[ el[0] for el in lats] #  Array 1D latitudes
            y=np.asarray(y)

            # RMSE table ini
            RMSE_file = open(workdir_path+"RMSE_stats.txt","w")
            RMSE_file.write('Root Mean Square Error (Amplitude differences)'+'\n')
            # Read Amp and Pha from MED24
            idx_comp = 0 
            for idx_2d in range(0,len(field_2d_name)):
                     var_2d=field_2d_name[idx_2d]
                     var_2d_udm=field_2d_units[idx_2d]
         
                     P_var_2d=P_field_2d_name[idx_2d]
                     P_var_2d_udm=P_field_2d_units[idx_2d]
         
                     # Build the path/name of the nc file and open it 
                     nc2open=model_path+model_fileprename+'2D_'+'0_'+'sossheig'+'_'+dates_label+'_'+model_postname+'.nc'
                     print ('Input file = ',nc2open)
                     model_med24 = NC.Dataset(nc2open,'r')
         
                     # Read field values 
                     vals = model_med24.variables[var_2d][:]*100 # Want cm not meters!
                     P_vals = model_med24.variables[P_var_2d][:]*(np.pi/180) # Want radians not deg
         
                     thresh = 0.0000
                     mask = np.abs(vals) == thresh
                     vals_ma = np.ma.masked_where(mask, vals)
                     # VALS_MA and P_VALS
                     

                     # Choose TPXO9 fields to be used and select the range  
                     print ('Reading Amplitude and Phase fields from TPXO9 model...')
                     # Build the pathname of the nc file and open it 
                     # Loop on tidal components
                     tidal_comp_name_tpxo=['m2','k1','o1','s2','p1','n2','q1','k2']
                     tidal_component=tidal_comp_name_tpxo[idx_comp]
                     if tidal_component == 'm2':
                        med_range=[-4.0,-3.5,-2.5,-1.5,-0.5,0.5,1.5,2.5,3.5,4.0]
                     elif tidal_component == 's2':
                        med_range=[-3.5,-2.5,-1.5,-0.5,0.5,1.5,2.5,3.5]
                     elif tidal_component == 'k1':
                        med_range=[-1.4,-1.0,-0.6,-0.2,0.2,0.6,1.0,1.4]
                     elif tidal_component == 'k2':
                        med_range=[-1.4,-1.0,-0.6,-0.2,0.2,0.6,1.0,1.4]
                     else:
                        med_range=[-1.4,-1.0,-0.6,-0.2,0.2,0.6,1.0,1.4]
     
                     nc2open=tpxo9_path+tpxo9_fileprename+tidal_component+tpxo9_postname+'.nc'
                     print ('Input file = ',nc2open)
                     model = NC.Dataset(nc2open,'r')
     
                     # Read lat, lon fields only from the first file 
                     if idx_comp == 0:
                        lons_tpxo = model.variables['lon_z'][:]
                        lats_tpxo = model.variables['lat_z'][:]
                        lons_tpxo=np.squeeze(lons_tpxo)
                        lats_tpxo=np.squeeze(lats_tpxo)
     
                     # Compute Amplitudes and Phases from Re and Im parts             
                     name_var='vals_tpxo_A'+tidal_component
                     globals()[name_var] = (np.sqrt(model.variables['hRe'][:]**2+model.variables['hIm'][:]**2))/10 # Want cm not millimiters!
                     name_pha='vals_tpxo_P'+tidal_component
                     globals()[name_pha] = np.arctan2(-model.variables['hIm'][:],model.variables['hRe'][:]) #*180/np.pi
     
                     field2interp=globals()[name_var]
                     P_vals_tpxo=globals()[name_pha]
     
                     # Sea over land
                     thresh2interp = 0.0000
                     mask2interp = np.abs(field2interp) == thresh2interp
                     field2interp_ma = np.ma.masked_where(mask2interp,field2interp)
                     input_matrix=field2interp_ma 
                     nloop=1
                     infill_value = input_matrix.fill_value
                     output_matrix = ma.copy(input_matrix)  # using ma.copy to prevent future field modifications
                     if np.sum(output_matrix.mask) == 0:  # nothing to fill
                        print('WARNING. Field does not have any land points or wrong field type. Exiting.', file=sys.stderr)
                     else:
                        # iterations loop
                        for loop in range(nloop):
                           # Create a nD x 8 matrix in which, last dimension fixed, the other dimensions
                           #  contains values that are shifted in one of the 8 possible directions
                           # of the last two axes compared to the original matrix
                           shift_matrix = ma.array(np.empty(shape=((8,) + output_matrix.shape)),
                                                   mask=True, fill_value=infill_value, dtype=float)
                           # up shift
                           shift_matrix[0, ..., : -1,:] = output_matrix[..., 1:,:]
                           # down shift
                           shift_matrix[1, ..., 1:, :] = output_matrix[..., : -1, :]
                           # left shift
                           shift_matrix[2, ..., :, : -1] = output_matrix[..., :, 1:]
                           # right shift
                           shift_matrix[3, ..., :, 1:] = output_matrix[..., :, : -1]
                           # up-left shift
                           shift_matrix[4, ..., : -1, : -1] = output_matrix[..., 1:, 1:]
                           # up-right shift
                           shift_matrix[5, ..., : -1, 1:] = output_matrix[..., 1:, : -1]
                           # down-left shift
                           shift_matrix[6, ..., 1:, : -1] = output_matrix[..., : -1, 1:]
                           # down-right shift
                           shift_matrix[7, ..., 1:, 1:] = output_matrix[..., : -1, : -1]
                           # Mediate the shift matrix among the third dimension
                           mean_matrix = ma.mean(shift_matrix, 0)
                           # Replace input masked points with new ones belonging to the mean matrix
                           output_matrix = ma.array(np.where(mean_matrix.mask + output_matrix.mask, mean_matrix, output_matrix),
                                                    mask=mean_matrix.mask, fill_value=infill_value, dtype=float)
                           output_matrix = ma.masked_where(mean_matrix.mask, output_matrix)
                           if np.sum(output_matrix.mask) == 0:  # nothing more to flood
                               print('WARNING. Field does not have anymore land points,', str(loop + 1),
                                     'steps were sufficient to flood it completely.', file=sys.stderr)
                               break
                     field2interp=output_matrix
                     #
                     field2interp=np.reshape(field2interp,len(lons_tpxo)*len(lats_tpxo))
                     P_vals_tpxo=np.reshape(P_vals_tpxo,len(lons_tpxo)*len(lats_tpxo))

                     # Interpolation from the old to the new grid 
                     f_tmp = interpolate.interp2d(lons_tpxo,lats_tpxo,field2interp)
                     p_tmp = interpolate.interp2d(lons_tpxo,lats_tpxo,P_vals_tpxo)
                     var_new = f_tmp(x+360,y)
                     pha_new = p_tmp(x+360,y)

                     # Plot the map and save in the path/name
                     plotname=workdir_path+'diff_tpxomed24_'+tidal_component+'.jpg'
                     print ('Plot path/name: ',plotname)
      
                     fig, ax = plt.subplots(1, 3, sharey=True,figsize=(14,4))
                     fig.subplots_adjust(wspace=0)
                     plt.rc('font', size=12)
                     # Plot Title
                     plt.suptitle ('Amplitude Diff: '+model_postname+' - TPXO9  --- MED24 grid --- '+tidal_component)

                     # Plot the frame to the map
                     plt.rcParams["axes.linewidth"]  = 1.25
                     # Plot the map and add the amphidromes
                     plt.subplot(1,3,1)
                     plt.ylim(30, 46)
                     plt.xlim(-20, 0)
                     cs = plt.contourf(x,y,np.abs(np.squeeze(vals))-np.abs(np.squeeze(var_new)),med_range,cmap='bwr',extend='both')
                     contourf1 = plt.contourf(x,y,np.abs(np.squeeze(vals_bathy)),[0.00000,0.00000001], colors='gray') 
                     contour2 = plt.contour(x,y,np.squeeze(vals_bathy),0.00000,colors='black')
                     ax[0].set_xticks(np.arange(340,360,10))
                     plt.tick_params(labelcolor='none', top=False, bottom=True, left=False, right=False)
                     plt.grid(color='black',linestyle='--')
                     max_n=np.max(np.abs(np.abs(np.squeeze(vals))-np.abs(np.squeeze(var_new))))

                     # RMSE diff Amp 4 bathy > 1000m
                     sq_vals_bathy=vals_bathy[0,:,:] 
                     rmse_diffamp_1000=[]
                     rmse_count=0
                     for lo_idx in range (0,len(x)):
                      if x[lo_idx]>-5 and x[lo_idx]<0:
                         for la_idx in range (0,len(y)):
                           if y[la_idx]>30 and y[la_idx]<40:
                             if sq_vals_bathy[la_idx][lo_idx]>bathylim4RMSE:
                                rmse_diffamp_1000.append(abs(vals[la_idx][lo_idx])-abs(var_new[la_idx][lo_idx]))
                                rmse_count=rmse_count+1

                     #
                     var_new2 = f_tmp(x,y)
                     pha_new2 = p_tmp(x,y)
                     #
                     plt.subplot(1,3,(2,3))
                     plt.ylim(30, 46)
                     plt.xlim(0, 40)
                     plt.tick_params(labelcolor='none', top=False, bottom=True, left=False, right=False)
                     cs = plt.contourf(x,y,np.abs(np.squeeze(vals))-np.abs(np.squeeze(var_new2)),med_range,cmap='bwr',extend='both')
                     contourf1 = plt.contourf(x,y,np.abs(np.squeeze(vals_bathy)),[0.00000,0.00000001], colors='gray') # 0.00000001
                     contour2 = plt.contour(x,y,np.squeeze(vals_bathy),0.00000,colors='black')
                     ax[1].get_yaxis().set_visible(False)
                     plt.grid(color='black',linestyle='--')
                     max_p=np.max(np.abs(np.abs(np.squeeze(vals))-np.abs(np.squeeze(var_new2))))
                     cbar = plt.colorbar(cs,extend='both')
                     cbar.set_label('Amplitude diff [cm]')

                     # RMSE
                     for lo_idx in range (0,len(x)):
                      if x[lo_idx]>0 and x[lo_idx]<40: 
                         for la_idx in range (0,len(y)):
                          if y[la_idx]>30 and y[la_idx]<46:
                             if sq_vals_bathy[la_idx][lo_idx]>bathylim4RMSE:
                                rmse_diffamp_1000.append(abs(vals[la_idx][lo_idx])-abs(var_new2[la_idx][lo_idx]))
                                rmse_count=rmse_count+1
                     rmse_val_ar=np.power(rmse_diffamp_1000,2)
                     rmse_val=np.sum(rmse_val_ar)
                     rmse_val=np.sqrt(rmse_val/rmse_count)
                     RMSE_file.write(tidal_component+' RMSE_TOT: '+str(rmse_val)+'\n')


                     # Print Max Amp diff computed on the whole domain
                     #mask_val=np.ma.masked_where(mask, vals)
                     #mask_var=np.ma.masked_where(mask, var_new2)
                     #vals_max=np.percentile(abs(np.abs(np.squeeze(mask_val))-np.abs(np.squeeze(mask_var))),98)
                     
                     #print ('Max diff ',tidal_component,vals_max)

                     # Save and close 
                     plt.savefig(plotname)
                     plt.clf()
 
                     idx_comp=idx_comp+1



# Intepolation of global TPXO9 fields on MED24 grid and vectorial distances  (1/30 --> 1/24)
if vectorial_dist_flag ==  1:
            # Read new grid structure from MED24 mesh_mask files
            # Open the mesh mask file and read lat/lon 
            nc2open3=model_bathy
            model3 = NC.Dataset(nc2open3,'r')
            nc2open4=model_meshmask # mesh mask
            model4 = NC.Dataset(nc2open4,'r')

            lons = model3.variables[longitude_name][:]
            lats = model3.variables[latitude_name][:]
            vals_bathy=model3.variables['Bathymetry'][:]
            vals_land=model4.variables['tmask'][0,0,:,:]
            vals_land=np.squeeze(vals_land)

            x=lons[0] # Array 1D longitudes
            y=[ el[0] for el in lats] #  Array 1D latitudes
            y=np.asarray(y)

            # RMSd table ini
            RMSd_file = open(workdir_path+"RMSd_stats.txt","w")
            RMSd_file.write('Root Mean Square Distances'+'\n')

            # Read Amp and Pha from MED24
            idx_comp = 0
            for idx_2d in range(0,len(field_2d_name)): 
                     var_2d=field_2d_name[idx_2d]
                     var_2d_udm=field_2d_units[idx_2d]

                     P_var_2d=P_field_2d_name[idx_2d]
                     P_var_2d_udm=P_field_2d_units[idx_2d]

                     # Build the path/name of the nc file and open it 
                     nc2open=model_path+model_fileprename+'2D_'+'0_'+'sossheig'+'_'+dates_label+'_'+model_postname+'.nc'
                     print ('Input file = ',nc2open)
                     model_med24 = NC.Dataset(nc2open,'r')

                     # Read field values 
                     vals = model_med24.variables[var_2d][:]*100
                     P_vals = model_med24.variables[P_var_2d][:]*(np.pi/180)

                     thresh = 0.0000
                     mask = np.abs(vals) == thresh
                     vals_ma = np.ma.masked_where(mask, vals)


                     # Choose TPXO9 fields to be used and select the range  
                     print ('Reading Amplitude and Phase fields from TPXO9 model...')
                     # Build the pathname of the nc file and open it 
                     # Loop on tidal components
                     tidal_comp_name_tpxo=['m2','k1','o1','s2','p1','n2','q1','k2']
                     tidal_component=tidal_comp_name_tpxo[idx_comp] 
                     if tidal_component == 'm2':
                        med_range=[0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5] 
                     elif tidal_component == 's2':
                        med_range=[0,0.5,1,1.5,2,2.5,3]
                     else:
                        med_range=[0,0.25,0.5,0.75,1,1.25,1.5]
                     nc2open=tpxo9_path+tpxo9_fileprename+tidal_component+tpxo9_postname+'.nc'
                     print ('Input file = ',nc2open)
                     model = NC.Dataset(nc2open,'r')

                     # Read lat, lon fields only from the first file 
                     if idx_comp == 0:
                        lons_tpxo = model.variables['lon_z'][:]
                        lats_tpxo = model.variables['lat_z'][:]
                        lons_tpxo=np.squeeze(lons_tpxo)
                        lats_tpxo=np.squeeze(lats_tpxo)

                     # Compute Amplitudes and Phases from Re and Im parts             
                     name_var='vals_tpxo_A'+tidal_component
                     globals()[name_var] = (np.sqrt(model.variables['hRe'][:]**2+model.variables['hIm'][:]**2))/10 # Want cm not millimiters!
                     name_pha='vals_tpxo_P'+tidal_component
                     globals()[name_pha] = np.arctan2(-model.variables['hIm'][:],model.variables['hRe'][:])

                     field2interp=globals()[name_var]
                     P_vals_tpxo=globals()[name_pha]
                     
                     # Sea over land Amplitude
                     thresh2interp = 0.0000
                     mask2interp = np.abs(field2interp) == thresh2interp
                     field2interp_ma = np.ma.masked_where(mask2interp,field2interp)
                     input_matrix=field2interp_ma

                     nloop=1
                     infill_value = input_matrix.fill_value
                     output_matrix = ma.copy(input_matrix)  # using ma.copy to prevent future field modifications
                     if np.sum(output_matrix.mask) == 0:  # nothing to fill
                        print('WARNING. Field does not have any land points or wrong field type. Exiting.', file=sys.stderr)
                     else:
                        # iterations loop
                        for loop in range(nloop):
                           # Create a nD x 8 matrix in which, last dimension fixed, the other dimensions
                           #  contains values that are shifted in one of the 8 possible directions
                           # of the last two axes compared to the original matrix
                           shift_matrix = ma.array(np.empty(shape=((8,) + output_matrix.shape)),
                                                   mask=True, fill_value=infill_value, dtype=float)
                           # up shift
                           shift_matrix[0, ..., : -1,:] = output_matrix[..., 1:,:]
                           # down shift
                           shift_matrix[1, ..., 1:, :] = output_matrix[..., : -1, :]
                           # left shift
                           shift_matrix[2, ..., :, : -1] = output_matrix[..., :, 1:]
                           # right shift
                           shift_matrix[3, ..., :, 1:] = output_matrix[..., :, : -1]
                           # up-left shift
                           shift_matrix[4, ..., : -1, : -1] = output_matrix[..., 1:, 1:]
                           # up-right shift
                           shift_matrix[5, ..., : -1, 1:] = output_matrix[..., 1:, : -1]
                           # down-left shift
                           shift_matrix[6, ..., 1:, : -1] = output_matrix[..., : -1, 1:]
                           # down-right shift
                           shift_matrix[7, ..., 1:, 1:] = output_matrix[..., : -1, : -1]
                           # Mediate the shift matrix among the third dimension
                           mean_matrix = ma.mean(shift_matrix, 0)
                           # Replace input masked points with new ones belonging to the mean matrix
                           output_matrix = ma.array(np.where(mean_matrix.mask + output_matrix.mask, mean_matrix, output_matrix),
                                                    mask=mean_matrix.mask, fill_value=infill_value, dtype=float)
                           output_matrix = ma.masked_where(mean_matrix.mask, output_matrix)
                           if np.sum(output_matrix.mask) == 0:  # nothing more to flood
                               print('WARNING. Field does not have anymore land points,', str(loop + 1),
                                     'steps were sufficient to flood it completely.', file=sys.stderr)
                               break
                     field2interp=output_matrix
                    
                     #
                     field2interp=np.reshape(field2interp,len(lons_tpxo)*len(lats_tpxo))
                     P_vals_tpxo=np.reshape(P_vals_tpxo,len(lons_tpxo)*len(lats_tpxo))


                     # Interpolation from the old to the new grid 
                     f_tmp = interpolate.interp2d(lons_tpxo,lats_tpxo,field2interp)
                     p_tmp = interpolate.interp2d(lons_tpxo,lats_tpxo,P_vals_tpxo)
                     var_new = f_tmp(x+360,y)
                     pha_new = p_tmp(x+360,y)

                     # Plot the diffs

                     # Plot the map and save in the path/name
                     plotname=workdir_path+'vecdist_tpxomed24_'+tidal_component+'.jpg'
                     print ('Plot path/name: ',plotname)

                     fig, ax = plt.subplots(1, 3, sharey=True,figsize=(14,4)) 
                     fig.subplots_adjust(wspace=0)
                     plt.rc('font', size=12)
                     # Plot Title
                     plt.suptitle ('Vectorial Distances: '+model_postname+' - TPXO9  --- MED24 grid --- '+tidal_component)

                     # Plot the frame to the map
                     plt.rcParams["axes.linewidth"]  = 1.25
                     # Plot the map and add the amphidromes
                     plt.subplot(1,3,1)
                     plt.ylim(30, 46)
                     plt.xlim(-20, 0)
                     ved=np.sqrt(np.power(np.squeeze(vals)*np.cos(np.squeeze(P_vals))-(np.squeeze(var_new)*np.cos(np.squeeze(pha_new))),2)+np.power(np.squeeze(vals)*np.sin(np.squeeze(P_vals))-(np.squeeze(var_new)*np.sin(np.squeeze(pha_new))),2))
                     cs = plt.contourf(x,y,ved,med_range,cmap='jet',extend='max')
                     contourf1 = plt.contourf(x,y,np.abs(np.squeeze(vals_bathy)),[0.00000,0.00000001], colors='gray') # 0.00000001
                     contour1 = plt.contour(x,y,np.abs(np.squeeze(vals_bathy)),0.00000,colors='black')
                     ax[0].set_xticks(np.arange(340,360,10))
                     plt.tick_params(labelcolor='none', top=False, bottom=True, left=False, right=False)
                     plt.grid(color='black',linestyle='--')
                     max_n=np.max(np.abs(np.abs(np.squeeze(vals))-np.abs(np.squeeze(var_new))))

                     # RMSd
                     rmsd_diffvecd=[]
                     rmsd_count=0
                     sq_vals_bathy=vals_bathy[0,:,:]
                     for lo_idx in range (0,len(x)):
                      if x[lo_idx]>-5 and x[lo_idx]<0:
                         for la_idx in range (0,len(y)):
                           if y[la_idx]>30 and y[la_idx]<40 and sq_vals_bathy[la_idx][lo_idx] != 0:
                                rmsd_diffvecd.append(ved[la_idx][lo_idx])
                                rmsd_count=rmsd_count+1

                     #
                     var_new2 = f_tmp(x,y)
                     pha_new2 = p_tmp(x,y)
                     #
                     plt.subplot(1,3,(2,3))
                     plt.ylim(30, 46)
                     plt.xlim(0, 40)
                     plt.tick_params(labelcolor='none', top=False, bottom=True, left=False, right=False)
                     ved=np.sqrt(np.power(np.squeeze(vals)*np.cos(np.squeeze(P_vals))-(np.squeeze(var_new2)*np.cos(np.squeeze(pha_new2))),2)+np.power(np.squeeze(vals)*np.sin(np.squeeze(P_vals))-(np.squeeze(var_new2)*np.sin(np.squeeze(pha_new2))),2))
                     cs = plt.contourf(x,y,ved,med_range,cmap='jet',extend='max')
                     contourf1 = plt.contourf(x,y,np.abs(np.squeeze(vals_bathy)),[0.00000,0.00000001], colors='gray') # 0.00000001
                     contour1 = plt.contour(x,y,np.abs(np.squeeze(vals_bathy)),0.00000,colors='black')
                     ax[1].get_yaxis().set_visible(False)
                     plt.grid(color='black',linestyle='--')
                     cbar = plt.colorbar(cs,extend='both')
                     cbar.set_label('Vectorial distances [cm]')

                     # RMSd
                     for lo_idx in range (0,len(x)):
                      if x[lo_idx]>0 and x[lo_idx]<40:
                         for la_idx in range (0,len(y)):
                          if y[la_idx]>30 and y[la_idx]<46 and sq_vals_bathy[la_idx][lo_idx] != 0:
                                rmsd_diffvecd.append(ved[la_idx][lo_idx])
                                rmsd_count=rmsd_count+1
                     rmsd_val_ar=np.power(rmsd_diffvecd,2)
                     rmsd_val=np.sum(rmsd_val_ar)
                     rmsd_val=np.sqrt(rmsd_val/(2*rmsd_count))
                     RMSd_file.write(tidal_component+' RMSd_TOT: '+str(rmsd_val)+'\n')

                     # Save and close 
                     plt.savefig(plotname)
                     plt.clf()

                     idx_comp=idx_comp+1



######### BATHY DIFFS
# Intepolation of global TPXO9 fields on MED24 grid and differences  (1/30 --> 1/24)
if bathy_diff_flag ==  1:
            # Read new grid structure from MED24 mesh_mask files
            # Open the mesh mask file and read lat/lon 
            nc2open3=model_bathy
            print ('Input file MED24 = ',nc2open3)
            model3 = NC.Dataset(nc2open3,'r')

            nc2open4=model_meshmask # mesh mask
            model4 = NC.Dataset(nc2open4,'r')

            lons = model3.variables[longitude_name][:]
            lats = model3.variables[latitude_name][:]
            vals_bathy=model3.variables['Bathymetry'][:]

            vals_land=model4.variables['tmask'][0,0,:,:]
            vals_land=np.squeeze(vals_land)

            x=lons[0] # Array 1D longitudes
            y=[ el[0] for el in lats] #  Array 1D latitudes
            y=np.asarray(y)

            vals=vals_bathy
            var_2d_udm='m'

            thresh = 0.0000
            mask = np.abs(vals) == thresh
            vals_ma = np.ma.masked_where(mask, vals)

            rel_bathy_range=[-20,-15,-10,-5,5,10,15,20]
            abs_bathy_range=[-180,-140,-100,-60,-20,20,60,100,140,180]

            # Do the same for TPXO grid files
            nc2open=tpxo9_path+'grid_tpxo9_atlas_30_v2.nc'
            print ('Input file TPXO9 = ',nc2open)
            model = NC.Dataset(nc2open,'r')

            # Read lat, lon fields 
            lons_tpxo = model.variables['lon_z'][:]
            lats_tpxo = model.variables['lat_z'][:]
            lons_tpxo=np.squeeze(lons_tpxo)
            lats_tpxo=np.squeeze(lats_tpxo)

            # Read TPXO bathy 
            field2interp = model.variables['hz'][:]

            # Sea over land
            thresh2interp = 0.0000
            mask2interp = np.abs(field2interp) == thresh2interp
            field2interp_ma = np.ma.masked_where(mask2interp,field2interp)
            input_matrix=field2interp_ma

            nloop=1
            infill_value = input_matrix.fill_value
            output_matrix = ma.copy(input_matrix)  # using ma.copy to prevent future field modifications
            if np.sum(output_matrix.mask) == 0:  # nothing to fill
               print('WARNING. Field does not have any land points or wrong field type. Exiting.', file=sys.stderr)
            else:
               # iterations loop
               for loop in range(nloop):
                  # Create a nD x 8 matrix in which, last dimension fixed, the other dimensions
                  #  contains values that are shifted in one of the 8 possible directions
                  # of the last two axes compared to the original matrix
                  shift_matrix = ma.array(np.empty(shape=((8,) + output_matrix.shape)),
                                          mask=True, fill_value=infill_value, dtype=float)
                  print ('shift_matrix: ',shift_matrix.shape)
                  # up shift
                  shift_matrix[0, ..., : -1,:] = output_matrix[..., 1:,:]
                  # down shift
                  shift_matrix[1, ..., 1:, :] = output_matrix[..., : -1, :]
                  # left shift
                  shift_matrix[2, ..., :, : -1] = output_matrix[..., :, 1:]
                  # right shift
                  shift_matrix[3, ..., :, 1:] = output_matrix[..., :, : -1]
                  # up-left shift
                  shift_matrix[4, ..., : -1, : -1] = output_matrix[..., 1:, 1:]
                  # up-right shift
                  shift_matrix[5, ..., : -1, 1:] = output_matrix[..., 1:, : -1]
                  # down-left shift
                  shift_matrix[6, ..., 1:, : -1] = output_matrix[..., : -1, 1:]
                  # down-right shift
                  shift_matrix[7, ..., 1:, 1:] = output_matrix[..., : -1, : -1]
                  # Mediate the shift matrix among the third dimension
                  mean_matrix = ma.mean(shift_matrix, 0)
                  # Replace input masked points with new ones belonging to the mean matrix
                  output_matrix = ma.array(np.where(mean_matrix.mask + output_matrix.mask, mean_matrix, output_matrix),
                                           mask=mean_matrix.mask, fill_value=infill_value, dtype=float)
                  output_matrix = ma.masked_where(mean_matrix.mask, output_matrix)
                  if np.sum(output_matrix.mask) == 0:  # nothing more to flood
                      print('WARNING. Field does not have anymore land points,', str(loop + 1),
                            'steps were sufficient to flood it completely.', file=sys.stderr)
                      break
            field2interp=output_matrix
            #
            field2interp=np.reshape(field2interp,len(lons_tpxo)*len(lats_tpxo))

            # Interpolation from the old to the new grid 
            f_tmp = interpolate.interp2d(lons_tpxo,lats_tpxo,field2interp)
            var_new = f_tmp(x+360,y)

            # Plot the diffs

            # Plot the absolute diff map and save in the path/name
            plotname=workdir_path+'diff_abs_bathy_tpxomed24.jpg'
            print ('Plot path/name: ',plotname)

            fig, ax = plt.subplots(1, 3, sharey=True,figsize=(14,4)) 
            fig.subplots_adjust(wspace=0)
            plt.rc('font', size=12)
            # Plot Title
            plt.suptitle ('Bathymetry Abs Diff: MED24 - TPXO9  --- MED24 grid ')

            plt.rcParams["axes.linewidth"]  = 1.25
            # Plot the map and add the amphidromes
            plt.subplot(1,3,1)
            plt.ylim(30, 46)
            plt.xlim(-20, 0)
            cs = plt.contourf(x,y,np.squeeze(vals)-np.squeeze(var_new),abs_bathy_range,cmap='bwr',extend='both') 
            contourf1 = plt.contourf(x,y,np.abs(np.squeeze(vals_bathy)),[0.00000,0.00000001], colors='gray') 
            contour2 = plt.contour(x,y,np.squeeze(vals_bathy),[0.00000,0.00000001],colors='black')
            ax[0].set_xticks(np.arange(340,360,10))
            plt.tick_params(labelcolor='none', top=False, bottom=True, left=False, right=False)
            plt.grid(color='black',linestyle='--')
            max_n=np.max(np.abs(np.abs(np.squeeze(vals))-np.abs(np.squeeze(var_new))))
            #
            var_new2 = f_tmp(x,y)
            #
            plt.subplot(1,3,(2,3))
            plt.ylim(30, 46)
            plt.xlim(0, 40)
            plt.tick_params(labelcolor='none', top=False, bottom=True, left=False, right=False)
            cs = plt.contourf(x,y,np.squeeze(vals)-np.squeeze(var_new2),abs_bathy_range,cmap='bwr',extend='both') 
            contourf1 = plt.contourf(x,y,np.abs(np.squeeze(vals_bathy)),[0.00000,0.00000001], colors='gray') 
            contour2 = plt.contour(x,y,np.squeeze(vals_bathy),[0,0.0000001],colors='black')
            ax[1].get_yaxis().set_visible(False)
            plt.grid(color='black',linestyle='--')
            max_p=np.max(np.abs(np.abs(np.squeeze(vals))-np.abs(np.squeeze(var_new2))))
            cbar = plt.colorbar(cs,extend='both')
            cbar.set_label('Bathymetry abs diff [m]')

            # Save and close 
            plt.savefig(plotname)
            plt.clf()

            # Plot the RELATIVE DIFF map and save in the path/name
            plotname=workdir_path+'diff_rel_bathy_tpxomed24.jpg'
            print ('Plot path/name: ',plotname)

            fig, ax = plt.subplots(1, 3, sharey=True,figsize=(14,4))
            fig.subplots_adjust(wspace=0)
            plt.rc('font', size=12)
            # Plot Title
            plt.suptitle ('Bathymetry Relative Diff: (MED24 - TPXO9)*100/(50+MED24)  --- MED24 grid ')

            plt.rcParams["axes.linewidth"]  = 1.25
            # Plot the map and add the amphidromes
            plt.subplot(1,3,1)
            plt.ylim(30, 46)
            plt.xlim(-20, 0)
            cs = plt.contourf(x,y,(np.squeeze(vals)-np.squeeze(var_new))*100.0/(50+np.squeeze(vals)),rel_bathy_range,cmap='bwr',extend='both') # do not consider diffs below 50 m
            contourf1 = plt.contourf(x,y,np.abs(np.squeeze(vals_bathy)),[0.00000,0.00000001], colors='gray')
            contour2 = plt.contour(x,y,np.squeeze(vals_bathy),[0.00000,0.00000001],colors='black')
            ax[0].set_xticks(np.arange(340,360,10))
            plt.tick_params(labelcolor='none', top=False, bottom=True, left=False, right=False)
            plt.grid(color='black',linestyle='--')
            max_n=np.max(np.abs(np.abs(np.squeeze(vals))-np.abs(np.squeeze(var_new))))
            #
            var_new2 = f_tmp(x,y)
            #
            plt.subplot(1,3,(2,3))
            plt.ylim(30, 46)
            plt.xlim(0, 40)
            plt.tick_params(labelcolor='none', top=False, bottom=True, left=False, right=False)
            cs = plt.contourf(x,y,(np.squeeze(vals)-np.squeeze(var_new2))*100.0/(50+np.squeeze(vals)),rel_bathy_range,cmap='bwr',extend='both') # do not consider diffs below 50 m
            contourf1 = plt.contourf(x,y,np.abs(np.squeeze(vals_bathy)),[0.00000,0.00000001], colors='gray')
            contour2 = plt.contour(x,y,np.squeeze(vals_bathy),[0,0.0000001],colors='black')
            ax[1].get_yaxis().set_visible(False)
            plt.grid(color='black',linestyle='--')
            max_p=np.max(np.abs(np.abs(np.squeeze(vals))-np.abs(np.squeeze(var_new2))))
            cbar = plt.colorbar(cs,extend='both')
            cbar.set_label('Bathymetry rel diff [%]')

            # Save and close 
            plt.savefig(plotname)
            plt.clf()

# tidal factors maps [Do-Seong] (1 map per factor)
if doseong_flag == 1:

            print ('# 2D VAR: Do-Seong factors')

            # Build the path/name of the nc file and open it 
            nc2open=model_path+model_fileprename+'2D_'+'0_'+'sossheig'+'_'+dates_label+'_'+model_postname+'.nc'
            print ('Input file = ',nc2open)
            model = NC.Dataset(nc2open,'r')

            # Open the mesh mask file and read lat/lon 
            nc2open3=model_bathy
            print ('Input file MED24 = ',nc2open3)
            model3 = NC.Dataset(nc2open3,'r')
            vals_bathy=model3.variables['Bathymetry'][:]

            vals_land=model4.variables['tmask'][0,0,:,:]
            vals_land=np.squeeze(vals_land)

            # Read lat, lon and fild values 
            lons = model3.variables[longitude_name][:]
            lats = model3.variables[latitude_name][:]

            vals_AM2 = model.variables['M2_Amp'][:]*100 # Want cm not meters!
            vals_AS2 = model.variables['S2_Amp'][:]*100 # Want cm not meters!
            vals_AK1 = model.variables['K1_Amp'][:]*100 # Want cm not meters!
            vals_AO1 = model.variables['O1_Amp'][:]*100 # Want cm not meters!
            vals_AN2 = model.variables['N2_Amp'][:]*100 # Want cm not meters!         

            vals_PM2 = model.variables['M2_Pha'][:]*np.pi/180 # Want rad not deg
            vals_PK1 = model.variables['K1_Pha'][:]*np.pi/180
            vals_PO1 = model.variables['O1_Pha'][:]*np.pi/180

            vals_PM2=np.where(vals_PM2>2*np.pi,vals_PM2-2*np.pi,vals_PM2)
            vals_PK1=np.where(vals_PK1>2*np.pi,vals_PK1-2*np.pi,vals_PK1)
            vals_PO1=np.where(vals_PO1>2*np.pi,vals_PO1-2*np.pi,vals_PO1)

            vals_PM2=np.where(vals_PM2>2*np.pi,vals_PM2-2*np.pi,vals_PM2)
            vals_PK1=np.where(vals_PK1>2*np.pi,vals_PK1-2*np.pi,vals_PK1)
            vals_PO1=np.where(vals_PO1>2*np.pi,vals_PO1-2*np.pi,vals_PO1)

            vals_PM2=np.where(vals_PM2>2*np.pi,vals_PM2-2*np.pi,vals_PM2)
            vals_PK1=np.where(vals_PK1>2*np.pi,vals_PK1-2*np.pi,vals_PK1)
            vals_PO1=np.where(vals_PO1>2*np.pi,vals_PO1-2*np.pi,vals_PO1)

            vals_PM2=np.where(vals_PM2>2*np.pi,vals_PM2-2*np.pi,vals_PM2)
            vals_PK1=np.where(vals_PK1>2*np.pi,vals_PK1-2*np.pi,vals_PK1)
            vals_PO1=np.where(vals_PO1>2*np.pi,vals_PO1-2*np.pi,vals_PO1)


            vals_PM2=np.where(vals_PM2<0,vals_PM2+2*np.pi,vals_PM2)
            vals_PK1=np.where(vals_PK1<0,vals_PK1+2*np.pi,vals_PK1)
            vals_PO1=np.where(vals_PO1<0,vals_PO1+2*np.pi,vals_PO1)

            vals_PM2=np.where(vals_PM2<0,vals_PM2+2*np.pi,vals_PM2)
            vals_PK1=np.where(vals_PK1<0,vals_PK1+2*np.pi,vals_PK1)
            vals_PO1=np.where(vals_PO1<0,vals_PO1+2*np.pi,vals_PO1)

            vals_PM2=np.where(vals_PM2<0,vals_PM2+2*np.pi,vals_PM2)
            vals_PK1=np.where(vals_PK1<0,vals_PK1+2*np.pi,vals_PK1)
            vals_PO1=np.where(vals_PO1<0,vals_PO1+2*np.pi,vals_PO1)

            vals_PM2=np.where(vals_PM2<0,vals_PM2+2*np.pi,vals_PM2)
            vals_PK1=np.where(vals_PK1<0,vals_PK1+2*np.pi,vals_PK1)
            vals_PO1=np.where(vals_PO1<0,vals_PO1+2*np.pi,vals_PO1)

            vals_PM2=np.where(vals_PM2<0,vals_PM2+2*np.pi,vals_PM2)
            vals_PK1=np.where(vals_PK1<0,vals_PK1+2*np.pi,vals_PK1)
            vals_PO1=np.where(vals_PO1<0,vals_PO1+2*np.pi,vals_PO1)

            # Compute factor fields
            F_factor=(np.abs(vals_AK1[:])+np.abs(vals_AO1[:]))/(np.abs(vals_AM2[:])+np.abs(vals_AS2[:]))
            E_factor=(np.abs(vals_AM2[:])+np.abs(vals_AN2[:]))/(np.abs(vals_AM2[:])+np.abs(vals_AS2[:]))
            Ea_factor=np.cos(vals_PK1[:]+vals_PO1[:]-vals_PM2[:])

            thresh = 0.0000
            mask = np.abs(vals_AM2) == thresh
            vals_ma = np.ma.masked_where(mask, vals_AM2)

            # Plot the map and save in the path/name
            for factor_idx in ('F_factor','E_factor','Ea_factor'):
                plotname=workdir_path+model_fileprename+'2D_'+'0_'+factor_idx+'_'+dates_label+'_'+model_postname+'.jpg'
                print ('Plot path/name: ',plotname)
                # Read the field
                field2plot=np.array(globals()[factor_idx])
                thresh = 0.0
                mask_factor = np.abs(field2plot) == thresh
                factor_masked=np.ma.masked_where(mask, field2plot) 

                plt.figure(figsize=(20,10))
                plt.rc('font', size=14)
                # Plot Title
                plt.title (factor_idx+' --- Harmonic Analysis '+dates_label)
                # Read the coordinates for the plot 
                lon_0 = lons.mean() 
                llcrnrlon = lons.min() 
                urcrnrlon = lons.max()  
                lat_0 =  lats.mean()  
                llcrnrlat = lats.min() 
                urcrnrlat = lats.max() 
                # Create the map
                m = Basemap(llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat,urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,resolution='c',projection='merc',lat_0=lat_0,lon_0=lon_0)
                xi, yi = m(lons, lats)
                # Plot the frame to the map
                plt.rcParams["axes.linewidth"]  = 1.25
                m.drawmapboundary(fill_color='gray')

                vals_max=np.amax(field2plot)
                vals_min=np.amin(field2plot)
                upth=vals_max
                lowth=vals_min
                if factor_idx == 'F_factor':
                   upth=4.0
                   lowth=0.0
                elif factor_idx == 'E_factor':
                   upth=1.5
                   lowth=0.0
                elif factor_idx == 'Ea_factor':
                   upth=1.0
                   lowth=-1.0

                # Factors contour
                if factor_idx == 'F_factor':
                   print ('F factor')
                   factor_contour = plt.contourf(xi,yi,np.squeeze(field2plot),[0,0.25,1.5,3.0,3.25],extend='max') 
                elif factor_idx == 'E_factor':
                   print ('E factor')
                   factor_contour = plt.contourf(xi,yi,np.squeeze(field2plot),[0.0,0.8,1.0,1.15,1.5])
                elif factor_idx == 'Ea_factor':
                   print ('Ea factor')
                   factor_contour = plt.contourf(xi,yi,np.squeeze(field2plot),[-1.0,-0.9,-0.1,0.0,0.1,0.9,1.0])
                # Land from bathymetry file
                contourf1 = plt.contourf(xi,yi,np.abs(np.squeeze(vals_bathy)),[0.00000,0.00000001], colors='gray')
                contour1 = plt.contour(xi,yi,np.abs(np.squeeze(vals_bathy)),0.00000,colors='black')
                # Add the grid
                m.drawparallels(np.arange(30., 46., 5.), labels=[1,0,0,0], fontsize=14) 
                m.drawmeridians(np.arange(-20., 40., 10.), labels=[0,0,0,1], fontsize=14) 
                # Plot the legend and its label
                cbar = plt.colorbar(factor_contour,orientation="horizontal",aspect=30)
                
                text_min_x,text_min_y= m(-17,29.0)
                text_max_x,text_max_y= m(32,29.0)
                plt.text(text_max_x,text_max_y,'max='+str(round(vals_max,1)), fontsize=14)
                plt.text(text_min_x,text_min_y,'min='+str(round(vals_min,2)), fontsize=14)
                #

                #
                # Save and close 
                plt.savefig(plotname)
                plt.clf()

# Doseong factors computed from TPXO model 
if doseong_tpxo == 1:

            print ('# 2D VAR: Do-Seong factors from TPXO9 model')
            # Build the path/name of the nc file and open it 
            for tidal_component in ('m2','s2','k1','o1','n2'):
                nc2open=tpxo9_path+tpxo9_fileprename+tidal_component+tpxo9_postname+'.nc'
                print ('Input file = ',nc2open)
                model = NC.Dataset(nc2open,'r')

                # Read lat, lon and fild values 
                lons = model.variables['lon_z'][:]
                lats = model.variables['lat_z'][:]
                lons=np.squeeze(lons)
                lats=np.squeeze(lats)
                name_var='vals_'+tidal_component
                globals()[name_var] = (np.sqrt(model.variables['hRe'][:]**2+model.variables['hIm'][:]**2))/10 # Want cm not millimiters!
                name_pha='vals_P'+tidal_component
                globals()[name_pha] = np.arctan2(-model.variables['hIm'][:],model.variables['hRe'][:]) #*180/np.pi

            vals_Pm2=np.where(vals_Pm2<0,vals_Pm2+2*np.pi,vals_Pm2)
            vals_Pk1=np.where(vals_Pk1<0,vals_Pk1+2*np.pi,vals_Pk1)
            vals_Po1=np.where(vals_Po1<0,vals_Po1+2*np.pi,vals_Po1)

            # Compute factor fields
            F_factor=(np.abs(vals_k1[:])+np.abs(vals_o1[:]))/(np.abs(vals_m2[:])+np.abs(vals_s2[:]))
            E_factor=(np.abs(vals_m2[:])+np.abs(vals_n2[:]))/(np.abs(vals_m2[:])+np.abs(vals_s2[:]))
            Ea_factor=np.cos(vals_Pk1[:]+vals_Po1[:]-vals_Pm2[:])


            thresh = 0.0000
            vals_m2=np.transpose(vals_m2)
            mask = np.abs(vals_m2) == thresh
            vals_ma = np.ma.masked_where(mask, vals_m2)

            # Plot the map and save in the path/name
            for factor_idx in ('F_factor','E_factor','Ea_factor'):
                plotname=workdir_path+model_fileprename+'2D_'+'0_'+factor_idx+'_'+dates_label+tpxo9_postname+'.jpg'
                print ('Plot path/name: ',plotname)
                # Read the field
                field2plot=np.array(globals()[factor_idx])
                field2plot=np.transpose(field2plot)
                thresh = 0.0
                mask_factor = np.abs(field2plot) == thresh
                factor_masked=np.ma.masked_where(mask, field2plot)

                fig, ax = plt.subplots(1, 3, sharey=True,figsize=(14,5))
                fig.subplots_adjust(wspace=0)
                plt.rc('font', size=12)
                # Plot Title
                plt.suptitle (factor_idx+' --- TPXO9 '+dates_label)
                # Read the coordinates for the plot 
                lon_0 = lons.mean()
                llcrnrlon = lons.min()
                urcrnrlon = lons.max()
                lat_0 = lats.mean()
                llcrnrlat = lats.min()
                urcrnrlat = lats.max()
                ## Plot the frame to the map
                plt.rcParams["axes.linewidth"]  = 1.25
                vals_max=np.amax(field2plot)
                vals_min=np.amin(field2plot)
                upth=vals_max
                lowth=vals_min
                if factor_idx == 'F_factor':
                   upth=4.0
                   lowth=0.0
                elif factor_idx == 'E_factor':
                   upth=1.5
                   lowth=0.0
                elif factor_idx == 'Ea_factor':
                   upth=1.0
                   lowth=-1.0
                # Factors contour
                land_value = 0.00000
                if factor_idx == 'F_factor':
                   print ('F factor')
                   plt.subplot(1,3,1)
                   plt.ylim(30, 46)
                   plt.xlim(340, 360)
                   factor_contour1 = plt.contourf(lons,lats,np.squeeze(field2plot),[0,0.25,1.5,3.0,3.25],extend='max') 
                   contourf1 = plt.contourf(lons,lats,np.abs(np.squeeze(vals_m2)),[land_value,0.00000001], colors='gray')
                   ax[0].set_xticks(np.arange(340,360,10))
                   plt.tick_params(labelcolor='none', top=False, bottom=True, left=False, right=False)
                   plt.grid(color='black',linestyle='--')
                   #
                   plt.subplot(1,3,(2,3))
                   plt.ylim(30, 46)
                   plt.xlim(0, 40)
                   plt.tick_params(labelcolor='none', top=False, bottom=True, left=False, right=False)
                   factor_contour2 = plt.contourf(lons,lats,np.squeeze(field2plot),[0,0.25,1.5,3.0,3.25],extend='max') 
                   contourf2 = plt.contourf(lons,lats,np.abs(np.squeeze(vals_m2)),[land_value,0.00000001], colors='gray')
                   contour2 = plt.contour(lons,lats,np.abs(np.squeeze(vals_m2)),land_value,colors='black')
                   ax[1].get_yaxis().set_visible(False)
                   plt.grid(color='black',linestyle='--')
                elif factor_idx == 'E_factor':
                   print ('E factor')
                   plt.subplot(1,3,1)
                   plt.ylim(30, 46)
                   plt.xlim(340, 360)
                   factor_contour1 = plt.contourf(lons,lats,np.squeeze(field2plot),[0.0,0.8,1.0,1.15,1.5])
                   contourf1 = plt.contourf(lons,lats,np.abs(np.squeeze(vals_m2)),[land_value,0.00000001], colors='gray')
                   contour1 = plt.contour(lons,lats,np.abs(np.squeeze(vals_m2)),land_value,colors='black')
                   ax[0].set_xticks(np.arange(340,360,10))
                   plt.tick_params(labelcolor='none', top=False, bottom=True, left=False, right=False)
                   plt.grid(color='black',linestyle='--')
                   #
                   plt.subplot(1,3,(2,3))
                   plt.ylim(30, 46)
                   plt.xlim(0, 40)
                   plt.tick_params(labelcolor='black', top=False, bottom=True, left=False, right=False)
                   factor_contour2 = plt.contourf(lons,lats,np.squeeze(field2plot),[0.0,0.8,1.0,1.15,1.5])
                   contourf2 = plt.contourf(lons,lats,np.abs(np.squeeze(vals_m2)),[land_value,0.00000001], colors='gray')
                   contour2 = plt.contour(lons,lats,np.abs(np.squeeze(vals_m2)),land_value,colors='black')
                   plt.tick_params(labelcolor='none', top=False, bottom=True, left=False, right=False)
                   plt.grid(color='black',linestyle='--')
                elif factor_idx == 'Ea_factor':
                   print ('Ea factor')
                   plt.subplot(1,3,1)
                   plt.ylim(30, 46)
                   plt.xlim(340, 360)
                   factor_contour1 = plt.contourf(lons,lats,np.squeeze(field2plot),[-1.0,-0.9,-0.1,0.0,0.1,0.9,1.0])
                   contourf1 = plt.contourf(lons,lats,np.abs(np.squeeze(vals_m2)),[land_value,0.00000001], colors='gray')
                   contour1 = plt.contour(lons,lats,np.abs(np.squeeze(vals_m2)),land_value,colors='black')
                   ax[0].set_xticks(np.arange(340,360,10))
                   plt.tick_params(labelcolor='none', top=False, bottom=True, left=False, right=False)
                   plt.grid(color='black',linestyle='--')
                   #
                   plt.subplot(1,3,(2,3))
                   plt.ylim(30, 46)
                   plt.xlim(0, 40)
                   plt.tick_params(labelcolor='black', top=False, bottom=True, left=False, right=False)
                   factor_contour2 = plt.contourf(lons,lats,np.squeeze(field2plot),[-1.0,-0.9,-0.1,0.0,0.1,0.9,1.0])
                   contourf2 = plt.contourf(lons,lats,np.abs(np.squeeze(vals_m2)),[land_value,0.00000001], colors='gray')
                   contour2 = plt.contour(lons,lats,np.abs(np.squeeze(vals_m2)),land_value,colors='black')
                   plt.tick_params(labelcolor='none', top=False, bottom=True, left=False, right=False)
                   plt.grid(color='black',linestyle='--')
                #
                # Save and close 
                plt.savefig(plotname)
                plt.clf()

print ('Output path: ',workdir_path)
