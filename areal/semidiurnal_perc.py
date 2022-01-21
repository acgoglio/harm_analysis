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
# Written: 10/12/2021
# Modified: 10/12/2021
#
# HOW TO BUILD THE INPUT FILE:   
# cdo expr,'ST_Amp=M2_Amp+S2_Amp+N2_Amp+K2_Amp;TOT_Amp=M2_Amp+S2_Amp+N2_Amp+K2_Amp+K1_Amp+O1_Amp+P1_Amp+Q1_Amp;ST_perc=100.0*ST_Amp/TOT_Amp' input.nc Amp_perc_tmp.nc
# ncatted -a units,'ST_perc',c,c,'%' /work/oda/ag15419/tmp/tidal_bdy_HA/a_drpo_2020/Amp_perc_tmp.nc
# cdo setmisstoc,0 /work/oda/ag15419/tmp/tidal_bdy_HA/a_drpo_2020/Amp_perc_tmp.nc /work/oda/ag15419/tmp/tidal_bdy_HA/a_drpo_2020/Amp_perc.nc
# 
# TPXO:
# for TOL in m2 s2 k1 o1 n2 p1 q1 k2 ; do cdo expr,"Amp_${TOL}=sqrt(hRe*hRe+hIm*hIm);lon_z=lon_z*1;lat_z=lat_z*1" /work/oda/ag15419/tmp/HA_twd/tpxo_semidiurnal/h_${TOL}_tpxo9_atlas_30_v2.nc /work/oda/ag15419/tmp/HA_twd/tpxo_semidiurnal/A_${TOL}_tpxo9_atlas_30_v2.nc ; done
# cdo merge /work/oda/ag15419/tmp/HA_twd/tpxo_semidiurnal/A_*_tpxo9_atlas_30_v2.nc /work/oda/ag15419/tmp/HA_twd/tpxo_semidiurnal/A_all_tpxo9_atlas_30_v2.nc
# cdo expr,'ST_Amp=Amp_m2+Amp_s2+Amp_n2+Amp_k2;TOT_Amp=Amp_m2+Amp_s2+Amp_n2+Amp_k2+Amp_k1+Amp_o1+Amp_p1+Amp_q1;ST_perc=100.0*ST_Amp/TOT_Amp;lon_z=lon_z*1;lat_z=lat_z*1' /work/oda/ag15419/tmp/HA_twd/tpxo_semidiurnal/A_all_tpxo9_atlas_30_v2.nc /work/oda/ag15419/tmp/HA_twd/tpxo_semidiurnal/tpxo_semid_perc.nc
#
#################################################################
# The user should modify the following lines to set his run
#################################################################
# General run parameters: 
#---------------------
# work dir path and bathymetry path/name
workdir_path = '/work/oda/ag15419/tmp/HA_twd/a_twd/'
model_bathy='/work/oda/ag15419/PHYSW24_DATA/TIDES/DATA0/bathy_meter.nc'
model_meshmask='/work/oda/ag15419/PHYSW24_DATA/TIDES/DATA0/mesh_mask.nc'
outplot=workdir_path+'/AmpST_perc.png'

flag_tpxo=1
workdir_path2 = '/work/oda/ag15419/tmp/HA_twd/tpxo_semidiurnal/'
nc2open_tpxo=workdir_path2+'/tpxo_semid_perc.nc'
outplot_tpxo=workdir_path2+'/tpxo_AmpST_perc.png'
outplot_glo=workdir_path2+'/tpxo_AmpST_perc_glo.png'

# Build the path/name of the nc file and open it 
nc2open=workdir_path+'/Amp_perc_tmp.nc'
print ('Input file = ',nc2open)
model = NC.Dataset(nc2open,'r')
ST_perc=model.variables['ST_perc'][:]
vals=ST_perc

# Read lat, lon and fild values 
nc2open3=model_bathy # tidal bathimetry
model3 = NC.Dataset(nc2open3,'r')
nc2open4=model_meshmask # mesh mask
model4 = NC.Dataset(nc2open4,'r')

vals_bathy=model3.variables['Bathymetry'][:]
lons = model3.variables['nav_lon'][:]
lats = model3.variables['nav_lat'][:]

vals_land=model4.variables['tmask'][0,0,:,:]
vals_land=np.squeeze(vals_land)

vals_ma = np.ma.masked_where(vals_land, vals)

plt.figure(figsize=(20,10))
plt.rc('font', size=16)

plt.title ('Semidiurnal Component Amplitude Percentages')

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

# Amp plot:
cs = plt.contourf(xi,yi,vals,levels=[0,10,20,30,40,50,60,70,80,90,100],cmap='coolwarm')
contour50 = plt.contour(xi,yi,vals,[50],colors='blue')

# Add the grid
m.drawparallels(np.arange(30., 46., 5.), labels=[1,0,0,0], fontsize=10)
m.drawmeridians(np.arange(-20., 40., 10.), labels=[0,0,0,1], fontsize=10)

# Land from bathymetry or mesh mask file
contourf1 = plt.contourf(xi,yi,np.abs(np.squeeze(vals_bathy)),[0.00000,0.00000001], colors='gray')
contour1 = plt.contour(xi,yi,np.abs(np.squeeze(vals_bathy)),[0.00000],colors='black')
#contourf1 = plt.contourf(xi,yi,np.abs(np.squeeze(vals_bathy)),[0.00000,500.0], colors='gray')
#contour1 = plt.contour(xi,yi,np.abs(np.squeeze(vals_bathy)),[500.0],colors='black')

# Plot the legend 
cbar = m.colorbar(cs, location='bottom', pad="10%")
bar_label_string=' Semidiurnal component [%]'
cbar.set_label(bar_label_string)

plt.savefig(outplot)
print ('Outfile ',outplot)
plt.clf()

#####################################
# TPXO
#####################################
# TPXO Phase maps
if flag_tpxo == 1:
   print ('Maps from TPXO9 model')
   print ('Input file = ',nc2open_tpxo)

   model = NC.Dataset(nc2open_tpxo,'r')

   # Read lat, lon and fild values 
   lons_tpxo = model.variables['lon_z'][:]
   lats_tpxo = model.variables['lat_z'][:]
   lons_tpxo=np.squeeze(lons_tpxo)
   lats_tpxo=np.squeeze(lats_tpxo)
   #
   field2plot=model.variables['ST_perc'][:]
   field2plot=np.transpose(field2plot)
   #
   land2plot=model.variables['TOT_Amp'][:]
   land2plot=np.transpose(land2plot)

   # Set the land 
   land_value = 0.00000
   thresh = -2.147484e+09
   mask = np.abs(field2plot) == thresh
   vals_tpxo_ma = np.ma.masked_where(mask, field2plot)
   current_cmap = plt.cm.get_cmap('coolwarm')
   current_cmap.set_bad(color='gray')

   # Plot the map and save in the path/name
   plotname=outplot_tpxo
   print ('Plot path/name: ',plotname)

   fig, ax = plt.subplots(1, 3, sharey=True,figsize=(14,4))
   fig.subplots_adjust(wspace=0)
   plt.rc('font', size=12)
   # Plot Title
   plt.suptitle (' TPXO9 Semidiurnal Component Amplitude Percentages')
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
   plt.ylim((30.188,45.979))
   plt.xlim(340,360) #(360-18.125, 360)
   factor_contour1 = plt.contourf(lons_tpxo,lats_tpxo,field2plot,levels=[0,10,20,30,40,50,60,70,80,90,100],cmap=current_cmap)
   contour50 = plt.contour(lons_tpxo,lats_tpxo,vals_tpxo_ma,[50],colors='blue')
   contourf_land = plt.contourf(lons_tpxo,lats_tpxo,np.abs(np.squeeze(land2plot)),[0.00000,0.00000001], colors='gray')
   contour_land = plt.contour(lons_tpxo,lats_tpxo,np.abs(np.squeeze(land2plot)),0.00000,colors='black')

   ax[0].set_xticks(np.arange(340,360,10))
   plt.tick_params(labelcolor='none', top=False, bottom=True, left=False, right=False)
   plt.grid(color='black',linestyle='--')
   #
   plt.subplot(1,3,(2,3))
   plt.ylim(30.188,45.979)
   plt.xlim(0,40) #(0,36.292 )
   plt.tick_params(labelcolor='none', top=False, bottom=True, left=False, right=False)
   factor_contour2 = plt.contourf(lons_tpxo,lats_tpxo,field2plot,levels=[0,10,20,30,40,50,60,70,80,90,100],cmap=current_cmap)
   contour50 = plt.contour(lons_tpxo,lats_tpxo,vals_tpxo_ma,[50],colors='blue')
   contourf_land2 = plt.contourf(lons_tpxo,lats_tpxo,np.abs(np.squeeze(land2plot)),[0.00000,0.00000001], colors='gray')
   contour_land2 = plt.contour(lons_tpxo,lats_tpxo,np.abs(np.squeeze(land2plot)),0.00000,colors='black')

   ax[1].get_yaxis().set_visible(False)
   plt.grid(color='black',linestyle='--')
   cbar = plt.colorbar(factor_contour2, extend='max',label=' Semidiurnal component [%]')

   # Save and close 
   plt.savefig(plotname)
   plt.clf()

   # GLOBAL MAP
   # Plot the map and save in the path/name
   plotname=outplot_glo
   print ('Plot path/name: ',plotname)
   plt.figure(figsize=(20,10))
   #fig, ax = plt.subplots(1, 3, sharey=True,figsize=(14,4))
   fig.subplots_adjust(wspace=0)
   plt.rc('font', size=12)
   # Plot Title
   plt.suptitle (' TPXO9 Semidiurnal Component Amplitude Percentages')
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

   plt.ylim((llcrnrlat,urcrnrlat))
   plt.xlim(llcrnrlon,urcrnrlon)
   factor_contour1 = plt.contourf(lons_tpxo,lats_tpxo,vals_tpxo_ma,levels=[0,10,20,30,40,50,60,70,80,90,100],cmap='coolwarm')
   contour50 = plt.contour(lons_tpxo,lats_tpxo,vals_tpxo_ma,[50],colors='blue')
   landcf=plt.contourf(lons_tpxo,lats_tpxo,field2plot,[-1e+10,-2.147484e+09,0],colors='gray')
   landc=plt.contour(lons_tpxo,lats_tpxo,field2plot,[-2.147484e+09],colors='black')

   plt.grid(color='black',linestyle='--')
   #
   plt.grid(color='black',linestyle='--')
   cbar = plt.colorbar(factor_contour2, extend='max',label=' Semidiurnal component [%]')

   # Save and close 
   plt.savefig(plotname)
   plt.clf()
