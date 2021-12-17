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
#################################################################
# The user should modify the following lines to set his run
#################################################################
# General run parameters: 
#---------------------
# work dir path and bathymetry path/name
workdir_path = '/work/oda/ag15419/tmp/HA_twd/a_twd_bf/'
model_bathy='/work/oda/ag15419/PHYSW24_DATA/TIDES/DATA0/bathy_meter.nc'
model_meshmask='/work/oda/ag15419/PHYSW24_DATA/TIDES/DATA0/mesh_mask.nc'
outplot='/work/oda/ag15419/tmp/tidal_bdy_HA/a_drpo_2020/AmpST_perc.png'

# Build the path/name of the nc file and open it 
nc2open='/work/oda/ag15419/tmp/tidal_bdy_HA/a_drpo_2020/Amp_perc.nc'
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
#contourf1 = plt.contourf(xi,yi,np.abs(np.squeeze(vals_bathy)),[0.00000,0.00000001], colors='gray')
#contour1 = plt.contour(xi,yi,np.abs(np.squeeze(vals_bathy)),[0.00000],colors='black')
contourf1 = plt.contourf(xi,yi,np.abs(np.squeeze(vals_bathy)),[0.00000,500.0], colors='gray')
contour1 = plt.contour(xi,yi,np.abs(np.squeeze(vals_bathy)),[500.0],colors='black')

# Plot the legend 
cbar = m.colorbar(cs, location='bottom', pad="10%")
bar_label_string=' Semidiurnal component [%]'
cbar.set_label(bar_label_string)

plt.savefig(outplot)
print ('Outfile ',outplot)
plt.clf()

