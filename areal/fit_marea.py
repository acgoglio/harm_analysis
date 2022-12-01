#
# imports
import sys
import warnings
#
warnings.filterwarnings("ignore") # Avoid warnings
#
import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as NC
import shutil
from scipy.optimize import curve_fit
from scipy import stats
import collections
import pandas as pd
import csv
import math
from datetime import datetime
from operator import itemgetter 
import ttide
#
# by AC Goglio (CMCC)
# annachiara.goglio@cmcc.it
#
# Written: 01/02/2021
#
# Script to fit sossheig filed by means of Foreman methodology
#
#################################################################
# The user should modify the following lines to set his run
#################################################################

# General run parameters:
workdir_path = '/work/oda/med_dev/EAS7/harm_ana/area_2020/'

# Dates
# Choose start and end dates of the period (format dd/mm/yyyy)
inidate = '01/01/2020'
enddate = '30/06/2020'

########################################################
# DO NOT CHANGE THE CODE BELOW THIS LINES
########################################################

# Parameter setting

# MODEL DATASET
# Path and name of inputs datasets
# Currently the extraction of nc is done externally by another script
model_fileprename='map' # DO NOT change this
model_postname='EAS7' # WARNING: Use the same string as in map_extr.ini (ANA_INTAG var)
model_path=workdir_path

# Dates
dates_label=inidate[6:10]+inidate[3:5]+inidate[0:2]+'_'+enddate[6:10]+enddate[3:5]+enddate[0:2]
print ('Whole time interval: ',dates_label)
print ('Model infile template: ',model_path,model_fileprename,'?D_','%vlev%_','%field%',dates_label,'_mod_',model_postname,'.nc')

# Fields to be analyzed
grid = 'T' # Choose T, U, V or uv2t grid

if grid == 'T':
   # 2D fields:
   var_2d='sossheig'
   field_2d_units='m'

# Outfile 
new_path=workdir_path
new_fileprename='amppha' # DO NOT change this
#--------------------------

# Open input and ini files and read infos
# Input file:
nc2open=model_path+model_fileprename+'2D_'+'0_'+var_2d+'_'+dates_label+'_mod_'+model_postname+'.nc'
print ('Input file = ',nc2open)
model = NC.Dataset(nc2open,'r')

# Ini and Out file:
nc2create=new_path+new_fileprename+'2D_'+'0_'+var_2d+'_'+dates_label+'_mod_'+model_postname+'.nc'
ncinifile=nc2create+'_ini'
shutil.copyfile(ncinifile,nc2create)
amppha=NC.Dataset(nc2create,'a') 
print ('Output file = ',nc2create)

# Read lat, lon and field values 
try:
   lat_2d='lat'  
   lon_2d='lon'
   lons = model.variables[lon_2d][:]
   lats = model.variables[lat_2d][:]
   lat_idx='lat' 
   lon_idx='lon' 
   grid_coo_type='new'
except:
   lat_2d='nav_lat'      
   lon_2d='nav_lon'
   lons = model.variables[lon_2d][:]
   lats = model.variables[lat_2d][:]
   lat_idx='y'
   lon_idx='x'
   grid_coo_type='old'

# Read start time and field
time=datetime(int(inidate[6:10]), int(inidate[3:5]), int(inidate[0:2]), 0, 30, 0)
vals = (model.variables[var_2d][:])

#--------------------
# Create the new fields

if grid_coo_type == 'new':
   # Amplitudes:
   M2_Amp=amppha.createVariable('M2_Amp',np.float64,(lat_idx,lon_idx))
   M2_Amp.units = 'm' 
   M2_Amp.standard_name = 'Amplitude'
   #
   K1_Amp=amppha.createVariable('K1_Amp',np.float64,(lat_idx,lon_idx))
   K1_Amp.units = 'm' 
   K1_Amp.standard_name = 'Amplitude'
   #
   O1_Amp=amppha.createVariable('O1_Amp',np.float64,(lat_idx,lon_idx))
   O1_Amp.units = 'm'
   O1_Amp.standard_name = 'Amplitude'
   #
   S2_Amp=amppha.createVariable('S2_Amp',np.float64,(lat_idx,lon_idx))
   S2_Amp.units = 'm'
   S2_Amp.standard_name = 'Amplitude'
   #
   P1_Amp=amppha.createVariable('P1_Amp',np.float64,(lat_idx,lon_idx))
   P1_Amp.units = 'm'
   P1_Amp.standard_name = 'Amplitude'
   #
   N2_Amp=amppha.createVariable('N2_Amp',np.float64,(lat_idx,lon_idx))
   N2_Amp.units = 'm'
   N2_Amp.standard_name = 'Amplitude'
   #
   Q1_Amp=amppha.createVariable('Q1_Amp',np.float64,(lat_idx,lon_idx))
   Q1_Amp.units = 'm'
   Q1_Amp.standard_name = 'Amplitude'
   #
   K2_Amp=amppha.createVariable('K2_Amp',np.float64,(lat_idx,lon_idx))
   K2_Amp.units = 'm'
   K2_Amp.standard_name = 'Amplitude'
   #
   # Phases:
   M2_Pha=amppha.createVariable('M2_Pha',np.float64,(lat_idx,lon_idx))
   M2_Pha.units = 'deg'
   M2_Pha.standard_name = 'Phase'
   #
   K1_Pha=amppha.createVariable('K1_Pha',np.float64,(lat_idx,lon_idx))
   K1_Pha.units = 'deg'
   K1_Pha.standard_name = 'Phase'
   #
   O1_Pha=amppha.createVariable('O1_Pha',np.float64,(lat_idx,lon_idx))
   O1_Pha.units = 'deg'
   O1_Pha.standard_name = 'Phase'
   #
   S2_Pha=amppha.createVariable('S2_Pha',np.float64,(lat_idx,lon_idx))
   S2_Pha.units = 'deg'
   S2_Pha.standard_name = 'Phase'
   #
   P1_Pha=amppha.createVariable('P1_Pha',np.float64,(lat_idx,lon_idx))
   P1_Pha.units = 'deg'
   P1_Pha.standard_name = 'Phase'
   #
   N2_Pha=amppha.createVariable('N2_Pha',np.float64,(lat_idx,lon_idx))
   N2_Pha.units = 'deg'
   N2_Pha.standard_name = 'Phase'
   #
   Q1_Pha=amppha.createVariable('Q1_Pha',np.float64,(lat_idx,lon_idx))
   Q1_Pha.units = 'deg'
   Q1_Pha.standard_name = 'Phase'
   #
   K2_Pha=amppha.createVariable('K2_Pha',np.float64,(lat_idx,lon_idx))
   K2_Pha.units = 'deg'
   K2_Pha.standard_name = 'Phase'
   #
   
   #--------------------
   #Loop on grid points:
   print ('len(lons): ',len(lons),' len(lats): ',len(lats))
   for idx_lon in range (0,len(lons)):
       for idx_lat in range (0,len(lats)): 
         print ('Indexes: ', idx_lon, idx_lat )
   
         # Read the ssh values
         ssh_mod = vals[:,idx_lat,idx_lon]
           
         # Extract Amp and Pha
         xin=np.array(ssh_mod)
         latitudes=lats[idx_lat]
         tidal_comp=['M2','S2','K1','O1','N2','P1','Q1','K2']
   
         # run ttide script
         harmanOUT = ttide.t_tide(np.array(xin), dt=1, stime=time, lat=latitudes, constitnames=tidal_comp, out_style=None, outfile=None)
         tideconout=harmanOUT['tidecon']
   
         # Write values in the arrays
         M2_Amp[idx_lat,idx_lon]=tideconout[5][0]
         S2_Amp[idx_lat,idx_lon]=tideconout[6][0]
         K1_Amp[idx_lat,idx_lon]=tideconout[3][0]
         O1_Amp[idx_lat,idx_lon]=tideconout[1][0]
         N2_Amp[idx_lat,idx_lon]=tideconout[4][0]
         P1_Amp[idx_lat,idx_lon]=tideconout[2][0]
         Q1_Amp[idx_lat,idx_lon]=tideconout[0][0]
         K2_Amp[idx_lat,idx_lon]=tideconout[7][0]
   
         M2_Pha[idx_lat,idx_lon]=tideconout[5][2]
         S2_Pha[idx_lat,idx_lon]=tideconout[6][2]
         K1_Pha[idx_lat,idx_lon]=tideconout[3][2]
         O1_Pha[idx_lat,idx_lon]=tideconout[1][2]
         N2_Pha[idx_lat,idx_lon]=tideconout[4][2]
         P1_Pha[idx_lat,idx_lon]=tideconout[2][2]
         Q1_Pha[idx_lat,idx_lon]=tideconout[0][2]
         K2_Pha[idx_lat,idx_lon]=tideconout[7][2]
   
   amppha.close()

elif grid_coo_type == 'old':

   # Amplitudes:
   M2_Amp=amppha.createVariable('M2_Amp',np.float64,(lat_idx,lon_idx))
   M2_Amp.units = 'm'
   M2_Amp.standard_name = 'Amplitude'
   #
   K1_Amp=amppha.createVariable('K1_Amp',np.float64,(lat_idx,lon_idx))
   K1_Amp.units = 'm'
   K1_Amp.standard_name = 'Amplitude'
   #
   O1_Amp=amppha.createVariable('O1_Amp',np.float64,(lat_idx,lon_idx))
   O1_Amp.units = 'm'
   O1_Amp.standard_name = 'Amplitude'
   #
   S2_Amp=amppha.createVariable('S2_Amp',np.float64,(lat_idx,lon_idx))
   S2_Amp.units = 'm'
   S2_Amp.standard_name = 'Amplitude'
   #
   P1_Amp=amppha.createVariable('P1_Amp',np.float64,(lat_idx,lon_idx))
   P1_Amp.units = 'm'
   P1_Amp.standard_name = 'Amplitude'
   #
   N2_Amp=amppha.createVariable('N2_Amp',np.float64,(lat_idx,lon_idx))
   N2_Amp.units = 'm'
   N2_Amp.standard_name = 'Amplitude'
   #
   Q1_Amp=amppha.createVariable('Q1_Amp',np.float64,(lat_idx,lon_idx))
   Q1_Amp.units = 'm'
   Q1_Amp.standard_name = 'Amplitude'
   #
   K2_Amp=amppha.createVariable('K2_Amp',np.float64,(lat_idx,lon_idx))
   K2_Amp.units = 'm'
   K2_Amp.standard_name = 'Amplitude'
   #
   # Phases:
   M2_Pha=amppha.createVariable('M2_Pha',np.float64,(lat_idx,lon_idx))
   M2_Pha.units = 'deg'
   M2_Pha.standard_name = 'Phase'
   #
   K1_Pha=amppha.createVariable('K1_Pha',np.float64,(lat_idx,lon_idx))
   K1_Pha.units = 'deg'
   K1_Pha.standard_name = 'Phase'
   #
   O1_Pha=amppha.createVariable('O1_Pha',np.float64,(lat_idx,lon_idx))
   O1_Pha.units = 'deg'
   O1_Pha.standard_name = 'Phase'
   #
   S2_Pha=amppha.createVariable('S2_Pha',np.float64,(lat_idx,lon_idx))
   S2_Pha.units = 'deg'
   S2_Pha.standard_name = 'Phase'
   #
   P1_Pha=amppha.createVariable('P1_Pha',np.float64,(lat_idx,lon_idx))
   P1_Pha.units = 'deg'
   P1_Pha.standard_name = 'Phase'
   #
   N2_Pha=amppha.createVariable('N2_Pha',np.float64,(lat_idx,lon_idx))
   N2_Pha.units = 'deg'
   N2_Pha.standard_name = 'Phase'
   #
   Q1_Pha=amppha.createVariable('Q1_Pha',np.float64,(lat_idx,lon_idx))
   Q1_Pha.units = 'deg'
   Q1_Pha.standard_name = 'Phase'
   #
   K2_Pha=amppha.createVariable('K2_Pha',np.float64,(lat_idx,lon_idx))
   K2_Pha.units = 'deg'
   K2_Pha.standard_name = 'Phase'
   #
    
   #--------------------
   #Loop on grid points:
   print ('len(lons): ',len(lons[1]),' len(lats): ',len(lats[:,0]))
   for idx_lon in range (0,len(lons[0])):
       for idx_lat in range (0,len(lats)):
         print ('Indexes: ', idx_lon, idx_lat )
    
         # Read the ssh values
         ssh_mod = vals[:,idx_lat,idx_lon]
    
         # Extract Amp and Pha
         xin=np.array(ssh_mod)
         latitudes=lats[idx_lat][0]
         tidal_comp=['M2','S2','K1','O1','N2','P1','Q1','K2']
    
         # run ttide script
         harmanOUT = ttide.t_tide(np.array(xin), dt=1, stime=time, lat=latitudes, constitnames=tidal_comp, out_style=None, outfile=None)
         tideconout=harmanOUT['tidecon']
   
         # Write values in the arrays
         M2_Amp[idx_lat,idx_lon]=tideconout[5][0]
         S2_Amp[idx_lat,idx_lon]=tideconout[6][0]
         K1_Amp[idx_lat,idx_lon]=tideconout[3][0]
         O1_Amp[idx_lat,idx_lon]=tideconout[1][0]
         N2_Amp[idx_lat,idx_lon]=tideconout[4][0]
         P1_Amp[idx_lat,idx_lon]=tideconout[2][0]
         Q1_Amp[idx_lat,idx_lon]=tideconout[0][0]
         K2_Amp[idx_lat,idx_lon]=tideconout[7][0]
   
         M2_Pha[idx_lat,idx_lon]=tideconout[5][2]
         S2_Pha[idx_lat,idx_lon]=tideconout[6][2]
         K1_Pha[idx_lat,idx_lon]=tideconout[3][2]
         O1_Pha[idx_lat,idx_lon]=tideconout[1][2]
         N2_Pha[idx_lat,idx_lon]=tideconout[4][2]
         P1_Pha[idx_lat,idx_lon]=tideconout[2][2]
         Q1_Pha[idx_lat,idx_lon]=tideconout[0][2]
         K2_Pha[idx_lat,idx_lon]=tideconout[7][2]
   
   amppha.close()
