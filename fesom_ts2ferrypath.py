# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 10:04:47 2017

@author: kuznetso
"""
# before applying, cut input file, leave only surface and variable of interest
# cdo daymean ts.nc daymean_ts.nc
# ncks -d x_2,0 -v salinity,lon,lat,lon_elem,lat_elem,nv ts.nc ts_1.nc

from netCDF4 import Dataset
import datetime
import os
import time 
from scipy.interpolate import griddata
import numpy as np
from netCDF4 import num2date

#input file read lon,lat for interpolation
fexample = Dataset("sns144-M161117n-P161118-bdyi3-mergedextract_phys_FB_yoana-buhel_salt_2012-2013.nc")
path_lon = fexample.variables['longitude'][:]
path_lat = fexample.variables['latitude'][:]
fexample.close()

#model input
fname_model = 'ts_ferry.nc'

#output files
fname_out = 'fesom.nc'

# open model, read model time,lon,lat
fin = Dataset(fname_model,'r', format="NETCDF4")
mtime_raw = fin.variables['time'][:]
a = fin.variables['time'].getncattr('units')
mtime = num2date(mtime_raw,a)  
mlon = fin.variables['lon'][:]
mlat = fin.variables['lat'][:]

#create output file, define dim and var
if os.path.exists(fname_out):
    os.remove(fname_out)
nc = Dataset(fname_out,'w', format="NETCDF4")
nc.description = 'File containing data from FESOM-coastal model interpolated on ferryline path'
nc.history     = 'Created ' + time.ctime(time.time())
nc.source      = fname_model
# define time dimension and variables
dim_time  = nc.createDimension('time', None)
dim_lon   = nc.createDimension('longitude', len(path_lon))
var_time  = nc.createVariable('time', 'd',('time'))
var_lon   = nc.createVariable('longitude','d',('longitude'))
var_lat   = nc.createVariable('latitude', 'd',('longitude'))
var_data  = nc.createVariable('data', 'd',('time','longitude'))
var_time.units    = 'seconds since 2000-01-01T00:00:00Z'
var_data.units  = 'PSU'
#put lon,lat,time to output  
var_lon[:]   = path_lon
var_lat[:]   = path_lat
times = np.zeros(len(mtime))
for j,item in enumerate(mtime):
   times[j] = item.timestamp() - datetime.datetime(2000,1,1).timestamp()
var_time[:] = times;

mod_data = np.zeros([len(mtime),len(path_lon)])

#loop over time steps in infile (read every time step in loop, file could be very big)
for it in range(len(mtime)):
    print(it)
    mdata = fin.variables['salinity'][it,:,:]
    d = griddata((mlon,mlat),mdata,(path_lon,path_lat), method = 'linear')
    mod_data[it,:] = d.reshape(1,len(d))
    
    
var_data[:] = mod_data;    

nc.close()    
fin.close()    
    
    