# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 16:06:57 2017

@author: kuznetso
"""
import os
import pickle
import numpy as np
from netCDF4 import num2date
from   netCDF4 import Dataset as open_ncfile

def readobs(timeint,pickledobsfile,depthints,path2data_emodnet,stations=[],):
    print ('reading ovservations:')
    
    # paths 2 datasets
    #path2data_emodnet   = '/workm/data/sns_ts_emodnet/'

    #check if the pickledobs exist
    #return obs
    #if os.path.exists(pickledobsfile):
    #    print('opening pickled obs file')
    #    (obs,) = np.load(pickledobsfile)
    #    return obs
    print(path2data_emodnet)    
    #ifnotexist:    
    if len(stations)==0:
        # get lists of availeble data 
        #EMODNET
        emodnet_files = [f for f in os.listdir(path2data_emodnet) if f.endswith('.nc')]      
        stats2read = emodnet_files[:]
        
    obs={}
    for skey in stats2read:
        print(skey)
        sdata = fill_stationdata(skey,depthints,path2data_emodnet)
        obs[skey] = sdata
    
    #pickle the obs file
    #f=open(pickledobsfile,'wb')
    #pickle.dump((obs,),f) #,protocol=-1
    #f.close()
    
    return obs

def fill_stationdata(skey,depthints,path2data_emodnet):
    tempdata={}; saltdata={}; sshdata={}
    # open net cdf file from emodnet    
    ncf = open_ncfile(path2data_emodnet+skey,'r')
    
    # reading metadata of station
    lon = ncf.variables['LONGITUDE'][:][0]
    lat = ncf.variables['LATITUDE'][:][0]
    depth = ncf.variables['DEPH'][:][0]
    tempfound = False; saltfound = False; sshfound = False
    #check if dates are in 
    #check if variables are in file
    if 'TEMP' in ncf.variables: tempfound = True
    if 'PSAL' in ncf.variables: saltfound = True
    if 'SLEV' in ncf.variables: sshfound  = True
        
    #a = ncf.variables['TIME'].getncattr('units')
    #time_emodnet.append([num2date(ncf.variables['TIME'][:][0],a),num2date(ncf.variables['TIME'][:][-1],a)])
    #ncf.close()
    
 #          if var in var_emodnet[i][:]: 
 #              if (time_emodnet[i][0] < t2 or time_emodnet[i][1] > t2):
 #                  # get data from file (first match)
 #                  ncf = open_ncfile(path2data_emodnet+item,'r')
#                   obs_z    = ncf.variables['DEPH'][:][0]                   
#                   obs_time = ncf.variables['TIME'][:]
#                   a = ncf.variables['TIME'].getncattr('units')
#                   obs_time = num2date(obs_time,a)
#                   if np.any(obs_time > t1) and np.any(obs_time < t2):
#                       time_ind1 = np.where(obs_time > t1)[0][0]
#                       time_ind2 = np.where(obs_time < t2)[0][-1]
#                   else:
#                       continue
#                   data = ncf.variables[v2v[var]][:][time_ind1:time_ind2,:]
#                   ncf.close()
#                   ind_time = obs_time[time_ind1:time_ind2]
#                   for iz,z in enumerate(obs_z):
#                       #check if data array is masked, if not make new masked array
#                       if isinstance(data[:,iz],np.ma.MaskedArray):
#                           dd = data[:,iz]
#                       else:
#                           dd = np.ma.array(data[:,iz])
#                       if dd.count():
#                           dset={'fname':item,'time':ind_time,'data':dd,'z':z}
#                           dset_list.append(dset)
    # get time    
    obs_time = ncf.variables['TIME'][:]                       
    a = ncf.variables['TIME'].getncattr('units')
    obs_time = num2date(obs_time,a)
    #handle temp
    if tempfound:
        tempdata['presence']=True
        for iz, z  in enumerate(depth):
            data = ncf.variables['TEMP'][:][:,iz]    
            tempdata[z]={'time':obs_time, 'value':data, 'depth_interval':z}
    else:
        tempdata['presence']=False
    
    #handle salt
    if saltfound:
        saltdata['presence']=True
        for iz, z in enumerate(depth):
            data = ncf.variables['PSAL'][:][:,iz]   
            saltdata[z]={'time':obs_time, 'value':data, 'depth_interval':z}
    else:
        saltdata['presence']=False
    
    #handle ssh 
    if sshfound:
        sshdata['presence']=True
        for iz, z in enumerate(depth):
            data = ncf.variables['SLEV'][:][:,iz]   
            sshdata[z]={'time':obs_time, 'value':data, 'depth_interval':z}
    else:
        sshdata['presence']=False
    
    #put all data in
    sdata ={'longname':'descriptive name of the station',
            'lon':lon,
            'lat':lat,
            'max_depth':0.0,
            'temp' : tempdata,
            'salt' : saltdata,
            'ssh': sshdata
            }
            
    ncf.close()
            
    return sdata