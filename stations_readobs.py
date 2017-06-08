# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 16:06:57 2017

@author: kuznetso
"""
import os
import pickle
import numpy as np

def readobs(timeint,pickledobsfile,depthints,stations=[],):
    print ('reading ovservations:')
    
    #check if the pickledobs exist
    #return obs
    #if os.path.exists(pickledobsfile):
    #    print('opening pickled obs file')
    #    (obs,) = np.load(pickledobsfile)
    #    return obs
        
    #ifnotexist:    
    if len(stations)==0:
        stats2read=['s1']
    
    obs={}
    for skey in stats2read:
        sdata = fill_stationdata(skey,depthints)
        obs[skey] = sdata
    
    #pickle the obs file
    #f=open(pickledobsfile,'wb')
    #pickle.dump((obs,),f) #,protocol=-1
    #f.close()
    
    return obs

def fill_stationdata(skey,depthints):
    tempdata={}; saltdata={}; sshdata={}
    
    #this info is coming from a function probably
    tempfound=True
    saltfound=True
    sshfound=False
    
    #handle temp
    if tempfound:
        tempdata['presence']=True
        for layername,depthint in depthints.items():    
            tempdata[layername]={'time':[], 'value':[], 'depth_interval':depthint}
    else:
        tempdata['presence']=False
    
    #handle salt
    if saltfound:
        saltdata['presence']=True
        for layername,depthint in depthints.items():    
            saltdata[layername]={'time':[], 'value':[], 'depth_interval':depthint}
    else:
        saltdata['presence']=False
    
    #handle ssh 
    if sshfound:
        sshdata['presence']=True
        sshdata['z1']={'time':[], 'value':[], 'depth_interval':depthint}
    else:
        sshdata['presence']=False
    
    #put all data in
    sdata ={'longname':'descriptive name of the station',
            'lon':0.0,
            'lat':0.0,
            'max_depth':0.0,
            'temp' : tempdata,
            'salt' : saltdata,
            'ssh': sshdata
            }
    return sdata