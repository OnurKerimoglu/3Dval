# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 16:06:57 2017

@author: kuznetso
"""
import os

def readobs(timeint,pickledobsfile,depthints,stations=[],):
    print ('reading ovservations:')
    
    #check if the pickledobs exist
    #return obs
    
    #ifnotexist:    
    if len(stations)==0:
        stats2read=['s1']
    
    obs={}
    for skey in stats2read:
        sdata = fill_stationdata(skey,depthints)
        obs[skey] = sdata
    
    #pickle the obs file
    
    return obs

def fill_stationdata(skey,depthints):
    tempdata={}; saltdata={}; sshdata={}    
    for layername,depthint in depthints.items():    
        tempdata[layername]={'time':[], 'value':[], 'depth_interval':depthint}
        saltdata[layername]={'time':[], 'value':[], 'depth_interval':depthint}
    sshdata['z0']={'time':[], 'value':[], 'depth_interval':depthint}
    sdata ={'longname':'descriptive name of the station',
            'lon':0.0,
            'lat':0.0,
            'max_depth':0.0,
            'temp' : tempdata,
            'salt' : saltdata,
            'ssh': sshdata,
            }
    return sdata