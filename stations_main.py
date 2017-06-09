# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 15:02:46 2017

Main file for stations validation.
Control of structre, call reading obs. model and plot.

@author: 

"""
import os,sys
import datetime
import stations_readobs as SRO

#to reload modules
import importlib
importlib.reload(SRO)

def main():
    
    rootpath='./'
    simpaths={'sim1':'some_path'}
    simnames={'sim1':'some_run_with_some_model'}
    stats2include=['s1','s2']
    sims2plot=['sim1']
    depthints={'surface':[0,5]}
    pickledobsfile=os.path.join(rootpath,'obs_all.pickle')
    pickledsimfile='bla'
    # path to observations data sets (now only EMODNET)    
    path2data_emodnet   = '/workm/data/INSITU_NWS_NRT_OBSERVATIONS_013_036/history/mooring/'

    #some common parameters
    timeint=[datetime.date(2011,1,1),datetime.date(2013,12,31)]
    
    #read the observations
    obs=SRO.readobs(timeint,pickledobsfile,depthints,path2data_emodnet)
    print(obs)
    print(len(obs))
    #read the simulations
    simset={}    
    for simno, modid in enumerate(sims2plot):
        print ('reading simulation:%s'%modid)
        #sim=stations_readobs(simpaths[modid],timeint,obs)
        sim=0        
        simset[modid]=sim
    
    #do the time series plots
    print ('doing the ts plots')    
    #stations_tsplots(obs,simset,depths)
    #stations_scatterplots()

if __name__=='__main__':
    main()
    