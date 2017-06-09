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

pathreg = {'onur': {'simroot': '',
                    'emodnet': '?',
                    'cosyna': '/home/onur/WORK/projects/GB/data/stations/stations_COSYNA',
                    'marnet': '/home/onur/WORK/projects/GB/data/stations/stations_MARNET'},
           'ivan': {'simroot': '',
                    'emodnet': '/workm/data/INSITU_NWS_NRT_OBSERVATIONS_013_036/history/mooring/',
                    'cosyna': '?',
                    'marnet': '?',}
           }

def main():
    #PARAMETERS:
    # common
    user='onur' #ivan,onur
    rootpath='./'
    depthints={'surface':[0,5]} # 'bottom':[5,0] #for bottom, depthint is relative to bottom depth
    #timeint = [datetime.datetime(2010, 1, 1,0,0,0), datetime.datetime(2010, 12, 31,23,59,59)]
    timeint = [datetime.datetime(2011, 1, 1, 0, 0, 0), datetime.datetime(2014, 12, 31, 23, 59, 59)]
    # regarding observations
    statsets = ['cosyna'] #,'marnet'
    stations = []
    # regarding simulations
    simpaths = {'sim1': 'some_path'}
    simnames = {'sim1': 'some_run_with_some_model'}
    # regarding plotting parameters
    sims2plot = ['sim1']
    # derived:
    pickledobsfile=os.path.join(rootpath,'obs_%s_%s_%s-%s.pickle'%('-'.join(statsets),'-'.join(depthints.keys()),timeint[0].year,timeint[1].year))
    pickledsimfile='bla'

    #READ OBSERVATIONS
    obs=SRO.readobs(pathreg[user],statsets,pickledobsfile,stations,timeint,depthints)
    print(obs)
    print(len(obs))
    return

    #READ SIMULATIONS
    simset={}    
    for simno, modid in enumerate(sims2plot):
        print ('reading simulation:%s'%modid)
        #sim=stations_readobs(simpaths[modid],timeint,obs)
        sim=0        
        simset[modid]=sim
    
    #PLOTS
    print ('doing the ts plots')    
    #stations_tsplots(obs,simset,depths)
    #stations_scatterplots()

if __name__=='__main__':
    main()
    