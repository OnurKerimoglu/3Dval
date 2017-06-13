# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 15:02:46 2017

Main file for stations validation.
Control of structre, call reading obs. model and plot.

@author: kerimoglu.o@gmail.com, ivan.kuznetsov@gmail.com

"""
import os,sys
import datetime
import stations_readobs as SRO
import stations_readsim as SRS
import stations_plots as SP

#to reload modules
import importlib
importlib.reload(SRO)
importlib.reload(SRS)
importlib.reload(SP)

pathreg = {'onur': {'GETM-SNS-new': '/home/onur/WORK/projects/GB/maecs/3d/sns144-M161117n-P161118-bdyi3-z01mm-wAtmN/sns144-M161117n-P161118-bdyi3-z01mm-wAtmN-mergedextract_phys_2006-2010_zSB.nc',
                    'GETM-SNS-old': '/home/onur/WORK/projects/GB/maecs/3d/sns144-M161117-P161118/sns144-M161117-P161118-mergedextract_phys_zSB_2006-2010.nc',
                    'plotrootpath':'/home/onur/WORK/projects/GB/maecs/3d/sns144-M161117n-P161118-bdyi3-z01mm-wAtmN',
                    'rootpath': './',
                    'emodnet': '?',
                    'cosyna': '/home/onur/WORK/projects/GB/data/stations/stations_COSYNA',
                    'marnet': '/home/onur/WORK/projects/GB/data/stations/stations_MARNET'},
           'ivan': {'rootpath': './',
                    'emodnet': '/workm/data/INSITU_NWS_NRT_OBSERVATIONS_013_036/history/mooring/',
                    'cosyna': '?',
                    'marnet': '?',}
           }

def main():
    #PARAMETERS:
    # general
    user='onur' #ivan,onur
    depthints={'surface':[0,5]} #,'bottom':[5,0]} #for bottom, depthint is relative to bottom depth
    timeint = [datetime.datetime(2006, 1, 1,0,0,0), datetime.datetime(2010, 12, 31,23,59,59)]
    ##timeint = [datetime.datetime(2000, 1, 1, 0, 0, 0), datetime.datetime(2010, 12, 31, 23, 59, 59)]
    # regarding observations. todo: place in a dictionary, specify variables to extract
    statsets = ['cosyna'] #,'marnet','cosyna'
    stations = []
    readobsraw=False #i.e., if the pickle file should be ignored
    # regarding simulations. todo: place in a dictionary, specify variables to extract
    sims2plot = ['GETM-SNS-new','GETM-SNS-old']
    readsimraw=False #i.e., if the pickle file should be ignored
    simdomain=''
    meth2D='int_tree'
    #regarding plots:
    plotopts={'TS':True,'TSstyle':'TSdefault','varns':['temp','salt']}

    #READ OBSERVATIONS
    obs=SRO.readobs(pathreg[user],readobsraw,statsets,stations,timeint,depthints)

    #READ SIMULATIONS
    simset={}    
    for simno, simname in enumerate(sims2plot):
        sim=SRS.readsim(pathreg[user],simname,readsimraw,simdomain,meth2D,statsets,timeint,depthints,obs)
        simset[simname]=sim

    #print(simset)

    #PLOTS
    SP.stations_plots(plotopts, obs, simset, pathreg[user]['plotrootpath'], stations, timeint, depthints)

if __name__=='__main__':
    main()