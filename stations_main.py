# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 15:02:46 2017

Main file for stations validation.
Control of structre, call reading obs. model and plot.

@author: kerimoglu.o@gmail.com

"""
import os,sys
import warnings
import datetime
import stations_readobs as SRO
import stations_readsim as SRS
import stations_plots as SP

#to reload modules
import importlib
importlib.reload(SRO)
importlib.reload(SRS)
importlib.reload(SP)

pathreg = {'onur': {#'GF': '/home/onur/WORK/projects/2013/maecs/sns144-M180109-nBF-Pbg2017-B180106-vsdetp4b1/sns144-M180109-nBF-Pbg2017-B180106-vsdetp4b1-1113_phys_zSB.nc',
                    #'plotrootpath':'/home/onur/WORK/projects/2013/maecs/sns144-M180106-nBF-Pbg2017-B180106-vsdetp4b1/','
                    #'GF': '/home/onur/WORK/projects/2013/maecs/sns144-M180109-nFpB-Pbg2017-B180106-vsdetp4b1/sns144-M180109-nFpB-Pbg2017-B180106-vsdetp4b1-1114_phys_zSB.nc',
                    'GF': '/home/onur/WORK/projects/2013/maecs/sns144-M180109-nFpB-Pbg2017-B180106-vsdetp4b1/sns144-M180109-nFpB-Pbg2017-B180106-vsdetp4b1-1114_chl_zSB.nc',
                    'plotrootpath':'/home/onur/WORK/projects/2013/maecs/sns144-M180109-nFpB-Pbg2017-B180106-vsdetp4b1/',
                    'pickledobspath': './',
                    'BGC':    '/home/onur/WORK/projects/GB/data/stations/individual/BGC/',
                    'BSH':    '/home/onur/WORK/projects/GB/data/stations/individual/BSH/',
                    'cosyna': '/home/onur/WORK/projects/GB/data/stations/COSYNA/proc/nc'},
           'newuser': {}
           }

def main():
    #PARAMETERS:
    # general
    user='onur'
    #vars=['temp','salt','DOs']
    vars = ['DIN','DIP','Chl']
    depthints={'surface':[0,10]} #,'bottom':[10,0]} #for bottom, depthint is relative to bottom depth
    timeint = [datetime.datetime(2012, 1, 1,0,0,0), datetime.datetime(2014, 12, 31,23,59,59)]
    ##timeint = [datetime.datetime(2000, 1, 1, 0, 0, 0), datetime.datetime(2010, 12, 31, 23, 59, 59)]
    # regarding observations.
    statsets = ['BGC'] #'cosyna', 'BSH', 'BGC'
    stations = []
    readobsraw=False #i.e., if the pickle file should be ignored
    # regarding simulations.
    sims2plot = ['GF']
    readsimraw=True #i.e., if the pickle file should be ignored
    simdomain=''
    meth2D='pretree'
    #regarding plots:
    plotopts={'TS':True,'TSstyle':'TSdefault','varns':vars,'sims2plot':sims2plot}

    #READ OBSERVATIONS
    obs=SRO.readobs(pathreg[user],readobsraw,statsets,stations,timeint,depthints,vars)

    #READ SIMULATIONS
    simset={}    
    for simno, simname in enumerate(sims2plot):
        sim=SRS.readsim(pathreg[user],simname,readsimraw,simdomain,meth2D,statsets,timeint,depthints,obs,vars)
        simset[simname]=sim

    #print(simset)

    #PLOTS
    SP.stations_plots(plotopts, obs, simset, pathreg[user]['plotrootpath'], statsets, stations, timeint, depthints)

if __name__=='__main__':
    main()