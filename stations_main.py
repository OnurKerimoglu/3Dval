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

pathreg = {'onur': {#'GF-Mnm': '/home/onur/WORK/projects/2013/maecs/sns144-M180109-nFpBpr-Pbg2017-B180106-vsdetp4b1/sns144-M180109-nFpBpr-Pbg2017-B180106-vsdetp4b1_TS_12-13_zSB.nc',
                    #'GF-Mfc': '/home/onur/WORK/projects/2013/maecs/sns144-M180109-nFpBpr-Pbg2017-B180106-vsdetp4b1-AHm2-c06/extract_Mphysred_sns144-M180109-nFpBpr-Pbg2017-B180106-vsdetp4b1-AHm2-c06.12-13_zSB.nc',
                    #'GF-Mvc': '/home/onur/WORK/projects/2013/maecs/sns144-M180109-nFpBpr-Pbg2017-B180106-vsdetp4b1-AHm2-c06-Ec01/extract_MphysTS_sns144-M180109-nFpBpr-Pbg2017-B180106-vsdetp4b1-AHm2-c06-Ec01.12-13_zSB.nc',
                    #'GF-Mnm': '/home/onur/WORK/projects/2013/maecs/sns144-M180109-nFpBpr-Pbg2017-B180106-vsdetp4b1/sns144-M180109-nFpBpr-Pbg2017-B180106-vsdetp4b1-chl_12-13_zSB.nc',
                    #'GF-Mfc': '/home/onur/WORK/projects/2013/maecs/sns144-M180109-nFpBpr-Pbg2017-B180106-vsdetp4b1-AHm2-c06/extract_chlred_sns144-M180109-nFpBpr-Pbg2017-B180106-vsdetp4b1-AHm2-c06.12-13_zSB.nc',
                    #'GF-Mvc': '/home/onur/WORK/projects/2013/maecs/sns144-M180109-nFpBpr-Pbg2017-B180106-vsdetp4b1-AHm2-c06-Ec01/extract_chlred_sns144-M180109-nFpBpr-Pbg2017-B180106-vsdetp4b1-AHm2-c06-Ec01.12-13_zSB.nc',
                    #'GF-M13R12': '/home/onur/WORK/projects/2013/maecs/sns144-M180109-nFpBpr-Pbg2017-B180106-vsdetp4b1-M13R12/sns144-M180109-nFpBpr-Pbg2017-B180106-vsdetp4b1-M13R12_TS_13_zSB.nc',
                    #'GF-M12R13': '/home/onur/WORK/projects/2013/maecs/sns144-M180109-nFpBpr-Pbg2017-B180106-vsdetp4b1-M12R13/sns144-M180109-nFpBpr-Pbg2017-B180106-vsdetp4b1-M12R13_TS_13_zSB.nc',
                    #'GF-ref': '/home/onur/WORK/projects/2013/maecs/sns144-M180109-nFpBpr-Pbg2017-B180106-vsdetp4b1/sns144-M180109-nFpBpr-Pbg2017-B180106-vsdetp4b1_BGC_12-13_S.nc',
                    #'plotrootpath':'/home/onur/WORK/projects/2013/maecs/sns144-M180109-nFpBpr-Pbg2017-B180106-vsdetp4b1-AHm2-c06-Ec01-comp/',
                    'GF-PPZZ': '/home/onur/WORK/projects/2013/gpmeh/sns144-GPMEH-P190529-fSG97dChl/extract_skillC_sns144-GPMEH-P190529-fSG97dChl.2012-2013_zSB.nc',
                    'plotrootpath':'/home/onur/WORK/projects/2013/gpmeh/sns144-GPMEH-P190529-fSG97dChl/',
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
    #vars=['temp','salt'] #,'DOs']
    vars = ['DIN','DIP','Chl']
    depthints={'surface':[0,10]} #,'bottom':[10,0]} #for bottom, depthint is relative to bottom depth
    timeint = [datetime.datetime(2012, 1, 1,0,0,0), datetime.datetime(2013, 12, 31,23,59,59)]
    ##timeint = [datetime.datetime(2000, 1, 1, 0, 0, 0), datetime.datetime(2010, 12, 31, 23, 59, 59)]
    # regarding observations.
    statsets = ['BGC'] #'cosyna', 'BSH', 'BGC'
    #stations = ['Ems', 'Deutsche Bucht','NBII']
    #stations = ['Cuxhaven','HPA-Elbe']
    stations=[]
    readobsraw=False #i.e., if the pickle file should be ignored
    # regarding simulations.
    sims2plot= ['GF-PPZZ']
    #sims2plot = ['GF-Mnm','GF-Mfc','GF-Mvc'] #'GF-c100','GF-ref'] #,'GF-M13R12','GF-M12R13']
    readsimraw=False #i.e., if the pickle file should be ignored
    simdomain=''
    meth2D='pretree'
    #regarding plots:
    plotopts={'TS':True,'TSstyle':'TSdefault','varns':vars,'sims2plot':sims2plot}
    olf=4.0 #
    getmv='mean'
    fabmv='GPMEH'

    #READ OBSERVATIONS
    obs=SRO.readobs(pathreg[user],readobsraw,statsets,stations,timeint,depthints,vars,olf)

    #READ SIMULATIONS
    simset={}    
    for simno, simname in enumerate(sims2plot):
        sim=SRS.readsim(pathreg[user],simname,readsimraw,simdomain,meth2D,statsets,timeint,depthints,obs,vars,getmv,fabmv)
        simset[simname]=sim

    #print(simset)

    #PLOTS
    SP.stations_plots(plotopts, obs, simset, pathreg[user]['plotrootpath'], statsets, stations, timeint, depthints)

if __name__=='__main__':
    main()