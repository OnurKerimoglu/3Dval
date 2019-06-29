# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 15:02:46 2017

Main file for stations validation.
Control of structre, call reading obs. model and plot.

stations_main.py modtype simfname
example call from shell:
python3.5 stations_main.py GF-PPZZ simfname.nc

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
                    #'GF-PPZZ': '/home/onur/WORK/projects/2013/gpmeh/sns144-GPMEH-P190529-fSG97dChl/extract_skillC_sns144-GPMEH-P190529-fSG97dChl.2012-2013_zSB.nc',
                    #'plotrootpath':'/home/onur/WORK/projects/2013/gpmeh/sns144-GPMEH-P190529-fSG97dChl/',
                    'pickledobspath': './',
                    'BGC':    '/home/onur/WORK/projects/GB/data/stations/individual/BGC/',
                    'BSH':    '/home/onur/WORK/projects/GB/data/stations/individual/BSH/',
                    'cosyna': '/home/onur/WORK/projects/GB/data/stations/COSYNA/proc/nc'
                   },
        'g260108': {
                    #'GF-PPZZ': '/work/gg0877/onur/simout-gpmeh/sns144-GPMEH-PPZZ-P190529-fSG97dChl/extract_skillC_sns144-GPMEH-PPZZ-P190529-fSG97dChl.2012-2013_zSB.nc',
                    #'plotrootpath':'/work/gg0877/onur/simout-gpmeh/sns144-GPMEH-PPZZ-P190529-fSG97dChl/',
                    'pickledobspath': './',
                    'BGC':    '/work/gg0877/onur/obsdata/stations/individual/BGC/',
                    'BSH':    '/work/gg0877/onur/obsdata/stations/individual/BSH/',
                    'cosyna': '/work/gg0877/onur/obsdata/stations/COSYNA/proc/nc'
                    },
         'newuser': {}
           }

def main(modtype,modfname):
    #PARAMETERS:
    # general
    user='g260108'
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
    #sims2plot= ['GF-PPZZ']
    sims2plot = ['GF-Mnm','GF-Mfc','GF-Mvc'] #'GF-c100','GF-ref'] #,'GF-M13R12','GF-M12R13']
    readsimraw=False #i.e., if the pickle file should be ignored
    simdomain=''
    meth2D='pretree'
    #regarding plots:
    olf=4.0 #
    getmv='mean'
    fabmv='GPMEH'
    
    pathreg_u=pathreg[user]
    
    #READ OBSERVATIONS
    obs=SRO.readobs(pathreg_u,readobsraw,statsets,stations,timeint,depthints,vars,olf)
    
    if not modfname == '':
        pathreg_u[modtype]=modfname
        sims2plot=[modtype]
        pathreg_u['plotrootpath']=os.path.dirname(modfname)

    #READ SIMULATIONS
    simset={}    
    for simno, simname in enumerate(sims2plot):
        sim=SRS.readsim(pathreg_u,simname,readsimraw,simdomain,meth2D,statsets,timeint,depthints,obs,vars,getmv,fabmv)
        simset[simname]=sim
    plotopts={'TS':True,'TSstyle':'TSdefault','varns':vars,'sims2plot':sims2plot}

    #print(simset)

    #PLOTS
    SP.stations_plots(plotopts, obs, simset, pathreg_u['plotrootpath'], statsets, stations, timeint, depthints)

if __name__=='__main__':
    if len(sys.argv)>1:
       modtype=sys.argv[1]
    else:
       modtype=''

    if len(sys.argv)>2:
       modfname=sys.argv[2]
    else:
       modfname=''
    
    main(modtype,modfname)
