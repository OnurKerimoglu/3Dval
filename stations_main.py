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
# import cftime
import numpy as np
import numpy.ma as ma
#home=os.getcwd()
#sys.path.insert(1, home+'/3Dsetups/postprocess')
#sys.path.insert(1, home+'/3Dval')
import stations_readobs as SRO
import stations_readsim as SRS
import stations_plots as SP

#to reload modules
import importlib
importlib.reload(SRO)
importlib.reload(SRS)
importlib.reload(SP)

only_use_cftime_dateimes=False
only_use_python_datetime=True
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
                    'GF-3DFnew': '/home/onur/WORK/projects/2013/gpmeh/sns144-GPMEH-G191216-Fnew3-PPZZSi-P191220-vS/extract_skillC_sns144-GPMEH-G191216-Fnew3-PPZZSi-P191220-vS.2011_zSB.nc',
                    'plotrootpath':'/home/onur/WORK/projects/2013/gpmeh/sns144-GPMEH-G191216-Fnew3-PPZZSi-P191220-vS/',
                    'DCSM':'/home/onur/WORK/projects/InterregA/model-comparison/DCSM-FM_2021-01-27/DCSM-FM_0_5nm_waq_0000_2012-2014_his.nc',
                    'plotrootpath':'/home/onur/WORK/projects/InterregA/model-comparison/',
                    'pickledobspath': './',
                    'BGC':    '/home/onur/WORK/projects/GB/data/stations/individual/BGC/',
                    'BSH':    '/home/onur/WORK/projects/GB/data/stations/individual/BSH/',
                    'cosyna': '/home/onur/WORK/projects/GB/data/stations/COSYNA/proc/nc',
                    'InterReg':'/home/onur/WORK/projects/GB/data/stations/InterReg/allNC'
                   },
        'g260108': {
                    'GF-PPZZ-fS': '/work/gg0877/onur/simout-gpmeh/sns144-GPMEH-PPZZ-P190628-fSG97dChl/extract_skillC_sns144-GPMEH-PPZZ-P190628-fSG97dChl.2012_zSB.nc',
                    'GF-PPZZ-vS': '/work/gg0877/onur/simout-gpmeh/sns144-GPMEH-PPZZ-P190628-vSG97dChl/extract_skillC_sns144-GPMEH-PPZZ-P190628-vSG97dChl.2012_zSB.nc',
                    'GF-ref': '/work/ku0646/UBA/simout-sns/sns144-GPMEH-G200124-Fnew3-PPZZSi-vS-P191223-OREF-IR-BCc/extract_skillphysC_sns144-GPMEH-G200124-Fnew3-PPZZSi-vS-P191223-OREF-IR-BCc.2015-2017_zSB.nc',
                    'GF-old': '/work/gg0877/onur/simout-gpmeh/sns144-GPMEH-G191216-Fnew3-PPZZSi-vS-P191223/extract_skillC_sns144-GPMEH-G191216-Fnew3-PPZZSi-vS-P191223.2010-2014_zSB.nc',
                    'plotrootpath':'/work/gg0877/onur/simout-gpmeh/sns144-GPMEH-G191216-Fnew3-PPZZSi-vS-P191223/3Dval_stations_2012-2013_scens',
                    'pickledobspath': './',
                    'BGC':    '/work/gg0877/onur/obsdata/stations/individual/BGC/',
                    'BSH':    '/work/gg0877/onur/obsdata/stations/individual/BSH/',
                    'cosyna': '/work/gg0877/onur/obsdata/stations/COSYNA/proc/nc',
                    'InterReg':'/work/ku0646/UBA/obsdata/stations/InterReg/allNC'
                    },
        'g260105': {
                    'GF-PPZZ-fS': '/home/daniel/IR/sns144-GPMEH-G200124-Fnew3-PPZZSi-vS-P191223-OREF-IR-BCc/extract_skillC_sns144-GPMEH-PPZZ-P190628-fSG97dChl.2017_zSB.nc',
                    'GF-PPZZ-vS': '/home/daniel/IR/sns144-GPMEH-G200124-Fnew3-PPZZSi-vS-P191223-OREF-IR-BCc/extract_skillC_sns144-GPMEH-PPZZ-P190628-vSG97dChl.2017_zSB.nc',
                    'GF-PPZZ': '/home/daniel/IR/sns144-GPMEH-G200124-Fnew3-PPZZSi-vS-P191223-OREF-IR-BCc/extract_skillphysC_sns144-GPMEH-G200124-Fnew3-PPZZSi-vS-P191223-OREF-IR-BCc.2015-2017_zSB.nc',
                    'GF-ref': '/home/daniel/IR/sns144-GPMEH-G200124-Fnew3-PPZZSi-vS-P191223-OREF-IR-BCc/extract_skillphysC_sns144-GPMEH-G200124-Fnew3-PPZZSi-vS-P191223-OREF-IR-BCc.2015-2017_zSB.nc',
                    'GF-old': '/home/daniel/IR/sns144-GPMEH-G200124-Fnew3-PPZZSi-vS-P191223-OREF-IR-BCc/extract_skillC_sns144-GPMEH-G191216-Fnew3-PPZZSi-vS-P191223.2015-2017_zSB.nc',
                    'plotrootpath':'/home/daniel/IR/3Dval',
                    'pickledobspath': './',
                    'BGC':    '/home/daniel/IR/stations/individual/BGC/',
                    'BSH':    '/home/daniel/IR/stations/individual/BSH/',
                    'cosyna': '/home/daniel/IR/stations/COSYNA/proc/nc',
                    'InterReg':'/home/daniel/IR/stations/InterReg/allNC'
                    },
        
         'newuser': {}
           }

def main(modtype,modfname,statsets,yint):
    #PARAMETERS:
    # general
    user='onur' #g260105,g260108
    vars_default=['temp','salt','DOs','DIN','DIP','Chl']
    depthints={'surface':[0,10]} #,'bottom':[10,0]} #for bottom, depthint is relative to bottom depth
    #timeint = [datetime.datetime(2012, 1, 1,0,0,0), datetime.datetime(2013, 12, 31,23,59,59)]
    timeint = [datetime.datetime(yint[0], 1, 1,0,0,0), datetime.datetime(yint[1], 12, 31,23,59,59)]
    # regarding observations.
    if len(statsets)==0:
       #statsets = ['cosyna', 'BSH', 'BGC']
       statsets = ['InterReg']
       #stations = ['Ems', 'Deutsche Bucht','NBII']
       #stations = ['Cuxhaven','HPA-Elbe']
    stations=[]
    if len(statsets)==1 and statsets[0]=='cosyna':
       #vars=['DOs']
       vars=['temp','salt'] #,'DOs']
    elif len(statsets)==1 and statsets[0]=='BGC':
       vars=['DIN','DIP','Si','Chl']
    elif len(statsets)==1 and statsets[0]=='InterReg':
       vars=['salt','DIN','DIP','Si','Chl']
    else:
       vars=vars_default
 
    # regarding simulations.
    #only if modfname==''
    #sims2plot=['GF-ref']
    sims2plot = ['DCSM']
    #sims2plot=['GF-ref', 'GF-R12', 'GF-M12', 'GF-W12']
    #sims2plot= ['GF-PPZZ-fS', 'GF-PPZZ-vS']
    #sims2plot = ['GF-Mnm','GF-Mfc','GF-Mvc'] #'GF-c100','GF-ref'] #,'GF-M13R12','GF-M12R13']
    readsimraw=True #i.e., if the pickle file should be ignored
    readobsraw=True #i.e., if the pickle file should be ignored
    simdomain=''
    meth2D='pretree'
    #regarding plots:
    olf=4.0 #
    getmv='mean' #mean,3d
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
        sim=SRS.readsim(pathreg_u,simname,readsimraw,simdomain,meth2D,statsets,timeint,depthints,obs,vars,getmv,modtype)
        simset[simname]=sim
    plotopts={'TS':True,'TSstyle':'TSdefault','varns':vars,'sims2plot':sims2plot}

    #print(simset)

    #PLOTS
    SP.stations_plots(plotopts, obs, simset, pathreg_u['plotrootpath'], statsets, stations, timeint, depthints)

if __name__=='__main__':
    if len(sys.argv)>1:
       modtype=sys.argv[1]
    else:
       modtype='GF-PPZZ'
    if len(sys.argv)>2:
       modfname=sys.argv[2]
    else:
       modfname=''
    
    if len(sys.argv)>3:
       statsets=sys.argv[3].split(',')
    else:
       statsets=[]
    
    if len(sys.argv)>4:
       yints=sys.argv[4].split(',')
       yint=[np.int(y) for y in yints]
    else:
       yint=[2012,2013]

    main(modtype,modfname,statsets,yint)
