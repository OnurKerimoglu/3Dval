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
from general_funcs import grab_data

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
                    'GF-IR': '/work/ku0646/UBA/simout-sns/sns144-GPMEH-G200124-Fnew3-PPZZSi-vS-P191223-OREF-IR-BCc/extract_skillphysC_sns144-GPMEH-G200124-Fnew3-PPZZSi-vS-P191223-OREF-IR-BCc.2015-2017_zSB.nc',
                    'GF-ref': '/work/ku0646/UBA/simout-sns/sns144-GPMEH-G200124-Fnew3-PPZZSi-vS-P191223-ICGEMO-CS/extract_skillMphysC_sns144-GPMEH-G200124-Fnew3-PPZZSi-vS-P191223-ICGEMO-CS.2012-2014_S10.nc',
                    'DCSM':'/work/ku0646/UBA/simout-deltares/DCSM-FM_2021-01-27/DCSM-FM_0_5nm_waq_0000_2012-2014_his.nc',
                    'plotrootpath':'/work/ku0646/UBA/simout-deltares/DCSM-FM_2021-01-27/',
                    #'plotrootpath': '/work/ku0646/UBA/simout-sns/sns144-GPMEH-G200124-Fnew3-PPZZSi-vS-P191223-ICGEMO-CS/',
                    'pickledobspath': './',
                    'BGC':    '/work/gg0877/onur/obsdata/stations/individual/BGC/',
                    'BSH':    '/work/gg0877/onur/obsdata/stations/individual/BSH/',
                    'cosyna': '/work/gg0877/onur/obsdata/stations/COSYNA/proc/nc',
                    'InterReg':'/work/ku0646/UBA/obsdata/stations/InterReg/allNC'
                    },
        'g260105': {
                    #'GF-PPZZ-fS': '/home/daniel/IR/sns144-GPMEH-G200124-Fnew3-PPZZSi-vS-P191223-OREF-IR-BCc/extract_skillC_sns144-GPMEH-PPZZ-P190628-fSG97dChl.2017_zSB.nc',
                    #'GF-PPZZ-vS': '/home/daniel/IR/sns144-GPMEH-G200124-Fnew3-PPZZSi-vS-P191223-OREF-IR-BCc/extract_skillC_sns144-GPMEH-PPZZ-P190628-vSG97dChl.2017_zSB.nc',
                    #'GF-PPZZ': '/home/daniel/IR/sns144-GPMEH-G200124-Fnew3-PPZZSi-vS-P191223-OREF-IR-BCc/extract_skillphysC_sns144-GPMEH-G200124-Fnew3-PPZZSi-vS-P191223-OREF-IR-BCc.2015-2017_zSB.nc',
                    'GF-v0': '/work/ku0646/g260105/IR/Harmonization/v0/extract_skillMphysC_sns144-GPMEH-G200124-Fnew3-PPZZSi-vS-P191223-ICGEMO-CS.2014_S10.nc',
                    'GF-v1': '/work/ku0646/g260105/IR/Harmonization/v1/extract_skillMphysC_sns144-GPMEH-G200124-Fnew3-PPZZSi-vS-P191223-ICGEMO-CS-rivWS.2014_S10.nc',
                    'GF-v2': '/work/ku0646/g260105/IR/Harmonization/v2/extract_skillMphysC_sns144-GPMEH-G200124-Fnew3-PPZZSi-vS-ICGEMO-CS-BCdcsmP-rivWS-DR.2014_S10.nc',
                    'GF-v3-sat': '/work/ku0646/g260105/IR/Harmonization/v3/extract_skillMphysC_sns144-GPMEH-G200124-Fnew3-PPZZSi-vS-ICGEMO-CS-BCdcsmP-rivWS.2014_S10.nc',
                    'GF-v3-bat': '/work/ku0646/g260105/IR/Harmonization/v3_ref/extract_skillMphysC_sns144-GPMEH-G200124-Fnew3-PPZZSi-vS-ICGEMO-CS-BCdcsmP-rivWS-old.2014_S10.nc',
                    'GF-v3-jt3': '/work/ku0646/g260105/IR/Harmonization/v3_con/extract_skillMphysC_sns144-GPMEH-G200124-Fnew3-PPZZSi-vS-ICGEMO-CS-BCdcsmP-rivWS-const.2014_S10.nc',
                    'GF-v4.1': '/work/ku0646/g260105/IR/Harmonization/v4.1/extract_skillMphysC_sns144-GPMEH-G200124-Fnew3-PPPMZZ-vS-ICGEMO-CS-BCdcsmP-rivWS.2014-2016_S10.nc',
                    'GF-v4.2': '/work/ku0646/g260105/IR/Harmonization/v4.2/extract_skillMphysC_sns144-GPMEH-G200124-Fnew3-PPPMZZ-vS-ICGEMO-CS-BCdcsmP-rivWS-kcl.2014-2016_S10.nc',
                    'GF-v4.3': '/work/ku0646/g260105/IR/Harmonization/v4.3/extract_skillMphysC_sns144-GPMEH-G200124-Fnew3-PPPMZZ-vS-ICGEMO-CS-BCdcsmP-rivWS-zmh3.2014_S10.nc',
                    'GF-v4.4': '/work/ku0646/g260105/IR/Harmonization/v4.4/extract_skillMphysC_sns144-GPMEH-G200124-Fnew3-PPPMZZ-vS-ICGEMO-CS-BCdcsmP-rivWS-zmh4.2014_S10.nc',
                    'GF-v4.7.1': '/work/ku0646/g260105/IR/Harmonization/v4.7.1/extract_skillMphysC_sns144-GPMEH-G200124-Fnew3-PPPMZZ-vS-ICGEMO-CS-BCdcsmP-rivWS-2.2014-2017_S10.nc',
                    'GF-v4.7.2': '/work/ku0646/g260105/IR/Harmonization/v4.7.2/extract_skillMphysC_sns144-GPMEH-G200124-Fnew3-PPPMZZ-vS-ICGEMO-CS-BCdcsmP-rivWS-2t.2014-2017_S10.nc',
                    'GF-CS': '/work/ku0646/g260105/IR/Harmonization/CS/extract_skillMphysC_sns144-GPMEH-G200124-Fnew3-PPPMZZ-vS-ICGEMO-CS-BCdcsmP-rivWS-2t.2014-2017_S10.nc',
                    'GF-28m': '/work/ku0646/g260105/IR/Harmonization/28/extract_skillMphysC_sns144-GPMEH-G200124-Fnew3-PPPMZZ-vS-ICGEMO-CS-BCdcsmP-rivWS-28.2014-2017_S10.nc',
                    'GF-28o': '/work/ku0646/g260105/IR/Harmonization/28M/extract_skillMphysC_sns144-GPMEH-G200124-Fnew3-PPPMZZ-vS-ICGEMO-CS-BCdcsmP-rivWS-28M.2014-2017_S10.nc',
                    'GF-HS1': '/work/ku0646/g260105/IR/Harmonization/HS1/extract_skillMphysC_sns144-GPMEH-G200124-Fnew3-PPPMZZ-vS-ICGEMO-HS1-BCdcsmP-rivWS.2014-2017_S10.nc',
                    'GF-HS2': '/work/ku0646/g260105/IR/Harmonization/HS2/extract_skillMphysC_sns144-GPMEH-G200124-Fnew3-PPPMZZ-vS-ICGEMO-HS2-BCdcsmP-rivWS.2014-2017_S10.nc',
                    'SNS-GPM': '/work/ku0646/g260105/IR/Harmonization/v3/extract_skillMphysC_sns144-GPMEH-G200124-Fnew3-PPZZSi-vS-ICGEMO-CS-BCdcsmP-rivWS.2014-2016_S10.nc',                    
                    'DCSM-CS': '/work/ku0646/g260105/IR/Harmonization/DCSM/DCSM-FM_0_5nm_waq_0000_2014-2017_CS.nc',
                    'DCSM-28m': '/work/ku0646/g260105/IR/Harmonization/DCSM/DCSM-FM_0_5nm_waq_0000_2014-2017_28.nc',
                    #'GF-old': '/home/daniel/IR/sns144-GPMEH-G200124-Fnew3-PPZZSi-vS-P191223-OREF-IR-BCc/extract_skillC_sns144-GPMEH-G191216-Fnew3-PPZZSi-vS-P191223.2015-2017_zSB.nc',
                    'plotrootpath':'/work/ku0646/g260105/IR/Harmonization/',
                    'pickledobspath': './',
                    'BGC':    '/home/daniel/IR/stations/individual/BGC/',
                    'BSH':    '/home/daniel/IR/stations/individual/BSH/',
                    'cosyna': '/home/daniel/IR/stations/COSYNA/proc/nc',
                    #'InterReg':'/work/ku0646/g260105/IR/stations/InterReg/allNC'
                    'InterReg':'/work/ku0646/g260105/IR/Harmonization/allNC',
                    'InterRegFG':'/work/ku0646/g260105/IR/Harmonization/allNC'
                    },
        
         'daniel': {
                    'GF-CS': '/home/daniel/levante_work/IR/Harmonization/CS/extract_skillMphysC_sns144-GPMEH-G200124-Fnew3-PPPMZZ-vS-ICGEMO-CS-BCdcsmP-rivWS-2t.2014-2017_S10.nc',
                    'GF-28m': '/home/daniel/levante_work/IR/Harmonization/28/extract_skillMphysC_sns144-GPMEH-G200124-Fnew3-PPPMZZ-vS-ICGEMO-CS-BCdcsmP-rivWS-28.2014-2017_S10.nc',
                    'GF-28o': '/home/daniel/levante_work/IR/Harmonization/28M/extract_skillMphysC_sns144-GPMEH-G200124-Fnew3-PPPMZZ-vS-ICGEMO-CS-BCdcsmP-rivWS-28M.2014-2017_S10.nc',
                    'GF-HS1': '/home/daniel/levante_work/IR/Harmonization/HS1/extract_skillMphysC_sns144-GPMEH-G200124-Fnew3-PPPMZZ-vS-ICGEMO-HS1-BCdcsmP-rivWS.2014-2017_S10.nc',
                    'GF-HS2': '/home/daniel/levante_work/IR/Harmonization/HS2/extract_skillMphysC_sns144-GPMEH-G200124-Fnew3-PPPMZZ-vS-ICGEMO-HS2-BCdcsmP-rivWS.2014-2017_S10.nc',
                    'SNS-GPM': '/home/daniel/levante_work/IR/Harmonization/v3/extract_skillMphysC_sns144-GPMEH-G200124-Fnew3-PPZZSi-vS-ICGEMO-CS-BCdcsmP-rivWS.2014-2016_S10.nc',                    
                    'DCSM-CS': '/home/daniel/levante_work/IR/Harmonization/DCSM/DCSM-FM_0_5nm_waq_0000_2014-2017_CS.nc',
                    'DCSM-28m': '/home/daniel/levante_work/IR/Harmonization/DCSM/DCSM-FM_0_5nm_waq_0000_2014-2017_28.nc',
                    'FSK-CS': '/home/daniel/levante_work/IR/Harmonization/FSK/FSK_waq_2017_CS_daily.nc',
                    'FSK-28m': '/home/daniel/levante_work/IR/Harmonization/FSK/FSK_waq_2017_2o8_daily.nc',
                    'FSK-CSwPr': '/home/daniel/levante_work/IR/Harmonization/FSK/FSK_waq_2017_CSwPr_daily.nc',
                    'FSK-28mwPr': '/home/daniel/levante_work/IR/Harmonization/FSK/FSK_waq_2017_2o8wPr_daily.nc',
                    #'GF-old': '/home/daniel/IR/sns144-GPMEH-G200124-Fnew3-PPZZSi-vS-P191223-OREF-IR-BCc/extract_skillC_sns144-GPMEH-G191216-Fnew3-PPZZSi-vS-P191223.2015-2017_zSB.nc',
                    'plotrootpath':'/home/daniel/levante_work/IR/Harmonization/',
                    'pickledobspath': '/home/daniel/levante_work/IR/Harmonization/allNC',
                    'GF-CSl':  '/home/daniel/levante_work2/sns144-GPMEH-G200124-Fnew3-PPPMZZ-vS-ICGEMO-CS-BCdcsmP-rivWS-2t/extract_var2C_sns144-GPMEH-G200124-Fnew3-PPPMZZ-vS-ICGEMO-CS-BCdcsmP-rivWS-2t.2017.nc',
                    'GF-28ml': '/home/daniel/levante_work2/sns144-GPMEH-G200124-Fnew3-PPPMZZ-vS-ICGEMO-CS-BCdcsmP-rivWS-28/extract_var2C_sns144-GPMEH-G200124-Fnew3-PPPMZZ-vS-ICGEMO-CS-BCdcsmP-rivWS-28.2017.nc',
                    'GF-HS1l': '/home/daniel/levante_work2/sns144-GPMEH-G200124-Fnew3-PPPMZZ-vS-ICGEMO-HS1-BCdcsmP-rivWS/extract_var2C_sns144-GPMEH-G200124-Fnew3-PPPMZZ-vS-ICGEMO-HS1-BCdcsmP-rivWS.2017.nc',
                    'GF-HS2l': '/home/daniel/levante_work2/sns144-GPMEH-G200124-Fnew3-PPPMZZ-vS-ICGEMO-HS2-BCdcsmP-rivWS/extract_var2C_sns144-GPMEH-G200124-Fnew3-PPPMZZ-vS-ICGEMO-HS2-BCdcsmP-rivWS.2017.nc',
                    'cosyna': '/home/daniel/IR/stations/COSYNA/proc/nc',
                    #'InterReg':'/home/daniel/levante_work/IR/stations/InterReg/allNC'
                    'InterReg':'/home/daniel/levante_work/IR/Harmonization/allNC',
                    'InterRegFG':'/home/daniel/levante_work/IR/Harmonization/allNC',
                    'InterRegFGlim-diat':'/home/daniel/levante_work/IR/Harmonization/allNC',
                    # 'InterReg':'/home/daniel/levante_work/IR/Harmonization/NWDM',
                    # 'InterRegFG':'/home/daniel/levante_work/IR/Harmonization/NWDM'
                    'SNS-r': '/home/daniel/levante_work/IR/sns-oe/simout-gpmeh/sns144-r/extract_skillMphysC_sns144-r.2014-2017_S10.nc',
                    'SNS-t2': '/home/daniel/levante_work/IR/sns-oe/simout-gpmeh/sns144-t2/extract_skillMphysC_sns144-t2.2014-2017_S10.nc',
                    'SNS-t4': '/home/daniel/levante_work/IR/sns-oe/simout-gpmeh/sns144-t4/extract_skillMphysC_sns144-t4.2014-2017_S10.nc',
                    'SNS-2g-CS': '/home/daniel/levante_work/IR/sns-oe/simout-gpmeh/sns144-2g-CS/extract_skillMphysC_sns144-2g-CS.2014-2017_S10.nc',
                    'SNS-4g-CS': '/home/daniel/levante_work/IR/sns-oe/simout-gpmeh/sns144-4g-CS/extract_skillMphysC_sns144-4g-CS.2014-2017_S10.nc',
                    'SNS-2g-CS-NEC': '/home/daniel/levante_work/IR/sns-oe/simout-gpmeh/sns144-2g-CS-NEC/extract_skillMphysC_sns144-2g-CS-NEC.2014-2017_S10.nc',
                    'SNS-4g-CS-NEC': '/home/daniel/levante_work/IR/sns-oe/simout-gpmeh/sns144-4g-CS-NEC/extract_skillMphysC_sns144-4g-CS-NEC.2014-2017_S10.nc'
                    }
           }

def main(modtype,modfname,statsets,yint):
    #PARAMETERS:
    # general
    # user='g260105' #daniel,g260105,g260108,onur
    user = 'daniel'
    vars_default = ['temp', 'salt', 'DOs', 'DIN', 'DIP', 'Chl']
    depthints={'surface': [0, 10]} #,'bottom':[10,0]} #for bottom, depthint is relative to bottom depth
    #timeint = [datetime.datetime(2012, 1, 1,0,0,0), datetime.datetime(2013, 12, 31,23,59,59)]
    timeint = [datetime.datetime(yint[0], 1, 1,0,0,0), datetime.datetime(yint[1], 12, 31,23,59,59)]
    # regarding observations.
    if len(statsets)==0:
       #statsets = ['cosyna', 'BSH', 'BGC']
       statsets = ['InterRegFG']
       # statsets = ['InterReg']
       #statsets = ['InterRegFGlim-diat']
       #stations = ['Ems', 'Deutsche Bucht','NBII']
       #stations = ['Cuxhaven','HPA-Elbe']
    stations = []
    #     stations = ['BOCHTVWTM','BOOMKDP','Bork_W_1','DANTZGT','DOOVBWT','GROOTGND','HUIBGOT','JaBu_W_1','MARSDND',
    #                    'Nney_W_1','Nney_W_2', 'Nney_W_3','ROTTMPT3','ROTTMPT70','TERSLG4','TERSLG10','WeMu_W_1',
    #                    'WeMu_W2', 'Wu_KU-W1','ZOUTKPLZGT','ZUIDOLWOT','220006', '220017', '220052_B', '220052_S',
    #                    '220057', '220065']
    if len(statsets)==1 and statsets[0]=='cosyna':
       #vars=['DOs']
       vars=['temp', 'salt'] #,'DOs']
    elif len(statsets)==1 and statsets[0]=='BGC':
       vars=['DIN', 'DIP', 'Si', 'Chl']
    elif len(statsets)==1 and statsets[0]=='InterReg':
       #vars=['salt','DIN','DIP','Si','Chl']
       vars=['salt', 'DIN', 'DIP', 'Si', 'Chl']
    elif len(statsets)==1 and statsets[0]=='InterRegFG':
       vars=['Chl', 'Diatoms', 'Flagellates', 'Dinoflagellates', 'Phaeocystis']#'Cyanobacteria','other']
       # vars=['Chl', 'Diatoms', 'Flagellates']#'Cyanobacteria','other']
    elif len(statsets) == 1 and statsets[0] == 'InterRegFGlim-diat':
       vars = ['diat_limI', 'diat_limN', 'diat_limP', 'diat_limSi']  # 'Cyanobacteria','other']
    #elif len(statsets)==1 and statsets[0]=='InterRegFG':
       #vars=['Chl',]#'Cyanobacteria','other']
    else:
       vars=vars_default
 
    # regarding simulations.
    #only if modfname==''
    #sims2plot=['DCSM','SNS-GPM']
    # sims2plot = ['SNS-2g-CS', 'SNS-4g-CS']
    # sims2plot = ['SNS-2g-CS-NEC', 'SNS-4g-CS-NEC']
    # sims2plot = ['SNS-2g-CS', 'SNS-2g-CS-NEC']
    # sims2plot = ['SNS-4g-CS', 'SNS-4g-CS-NEC']
    sims2plot = ['SNS-2g-CS', 'SNS-2g-CS-NEC', 'SNS-4g-CS', 'SNS-4g-CS-NEC']
    # sims2plot=['DCSM-CS', 'GF-CS']
    # sims2plot=['GF-CS', 'FSK-CS', 'FSK-CSwPr']
    # sims2plot=['DCSM-28m','GF-28m', 'GF-28o']
    # sims2plot=['DCSM-CS', 'GF-CS', 'FSK-CS', 'FSK-CSwPr']
    #sims2plot = ['GF-CSl', 'GF-28ml', 'GF-HS2l']
    # sims2plot=['GF-CS','GF-HS1','GF-HS2']
    # sims2plot=['DCSM-CS','DCSM-28m']
    # sims2plot = ['SNS-r', 'SNS-t4']
    #sims2plot = ['GF-v3-jt3','GF-v3-bat','GF-v3-sat']
    #sims2plot=['GF-ref', 'GF-R12', 'GF-M12', 'GF-W12']
    #sims2plot= ['GF-PPZZ-fS', 'GF-PPZZ-vS']
    #sims2plot = ['GF-Mnm','GF-Mfc','GF-Mvc'] #'GF-c100','GF-ref'] #,'GF-M13R12','GF-M12R13']
    readsimraw = False #i.e., if the pickle file should be ignored
    readobsraw = False #i.e., if the pickle file should be ignored
    simdomain = ''
    meth2D = 'pretree'
    #regarding plots:
    olf = 4.0  #
    getmv = 'mean' #mean,3d
    pathreg_u = pathreg[user]
    
    #READ OBSERVATIONS
    use_NWDM = False
    obs = SRO.readobs(pathreg_u, readobsraw, statsets, stations, timeint, depthints, vars, olf, use_NWDM)
    
    if not modfname == '':
        pathreg_u[modtype] = modfname
        sims2plot = [modtype]
        pathreg_u['plotrootpath'] = os.path.dirname(modfname)

    #READ SIMULATIONS
    simset={}    
    for simno, simname in enumerate(sims2plot):
        sim=SRS.readsim(pathreg_u, simname, readsimraw, simdomain, meth2D, statsets, timeint, depthints, obs, vars, getmv, modtype)
        simset[simname] = sim
    plotopts={'TS': True, 'TSstyle':'TSdefault','varns':vars,'sims2plot': sims2plot, 'plot_obs_override': True, 'plot_obs_force': True}

    # grab data, if necessary
    if grab_data_switch:
        varlist = ['DIN', 'DIP', 'Si', 'Chl']
        # statoutlist = ['BOCHTVWTM', 'BOOMKDP', 'Bork_W_1', 'DANTZGT', 'DOOVBWT', 'GROOTGND', 'HUIBGOT', 'JaBu_W_1', 'MARSDND',
        #                'Nney_W_1','Nney_W_2', 'Nney_W_3','ROTTMPT3','ROTTMPT70','TERSLG4','TERSLG10', 'WeMu_W_1',
        #                'WeMu_W2', 'Wu_KU-W1','ZOUTKPLZGT','ZUIDOLWOT','220006', '220017', '220052_B', '220052_S',
        #                '220057', '220065']
        statoutlist = ['BOCHTVWTM', 'BOOMKDP', 'Bork_W_1', 'DANTZGT', 'DOOVBWT', 'GROOTGND', 'HUIBGOT', 'JaBu_W_1', 'MARSDND',
                       'Nney_W_2', 'ROTTMPT3', 'ROTTMPT50', 'ROTTMPT70', 'TERSLG4', 'TERSLG10', 'TERSLG100', 'TERSLG135', 'TERSLG175',
                       'TERSLG235', 'TERSLG50', 'WeMu_W_1']
        grab_data(simset, varlist, pathreg_u['plotrootpath'], statoutlist)

    #PLOTS
    SP.stations_plots(plotopts, obs, simset, pathreg_u['plotrootpath'], statsets, stations, timeint, depthints)

if __name__=='__main__':
    if len(sys.argv)>1:
       modtype=sys.argv[1]
    else:
       modtype='GF-PPPMZZ'
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
       # yint=[2017, 2017]
       yint=[2014, 2017]

    if len(sys.argv) > 5:
       grab_data_switch = sys.argv[5]
    else:
       grab_data_switch = False

    main(modtype, modfname, statsets, yint)
