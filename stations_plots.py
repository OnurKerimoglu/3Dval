# -*- coding: utf-8 -*-
"""
Created on Mon Jun  12 12:10 2017

@author: kerimoglu.o@gmail.com

"""
import os

def stations_plots(plotopts,obs,sim,plotrootpath,stations,timeint,depthints):
    code = '_%s_%s-%s' % ('-'.join(depthints.keys()), timeint[0].year, timeint[1].year)

    if plotopts['TS']==True:
        plotpath=os.path.join(plotrootpath,'3Dval_stations'+code)
        stations_plots_ts(plotopts, obs, sim, plotpath, stations, timeint, depthints, code)


def stations_plots_ts(plotopts,obs,sim,plotpath,stations,timeint,depthints,code):
    print('Doing the time series plots')

    #if the plotpath doesn't exist, create it
    #if
    print ('plot path:%s'%plotpath)
    return