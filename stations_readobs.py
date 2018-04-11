# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 16:06:57 2017

@author: kerimoglu.o@gmail.com
"""
import os
import pickle
import numpy as np
import netCDF4
from general_funcs import get_botdepth

def readobs(paths,readraw,statsets,stations,timeint,depthints,vars=[]):
    print ('Reading observations:')

    pickleobs=True
    if pickleobs:
        picklecode = '_%s_%s_%s-%s' % ('-'.join(statsets), '-'.join(depthints.keys()), timeint[0].year, timeint[1].year)
        pickledobsfile = os.path.join(paths['pickledobspath'], 'obs' + picklecode + '.pickle')
        #pickledobsfile = os.path.join(paths[statset], 'obs' + picklecode + '.pickle')

        #check if the pickledobs exist
        if (not readraw) and os.path.exists(pickledobsfile):
           print('Opening pickled obs file')
           (obs,) = np.load(pickledobsfile)
           # filter requested stations?
           return obs

        # if pickledfile does not exist:
        obs={}
        for statset in statsets:
            print('Filling obs. dataset:%s'%statset)
            obs=looppath_fill_stationdata_obs(obs,paths,statset,stations,timeint,depthints,vars)

        #pickle the obs file
        f=open(pickledobsfile,'wb')
        pickle.dump((obs,),f) #,protocol=-1
        print('Pickled obs file for later use:' + pickledobsfile)
        f.close()
    else:
        for statset in statsets:
            print('Filling obs. dataset:%s' % statset)
            obs = looppath_fill_stationdata_obs({}, paths, statset, stations, timeint, depthints,vars)

    return obs

def looppath_fill_stationdata_obs(obs,paths,statset,stations,timeint,depthints,vars):

    obspath = paths[statset]

    # if stations to include are not specified,
    if len(stations) == 0:
        # get lists of available data
        files = [f for f in os.listdir(obspath) if f.endswith('.nc')]
        sfiles = files[:]
    else:
        #TODO: check the station attribute of each file, return the needed file names
        raise(Exception('sfile generation based on stations not yet implemented. leave the station list empty'))

    for sfile in sfiles:
        # station=sfile.split('.nc')[0]
        print('  ' + sfile)
        if len(vars)==0:
            vars = ['temp', 'salt', 'DOs', 'ssh']
        sdata,station = fill_stationdata_obs(os.path.join(obspath,sfile),statset,vars,timeint,depthints)
        obs[station] = sdata

    return obs

def fill_stationdata_obs(file,statset,vars,timeint,depthints0):

    if statset in ['cosyna']:
        vlib = {'t': 'time', 'x': 'lon', 'y': 'lat', 'z': 'depth', 'temp': 'temp', 'salt': 'salt','DOs': 'DOsat', 'ssh':'ssh'}
    elif statset in ['BSH']:
        vlib = {'t': 'time', 'x': 'lon', 'y': 'lat', 'z': 'depth', 'temp': 'temp', 'salt': 'sal','DOs': 'DOsat'}
    elif statset in ['BGC']:
        vlib = {'t': 'time', 'x': 'lon', 'y': 'lat', 'z': 'depth', 'Chl': 'chl', 'DIN': 'DIN', 'DIP': 'DIP'}

    #variables without depth dimension
    noZDvars = ['ssh']

    #open the data set
    ncf = netCDF4.Dataset(file,'r')
    station=ncf.station

    # reading metadata of station
    lon = ncf.variables[vlib['x']][:][0]
    lat = ncf.variables[vlib['y']][:][0]
    if vlib['z'] in ncf.variables: #if a depth dimension exists, extract it
        depth = ncf.variables[vlib['z']][:]
    else: #assume that all are surface depths
        depth=-9999*np.ones(1)
    time_num = ncf.variables[vlib['t']][:]
    time = netCDF4.num2date(time_num, ncf.variables[vlib['t']].getncattr('units'))

    # find the max_depth, if necessary
    try:
        maxz=ncf.bottom_depth
    except:
        maxz = get_botdepth(lon, lat)

    # update the bottom depthint
    depthints = depthints0.copy()
    if 'bottom' in depthints.keys() and not np.isnan(maxz):
        depthints['bottom']=[maxz-depthints0['bottom'][0], maxz-depthints0['bottom'][1]]
        print ('updated depthints for bottom: %s-%s'%(depthints['bottom'][0],depthints['bottom'][1]))

    # put all data in
    sdata = {'longname': 'descriptive name of the station',
             'lon': lon,
             'lat': lat,
             'bottom_depth': maxz,
             }

    # check if dates are in
    tind=np.where((time>=timeint[0]) * (time<=timeint[1]))[0]

    # check if depths are in
    if list(depth) != [-9999]:
        depthintmin=np.nanmin([dint[0] for dint in depthints.values()]) #find the minimum lower lim of depthints
        depthintmax = np.nanmax([dint[1] for dint in depthints.values()]) #find the maximum upper lim of depthints
        zind=np.where((depth>=depthintmin) * (depth<=depthintmax))[0]
        depthsin=True
    else:
        depthsin = False
        #zind=-1*np.ones(1)

    # for each variable, fill in the data, if exists
    for var in vars:
        sdata[var]={}
        if len(tind)>0 and depthsin and vlib[var] in ncf.variables:
            sdata[var]['presence'] = True
            if (var in noZDvars) or (list(zind) == [-1]): #if a variable with no vertical dimension:
                sdata[var]['surface']={'time':time, 'value':ncf.variables[vlib[var]][tind,0,0], 'depth_interval':[0,0]}
            else: #if vertical dimension (may) exist
                for layername, depthint in depthints.items():
                    zind = np.where((depth >= depthint[0]) * (depth <= depthint[1]))[0]
                    if len(zind) > 0:
                        #calculate average over zind
                        vals=np.nanmean(ncf.variables[vlib[var]][tind,zind,0,0],axis=1)
                        sdata[var][layername] = {'time': time[tind], 'value': vals, 'depth_interval': depthint}
                    else:
                        sdata[var][layername] = {'time': np.array([]), 'value': np.array([]), 'depth_interval': depthint}
        else:
            sdata[var]['presence']=False

    ncf.close()

    return sdata,station