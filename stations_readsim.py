# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 18:55:00 2017

@author: kerimoglu.o@gmail.com

"""
import os
import pickle
import numpy as np
from netCDF4 import num2date
from netCDF4 import Dataset

def readsim(paths,simname,statsets,timeint,depthints,obs):
    print('Reading simulation:'+simname)

    simf=paths[simname]
    simfb=os.path.basename(simf)
    simfp=os.path.dirname(simf)
    picklecode = '_%s_%s_%s-%s' % ('-'.join(statsets), '-'.join(depthints.keys()), timeint[0].year, timeint[1].year)
    pickledsimfile = os.path.join(simfp, simfb.split('.nc')[0] + picklecode + '.pickle')

    # check if the pickledsim exist
    if os.path.exists(pickledsimfile):
        print('Opening pickled sim file')
        (sim,) = np.load(pickledsimfile)
        return sim

    # if pickledsim does not exist:
    sim = {}
    stations=obs.keys()
    for station in stations:
        lon = obs[station]['lon']
        lat = obs[station]['lat']
        maxz_obs= obs[station]['bottom_depth']
        sim = fill_stationdata_sim(sim,simf, simname, station, lon, lat, maxz_obs, timeint, depthints)

    # pickle the sim file
    f = open(pickledsimfile, 'wb')
    pickle.dump((sim,), f)  # ,protocol=-1
    print('Pickled sim file for later use:'+pickledsimfile)
    f.close()

    return sim

def fill_stationdata_sim(sim,simf, simname, station, lon, lat, maxz_obs, timeint, depthints):
    print('  ' + station)
    if simname[0:4]=='GETM':
        ncf = Dataset(simf)
        sdata=get_station_data_getm(ncf,lon,lat,maxz_obs,timeint,depthints)
        sim[station]=sdata
        ncf.close()
    elif simname[0:5]=='FVCOM':
        #sdata=get_station_data_fvcom() #TODO:IVAN
        sim[station] = sdata
    return sim

def get_station_data_getm(ncf,lon,lat,maxz_obs,timeint,depthints):

    vlib = {'t': 'time', 'z': 'depth', 'temp': 'temp', 'salt': 'salt', 'ssh': 'elev'}

    #get lons,lats,z_max
    #topo=get_topo()

    tempfound = False
    saltfound = False
    sshfound = False

    # reading metadata of station
    depth = ncf.variables[vlib['z']][:][0]
    time_num = ncf.variables[vlib['t']][:]
    time = num2date(time_num, ncf.variables[vlib['t']].getncattr('units'))

    # find lon lat indices: TODO
    #topo=
    #loni,lati=
    loni = 50; lati = 50

    # find the max_depth, if necessary
    if 'bottom' in depthints.keys():
        #maxz =  topo['H'][loni,lati] #find maxz: todo
        maxz=100
        #todo: perhaps compare with maxz_obs, and do something about it?
        # update the bottom depthint
        depthints['bottom'] = [maxz - depthints['bottom'][0], maxz - depthints['bottom'][1]]
    else:
        maxz=np.nan

    # check if dates are in
    tind = np.where((time >= timeint[0]) * (time <= timeint[1]))[0]

    # check if depths are in
    depthintmin = np.min([dint[0] for dint in depthints.values()])  # find the minimum lower lim of depthints
    depthintmax = np.max([dint[1] for dint in depthints.values()])  # find the maximum upper lim of depthints
    zind = np.where((depth >= depthintmin) * (depth <= depthintmax))[0]

    # check if variables are in, decide if anything relevant found
    if len(tind) > 0 and len(zind) > 0 and vlib['temp'] in ncf.variables: tempfound = True
    if len(tind) > 0 and len(zind) > 0 and vlib['salt'] in ncf.variables: saltfound = True
    if len(tind) > 0 and len(zind) > 0 and vlib['ssh'] in ncf.variables: sshfound = True

    #fill in the data:
    tempdata = {}; saltdata = {}; sshdata = {}

    # handle temp
    if tempfound:
        tempdata['presence'] = True
        for layername, depthint in depthints.items():
            zind = np.where((depth >= depthint[1]) * (depth <= depthint[1]))[0]
            if len(zind) > 0:
                data = ncf.variables[vlib['temp']][tind, zind,loni,lati]
                tempdata[layername] = {'time': time[tind], 'value': data, 'depth_interval': depthint}
            else:
                tempdata[layername] = {'time': np.array([]), 'value': np.array([]), 'depth_interval': depthint}
    else:
        tempdata['presence'] = False

    # handle salt
    if saltfound:
        saltdata['presence'] = True
        for layername, depthint in depthints.items():
            zind = np.where((depth >= depthint[0]) * (depth <= depthint[1]))[0]
            if len(zind) > 0:
                data = ncf.variables[vlib['salt']][tind, zind,loni,lati]
                saltdata[layername] = {'time': time[tind], 'value': data, 'depth_interval': depthint}
            else:
                saltdata[layername] = {'time': np.array([]), 'value': np.array([]), 'depth_interval': depthint}
    else:
        saltdata['presence'] = False

    # handle ssh
    if sshfound:
        sshdata['presence'] = True
        data = ncf.variables[vlib['ssh']][tind,loni,lati]
        sshdata['z0'] = {'time': time, 'value': data, 'depth_interval': [0, 0]}
    else:
        sshdata['presence'] = False


    sdata = {'longname': 'descriptive name of the station',
             'lon': lon,
             'lat': lat,
             'max_depth': maxz,
             'temp': tempdata,
             'salt': saltdata,
             'ssh': sshdata
             }
    return sdata
