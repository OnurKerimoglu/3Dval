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

from getm_funcs import get_getm_dom_vars,get_getm_dataF
from general_funcs import interp_2d_tree,get_2Dtree_p3

def readsim(paths,simname,readraw,simdomain,meth2D,statsets,timeint,depthints,obs):
    print('Reading simulation:'+simname)

    simf=paths[simname]
    simfb=os.path.basename(simf)
    simfp=os.path.dirname(simf)
    picklecode = '_%s_%s_%s-%s' % ('-'.join(statsets), '-'.join(depthints.keys()), timeint[0].year, timeint[1].year)
    pickledsimfile = os.path.join(simfp, simfb.split('.nc')[0] + picklecode + '.pickle')

    # check if the pickledsim exist
    if (not readraw) and os.path.exists(pickledsimfile):
        print('Opening pickled sim file')
        (sim,) = np.load(pickledsimfile)
        return sim

    # if pickledsim does not exist:
    #collect all info that's needed to extract data from each station in the loop
    if simname[0:4] == 'GETM':
        print('Accessing getm data')
        lons,lats,bat,ysl,xsl=get_getm_dom_vars(simdomain)
        if meth2D == 'int_tree':
            domaintree = get_2Dtree_p3(lons,lats)
        else:
            raise (Exception('unknown spatial method for extracting values from GETM'))
        simdata,simtime = get_getm_dataF(simf, ['temp','salt','ssh'],ysl,xsl)

    #fill the data in correct structure
    sim = {}
    stations=obs.keys()
    for station in stations:
        print('  ' + station)
        lon = obs[station]['lon']
        lat = obs[station]['lat']
        maxz_obs= obs[station]['bottom_depth']
        if (simname[0:4] == 'GETM') and (meth2D == 'int_tree'):
            sdata = get_station_data_getm_inttree(station,simdata,simtime,domaintree,bat,lon,lat,timeint,depthints,maxz_obs)
        elif simname[0:5] == 'FVCOM':
            #sdata =get_station_data_fvcom() #TODO:IVAN
            raise(Exception('not yet implemented'))
        else:
            raise (Exception('unknown spatial method for extracting values from simulation'))
        sim[station] = sdata

    # pickle the sim file
    f = open(pickledsimfile, 'wb')
    pickle.dump((sim,), f)  # ,protocol=-1
    print('Pickled sim file for later use:'+pickledsimfile)
    f.close()

    return sim

def get_station_data_getm_inttree(station,simdata,time,domaintree,bat,lon,lat,timeint,depthints,maxz_obs):

    varns={'temp':'3D','salt':'3D','ssh':'2D'}
    vlib = {'t': 'time', 'z': 'depth', 'temp': 'temp', 'salt': 'salt', 'ssh': 'elev'}

    # maybe no need to check if other conditions are not satisfied
    XY_in=False
    z_in=True #assume z_in is ok by default
    t_in=False

    #first check if the data is available at all in the sim file
    varfound={}
    for varn in varns.keys():
        varfound[varn] = True if varn in simdata.keys() else False

    # get zmax, see if it's a finite value (np.nan) means it's outside the domain, or interpolation can't be done
    if any(varfound.values()):
        maxz = interp_2d_tree(bat, domaintree, lon, lat)
        if not np.isnan(maxz): XY_in = True
        # If the diff between  maxz with maxz_obs too large (=?), throw a warning and skip the station
        if abs(maxz - maxz_obs) > 5:
            raise (Warning('For station %s, abs(maxz(sim=%s)-maxz(obs=%s))>5' %(station, maxz, maxz_obs)))
            XY_in = False

    # quick scan if relevant depths are definitely not available
    if XY_in:
        if 'z' in simdata:
            if 'surface' in depthints.keys():
                if np.nanmin(simdata['z'])>depthints['surface'][1]:
                    z_in = False
            elif 'bottom' in dpethints.keys():
                if np.nanmax(simdata['z'])<20:
                    z_in = False
            numz = simdata['z'].shape[1] #will be later used

    # check if dates are in
    if XY_in and z_in:
        tind = np.where((time >= timeint[0]) * (time <= timeint[1]))[0]
        if len(tind) > 0: t_in = True

    # whether the data is available relevant for the station, depthinterval and time interval
    for varn in varns.keys():
        varfound[varn] = True if varfound[varn] and t_in else False #t_in implies XY_in=z_in=True

    # update the bottom depth interval if necessary
    if 'bottom' in depthints.keys():
        depthints['bottom'] = [maxz - depthints['bottom'][0], maxz - depthints['bottom'][1]]

    #fill in the data:
    sdata = {  # 'longname': '',
        'lon': lon,
        'lat': lat,
        'bottom_depth': maxz,
    }

    for varn in varns.keys():
        vdata = {}
        if varfound[varn]:
            if varns[varn] == '2D':
                vdata['presence'] = True
                # value from a single cell:
                # data = ncf.variables[simdata[varn][tind,loni,lati]
                # values interpolated from nearest 4 cells
                data = interp_2d_tree(simdata[varn][tind, :, :], domaintree, lon, lat, k=4)
                vdata['z0'] = {'time': time, 'value': data, 'depth_interval': [0, 0]}
            elif varns[varn]=='3D': #handle 3-D (t,z,x,y) vars
                vdata['presence'] = True
                for layername, depthint in depthints.items():
                    #for all layers in sim file, calculate average depth and find the one in the correct interval
                    zsimstat=np.zeros(numz)
                    for zi in range(numz):
                        zsimstat[zi]=interp_2d_tree(simdata['z'][tind,zi,:,:],domaintree,lon,lat,k=4)
                    zind = np.where((zsimstat >= depthint[0]) * (zsimstat <= depthint[1]))[0]

                    if len(zind)>0:
                        #interpolate all values from the layers in the correct depth interval, and calculate the average
                        data=np.zeros(len(tind))*np.nan
                        for tii,ti in enumerate(tind):
                            vsimstat = np.zeros(len(zind))*np.nan
                            for zi in zind:
                                # value from a single cell
                                # data = ncf.variables[simdata[varn][tind, zind,loni,lati]
                                # values interpolated from nearest 4 cells
                                vsimstat[zi] = interp_2d_tree(simdata[varn][ti, zi, :, :], domaintree, lon, lat, k=4)
                            data[tii]=np.mean(vsimstat)
                        vdata[layername] = {'time': time[tind], 'value': data, 'depth_interval': depthint}
                    else:
                        vdata['presence'] = False
        else:
            vdata['presence'] = False

        sdata[varn]=vdata

    return sdata