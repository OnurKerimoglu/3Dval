# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 18:55:00 2017

@author: kerimoglu.o@gmail.com

"""
import os
import pickle
import numpy as np
import warnings
import datetime
from getm_funcs import get_getm_dom_vars,get_getm_dataF
from dcsm_funcs import get_dcsm_dataF,structure_dcsm_data
from general_funcs import interpval2D,get_2Dtree,getproj

def readsim(paths,simname,readraw,simdomain,meth2D,statsets,timeint,depthints,obs,vars,getmv,modtype):
    print('Reading simulation:'+simname)

    simf=paths[simname]
    simfb=os.path.basename(simf)
    simfp=os.path.dirname(simf)
    picklecode = '_%s_%s_%s-%s' % ('-'.join(statsets), '-'.join(depthints.keys()), timeint[0].year, timeint[1].year)
    pickledsimfile = os.path.join(simfp, simfb.split('.nc')[0] + picklecode + '.pickle')

    # check if the pickledsim exist
    if (not readraw) and os.path.exists(pickledsimfile):
        print('Opening pickled sim file')
        (sim,) = np.load(pickledsimfile,allow_pickle=True)
        return sim

    # if pickledsim does not exist:
    #collect all info that's needed to extract data from each station in the loop
    if simname[0:2] == 'GF' or simname[0:3] == 'SNS':
        print('Accessing getm data')
        lons,lats,bat,ysl,xsl=get_getm_dom_vars(simf,simdomain)
        if meth2D == 'pretree':
            proj = getproj(setup = 'SNSfull', projpath = os.path.dirname(os.path.realpath(__file__)))
            # proj = getproj(setup = 'NS', projpath = os.path.dirname(os.path.realpath(__file__)))
            domaintree = get_2Dtree(lons,lats,proj)
        else:
            raise (Exception('unknown spatial method for extracting values from GETM'))
        simdata,simtime = get_getm_dataF(simf,vars,ysl,xsl,getmv=getmv,modtype=modtype)
    elif simname[0:4] == 'DCSM':
        print('Accessing DCSM data')
        simdata, simtime, lons, lats, maxz_sim, StInds = get_dcsm_dataF(simf, vars)

    #fill the data in correct structure
    sim = {}
    stations=obs.keys()
    for station in stations:
        print('  ' + station)
        lon = obs[station]['lon']
        lat = obs[station]['lat']
        maxz_obs= obs[station]['bottom_depth']
        if (simname[0:2] == 'GF' or simname[0:3] == 'SNS') and (meth2D == 'pretree'):
            sdata = interp_simdata_on_station(station,simdata,simtime,proj,domaintree,bat,lon,lat,maxz_obs,timeint,depthints,vars)
        elif simname[0:6] == 'DCSM':
            #sdata = structure_dcsm_data(station,simf,lon,lat,maxz_obs,timeint,depthints,vars)
            sdata = structure_dcsm_data(station,lon,lat,maxz_obs,timeint,depthints,vars,simdata,simtime,lons,lats,maxz_sim,StInds)
        else:
            raise (Exception('unknown spatial method for extracting values from simulation'))
        sim[station] = sdata

    # pickle the sim file
    f = open(pickledsimfile, 'wb')
    pickle.dump((sim,), f)  # ,protocol=-1
    print('Pickled sim file for later use:'+pickledsimfile)
    f.close()

    return sim

def interp_simdata_on_station(station,simdata,time,proj,domaintree,bat,lon,lat,maxz_obs,timeint,depthints,vars,quickzfind=True):

    vardims={'ssh':'2D','temp':'3D','salt':'3D','DO':'3D','DOs':'3D','DIN':'3D','DIP':'3D','Si':'3D','NO3':'3D','NH4':'NH4','Chl':'3D',
             'Cyanobacteria':'3D','Diatoms':'3D','Dinoflagellates':'3D','Flagellates':'3D','Phaeocystis':'3D','other':'3D'}
    # maybe no need to check if other conditions are not satisfied
    XY_in=False
    z_in=True #assume z_in is ok by default
    t_in=False

    #first check if the data is available at all in the sim file
    varfound={}
    for varn in vars:
        varfound[varn] = True if varn in simdata.keys() else False

    # get zmax, see if it's a finite value (np.nan) means it's outside the domain, or interpolation can't be done
    if any(varfound.values()):
        #maxz = interp_2d_tree(bat, domaintree, lon, lat)
        maxz = interpval2D(0, 0, bat, lat, lon, 'pretree', proj, domaintree)
        if not np.isnan(maxz): XY_in = True
        # If the diff between  maxz with maxz_obs too large (=?), throw a warning
        if (not np.isnan(maxz_obs)) and (abs(maxz - maxz_obs) > 10):
            warnings.warn('For station %s, abs(maxz(sim=%s)-maxz(obs=%s))>5' %(station, maxz, maxz_obs))
            #XY_in = False
    else:
        raise(Exception('Requested variables were not found in the simulation file'))

    # quick scan if relevant depths are definitely not available
    if XY_in:
        if 'z' in simdata:
            if 'surface' in depthints.keys():
                z_inS = False if np.nanmin(simdata['z'])>depthints['surface'][1] else True
            if'bottom' in depthints.keys():
                z_inB = False if np.nanmax(simdata['z'])<20 else True
            z_in=False if ('bottom' in depthints.keys() and (not z_inB)) or ('surface' in depthints.keys() and (not z_inS)) else True
            numz = simdata['z'].shape[1] #will be later used

    # check if dates are in
    if XY_in and z_in:
        tind = np.where((time >= timeint[0]) * (time <= timeint[1]))[0]
        if len(tind) > 0: t_in = True

    # whether the data is available relevant for the station, depthinterval and time interval
    for varn in vars:
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

    for varn in vars:
        vdata = {}
        if varfound[varn]:
            if vardims[varn] == '2D':
                vdata['presence'] = True
                # value from a single cell:
                # data = ncf.variables[simdata[varn][tind,loni,lati]
                # values interpolated from nearest 4 cells
                data = np.zeros(len(tind)) * np.nan
                for tii, ti in enumerate(tind):
                    data[tii]= interpval2D(0, 0, simdata[varn][ti, :, :], lat, lon, 'pretree', proj, domaintree)
                vdata['z0'] = {'time': time, 'value': data, 'depth_interval': [0, 0]}
            elif vardims[varn]=='3D': #handle 3-D (t,z,x,y) vars
                vdata['presence'] = True
                for layername, depthint in depthints.items():
                    if len(simdata[varn].shape)==3: #no z dimension
                        data = np.zeros(len(tind)) * np.nan
                        for tii, ti in enumerate(tind):
                            # value from a single cell
                            # data = ncf.variables[simdata[varn][tind, zind,loni,lati]
                            # values interpolated from nearest 4 cells
                            data[tii] = interpval2D(0, 0, simdata[varn][ti, :], lat, lon, 'pretree', proj, domaintree)
                        vdata[layername] = {'time': time[tind], 'value': data, 'depth_interval': depthint}
                    else:
                        if quickzfind:
                            if ('z' not in simdata) and simdata[varn].shape[1]==1:
                                if layername=='surface':
                                    zi=0 #assume that the only provided layer is the surface
                                elif layername=='bottom':
                                    Warning('Skipping the bottom layer, as only 1 layer is found in the data set which is assumed to be surface')
                                    continue
                            elif simdata['z'].shape[1]==2:
                                if layername=='surface':
                                    zi=1
                                elif layername=='bottom':
                                    zi=0
                                else:
                                    #raise(Exception('simdata has 2 z levels and quickzfind was requested, but the layer to search (%s) is neither surface notr bottom'%layername))
                                    # if other layers are requested, recurse, but disable the quickzfind
                                    interp_simdata_on_station(station, simdata, time, proj, domaintree, bat, lon, lat,
                                                              maxz_obs,timeint, depthints, vars, quickzfind=False)
                            else:
                                #if more levels are available, recurse, but disable the quickzfind
                                interp_simdata_on_station(station, simdata, time, proj, domaintree, bat, lon, lat, maxz_obs,
                                                          timeint, depthints, vars, quickzfind=False)
                            data = np.zeros(len(tind)) * np.nan
                            for tii, ti in enumerate(tind):
                                # value from a single cell
                                # data = ncf.variables[simdata[varn][tind, zind,loni,lati]
                                # values interpolated from nearest 4 cells
                                data[tii] = interpval2D(0, 0, simdata[varn][ti, zi, :], lat, lon, 'pretree', proj, domaintree)
                            vdata[layername] = {'time': time[tind], 'value': data, 'depth_interval': depthint}
                        else: #exact search
                            # for all layers in sim file, calculate average depth and find the one in the correct interval
                            #doesn't work somehow?
                            zsimstat=np.zeros(numz)
                            for zi in range(numz):
                                #zsimstat[zi]=interp_2d_tree(simdata['z'][0,zi,:,:],domaintree,lon,lat,k=4)
                                zsimstat[zi] = interpval2D(0, 0, simdata['z'][0,zi,:,:], lat, lon, 'pretree', proj, domaintree)
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
                                        zsimstat[zi] = interpval2D(0, 0, simdata['z'][ti, zi, :, :], lat, lon, 'pretree', proj, domaintree)
                                    data[tii]=np.mean(vsimstat)
                                vdata[layername] = {'time': time[tind], 'value': data, 'depth_interval': depthint}
                            else:
                                vdata[layername] = {'time': np.array([]), 'value': np.array([]), 'depth_interval': depthint}
        else:
            vdata['presence'] = False

        sdata[varn]=vdata

    return sdata
