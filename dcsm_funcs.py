"""
Created on 27 Jan 2021
@authors: onur.kerimoglu@uol.de
provides functions relevant to dcsm (deltares) simulations
"""

import netCDF4
import numpy as np
import datetime
from general_funcs import get_var_from_ncf
#import warnings

def get_dcsm_dataF(simf,vars,dcsmv='surface'):

    #if dcsmv == 'surface':
    vlib = {'t': 'time', 'zmax': 'TotalDepth',
            'temp': 'temperature', 'salt': 'salinity',
            'DO': 'OXY', 'Chl': 'Chlfa',
            'DIN': 'DIN','DIP': 'PO4', 'Si': 'Si',
            'Diatoms':'MDIATOMS','Flagellates':'MFLAGELA',
            'Phaeocystis':'PHAEOCYS','Dinoflagellates':'DINOFLAG',
            'SPM':'IM1'
            }

    try:
        ncf = netCDF4.Dataset(simf)
    except:
        raise (Exception('File not found:%s' % simf))
    simdata = {}

    for varn in vars:
        # attempt to retrieve the variable
        varF, success = get_var_from_ncf(vlib[varn], ncf)
        conv_factor=1.0
        if success:
            varF=ncf.variables[vlib[varn]][:]
            units=ncf.variables[vlib[varn]].units
            if len(varF.shape) == 2:
                var = varF[:, :]
            if np.ma.is_masked(var):
                var[var.mask] = np.nan  # transform the masked values to nan, such that intp will result in nan if any cell is nan
            if varn=='DIN' and units=='(gN/m3)':
                conv_factor=1000/14 #1000mg/g*1mmolN/14gN -> mmolN/m3
            elif varn=='DIP' and units=='(gP/m3)':
                conv_factor = 1000 / 31  # 1000mg/g*1mmolP/31gP -> mmolP/m3
            elif varn=='Si' and units=='(gSi/m3)':
                conv_factor = 1000 / 28  # 1000mg/g*1mmolP/28gP -> mmolSi/m3
            elif varn=='OXY' and units=='(g/m3)':
                #!! is it O (atomic weight 16) or O2 (molecular weight: 32)? assume O (16)
                conv_factor = 1000 / 16  # 1000mg/g*1mmolO/16gO -> mmolO2/m3
            elif units=='(gC/m3)' and varn in ['Flagellates','Phaeocystis','Diatoms','Dinoflagellates']:
                conv_factor = 1000/12
            simdata[varn] = var*conv_factor
            # multiply depth with -1
            if varn == 'z':
                simdata[varn] = -1 * simdata[varn]
    # add time
    time_num = ncf.variables[vlib['t']][:]
    #default netCDF4
    units=ncf.variables[vlib['t']].getncattr('units').replace('minutes','seconds')
    # simtime = netCDF4.num2date(time_num, units)
    #to use in combinaton with cftime lib:
    simtime = netCDF4.num2date(time_num, units,
                              only_use_cftime_datetimes=False,
                              only_use_python_datetimes=True)

    lons=ncf.variables['station_x_coordinate'][:]
    lats = ncf.variables['station_y_coordinate'][:]
    try:
        TotDepths=ncf.variables['TotalDepth'][:, :]
    except:
        TotDepths = ncf.variables['salinity'][:, :]*0
    
    if np.ndim(lons)>1:
        print('DCSM coordinates are given in 2D format! Eliminating time dimension')
        lonsxy=lons
        latsxy=lats
        lons=np.squeeze(lonsxy[0,:])
        lats=np.squeeze(latsxy[0,:])
    
    #construct a station index, i.e., the index that points to the respective station
    StInd = {}
    bla = [di for di, dim in enumerate(ncf.variables['station_name'].dimensions) if dim == 'stations']
    stationrange = range(0,ncf.variables['station_name'].shape[bla[0]])

    for stno in stationrange:
        if bla[0] == 0:
            ststr=ncf.variables['station_name'][stno, :].tostring()
        else:
            ststr = ncf.variables['station_name'][:,stno].tostring()
        stname=str(ststr.decode('ascii').strip().strip('\x00'))
        # if stname=='Nney_W_2' or stname=='NneyW21' or stname=='NneyW23':
        #     stname='Nney'
        # elif stname=='Bork_W_1' or stname=='BorkW1':
        #     stname='Bork'
        # elif stname=='JaBu_W_1':
            # stname='JaBu'
        # elif stname=='WeMu_W_1':
            #stname='WeMu'
            # stname='Wesermuendung2'
        # elif stname=='WeMu_W_2':
            #stname='WeMu'
            # stname='Wesermuendung'
        if stname=='220065':
            stname='GE_Norderelbe'
        elif stname=='220052_S':
            stname='GE_Sderpiep'
        elif stname=='220052_B':
            stname='GE_Bsum'
        elif stname=='220057':
            stname='GE_Hrnum_Vortrapptief'
        elif stname=='220006':
            stname='GE_Sdl_Amrum'
            
        StInd[stname]=stno
        StInd[stname.encode()]=stno #convert to byte class

    ncf.close()
    return (simdata,simtime,lons,lats,TotDepths,StInd)


def structure_dcsm_data(station,lon,lat,maxz,timeint,depthints,vars,simdata, time, lons, lats, TotDepths, StInds):
    # retrieve data from file
    #simdata, time, lons, lats, TotDepths, StInds = get_dcsm_dataF(simf, vars)

    vardims = {'ssh': '2D', 'temp': '3D', 'salt': '3D', 'DO': '3D', 'DOs': '3D', 'DIN': '3D', 'DIP': '3D', 'Si': '3D',
               'NO3': '3D', 'NH4': 'NH4', 'Chl': '3D','Diatoms':'3D','Dinoflagellates':'3D','Flagellates':'3D','Phaeocystis':'3D'}

    t_in = False
    st_in = False
    consistent=False

    # first check if the data is available at all in the sim file
    varfound = {}
    for varn in vars:
        varfound[varn] = True if varn in simdata.keys() else False

    # check if dates are in
    tind = np.where((time >= timeint[0]) * (time <= timeint[1]))[0]
    if len(tind) > 0: t_in = True

    #find the station index
    if station in StInds.keys():
        Sind = StInds[station]
        st_in = True
    else:
        print(station,StInds)
        
    #check whether the coordinates and max depth (roughly) match
    if t_in and st_in:
        if abs(lon - lons[Sind])>0.1:
            print('for station: %s inconsistent longitude. Expected: %s, Found: %s'%(station,lon,lons[Sind]))
        elif abs(lat - lats[Sind])>0.1:
            print('for station: %s inconsistent latitude. Expected: %s, Found: %s' % (station, lat, lats[Sind]))
        else:
            consistent=True

        maxz_sim=np.max(TotDepths[:,Sind])
        if abs(maxz - maxz_sim) > 5.0:
            #there are some differnces, but this may be due to interpolated bat, just report but don't mark as inconsistent
            print('for station: %s inconsistent max depth. Expected: %s, Found: %s' %(station,maxz,maxz_sim))

    # whether the data is available relevant for the station, time interval, and coordinates/maxdepths are consistent
    for varn in vars:
        varfound[varn] = True if varfound[varn] and st_in and t_in and consistent else False

    # fill in the data:
    sdata = {  # 'longname': '',
        'lon': lon,
        'lat': lat,
        'bottom_depth': maxz,
    }

    for varn in vars:
        vdata = {}
        if varfound[varn]:
            if vardims[varn] == '3D':  # handle 3-D (t,z,x,y) vars
                vdata['presence'] = True
                for layername, depthint in depthints.items():
                    if len(simdata[varn].shape) == 2:  # i.e., time and station, no z dimension
                        data = np.zeros(len(tind)) * np.nan
                        for tii, ti in enumerate(tind):
                            data[tii] = simdata[varn][ti, Sind]
                        vdata[layername] = {'time': time[tind], 'value': data, 'depth_interval': depthint}
        else:
            vdata['presence'] = False

        sdata[varn] = vdata

    return sdata
