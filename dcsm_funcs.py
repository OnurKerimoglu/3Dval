"""
Created on 27 Jan 2021
@authors: onur.kerimoglu@uol.de
provides functions relevant to dcsm (deltares) simulations
"""

import netCDF4
import numpy as np

from general_funcs import get_var_from_ncf

def get_dcsm_dataF(simf,vars,dcsmv='surface'):

    #if dcsmv == 'surface':
    vlib = {'t': 'time', 'zmax': 'TotalDepth',
            'temp': 'temperature', 'salt': 'salinity',
            'DO': 'OXY', 'Chl': 'Chlfa',
            'DIN': 'DIN','DIP': 'PO4', 'Si': 'Si'
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
            simdata[varn] = var*conv_factor
            # multiply depth with -1
            if varn == 'z':
                simdata[varn] = -1 * simdata[varn]
    # add time
    time_num = ncf.variables[vlib['t']][:]
    #default netCDF4
    units=ncf.variables[vlib['t']].getncattr('units').replace('minutes','seconds')
    simtime = netCDF4.num2date(time_num, units)
    #to use in combinaton with cftime lib:
    #simtime = netCDF4.num2date(time_num, units,
    #                           only_use_cftime_datetimes=False,
    #                           only_use_python_datetimes=True)

    lons=ncf.variables['station_x_coordinate'][:]
    lats = ncf.variables['station_y_coordinate'][:]
    TotDepths=ncf.variables['TotalDepth'][:,:]

    #construct a station index, i.e., the index that points to the respective station
    StInd = {}
    for stno in range(0,len(ncf.variables['station_name'])):
        ststr=ncf.variables['station_name'][stno, :].tostring()
        stname=str(ststr.decode('ascii').strip().strip('\x00'))
        StInd[stname]=stno
        StInd[stname.encode()]=stno #convert to byte class

    ncf.close()
    return (simdata,simtime,lons,lats,TotDepths,StInd)


def structure_dcsm_data(station,lon,lat,maxz,timeint,depthints,vars,simdata, time, lons, lats, TotDepths, StInds):
    # retrieve data from file
    #simdata, time, lons, lats, TotDepths, StInds = get_dcsm_dataF(simf, vars)

    vardims = {'ssh': '2D', 'temp': '3D', 'salt': '3D', 'DO': '3D', 'DOs': '3D', 'DIN': '3D', 'DIP': '3D', 'Si': '3D',
               'NO3': '3D', 'NH4': 'NH4', 'Chl': '3D'}

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