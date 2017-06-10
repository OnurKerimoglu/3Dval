"""
Created on 10 June 2017
@authors: onur.kerimoglu@hzg.de
provides functions relevant to getm simulations
"""

import netCDF4
import numpy as np

def get_getm_dataF(simf,varns,ysl,xsl):
    vlib = {'t': 'time', 'z': 'depth', 'temp': 'temp', 'salt': 'salt', 'ssh': 'elev'}
    ncf = netCDF4.Dataset(simf)
    simdata={}
    #add depth to the varlist
    varns.append('z')
    for varn in varns:
        if vlib[varn] in ncf.variables:
            varF=ncf.variables[vlib[varn]][:]
            if len(varF.shape)==2:
                var = varF[ysl,xsl]
            elif len(varF.shape)==3:
                var = varF[:,ysl, xsl]
            elif len(varF.shape)==4:
                var = varF[:,:, ysl, xsl]
            if np.ma.is_masked(var):
                var[var.mask] = np.nan  # transform the masked values to nan, such that intp will result in nan if any cell is nan
            simdata[varn] = var
        #multiply depth with -1
        if varn=='z':
            simdata[varn]=-1*simdata[varn]
    #add time
    time_num = ncf.variables[vlib['t']][:]
    simtime = netCDF4.num2date(time_num, ncf.variables[vlib['t']].getncattr('units'))
    ncf.close()
    return (simdata,simtime)

def get_getm_dom_vars(simdomain):

    # read getm- topo file for coordinates and bathymetry
    topo = get_getm_bathymetry_cropped()

    if simdomain == 'SNSe':
        ysl = slice(10, 75)
        xsl = slice(60, 137)
    else:
        ysl = slice(0, topo['lonc'].shape[0])
        xsl = slice(0, topo['lonc'].shape[1])

    lons=topo['lonc'][ysl, xsl]
    lats=topo['latc'][ysl, xsl]
    bat=topo['H'][ysl, xsl]
    if np.ma.is_masked(bat):
        bat[bat.mask] = np.nan  # transform the masked values to nan

    return (lons,lats,bat,ysl,xsl)

def get_getm_bathymetry_cropped(fname='/home/onur/WORK/projects/GB/data/topo/topo_area_sns.nc'):
    ncB=netCDF4.Dataset(fname)
    ncBv=ncB.variables
    #bathymetry from topo_sns.nc.
    lonx=ncBv['lonx'][4:-1,1:-1] #this should be [95,138]
    latx=ncBv['latx'][4:-1,1:-1] #this should be [95,138]
    lonc=0.25*(lonx[:-1,:-1]+lonx[:-1,1:]+lonx[1:,:-1]+lonx[1:,1:]) #this should be [94,137]
    latc=0.25*(latx[:-1,:-1]+latx[:-1,1:]+latx[1:,:-1]+latx[1:,1:])
    H = ncBv['bathymetry'][4:-1,1:-1] #this should be [94,137])
    A= ncBv['A'][3:-1,1:-2] #this should be [94,137])
    topo={'H':H,'latc':latc, 'lonc':lonc,'latx':latx, 'lonx':lonx,'Hunit':ncBv['bathymetry'].units,'A':A}
    ncB.close()
    return(topo)