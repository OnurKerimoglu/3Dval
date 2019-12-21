"""
Created on 10 June 2017
@authors: onur.kerimoglu@hzg.de
provides functions relevant to getm simulations
"""

import netCDF4
import numpy as np

def get_getm_dataF(simf,varns,ysl,xsl,getmv='mean',modtype='GF-PPZZ'):
    vlib = {'t': 'time', 'z': 'depth'}

    if getmv == '3d':
        vlibgetm = {'temp': 'temp', 'salt': 'salt', 'ssh': 'elev'}
    elif getmv == 'mean':
        vlibgetm={'temp': 'tempmean', 'salt': 'saltmean', 'ssh': 'elevmean'}
    vlib.update(vlibgetm)

    if modtype == 'GF-maecs':
        vlibfabm={'DOs': 'hzg_maecs_O2_percSat','DIN': 'hzg_maecs_nutN', 'DIP': 'hzg_maecs_nutP', 'Chl': 'hzg_maecs_chl'}
    elif modtype == 'GF-PZ':
        vlibfabm={'DOs':'EH_abioP_O2_percSat','DIN':'EH_abioP_DINO3+EH_abioP_DINH4','NO3':'EH_abioP_DINO3','NH4':'EH_abioP_DINH4','DIP':'EH_abioP_DIP','Si':'EH_abioP_DISi','Chl':'GPM_phy_Chl'}
    elif modtype=='GF-PPZZ':
        vlibfabm={'DOs':'EH_abioP_O2_percSat','DIN':'EH_abioP_DINO3+EH_abioP_DINH4','NO3':'EH_abioP_DINO3','NH4':'EH_abioP_DINH4','DIP':'EH_abioP_DIP','Si':'EH_abioP_DISi','Chl':'total_chlorophyll_calculator_result'}
    vlib.update(vlibfabm)

    try:
        ncf = netCDF4.Dataset(simf)
    except:
        raise(Exception('File not found:%s'%simf))
    simdata={}
    #add depth to the varlist
    #varns.append('z')
    varnsnew=varns+['z']
    for varn in varnsnew:
        #attempt to retrieve the variable
        varF,success = get_var_from_ncf(vlib[varn], ncf)
        if success:
            #varF=ncf.variables[vlib[varn]][:]
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

def get_var_from_ncf(varn_vl,ncf):
    #if '+' exists (eg., DIN=DINO3+DINH4); operate

    if '+' in varn_vl:
        varn_vl_list=varn_vl.split('+')
        varF=0
        success = True
        for varn_vl_i in varn_vl_list:
            varFi,success_i=get_var_from_ncf(varn_vl_i,ncf)
            varF=varF+varFi
            if not success_i: #if any of the additive variables can not be retrieved, assume failure
                success=False
    else:
        if varn_vl in ncf.variables:
            varF = ncf.variables[varn_vl][:]
            success=True
        else:
            varF=0
            success = False
            raise(Warning('requested variable not found:'+varn))

    return(varF,success)

def get_getm_dom_vars(simf,simdomain=''):
    dominfo_found=False
    #see if the domain info is provided in the simulation file
    print('attempting to extract domain info from:'+simf) 
    ncf=netCDF4.Dataset(simf)    
    ncv=ncf.variables.keys()
    if 'bathymetry' in ncv and 'lon' in ncv and 'lat' in ncv:
        topo={'lonc':ncf.variables['lon'][:], 
              'latc':ncf.variables['lat'][:], 
              'H':ncf.variables['bathymetry'][:]} 
        dominfo_found=True
    ncf.close()
    
    if not dominfo_found:
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
        bat=bat.filled(np.nan) # transform the masked values to nan
    return (lons,lats,bat,ysl,xsl)

def get_getm_bathymetry_cropped(fname='/home/onur/WORK/projects/GB/data/topo/topo_area_sns.nc',setup='SNS'):
    ncB=netCDF4.Dataset(fname)
    ncBv=ncB.variables
    #bathymetry from a topo file
    if setup=='SNS':
        lonx=ncBv['lonx'][4:-1,1:-1] #this should be [95,138]
        latx=ncBv['latx'][4:-1,1:-1] #this should be [95,138]
        H = ncBv['bathymetry'][4:-1,1:-1] #this should be [94,137])
    elif setup=='GB300':
        lonx = ncBv['lonx'][:, :]  # this should be [732,562]
        latx = ncBv['latx'][:, :]  # this should be [732,562]
        H = ncBv['bathymetry'][:, :]  # this should be [731,561])

    lonc = 0.25 * (lonx[:-1, :-1] + lonx[:-1, 1:] + lonx[1:, :-1] + lonx[1:, 1:])
    latc = 0.25 * (latx[:-1, :-1] + latx[:-1, 1:] + latx[1:, :-1] + latx[1:, 1:])
    topo={'H':H,'latc':latc, 'lonc':lonc,'latx':latx, 'lonx':lonx,'Hunit':ncBv['bathymetry'].units} #,'A':A}
    ncB.close()
    return(topo)
