"""
Created on 10 June 2017
@authors: onur.kerimoglu@hzg.de
provides general functions useful for multiple modules
"""

import os,sys
import numpy as np
import datetime
import netCDF4
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from scipy import interpolate
#import time

def cm2inch(*tupl):
    inch = 2.54
    if type(tupl[0]) == tuple:
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)

def date_matchup(obs0,sim0):
    # get rid of invalid entries
    ovali=np.isfinite(obs0['value'])
    svali=np.isfinite(sim0['value'])
    obs= {'value':obs0['value'][ovali], 'time':obs0['time'][ovali]}
    sim = {'value': sim0['value'][svali], 'time': sim0['time'][svali]}
    #datetime to date
    obsd = [dt.date() for dt in obs['time']]
    simd = [dt.date() for dt in sim['time']]
    # find the common dates
    comd = np.sort(list(set(obsd).intersection(simd)))
    # find the indices
    oi = [obsd.index(d) for d in comd]
    si = [simd.index(d) for d in comd]
    #return the (unique) dates, and corresponding obs & sim pairs as a dict
    return({'dates':obs['time'][oi],'ref':obs['value'][oi], 'model':sim['value'][si]})

def discrete_cmap(N, base_cmap=None):
    # By Jake VanderPlas
    # License: BSD-style
    """Create an N-bin discrete colormap from the specified input map"""

    # Note that if base_cmap is a string or None, you can simply do
    #    return plt.cm.get_cmap(base_cmap, N)
    # The following works for string, None, or a colormap instance:

    base = plt.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, N))
    cmap_name = base.name + str(N)
    try:
        #standard colormaps, eg., jet
        new_cmap=base.from_list(cmap_name, color_list, N)
    except:
        #eg., for extra colormaps implemented as a ListedColormap like viridis
        new_cmap=mpl.colors.LinearSegmentedColormap.from_list(cmap_name, color_list, N)
    return new_cmap

def fnclean(fname):
    import collections
    #print 'cleaning fname.\nin:'+fname
    extlib=['.png','.pdf']
    substdict=collections.OrderedDict()
    substdict['nc']=''
    substdict['.'] = ''
    substdict[','] = '_'

    #assume that there's no extension
    bn=fname
    ext=''
    #strip the extension away if detected
    for i,exti in enumerate(extlib):
        if exti in fname:
            ext=exti
            bn=fname.split(ext)[0]
            break

    #substitute all characters in the dict
    for rem,add in substdict.items():
        bn=bn.replace(rem,add)

    fname=bn+ext
    #print 'out:'+fname
    return fname

def discrete_cmap_tuner(clim,vallims,Nlev,colmap,nancol='white'):

    cbt = np.linspace(clim[0], clim[-1], Nlev)
    cbtstep = cbt[-1] - cbt[0]
    intbounds = list(cbt)

    cmap = discrete_cmap(Nlev - 1, colmap)

    extover = False;
    extunder = False
    allbounds = intbounds
    if vallims[0] < clim[0]:  # and clim[0]!=0:
        extunder = True
        allbounds = [allbounds[0] - cbtstep] + allbounds
        # Nlev = Nlev + 1
        if colmap in ['Blues_r', 'bwr', 'Reds_r']:
            cmap.set_under((0., 0., 0.))  # black
        elif colmap == 'Reds':
            cmap.set_under((0.5, 0.5, 1.0))  # light blue
        elif colmap == 'viridis':
            cmap.set_under((0.0, 0.0, 0.0))  # black
        elif colmap == 'viridis_r':
            cmap.set_under((1.0, 0.8, 0.2))  # light orange
    if vallims[1] > clim[-1] * 1.05:  # and clim[1]!=0:
        extover = True
        allbounds = allbounds + [allbounds[-1] + cbtstep]
        # Nlev = Nlev + 1
        if colmap == 'Blues_r':
            cmap.set_over((1.0, 0.8, 0.8))  # light red
        elif colmap == 'Reds_r':
            cmap.set_over((0.8, 0.8, 1.0))  # light blue
        elif colmap in ['Reds', 'bwr']:
            cmap.set_over((0.0, 0.0, 0.0))  # black
        elif colmap == 'viridis':
            cmap.set_over((1.0, 0.8, 0.2))  # light orange
        elif colmap == 'viridis_r':
            cmap.set_over((0.0, 0.0, 0.0))  # black

    cmap.set_bad(str(nancol))

    if (extunder==False and extover==False):
        extopt='neither'
    elif extunder==False and extover==True:
        extopt='max'
    elif extunder==True and extover==False:
        extopt='min'
    else:
        extopt = 'both'

    return(cmap,extopt,cbt,intbounds,allbounds)

def format_date_axis(ax,tspan):
    #ax.set_xlim(datetime.date(tspan[0].year,1,1), datetime.date(tspan[1].year,tspan[1].12,31))
    ax.set_xlim(datetime.date(tspan[0].year,1,1), tspan[1])
    if np.diff(tspan)[0].days<63:
        ax.xaxis.set_major_locator(mpl.dates.WeekdayLocator(byweekday=mpl.dates.MO) )
        ax.xaxis.set_major_formatter(mpl.dates.DateFormatter('%b\n%d'))
    elif np.diff(tspan)[0].days<367:
        ax.xaxis.set_major_locator(mpl.dates.MonthLocator(bymonthday=1, interval=3) )
        ax.xaxis.set_major_formatter(mpl.dates.DateFormatter('%b'))
        ax.xaxis.set_minor_locator(mpl.dates.MonthLocator(bymonthday=1, interval=1))
    elif np.diff(tspan)[0].days<732:
        ax.xaxis.set_major_locator(mpl.dates.MonthLocator(bymonth=[1],interval=1))
        ax.xaxis.set_major_formatter(mpl.dates.DateFormatter('%Y')) #%m %y
        ax.xaxis.set_minor_locator(mpl.dates.MonthLocator(bymonthday=[1], interval=1))
        ax.xaxis.set_minor_formatter(mpl.dates.DateFormatter(''))
    elif np.diff(tspan)[0].days<1466:
        ax.xaxis.set_major_locator(mpl.dates.MonthLocator(bymonth=[1], interval=1) )
        ax.xaxis.set_major_formatter(mpl.dates.DateFormatter('%Y'))
        ax.xaxis.set_minor_locator(mpl.dates.MonthLocator(bymonthday=1, interval=1))
        ax.xaxis.set_minor_formatter(mpl.dates.DateFormatter(''))
    elif np.diff(tspan)[0].days<3655:
        ax.xaxis.set_major_locator(mpl.dates.MonthLocator(bymonth=1,bymonthday=1, interval=1) )
        ax.xaxis.set_minor_locator(mpl.dates.MonthLocator(bymonth=7,bymonthday=1, interval=1))
        ax.xaxis.set_minor_formatter(mpl.dates.DateFormatter('%Y'))
    elif np.diff(tspan)[0].days<9130: #25*365=9125
        ax.xaxis.set_major_locator(mpl.dates.MonthLocator(bymonth=1,bymonthday=1, interval=1) )
        ax.xaxis.set_minor_locator(mpl.dates.MonthLocator(bymonth=7,bymonthday=1, interval=1))
        ax.xaxis.set_minor_formatter(mpl.dates.DateFormatter('%Y'))
    else:
        ax.xaxis.set_major_locator(mpl.dates.MonthLocator(bymonthday=1, interval=1))
        ax.xaxis.set_minor_locator(mpl.dates.MonthLocator(bymonthday=1, interval=1))
        ax.xaxis.set_minor_formatter(mpl.dates.DateFormatter('%y'))

    if np.diff(tspan)[0].days<367:
        ax.tick_params(axis='x', which='major', direction='out',labelsize=9)
        ax.tick_params(axis='x', which='minor', direction='out')
    elif np.diff(tspan)[0].days< 1466:
        ax.tick_params(axis='x', which='major', length=2, direction='out', labelsize=9, pad=1)
        ax.tick_params(axis='x', which='minor', length=2, direction='in', labelsize=9, pad=1)  # labelbottom='off')
        #ax.grid(axis='x', which='major', color='0.5', linestyle='-', linewidth=1.5)
        ax.grid(axis='x', which='major', color='0.5', linestyle='-', linewidth=.5)
    else:
        ax.tick_params(axis='x', which='major',direction='out',labelbottom='off')
        ax.tick_params(axis='x', which='minor', length=0, labelsize=9)
        ax.grid(axis='x', which='major', color='.5', linestyle='-', linewidth=1)
        #ax.grid(axis='x', which='minor', color='k', linestyle='-', linewidth=.5)

def get_2Dtree(lons,lats,proj,Vpy=3):
    if len(lons.shape) == 2:
        lonsM, latsM = lons, lats
    else:
        lonsM, latsM = np.meshgrid(lons, lats)
    x1, y1 = proj(lonsM, latsM)
    if Vpy==2:
        xypairs = zip(x1.flat, y1.flat)
    else:
        xypairs = list(zip(x1.flat, y1.flat))  # p3: must be placed in list

    domaintree = cKDTree(xypairs)
    return domaintree

def get_botdepth(lon,lat,method='tree'):
    proj = getproj(setup='SNSfull', projpath=os.path.dirname(os.path.realpath(__file__)))
    # proj = getproj(setup = 'NS', projpath = os.path.dirname(os.path.realpath(__file__)))
    topofile = '/home/onur/WORK/projects/GB/data/topo/topo_HR2.nc'
    #topofile='/home/onur/WORK/projects/GB/data/topo/topo_GreaterNorthSea.nc'
    if os.path.exists(topofile):
        nc=netCDF4.Dataset(topofile)
        # depths=nc.variables['depth_average'][:]
        # lats = nc.variables['lat'][:]
        # lons = nc.variables['lon'][:]
        depths = nc.variables['bathymetry'][:]
        lats = nc.variables['latx'][:-1,:-1]
        lons = nc.variables['lonx'][:-1,:-1]
        nc.close()
        depth = interpval2D(lats, lons, depths, lat, lon, method,proj)
    else:
        warnings.warn('Topo file not found in the provided location (%s). Filling nan for station depth:'%topofile)
        depth=np.nan
    return depth

def get_skills_atstations(obs,sim,vars,layer):
    #do the matchups and combine the obs&sim in a single structure:
    #V[var] = { 'ref', 'model'} and optionally {'dates', 'lats', 'lons',}
    V={}
    for var in vars:
        if obs[var]['presence']:
            V[var]=date_matchup(obs[var][layer],sim[var][layer])
        else:
            print ('%s not found'%var)
            V[var]={'dates':[np.nan],'ref':[np.nan],'model':[np.nan]}
    stats=get_stats(V)
    return(stats)

def get_stats(V):
    #copied from 'do_3dstats'
    #expected structure: V[var]={'dates','lats','lons','ref','model'}
    Vstats={}
    for v in V.keys():
        ref=V[v]['ref']
        model=V[v]['model']
        if len(ref)!=1:
            stats = {}
            stats['Rstd'] = np.std(ref)
            stats['Mstd'] = np.std(model)
            stats['n'] = len(ref)
            stats['rho']=np.corrcoef([ref,model])[0,1]
            if np.mean(ref)>-0.0001 and np.mean(ref)<0.0001: #eg, because the values are anomolies
                stats['B']=np.nan
                stats['B*']=np.nan
            else:
                stats['B'] = np.mean(model) - np.mean(ref)
                stats['B*'] = (np.mean(model) - np.mean(ref))/(np.mean(ref))
            #todo: other statistics?
        else:
            stats={'Rstd':np.nan,'Mstd':np.nan,'n':0, 'rho':np.nan, 'B':np.nan,'B*':np.nan}
        Vstats[v]=stats
    return Vstats

def getproj(setup,projpath=os.path.dirname(os.path.realpath(__file__))):
    pickledproj=os.path.join(projpath,'proj_'+setup+'.pickle')
    if os.path.isfile(pickledproj):
        print ('opening an existing projection: '+ pickledproj)
        #if a projection exists, just load it (fast)
        (proj,) = np.load(pickledproj, encoding='latin1')
    else:
        from mpl_toolkits.basemap import Basemap
        import pickle
        print ('projecting for: '+ setup)
        if setup=='SNSfull':
            proj=Basemap(projection='lcc',
                       resolution='i',
                       llcrnrlon=0.0,
                       llcrnrlat=51.0, #51.2,
                       urcrnrlon=9.5,
                       urcrnrlat=56.0, #55.8,
                       lat_0=52.0,
                       lon_0=5.)
        elif setup=='WadSea':
            proj = Basemap(projection='lcc',
                           resolution='i',
                           llcrnrlon=6.2,
                           llcrnrlat=53.2,
                           urcrnrlon=9.3,
                           urcrnrlat=55.1,
                           lat_0=52.0,
                           lon_0=5.)
        elif setup=='Ems':
            proj = Basemap(projection='lcc',
                           resolution='f',
                           llcrnrlon=6.0296,
                           llcrnrlat=52.8525,
                           urcrnrlon=7.7006,
                           urcrnrlat=53.9622,
                           lat_0=52.0,
                           lon_0=5.)
        elif setup=='GBight':
            proj = Basemap(projection='lcc',
                           resolution='i',
                           llcrnrlon=4,
                           llcrnrlat=52.9,
                           urcrnrlon=9.7,
                           urcrnrlat=55.6,
                           lat_0=52.0,
                           lon_0=5.)
        elif setup=='SNSext':
            proj=Basemap(projection='lcc',
                       resolution='i',
                       llcrnrlon=0.0,
                       llcrnrlat=51.0, #51.2,
                       urcrnrlon=9.5,
                       urcrnrlat=58.0, #55.8,
                       lat_0=52.0,
                       lon_0=5.)
        elif setup=='SENS':
            proj=Basemap(projection='lcc',
                       resolution='i',
                       llcrnrlon=4.3,
                       llcrnrlat=52.4, #51.2,
                       urcrnrlon=9.3,
                       urcrnrlat=56, #55.8,
                       lat_0=52.0,
                       lon_0=5.)
        elif setup == 'WadSea':
            proj = Basemap(projection='lcc',
                           resolution='i',
                           llcrnrlon=6.2,
                           llcrnrlat=53.2,  # 51
                           urcrnrlon=9.3,
                           urcrnrlat=55.1,
                           lat_0=52.0,
                           lon_0=5.)
        elif setup=='CuxBusHel':
            proj=Basemap(projection='lcc',
                       resolution='i',
                       llcrnrlon=7.7,
                       llcrnrlat=53.5, #51
                       urcrnrlon=9.2,
                       urcrnrlat=54.3,
                       lat_0=52.0,
                       lon_0=5.)
        else:
            raise (Exception('unknown setup:' + setup + '. cannot produce projection.'))

        #pickle for later use:
        f=open(pickledproj,'wb')
        pickle.dump((proj,),f) #,protocol=-1
        f.close()
    return proj

def interp_2d_tree(vals,domaintree,lon,lat,k=4,kmin=1):
    #vals: 2d field of values to be interpolated
    #domaintree: query tree generated with scipy.spatial.cKDTree
    #lon,lat:  coordinates to be used for interpolation
    #k: number of cells to use for interpolation

    dLi, gridindsLi = domaintree.query((lon, lat), k)
    wLi = 1.0 / dLi ** 2

    #remove the mask, if necessary
    if isinstance(vals,np.ma.MaskedArray):
        #vals[vals.mask] = np.nan
        #vals.Mask=False
        vals = vals.filled(np.nan)

    if len(vals.shape)==2:
        vals2int=vals[:, :].flatten()[gridindsLi]
        #find the nan-indices and remove them (by nullifying their weights)
        nani=np.where(np.isnan(vals2int))
        wLi[nani]=0
        vals2int[nani]=0
        intval = np.sum(wLi * vals2int) / np.sum(wLi)
    else:
        raise(ValueError('shape of the data to be interpolated is more than 2-dimensional'))
    return intval


def interpval2D(lats,lons,vals,lat,lon,method,proj,domaintree=0):

    # convert lat-lons to cartesian coordinates
    x2,y2 = proj(lon,lat)

    if method=='pretree':
        if not isinstance(domaintree,int):  #p2: istype(domaintree,'integer'):
            domaintree=domaintree
        else:
            warnings.WarningMessage('no pre-generated tree provided. Generating a new one')
            method='tree'
    else:
        # meshscipy.spatial.ckdtree.cKDTree.query
        if len(lons.shape) == 2:
            lonsM, latsM = lons, lats
        else:
            lonsM, latsM = np.meshgrid(lons, lats)
        x1, y1 = proj(lonsM, latsM)

    if method == 'tree':
        xypairs = list(zip(x1.flat, y1.flat)) #p3: must be placed in list
        # ts=time.time()
        domaintree = cKDTree(xypairs)
        # print 'tree generated in:%.1f' %(time.time() - ts)

    if method in ['NN','linear']:
        coords = np.asarray(zip(np.reshape(x1,x1.size),np.reshape(y1,y1.size)))
        if method=='NN':
            f = interpolate.NearestNDInterpolator(coords, np.reshape(vals, vals.size))
        elif method=='linear':
            f = interpolate.LinearNDInterpolator(coords, np.reshape(vals,vals.size), fill_value=np.nan)
        val = f(x2, y2)
    elif method in ['tree','pretree']:
        val=interp_2d_tree(vals, domaintree, x2, y2, k=4, kmin=1)

    return val
