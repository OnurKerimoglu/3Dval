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

    if base_cmap in ['viridis','viridis_r']:
        reg_viridis()

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

def discrete_cmap_tuner(clim,vallims,Nlev,colmap,nancol='white'):

    cbt = np.linspace(clim[0], clim[1], Nlev)
    cbtstep = cbt[1] - cbt[0]
    intbounds = list(cbt)

    cmap = discrete_cmap(Nlev + 1, colmap)

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
    if vallims[1] > clim[1] * 1.05:  # and clim[1]!=0:
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
        ax.xaxis.set_major_locator(mpl.dates.MonthLocator(bymonth=[1,7],interval=1))
        ax.xaxis.set_major_formatter(mpl.dates.DateFormatter('%m/%y')) #%m %y
        ax.xaxis.set_minor_locator(mpl.dates.MonthLocator(bymonthday=1, interval=1))
        ax.xaxis.set_minor_formatter(mpl.dates.DateFormatter(''))
    elif np.diff(tspan)[0].days<1466:
        ax.xaxis.set_major_locator(mpl.dates.MonthLocator(bymonth=7, bymonthday=1, interval=1) )
        ax.xaxis.set_minor_locator(mpl.dates.MonthLocator(bymonth=7,bymonthday=1, interval=1))
        ax.xaxis.set_minor_formatter(mpl.dates.DateFormatter('%Y'))
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
    elif np.diff(tspan)[0].days< 732:
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

def getproj(setup,projpath):
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
        elif setup=='SNSext':
            proj=Basemap(projection='lcc',
                       resolution='i',
                       llcrnrlon=0.0,
                       llcrnrlat=51.0, #51.2,
                       urcrnrlon=9.5,
                       urcrnrlat=58.0, #55.8,
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

def reg_viridis():
    # Onur Kerimoglu: Following code is adapted from the mpl-colormaps package by Nathaniel Smith & Stefan van der Walt
    # https://github.com/BIDS/colormap.git
    #
    # New matplotlib colormaps by Nathaniel J. Smith, Stefan van der Walt,
    # and (in the case of viridis) Eric Firing.
    #
    # This file and the colormaps in it are released under the CC0 license /
    # public domain dedication. We would appreciate credit if you use or
    # redistribute these colormaps, but do not impose any legal restrictions.
    #
    # To the extent possible under law, the persons who associated CC0 with
    # mpl-colormaps have waived all copyright and related or neighboring rights
    # to mpl-colormaps.
    #
    # You should have received a copy of the CC0 legalcode along with this
    # work.  If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.

    __all__ = ['viridis']

    _viridis_data = [[0.267004, 0.004874, 0.329415],
                     [0.268510, 0.009605, 0.335427],
                     [0.269944, 0.014625, 0.341379],
                     [0.271305, 0.019942, 0.347269],
                     [0.272594, 0.025563, 0.353093],
                     [0.273809, 0.031497, 0.358853],
                     [0.274952, 0.037752, 0.364543],
                     [0.276022, 0.044167, 0.370164],
                     [0.277018, 0.050344, 0.375715],
                     [0.277941, 0.056324, 0.381191],
                     [0.278791, 0.062145, 0.386592],
                     [0.279566, 0.067836, 0.391917],
                     [0.280267, 0.073417, 0.397163],
                     [0.280894, 0.078907, 0.402329],
                     [0.281446, 0.084320, 0.407414],
                     [0.281924, 0.089666, 0.412415],
                     [0.282327, 0.094955, 0.417331],
                     [0.282656, 0.100196, 0.422160],
                     [0.282910, 0.105393, 0.426902],
                     [0.283091, 0.110553, 0.431554],
                     [0.283197, 0.115680, 0.436115],
                     [0.283229, 0.120777, 0.440584],
                     [0.283187, 0.125848, 0.444960],
                     [0.283072, 0.130895, 0.449241],
                     [0.282884, 0.135920, 0.453427],
                     [0.282623, 0.140926, 0.457517],
                     [0.282290, 0.145912, 0.461510],
                     [0.281887, 0.150881, 0.465405],
                     [0.281412, 0.155834, 0.469201],
                     [0.280868, 0.160771, 0.472899],
                     [0.280255, 0.165693, 0.476498],
                     [0.279574, 0.170599, 0.479997],
                     [0.278826, 0.175490, 0.483397],
                     [0.278012, 0.180367, 0.486697],
                     [0.277134, 0.185228, 0.489898],
                     [0.276194, 0.190074, 0.493001],
                     [0.275191, 0.194905, 0.496005],
                     [0.274128, 0.199721, 0.498911],
                     [0.273006, 0.204520, 0.501721],
                     [0.271828, 0.209303, 0.504434],
                     [0.270595, 0.214069, 0.507052],
                     [0.269308, 0.218818, 0.509577],
                     [0.267968, 0.223549, 0.512008],
                     [0.266580, 0.228262, 0.514349],
                     [0.265145, 0.232956, 0.516599],
                     [0.263663, 0.237631, 0.518762],
                     [0.262138, 0.242286, 0.520837],
                     [0.260571, 0.246922, 0.522828],
                     [0.258965, 0.251537, 0.524736],
                     [0.257322, 0.256130, 0.526563],
                     [0.255645, 0.260703, 0.528312],
                     [0.253935, 0.265254, 0.529983],
                     [0.252194, 0.269783, 0.531579],
                     [0.250425, 0.274290, 0.533103],
                     [0.248629, 0.278775, 0.534556],
                     [0.246811, 0.283237, 0.535941],
                     [0.244972, 0.287675, 0.537260],
                     [0.243113, 0.292092, 0.538516],
                     [0.241237, 0.296485, 0.539709],
                     [0.239346, 0.300855, 0.540844],
                     [0.237441, 0.305202, 0.541921],
                     [0.235526, 0.309527, 0.542944],
                     [0.233603, 0.313828, 0.543914],
                     [0.231674, 0.318106, 0.544834],
                     [0.229739, 0.322361, 0.545706],
                     [0.227802, 0.326594, 0.546532],
                     [0.225863, 0.330805, 0.547314],
                     [0.223925, 0.334994, 0.548053],
                     [0.221989, 0.339161, 0.548752],
                     [0.220057, 0.343307, 0.549413],
                     [0.218130, 0.347432, 0.550038],
                     [0.216210, 0.351535, 0.550627],
                     [0.214298, 0.355619, 0.551184],
                     [0.212395, 0.359683, 0.551710],
                     [0.210503, 0.363727, 0.552206],
                     [0.208623, 0.367752, 0.552675],
                     [0.206756, 0.371758, 0.553117],
                     [0.204903, 0.375746, 0.553533],
                     [0.203063, 0.379716, 0.553925],
                     [0.201239, 0.383670, 0.554294],
                     [0.199430, 0.387607, 0.554642],
                     [0.197636, 0.391528, 0.554969],
                     [0.195860, 0.395433, 0.555276],
                     [0.194100, 0.399323, 0.555565],
                     [0.192357, 0.403199, 0.555836],
                     [0.190631, 0.407061, 0.556089],
                     [0.188923, 0.410910, 0.556326],
                     [0.187231, 0.414746, 0.556547],
                     [0.185556, 0.418570, 0.556753],
                     [0.183898, 0.422383, 0.556944],
                     [0.182256, 0.426184, 0.557120],
                     [0.180629, 0.429975, 0.557282],
                     [0.179019, 0.433756, 0.557430],
                     [0.177423, 0.437527, 0.557565],
                     [0.175841, 0.441290, 0.557685],
                     [0.174274, 0.445044, 0.557792],
                     [0.172719, 0.448791, 0.557885],
                     [0.171176, 0.452530, 0.557965],
                     [0.169646, 0.456262, 0.558030],
                     [0.168126, 0.459988, 0.558082],
                     [0.166617, 0.463708, 0.558119],
                     [0.165117, 0.467423, 0.558141],
                     [0.163625, 0.471133, 0.558148],
                     [0.162142, 0.474838, 0.558140],
                     [0.160665, 0.478540, 0.558115],
                     [0.159194, 0.482237, 0.558073],
                     [0.157729, 0.485932, 0.558013],
                     [0.156270, 0.489624, 0.557936],
                     [0.154815, 0.493313, 0.557840],
                     [0.153364, 0.497000, 0.557724],
                     [0.151918, 0.500685, 0.557587],
                     [0.150476, 0.504369, 0.557430],
                     [0.149039, 0.508051, 0.557250],
                     [0.147607, 0.511733, 0.557049],
                     [0.146180, 0.515413, 0.556823],
                     [0.144759, 0.519093, 0.556572],
                     [0.143343, 0.522773, 0.556295],
                     [0.141935, 0.526453, 0.555991],
                     [0.140536, 0.530132, 0.555659],
                     [0.139147, 0.533812, 0.555298],
                     [0.137770, 0.537492, 0.554906],
                     [0.136408, 0.541173, 0.554483],
                     [0.135066, 0.544853, 0.554029],
                     [0.133743, 0.548535, 0.553541],
                     [0.132444, 0.552216, 0.553018],
                     [0.131172, 0.555899, 0.552459],
                     [0.129933, 0.559582, 0.551864],
                     [0.128729, 0.563265, 0.551229],
                     [0.127568, 0.566949, 0.550556],
                     [0.126453, 0.570633, 0.549841],
                     [0.125394, 0.574318, 0.549086],
                     [0.124395, 0.578002, 0.548287],
                     [0.123463, 0.581687, 0.547445],
                     [0.122606, 0.585371, 0.546557],
                     [0.121831, 0.589055, 0.545623],
                     [0.121148, 0.592739, 0.544641],
                     [0.120565, 0.596422, 0.543611],
                     [0.120092, 0.600104, 0.542530],
                     [0.119738, 0.603785, 0.541400],
                     [0.119512, 0.607464, 0.540218],
                     [0.119423, 0.611141, 0.538982],
                     [0.119483, 0.614817, 0.537692],
                     [0.119699, 0.618490, 0.536347],
                     [0.120081, 0.622161, 0.534946],
                     [0.120638, 0.625828, 0.533488],
                     [0.121380, 0.629492, 0.531973],
                     [0.122312, 0.633153, 0.530398],
                     [0.123444, 0.636809, 0.528763],
                     [0.124780, 0.640461, 0.527068],
                     [0.126326, 0.644107, 0.525311],
                     [0.128087, 0.647749, 0.523491],
                     [0.130067, 0.651384, 0.521608],
                     [0.132268, 0.655014, 0.519661],
                     [0.134692, 0.658636, 0.517649],
                     [0.137339, 0.662252, 0.515571],
                     [0.140210, 0.665859, 0.513427],
                     [0.143303, 0.669459, 0.511215],
                     [0.146616, 0.673050, 0.508936],
                     [0.150148, 0.676631, 0.506589],
                     [0.153894, 0.680203, 0.504172],
                     [0.157851, 0.683765, 0.501686],
                     [0.162016, 0.687316, 0.499129],
                     [0.166383, 0.690856, 0.496502],
                     [0.170948, 0.694384, 0.493803],
                     [0.175707, 0.697900, 0.491033],
                     [0.180653, 0.701402, 0.488189],
                     [0.185783, 0.704891, 0.485273],
                     [0.191090, 0.708366, 0.482284],
                     [0.196571, 0.711827, 0.479221],
                     [0.202219, 0.715272, 0.476084],
                     [0.208030, 0.718701, 0.472873],
                     [0.214000, 0.722114, 0.469588],
                     [0.220124, 0.725509, 0.466226],
                     [0.226397, 0.728888, 0.462789],
                     [0.232815, 0.732247, 0.459277],
                     [0.239374, 0.735588, 0.455688],
                     [0.246070, 0.738910, 0.452024],
                     [0.252899, 0.742211, 0.448284],
                     [0.259857, 0.745492, 0.444467],
                     [0.266941, 0.748751, 0.440573],
                     [0.274149, 0.751988, 0.436601],
                     [0.281477, 0.755203, 0.432552],
                     [0.288921, 0.758394, 0.428426],
                     [0.296479, 0.761561, 0.424223],
                     [0.304148, 0.764704, 0.419943],
                     [0.311925, 0.767822, 0.415586],
                     [0.319809, 0.770914, 0.411152],
                     [0.327796, 0.773980, 0.406640],
                     [0.335885, 0.777018, 0.402049],
                     [0.344074, 0.780029, 0.397381],
                     [0.352360, 0.783011, 0.392636],
                     [0.360741, 0.785964, 0.387814],
                     [0.369214, 0.788888, 0.382914],
                     [0.377779, 0.791781, 0.377939],
                     [0.386433, 0.794644, 0.372886],
                     [0.395174, 0.797475, 0.367757],
                     [0.404001, 0.800275, 0.362552],
                     [0.412913, 0.803041, 0.357269],
                     [0.421908, 0.805774, 0.351910],
                     [0.430983, 0.808473, 0.346476],
                     [0.440137, 0.811138, 0.340967],
                     [0.449368, 0.813768, 0.335384],
                     [0.458674, 0.816363, 0.329727],
                     [0.468053, 0.818921, 0.323998],
                     [0.477504, 0.821444, 0.318195],
                     [0.487026, 0.823929, 0.312321],
                     [0.496615, 0.826376, 0.306377],
                     [0.506271, 0.828786, 0.300362],
                     [0.515992, 0.831158, 0.294279],
                     [0.525776, 0.833491, 0.288127],
                     [0.535621, 0.835785, 0.281908],
                     [0.545524, 0.838039, 0.275626],
                     [0.555484, 0.840254, 0.269281],
                     [0.565498, 0.842430, 0.262877],
                     [0.575563, 0.844566, 0.256415],
                     [0.585678, 0.846661, 0.249897],
                     [0.595839, 0.848717, 0.243329],
                     [0.606045, 0.850733, 0.236712],
                     [0.616293, 0.852709, 0.230052],
                     [0.626579, 0.854645, 0.223353],
                     [0.636902, 0.856542, 0.216620],
                     [0.647257, 0.858400, 0.209861],
                     [0.657642, 0.860219, 0.203082],
                     [0.668054, 0.861999, 0.196293],
                     [0.678489, 0.863742, 0.189503],
                     [0.688944, 0.865448, 0.182725],
                     [0.699415, 0.867117, 0.175971],
                     [0.709898, 0.868751, 0.169257],
                     [0.720391, 0.870350, 0.162603],
                     [0.730889, 0.871916, 0.156029],
                     [0.741388, 0.873449, 0.149561],
                     [0.751884, 0.874951, 0.143228],
                     [0.762373, 0.876424, 0.137064],
                     [0.772852, 0.877868, 0.131109],
                     [0.783315, 0.879285, 0.125405],
                     [0.793760, 0.880678, 0.120005],
                     [0.804182, 0.882046, 0.114965],
                     [0.814576, 0.883393, 0.110347],
                     [0.824940, 0.884720, 0.106217],
                     [0.835270, 0.886029, 0.102646],
                     [0.845561, 0.887322, 0.099702],
                     [0.855810, 0.888601, 0.097452],
                     [0.866013, 0.889868, 0.095953],
                     [0.876168, 0.891125, 0.095250],
                     [0.886271, 0.892374, 0.095374],
                     [0.896320, 0.893616, 0.096335],
                     [0.906311, 0.894855, 0.098125],
                     [0.916242, 0.896091, 0.100717],
                     [0.926106, 0.897330, 0.104071],
                     [0.935904, 0.898570, 0.108131],
                     [0.945636, 0.899815, 0.112838],
                     [0.955300, 0.901065, 0.118128],
                     [0.964894, 0.902323, 0.123941],
                     [0.974417, 0.903590, 0.130215],
                     [0.983868, 0.904867, 0.136897],
                     [0.993248, 0.906157, 0.143936]]

    #from mpl.colors import ListedColormap
    plt.register_cmap(name='viridis', cmap=mpl.colors.ListedColormap(_viridis_data, name='viridis'))
    plt.register_cmap(name='viridis_r', cmap=mpl.colors.ListedColormap(np.flipud(_viridis_data), name='viridis_r'))