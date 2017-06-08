"""
Created on 30 May 2017
@authors: onur.kerimoglu@hzg.de, richard.hofmeister@hzg.de
example shell call:
export simf=/home/onur/WORK/projects/GB/maecs/3d/sns144-M161117n-P161118-bdyi3/sns144-M161117n-P161118-bdyi3-mergedextract_phys_FB_yoana-buhel_salt_2012-2012.nc
python plot3d_ferrybox.py $simf salt /home/onur/WORK/projects/GB/data/ferrybox yoana-buhel 2012,2013 1
"""

import sys,os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import netCDF4, netcdftime
from datetime import date,datetime,timedelta

def main(simfile,preint,gridtype,varn,fbrootpath,fbname,yrint):

    if fbname in ['cosyna-tordania', 'richard-cuximm']:
        lonlims = [0.1, 8.65]
    elif fbname in ['cosyna-funnygirl']:
        lonlims = [7.94, 8.65]
    elif fbname in ['yoana-buhel']:
        lonlims = [7.9, 8.86]

    #Get the ferrybox(fb) data
    date_fb, lon_fb, lat_fb, data_fb=get_ferrybox(fbrootpath,fbname,lonlims,varn,yrint)

    # Get model data, interpolated on a grid across fb track
    date_mt, lon_mt, lat_mt, data_mt, bg_contour, lon_FMD, lat_FMD, var_FMD= get_model_intonfb(simfile,preint,gridtype,fbname,varn,yrint,lonlims,date_fb, lat_fb, lon_fb)

    # Do the Plots
    colmap='viridis_r'
    figw,figh=(10., 5.)

    if varn == 'salt':
        cbartitle = 'Salinity [psu]'
        if fbname in ['cosyna-funnygirl', 'yoana-buhel']:
            clims=[15, 33]
            Nlev=7
        elif fbname in ['cosyna-tordania','richard-cuximm']:
            clims=[28, 34]
            Nlev=7

    vallims = [min(data_mt.min(), data_fb.min()), max(data_mt.max(), data_fb.max())]

    # cmap = cm.RdYlBu_r
    #Discrete colmap/colbar
    cmap, extopt, cbt, intbounds, allbounds = discrete_cmap_tuner(clims, vallims, Nlev=Nlev, colmap=colmap, nancol='white')

    plt.figure(figsize=(figw, figh))
    #axes positions
    ax_data = plt.axes([0.12, 0.15, 0.32, 0.8])
    ax_model = plt.axes([0.46, 0.15, 0.32, 0.8])
    cbarax = plt.axes([0.84, 0.15, 0.03, 0.5])  # common colorbor
    ax_map = plt.axes([0.78, 0.75, 0.15, 0.20])  # transect track

    dsc = ax_data.scatter(lon_fb, date_fb, s=1.0, c=data_fb, edgecolor='', cmap=cmap)
    dsc.set_clim(clims[0], clims[1])
    ax_data.set_xlabel('Longitude [degE]\nFerryBox')
    # ax_data.invert_xaxis()
    ax_data.grid(axis='y', lw=0.3)

    mpc = ax_model.pcolor(lon_mt, date_mt, data_mt, cmap=cmap)
    mpc.set_clim(clims[0], clims[1])
    ax_model.set_xlabel('Longitude [degE]\nModel Result')
    ax_model.grid(axis='y', lw=0.3)

    #handle date axis:
    if yrint[1]-yrint[0]==0:
        ax_data.yaxis.set_major_locator(mpl.dates.MonthLocator(bymonth=range(1, 13), interval=2))
        ax_model.yaxis.set_major_locator(mpl.dates.MonthLocator(bymonth=range(1, 13), interval=2))
        ax_data.yaxis.set_minor_locator(mpl.dates.MonthLocator(bymonth=range(1, 13), interval=1))
        ax_model.yaxis.set_minor_locator(mpl.dates.MonthLocator(bymonth=range(1, 13), interval=1))
        #ax_data.tick_params(axis='y', which='minor', length=0, labelsize=9)
    elif yrint[1]-yrint[0]==1:
        ax_data.yaxis.set_major_locator(mpl.dates.MonthLocator(bymonth=range(1,13), interval=3))
        ax_model.yaxis.set_major_locator(mpl.dates.MonthLocator(bymonth=range(1, 13), interval=3))
        ax_data.yaxis.set_minor_locator(mpl.dates.MonthLocator(bymonth=range(1,13), interval=1))
        ax_model.yaxis.set_minor_locator(mpl.dates.MonthLocator(bymonth=range(1,13), interval=1))
    ax_data.yaxis.set_major_formatter(mpl.dates.DateFormatter('%d.%m.%Y'))
    ax_data.tick_params(axis='y', which='major', direction='out', labelsize=10)
    ax_model.tick_params(axis='y', which='major', direction='out', labelsize=10)
    ax_data.tick_params(axis='y', which='minor', direction='out', labelsize=0)
    ax_model.tick_params(axis='y', which='minor', direction='out', labelsize=0)
    ax_model.set_yticklabels([])

    #set the xlims and ylims for both plots equally
    ax_model.set_ylim(date(yrint[0],1,1), date(yrint[1],12,31))
    ax_data.set_ylim(date(yrint[0],1,1), date(yrint[1],12,31))
    lonlims = [min(lon_mt[0], lon_fb[0]), max(lon_mt[-1], lon_fb[-1])]
    ax_data.set_xlim(lonlims[0], lonlims[1])
    ax_model.set_xlim(lonlims[0], lonlims[1])

    #discrete colorbar
    norm = mpl.colors.BoundaryNorm(intbounds, cmap.N)
    cb = mpl.colorbar.ColorbarBase(cbarax, cmap=cmap,
                                   norm=norm,
                                   boundaries=allbounds,
                                   extend=extopt,
                                   ticks=cbt,
                                   spacing='proportional')
    cb.ax.set_title(cbartitle, size=12.)

    #Plot the transect track
    if fbname in ['cosyna-funnygirl','yoana-buhel']:
        proj = getproj(setup='CuxBusHel',projpath=os.path.dirname(os.path.realpath(__file__)))
    elif fbname in ['cosyna-tordania','richard-cuximm']:
        proj = getproj(setup='SNSfull',projpath=os.path.dirname(os.path.realpath(__file__)))

    #background as a salinity contour plot
    if bg_contour==True:
        lonx,laty = proj(lon_FMD,lat_FMD)
        cmap.set_bad(str(.8))
        varavg=var_FMD[:, :-1, :-1].mean(axis=0)
        mpc=proj.pcolormesh(lonx,laty,varavg,cmap=cmap)
        mpc.set_clim(clims[0],clims[1])
        #xlabel('longitude [degrees east]')
        #ylabel('latitude [degrees north]')

    tx,ty=proj(lon_mt,lat_mt)
    proj.plot(tx, ty, 'k.',markersize=3)
    ax_map.set_yticklabels([])
    ax_map.set_xticklabels([])
    #proj.drawlsmask(land_color='.8', ocean_color='w',resolution='f')
    proj.drawcoastlines(color=(0.3, 0.3, 0.3), linewidth=0.5)
    proj.fillcontinents((.8, .8, .8), lake_color=(0.6, 0.6, 1.0))

    fname = simfile.replace('.nc', '_ferrybox_%s_%s.png' %(varn,fbname))
    plt.savefig(fname, dpi=300)
    plt.close()
    # show()
    print 'figure saved:%s' % fname


def get_model_intonfb(simfile,preint,gridtype,fbname,varn,yrint,lonlims,date_fb, lat_fb, lon_fb):
    if preint:
        nc = netCDF4.Dataset(simfile)
        ncv = nc.variables
        nct = ncv['time']
        utime = netcdftime.utime(nct.units)
        date_mt = np.array([t.date() for t in utime.num2date(nct[:])])
        lon_mt=ncv['longitude'][:]
        lat_mt=ncv['latitude'][:]
        data_mt=ncv['data'][:]
        nc.close()
        bg_contour=False;lon_FMD=0;lat_FMD=0;var_FMD=0
    else:
        if gridtype == 'getm-sns':
            from scipy.spatial import cKDTree
            if fbname in ['cosyna-tordania', 'richard-cuximm']:
                pnum = 300  # number of grid points to which the model will be interpolated
                xsl = slice(5, 135)  # getm-sns slice
                ysl = slice(26, 90)  # getm-sns slice
            elif fbname in ['cosyna-funnygirl']:
                pnum = 100  # number of grid points along the transect to which the model will be interpolated
                xsl = slice(98, 137)  # getm-sns slice
                ysl = slice(26, 49)  # getm-sns slice
            elif fbname in ['yoana-buhel']:
                pnum = 100  # number of grid points along the transect to which the model will be interpolated
                xsl = slice(98, 137)  # getm-sns slice
                ysl = slice(26, 49)  # getm-sns slice

            nc = netCDF4.Dataset(simfile)
            ncv = nc.variables
            nct = ncv['time']
            utime = netcdftime.utime(nct.units)
            # reduce to a daily resolution
            date_mt_all = np.array([t.date() for t in utime.num2date(nct[:])])
            tind = np.where((date_mt_all >= date(yrint[0], 1, 1)) * (date_mt_all <= date(yrint[1], 12, 31)))[0]
            # days = (time - time[0])/86400. + 0.5
            date_mt = date_mt_all[tind]
            var_FMD = ncv[varn][tind, -1, :, :].squeeze()  # FMD=full model domain
            varunit=ncv[varn].units
            nc.close()

            # read coordinates
            topo = get_getm_bathymetry_cropped()
            lat_FMD = topo['latc']
            lon_FMD = topo['lonc']

            # extract the spatially relevant region from model data
            lat = lat_FMD[ysl, xsl]
            lon = lon_FMD[ysl, xsl]
            var = var_FMD[:, ysl, xsl]
            if np.ma.is_masked(var):
                var[var.mask]=np.nan #transform the masked values to nan

            # del lon_FMD,lat_FMD,var_FMD #needed later for contourplot as the background for fb track plot

            # define a transect grid based on longitudes
            lon_mt = np.linspace(lonlims[0], lonlims[1], pnum)
            dlon = float(np.diff(lon_mt[:2]))  # this is the bin size

            # average lat within the lon bin
            # lat_mt = asarray([ missmean(lat_fb[where(abs(lon_fb - ll) < dlon)]) for ll in lon_mt])
            lat_mt = np.ones(len(lon_mt)) * -99
            for i, ll in enumerate(lon_mt):
                binlats = lat_fb[np.where(abs(lon_fb - ll) < dlon)]
                if len(binlats) > 0:
                    lat_mt[i] = binlats.mean()

            # interpolate model data to transect grid
            # choose the dates from the model and do the spatial interpolation along a grid across the transect
            tnum, ynum, xnum = var.shape
            tree = cKDTree(zip(lon.flat, lat.flat))
            #calculate the weights all at once: probably faster but debugging is difficult
            #dL, gridindsL = tree.query(zip(lon_mt, lat_mt), k=4)
            #wL = 1.0 / dL ** 2

            data_mt = -99. * np.ones((tnum, pnum))
            for tidx in range(tnum):
                if date_mt[tidx] in date_fb:
                    # all at once: probably fastest
                    # data_mt[tidx,:] = np.sum(wL * var[tidx,:,:].flatten()[gridindsL], axis=1) / np.sum(wL, axis=1)
                    for pi in range(pnum):
                        #use the pre-calculated weights: probably faster but debugging is difficult
                        #intval = np.sum(wL[pi] * var[tidx, :, :].flatten()[gridindsL[pi]]) / np.sum(wL[pi])
                        #calculate the weights here: probably slow, but
                        dLi, gridindsLi = tree.query((lon_mt[pi],lat_mt[pi]), k=4)
                        wLi = 1.0 / dLi ** 2
                        intvali=np.sum(wLi * var[tidx, :, :].flatten()[gridindsLi]) / np.sum(wLi)
                        data_mt[tidx, pi] = intvali
                    data_mt[tidx][np.isnan(data_mt[tidx])] = -99.
                    data_mt[tidx][np.where(lat_mt == -99.)] = -99.

            data_mt = np.ma.masked_equal(data_mt, -99.)
            lat_mt = np.ma.masked_equal(lat_mt, -99.)

            #write as a nc file:
            fname=simfile.replace('.nc', '_FB_%s_%s_%s-%s.nc'%(fbname,varn,yrint[0],yrint[1]))
            write_modintnc(fname,date_mt,lon_mt,lat_mt,data_mt,varunit)
            bg_contour=True
    return (date_mt, lon_mt, lat_mt, data_mt, bg_contour, lon_FMD, lat_FMD, var_FMD)

def write_modintnc(fname,dates,lons,lats,data,varunit):
    import time
    refdate=datetime(2000, 1, 1,0,0,0)
    nc = netCDF4.Dataset(fname, 'w', format="NETCDF4")
    nc.history = 'Created ' + time.ctime(time.time())
    nc.source = fname.split('/')[-1]

    # define time dimension and variables
    dim_time = nc.createDimension('time', dates.size)
    dim_lon   = nc.createDimension('longitude', lons.size)
    var_time = nc.createVariable('time', 'd', ('time'))
    var_time.units = 'seconds since ' + str(refdate)
    var_lon = nc.createVariable('longitude', 'd', ('longitude'))
    var_lat = nc.createVariable('latitude', 'd', ('longitude'))
    var_data = nc.createVariable('data', 'd', ('time','longitude'))
    var_data.units = varunit

    #write variables
    var_lon[:] = lons
    var_lat[:] = lats
    var_data[:] = data
    #time: convert dates to seconds
    dts = [datetime(d.year, d.month, d.day, 12, 0, 0) for d in dates] # first dates to datetimes, assuming 12:00 as H:M
    deldt_vec = [date - refdate for date in dts]
    var_time[:] = [deldt.seconds + 86400 * deldt.days for deldt in deldt_vec];

    nc.sync()
    nc.close()
    print('For the model data interpolated on ferry transect, NC file created:%s'%fname)

def get_ferrybox(fbrootpath, fbname, lonlims,varn, yrint):

    if fbname[0:7] == 'richard':
        nc = netCDF4.Dataset(os.path.join(fbrootpath, fbname, fbname + '_' + varn + '.nc'))
        nct = nc['doy']
        utime = netcdftime.utime(nct.units)
        tvec_all = np.array([t.date() for t in utime.num2date(nct[:])])
        tinds = np.where((tvec_all >= date(yrint[0], 1, 1)) * (tvec_all <= date(yrint[1], 12, 31)))[0]
        tvecD = tvec_all[tinds]  # it's a dimension variable
        lonsD = nc['lon'][:]  # it's a dimension variable
        loninds=(lonsD>=lonlims[0]) & (lonsD<=lonlims[1])
        lonsD= lonsD[loninds]
        # vectorize the matrices
        date_fb = np.repeat(tvecD, len(lonsD))
        lon_fb = np.tile(lonsD, len(tinds))
        lat_fb = nc['lat'][tinds, loninds].flatten('C')  # 'F' would be column-wise
        data_fb = nc[varn][tinds, loninds].flatten('C')
        nc.close()
    elif fbname[0:6] == 'cosyna' or fbname[0:5] == 'yoana':  # in these cases, data from different years are in different files
        # define empty arrays to be appended with data from each file
        date_fb = np.array([]);
        lon_fb = np.array([]);
        lat_fb = np.array([]);
        data_fb = np.array([])
        for y in range(yrint[0], yrint[1] + 1):
            fname = os.path.join(fbrootpath, fbname, fbname + '_' + varn + '_' + str(y))
            date_fb, lon_fb, lat_fb, data_fb = append_yearlydata(fname, fbname, y, lonlims, date_fb, lon_fb, lat_fb, data_fb)
    else:
        raise (ValueError('Unregistered file type:%s' % fbname))
    return (date_fb, lon_fb, lat_fb, data_fb)

def append_yearlydata(fname, fbname, y, lonlims,dates, lons, lats, values):
    if fbname[0:6] == 'cosyna':
        nc = netCDF4.Dataset(fname + '.nc')
        nct = nc['time']
        utime = netcdftime.utime(nct.units)
        datesC = np.array([t.date() for t in utime.num2date(nct[:])])
        dataC = np.hstack(nc['latitude'][:], nc['longitude'][:], nc['data'][:])
        nc.close()
    elif fbname[0:5] == 'yoana':
        hl = 1  # number of header lines
        valsep = ','  # character separating the values in each line
        f = open(fname + '.dat')
        lines = f.readlines()[hl:]
        f.close()
        # parse the current set
        dataC = np.nan * np.ones((len(lines), 3))
        datesC = [None] * len(lines)
        for lno, l in enumerate(lines):
            dat = [float(val) for val in l.split(valsep)]
            if len(dat) >= 4:
                ldate = datetime.fromordinal(int(dat[0])) + timedelta(days=dat[0] % 1) - timedelta(days=366)
                datesC[lno] = ldate.date()
                dataC[lno, :] = np.array([dat[1], dat[2], dat[3]]) #lat,lon,value
    datesC = np.asarray(datesC)  # convert date list to an indexible array
    #filter lons
    loninds=(dataC[:,1]>=lonlims[0])&(dataC[:,1]<=lonlims[1])
    dataC=dataC[loninds,:]
    datesC=datesC[loninds]
    # remove nan entries
    nani = np.where(np.isnan(dataC[:, 2]))[0]
    if len(nani) > 0:
        dataC = np.delete(dataC, nani, axis=0)
        datesC = np.delete(datesC, nani, axis=0)
    # filter only the requested year
    tind = np.where((datesC >= date(y, 1, 1)) * (datesC <= date(y, 12, 31)))[0]
    # append the current set to the existing
    dates = np.hstack((dates, datesC[tind]))
    lats = np.hstack((lats, dataC[tind, 0]))
    lons = np.hstack((lons, dataC[tind, 1]))
    values = np.hstack((values, dataC[tind, 2]))
    return (dates, lons, lats, values)

def discrete_cmap_tuner(clim,vallims,Nlev,colmap,nancol='white'):

    cbt = np.linspace(clim[0], clim[1], Nlev)
    cbtstep = cbt[1] - cbt[0]
    intbounds = list(cbt)

    cmap = discrete_cmap(Nlev + 1, colmap)
    # cmap = discrete_cmap(Nlev, colmap)

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


def getproj(setup,projpath):
    pickledproj=os.path.join(projpath,'proj.'+setup+'.pickle')
    if os.path.isfile(pickledproj):
        print 'opening an existing projection: '+ pickledproj
        #if a projection exists, just load it (fast)
        (proj,) = np.load(pickledproj)
    else:
        from mpl_toolkits.basemap import Basemap
        import pickle
        print 'projecting for: '+ setup
        if setup=='SNSfull':
            proj=Basemap(projection='lcc',
                       resolution='i',
                       llcrnrlon=0.0,
                       llcrnrlat=51.0, #51.2,
                       urcrnrlon=9.5,
                       urcrnrlat=56.0, #55.8,
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

if __name__=='__main__':
    if len(sys.argv) > 1:
        simfile = sys.argv[1]
    else:
        #simfile = '/home/onur/WORK/projects/GB/maecs/3d/sns144-M161117-P161118/sns144-M161117-P161118-mergedextract_phys_zSB_2009-2010.nc'
        #simfile = '/home/onur/WORK/projects/GB/maecs/3d/sns144-M161117n-P161118-bdyi3/sns144-M161117n-P161118-bdyi3-mergedextract_phys.nc'
        #simfile='/home/onur/WORK/projects/GB/maecs/3d/sns144-M161117n-P161118-bdyi3-z01mm/sns144-M161117n-P161118-bdyi3-z01mm-mergedextract_phys.nc'
        #simfile= '/home/onur/WORK/projects/GB/maecs/3d/sns144-M161117n-P161118-bdyi3/sns144-M161117n-P161118-bdyi3-mergedextract_phys_FB_yoana-buhel_salt_2012-2012.nc'
        simfile= '/workm/sns/validation/ferrybox/model/fesom.nc'

    if len(sys.argv) > 2:
        varn = sys.argv[2]
    else:
        varn = 'salt'

    if len(sys.argv) > 3:
        fbrootpath = sys.argv[3]
    else:
        #fbrootpath = '/home/onur/WORK/projects/GB/data/ferrybox'
	fbrootpath = '/workm/sns/validation/ferrybox/data'

    if len(sys.argv) > 4:
        fbname = sys.argv[4]
    else:
        #fbname = 'richard-cuximm'  # 'tordania' #funnygirl 'richard-cuximm'
        fbname = 'yoana-buhel'
        #fbname= 'cosyna-funnygirl'

    if len(sys.argv) > 5:
        yrint = [int(y) for y in sys.argv[5].split(',')]
        print yrint
    else:
        #yrint = [2010, 2010]
        yrint=[2012,2013]

    if len(sys.argv) > 6:
        preint = True if int(sys.argv[6])==1 else False
    else:
        preint = True

    if len(sys.argv) > 7:
        gridtype = sys.argv[7]
    else:
        gridtype = 'getm-sns'

    main(simfile,preint,gridtype,varn,fbrootpath,fbname,yrint)