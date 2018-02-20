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

from getm_funcs import get_getm_bathymetry_cropped
from general_funcs import discrete_cmap_tuner,getproj

def main(simfile,preint,gridtype,varn,fbrootpath,fbname,yrint):

    if fbname in ['cosyna-tordania', 'richard-cuximm']:
        lonlims = [0.1, 8.65]
    elif fbname in ['cosyna-funnygirl']:
        lonlims = [7.94, 8.65]
    elif fbname in ['yoana-buhel']:
        lonlims = [7.9, 8.86]

    #Get the ferrybox(fb) data
    date_fb, lon_fb, lat_fb, data_fb=get_ferrybox(fbrootpath,fbname,lonlims,varn,yrint)

    #Write the fb data on a fixed grid (although this won't be used in the remainder of this script)
    store_interpfb(fbrootpath, fbname, varn, date_fb, lon_fb, lat_fb, data_fb, lonlims, yrint)
    return

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

def store_interpfb(fbrootpath,fbname,varn,date_fb,lon_fb,lat_fb,data_fb,lonlims,yrint=[2012,2013]):
    unitdict={'salt':'psu'}
    if fbname in ['cosyna-tordania', 'richard-cuximm']:
        pnum = 300  # number of grid points to which the model will be interpolated
    elif fbname in ['cosyna-funnygirl']:
        pnum = 100  # number of grid points along the transect to which the model will be interpolated
    elif fbname in ['yoana-buhel']:
        pnum = 100  # number of grid points along the transect to which the model will be interpolated

    # define an average transect grid based on longitudes
    lon_at = np.linspace(lonlims[0], lonlims[1], pnum)
    dlonlat = float(np.diff(lon_at[:2]))  # this is the bin size (of both lon and lat

    # average lat within the lon bin
    #lat_at = np.ones(len(lon_at)) * -99
    #for i, ll in enumerate(lon_at):
    #    binlats = lat_fb[np.where(abs(lon_fb - ll) < dlonlat)]
    #    if len(binlats) > 0:
    #        lat_at[i] = binlats.mean()

    #dates to fill:
    delta=datetime(yrint[1],12,31,0,0,0)-datetime(yrint[0], 1, 1, 0, 0, 0)
    delta_sec = 1 * 24 * 60 * 60*delta.days+delta.seconds
    dates_at = np.array([date(yrint[0], 1, 1) + timedelta(0, t) for t in range(0, delta_sec, 86400)])
    tnum= len(dates_at)

    #bin the data on the grid
    data_at = -99. * np.ones((tnum, pnum))
    lat_at = -99. * np.ones((tnum,pnum))
    for di,d in enumerate(dates_at):
        if d in date_fb:
            # reduce the data to the relevant date
            ti = np.where(date_fb == d)[0]
            lon_fb_t = lon_fb[ti]
            lat_fb_t = lat_fb[ti]
            data_fb_t = data_fb[ti]
            #search relevant data for each lon bin
            for i in range(len(lon_at)):
                loni=np.where(lon_fb_t>=(lon_at[i]-dlonlat/2)) * (lon_fb_t<=(lon_at[i]+dlonlat/2.))[0]
                if len(loni)>0:
                    try:
                        data_at[di,i]=data_fb_t[loni].mean()
                        lat_at[di,i]=lat_fb_t[loni].mean()
                    except:
                        print('?')

    fname = fbname+'_binned_%s_%s-%s.nc' %(varn,yrint[0],yrint[1])
    write_fbint2nc(os.path.join(fbrootpath,fbname,fname), dates_at, lon_at, lat_at, data_at, unitdict[varn])

def write_fbint2nc(fname,dates,lons,lats,data,varunit):
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
    var_lat = nc.createVariable('latitude', 'd', ('time','longitude'))
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
            lonlatpairs=zip(lon.flat, lat.flat) #(python2 zip)
            tree = cKDTree(lonlatpairs)

            #calculate the weights all at once: probably faster but debugging is difficult
            #dL, gridindsL = tree.query(zip(lon_mt, lat_mt), k=4)
            #wL = 1.0 / dL ** 2

            data_mt = -99. * np.ones((tnum, pnum))
            for tidx in range(tnum):
                if True: #date_mt[tidx] in date_fb:
                    # all at once: probably fastest
                    # data_mt[tidx,:] = np.sum(wL * var[tidx,:,:].flatten()[gridindsL], axis=1) / np.sum(wL, axis=1)
                    for pi in range(pnum):
                        #use the pre-calculated weights: probably faster but debugging is difficult
                        #intval = np.sum(wL[pi] * var[tidx, :, :].flatten()[gridindsL[pi]]) / np.sum(wL[pi])
                        #calculate the weights here: probably slow, but
                        dLi, gridindsLi = tree.query((lon_mt[pi],lat_mt[pi]))
                        wLi = 1.0 / dLi ** 2
                        intvali=np.sum(wLi * var[tidx, :, :].flatten()[gridindsLi]) / np.sum(wLi)
                        data_mt[tidx, pi] = intvali
                    data_mt[tidx][np.isnan(data_mt[tidx])] = -99.
                    data_mt[tidx][np.where(lat_mt == -99.)] = -99.

            data_mt = np.ma.masked_equal(data_mt, -99.)
            lat_mt = np.ma.masked_equal(lat_mt, -99.)

            #write as a nc file:
            fname=simfile.replace('.nc', '_FB_%s_%s_%s-%s.nc'%(fbname,varn,yrint[0],yrint[1]))
            write_modint2nc(fname,date_mt,lon_mt,lat_mt,data_mt,varunit)
            bg_contour=True
    return (date_mt, lon_mt, lat_mt, data_mt, bg_contour, lon_FMD, lat_FMD, var_FMD)

def write_modint2nc(fname,dates,lons,lats,data,varunit):
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

if __name__=='__main__':
    if len(sys.argv) > 1:
        simfile = sys.argv[1]
    else:
        #simfile = '/home/onur/WORK/projects/GB/maecs/3d/sns144-M161117-P161118/sns144-M161117-P161118-mergedextract_phys_zSB_2009-2010.nc'
        #simfile = '/home/onur/WORK/projects/GB/maecs/3d/sns144-M161117n-P161118-bdyi3/sns144-M161117n-P161118-bdyi3-mergedextract_phys.nc'
        simfile='/home/onur/WORK/projects/GB/maecs/3d/sns144-M161117n-P161118-bdyi3-z01mm/sns144-M161117n-P161118-bdyi3-z01mm-mergedextract_phys.nc'
        #simfile= '/home/onur/WORK/projects/GB/maecs/3d/sns144-M161117n-P161118-bdyi3/sns144-M161117n-P161118-bdyi3-mergedextract_phys_FB_yoana-buhel_salt_2012-2012.nc'
        #simfile= '/workm/sns/validation/ferrybox/model/fesom.nc'

    if len(sys.argv) > 2:
        varn = sys.argv[2]
    else:
        varn = 'salt'

    if len(sys.argv) > 3:
        fbrootpath = sys.argv[3]
    else:
        fbrootpath = '/home/onur/WORK/projects/GB/data/ferrybox'
	    #fbrootpath = '/workm/sns/validation/ferrybox/data'

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
        preint = False

    if len(sys.argv) > 7:
        gridtype = sys.argv[7]
    else:
        gridtype = 'getm-sns'

    main(simfile,preint,gridtype,varn,fbrootpath,fbname,yrint)