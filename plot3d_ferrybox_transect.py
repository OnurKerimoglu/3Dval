
import os,sys
import netCDF4,netcdftime
import numpy as np
import datetime
from scipy import interpolate
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
import matplotlib
#local python modules
from plot3d_ferrybox import get_ferrybox,get_model_intonfb
from getm_funcs import get_getm_bathymetry_cropped

def main(simfile, preint, gridtype, varn, fbrootpath, fbname, yrint, fname_topo, fbdata='raw'):

    plotfb=False
    plotsim=True
    M = 7;
    #Yvec=[2013]
    Yvec=range(yrint[0],yrint[1]+1);
    # Yvec=[-1]; M=-1
    tempsuf= 'July-2012-2013' #''%s-%s'%(Y,M)
    labels=[str(Y) for Y in Yvec]

    if fbname in ['cosyna-tordania', 'richard-cuximm']:
        lonlims = [0.1, 8.65]
    elif fbname in ['cosyna-funnygirl']:
        lonlims = [7.94, 8.65]
    elif fbname in ['yoana-buhel']:
        lonlims = [7.9, 8.86]

    # Get the ferrybox(fb) data
    # Get the ferrybox(fb) data
    if fbdata=='raw':
        date_fb, lon_fb, lat_fb, data_fb = get_ferrybox(fbrootpath, fbname, lonlims, varn, yrint)
    else:
        date_fb,lon_fb,lat_fb,data_fb=get_FB_data(fnamefb)

    #interpolate sim on fb transect
    print ('interpolating')
    data_mt_vec=[None]*len(Yvec)
    datafb_mt_vec = [None] * len(Yvec)
    for Yno,Y in enumerate(Yvec):
        date_mt, lon_mt, lat_mt, data_mt, bg_contour, lon_FMD, lat_FMD, var_FMD = get_model_intonfb(simfile, preint,
                                                                                                gridtype, fbname, varn,
                                                                                                yrint, lonlims, date_fb,
                                                                                                lat_fb, lon_fb,
                                                                                                fname_topo,Y,M)
        data_mt_vec[Yno]=data_mt
        if plotfb == True:
            datafb_mt_vec[Yno] = regrid_fb(date_fb, lon_fb, lat_fb, data_fb, date_mt, lon_mt, lat_mt, M, Y)

    # plot
    print ('plotting')
    if plotsim==True:
        plotfname = simfile.replace('.nc', '_ferrybox_%s_%s_prof_%s.png' % (varn, fbname, tempsuf))
        plot_transects_loop(plotfname,lon_mt,lat_mt,data_mt_vec,varn,labels)

    if plotfb==True:
        plotfname=os.path.join(fbrootpath, fbname, fbname + '_%s_prof_%s.png' % (varn, tempsuf))
        plot_transects_loop(plotfname, lon_mt, lat_mt, datafb_mt_vec, varn, labels)

def plot_transects_loop(plotfname,lon,lat,data_mt_vec,var,labels):
    udict={'salt':'g/kg'}
    lims={'salt':[12, 33]}
    lcols={'2012':'darkblue','2013':'darkorange','2014':'g'}
    fcols={'2012':'lightskyblue','2013':'orange','2014':'yellowgreen'}

    fig = plt.figure(figsize=(4, 2.5), dpi=300)
    fig.subplots_adjust(hspace=.4, wspace=.5, left=0.12, right=0.98, bottom=0.17,top=0.88)
    ax = plt.subplot(1, 1, 1)
    ax.tick_params(axis='both', which='both', direction='in')
    matplotlib.rcParams.update({'font.size': 10})
    matplotlib.rcParams['xtick.direction']='in'
    matplotlib.rcParams['ytick.direction'] = 'in'
    for dno,data_mt in enumerate(data_mt_vec):
        data_tmean=data_mt.mean(axis=0)
        data_tstd = data_mt.std(axis=0)
        cumdist= np.flipud(lonlat_to_cumdist(np.flipud(lon),np.flipud(lat)))
        ax.plot(cumdist,data_tmean, lcols[labels[dno]],lw=2,label=labels[dno])
        #data_p25 = np.percentile(data_mt.filled(np.nan),25,axis=0)
        #data_p75 = np.percentile(data_mt.filled(np.nan), 75, axis=0)
        #ax.fill_between(cumdist, data_p25, data_p75, facecolor=fcols[dno], edgecolor=fcols[dno], alpha='0.5')
        ax.fill_between(cumdist, data_tmean-data_tstd, data_tmean+data_tstd,
                        facecolor=fcols[labels[dno]], edgecolor=fcols[labels[dno]], alpha='0.5')

        #disti=np.where(cumdist<13)[0]
        disti = np.where(cumdist < 30)[0]
        y = data_tmean[disti[0]]
        if dno==0:
            ax.text(36, y+1., labels[dno], color=lcols[labels[dno]], size=12, verticalalignment='bottom')
        elif dno==1:
            y=data_tmean[disti[0]] - 0.5
            ax.text(36, y-1.5, labels[dno], color=lcols[labels[dno]], size=12, verticalalignment='top')

    #ax.set_title('Ferry Data')
    ax.set_title('Simulation')
    ax.title.set_position([.5,1.05])
    #plt.legend(labels,loc='lower left')
    ax.invert_xaxis()

    plt.xlim(61, 4)
    plt.ylim(lims[var][0], lims[var][1])
    plt.xlabel('Distance from Coast [km]')
    if var=='salt':
        plt.ylabel('%s [%s]'%('Salinity',udict[var]))
    else:
        plt.ylabel('%s [%s]' % (var, udict[var]))
    #plt.show()
    plt.savefig(plotfname, dpi=300)
    print ('figure saved:%s'%plotfname)


def regrid_fb(dates,lon,lat,data,date_mt,lon_mt,lat_mt,M,Y):
    #filter time
    dates,lon,lat,data=timesel_fb(dates,lon,lat,data,M,Y)

    tnum=len(date_mt); pnum=len(lon_mt)
    data_mt = -99. * np.ones((tnum, pnum))
    for tit,date in enumerate(date_mt):

        # create a tree
        ti=np.where(dates==date)[0]
        if len(ti)>0:
            lonlatpairs = zip(lon[ti], lat[ti])  # (python2 zip)
            try:
                tree = cKDTree(lonlatpairs)
            except:
                print ('?')
            datad=data[ti]

            for pit in range(pnum):
                dLi, gridindsLi = tree.query((lon_mt[pit], lat_mt[pit]),k=4)
                wLi = 1.0 / dLi ** 2
                intvali = np.sum(wLi * datad[gridindsLi]) / np.sum(wLi)
                #intvali = np.sum(wLi * data[tidx, :, :].flatten()[gridindsLi]) / np.sum(wLi)
                data_mt[tit, pit] = intvali
        data_mt[tit,np.isnan(data_mt[tit,:])] = -99.
        data_mt[tit,np.where(lat_mt == -99.)] = -99.

    data_mt = np.ma.masked_equal(data_mt, -99.)
    return(data_mt)

def timesel_fb(dates,lon,lat,data,M,Y,yop='inc'):
    mons=np.array([d.month for d in dates])
    years = np.array([d.year for d in dates])
    if yop=='inc':
        if type(M)==list:
            ti = (mons>=M[0]) * (mons<=M[1]) * (years == Y)
        else:
            ti=(mons==M)*(years==Y)
        tempsuf = '%s-%s' % (Y, M)
    else:
        if type(M)==list:
            ti = (mons>=M[0]) * (mons<=M[1]) * (years != Y)
        else:
            ti=(mons==M)*(years!=Y)
        tempsuf = '%s-%s' % (Ystr, M)
    #return (dates[ti],lat[ti,:],data[ti,:])
    return (dates[ti],lon[ti],lat[ti], data[ti])

def get_sim(fname,var):
    nc = netCDF4.Dataset(fname)
    tv = nc.variables['time']
    utime = netcdftime.utime(tv.units)
    dates = utime.num2date(tv[:])
    valsS=nc.variables[var][:]
    topo = get_getm_bathymetry_cropped()  # fname_topo
    lons,lats=(topo['lons'], topo['lats'])
    return (dates,lons,lats,valsS)

def get_FB_data(fname):
    nc=netCDF4.Dataset(fname)
    tv=nc.variables['time']
    utime=netcdftime.utime(tv.units)
    dates=utime.num2date(tv[:])
    lat = nc.variables['latitude'][:]
    lat = np.ma.masked_array(lat,mask=lat==-99)
    lon = nc.variables['longitude'][:]
    valsfb = nc.variables['data'][:]
    valsfb=np.ma.masked_array(valsfb,mask=valsfb==-99)
    nc.close()
    return(dates,lon,lat,valsfb)

def lonlat_to_cumdist(lon, lat):
    dist = np.zeros(len(lon))
    for il in range(0, len(lon) - 1):
        dist[il + 1] = distance([lat[il], lon[il]], [lat[il + 1], lon[il + 1]])

    cumdist = np.cumsum(dist)
    return cumdist

def distance(origin, destination):
#calculates the length of an arch on a sphere defined by start and end coordinates

    import math
    lat1, lon1 = origin
    lat2, lon2 = destination
    radius = 6371 # km
    dlat = math.radians(lat2-lat1)
    dlon = math.radians(lon2-lon1)
    a = math.sin(dlat/2) * math.sin(dlat/2) + math.cos(math.radians(lat1)) \
    * math.cos(math.radians(lat2)) * math.sin(dlon/2) * math.sin(dlon/2)
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
    d = radius * c
    return d


if __name__=="__main__":

    if len(sys.argv) > 1:
        simfile = sys.argv[1]
    else:
        #simfile = '/home/onur/WORK/projects/GB/maecs/3d/sns144-M161117n-P161118-bdyi3-z01mm/sns144-M161117n-P161118-bdyi3-z01mm-mergedextract_phys.nc'
        # simfile = '/home/onur/WORK/projects/GB/maecs/3d/sns144-M161117-P161118/sns144-M161117-P161118-mergedextract_phys_zSB_2009-2010.nc'
        #simfile = '/home/onur/WORK/projects/2013/maecs/sns144-M180109-nFpBpr-Pbg2017-B180106-vsdetp4b1/sns144-M180109-nFpBpr-Pbg2017-B180106-vsdetp4b1_TS_12-13_zSB.nc'
        #simfile = '/home/onur/WORK/projects/2013/maecs/sns144-M180109-nFpBpr-Pbg2017-B180106-vsdetp4b1-AHm1/extract_Mphysred_sns144-M180109-nFpBpr-Pbg2017-B180106-vsdetp4b1-AHm1.12-13_S3.nc'
        #simfile = '/home/onur/WORK/projects/2013/maecs/sns144-M180109-nFpB-Pbg2017-B180106-vsdetp4b1-M12R13/sns144-M180109-nFpB-Pbg2017-B180106-vsdetp4b1-M12R13-1314_phys_zSB.nc'
        #simfile = '/home/onur/WORK/projects/2013/maecs/sns144-M180109-nFpB-Pbg2017-B180106-vsdetp4b1/extract_phys_sns144-M180109-nFpB-Pbg2017-B180106-vsdetp4b1.2013_zSB.nc'
        #simfile = '/home/onur/WORK/projects/2013/maecs/sns144-M180109-nFpB-Pbg2017-B180106-vsdetp4b1-M13R12/extract_phys_sns144-M180109-nFpB-Pbg2017-B180106-vsdetp4b1-M13R12.2013_zSB.nc'
        #simfile = '/home/onur/WORK/projects/2013/maecs/sns144-M180109-nFpBpr-Pbg2017-B180106-vsdetp4b1-AHm2-c06-Ec01/extract_Mphysred_sns144-M180109-nFpBpr-Pbg2017-B180106-vsdetp4b1-AHm2-c06-Ec01.12-13_S3.nc'
        #simfile = '/home/onur/WORK/projects/2013/gb300/GB300_2013-0608_dmean.nc'
        #simfile = '/home/onur/WORK/projects/2013/gb300/GB300_2012-06_2013-0608_dmean.nc'
        #simfile = '/home/onur/WORK/projects/2013/ecosmoHR/2019-04-09/ef1ssv_12-13_GB_T.nc'
        simfile = '/home/onur/WORK/projects/2013/gpmeh/sns144-GPMEH-G191007-Fnew3-PPZZSi-PinR-P191010-vS/extract_MphysC_sns144-GPMEH-G191007-Fnew3-PPZZSi-PinR-P191010-vS-L1719avg.2012-2013.nc'

    #fname_topo = '/home/onur/WORK/projects/2013/gb300/topo_german_bight_300m_v06.0_curv.nc'
    fname_topo = '/home/onur/WORK/projects/GB/data/topo/topo_area_sns.nc'

    if len(sys.argv) > 2:
        varn = sys.argv[2]
    else:
        varn = 'salt'

    if len(sys.argv) > 3:
        fbrootpath = sys.argv[3]
    else:
        fbrootpath = '/home/onur/WORK/projects/GB/data/ferrybox'

    if len(sys.argv) > 4:
        fbname = sys.argv[4]
    else:
        fbname = 'yoana-buhel'  # for 2013

    if len(sys.argv) > 5:
        yrint = [int(y) for y in sys.argv[5].split(',')]
        print (yrint)
    else:
        #yrint = [2013, 2013]
        yrint = [2012, 2013]

    if len(sys.argv) > 6:
        preint = True if int(sys.argv[6]) == 1 else False
    else:
        preint = False

    if len(sys.argv) > 7:
        gridtype = sys.argv[7]
    else:
        gridtype = 'getm-sns'   # 'ecosmoHR' 'getm-gb300' #'getm-sns'

    main(simfile, preint, gridtype, varn, fbrootpath, fbname, yrint, fname_topo)
