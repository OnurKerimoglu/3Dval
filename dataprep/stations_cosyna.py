__author__ = 'onur'

import os,sys
import collections
import datetime
import pickle
import numpy as np
import netCDF4
import netcdftime

#sys.path.insert(1, '/home/onur/WORK/codes/python/')
from data_tools import get_botdepth,create_ncfile,create_ncfile_0D,mapon_container_tvec,mapon_container_tzvec

refdate=datetime.datetime(2000,1,1)
varlims={'temp':[-2.,30.], 'salt':[0.0,36.0],'DOsat':[0.0,200.0], 'Gauge':[1e-2,10.0]}

def main():
    #merge=True;read=False
    merge = False; read = True
    #stations=['marnet_twems', 'marnet_NBII','marnet_BSHdb','LS_cuxhaven','pile_HPAelbe']
    #stations=['pileG_hoernum']
    stations = ['LS_helgoland'] #['wsv_norderney','wsv_list','wsv_buesum'] #['wsv_bake-z'] 'wsv_helgoland'
    rootpath='/home/onur/WORK/projects/GB/data/stations/COSYNA/proc/'
    fstats=readstats(rootpath,stations,read)

    if merge:
        mergestats(rootpath,fstats)

def mergestats(rootpath,fstats):
    print 'Merging'
    #vars={'phyC':'mgC/m3','zooC':'mgC/m3','chl':'mg/m3','DIN':'mmolN/m3','DIP':'mmolP/m3'}
    #vars={'chl':'mg m-3','SPM':'g m-3','DOC':'mmolC m-3', 'POC':'mmolC m-3','SecchiDepth':'m', 'SiO2':'mmolSi m-3',
    #      'NH4': 'mmolN m-3','NO3': 'mmolN m-3','NO2':'mmolN m-3','PO4':'mmolP m-3'}

    yearint=[1960,2015]
    #yearint=[1970,2015]
    zint=[0,10] #only the surface values

    #build the fout
    fout=''
    for i,station in enumerate(fstats.keys()):
        if i==0:
            fout=station
        else:
            fout=fout+'_'+station
    fout_abs=os.path.join(rootpath,'stations_'+fout+'_'+str(yearint[0])+'-'+str(yearint[1])+'.nc')
    #print('fout:'+fout_abs)

    #collect variables from files
    M=[]; V=[]; U=[]; lats=[]; lons=[]
    for station,fin in fstats.iteritems():
        print fin
        nc=netCDF4.Dataset(fin)
        ncv=nc.variables
        #procedures specific to stations
        print 'Reading:'+station,

        if station== 'Norderney':
            varkeys={'phyC':'','zooC':'','chl':'CPHL','DIN':'DIN','DIP':'PHOS'}
        elif station in ['Norderelbe','SAmrum']:
            #measurements are NH4-N, NO3-N, NO2-N and PO4-P
            varkeys={'phyC':'','zooC':'','chl':'Chlor','DIN':'71.4*Ammon+71.4*Nitra+71.4*Nitri','DIP':'32.285*o-Pho'}
        elif station=='Sylt':
            varkeys={'phyC':'','zooC':'','chl':'chl','DIN':'nh4+ntri+ntra','DIP':'phos'}
        elif station=='Helgoland':
            varkeys={'phyC':'phytoplankton','zooC':'copepods','chl':'chl-a','DIN':'DIN','DIP':'PO4'}
        elif station in ['T36','T26','T41','T8','T2','T22','T5','T12','T11']:
            varkeys={'phyC':'','zooC':'','chl':'chl','DIN':'DIN','DIP':'DIP'}
        elif station[0:6] in ['NOORDW', 'TERSLG']:
            # measurements are NH4-N, NO3-N, NO2-N and PO4-P
            varkeys={'phyC':'','zooC':'','chl':'CHL','DIN':'71.4*NH4+71.4*NO3+71*NO2','DIP':'32.285*PO4',
                     'SPM':'0.001*SPM','NH4':'71.4*NH4','NO3':'71.4*NO3','NO2':'71.4*NO2','PO4':'32.285*PO4',
                     'SiO2':'35.6*SiO2', 'SecchiDepth':'SD', 'DOC':'83.26*DOC','POC':'83.26*POC'}
            #measurements are NH4,NO3,NO2,PO4
            #varkeys = {'phyC': '', 'zooC': '', 'chl': 'CHL', 'DIN': '55.43*NH4+16.13*NO3+21.74*NO2', 'DIP':'10.56*PO4'}
        M,V,U,lats,lons=getvar_fromstat(M,V,U,lats,lons,ncv,vars,varkeys,yearint,zint,station)
        nc.close()
        print ': done.'

    #write in a ncdf file:
    Mv,tvec=mapon_container_tvec(M) #map all vars on along a container time vector
    #transform tvec to dates
    datevec = [datetime.datetime(refdate.year, refdate.month, refdate.day, 0, 0, 0) + datetime.timedelta(0, t_s) for t_s in tvec]
    create_ncfile_0D(fout_abs,datevec,-1,Mv,V,V,U,lats,lons,dims={'t':'time'},refdate=refdate,missing_value=-99,notify=True)

    #print 'Written:'+fout_abs

def getvar_fromstat(M,V,U,lat,lon,ncv,vars,varkeys,yearint,zint,station):

    for globvar in vars.keys():
    #for globvar,statvar in varkeys.iteritems():
        print globvar,
        statvar=varkeys[globvar]
        #get values
        v=getaddvar(ncv,statvar)
        if v==[]:
            continue
        else:
            #filter time
            tv=ncv['time']
            #dates=[refdate+datetime.timedelta(0,t[r]) for r in range(t.shape[0])]
            if tv.units.startswith('seconds'):
                utime=netcdftime.utime(tv.units)
                dates = utime.num2date(tv[:])
            years=np.array([dates[r].year for r in range(len(dates))])
            tind=np.where((years>=yearint[0]) * (years<=yearint[1]))[0]

            #reconstruct time
            datesnew=dates[tind]
            deltadate=datesnew - datetime.datetime(refdate.year,refdate.month,refdate.day,0,0)
            tvnew=[86400*deltadate[r].days +deltadate[r].seconds for r in range(len(deltadate))]

            #use the original reference date
            #tvnew = utime.date2num(dates[tind])

            #merge t and v in an nd array
            Mv=np.ndarray((len(tind),2))
            Mv[:,0]=tvnew
            if len(v.shape)>1:
                if 'depth' in ncv.keys():
                    zv=ncv['depth'][:]
                    zind=(zv>zint[0]) * (zv<=zint[1])
                    for ti in tind:
                        tiB=np.where(datesnew==dates[ti])
                        Mv[tiB,1]=np.nanmean(v[ti,zind])
                else:
                    raise(Exception('unknown dimensions encountered'))
            else:
                Mv[:,1]=v[tind]
            #get rid of empty lines
            vali=np.invert(np.isnan(Mv[:,1]))
            Mv=Mv[vali,:]
            M.append(Mv)
            V.append(station+'-'+globvar)
            U.append(vars[globvar])
            if 'lat' in ncv.keys():
                lat.append(ncv['lat'][:])
            elif 'latitude' in ncv.keys():
                lat.append(ncv['latitude'][:])
            else:
                print 'no lat. variable found. Skipping.'

            if 'lon' in ncv.keys():
                lon.append(ncv['lon'][:])
            elif 'longitude' in ncv.keys():
                lon.append(ncv['longitude'][:])
            else:
                print 'no lat. variable found. Skipping.'

    return (M,V,U,lat,lon)

def getaddvar(ncv,var):
    if '+' in var:
        vars2add=var.split('+')
        v=0
        for var in vars2add:
            v=v+getfactvar(ncv,var)
    else:
        v=getfactvar(ncv,var)

    return v

def getfactvar(ncv,var):
    if '*' in var:
        f=float(var.split('*')[0])
        realvar=var.split('*')[1]
    else:
        f=1.0
        realvar=var
    if realvar in ncv.keys():
        v=np.squeeze(ncv[realvar][:])
        try:
            v.mask=False
        except:
            masked=False
        nonvali=(v==ncv[realvar]._FillValue) + (np.isnan(v))
        v[nonvali]=-99.0
        nonvali=(v <= -99.)
        v[nonvali]=np.nan
        vf=v*f
    else:
        print '(not found !)',
        vf=[]
    return vf

def perdelta(start, end, delta):
    curr = start
    while curr < end:
        yield curr
        curr += delta

def temp_avg(data,interval):
    timewindow='none'

    if 'secs' in interval:
        secs=int(interval.split('secs')[0])

    #define a constant time vector
    ds=refdate + datetime.timedelta(seconds=data[0, 0])
    ds=datetime.datetime(ds.year,ds.month,ds.day,ds.hour,30,0)
    de=refdate + datetime.timedelta(seconds=data[-1, 0])
    de=datetime.datetime(de.year,de.month,de.day,de.hour,00,0)
    tvec=[t for t in perdelta(ds, de, datetime.timedelta(seconds=secs))]

    dataA=np.nan * np.ones((len(tvec), 2))
    for tno,t in enumerate(tvec):
        deltadate = t - refdate
        t_sec = 86400 * deltadate.days + deltadate.seconds
        dataA[tno,0]=t_sec
        #averaging
        if timewindow=='nogap':
            tw = [t_sec-secs/2,t_sec+secs/2] #time window
        elif timewindow=='none':
            tw = [t_sec, t_sec]  # time window
        else :
            raise(ValueError('unknown timewindow option:'+timewindow))
        ti = (data[:,0]>=tw[0]) & (data[:,0]<=tw[1]) #time index for averaging
        dataA[tno,1]=np.nanmean(data[ti,1])

    return dataA

def read_wsv_gauge(fin,fout,dims,vars,vardict,units,lat,lon,tlim,station='-'):
    # find bot. depth
    print ' finding bottom depth..',
    # maxz=np.nan
    maxz = get_botdepth(lon[0], lat[0], 'tree')
    print (' found: %sN, %sE: %.1f m' % (lat[0], lon[0], maxz))

    M = [None] * len(vars)

    with open(fin, 'rb') as csvfile:
        RDraw = csvfile.readlines()

    #headers = RDraw[5].split(',')
    #headers = [header.split('\n')[0] for header in headers]
    headers=['timestamp','Gauge']
    del RDraw[0:6]  # delete the header rows

    numcolO = 2 #len(headers)
    numrow = len(RDraw)
    print str(numrow) + ' lines,' + str(numcolO) + ' columns. '

    # read time:
    secs = np.ndarray(numrow)  # seconds
    dates= np.ndarray(numrow,dtype=datetime.datetime)  # dates
    # depths=np.ndarray((numrow,1))
    for r in np.arange(numrow):
        years = int(RDraw[r][0:4])
        months = int(RDraw[r][4:6])
        days = int(RDraw[r][6:8])
        hours = np.mod(int(RDraw[r][8:10]),24)
        mins = int(RDraw[r][10:12])
        dates[r]=datetime.datetime(years, months, days, hours, mins)
        deltadate = dates[r] - refdate
        secs[r] = 86400 * deltadate.days + deltadate.seconds

    # read each variable into a time-value matrix
    for varno, var in enumerate(vars):
        print var + ', z:',
        M[varno] = np.nan * np.ones((numrow, 2))  # seconds, var

        coln = vardict[var]
        colno = headers.index(coln)

        for r0 in range(numrow):  # row in original data
            r = r0  # row in M #+numrow * zno
            M[varno][r, 0] = secs[r0]
            #M[varno][r, 1] = z

            #val = RDraw[r0].split(',')[colno].split('\n')[0]
            val=RDraw[r0].split(' ')[colno].split('\n')[0]
            if val == '':
                M[varno][r, 1] = np.nan
            else:
                M[varno][r, 1] = float(val)/100 #cm to m

        # get rid of the values outside the specified time interval
        ti=(dates>=tlim[0]) * (dates<=tlim[1])
        M[varno]=M[varno][ti,:]

        # remove the implausible values
        if station == 'wsv_bake-z':
            #maxval=5.6
            impi = (M[varno][:, 1] < varlims[var][0]) + (M[varno][:, 1] > varlims[var][1])
            M[varno][impi, 1] = np.nan
            #get rid of the values outside

        # get rid of empty lines
        M[varno] = M[varno][np.invert(np.isnan(M[varno][:, 1])), :]
        print 'done.'

        #temporal average
        M[varno]=temp_avg(M[varno],interval='1800secs')

        #convert 'Gauge'-> 'Zeta'
        if var=='Gauge':
            # need a ref value
            refSSH=np.nanmean(M[varno][:,1]) #?
            M[varno][:,1] = M[varno][:,1] - refSSH
            vars[varno]='zeta'

    # write in a ncdf file:
    Mv, tvec_s = mapon_container_tvec(M)  # map all vars along a container time/z vector
    tvec = [datetime.datetime(refdate.year, refdate.month, refdate.day, 0, 0, 0) + datetime.timedelta(0, t_s) for t_s in
            tvec_s]

    create_ncfile_0D(fout, tvec, -1, Mv, vars, vars, units, lat, lon, dims=dims, climatology=False,
                     refdate=refdate, missing_value=-99, notify=True)
    # print 'file written:'+fout
    addncatt(fout, {'station': station, 'bottom_depth': maxz})

def read_cosyna_gauge(fin,fout,dims,vars,vardict,units,lat,lon,tlim,station='-'):
    # find bot. depth
    print ' finding bottom depth..',
    # maxz=np.nan
    maxz = get_botdepth(lon[0], lat[0], 'tree')
    print (' found: %sN, %sE: %.1f m' % (lat[0], lon[0], maxz))

    M = [None] * len(vars)

    with open(fin, 'rb') as csvfile:
        RDraw = csvfile.readlines()

    headers = RDraw[5].split(',')
    headers = [header.split('\n')[0] for header in headers]
    del RDraw[0:7]  # delete the header rows

    numcolO = len(headers)
    numrow = len(RDraw)
    print str(numrow) + ' lines,' + str(numcolO) + ' columns. '

    # read time:
    secs = np.ndarray((numrow, 1))  # seconds
    secs = np.ndarray(numrow)  # seconds
    dates= np.ndarray(numrow,dtype=datetime.datetime)  # dates
    # depths=np.ndarray((numrow,1))
    for r in np.arange(numrow):
        years = int(RDraw[r].split(' ')[0].split('-')[0])
        months = int(RDraw[r].split(' ')[0].split('-')[1])
        days = int(RDraw[r].split(' ')[0].split('-')[2])
        hours = int(RDraw[r].split(' ')[1].split(':')[0])
        mins = int(RDraw[r].split(' ')[1].split(':')[1])
        dates[r]=datetime.datetime(years, months, days, hours, mins)
        deltadate = dates[r] - refdate
        secs[r] = 86400 * deltadate.days + deltadate.seconds

    # read each variable into a time-value matrix
    for varno, var in enumerate(vars):
        print var + ', z:',
        M[varno] = np.nan * np.ones((numrow, 2))  # seconds, var

        coln = vardict[var]
        colno = headers.index(coln)

        for r0 in range(numrow):  # row in original data
            r = r0  # row in M #+numrow * zno
            M[varno][r, 0] = secs[r0]
            #M[varno][r, 1] = z

            val = RDraw[r0].split(',')[colno].split('\n')[0]
            if val == '':
                M[varno][r, 1] = np.nan
            else:
                if val[0] == '<':
                    M[varno][r, 1] = 0.0
                elif val[0] > 1000.:
                    M[varno][r, 1]= float(val)/1000.
                else:
                    M[varno][r, 1] = float(val)

        # get rid of the values outside the specified time interval
        ti=(dates>=tlim[0]) * (dates<=tlim[1])
        M[varno]=M[varno][ti,:]

        # remove the implausible values
        if station == 'pileG_hoernum':
            maxval=5.6
            impi = (M[varno][:, 1] < varlims[var][0]) + (M[varno][:, 1] > maxval)
            M[varno][impi, 1] = np.nan
            #get rid of the values outside

        # get rid of empty lines
        M[varno] = M[varno][np.invert(np.isnan(M[varno][:, 1])), :]
        print 'done.'

        #convert 'Gauge'-> 'Zeta'
        if var=='Gauge':
            # need a ref value
            refSSH=np.mean(M[varno][:,1]) #?
            M[varno][:,1] = M[varno][:,1] - refSSH
            vars[varno]='zeta'

    # write in a ncdf file:
    Mv, tvec_s = mapon_container_tvec(M)  # map all vars along a container time/z vector
    tvec = [datetime.datetime(refdate.year, refdate.month, refdate.day, 0, 0, 0) + datetime.timedelta(0, t_s) for t_s in
            tvec_s]

    create_ncfile_0D(fout, tvec, -1, Mv, vars, vars, units, lat, lon, dims=dims, climatology=False,
                     refdate=refdate, missing_value=-99, notify=True)
    # print 'file written:'+fout
    addncatt(fout, {'station': station, 'bottom_depth': maxz})

def read_cosyna(fin,fout,dims,vars,vardict,units,zin,lat,lon,station='-'):

    #find bot. depth
    print ' finding bottom depth..',
    #maxz=np.nan
    maxz = get_botdepth(lon[0],lat[0],'tree')
    print (' found: %sN, %sE: %.1f m'%(lat[0],lon[0],maxz))

    M = [None] * len(vars)

    with open(fin, 'rb') as csvfile:
        RDraw = csvfile.readlines()

    headers = RDraw[5].split(',')
    headers=[header.split('\n')[0] for header in headers]
    del RDraw[0:7]  # delete the header rows

    numcolO=len(headers)
    numrow=len(RDraw)
    print str(numrow)+' lines,'+str(numcolO)+' columns. '

    #read time:
    secs=np.ndarray((numrow,1)) #seconds
    #depths=np.ndarray((numrow,1))
    for  r in np.arange(numrow):
        years=int(RDraw[r].split(' ')[0].split('-')[0])
        months=int(RDraw[r].split(' ')[0].split('-')[1])
        days=int(RDraw[r].split(' ')[0].split('-')[2])
        hours=int(RDraw[r].split(' ')[1].split(':')[0])
        mins=int(RDraw[r].split(' ')[1].split(':')[1])
        deltadate=datetime.datetime(years,months,days,hours,mins)- refdate
        secs[r]=86400*deltadate.days +deltadate.seconds

    #read each variable into a time-value matrix
    for varno,var in enumerate(vars):
        print var+', z:',
        M[varno] = np.nan*np.ones((numrow*len(zin), 3))  # seconds,depth, var
        for zno,z in enumerate(zin):
            print z,
            if len(zin)==1:
                coln=vardict[var]
            else:
                coln='%s_%3.1f'%(vardict[var],z)
            colno=headers.index(coln)

            for  r0 in range(numrow): #row in original data, for the current depth
                r=numrow*zno+r0 #real row in M
                M[varno][r,0] = secs[r0]
                M[varno][r,1] = z

                val=RDraw[r0].split(',')[colno].split('\n')[0]
                if val=='' or val=='0':
                    M[varno][r,2]=np.nan
                else:
                    if val[0]=='<':
                        M[varno][r,2]=0.0
                    else:
                        M[varno][r,2]=float(val)
        #remove the implausible values
        impval=(M[varno][:,2]<varlims[var][0]) +(M[varno][:,2]>varlims[var][1])
        M[varno][impval,2]=np.nan
        #get rid of empty lines
        M[varno]=M[varno][np.invert(np.isnan(M[varno][:,2])),:]
        print 'done.'

    #write in a ncdf file:
    Mv,tvec_s,zvec=mapon_container_tzvec(M) #map all vars along a container time/z vector
    tvec=[datetime.datetime(refdate.year,refdate.month,refdate.day,0,0,0)+datetime.timedelta(0,t_s) for t_s in tvec_s]
    create_ncfile(fout, lon, lat, Mv, vars, vars, units, dims, tvec, zvec, climatology=False,
                  refdate=refdate, missing_value=-99, notify=True)
    # print 'file written:'+fout
    addncatt(fout,{'station':station, 'bottom_depth':maxz})

def addncatt (fout,atts):
    nc=netCDF4.Dataset(fout,'a')
    for attkey,attval in atts.iteritems():
        if attkey=='station':
            nc.station=attval
        elif attkey=='bottom_depth':
            nc.bottom_depth=attval
    nc.sync()
    nc.close()

def readstats(rootpath,stations,read):

    coords={'pile_HPAelbe': [53.859, 8.944],
            'LS_cuxhaven':[53.87, 8.72],
            'LS_helgoland':[54.18,7.89],
            'marnet_BSHdb': [54.17, 7.45],
            'marnet_twems': [54.17, 6.34],
            'marnet_NBII': [55.00, 6.34],
            'pileG_hoernum': [54.70,8.28],
            'wsv_helgoland':[54.18,7.89], #54.10'33'',07.53'29
            'wsv_buesum':[54.12,8.86], #54.07'12'',08.51'35''
            'wsv_bake-z':[54.014,8.315], # 54.00'49'', 08.18'53''
            'wsv_list':[55.02,8.44], #55.00'60'', 08.26'31''
            'wsv_norderney':[53.7,7.16] # 53.41'47'', 07.09'21''
            }

    dims={'t':'time','x':'lon','y':'lat','z':'depth'}
    #udict={'PSAL': 'PSU', 'TEMP': 'Celcius', 'DOXY': 'percent']

    fstats=collections.OrderedDict({})
    for sno,station in enumerate(stations):

        print '\nReading:'+station
        if station in ['wsv_helgoland','wsv_buesum','wsv_list','wsv_bake-z','wsv_norderney']:
            rootpath='/home/onur/WORK/projects/GB/data/stations/WSV/'
            if station=='wsv_helgoland':
                fin='Helgoland-Binnenhafen-W1.zrx'; fout='Helgoland.nc'
            elif station=='wsv_buesum':
                fin='Buesum-W1.zrx'; fout='Buesum.nc'
            elif station=='wsv_list':
                fin='List-W1.zrx'; fout='List.nc'
            elif station=='wsv_bake-z':
                fin='Bake-Z-W1.zrx'; fout='Bake-Z.nc'
            elif station=='wsv_norderney':
                fin='Norderney-Riffgat-W1.zrx'; fout='Norderney.nc'
            pfin=os.path.join(rootpath,'raw',fin)
            pfout=os.path.join(rootpath, 'nc', fout)
            vars=['Gauge'];units=['m'];vardict={'Gauge':'Gauge'}
            if read: read_wsv_gauge(pfin,pfout,dims,vars,vardict,units,lat=np.array([coords[station][0]]),lon=np.array([coords[station][1]]),
                                       tlim=[datetime.datetime(2013,6,01),datetime.datetime(2013,8,30)],station=station)
        elif station == 'pileG_hoernum':
            fin=os.path.join(rootpath,'csv','gauge','PileGauge_hoernum_2012-2013_1800smean_tw900.csv')
            fout=os.path.join(rootpath, 'nc', 'gauge','Hoernum.nc')
            vars=['Gauge']
            units=['m']
            vardict={'Gauge':'Gauge'}
            if read: read_cosyna_gauge(fin,fout,dims,vars,vardict,units,lat=np.array([coords[station][0]]),lon=np.array([coords[station][1]]),
                                       tlim=[datetime.datetime(2012,01,01),datetime.datetime(2013,10,15)],station='pileG_hoernum')
        elif station=='HPAelbe':
            fin=os.path.join(rootpath,'csv','hpaelbe_2012-2014_dmean_tw3600.csv')
            fout=os.path.join(rootpath, 'nc', 'hpaelbe_2013-2014_dmean_tw3600.nc')
            vars=['temp','salt','DOsat']
            units = ['Celcius', 'PSU', '% sat']
            vardict={'temp':'WaterTemperature', 'salt':'Salinity','DOsat':'O2_Saturation'}
            zin=[0.5]
            if read: read_cosyna(fin,fout,dims,vars,vardict,units,zin,lat=np.array([coords[station][0]]),lon=np.array([coords[station][1]]),station='HPA-Elbe')
        elif station=='LS_helgoland':
            fin=os.path.join(rootpath,'csv','LS_helgoland_T_2011-2014_dmean_tw3600.csv')
            fout=os.path.join(rootpath, 'nc', 'LS_helgoland_T_2011-2014_dmean_tw3600.nc')
            vars=['temp'] #'salt'
            units = ['Celcius'] #'PSU'
            vardict={'salt':'Salinity', 'temp':'WaterTemp(FSI)'}
            zin=[0.5]
            if read: read_cosyna(fin, fout, dims, vars, vardict, units, zin, lat=np.array([coords[station][0]]),lon=np.array([coords[station][1]]), station='Helgoland')
        elif station=='LS_cuxhaven':
            fin=os.path.join(rootpath,'csv','LScuxhaven_2011-2014_dmean_tw3600.csv')
            fout=os.path.join(rootpath, 'nc', 'LScuxhaven_2011-2014_dmean_tw3600.nc')
            vars=['temp','salt','DOsat']
            units = ['Celcius', 'PSU', '% sat']
            vardict={'temp':'WaterTemp(FSI)', 'salt':'Salinity','DOsat':'DO_Saturation'}
            zin=[0.5]
            if read: read_cosyna(fin,fout,dims,vars,vardict,units,zin,lat=np.array([coords[station][0]]),lon=np.array([coords[station][1]]),station='Cuxhaven')
        elif station=='marnet_BSHdb':
            fin=os.path.join(rootpath,'csv','marnet_DB_2011-2014_dmean_tw3600.csv')
            fout=os.path.join(rootpath, 'nc', 'marnet_DB_2011-2014_dmean_tw3600.nc')
            vars=['temp','salt','DOsat']
            units = ['Celcius', 'PSU', '% sat']
            vardict={'temp':'TEMP', 'salt':'PSAL','DOsat':'DOXY'}
            zin=[6.0,30.0]
            if read: read_cosyna(fin,fout,dims,vars,vardict,units,zin,lat=np.array([coords[station][0]]),lon=np.array([coords[station][1]]),station='Deutsche Bucht')
        elif station=='marnet_twems':
            fin=os.path.join(rootpath,'csv','marnet_twems_2011-2014_dmean_tw3600.csv')
            fout=os.path.join(rootpath, 'nc', 'marnet_twems_2011-2014_dmean_tw3600.nc')
            vars=['temp','salt','DOsat']
            units = ['Celcius', 'PSU', '% sat']
            vardict={'temp':'TEMP', 'salt':'PSAL','DOsat':'DOXY'}
            zin=[6.0,30.0]
            if read: read_cosyna(fin,fout,dims,vars,vardict,units,zin,lat=np.array([coords[station][0]]),lon=np.array([coords[station][1]]),station='Ems')
        elif station=='marnet_NBII':
            fin=os.path.join(rootpath,'csv','marnet_NBII_2011-2014_dmean_tw3600.csv')
            fout=os.path.join(rootpath, 'nc', 'marnet_NBII_2011-2014_dmean_tw3600.nc')
            vars=['temp','salt','DOsat']
            units = ['Celcius', 'PSU', '% sat']
            vardict={'temp':'TEMP', 'salt':'PSAL','DOsat':'DOXY'}
            zin=[6.0,35.0]
            if read: read_cosyna(fin,fout,dims,vars,vardict,units,zin,lat=np.array([coords[station][0]]),lon=np.array([coords[station][1]]),station='NBII')
        fstats[station]=fout

    return fstats

if __name__=='__main__':
    main()
