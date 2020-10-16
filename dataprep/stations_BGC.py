__author__ = 'onur'

import os,sys
import collections
import datetime
import pickle
import numpy as np
import netCDF4
import netcdftime

#sys.path.insert(1, '/home/onur/WORK/codes/python/')
from data_tools import calc_osol,get_botdepth,create_ncfile,create_ncfile_0D,mapon_container_tvec,mapon_container_tzvec

refdate=datetime.datetime(1960,1,1)

#statnamelib={'ROTTMPT70'}

def main():
    merge=False;isoBGC=True;read=False
    #merge = True; read = False
    #stations=['JMA_Norderelbe','JMA_Suederpiep','JMA_Westerhever']
    #stations=['Helgoland', 'SAmrum','Norderelbe','Norderney'] # 'Sylt',
    #stations=['T36','T26','T41','T8','T2','T22','T5','T12','T11']
    #stations=['NOORDWK2','NOORDWK10','NOORDWK50','NOORDWK70','TERSLG4','TERSLG50','TERSLG70']
    stations=['ROTTMPT50','ROTTMPT70']#['ROTTMPT5', 'ROTTMPT15', 'ROTTMPT30','ROTTMPT50']
    #stations = ['Helgoland', 'NOORDWK2', 'NOORDWK10', 'NOORDWK50', 'NOORDWK70', 'TERSLG4', 'TERSLG50', 'TERSLG70']
    # stations=['Helgoland','Sylt','SAmrum','Norderelbe','Norderney',
    #             'T36','T26','T41','T8','T2','T22','T5','T12','T11',
    #             'NOORDWK2','NOORDWK10','NOORDWK50','NOORDWK70',
    #             'TERSLG4','TERSLG50','TERSLG70',
    #             'ROTTMPT50','ROTTMPT70']
    #stations=['Norderney']

    rootpath='/home/onur/WORK/projects/GB/data/stations/individual'
    fstats=readstats(rootpath,stations,read)

    if merge:
        mergestats(rootpath,fstats)

    if isoBGC:
        isolateBGC(rootpath,fstats)

def isolateBGC(rootpath,fstats):
    print 'Isolating'
    vars = {'DIN': 'mmolN/m3', 'NO3': 'mmolN/m3', 'NH4': 'mmolN/m3',
            'chl': 'mg/m3','DIP': 'mmolP/m3','Si': 'mmolSi/m3','zooC':'mmolC/m3'}

    # collect variables from files
    for station, fin in fstats.iteritems():
        print fin
        fout = os.path.join(rootpath,'BGCfull',station+'.nc')
        nc = netCDF4.Dataset(fin)
        ncv = nc.variables
        # procedures specific to stations
        print 'Reading:' + station,

        varkeys = get_stat_varkeys(station)
        M = []; V = []; U = []; lats = []; lons = [];
        M, V, U, lats, lons = getvar_fromstat(M, V, U, lats, lons, station, ncv, vars, varkeys,merge_statvar=False,keepZD=True)
        if hasattr(nc, 'station'):
            stationstr=nc.station
        else:
            stationstr=station
        if hasattr(nc, 'bottom_depth'):
            botdepth=nc.bottom_depth
        else:
            botdepth = get_botdepth(lons[0][0], lats[0][0], 'tree')
        nc.close()

        print ': done. Writing:'

        # if 0-d, convert to 1-d, assuming that all measurements are from surface
        if M[0].shape[1] == 2:
            for vno in range(len(M)):
                Mnew=np.ndarray((len(M[vno][:,0]), 3))
                Mnew[:,0]=M[vno][:,0]
                Mnew[:,1]=1.0
                Mnew[:,2]=M[vno][:,1]
                M[vno]=Mnew

        # write in a ncdf file:
        Mv, tvec_s, zvec = mapon_container_tzvec(M)  # map all vars along a container time/z vector
        tvec = [datetime.datetime(refdate.year, refdate.month, refdate.day, 0, 0, 0) + datetime.timedelta(0, t_s) for t_s in tvec_s]
        dims = {'t': 'time', 'z': 'depth', 'x': 'lon', 'y': 'lat'}
        create_ncfile(fout, lons[0], lats[0], Mv, V, V, U, dims, tvec, zvec, climatology=False,
                      refdate=refdate, missing_value=-99, notify=True)
        addncatt(fout, {'station': stationstr, 'bottom_depth': botdepth})

def mergestats(rootpath,fstats):
    print 'Merging'
    vars = {'DIN': 'mmolN/m3', 'DIP': 'mmolP/m3'}
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

        varkeys=get_stat_varkeys(station)
        M,V,U,lats,lons=getvar_fromstat(M,V,U,lats,lons,station,ncv,vars,varkeys,merge_statvar=True,keepZD=False,yearint=yearint,zint=zint)
        nc.close()
        print ': done.'

    #write in a ncdf file:
    Mv,tvec=mapon_container_tvec(M) #map all vars on along a container time vector
    #transform tvec to dates
    datevec = [datetime.datetime(refdate.year, refdate.month, refdate.day, 0, 0, 0) + datetime.timedelta(0, t_s) for t_s in tvec]
    create_ncfile_0D(fout_abs,datevec,-1,Mv,V,V,U,lats,lons,dims={'t':'time'},refdate=refdate,missing_value=-99,notify=True)

    #print 'Written:'+fout_abs

def get_stat_varkeys(station):
    if station== 'Norderney':
        varkeys={'phyC':'','zooC':'','chl':'CPHL','DIN':'DIN','DIP':'PHOS','Si':'SLCA', 'NH4':'AMON','NO3':'NTRA'}
    elif station in ['JMA_Norderelbe', 'JMA_Suederpiep', 'JMA_Westerhever']:
        varkeys={'DIN':'NO3','DIP':'PO4'}
    elif station in ['Norderelbe','SAmrum']:
        #measurements are NH4-N, NO3-N, NO2-N and PO4-P
        varkeys={'phyC':'','zooC':'','chl':'Chlor','NH4':'71.4*Ammon','Si':'','NO3':'71.4*Nitra','DIN':'71.4*Ammon+71.4*Nitra+71.4*Nitri','DIP':'32.285*o-Pho'}
    elif station=='Sylt':
        varkeys={'phyC':'','zooC':'','chl':'chl','DIN':'nh4+ntri+ntra','DIP':'phos','NO3':'ntra','NH4':'nh4', 'Si':'si'}
    elif station=='Helgoland':
        varkeys={'phyC':'phytoplankton','zooC':'copepods','chl':'chl-a','DIN':'DIN','DIP':'PO4','Si':'SiO4','NO3':'NO3','NH4':'NH4'}
    elif station in ['T36','T26','T41','T8','T2','T22','T5','T12','T11']:
        varkeys={'phyC':'','zooC':'','chl':'chl','DIN':'DIN','DIP':'DIP'}
    elif station[0:6] in ['NOORDW', 'TERSLG','ROTTMP']:
        # measurements are NH4-N, NO3-N, NO2-N and PO4-P
        varkeys={'phyC':'','zooC':'','chl':'CHL','DIN':'71.4*NH4+71.4*NO3+71*NO2','DIP':'32.285*PO4',
                 'SPM':'0.001*SPM','NH4':'71.4*NH4','NO3':'71.4*NO3','NO2':'71.4*NO2','PO4':'32.285*PO4',
                 'Si':'35.6*SiO2', 'SecchiDepth':'SD', 'DOC':'83.26*DOC','POC':'83.26*POC'}
        #measurements are NH4,NO3,NO2,PO4
        #varkeys = {'phyC': '', 'zooC': '', 'chl': 'CHL', 'DIN': '55.43*NH4+16.13*NO3+21.74*NO2', 'DIP':'10.56*PO4'}
    return varkeys

def getvar_fromstat(M,V,U,lat,lon,station,ncv,vars,varkeys,merge_statvar=False,keepZD=False,yearint=[0,0],zint=[0,0]):

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
            if yearint==[0,0]: #take it all
                tind = np.arange(len(dates))
            else:
                tind=np.where((years>=yearint[0]) * (years<=yearint[1]))[0]

            #reconstruct time
            datesnew=dates[tind]
            deltadate=datesnew - datetime.datetime(refdate.year,refdate.month,refdate.day,0,0)
            tvnew=[86400*deltadate[r].days +deltadate[r].seconds for r in range(len(deltadate))]

            #use the original reference date
            #tvnew = utime.date2num(dates[tind])

            if 'depth' in ncv.keys():
                zv = ncv['depth'][:]
                if zint == [0, 0]:  # take it all
                    zind = np.arange(len(v[0, :]))
                else:
                    zind = (zv > zint[0]) * (zv <= zint[1])

            #merge t, v and z in an nd array
            if keepZD and ('depth' in ncv.keys()):
                Mv = np.ndarray((len(tind)*len(zind), 3))
            else:
                Mv = np.ndarray((len(tind), 2))

            if len(v.shape)>1:
                if 'depth' in ncv.keys():
                    if keepZD:
                        for zno,zi in enumerate(zind):
                            for tno,ti in enumerate(tind):
                                #tiB=np.where(datesnew==dates[ti])
                                tiB = len(tind)*zno + tno
                                Mv[tiB, 0] = tvnew[ti]
                                Mv[tiB, 1] = zv[zi]
                                Mv[tiB, 2] = v[ti,zi]
                    else:
                        for tno,ti in enumerate(tind):
                            tiB = tno
                            Mv[tiB, 0] = tvnew[ti]
                            tiB=np.where(datesnew==dates[ti])
                            Mv[tiB,1]=np.nanmean(v[ti,zind])
                else:
                    raise(Exception('unknown dimensions encountered'))
            else:
                Mv[:,0]=tvnew
                Mv[:,1]=v[tind]

            #get rid of empty lines
            if keepZD and ('depth' in ncv.keys()):
                nanind = np.isnan(Mv[:, 2])
            else:
                nanind = np.isnan(Mv[:, 1])
            vali=np.invert(nanind)
            Mv=Mv[vali,:]

            #append
            M.append(Mv)
            if merge_statvar:
                if 'JMA' in station:
                    statname=station.split('JMA_')[1]
                else:
                    statname=station
                V.append(statname + '-' + globvar)
            else:
                V.append(globvar)
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

def read_bsh(fin,fout,dims,lat,lon,station='-'):

    print ' finding bottom depth..',
    #maxz=np.nan
    maxz = get_botdepth(lon[0], lat[0], 'tree')
    print (' found: %sN, %sE: %.1f m' % (lat[0], lon[0], maxz))

    M=[None]*9
    cols=[6,7,8,9,10,11,12,13,14]
    U=['PSU','Celcius', 'mmol/m3','mmol/m3','mmol/m3','mmol/m3','mmol/m3','ml/l','mg/m3']
    V=['sal','temp',    'nh4',    'ntri',   'ntra',   'DIN',    'DIP',    'DO',  'chl']

    with open(fin, 'rb') as csvfile:
        RDraw=csvfile.readlines()

    headers=RDraw[0].split(',')
    del RDraw[0] #delete the header
    del RDraw[0] #delete the units

    numcolO=len(headers)
    numrow=len(RDraw)
    print str(numrow)+' lines,'+str(numcolO)+' columns. ',

    #read time:
    secs=np.ndarray((numrow,1)) #seconds
    depths=np.ndarray((numrow,1))
    for  r in np.arange(numrow):
        years=int(RDraw[r].split(',')[3].split('/')[2])
        months=int(RDraw[r].split(',')[3].split('/')[0])
        days=int(RDraw[r].split(',')[3].split('/')[1])
        hours=int(RDraw[r].split(',')[4].split(':')[0])
        mins=int(RDraw[r].split(',')[4].split(':')[1])
        deltadate=datetime.datetime(years,months,days,hours,mins)- refdate
        secs[r]=86400*deltadate.days +deltadate.seconds
        depths[r]=float(RDraw[r].split(',')[5])

    #read each variable into a time-value matrix
    for varno,var in enumerate(V):
        print var,

        M[varno]=np.ndarray((numrow,3)) #seconds,depth, var
        for  r in range(numrow):
            M[varno][r,0] = secs[r]
            M[varno][r,1] = depths[r]

            val=RDraw[r].split(',')[cols[varno]].split('\n')[0]
            if val=='':
                M[varno][r,2]=np.nan
            else:
                if val[0]=='<':
                    M[varno][r,2]=0.0
                else:
                    M[varno][r,2]=float(val)

        #if O2 (ml/l) convert to O2sat (%)
        if (var=='DO') and (U[varno]=='ml/l'):
            #M[varno][:,2]=M[varno][:,2]*1000/22.391
            temp=M[V.index('temp')][:,2]
            sal=M[V.index('sal')][:,2]
            osol=calc_osol(temp,sal) #o2 solubility as f(t,s) in ml/l
            M[varno][:, 2]=M[varno][:,2]/osol*100.
            V[varno]='DOsat'
            U[varno]='% sat'
        #get rid of empty lines
        M[varno]=M[varno][np.invert(np.isnan(M[varno][:,1])),:]

    #write in a ncdf file:
    Mv,tvec_s,zvec=mapon_container_tzvec(M) #map all vars along a container time/z vector
    tvec=[datetime.datetime(refdate.year,refdate.month,refdate.day,0,0,0)+datetime.timedelta(0,t_s) for t_s in tvec_s]
    create_ncfile(fout, lon, lat, Mv, V, V, U, dims, tvec, zvec, climatology=False,
                  refdate=refdate, missing_value=-99, notify=True)
    addncatt(fout, {'station': station, 'bottom_depth': maxz})
    #print 'file written:'+fout

def read_sylt(fin,fout,dims,lat,lon,fno,station='-'):

    print ' finding bottom depth..',
    #maxz=np.nan
    maxz = get_botdepth(lon[0], lat[0], 'tree')
    print (' found: %sN, %sE: %.1f m' % (lat[0], lon[0], maxz))

    if fno==1: #f1 (1976-2008)
        M = [None] * 11  # sal,temp,ph,chl,nh4,si,ntri,ntra,phos,dop,don
        cols=[3,4,5,6,7,8,9,10,11,13,14]
        U=['PSU','Celcius','-','mg/m3','mmol/m3','mmol/m3','mmol/m3','mmol/m3','mmol/m3','mmol/m3','mmol/m3']
        V=['sal','temp','ph','chl','nh4','si','ntri','ntra','phos','dop','don']
    elif fno==2: # f1 (2008-2013)
        M = [None] * 8  # sal,temp,phos,si,nh4,ntri,ntra,chl
        cols = [2, 3, 4, 5, 6, 7, 8, 9]
        U = ['PSU', 'Celcius', 'mmol/m3', 'mmol/m3', 'mmol/m3', 'mmol/m3', 'mmol/m3', 'mg/m3']
        V = ['sal', 'temp',    'phos',    'si',      'nh4',     'ntri',    'ntra',    'chl']

    with open(fin, 'rb') as csvfile:
        RDraw=csvfile.readlines()

    headers=RDraw[0].split(',')
    del RDraw[0] #delete the header

    numcolO=len(headers)
    numrow=len(RDraw)
    print str(numrow)+' lines,'+str(numcolO)+' columns. ',

    #read time:
    secs=np.ndarray((numrow,1)) #seconds
    for  r in np.arange(numrow):
        years=int(RDraw[r].split(',')[0])
        months=int(RDraw[r].split(',')[1])
        deltadate=datetime.date(years,months,15) - datetime.date(refdate.year,refdate.month,refdate.day)
        secs[r]=86400*deltadate.days +deltadate.seconds

    #read each variable into a time-value matrix
    for varno,var in enumerate(V):
        print var,

        M[varno]=np.ndarray((numrow,2)) #seconds,var
        for  r in range(numrow):
            M[varno][r,0] = secs[r]
            M[varno][r,1]=float(RDraw[r].split(',')[cols[varno]]) if RDraw[r].split(',')[cols[varno]]!='' else np.nan
        #get rid of empty lines
        M[varno]=M[varno][np.invert(np.isnan(M[varno][:,1])),:]

    #write in a ncdf file:
    Mv,tvec_s=mapon_container_tvec(M) #map all vars on along a container time vector
    tvec=[datetime.datetime(refdate.year,refdate.month,refdate.day,0,0,0)+datetime.timedelta(0,t_s) for t_s in tvec_s]
    create_ncfile(fout, lon, lat, Mv, V, V, U, dims, tvec, zvec=-1, climatology=False,
                  refdate=refdate, missing_value=-99, notify=True)
    #print 'file written:'+fout
    addncatt(fout, {'station': station, 'bottom_depth': maxz})

def read_justusMA(fin,fout,dims,lat,lon,station='-'):

    print ' finding bottom depth..',
    #maxz=np.nan
    maxz = get_botdepth(lon[0], lat[0], 'tree')
    print (' found: %sN, %sE: %.1f m' % (lat[0], lon[0], maxz))


    M = [None] * 2  # sal,temp,phos,si,nh4,ntri,ntra,chl
    cols = [3,4]
    U = ['mmol/m3', 'mmol/m3']
    V = ['NO3',     'PO4']

    with open(fin, 'rb') as csvfile:
        RDraw=csvfile.readlines()

    headers=RDraw[0].split(',')
    del RDraw[0] #delete the header

    numcolO=len(headers)
    numrow=len(RDraw)
    print str(numrow)+' lines,'+str(numcolO)+' columns. ',

    #read time:
    secs=np.ndarray((numrow,1)) #seconds
    for  r in np.arange(numrow):
        years=int(RDraw[r].split(',')[0])
        months=int(RDraw[r].split(',')[1])
        deltadate=datetime.date(years,months,15) - datetime.date(refdate.year,refdate.month,refdate.day)
        secs[r]=86400*deltadate.days +deltadate.seconds

    #read each variable into a time-value matrix
    for varno,var in enumerate(V):
        print var,

        M[varno]=np.ndarray((numrow,2)) #seconds,var
        for  r in range(numrow):
            M[varno][r,0] = secs[r]
            val=RDraw[r].split(',')[cols[varno]].split('\n')[0] if RDraw[r].split(',')[cols[varno]]!='' else np.nan
            M[varno][r,1]=float(val)
        #get rid of empty lines
        M[varno]=M[varno][np.invert(np.isnan(M[varno][:,1])),:]

    #write in a ncdf file:
    Mv,tvec_s=mapon_container_tvec(M) #map all vars on along a container time vector
    tvec=[datetime.datetime(refdate.year,refdate.month,refdate.day,0,0,0)+datetime.timedelta(0,t_s) for t_s in tvec_s]
    create_ncfile(fout, lon, lat, Mv, V, V, U, dims, tvec, zvec=-1, climatology=False,
                  refdate=refdate, missing_value=-99, notify=True)
    #print 'file written:'+fout
    addncatt(fout, {'station': station, 'bottom_depth': maxz})

def read_opendapnc(fin_base,fout,dims,station='-'):
    vars={'CHL':'concentration_of_chlorophyll_in_water',
          'DOC':'DOC',
          'NH4':'NH4',
          'NO2':'NO2',
          'NO3':'NO3',
          'PO4':'PO4',
          'POC':'POC',
          'SD': 'secchi_depth',
          'SiO2':'SiO2',
          'SPM':'concentration_of_suspended_matter_in_water'
          }

    print ' finding bottom depth..',
    # maxz=np.nan
    maxz = get_botdepth(lon[0], lat[0], 'tree')
    print (' found: %sN, %sE: %.1f m' % (lat[0], lon[0], maxz))

    lon =-999.0; lat =-999.0
    M=[];U=[];V=[]

    print 'variable:',
    for varshort,varlong in vars.iteritems():
        print varshort,

        fin = fin_base + '-' + varshort + '.nc'
        if not os.path.isfile(fin):
            print '(not available)',
            continue
        try:
            nc = netCDF4.Dataset(fin)
        except:
            print '(not available)',
            continue

        utime = netcdftime.utime(nc.variables['time'].units)
        dates = utime.num2date(nc.variables['time'][:])
        tdels= [date-datetime.datetime(refdate.year,refdate.month,refdate.day,0,0,0) for date in dates]

        secs=np.array([86400*tdel.days+tdel.seconds for tdel in tdels])
        depths=nc.variables['z'][:].squeeze()
        if depths.size==1:
            depths=np.array([depths])
        vals=nc.variables[varlong][:].squeeze()


        if type(vals)==np.ma.core.MaskedArray:
            vali=np.invert(vals.mask)
        else:
            vali=np.invert(np.isnan(vals))

        if sum(vali)==0:
            print 'sum(vali)=0'

        Mv = np.ndarray((sum(vali), 3))
        Mv[:,0] = secs[vali]

        if depths.shape[0]==1:
            Mv[:,1]=-1*depths[0]*np.ones(sum(vali))
        else:
            Mv[:,1] = -1* depths[vali]
        Mv[:,2] = vals[vali]

        #filter outliers
        SD = np.std(Mv[:,2])
        AVG=np.mean(Mv[:,2])
        OLi = abs(Mv[:,2]-AVG)>SD*2.0
        if varshort == 'PO4':
            print ''
        Mv[OLi,2] = np.nan


        M.append(Mv)
        U.append(nc.variables[varlong].units)
        V.append(varshort)

        newlon=nc.variables['lon'][:]
        newlat=nc.variables['lat'][:]
        if lon==-999.0:
            lon=newlon
            lat=newlat
        else:
            if lon!=newlon or lat!=newlat:
                raise(ValueError('previous lon-lat pair (%s-%s) do not match with the new pair(%s-%s)'%(lon,lat,nowlon,newlat)))

        nc.close()
    print ''
    #Mv,tvec_s=mapon_container_tvec(M) #map all vars on along a container time vector
    Mv, tvec_s, zvec = mapon_container_tzvec(M)  # map all vars along a container time/z vector
    tvec=[datetime.datetime(refdate.year,refdate.month,refdate.day,0,0,0)+datetime.timedelta(0,t_s) for t_s in tvec_s]
    create_ncfile(fout, lon, lat, Mv, V, V, U, dims, tvec,zvec,
                  climatology=False,refdate=refdate, missing_value=-99, notify=True)
    #print 'file written:'+fout
    addncatt(fout, {'station': station, 'bottom_depth': maxz})

def read_SchlHolst(fin,fout,dims,lat,lon,station='-'):
    print ' finding bottom depth..',
    # maxz=np.nan
    maxz = get_botdepth(lon[0], lat[0], 'tree')
    print (' found: %sN, %sE: %.1f m' % (lat[0], lon[0], maxz))

    M=[None]*5 #chl,amon,ntra,ntri,phos
    U=[None]*5
    V=[None]*5

    with open(fin, 'rb') as csvfile:
        RDraw=csvfile.readlines()

    headers=RDraw[0].split(',')
    del RDraw[0] #delete the header

    numcolO=len(headers)
    numrow=len(RDraw)
    #print str(numrow)+' lines,'+str(numcolO)+' columns. ',

    #read all variables:
    varar=np.chararray((numrow,1), itemsize=5)
    secs=np.ndarray((numrow,1)) #seconds
    for  r in np.arange(numrow):
        varar[r,0]=RDraw[r].split(',')[7]
        years=int(RDraw[r].split(',')[4].split('.')[2])
        months=int(RDraw[r].split(',')[4].split('.')[1])
        days=int(RDraw[r].split(',')[4].split('.')[0])
        #hours=int(RDraw[r].split(',')[2].split(':')[0])
        #mins=int(RDraw[r].split(',')[2].split(':')[1])
        #deltadate=datetime.datetime(years,months,days,hours,mins)- refdate
        deltadate=datetime.date(years,months,days)- datetime.date(refdate.year,refdate.month,refdate.day)
        secs[r]=86400*deltadate.days +deltadate.seconds

    #read each variable into a time-value matrix
    uvars=np.unique(varar)
    for varno,var in enumerate(uvars):
        vari=np.array([varar[r][0]==var for r in range(len(varar))])
        secsV=secs[vari]
        RD=[None]*len(secsV)
        for r0,r1 in enumerate(np.where(vari)[0]):
            RD[r0]=RDraw[r1] #i.e., skip the header

        numrow=len(RD)
        print var +':'+str(numrow)+' lines,'+str(numcolO)+' columns. ',

        U[varno]=RD[0].split(',')[11]
        V[varno]=var
        M[varno]=np.ndarray((numrow,2)) #seconds,var
        for  r in range(numrow):
            M[varno][r,0] = secsV[r]
            M[varno][r,1]=float(RD[r].split(',')[10]) if RD[r].split(',')[10]!='' else np.nan

        #get rid of empty lines
        M[varno]=M[varno][np.invert(np.isnan(M[varno][:,1])),:]

    #write in a ncdf file:
    Mv,tvec_s=mapon_container_tvec(M) #map all vars on along a container time vector
    tvec=[datetime.datetime(refdate.year,refdate.month,refdate.day,0,0,0)+datetime.timedelta(0,t_s) for t_s in tvec_s]
    create_ncfile(fout, lon, lat, Mv, V, V, U, dims, tvec, zvec=-1, climatology=False,
                  refdate=refdate, missing_value=-99, notify=True)
    #print 'file written:'+fout
    addncatt(fout, {'station': station, 'bottom_depth': maxz})


def read_norderney(rootpath,fin,fout,dims,lat,lon,station='-'):
    print ' finding bottom depth..',
    # maxz=np.nan
    maxz = get_botdepth(lon[0], lat[0], 'tree')
    print (' found: %sN, %sE: %.1f m' % (lat[0], lon[0], maxz))

    M=[None]*10 #chl,amon,din,ntra,ntri,phos,ntot,ptot,slca,doc
    U=[None]*10
    V=[None]*10

    #read chl
    with open(fin, 'rb') as csvfile:
        RD=csvfile.readlines()

    headers=RD[0].split(',')
    del RD[0]
    numrow=len(RD)

    numcolO=len(headers)
    print 'CHL:'+ str(numrow)+' lines,'+str(numcolO)+' columns. ',

    M[0]=np.ndarray((numrow,2)) #seconds,chl
    for  r in np.arange(0,numrow):
        years=int(RD[r].split(',')[1].split('.')[2])
        months=int(RD[r].split(',')[1].split('.')[1])
        days=int(RD[r].split(',')[1].split('.')[0])
        hours=int(RD[r].split(',')[2].split(':')[0])
        mins=int(RD[r].split(',')[2].split(':')[1])
        deltadate=datetime.datetime(years,months,days,hours,mins)- refdate
        M[0][r,0] = 86400*deltadate.days +deltadate.seconds
        M[0][r,1]=float(RD[r].split(',')[10]) if RD[r].split(',')[10]!='' else np.nan

    #get rid of empty lines
    M[0]=M[0][np.invert(np.isnan(M[0][:,1])),:]
    U[0]='mg/m3'
    V[0]='CPHL'

    #read nut
    with open(os.path.join(rootpath,'Norderney','Nney_W_2_1999-2014_nutrients.csv'), 'rb') as csvfile:
        RDraw=csvfile.readlines()

    headers=RDraw[0].split(',')
    del RDraw[0]
    numcolO=len(headers)
    numrow=len(RDraw)
    #print str(numrow)+' lines,'+str(numcolO)+' columns. ',

    #read all variables:
    varar=np.chararray((numrow,1), itemsize=5)
    secs=np.ndarray((numrow,1)) #seconds
    for  r in np.arange(numrow):
        varar[r,0]=RDraw[r].split(',')[11]
        years=int(RDraw[r].split(',')[1].split('.')[2])
        months=int(RDraw[r].split(',')[1].split('.')[1])
        days=int(RDraw[r].split(',')[1].split('.')[0])
        hours=int(RDraw[r].split(',')[2].split(':')[0])
        mins=int(RDraw[r].split(',')[2].split(':')[1])
        deltadate=datetime.datetime(years,months,days,hours,mins)- refdate
        secs[r]=86400*deltadate.days +deltadate.seconds

    #read each variable into a time-value matrix
    uvars=np.unique(varar)
    for varno,var in enumerate(uvars):
        vari=np.array([varar[r][0]==var for r in range(len(varar))])
        secsV=secs[vari]
        RD=[None]*len(secsV)
        for r0,r1 in enumerate(np.where(vari)[0]):
            RD[r0]=RDraw[r1] #i.e., skip the header

        numrow=len(RD)
        print var +':'+str(numrow)+' lines,'+str(numcolO)+' columns. ',

        U[varno+1]=RD[0].split(',')[13]
        V[varno+1]=var
        M[varno+1]=np.ndarray((numrow,2)) #seconds,var
        for  r in range(numrow):
            M[varno+1][r,0] = secsV[r]
            M[varno+1][r,1]=float(RD[r].split(',')[12]) if RD[r].split(',')[12]!='' else np.nan

        #get rid of empty lines
        M[varno+1]=M[varno+1][np.invert(np.isnan(M[varno+1][:,1])),:]

    #write in a ncdf file:
    Mv,tvec_s=mapon_container_tvec(M) #map all vars on along a container time vector
    tvec=[datetime.datetime(refdate.year,refdate.month,refdate.day,0,0,0)+datetime.timedelta(0,t_s) for t_s in tvec_s]
    #create_ncfile(fout,lon,lat,tvec,-1,Mv,V,V,U,dims,refdate=refdate,missing_value=-99,notify=True)
    create_ncfile(fout, lon, lat, Mv, V, V, U, dims, tvec, zvec=-1, climatology=False,
                  refdate=refdate, missing_value=-99, notify=True)
    #print 'file written:'+fout
    addncatt(fout, {'station': station, 'bottom_depth': maxz})

def read_helgoland(dims, rootpath, fout, lat, lon, station='-'):

    #read_newhelgoland()
    combine_helgoland(dims,rootpath, fout,lat,lon,station=station)

def combine_helgoland(dims,rootpath, fout,lat,lon,MV=-99,station='-'):
    print ' finding bottom depth..',
    # maxz=np.nan
    maxz = get_botdepth(lon[0], lat[0], 'tree')
    print (' found: %sN, %sE: %.1f m' % (lat[0], lon[0], maxz))

    f = [None] * 2
    yint = [None] * 2
    varkeys = [None] * 2

    #varkeys when reading the output:
    #varkeys = {'phyC': 'phy-molC', 'zooC': 'cop-molC', 'chl': 'chl-a', 'DIN': 'DIN', 'DIP': 'PO4'}

    f[0] = os.path.join(rootpath, 'Helgoland', 'Helgoland_1962-2005.nc')
    yint[0] = [1962, 2004]
    varkeys[0] = {'secchi_depth':'Sec','temperature':'Temp','salinity':'Sal',
                  'PO4':'PO4','NO2':'NO2','NO3':'NO3','NH4':'NH4','DIN':'DIN','SiO4':'SiO4',
                  'chl-a':'chl-a','diatoms_count':'cell','diatoms_lESD':'lESD',
                  'diatoms':'diat-molC','flagellates':'flag-molC','phytoplankton':'phy-molC','copepods':'cop-molC'}
                  #'diat-C':'diat-C','flag-C':'flag-C','phy-C':'phy-C','cop-C':'cop-C'

    f[1] = os.path.join(rootpath, 'Helgoland', 'Helgoland_2005-2014.nc')
    yint[1] = [2005, 2014]
    varkeys[1] = {'chl-a':'Chla_BBE', 'chl-a_hplc':'Chla_Hplc', 'DIN': 'DIN', 'PO4': 'PO4', 'SiO4': 'SiO4'}

    vars = {'secchi_depth':'m','temperature':'Celsius','salinity':'o/oo',
            'PO4': 'mmolP/m3','NO2':'mmolN/m3','NO3':'mmolN/m3','NH4':'mmolN/m3', 'DIN': 'mmolN/m3','SiO4': 'mmolSi/m3',
            'diatoms_count':'ind/l','diatoms_lESD':'log (\mu m)', 'chl-a': 'mg/m3', 'chl-a_hplc': 'mg/m3',
            'diatoms':'mmolC/m3','flagellates':'mmolC/m3','phytoplankton': 'mmolC/m3', 'copepods': 'mmolC/m3'}

    yearint = [yint[0][0], yint[1][1]]
    fout_abs = os.path.join(rootpath, 'Helgoland' + '_' + str(yearint[0]) + '-' + str(yearint[1]) + '.nc')

    #find time indices
    tinds=[None]*len(f)
    datevec=np.empty((0,))
    for fno, fin in enumerate(f):
        nc = netCDF4.Dataset(fin)
        #time
        tv = nc.variables['time']
        utime = netcdftime.utime(tv.units)
        dates = utime.num2date(tv[:])
        tinds[fno]=(dates>datetime.datetime(yint[fno][0]-1,12,31,23,0,0)) * (dates<datetime.datetime(yint[fno][1]+1,1,1,1,0,0))
        datevec=np.hstack((datevec, dates[tinds[fno]]))
        #tvec=np.vstack(tvec,tv[tinds[fno]])
        nc.close()

    #convert dates to secs
    #deltadate = [date - refdate for date  in datevec]
    #tvec_s = [86400 * deltadate[r].days + deltadate[r].seconds for r in range(len(deltadate))] #
    #tvec = [86400 * deltadate[r].days for r in range(len(deltadate))] #trim away the time info

    # collect variables from files and concatenate
    M = []; V = []; L=[]; U = [];
    for var,unit in vars.iteritems():
        V.append(var)
        U.append(unit)
        m=np.empty((0,))
        sname = '?'
        print var+':',
        for fno, fin in enumerate(f):
            if var in varkeys[fno]:
                varn=varkeys[fno][var]
                print 'f'+str(fno+1)+':'+varn,
                nc = netCDF4.Dataset(fin)
                #if varn in nc.variables:
                v0=nc.variables[varn][:].squeeze()
                v=v0[tinds[fno]]
                v.mask=False
                v[v==nc.variables[varn]._FillValue]=np.nan
                if fno==0:
                    sname=str(nc.variables[varn].standard_name)
                nc.close()
            else:
                print 'f'+str(fno+1)+': NA',
                v=np.nan*np.ones((sum(tinds[fno]),))
            m = np.hstack((m, v))
        print '.'
        m[np.isnan(m)]=MV
        M.append(m)
        L.append(sname)

    create_ncfile(fout, lon, lat, M, V, L, U, dims, datevec, -1, climatology=False,
                  refdate=refdate, missing_value=MV, notify=True)
    addncatt(fout, {'station': station, 'bottom_depth': maxz})

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

    coords={'Helgoland': [54.18,7.90],  #54 11'N 7.54'E
             'Norderney':[53.7,7.17],
             'SAmrum':[54.59,8.39],
             'Norderelbe':[54.0,8.67],
             'Sylt': [55.0, 8.45],  # 55.0'N 8.27'E
             'Suederpiep':[54.10, 8.45],
             'Westerhever':[54.37, 8.50],
             'T2':[55.2,5.0],
             'T5':[55.0,6.33],
             'T8':[55.0,8.0],
             'T11':[54.67,6.92],
             'T12':[54.67,7.42],
             'T22':[54.13,6.35],
             'T26':[54.18,7.46],
             'T36':[53.68,6.42],
             'T41':[54.0,8.11]
            }

    dims={'t':'time','x':'lon','y':'lat'}

    fstats=collections.OrderedDict({})
    for sno,station in enumerate(stations):

        print 'Reading:'+station
        if station=='Helgoland':
            #fin=
            fout=os.path.join(rootpath, 'Helgoland', 'Helgoland.nc')
            if read: read_helgoland(dims, rootpath, fout, lat=np.array([coords[station][0]]),lon=np.array([coords[station][1]]),station='Helgoland')
        elif station=='Norderney':
            fin=os.path.join(rootpath,'Norderney','Nney_W_2_1999-2014_chl.csv')
            fout=os.path.join(rootpath, 'Norderney', 'Norderney.nc')
            if read: read_norderney(rootpath,fin, fout,dims,lat=np.array([coords[station][0]]),lon=np.array([coords[station][1]]),station='Norderney')
        elif station== 'SAmrum':
            fin=os.path.join(rootpath,'SchlHolst','DIN-DIP-Chl-220006-220065-2000-2014_SAmrum.csv')
            fout=os.path.join(rootpath, 'SchlHolst', 'SAmrum.nc')
            if read: read_SchlHolst(fin,fout,dims,lat=np.array([coords[station][0]]),lon=np.array([coords[station][1]]),station='S.Amrum')
        elif station== 'Norderelbe':
            fin=os.path.join(rootpath,'SchlHolst','DIN-DIP-Chl-220006-220065-2000-2014_Norderelbe.csv')
            fout=os.path.join(rootpath, 'SchlHolst', 'Norderelbe.nc')
            if read: read_SchlHolst(fin, fout,dims,lat=np.array([coords[station][0]]),lon=np.array([coords[station][1]]),station='Norderelbe')
        elif station=='JMA_Sylt_1':
            fin=os.path.join(rootpath,'justusMA','SyltMonthlymeans(1973-2008).csv')
            fout=os.path.join(rootpath, 'justusMA', 'Sylt_1973-2008_all.nc')
            if read: read_sylt(fin, fout,dims,lat=np.array([coords['Sylt'][0]]),lon=np.array([coords['Sylt'][1]]),fno=1,station='Sylt')
        elif station=='JMA_Sylt_2':
            fin=os.path.join(rootpath,'justusMA','Sylt_2009-2013.csv')
            fout=os.path.join(rootpath, 'justusMA','Sylt_2009-2013.nc')
            if read: read_sylt(fin, fout,dims,lat=np.array([coords['Sylt'][0]]),lon=np.array([coords['Sylt'][1]]),fno=2,station='Sylt')
        elif station=='Sylt':
            fout=os.path.join(rootpath,'justusMA','Sylt_1973-2013.nc')
        elif station in ['JMA_Norderelbe','JMA_Suederpiep','JMA_Westerhever']:
            stname=station.split('JMA_')[1]
            fin = os.path.join(rootpath, 'justusMA', stname+'.csv')
            fout = os.path.join(rootpath, 'justusMA', stname+'.nc')
            if read: read_justusMA(fin, fout, dims, lat=np.array([coords[stname][0]]), lon=np.array([coords[stname][1]]),station=stname)
        elif station== 'Norderelbe':
            fin=os.path.join(rootpath,'SchlHolst','DIN-DIP-Chl-220006-220065-2000-2014_Norderelbe.csv')
            fout=os.path.join(rootpath, 'SchlHolst', 'Norderelbe.nc')
            if read: read_SchlHolst(fin, fout,dims,lat=np.array([coords[station][0]]),lon=np.array([coords[station][1]]),station='Norderelbe')
        elif station in ['T36','T26','T41','T8','T2','T22','T5','T12','T11']:
            dims={'t':'time','z':'depth','x':'lon','y':'lat'}
            fin=os.path.join(rootpath,'BSH',station+'.csv')
            fout=os.path.join(rootpath,'BSH',station+'.nc')
            if read: read_bsh(fin,fout,dims,lat=np.array([coords[station][0]]),lon=np.array([coords[station][1]]),station=station)
        elif station[0:6] in ['NOORDW','TERSLG','ROTTMP']:
            dims = {'t': 'time', 'z': 'depth', 'x': 'lon', 'y': 'lat'}
            fin=os.path.join(rootpath,'Dutch','opendap','raw',station)
            fout = os.path.join(rootpath, 'Dutch', 'opendap', 'merged', station+'.nc')
            if read: read_opendapnc(fin,fout,dims,station=station)

        fstats[station]=fout

    return fstats

if __name__=='__main__':
    main()
