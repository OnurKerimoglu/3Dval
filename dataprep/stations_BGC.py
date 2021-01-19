__author__ = 'onur'

import os,sys
import collections
import datetime
import pickle
import numpy as np
import netCDF4
import netcdftime

home=os.getcwd()
#sys.path.insert(1, home+'/3Dsetups/postprocess')
#sys.path.insert(1, home+'/3Dval/dataprep')
from data_tools import calc_osol,get_botdepth,create_ncfile,create_ncfile_0D,mapon_container_tvec,mapon_container_tzvec

refdate=datetime.datetime(1960,1,1)

#statnamelib={'ROTTMPT70'}

def main():
    merge=False;isoBGC=False;read=True
    #merge = True; read = False
    #stations=['JMA_Norderelbe','JMA_Suederpiep','JMA_Westerhever']
    #stations=['Helgoland', 'SAmrum','Norderelbe','Norderney'] # 'Sylt',
    #stations=['T36','T26','T41','T8','T2','T22','T5','T12','T11']
    #stations=['BOCHTVWTM','BOOMKDP','DANTZGT','DOOVBWT','GROOTGND'] #
    #stations=['HUIBGOT','MARSDND','ZOUTKPLZGT','ZUIDOLWOT']
    #stations=['NOORDWK2'] #,'NOORDWK10','NOORDWK50','NOORDWK70'],au
    #stations=['TERSLG4','TERSLG10'] #'TERSLG50','TERSLG70']
    #stations=['ROTTMPT3','ROTTMPT70'] #['ROTTMPT5', 'ROTTMPT15', 'ROTTMPT30','ROTTMPT50','ROTTMPT70']
    #stations = ['Helgoland', 'NOORDWK2', 'NOORDWK10', 'NOORDWK50', 'NOORDWK70', 'TERSLG4', 'TERSLG50', 'TERSLG70']
    # stations=['Helgoland','Sylt','SAmrum','Norderelbe','Norderney',
    #             'T36','T26','T41','T8','T2','T22','T5','T12','T11',
    #             'NOORDWK2','NOORDWK10','NOORDWK50','NOORDWK70',
    #             'TERSLG4','TERSLG50','TERSLG70',
    #             'ROTTMPT50','ROTTMPT70']
    #stations=['Norderney_2','Wesermuendung','JaBu','Bork']
    #stations=['BOCHTVWTM','BOOMKDP','DANTZGT','DOOVBWT','GROOTGND'] #
    #stations=['BOCHTVWTM','BOOMKDP','DANTZGT','DOOVBWT','GROOTGND','HUIBGOT','MARSDND','ROTTMPT3','ROTTMPT70','TERSLG10','TERSLG4','ZUIDOLWOT','ZOUTKPLZGT']
    stations=['SchlHolst2','Nney','Wesermuendung','JaBu','Bork','BOCHTVWTM','BOOMKDP','DANTZGT','DOOVBWT','GROOTGND','HUIBGOT',
              'MARSDND','ROTTMPT3','ROTTMPT70','TERSLG10','TERSLG4','ZUIDOLWOT','ZOUTKPLZGT']
    # stations=['BOCHTVWTM','BOOMKDP','DANTZGT','DOOVBWT','GROOTGND','HUIBGOT',
    #             'MARSDND','ROTTMPT3','ROTTMPT70','TERSLG10','TERSLG4','ZUIDOLWOT','ZOUTKPLZGT']
    # stations=['Nney','NOORDW','TERSLG','ROTTMP','BOCHTV','BOOMKD','DANTZG','DOOVBW','GROOTG','HUIBGO','MARSDN','ZOUTKP',
    #           'ZUIDOL']
    #rootpath='/home/onur/WORK/projects/GB/data/stations/individual'
    rootpath=home+'/stations/InterReg'
    dataset='InterReg'
    fstats=readstats(rootpath,stations,read)
    if merge:
        mergestats(rootpath,fstats)

    if isoBGC:
        isolate_convert_vars(rootpath,fstats,dataset)

def isolate_convert_vars(rootpath,fstats,dataset):
    print('Isolating and Converting Units')
    vars = {'DIN': 'mmolN/m3', 'NO3': 'mmolN/m3', 'NH4': 'mmolN/m3', 'KC':'1/m','SALT':'g/kg',
            'chl': 'mg/m3','DIP': 'mmolP/m3','Si': 'mmolSi/m3','zooC':'mmolC/m3'}

    # collect variables from files
    for station, fin in fstats.items():
        print(fin)
        fout = os.path.join(rootpath,dataset,station+'.nc')
        nc = netCDF4.Dataset(fin)
        ncv = nc.variables
        # procedures specific to stations
        print('Reading: ' + station,)

        varkeys = get_stat_varkeys(station)
        M = []; V = []; U = []; lats = []; lons = []
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

        print(': done. Writing:')

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
    print('Merging')
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
    for station,fin in fstats.items():
        print(fin)
        nc=netCDF4.Dataset(fin)
        ncv=nc.variables
        #procedures specific to stations
        print('Reading:'+station,)

        varkeys=get_stat_varkeys(station)
        M,V,U,lats,lons=getvar_fromstat(M,V,U,lats,lons,station,ncv,vars,varkeys,merge_statvar=True,keepZD=False,yearint=yearint,zint=zint)
        nc.close()
        print(': done.')

    #write in a ncdf file:
    Mv,tvec=mapon_container_tvec(M) #map all vars on along a container time vector
    #transform tvec to dates
    datevec = [datetime.datetime(refdate.year, refdate.month, refdate.day, 0, 0, 0) + datetime.timedelta(0, t_s) for t_s in tvec]
    create_ncfile_0D(fout_abs,datevec,-1,Mv,V,V,U,lats,lons,dims={'t':'time'},refdate=refdate,missing_value=-99,notify=True)

    #print 'Written:'+fout_abs

def get_stat_varkeys(station):
    if station== 'Bork':
        varkeys={'DOC':'DOC','Ntotal':'N','PO4':'DIP','Salz n.LF':'SALT','Seston':'TSM','Gluehverlust':'LOI','Pges':'P','Silikat':'Si','Nitrat':'NO3','Amon':'NH4','Nitrit':'NO2','Wassertemp.':'T'}
    elif station== 'Norderney':
        varkeys={'phyC':'','zooC':'','chl':'CPHL','DIN':'DIN','DIP':'PHOS','Si':'SLCA', 'NH4':'AMON','NO3':'NTRA'}
    elif station== 'Nney':
        varkeys={'CPHL':'chl','DIN':'DIN','Phosphat':'DIP','DOC':'DOC','Silikat':'Si', 'Ammonium':'NH4','Nitrat':'NO3','N':'N','P':'P','Temperatur':'T','Salzgeh. nach Leitfaehigkeit':'SALT','Seston Trockengew.':'TSM','Gluehverlust':'LOI','Secci-Sichttiefe':'zS'}
    elif station=='JaBu':
        varkeys={'Chlorophyll-a total Ethanol':'chl','DOC':'DOC','Gluehverlust Seston':'LOI','Abfiltrierbare_Stoffe (Seston)':'TSM','Salinitaet Sonde':'SALT',
                 'ortho P':'DIP','P total':'P','NH4N':'NH4','NO2N':'NO2','NO3N':'NO3','N total':'N',
                 'N total geloest':'DIN','SiO2':'SiO2','Si':'Si','TOC':'TOC','W_Sichttiefe':'zS','W_T_Wasser':'T'}
    elif station in ['JMA_Norderelbe', 'JMA_Suederpiep', 'JMA_Westerhever']:
        varkeys={'DIN':'NO3','DIP':'PO4'}
    elif station in ['Norderelbe','SAmrum']:
        #measurements are NH4-N, NO3-N, NO2-N and PO4-P
        varkeys={'phyC':'','zooC':'','chl':'Chlor','NH4':'71.4*Ammon','Si':'','NO3':'71.4*Nitra','DIN':'71.4*Ammon+71.4*Nitra+71.4*Nitri','DIP':'32.285*o-Pho'}
    elif station=='Sylt':
        varkeys={'phyC':'','zooC':'','chl':'chl','DIN':'nh4+ntri+ntra','DIP':'phos','NO3':'ntra','NH4':'nh4', 'Si':'si'}
    elif station=='Helgoland':
        varkeys={'phyC':'phytoplankton','zooC':'copepods','chl':'chl-a','DIN':'DIN','DIP':'PO4','Si':'SiO4','NO3':'NO3','NH4':'NH4'}
    elif station=='Wesermuendung':
        varkeys={'Chlorophyll-a gesamt Ethanol':'chl','DOC':'DOC','Gluehverlust d. abfiltrierbaren Stoffe':'LOI','Abfiltrierbare_Stoffe (Seston)':'TSM','W_Salzgehalt':'SALT',
                 'ortho P':'DIP','P total':'P','NH4N':'NH4','NO2N':'NO2','NO3N':'NO3','N total':'N',
                 'N total geloest':'DIN','SiO2':'SiO2','Si':'Si','TOC':'TOC','W_Sichttiefe':'zS','Temperatur Wasser':'T'}
    elif station in ['T36','T26','T41','T8','T2','T22','T5','T12','T11']:
        varkeys={'phyC':'','zooC':'','chl':'chl','DIN':'DIN','DIP':'DIP'}
    elif station[0:6] in ['NOORDW', 'TERSLG','ROTTMP','BOCHTV','BOOMKD','DANTZG','DOOVBW','GROOTG','HUIBGO','MARSDN','ZOUTKP','ZUIDOL']:
        # measurements are NH4-N, NO3-N, NO2-N and PO4-P
        varkeys={'phyC':'','zooC':'','chl':'CHL','DIN':'71.4*NH4+71.4*NO3','DIP':'32.285*PO4', 'KC':'KC','SALT':'SALT', #+71*NO2
                 'SPM':'0.001*SPM','NH4':'71.4*NH4','NO3':'71.4*NO3','NO2':'71.4*NO2','PO4':'32.285*PO4',
                 'Si':'35.6*SiO2', 'SecchiDepth':'SD', 'DOC':'83.26*DOC','POC':'83.26*POC'}
        #measurements are NH4,NO3,NO2,PO4
        #varkeys = {'phyC': '', 'zooC': '', 'chl': 'CHL', 'DIN': '55.43*NH4+16.13*NO3+21.74*NO2', 'DIP':'10.56*PO4'}
    else:
        varkeys={}
    return varkeys

def getvar_fromstat(M,V,U,lat,lon,station,ncv,vars,varkeys,merge_statvar=False,keepZD=False,yearint=[0,0],zint=[0,0]):

    for globvar in vars.keys():
    #for globvar,statvar in varkeys.items():
        print(globvar,)
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
                print('no lat. variable found. Skipping.')

            if 'lon' in ncv.keys():
                lon.append(ncv['lon'][:])
            elif 'longitude' in ncv.keys():
                lon.append(ncv['longitude'][:])
            else:
                print('no lat. variable found. Skipping.')

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
        #v=np.squeeze(ncv[realvar][:])
        v=ncv[realvar][:]
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
        print('(not found !)',)
        vf=[]
    return vf

def read_bsh(fin,fout,dims,lat,lon,station='-'):

    print(' finding bottom depth..',)
    #maxz=np.nan
    maxz = get_botdepth(lon[0], lat[0], 'tree')
    print(' found: %sN, %sE: %.1f m' % (lat[0], lon[0], maxz))

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
    print(str(numrow)+' lines,'+str(numcolO)+' columns. ',)

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
        print(var,)

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

    print(' finding bottom depth..',)
    #maxz=np.nan
    maxz = get_botdepth(lon[0], lat[0], 'tree')
    print(' found: %sN, %sE: %.1f m' % (lat[0], lon[0], maxz))

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
    print(str(numrow)+' lines,'+str(numcolO)+' columns. ',)

    #read time:
    secs=np.ndarray((numrow,1)) #seconds
    for  r in np.arange(numrow):
        years=int(RDraw[r].split(',')[0])
        months=int(RDraw[r].split(',')[1])
        deltadate=datetime.date(years,months,15) - datetime.date(refdate.year,refdate.month,refdate.day)
        secs[r]=86400*deltadate.days +deltadate.seconds

    #read each variable into a time-value matrix
    for varno,var in enumerate(V):
        print(var,)

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

    print(' finding bottom depth..',)
    #maxz=np.nan
    maxz = get_botdepth(lon[0], lat[0], 'tree')
    print(' found: %sN, %sE: %.1f m' % (lat[0], lon[0], maxz))


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
    print(str(numrow)+' lines,'+str(numcolO)+' columns. ',)

    #read time:
    secs=np.ndarray((numrow,1)) #seconds
    for  r in np.arange(numrow):
        years=int(RDraw[r].split(',')[0])
        months=int(RDraw[r].split(',')[1])
        deltadate=datetime.date(years,months,15) - datetime.date(refdate.year,refdate.month,refdate.day)
        secs[r]=86400*deltadate.days +deltadate.seconds

    #read each variable into a time-value matrix
    for varno,var in enumerate(V):
        print(var,)

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
          #'DOC':'DOC',
          'NH4':'NH4',
          #'NO2':'NO2',
          'NO3':'NO3',
          #'N': 'N',
          #'NO3NO2': 'NO3NO2',
          'PO4':'PO4',
          #'P': 'P',
          #'POC':'POC',
          #'SD': 'secchi_depth',
          'SiO2':'SiO2',
          #'SPM':'concentration_of_suspended_matter_in_water'
          'SALT':'sea_water_salinity',
          'KC': 'E'
          }

    lon =-999.0; lat =-999.0 #these are retrieved from the raw files
    M=[];U=[];V=[]

    print('variable:',)
    for varshort,varlong in vars.items():
        print(varshort,)

        fin = fin_base + '-' + varshort + '.nc'
        if not os.path.isfile(fin):
            print('(not available)',)
            continue
        try:
            nc = netCDF4.Dataset(fin)
        except:
            print('(not available)',)
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
            print('sum(vali)=0')

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
            print('')
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
            if abs(lon-newlon)>0.05 or abs(lat-newlat)>0.05:
                raise(ValueError('previous lon-lat pair (%s-%s) do not match with the new pair(%s-%s)'%(lon,lat,newlon,newlat)))

        nc.close()
    print('')
    #Mv,tvec_s=mapon_container_tvec(M) #map all vars on along a container time vector
    Mv, tvec_s, zvec = mapon_container_tzvec(M)  # map all vars along a container time/z vector
    tvec=[datetime.datetime(refdate.year,refdate.month,refdate.day,0,0,0)+datetime.timedelta(0,t_s) for t_s in tvec_s]
    create_ncfile(fout, lon, lat, Mv, V, V, U, dims, tvec,zvec,
                  climatology=False,refdate=refdate, missing_value=-99, notify=True)

    print(' appending bottom depth..',)
    # maxz=np.nan
    maxz = get_botdepth(lon, lat, 'tree')
    print(' found: %sN, %sE: %.1f m' % (lat, lon, maxz))
    addncatt(fout, {'station': station, 'bottom_depth': maxz})

    print('file written:'+fout)

def read_dutchcsv(rootpath,fin,fin2,fout,dims,station='-'):
    varkeys={
        'CHLFa':'chl',
        'DIN':'DIN',
        'DOC':'DOC',
        'PO4':'DIP',
        'E':'Kc',
        'NH4':'NH4',
        'NO2':'NO2',
        'NO3':'NO3',
        'POC':'POC',
        'SALNTT':'SALT',
        'SiO2':'Si',
        'T':'T',
        'TN':'N',
        'TP':'P',
        'ZICHT':'Zs',
        'ZS':'TSM'
            }
    
    with open(fin2, 'rb') as csvfile:
        RD=csvfile.readlines()
        
    headers=RD[0].split(';')
    del RD[0]
    numrows=len(RD)
    numcols=len(headers)
    
    varst=np.array([RD[r].split(';')[0]==station for r in range(numrows)])
    varst2=np.where(varst)[0]
    r=varst2[0]
    
    lat=np.array([RD[r].split(';')[2]])
    lon=np.array([RD[r].split(';')[3]])
    
    print ' finding bottom depth..',
    maxz = get_botdepth(lon[0], lat[0], 'tree')
    print ' found: %sN, %sE: %.1f m' % (lat[0], lon[0], maxz)

    with open(fin, 'rb') as csvfile:
        RD=csvfile.readlines()
        
    headers=RD[0].split(';')
    del RD[0]
    numrows=len(RD)
    numcols=len(headers)
    allvars=headers[5:] 
    
    for r in range(len(allvars)):
        allvars[r]=allvars[r].replace('\r\n','')
    
    for r in range(len(headers)):
        headers[r]=headers[r].replace('\r\n','')
    
    uvars=allvars   
    numvars=0
    TSM_switch=False
    LOI_switch=False
    
    for varno,var in enumerate(allvars):
        if var in varkeys:
            varaux=varkeys[var]
            if varaux=='TSM':
                print 'TSM_switch==True'
                TSM_switch=True
                TSM_index=numvars
            elif varaux=='LOI':
                print 'LOI_switch==True'
                LOI_switch=True
                LOI_index=numvars
            numvars+=1
        else:
            uvars[varno]='wzf'
            
    while 'wzf' in uvars: uvars.remove('wzf')
    
    if TSM_switch==True and LOI_switch==True:
        numvars=numvars+2
        
    M=[None]*numvars #chl,amon,din,ntra,ntri,phos,ntot,ptot,slca,doc,tsm,loi,zs
    U=[None]*numvars
    V=[None]*numvars
    D=[None]*numvars

    conversion_factor={
        'CHLFa':1,
        'DIN':7.14e1,
        'DOC':8.33e1,
        'DON':7.14e1,
        'DOP':3.23e1,
        'E':1,
        'NH4':7.14e1,
        'NO2':7.14e1,
        'NO3':7.14e1,
        'O2':3.13e1,
        'PO4':3.23e1,
        'POC':8.33e1,
        'SALNTT':1,
        'SiO2':3.56e1,
        'T':1,
        'TN':7.14e1,
        'TP':3.23e1,
        'ZICHT':0.1,
        'ZS':1
            }
    
    units={
        'CHLFa':'mg/m3',
        'DIN':'mmol N/m3',
        'DOC':'mmol C/m3',
        'DON':'mmol N/m3',
        'DOP':'mmol P/m3',
        'E':'1/m',
        'NH4':'mmol N/m3',
        'NO2':'mmol N/m3',
        'NO3':'mmol N/m3',
        'O2':'g/m3',
        'PO4':'mmol P/m3',
        'POC':'mmol C/m3',
        'SALNTT':'PSU',
        'SiO2':'mmol Si/m3',
        'T':'degC',
        'TN':'mmol N/m3',
        'TP':'mmol P/m3',
        'ZICHT':'m',
        'ZS':'g/m3'
            }

    depth_key={
        'BODM':maxz-1,
        'HALVWTKL':maxz*0.5,
        'NAP':0,
        'NVT': np.nan,
        'SPRONGLG':maxz*0.5,
        'WATSGL':3
            }

    secs=np.ndarray((numrows,1)) #seconds
    length_station=-1
    for r in np.arange(numrows):
        if RD[r].split(';')[0]==station:
            length_station+=1
       
    for varno,var in enumerate(uvars):
        if var in varkeys:
            varaux=varkeys[var]
            V[varno]=varaux
            U[varno]=units[var]
        
    
    varr=np.array([RD[r].split(';')[0]==station for r in range(numrows)])
    varr2=np.where(varr)[0]
    
    for  r in varr2:
        dateaux=datetime.datetime.strptime(RD[r].split(';')[2],'%d-%b-%y').date()
        deltadate=dateaux-datetime.date(refdate.year,refdate.month,refdate.day)
        secs[r]=86400*deltadate.days+deltadate.seconds
    
    rr=len(varr2)
    
    for varno,var in enumerate(uvars):
        vari=np.array([headers[r]==var for r in range(len(headers))])
        vari2=np.asscalar(np.where(vari)[0])
        M[varno]=np.ndarray((rr,3))
        D[varno]=np.ndarray((rr,1))
        conv_fac2=conversion_factor[var]
        #print 'bla 744: ',vari2,var,varkeys[var],conv_fac2
        kk=0
        for r in varr2:
            M[varno][kk,0] = secs[r]
            M[varno][kk,1] = depth_key[RD[r].split(';')[3]]
            # D[varno][kk,0] = depth_key[RD[r].split(';')[3]]
            val=RD[r].split(';')[vari2]
            
            val=val.replace('\n','')
                
            if val=='' or '\n' in val:
                M[varno][kk,2]=np.nan
            else:
                try:
                    if val[0]=='<':
                        M[varno][kk,2]=0.0
                    else:
                        M[varno][kk,2]=float(val)*conv_fac2
                except:
                        #print 'bla 757!'  
                        M[varno][kk,:]=np.nan
            kk+=1
        M[varno]=M[varno][np.invert(np.isnan(M[varno][:,1])),:]
        M[varno]=M[varno][np.invert(np.isnan(M[varno][:,2])),:]
                    
    if TSM_switch==True and LOI_switch==True:
        
        M[-2]=np.ndarray((numrows,3))
        V[-2]='PIM'
        U[-2]='g/m3'
        
        for r in range(numrows):
            try:
                M[-2][r,0]=M[TSM_index][r,0]
                M[-2][r,1]=M[TSM_index][r,1]
                M[-2][r,2]=M[TSM_index][r,2]*(1-0.01*M[LOI_index][r,2])
            except:
                M[-2][r,:]=np.nan

        M[-2]=M[-2][np.invert(np.isnan(M[-2][:,1])),:]
        M[-2]=M[-2][np.invert(np.isnan(M[-2][:,2])),:]
        
        M[-1]=np.ndarray((numrows,3))
        V[-1]='POM'
        U[-1]='g/m3'
        
        for r in range(numrows):
            try:
                M[-1][r,0]=M[TSM_index][r,0]
                M[-1][r,1]=M[TSM_index][r,1]
                M[-1][r,2]=M[TSM_index][r,2]*0.01*M[LOI_index][r,2]
            except:
                M[-1][r,:]=np.nan

        M[-1]=M[-1][np.invert(np.isnan(M[-1][:,1])),:]
        M[-1]=M[-1][np.invert(np.isnan(M[-1][:,2])),:]
    
    # Mv,tvec_s=mapon_container_tvec(M) #map all vars on along a container time vector
    Mv, tvec_s, zvec = mapon_container_tzvec(M)  # map all vars along a container time/z vector
    tvec=[datetime.datetime(refdate.year,refdate.month,refdate.day,0,0,0)+datetime.timedelta(0,t_s) for t_s in tvec_s]
    #create_ncfile(fout,lon,lat,tvec,-1,Mv,V,V,U,dims,refdate=refdate,missing_value=-99,notify=True)
    create_ncfile(fout, lon, lat, Mv, V, V, U, dims, tvec, zvec, climatology=False,
                  refdate=refdate, missing_value=-99, notify=True)
    #print 'file written:'+fout
    addncatt(fout, {'station': station, 'bottom_depth': maxz})
    
def read_SchlHolst(fin,fout,dims,lat,lon,station='-'):
    print(' finding bottom depth..',)
    # maxz=np.nan
    maxz = get_botdepth(lon[0], lat[0], 'tree')
    print(' found: %sN, %sE: %.1f m' % (lat[0], lon[0], maxz))

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
        print(var +':'+str(numrow)+' lines,'+str(numcolO)+' columns. ',)

        U[varno]=RD[0].split(',')[11]
        V[varno]=var
        M[varno]=np.ndarray((numrow,2)) #seconds,var
        for  r in range(numrow):
            M[varno][r,0] = secsV[r]
            M[varno][r,1]=float(RD[r].split(';')[10]) if RD[r].split(';')[10]!='' else np.nan

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
    print(' finding bottom depth..',)
    # maxz=np.nan
    maxz = get_botdepth(lon[0], lat[0], 'tree')
    print(' found: %sN, %sE: %.1f m' % (lat[0], lon[0], maxz))

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
    print('CHL:'+ str(numrow)+' lines,'+str(numcolO)+' columns. ',)

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
        print(var +':'+str(numrow)+' lines,'+str(numcolO)+' columns. ',)

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

def read_Nney(rootpath,fin,fin2,fin3,fout,dims,lat,lon,station='-'):
    print ' finding bottom depth..',
    # maxz=np.nan
    maxz = get_botdepth(lon[0], lat[0], 'tree')
    varkeys=get_stat_varkeys('Nney')
    print ' found: %sN, %sE: %.1f m' % (lat[0], lon[0], maxz)

    M=[None]*16 #chl,amon,din,ntra,ntri,phos,ntot,ptot,slca,doc,tsm,loi,zs
    U=[None]*16
    V=[None]*16

    #read chl
    with open(fin2, 'rb') as csvfile:
        RD=csvfile.readlines()

    headers=RD[0].split(',')
    del RD[0]
    del RD[0]
    numrow=len(RD)

    numcolO=len(headers)
    print 'CHL:'+ str(numrow)+' lines,'+str(numcolO)+' columns. ',

    M[0]=np.ndarray((numrow,3)) #seconds,chl
    for  r in np.arange(0,numrow):
        years=int(RD[r].split(';')[1].split('.')[2])
        months=int(RD[r].split(';')[1].split('.')[1])
        days=int(RD[r].split(';')[1].split('.')[0])
        hours=int(RD[r].split(';')[2].split(':')[0])
        mins=int(RD[r].split(';')[2].split(':')[1])
        deltadate=datetime.datetime(years,months,days,hours,mins)- refdate
        M[0][r,0] = 86400*deltadate.days +deltadate.seconds
        M[0][r,2]=float(RD[r].split(';')[4]) if RD[r].split(';')[4]!='' else np.nan

    #get rid of empty lines
    M[0][:,1]=0
    M[0]=M[0][np.invert(np.isnan(M[0][:,2])),:]
    U[0]='mg/m3'
    V[0]='chl'
    
    
    #read nut
    with open(fin, 'rb') as csvfile:
        RDraw=csvfile.readlines()

    headers=RDraw[0].split(';')
    varar=RDraw[1].split(';')
    del RDraw[0]
    del RDraw[0]
    numcolO=len(headers)
    numrow=len(RDraw)
    #print str(numrow)+' lines,'+str(numcolO)+' columns. ',

    #read all variables:
    #varar=np.chararray((numrow,1), itemsize=5)
    secs=np.ndarray((numrow,1)) #seconds
    for  r in np.arange(0,numrow):
        #print 'bla 786: ',RDraw[r].split(';')[11]
        #varar[r,0]=RDraw[r].split(';')[11]
        years=int(RDraw[r].split(';')[1].split('.')[2])
        months=int(RDraw[r].split(';')[1].split('.')[1])
        days=int(RDraw[r].split(';')[1].split('.')[0])
        hours=int(RDraw[r].split(';')[2].split(':')[0])
        mins=int(RDraw[r].split(';')[2].split(':')[1])
        deltadate=datetime.datetime(years,months,days,hours,mins)- refdate
        secs[r]=86400*deltadate.days +deltadate.seconds

    #read each variable into a time-value matrix
    uvars=np.unique(varar)
    uvars=np.delete(uvars,np.where(uvars==''))
    kk=0
    for varno,var in enumerate(uvars):
        #vari=np.array([varar[r][0]==var for r in range(len(varar))])
        vari=np.array([varar[r]==var for r in range(len(varar))])
        vari2=np.where(vari)[0]
        if len(vari2)==1:
            vari2=vari2[0]
        else:
            print 'bla 1038!'
        secsV=secs
        #RD=[None]*len(secsV)
        
        #for r0,r1 in enumerate(np.where(vari)[0]):
            #RD[r0]=RDraw[r1] #i.e., skip the header
        #print 'bla 813: ',len(RDraw)
        #for r0,r1 in enumerate(np.where(vari)[0]):
            #for r in range(len(RDraw)):
                #print 'bla 817: ',RDraw[r].split(';')[r1],r1,np.shape(RDraw[r].split(';'))
                #RD[r]=RDraw[r].split(';')[r1]
                #print 'bla 818: ',RD[r]
            
        numrow=len(RDraw)
        
        #print 'bla 822: ',varkeys[var.split('[')[0].split(' ')[0]]
        varn=var.split('[')[0].split(' ')[0]
        if varn in varkeys:
            kk=kk+1
            varaux = varkeys[varn]
            #U[kk]=var.split('[')[1][0:-1]
            U[kk]=var.split('[')[1].split(']')[0]
            V[kk]=varaux
            M[kk]=np.ndarray((numrow,3)) #seconds,var
            print varaux +':'+str(numrow)+' lines,'+str(numcolO)+' columns. '
            
            for  r in range(numrow):
                M[kk][r,0] = secsV[r]
                val=RDraw[r].split(';')[vari2]
                if val=='' or '\n' in val:
                    M[kk][r,2]=np.nan
                else:
                    try:
                        if val[0]=='<':
                            M[kk][r,2]=0.0
                        else:
                            M[kk][r,2]=float(val)
                    except:
                            print 'bla 854!'
            #get rid of empty lines
            M[kk][:,1]=0
            M[kk]=M[kk][np.invert(np.isnan(M[kk][:,2])),:]

    #read nut
    with open(fin3, 'rb') as csvfile:
        RDraw2=csvfile.readlines()

    headers=RDraw2[0].split(';')
    varar=RDraw2[1].split(';')
    del RDraw2[0]
    del RDraw2[0]
    numcolO=len(headers)
    numrow=len(RDraw2)
    
    secs2=np.ndarray((numrow,1)) #seconds
    for  r in np.arange(0,numrow):
        years=int(RDraw2[r].split(';')[1].split('.')[2])
        months=int(RDraw2[r].split(';')[1].split('.')[1])
        days=int(RDraw2[r].split(';')[1].split('.')[0])
        hours=int(RDraw2[r].split(';')[2].split(':')[0])
        mins=int(RDraw2[r].split(';')[2].split(':')[1])
        deltadate=datetime.datetime(years,months,days,hours,mins)- refdate
        secs2[r]=86400*deltadate.days +deltadate.seconds

    uvars=np.unique(varar)
    uvars=np.delete(uvars,np.where(uvars==''))
    secsV=secs2
    for varno,var in enumerate(uvars):
        vari=np.array([varar[r]==var for r in range(len(varar))])
        vari2=np.where(vari)[0]
        if len(vari2)==1:
            vari2=vari2[0]
        else:
            print 'bla 1110!'
        numrow=len(RDraw2)
        
        varn=var.split('[')[0][0:-1]
        if varn in varkeys:
            kk=kk+1
            varaux = varkeys[varn]
            U[kk]=var.split('[')[1].split(']')[0]
            V[kk]=varaux
            M[kk]=np.ndarray((numrow,3)) #seconds,var
            print varaux +':'+str(numrow)+' lines,'+str(numcolO)+' columns. '
            if varaux=='TSM':
                TSM_index=kk
                TSM_switch=True
            elif varaux=='LOI':
                LOI_index=kk
                LOI_switch=True
                
            for  r in range(numrow):
                M[kk][r,0] = secsV[r]
                val=RDraw2[r].split(';')[vari2]
                if val=='' or '\n' in val:
                    M[kk][r,2]=np.nan
                else:
                    try:
                        if val[0]=='<':
                            M[kk][r,2]=0.0
                        else:
                            M[kk][r,2]=float(val)
                    except:
                            print 'bla 915!'
            #get rid of empty lines
            M[kk][:,1]=0
            M[kk]=M[kk][np.invert(np.isnan(M[kk][:,2])),:]

    if TSM_switch==True and LOI_switch==True:
        M[kk+1]=np.ndarray((numrow,3))
        V[kk+1]='PIM'
        U[kk+1]='g/m3'
        
        for r in range(numrow):
            try:
                M[kk+1][r,0]=M[TSM_index][r,0]
                M[kk+1][r,2]=M[TSM_index][r,2]*(1-0.01*M[LOI_index][r,2])
            except:
                M[kk+1][r,:]=np.nan

        M[kk+1][:,1] = 0
        M[kk+1]=M[kk+1][np.invert(np.isnan(M[kk+1][:,2])),:]

        M[kk+2]=np.ndarray((numrow,3))
        V[kk+2]='POM'
        U[kk+2]='g/m3'
        
        for r in range(numrow):
            try:
                M[kk+2][r,0]=M[TSM_index][r,0]
                M[kk+2][r,2]=M[TSM_index][r,2]*0.01*M[LOI_index][r,2]
            except:
                M[kk+2][r,:]=np.nan

        M[kk+1][:,2] = 0
        M[kk+2]=M[kk+2][np.invert(np.isnan(M[kk+2][:,2])),:]
    
    #write in a ncdf file:
    # Mv,tvec_s=mapon_container_tvec(M) #map all vars on along a container time vector
    Mv, tvec_s, zvec = mapon_container_tzvec(M)  # map all vars along a container time/z vector
    tvec=[datetime.datetime(refdate.year,refdate.month,refdate.day,0,0,0)+datetime.timedelta(0,t_s) for t_s in tvec_s]
    #create_ncfile(fout,lon,lat,tvec,-1,Mv,V,V,U,dims,refdate=refdate,missing_value=-99,notify=True)
    create_ncfile(fout, lon, lat, Mv, V, V, U, dims, tvec, zvec, climatology=False,
                  refdate=refdate, missing_value=-99, notify=True)
    #print 'file written:'+fout
    addncatt(fout, {'station': station, 'bottom_depth': maxz})

def read_JaBu(rootpath,fin,fout,dims,lat,lon,station='-'):
    print ' finding bottom depth..',
    maxz = get_botdepth(lon[0], lat[0], 'tree')
    varkeys=get_stat_varkeys('JaBu')
    print ' found: %sN, %sE: %.1f m' % (lat[0], lon[0], maxz)

    numvars=len(varkeys)+2
    M=[None]*numvars
    U=[None]*numvars
    V=[None]*numvars

    #read chl
    with open(fin, 'rb') as csvfile:
        RD=csvfile.readlines()
    headers=RD[0].split(';')
    del RD[0]
    allunits=RD[0].split(';')
    del RD[0]
    numrow=len(RD)

    numcolO=len(headers)

    secs=np.ndarray((numrow,1)) #seconds
    for  r in np.arange(numrow):
        years=int(RD[r].split(';')[2].split('.')[2])
        months=int(RD[r].split(';')[2].split('.')[1])
        days=int(RD[r].split(';')[2].split('.')[0])
        #hours=int(RDraw[r].split(',')[2].split(':')[0])
        #mins=int(RDraw[r].split(',')[2].split(':')[1])
        #deltadate=datetime.datetime(years,months,days,hours,mins)- refdate
        deltadate=datetime.date(years,months,days)- datetime.date(refdate.year,refdate.month,refdate.day)
        secs[r]=86400*deltadate.days +deltadate.seconds

    #read each variable into a time-value matrix
    uvars=np.unique(headers)
    kk=0
    LOI_switch=False
    TSM_switch=False
    for varno,var in enumerate(uvars):
        if var in varkeys:
            varaux=varkeys[var]
            if varaux=='TSM':
                print 'TSM_switch==True'
                TSM_switch=True
                TSM_index=kk
            elif varaux=='LOI':
                print 'LOI_switch==True'
                LOI_switch=True
                LOI_index=kk
                
            vari=np.array([headers[r]==var for r in range(len(headers))])
            vari2=np.where(vari)[0]
            unitaux=[allunits[i] for i in vari2]
            if len(unitaux)>1:
                unitout=unitaux[1]
            else:
                unitout=unitaux[0]
            if unitout==' [umol/l]':
                unitout='umol/l'
            V[kk]=varaux
            U[kk]=unitout
            M[kk]=np.ndarray((numrow,3))
            M[kk][:,1]=0
            for r in range(numrow):
                M[kk][r,0] = secs[r]
                val=RD[r].split(';')[vari2[-1]]
                if val=='' or '\n' in val:
                    M[kk][r,2]=np.nan
                else:
                    try:
                        if val[0]=='<':
                            M[kk][r,2]=0.0
                        else:
                            M[kk][r,2]=float(val)
                    except:
                            print 'bla 932!'
            M[kk]=M[kk][np.invert(np.isnan(M[kk][:,2])),:]
            kk=kk+1
            
    M[kk]=np.ndarray((numrow,3))
    M[kk][:,1]=0
    V[kk]='PIM'
    U[kk]='g/m3'
    
    for r in range(numrow):
        try:
            M[kk][r,0]=M[TSM_index][r,0]
            M[kk][r,2]=M[TSM_index][r,2]*(1-0.01*M[LOI_index][r,2])
        except:
            M[kk][r,:]=np.nan
            
    M[kk]=M[kk][np.invert(np.isnan(M[kk][:,2])),:]

    
    M[kk+1]=np.ndarray((numrow,3))
    M[kk+1][:,1]=0
    V[kk+1]='POM'
    U[kk+1]='g/m3'
    
    for r in range(numrow):
        try:
            M[kk+1][r,0]=M[TSM_index][r,0]
            M[kk+1][r,2]=M[TSM_index][r,2]*0.01*M[LOI_index][r,2]
        except:
            M[kk+1][r,:]=np.nan
            
    M[kk+1]=M[kk+1][np.invert(np.isnan(M[kk+1][:,2])),:]
    
    # Mv,tvec_s=mapon_container_tvec(M) #map all vars on along a container time vector
    Mv, tvec_s, zvec = mapon_container_tzvec(M)  # map all vars along a container time/z vector
    tvec=[datetime.datetime(refdate.year,refdate.month,refdate.day,0,0,0)+datetime.timedelta(0,t_s) for t_s in tvec_s]
    #create_ncfile(fout,lon,lat,tvec,-1,Mv,V,V,U,dims,refdate=refdate,missing_value=-99,notify=True)
    create_ncfile(fout, lon, lat, Mv, V, V, U, dims, tvec, zvec, climatology=False,
                  refdate=refdate, missing_value=-99, notify=True)
    #print 'file written:'+fout
    addncatt(fout, {'station': station, 'bottom_depth': maxz})

def read_Bork(rootpath,fin,fout,dims,lat,lon,station='-'):
    print ' finding bottom depth..',
    maxz = get_botdepth(lon[0], lat[0], 'tree')
    varkeys=get_stat_varkeys('Bork')
    
    print ' found: %sN, %sE: %.1f m' % (lat[0], lon[0], maxz)

    M=[None]*14 #chl,amon,din,ntra,ntri,phos,ntot,ptot,slca,doc,tsm,loi,zs
    U=[None]*14
    V=[None]*14

    with open(fin, 'rb') as csvfile:
        RD=csvfile.readlines()
    headers=RD[0].split(';')
    del RD[0]
    numrow=len(RD)

    numcolO=len(headers)
    
    secs=np.ndarray((numrow,1)) #seconds
    for  r in np.arange(numrow):
        years=int(RD[r].split(';')[1].split('.')[2])
        months=int(RD[r].split(';')[1].split('.')[1])
        days=int(RD[r].split(';')[1].split('.')[0])
        #hours=int(RDraw[r].split(',')[2].split(':')[0])
        #mins=int(RDraw[r].split(',')[2].split(':')[1])
        #deltadate=datetime.datetime(years,months,days,hours,mins)- refdate
        deltadate=datetime.date(years,months,days)- datetime.date(refdate.year,refdate.month,refdate.day)
        secs[r]=86400*deltadate.days +deltadate.seconds

    uvarsraw=np.unique(headers)
    #uvarsraw=np.delete(uvarsraw,np.where(uvarsraw==''))
    #uvarsraw=np.delete(uvarsraw,0)
    uvars=len(uvarsraw)*[None]
    allunits=len(uvarsraw)*[None]
    for r in range(len(uvarsraw)):
        varaux=uvarsraw[r].split('[')[0][0:-1]
        uvars[r]=varaux
        try:
            unitaux=uvarsraw[r].split('[')[1][0:-1]
        except:
            unitaux=''
        if '\n' in unitaux:
            unitaux=''
        allunits[r]=unitaux
        
    kk=0
    LOI_switch=False
    TSM_switch=False
    for varno,var in enumerate(uvars):
        if var in varkeys:
            varaux=varkeys[var]
            if varaux=='TSM':
                print 'TSM_switch==True'
                TSM_switch=True
                TSM_index=kk
            elif varaux=='LOI':
                print 'LOI_switch==True'
                LOI_switch=True
                LOI_index=kk
    
            #vari=np.array([headers[r]==var for r in range(len(headers))])
            vari=np.array([var in headers[r] for r in range(len(headers))])
            vari2=np.where(vari)[0]
            unitout=allunits[varno]
            
            V[kk]=varaux
            U[kk]=unitout
            M[kk]=np.ndarray((numrow,3))
            M[kk][:,1]=0

            for r in range(numrow):
                M[kk][r,0] = secs[r]
                val=RD[r].split(';')[vari2[-1]]
                if val=='' or '\n' in val or val=='not measured':
                    M[kk][r,2]=np.nan
                else:
                    try:
                        if val[0]=='<':
                            M[kk][r,2]=0.0
                        else:
                            M[kk][r,2]=float(val)
                    except:
                            print 'bla 932: ',val
            M[kk]=M[kk][np.invert(np.isnan(M[kk][:,2])),:]
            kk=kk+1
            
    M[kk]=np.ndarray((numrow,3))
    M[kk][:,1]=0
    V[kk]='PIM'
    U[kk]='g/m3'
    
    for r in range(numrow):
        try:
            M[kk][r,0]=M[TSM_index][r,0]
            M[kk][r,2]=M[TSM_index][r,2]*(1-0.01*M[LOI_index][r,2])
        except:
            M[kk][r,:]=np.nan
            
    M[kk]=M[kk][np.invert(np.isnan(M[kk][:,2])),:]
    
    kk=kk+1
    M[kk]=np.ndarray((numrow,3))
    M[kk][:,1]=0
    V[kk]='P0M'
    U[kk]='g/m3'
    
    for r in range(numrow):
        try:
            M[kk][r,0]=M[TSM_index][r,0]
            M[kk][r,2]=M[TSM_index][r,2]*0.01*M[LOI_index][r,2]
        except:
            M[kk][r,:]=np.nan
            
    M[kk]=M[kk][np.invert(np.isnan(M[kk][:,2])),:]
            
    # Mv,tvec_s=mapon_container_tvec(M) #map all vars on along a container time vector
    Mv, tvec_s, zvec = mapon_container_tzvec(M)  # map all vars along a container time/z vector
    tvec=[datetime.datetime(refdate.year,refdate.month,refdate.day,0,0,0)+datetime.timedelta(0,t_s) for t_s in tvec_s]
    #create_ncfile(fout,lon,lat,tvec,-1,Mv,V,V,U,dims,refdate=refdate,missing_value=-99,notify=True)
    create_ncfile(fout, lon, lat, Mv, V, V, U, dims, tvec, zvec, climatology=False,
                  refdate=refdate, missing_value=-99, notify=True)
    #print 'file written:'+fout
    addncatt(fout, {'station': station, 'bottom_depth': maxz})

def read_Wesermuendung(rootpath,fin,fout,dims,lat,lon,station='-'):
    varkeys=get_stat_varkeys('Wesermuendung')
    
    numvars=len(varkeys)+1
    M=[None]*numvars #chl,amon,din,ntra,ntri,phos,ntot,ptot,slca,doc,tsm,loi,zs
    U=[None]*numvars
    V=[None]*numvars

    with open(fin, 'rb') as csvfile:
        RD=csvfile.readlines()
    headers=RD[0].split(';')
    del RD[0]
    allunits=RD[0].split(';')
    del RD[0]
    numrow=len(RD)

    numcolO=len(headers)
    
    secs=np.ndarray((numrow,1)) #seconds
    substatname=[None]*numrow
    for  r in np.arange(numrow):
        years=int(RD[r].split(';')[2].split('.')[2])
        months=int(RD[r].split(';')[2].split('.')[1])
        days=int(RD[r].split(';')[2].split('.')[0])
        #hours=int(RDraw[r].split(',')[2].split(':')[0])
        #mins=int(RDraw[r].split(',')[2].split(':')[1])
        #deltadate=datetime.datetime(years,months,days,hours,mins)- refdate
        deltadate=datetime.date(years,months,days)- datetime.date(refdate.year,refdate.month,refdate.day)
        secs[r]=86400*deltadate.days +deltadate.seconds
        substatname[r]=RD[r].split(';')[0]
    
    uvars=np.unique(headers)
    statid=['WeMu_W2','WeMu_W2','WuKu_W1']
    for r in range(len(substatname)):
        if substatname[r] in ['WeMu_W2 (Auto)','WeMu_W_2 (Schiff)','WeMu_W_2']:
            substatname[r]=statid[0]
        elif substatname[r]=='WeMu_W_1 (km  085 Aussenweser)':
            substatname[r]=statid[1]
        elif substatname[r]=='WuKu_W_1':
            substatname[r]=statid[2]
        else:
            substatname[r]='wzf'
    
    #while 'wzf' in substatname: substatname.remove('wzf')
    
    for si in range(0,3):
        print 'Substation: '+statid[si]
        print ' finding bottom depth..',
        maxz = get_botdepth(lon[si], lat[si], 'tree')
        
        print ' found: %sN, %sE: %.1f m' % (lat[si], lon[si], maxz)
        
        M2=[None]*numvars
            
        kk=0
        LOI_switch=False
        TSM_switch=False
        for varno,var in enumerate(uvars):
            if var in varkeys:
                varaux=varkeys[var]
                if varaux=='TSM':
                    print 'TSM_switch==True'
                    TSM_switch=True
                    TSM_index=kk
                elif varaux=='LOI':
                    print 'LOI_switch==True'
                    LOI_switch=True
                    LOI_index=kk
                    
                vari=np.array([headers[r]==var for r in range(len(headers))])
                vari2=np.where(vari)[0]
                unitaux=[allunits[i] for i in vari2]
                if len(unitaux)>1:
                    unitout=unitaux[1]
                else:
                    unitout=unitaux[0]
                if unitout==' [umol/l]':
                    unitout='umol/l'
                
                V[kk]=varaux
                U[kk]=unitout
                M2[kk]=np.ndarray((numrow,3))
                M2[kk][:,1]=0

                for r in range(numrow):
                    M2[kk][r,0] = secs[r]
                    if substatname[r]==statid[si]:
                        val=RD[r].split(';')[vari2[-1]]
                        if val=='' or '\n' in val:
                            M2[kk][r,2]=np.nan
                        else:
                            try:
                                if val[0]=='<':
                                    M2[kk][r,2]=0.0
                                else:
                                    M2[kk][r,2]=float(val)
                                    if kk==14 and float(val)>40:
                                        print 'bla'
                            except:
                                    print 'bla 932!'
                    else:
                        M2[kk][r,2]=np.nan
                            
                M2[kk]=M2[kk][np.invert(np.isnan(M2[kk][:,2])),:]
                kk=kk+1
                
        M2[kk]=np.ndarray((numrow,3))
        M2[kk][:,1]=0
        V[kk]='PIM'
        U[kk]='g/m3'
        
        for r in range(numrow):
            try:
                M2[kk][r,0]=M2[TSM_index][r,0]
                M2[kk][r,2]=M2[TSM_index][r,2]*(1-0.01*M2[LOI_index][r,2])
            except:
                M2[kk][r,:]=np.nan
                
        M2[kk]=M2[kk][np.invert(np.isnan(M2[kk][:,2])),:]
              
        #print 'bla 1190: ',np.shape(M2),kk+1
        M2[kk+1]=np.ndarray((numrow,3))
        M2[kk+1][:,1]=0
        V[kk+1]='POM'
        U[kk+1]='g/m3'
        
        for r in range(numrow):
            try:
                M2[kk+1][r,0]=M2[TSM_index][r,0]
                M2[kk+1][r,2]=M2[TSM_index][r,2]*0.01*M2[LOI_index][r,2]
            except:
                M2[kk+1][r,:]=np.nan
                
        M2[kk+1]=M2[kk+1][np.invert(np.isnan(M2[kk+1][:,2])),:]
        
        
        # Mv,tvec_s=mapon_container_tvec(M2) #map all vars on along a container time vector
        Mv, tvec_s, zvec = mapon_container_tzvec(M2)  # map all vars along a container time/z vector
        tvec=[datetime.datetime(refdate.year,refdate.month,refdate.day,0,0,0)+datetime.timedelta(0,t_s) for t_s in tvec_s]
        #create_ncfile(fout,lon,lat,tvec,-1,Mv,V,V,U,dims,refdate=refdate,missing_value=-99,notify=True)
        lon2=np.array([lon[si]])
        lat2=np.array([lat[si]])
        create_ncfile(fout[si], lon2, lat2, Mv, V, V, U, dims, tvec, zvec, climatology=False,
                    refdate=refdate, missing_value=-99, notify=True)
        #print 'file written:'+fout
        addncatt(fout[si], {'station': station, 'bottom_depth': maxz})

def read_SchlHolst2(fin,fout,dims,station='-'):
    #better not modify this variable:
    coords_SH={
        '220017':[54.23,8.6],
        '220057':[54.71,8.23],
        '220006':[54.59,8.39],
        '220051':[54.12,8.86],
        '220052':[54.11,8.46],
        '220065':[54.00,8.67]
             }
    
    stations_ln={
        '220017':'Eider_west_Fliegenplate',
        '220057':'Hoernum_Vortrapptief',
        '220006':'Suedl_Amrum',
        '220051':'Buesum',
        '220052':'Suederpiep',
        '220065':'Norderelbe'
        }
    
    #comment out the unwanted stations, add additional stations 
    stations_SH=[
        '220017', #Eider west. Fliegenplate
        '220057', #Hoernum/Vortrapptief
        '220006', #Suedl. Amrum
        '220051', #Buesum
        '220052', #Suederpiep
        '220065'  # Norderelbe
        ]
    
    varkeys={
        'Ammonium-N':'NH4',
        'Chlorophyll-A':'chl',
        'Nitrat-N':'NO3',
        'Nitrit-N':'NO2',
        'Salzgehalt':'SALT',
        'Sichttiefe':'zS',
        'Silikat-Si':'Si',
        'o-Phosphat-P':'PO4'
        }
    
    conversion_factor={
        'Ammonium-N':7.14e1,
        'Chlorophyll-A':1,
        'Nitrat-N':7.14e1,
        'Nitrit-N':7.14e1,
        'Salzgehalt':1,
        'Sichttiefe':1,
        'Silikat-Si':3.56e1,
        'o-Phosphat-P':3.23e1
            }
    
    varunits={
        'Ammonium-N':'mmol N/m3',
        'Chlorophyll-A':'mg/m3',
        'Nitrat-N':'mmol N/m3',
        'Nitrit-N':'mmol N/m3',
        'Salzgehalt':'PSU',
        'Sichttiefe':'m',
        'Silikat-Si':'mmol Si/m3',
        'o-Phosphat-P':'mmol P/m3'
            }
    #read nut
    with open(fin, 'rb') as csvfile:
        RDraw=csvfile.readlines()
    headers=RDraw[0].split(';')
    del RDraw[0]
    numcolO=len(headers)
    numrow=len(RDraw)

    
    Mraw=[None]*len(stations_SH)
    Varraw=[None]*len(stations_SH)
    Uraw=[None]*len(stations_SH)
    Vraw=[None]*len(stations_SH)
    
    secs=np.ndarray((numrow,1)) #seconds
    for r in np.arange(numrow):
        try:
            years=int(RDraw[r].split(';')[2].split('.')[2])
            months=int(RDraw[r].split(';')[2].split('.')[1])
            days=int(RDraw[r].split(';')[2].split('.')[0])
        except:
            years=int(RDraw[r].split(';')[2].split(' ')[2])
            months=int(RDraw[r].split(';')[2].split(' ')[1])
            days=int(RDraw[r].split(';')[2].split(' ')[0])
            
        deltadate=datetime.date(years,months,days)- datetime.date(refdate.year,refdate.month,refdate.day)
        secs[r]=86400*deltadate.days +deltadate.seconds
        
    for sic, si in enumerate(stations_SH):
        print 'Station: '+stations_ln[si]
        lat=np.array([coords_SH[si][0]])
        lon=np.array([coords_SH[si][1]])
        print ' finding bottom depth..',
        maxz = get_botdepth(lon[0], lat[0], 'tree')
        print ' found: %sN, %sE: %.1f m' % (lat[0], lon[0], maxz)
        
        fout2=os.path.join(fout,stations_ln[si]+'.nc')
        
        stati=np.array([RDraw[r].split(';')[0]==si for r in range(numrow)])
        stati2=np.where(stati)[0]
        
        if np.sum(stati)>0:
            Mraw[sic]=np.ndarray((numrow,3))
            Varraw[sic]=[None]*numrow
            for r0,r1 in enumerate(np.where(stati)[0]):
                Mraw[sic][r0,0]=secs[r1]
                try:
                    Mraw[sic][r0,1]=RDraw[r1].split(';')[5]
                except:
                    Mraw[sic][r0,1]=np.nan
                
                val=RDraw[r1].split(';')[9]
                Mraw[sic][r0,2]=val
                Varraw[sic][r0]=RDraw[r1].split(';')[6]

            Mraw[sic]=Mraw[sic][np.invert(np.isnan(Mraw[sic][:,1])),:]
            Mraw[sic]=Mraw[sic][np.invert(np.isnan(Mraw[sic][:,2])),:]
            
            varar=np.unique(Varraw[sic])
            varar=filter(None,varar)
            
            M=[None]*len(varar)
            U=[None]*len(varar)
            V=[None]*len(varar)
            for varc,var in enumerate(varar):
                vari=np.array([Varraw[sic][r]==var for r in range(len(Varraw[sic]))])
                laux=len(np.where(vari)[0])
                
                M[varc]=np.ndarray((len(vari),3))
                
                V[varc]=varkeys[var]
                U[varc]=varunits[var]
                
                conv_fac2=conversion_factor[var]
                for r0,r1 in enumerate(np.where(vari)[0]):
                    
                    val=Mraw[sic][r1,:]
                    #if val[0]<0:
                        #print 'bla 1687: ',val[0]
                    val[2]=val[2]*conv_fac2
                    M[varc][r0,:]=val
                
                M[varc]=M[varc][0:laux,:]
                
            # Mv,tvec_s=mapon_container_tvec(M) #map all vars on along a container time vector
            #if sic==2:
            Mv, tvec_s, zvec = mapon_container_tzvec(M)  # map all vars along a container time/z vector
            tvec=[datetime.datetime(refdate.year,refdate.month,refdate.day,0,0,0)+datetime.timedelta(0,t_s) for t_s in tvec_s]
            #create_ncfile(fout,lon,lat,tvec,-1,Mv,V,V,U,dims,refdate=refdate,missing_value=-99,notify=True)
            create_ncfile(fout2, lon, lat, Mv, V, V, U, dims, tvec, zvec, climatology=False,
                        refdate=refdate, missing_value=-99, notify=True)
            #print 'file written:'+fout
            addncatt(fout2, {'station': stations_ln[si], 'bottom_depth': maxz})
        else:
            print 'station not found in file!\n'

def read_helgoland(dims, rootpath, fout, lat, lon, station='-'):

    #read_newhelgoland()
    combine_helgoland(dims,rootpath, fout,lat,lon,station=station)

def combine_helgoland(dims,rootpath, fout,lat,lon,MV=-99,station='-'):
    print(' finding bottom depth..',)
    # maxz=np.nan
    maxz = get_botdepth(lon[0], lat[0], 'tree')
    print(' found: %sN, %sE: %.1f m' % (lat[0], lon[0], maxz))

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
    for var,unit in vars.items():
        V.append(var)
        U.append(unit)
        m=np.empty((0,))
        sname = '?'
        print(var+':',)
        for fno, fin in enumerate(f):
            if var in varkeys[fno]:
                varn=varkeys[fno][var]
                print('f'+str(fno+1)+':'+varn,)
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
                print('f'+str(fno+1)+': NA',)
                v=np.nan*np.ones((sum(tinds[fno]),))
            m = np.hstack((m, v))
        print('.')
        m[np.isnan(m)]=MV
        M.append(m)
        L.append(sname)

    create_ncfile(fout, lon, lat, M, V, L, U, dims, datevec, -1, climatology=False,
                  refdate=refdate, missing_value=MV, notify=True)
    addncatt(fout, {'station': station, 'bottom_depth': maxz})


def addncatt (fout,atts):
    nc=netCDF4.Dataset(fout,'a')
    for attkey,attval in atts.items():
        if attkey=='station':
            nc.station=attval
        elif attkey=='bottom_depth':
            nc.bottom_depth=attval
    nc.sync()
    nc.close()

def readstats(rootpath,stations,read):

    coords={'Bork': [53.48,6.92],
             'Helgoland': [54.18,7.90],  #54 11'N 7.54'E
             'JaBu':[53.51,8.15],
             'Norderney':[53.7,7.17],
             'Nney':[53.7,7.17],
             'SAmrum':[54.59,8.39],
             'Norderelbe':[54.0,8.67],
             'Sylt': [55.0, 8.45],  # 55.0'N 8.27'E
             'Suederpiep':[54.10, 8.45],
             'Wesermuendung':[53.67, 53.63, 53.78, 8.38, 8.31, 8.46],
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
        Germanpathin='German'
        Germanpathout='allNC'
        Dutchpathin='Dutch'
        Dutchpathout='allNC'

        print('Reading:'+station)
        if station=='Bork':
            dims={'t':'time','z':'depth','x':'lon','y':'lat'}
            fin=os.path.join(rootpath,Germanpathin,'Bork_W_1_1994-2019.CSV')
            fout=os.path.join(rootpath,Germanpathout, 'Bork.nc')
            if read: read_Bork(rootpath,fin, fout,dims,lat=np.array([coords[station][0]]),lon=np.array([coords[station][1]]),station='Bork')
        elif station=='Helgoland':
            #fin=
            fout=os.path.join(rootpath, 'Helgoland', 'Helgoland.nc')
            if read: read_helgoland(dims, rootpath, fout, lat=np.array([coords[station][0]]),lon=np.array([coords[station][1]]),station='Helgoland')
        elif station=='JaBu':
            dims={'t':'time','z':'depth','x':'lon','y':'lat'}
            fin=os.path.join(rootpath,Germanpathin,'JaBu_W_1_chem_f_Projekt_2009-19.CSV')
            fout=os.path.join(rootpath, Germanpathout, 'JaBu.nc')
            if read: read_JaBu(rootpath,fin, fout,dims,lat=np.array([coords[station][0]]),lon=np.array([coords[station][1]]),station='JaBu')
        elif station=='Norderney':
            fin=os.path.join(rootpath,'Norderney','Nney_W_2_1999-2014_chl.csv')
            fout=os.path.join(rootpath, 'Norderney', 'Norderney.nc')
            if read: read_norderney(rootpath,fin, fout,dims,lat=np.array([coords[station][0]]),lon=np.array([coords[station][1]]),station='Norderney')
        elif station=='Nney':
            dims={'t':'time','z':'depth','x':'lon','y':'lat'}
            fin=os.path.join(rootpath,Germanpathin,'Nney_W_2_chem_f_Projekt_1999-18.CSV')
            fin2=os.path.join(rootpath,Germanpathin,'Chlorophyll1_Nney_W_2_1999-19.CSV')
            fin3=os.path.join(rootpath,Germanpathin,'Nney_W_2_chem_f_Projekt_1999-18_tab2.CSV')
            fout=os.path.join(rootpath, Germanpathout, 'Nney.nc')
            if read: read_Nney(rootpath,fin, fin2, fin3, fout,dims,lat=np.array([coords[station][0]]),lon=np.array([coords[station][1]]),station='Nney')
        elif station== 'SAmrum':
            dims={'t':'time','x':'lon','y':'lat'}
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
        elif station=='SchlHolst2':
            dims={'t':'time','z':'depth','x':'lon','y':'lat'}
            fin=os.path.join(rootpath,Germanpathin,'JMP-EUNOSAT_Data_North_Sea_Schleswig-Holstein_2015-2017.CSV')
            fout=os.path.join(rootpath, Germanpathout)
            if read: read_SchlHolst2(fin,fout,dims,station='SchlHolst2')
        elif station=='Sylt':
            fout=os.path.join(rootpath,'justusMA','Sylt_1973-2013.nc')
        elif station in ['JMA_Norderelbe','JMA_Suederpiep','JMA_Westerhever']:
            stname=station.split('JMA_')[1]
            fin = os.path.join(rootpath, 'justusMA', stname+'.csv')
            fout = os.path.join(rootpath, 'justusMA', stname+'.nc')
            if read: read_justusMA(fin, fout, dims, lat=np.array([coords[stname][0]]), lon=np.array([coords[stname][1]]),station=stname)
        elif station=='Wesermuendung':
            dims={'t':'time','z':'depth','x':'lon','y':'lat'}
            fout=[None]*3
            fin=os.path.join(rootpath,Germanpathin,'Wesermuendung2008-2019_fuer_Projekt.CSV')
            fout[0]=os.path.join(rootpath, Germanpathout, 'WeMu_W1.nc')
            fout[1]=os.path.join(rootpath, Germanpathout, 'WeMu_W2.nc')
            fout[2]=os.path.join(rootpath, Germanpathout, 'Wu_KU-W1.nc')
            if read: read_Wesermuendung(rootpath,fin, fout, dims,lat=np.array([coords[station][i] for i in range(0,3)]),lon=np.array([coords[station][i] for i in range(3,6)]),station='Wesermuendung')
        elif station== 'Norderelbe':
            fin=os.path.join(rootpath,'SchlHolst','DIN-DIP-Chl-220006-220065-2000-2014_Norderelbe.csv')
            fout=os.path.join(rootpath, 'SchlHolst', 'Norderelbe.nc')
            if read: read_SchlHolst(fin, fout,dims,lat=np.array([coords[station][0]]),lon=np.array([coords[station][1]]),station='Norderelbe')
        elif station in ['T36','T26','T41','T8','T2','T22','T5','T12','T11']:
            dims={'t':'time','z':'depth','x':'lon','y':'lat'}
            fin=os.path.join(rootpath,'BSH',station+'.csv')
            fout=os.path.join(rootpath,'BSH',station+'.nc')
            if read: read_bsh(fin,fout,dims,lat=np.array([coords[station][0]]),lon=np.array([coords[station][1]]),station=station)
        elif station[0:6] in ['NOORDW','TERSLG','ROTTMP','BOCHTV','BOOMKD','DANTZG','DOOVBW','GROOTG','HUIBGO','MARSDN','ZOUTKP','ZUIDOL']:
            dims = {'t': 'time', 'z': 'depth', 'x': 'lon', 'y': 'lat'}
            #fin=os.path.join(rootpath,'Dutch','opendap_new','raw',station)
            fin=os.path.join(rootpath,Dutchpathin,'Monitoring_data_crosstable.CSV')
            fin2=os.path.join(rootpath,Dutchpathin,'Monitoring_station_description.CSV')
            fout = os.path.join(rootpath,Dutchpathout, station+'.nc')
            #if read: read_opendapnc(fin,fout,dims,station=station)
            if read: read_dutchcsv(rootpath,fin,fin2,fout,dims,station=station)

        fstats[station]=fout

    return fstats

if __name__=='__main__':
    main()
