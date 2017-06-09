# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 16:06:57 2017

@author: ivan.kuznetsov@gmail.com, kerimoglu.o@gmail.com
"""
import os
import pickle
import numpy as np
from netCDF4 import num2date
from netCDF4 import Dataset as open_ncfile

def readobs(paths,statsets,stations,timeint,depthints):
    print ('Reading observations:')

    picklecode = '_%s_%s_%s-%s' % ('-'.join(statsets), '-'.join(depthints.keys()), timeint[0].year, timeint[1].year)
    pickledobsfile = os.path.join(paths['rootpath'], 'obs' + picklecode + '.pickle')

    #check if the pickledobs exist
    if os.path.exists(pickledobsfile):
       print('Opening pickled obs file')
       (obs,) = np.load(pickledobsfile)
       # filter requested stations?
       return obs

    # if pickledfile does not exist:
    obs={}
    for statset in statsets:
        print('Filling obs. dataset:%s'%statset)
        obs=looppath_fill_stationdata_obs(obs,paths,statset,stations,timeint,depthints)
    
    #pickle the obs file
    f=open(pickledobsfile,'wb')
    pickle.dump((obs,),f) #,protocol=-1
    print('Pickled obs file for later use:' + pickledobsfile)
    f.close()

    return obs

def looppath_fill_stationdata_obs(obs,paths,statset,stations,timeint,depthints):

    #TODO: pickle/unpickle each station set?

    obspath = paths[statset]

    #if stations to include are not specified,
    if len(stations) == 0:
        # get lists of available data
        files = [f for f in os.listdir(obspath) if f.endswith('.nc')]
        stations= files[:]

    for station in stations:
        print('  '+station)
        sdata = fill_stationdata_obs(os.path.join(obspath,station),statset,timeint,depthints)
        obs[station] = sdata

    return obs

def fill_stationdata_obs(file,statset,timeint,depthints):

    tempfound = False
    saltfound = False
    sshfound = False

    if statset in ['marnet']:
        vlib={'t':'TIME','x':'LONGITUDE','y':'LATITUDE','z':'DEPH','temp':'TEMP','salt':'PSAL', 'ssh':'?'}
    elif statset in ['emodnet']:
        vlib = {'t':'TIME','x': 'LONGITUDE', 'y': 'LATITUDE', 'z': 'DEPH', 'temp': 'TEMP', 'salt': 'PSAL','ssh':'SLEV'}
    elif statset in ['cosyna']:
        vlib = {'t': 'TIME', 'x': 'LONGITUDE', 'y': 'LATITUDE', 'z': 'DEPTH', 'temp': 'TEMP', 'salt': 'PSAL','ssh': 'SLEV'}
    ncf = open_ncfile(file,'r')

    # reading metadata of station
    lon = ncf.variables[vlib['x']][:][0]
    lat = ncf.variables[vlib['y']][:][0]
    depth = ncf.variables[vlib['z']][:][0]
    time_num = ncf.variables[vlib['t']][:]
    time = num2date(time_num, ncf.variables[vlib['t']].getncattr('units'))

    # find the max_depth, if necessary
    if 'bottom' in depthints.keys():
        maxz=get_maxdepth_obs(lon,lat)
        #update the bottom depthint
        depthints['bottom']=[maxz-depthints['bottom'][0], maxz-depthints['bottom'][1]]

    # check if dates are in
    tind=np.where((time>=timeint[0]) * (time<=timeint[1]))[0]

    # check if depths are in
    depthintmin=np.min([dint[0] for dint in depthints.values()]) #find the minimum lower lim of depthints
    depthintmax = np.max([dint[1] for dint in depthints.values()]) #find the maximum upper lim of depthints
    zind=np.where((depth>=depthintmin) * (depth<=depthintmax))[0]

    # check if variables are in, decide if anything relevant found
    if len(tind)>0 and len(zind)>0 and vlib['temp'] in ncf.variables: tempfound = True
    if len(tind)>0 and len(zind)>0 and vlib['salt'] in ncf.variables: saltfound = True
    if len(tind)>0 and len(zind)>0 and vlib['ssh'] in ncf.variables: sshfound = True

    # a = ncf.variables['TIME'].getncattr('units')
    # time_emodnet.append([num2date(ncf.variables['TIME'][:][0],a),num2date(ncf.variables['TIME'][:][-1],a)])
    # ncf.close()

    #          if var in var_emodnet[i][:]:
    #              if (time_emodnet[i][0] < t2 or time_emodnet[i][1] > t2):
    #                  # get data from file (first match)
    #                  ncf = open_ncfile(path2data_emodnet+item,'r')
    #                   obs_z    = ncf.variables['DEPH'][:][0]
    #                   time = ncf.variables['TIME'][:]
    #                   a = ncf.variables['TIME'].getncattr('units')
    #                   time = num2date(time,a)
    #                   if np.any(time > t1) and np.any(time < t2):
    #                       time_ind1 = np.where(time > t1)[0][0]
    #                       time_ind2 = np.where(time < t2)[0][-1]
    #                   else:
    #                       continue
    #                   data = ncf.variables[v2v[var]][:][time_ind1:time_ind2,:]
    #                   ncf.close()
    #                   ind_time = time[time_ind1:time_ind2]
    #                   for iz,z in enumerate(obs_z):
    #                       #check if data array is masked, if not make new masked array
    #                       if isinstance(data[:,iz],np.ma.MaskedArray):
    #                           dd = data[:,iz]
    #                       else:
    #                           dd = np.ma.array(data[:,iz])
    #                       if dd.count():
    #                           dset={'fname':item,'time':ind_time,'data':dd,'z':z}
    #                           dset_list.append(dset)

    # fill in the data:
    tempdata = {}; saltdata = {}; sshdata = {}

    #handle temp
    if tempfound:
        tempdata['presence']=True
        for layername, depthint  in depthints.items():
            zind=np.where((depth>=depthint[1]) * (depth<=depthint[1]))[0]
            if len(zind)>0:
                data = ncf.variables[vlib['temp']][tind,zind]
                tempdata[layername] = {'time': time[tind], 'value': data, 'depth_interval': depthint}
            else:
                tempdata[layername] = {'time': np.array([]), 'value': np.array([]), 'depth_interval': depthint}


    else:
        tempdata['presence']=False

    #handle salt
    if saltfound:
        saltdata['presence']=True
        for layername, depthint  in depthints.items():
            zind=np.where((depth>=depthint[0]) * (depth<=depthint[1]))[0]
            if len(zind)>0:
                data= ncf.variables[vlib['salt']][tind,zind]
                saltdata[layername] = {'time': time[tind], 'value': data, 'depth_interval': depthint}
            else:
                tempdata[layername] = {'time': np.array([]), 'value': np.array([]), 'depth_interval': depthint}
    else:
        saltdata['presence']=False

    #handle ssh
    if sshfound:
        sshdata['presence']=True
        data = ncf.variables[vlib['ssh']][tind]
        sshdata['z0']={'time':time, 'value':data, 'depth_interval':[0,0]}
    else:
        sshdata['presence']=False

    #put all data in
    sdata ={'longname':'descriptive name of the station',
            'lon':lon,
            'lat':lat,
            'max_depth':0.0,
            'temp' : tempdata,
            'salt' : saltdata,
            'ssh': sshdata
            }

    ncf.close()
            
    return sdata

def get_maxdepth_obs(lon,lat):
    maxz=100 #TODO:find the depth using lon,lat
    return maxz