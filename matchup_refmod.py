"""
Created on 23 Nov 2016
@author: onur.kerimoglu@hzg.de
"""

import os,sys
import pickle
import numpy as np
import netCDF4
import netcdftime
import datetime
from scipy.spatial import cKDTree

from getm_funcs import get_getm_bathymetry_cropped

cordict={'mgN m-3':'$\mu$MN'}
corfacts={'mgN m-3':1./14.0067}

def get_matchups(files,ftypes,fnames,varnames,yrs):

    #for each var, save/load an individual match library, and combine them in one set
    matchupset = {};
    unitset = {}
    for v in varnames:
        # construct a pickle name from varname
        picklefname=files[1].split('.nc')[0] + '_matchup_' + fnames[0] + '_' + v + '_' + str(yrs[0]) + '-'+str(yrs[1]) +'.pickle'

        #check if the pickled data exists:
        if os.path.isfile(picklefname):
            matchupsetV,unitsetV = np.load(picklefname)
            print ('un-pickled:'+picklefname)
        else:
            #if it does not exist:
            print ('Collecting data:')
            matchupsetR = {}
            unitsetR = {}
            print ('  Ref:'+files[0])

            #retrieve the model as a data object, for the requested year interval
            timeint=[datetime.date(yrs[0], 1, 1), datetime.date(yrs[1], 12, 31)]
            model = data(ftypes[1], files[1], [v], timeint)

            # refine the timeint based on the available data
            timeint = [model.dates[0], model.dates[-1]]

            # retrieve observations as a data object, trim the data outside the timeint
            ref = data(ftypes[0], files[0], [v], timeint)

            #find the matchups
            matchupsetVR, unitsetVR=find_matchups(ref,model,[v])

            #convert the non-matching units
            matchupsetV,unitsetV=convert_units(matchupsetVR,unitsetVR)

            #save the data as pickle
            f=open(picklefname,'wb')
            pickle.dump([matchupsetV,unitsetV],f) #,protocol=-1
            f.close()
            print ('pickled:'+picklefname)

        #collect in one set
        matchupset[v]=matchupsetV[v]
        unitset[v]=unitsetV[v]
    return (matchupset,unitset)

def find_matchups(ref,model,varnames):
    print ('\nFinding match-ups:')
    #create a tree for lateral interpolation:
    #tree = cKDTree(zip(model.lons.flat, model.lats.flat))
    tree = cKDTree(list(zip(model.lons.flat, model.lats.flat)))#p3: must be placed in list
    # for each variable:
    matchupsetR={};unitsetR={}
    for v in varnames:
        print ('  '+v+':'),

        #if list(model.depths)!=[-1]:
        # raise(ValueError('Unable to handle variable with non-squeezable (len>1) depth dimension'))
        # calculate a space-time average depth vector, to have an idea about the approximate depth of each layer
        # this is going to be used for matching the depths in ref data
        # Mdepths=np.mean(np.mean(np.mean(model.depths[:],0),2),3)

        #do the time matching, if necessary
        if (len(ref.dates)==1) and (len(model.dates)==1):
            #special case: only one time frame in ref data.
            #assume that the times of the ref and model match
            if len(model.values[v].shape)==3:
                MT=model.values[v][-1,:,:] #i.e., at the surface
            elif len(model.values[v].shape)==2:
                MT=model.values[v][:,:]
            else:
                raise(ValueError('cannot handle model with shape'+str(model.values[v].shape)))
        else:
            #pick the time-matching frames:
            sti=np.nan*np.ones(len(ref.dates[v]))
            for r, date in enumerate(ref.dates[v]):
                m=np.where(model.dates==date)[0]
                if len(m)>0:
                    sti[r]=m[0]
            #remove the no-match cases
            nomatchi=np.isnan(sti)
            if any(nomatchi):
                remi=np.where(nomatchi)[0]
                sti=np.delete(sti,remi)
                ref.dates[v]=np.delete(ref.dates[v],remi)
                ref.lons[v]=np.delete(ref.lons[v],remi)
                ref.lats[v] = np.delete(ref.lats[v], remi)
                ref.values[v] = np.delete(ref.values[v],remi)
                ref.depths[v] = np.delete(ref.depths[v], remi)
                ref.maxdepths[v] = np.delete(ref.maxdepths[v], remi)
            #isolate the match
            if len(model.values[v].shape) == 3:
                MT=model.values[v][sti.astype(int),:,:]
            elif len(model.values[v].shape) == 4:
                #ref.depths[v] = np.delete(ref.depths[v], remi)
                MTZ = model.values[v][sti.astype(int), :, :, :]
                if list(model.depths)==[-1]:
                    if model.values[v].shape[1]==2: #assume that the 2 layers are the surface and bottom
                        MdepthsS=np.tile(1.0, (MTZ.shape[0],1,MTZ.shape[2],MTZ.shape[3]))
                        MdepthsB=np.tile(model.maxdepths-2.0,(len(sti),1,1,1))
                        Mdepths=np.concatenate((MdepthsS,MdepthsB),axis=1)
                        print ('.')
                else:
                    Mdepths=model.depths[sti.astype(int),:,:,:]

        #do the lateral interpolation, if necessary
        if (len(ref.dates[v])==1) and (ref.lons[v]==[-1]) and (MT.shape==ref.values[v].shape):
            reportskipstats = False
            #special case: no ref.lons was given, and the shapes of model and ref match
            #assume that they are on identical coords
            Vmodel=MT.flatten()
            if len(ref.values[v].shape) == 3:
                Vref=ref.values[v][-1,:,:].flatten() #i.e., at the surface
            elif len(ref.values[v].shape) == 2:
                Vref = ref.values[v][:,:].flatten()
            else:
                raise(ValueError('cannot handle ref with shape'+str(ref.values[v].shape)))
            dates=np.array([model.dates[0] for i in range(len(Vref))])
            lats=model.lats.flatten()
            lons=model.lons.flatten()
            depths = [-1]
            maxdepths=[-1] #todo:can be included?
        else:
            reportskipstats=True
            Vref = ref.values[v]
            dates = ref.dates[v]
            lats = ref.lats[v]
            lons = ref.lons[v]


            depths=ref.depths[v]
            maxdepths = ref.maxdepths[v]

            # interpolate the model results

            #vectorized method does not work:
            # if list(model.depths) == [-1]:
            #     # lateral indices to be used later for interpolation
            #     d, inds = tree.query(zip(ref.lons[v], ref.lats[v]), k=4) #,distance_upper_bound=0.15)
            #     w = 1.0 / d ** 2  # weights
            #     Vmodel=np.sum(w * MT.flatten()[inds], axis=1) / np.sum(w, axis=1)
            # else:
            #     raise(ValueError('No quick method when depth dimension is present'))

            #loop over each observation instance
            Vmodel=np.nan*Vref #allocate
            latdistskip=0
            vertdistskip=0
            maskskip=0
            for i,date in enumerate(dates):
                skip=False #by default, don't skip
                #find spatial weights
                d, inds = tree.query([ref.lons[v][i], ref.lats[v][i]], k = 4) #,distance_upper_bound=0.15)
                w = 1.0 / d ** 2  # weights
                Mlon=np.sum(w * model.lons.flatten()[inds]) / np.sum(w)
                Mlat = np.sum(w * model.lats.flatten()[inds]) / np.sum(w)
                #if the nearest cell in the model is too far from the ref, discard  it
                dist=distance((Mlat,Mlon),(ref.lats[v][i],ref.lons[v][i]))
                if dist>2.: #(km)
                    latdistskip=latdistskip+1
                    skip=True
                if list(model.depths) == [-1] and len(model.values[v].shape) == 3:
                    M=MT[i,:,:]
                else: #choose the correct layer for model
                    Zref=ref.depths[v][i]
                    #calculate an average depth grid
                    Zavg=np.nan*np.zeros(Mdepths.shape[1])
                    for mzii in range(Mdepths.shape[1]):
                        Zti=Mdepths[i,mzii,:,:]
                        Zti_latlon=Zti.flatten()[inds]
                        if np.ma.is_masked(Zti_latlon) or np.any(np.isnan(Zti_latlon)): #if data coords comeout masked, skip
                            maskskip = maskskip + 1
                            skip=True
                        Zavg[mzii] = np.sum(w * Zti_latlon) / np.sum(w)
                    #find the layer where the modeled depth is closest to the ref depth
                    Zabsdif=np.abs(Zref-Zavg)
                    if Zabsdif.min()>5: #if min.abs.dif.>5m, skip
                        vertdistskip=vertdistskip+1
                        skip=True
                    mzi=np.where(Zabsdif.min()==Zabsdif)[0]
                    M=MTZ[i,mzi,:,:]
                if not skip: #if skip, the instance will be left as a nan, which will be removed later
                    Vmodel[i] = np.sum(w * M.flatten()[inds]) / np.sum(w)
        #remove any potential nan instance
        #first convert the masked elements to nan
        if np.ma.is_masked(Vmodel):
            maski=np.where(Vmodel.mask)[0]
            Vmodel.mask=False
            Vmodel[maski]=np.nan
        if np.ma.is_masked(Vref):
            maski = np.where(Vref.mask)[0]
            Vref.mask = False
            Vref[maski] = np.nan
        vali=np.where(np.isfinite(Vmodel) * np.isfinite(Vref))[0]
        if reportskipstats:
            print (str(len(vali)) + ' data points out of %s (skipped %s, %s and %s instances due to latdist,vertdist,mask)'%(len(dates),latdistskip,vertdistskip,maskskip))

        if list(depths) != [-1]:
            depths=depths[vali]
        if list(maxdepths) != [-1]:
            maxdepths=maxdepths[vali]

        #store
        matchupsetR[v]={'dates': dates[vali],'depths': depths, 'maxdepths':maxdepths, 'lats':lats[vali] , 'lons': lons[vali], 'ref': Vref[vali], 'model':Vmodel[vali]}
        unitsetR[v] = {'ref': ref.units[v], 'model': model.units[v]}

    return (matchupsetR,unitsetR)

def distance(origin, destination):
    import math
    lat1, lon1 = origin
    lat2, lon2 = destination
    radius = 6371 # km
    dlat = math.radians(lat2-lat1)
    dlon = math.radians(lon2-lon1)
    a = math.sin(dlat/2) * math.sin(dlat/2) + math.cos(math.radians(lat1)) \
    * math.cos(math.radians(lat2)) * math.sin(dlon/2) * math.sin(dlon/2)#get the topography
    #topo=get_getm_bathymetry()
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
    d = radius * c
    return d

def filter_matchups_region(Vset0,r):
    Rdict = {'NW': [[53.5,56.0],[0,5.0]], 'NE': [[54.5,56.0],[5.0,10.0]],
             'SW': [[51.0,53.5],[0,5.0]], 'SE': [[51.0,54.5],[5.0,10.0]],
             'W': [[51.0,56.0],[0,5.0]], 'E': [[51.0,56.0],[5.0,10.0]],'all':[]}

    if r=='all':
        suf=''
        return (Vset0,suf)
    else:
        suf = '_' + r

        Vset = {}
        for v in Vset0.keys():
            latint = Rdict[r][0]
            lonint = Rdict[r][1]
            ri = (Vset0[v]['lats'] >= latint[0]) * (Vset0[v]['lats'] <= latint[1]) * \
                 (Vset0[v]['lons'] >= lonint[0]) * (Vset0[v]['lons'] <= lonint[1])

            Vset[v] = {'dates': Vset0[v]['dates'][ri], 'lats': Vset0[v]['lats'][ri], 'lons': Vset0[v]['lons'][ri],
                    'depths': Vset0[v]['depths'][ri], 'maxdepths':Vset0[v]['maxdepths'][ri],
                   'ref': Vset0[v]['ref'][ri], 'model': Vset0[v]['model'][ri]}

        return (Vset, suf)

def filter_matchups_vertloc(Vset0,vertloc):
    suf = '_V' + vertloc

    if vertloc=='all':
        return (Vset0,suf)
    else:
        Vset0keys=list(Vset0.keys())
        if list(Vset0[Vset0keys[0]]['depths'])==[-1]:
            raise(Warning('depths not available, cannot filter the vertloc'))
            return (Vset0,suf)
        print ('extracting %s:'%vertloc)
        Vset={}
        for v in Vset0keys:
            Vv=Vset0[v]

            #find the valid indices
            if vertloc=='surf':
                vali=np.where(Vv['depths']<= 5.0)[0]
            elif vertloc == 'bot':
                Zdif = Vv['maxdepths'] - Vv['depths']
                vali = np.where(Zdif <= 5.0)[0]
            else:
                raise(ValueError('vertloc %s unidentified for filtering'%vertloc))

            #extract the valid indices
            Vset[v]={}
            for par in Vv.keys():
                Vset[v][par] = Vv[par][vali]
            print ('  %s: %s values out of %s'%(v,len(vali),len(Vv[par][:])))

        return (Vset,suf)


def convert_units(matchupset,unitset):
    unitsetC={}
    print ('\nCross-checking units:')
    for v in unitset.keys():
        print ('  '+v+':'),
        if unitset[v]['ref']==unitset[v]['model']:
            unitsetC[v]=unitset[v]['ref']
            print ('OK')
        else:
            corrected=False
            if unitset[v]['ref'] in cordict.keys():
                unitsetC[v]=cordict[unitset[v]['ref']]
                matchupset[v]['ref']=matchupset[v]['ref']*corfacts[unitset[v]['ref']]
                corrected=True
                print (' ref:'+ unitset[v]['ref']+'->'+unitsetC[v]),
            if unitset[v]['model'] in cordict.keys():
                unitsetC[v]=cordict[unitset[v]['model']]
                matchupset[v]['model']=matchupset[v]['model']*corfacts[unitset[v]['model']]
                corrected=True
                print (' model:' + unitset[v]['model'] + '->' + unitsetC[v]),
            print ('') #to end the line

            #throw error if no correction could be made
            if not corrected:
                msg='Unable to correct the unit mismatch for %s. ref:%s, model:%s'%(v,unitset[v]['ref'],unitset[v]['model'])
                raise (ValueError(msg))

    return(matchupset,unitsetC)

class data(object):
    def __init__(self,type,file,varnames,timeint=np.nan):
        self.type=type
        self.file=file
        self.varnames=varnames
        self.Xdict={'T':'T','S':'S','Chl':'Chl','NO3':'N','NH4':'N','DIP':'P'}

        if 'GETM' in type:
            type_phys=type.split('-')[0]
            type_fabm = type.split('-')[1]
            vardict={}
            if type_phys == 'GETM':
                vlibphys = {'T': 'temp', 'S': 'salt'}
            elif type_phys == 'GETM.M':
                vlibphys={'T':'tempmean','S':'saltmean'}
            vardict.update(vlibphys)

            if type_fabm == 'MAECS':
                vlibfabm = {'Chl':'hzg_maecs_chl', 'NO3': 'hzg_maecs_nutN', 'DIP': 'hzg_maecs_nutP', 'Kd':'hzg_maecs_att'}

            elif type_fabm == 'GPMEH.PZ':
                vlibfabm = {'DOs': 'EH_abioP_O2_percSat', 'NO3': 'EH_abioP_DINO3', 'NH4': 'EH_abioP_DINH4', 'DIP': 'EH_abioP_DIP',
                            'Chl': 'GPM_phy_Chl'}
            elif type_fabm[:3] == 'GPM':
                vlibfabm = {'DOs': 'EH_abioP_O2_percSat', 'NO3': 'EH_abioP_DINO3', 'NH4': 'EH_abioP_DINH4', 'DIP': 'EH_abioP_DIP',
                            'Chl': 'total_chlorophyll_calculator_result'}
            vardict.update(vlibfabm)
            #unitdict = {'[d': '$^\circ$C', '[degC]': '$^\circ$C', '[psu]': 'g/kg', '[ug/l]': 'mgX m-3','[umol/l]': 'mmolX m-3'}
            unitdict = {'T':u'\N{DEGREE SIGN}C','S':'g/kg','Chl': 'mgChl/m$^3$', 'NO3': '$\mu$MN', 'NH4': '$\mu$MN', 'DIP': '$\mu$MP'}
            self.dates,self.depths,self.lons,self.lats,self.maxdepths,self.values,self.units = self.get_getmdata(vardict,unitdict,timeint)
        elif 'ICES' in type:
            vardict = {'T': 'TEMP', 'S': 'PSAL', 'Chl': 'CPHL', 'NO3': 'NTRA', 'NH4': 'AMON', 'DIP': 'PHOS'}
            unitdict = {'[degC]': u'\N{DEGREE SIGN}C', '[psu]\r\n': 'g/kg', '[ug/l]': 'mgX/m$^3$',
                        '[umol/l]': '$\mu$MX'}
            self.dates, self.depths, self.lons, self.lats,self.maxdepths,self.values, self.units = self.get_ices_data(timeint,vardict,unitdict)
        elif 'ESA-CCI' in type:
            vardict = {'Chl': 'chl', 'Kd': 'kd490'}
            unitdict = {'milligram m-3': 'mgX/m$^3$', 'm-1': 'm-1'}
            self.dates, self.depths, self.lons, self.lats, self.values, self.units = self.get_esacci_data(vardict,unitdict)
        else:
            raise(ValueError('unable to handle '+file+':unknown data type ('+type+')'))

    def get_getmdata(self,vardict,unitdict,timeint=np.nan):

        nc = netCDF4.Dataset(self.file)
        #dates
        tv=nc.variables['time']
        utime = netcdftime.utime(tv.units)
        datetimes = utime.num2date(tv[:])
        datesF=np.array([datetimes[r].date() for r in range(len(datetimes))])

        ti = (datesF >= timeint[0]) & (datesF <= timeint[1])
        dates = datesF[ti]

        #check if a non-squeezable (len>1) depth dimension exists:
        if 'level' not in nc.variables.keys():
            depthdim=False
        elif len(nc.variables['level'].shape)==0:
            depthdim = False
        elif nc.variables['level'].shape[0]==1:
            depthdim = False
        else:
            depthdim = True

        if not depthdim:
            depths=[-1]
        else:
            if ('depth' not in nc.variables.keys()) or ('bathymetry' not in nc.variables.keys()):
                raise(ValueError('data seems to have a depth dimension, but no depth or bathymetry provided'))
            #h=nc.variables['h'][:,1:,:,:]
            #bath=nc.variables['bathymetry'][:]
            #depths=np.nan*h[:]
            #for lidx in range(len(h)):
            #    if lidx == 0:
            #        depths[:,lidx,:,:] = -bath + 0.5 * h[:, lidx,:,:]
            #    else:  # lidx > 0
            #        depths[:, lidx, :, :] = depths[:,lidx-1,:,:] + 0.5*depths[:, lidx-1,:,:] + 0.5*depths[:, lidx,:,:]
            depths=-1*nc.variables['depth'][:]

        values={};  units={}
        for v in self.varnames:
            if len(nc.variables[vardict[v]][:].shape)==4:
                values[v]= nc.variables[vardict[v]][ti,:,:,:].squeeze() #squeeze is to get rid of len:1 dimension (eg., depth)
            else:
                values[v] = nc.variables[vardict[v]][ti, :, :]
            try:
                unitR=nc.variables[vardict[v]].units
            except:
                unitR='-'

            if v in unitdict.keys():
                if 'X' in unitdict[v]:
                    units[v] = unitdict[v].replace('X', self.Xdict[v])
                else:
                    units[v] = unitdict[v]
            else:
                units[v]=unitR

        ncvk=nc.variables.keys()
        if 'bathymetry' in ncvk and 'lon' in ncvk and 'lat' in ncvk:
            lons = nc.variables['lon'][:]
            lats = nc.variables['lat'][:]
            H = nc.variables['bathymetry'][:]
        else:
            topo = get_getm_bathymetry_cropped()
            lats = topo['lats']; lons = topo['lons'];H = topo['H']
        return (dates,depths,lons,lats,H,values,units)
        nc.close()

    def get_esacci_data(self,vardict,unitdict):

        nc=netCDF4.Dataset(self.file)

        # collect different variables in different structures
        dates = {};depths = {};lons = {}; lats = {};values = {}; units = {}
        for varno, v in enumerate(self.varnames):
            print ('  '+v),

            # dates
            tv = nc.variables['time']
            utime = netcdftime.utime(tv.units)
            datetimes = utime.num2date(tv[:])
            dates[v] = np.array([datetimes[r].date() for r in range(len(datetimes))])

            depths[v]=[-1]
            lons[v]=[-1]
            lats[v]=[-1]

            values[v]=nc.variables[vardict[v]][:].squeeze()
            unitR=nc.variables[vardict[v]].units
            if unitR in unitdict.keys():
                if 'X' in unitdict[unitR]:
                    units[v] = unitdict[unitR].replace('X', self.Xdict[v])
                else:
                    units[v] = unitdict[unitR]
            else:
                units[v]=unitR

        nc.close()
        return (dates, depths, lons, lats, values, units)

    # to solve compatibility problems  between py2 and py3
    def load_pickle(self,pickle_file):
        try:
            with open(pickle_file, 'rb') as f:
                pickle_data = pickle.load(f)
        except UnicodeDecodeError as e:
            with open(pickle_file, 'rb') as f:
                pickle_data = pickle.load(f, encoding='latin1')
        except Exception as e:
            print('Unable to load data ', pickle_file, ':', e)
            raise
        return pickle_data

    def get_ices_data(self,timeint,vardict,unitdict):
        # with open(self.file) as f:
        # dataT, varnsT, unitsT = pickle.load(f) #may lead to compatability problem with loading pickle saved with py2
        dataT, varnsT, unitsT = self.load_pickle(self.file)

        rd = datetime.date(2000, 1, 1)
        zmaxint=[0.,60.]
        zint = [0., 60.]
        lonint=[0.,9.]
        latint=[51.,55.5]

        # collect different variables in different structures
        dates={};depths={};lons={};lats={};maxdepths={};values={}; units={}
        for varno, v in enumerate(self.varnames):
            print ('  '+v),
            # find the corresponding variable
            vari = varnsT.index(vardict[v])
            if vari==[]:
                raise(ValueError('variable '+vardict[v]+' does not exist in '+self.file))
            dataV = dataT[vari]

            # if id in ['ices']:
            data= self.reduce_ices(dataV, rd, timeint, zint, zmaxint,lonint,latint)
            #data[varno]=[tcol,zcol,loncol,latcol,zmaxcol,varcol]

            dates[v]= np.array([rd + datetime.timedelta(0, data[r, 0]) for r in range(len(data[:, 0]))])  # reconstruct dates
            depths[v]=data[:,1]
            lons[v] = data[:, 2]
            lats[v] = data[:, 3]
            maxdepths[v] = data[:, 4]
            values[v] = data[:, 5]
            #store the unit, but in standard form if registered
            if unitsT[vari] in unitdict.keys():
                unitR = unitsT[vari]
                units[v] = unitdict[unitR].replace('X', self.Xdict[v])
            else:
                units[v]=unitsT[vari]

        return (dates,depths,lons,lats,maxdepths,values,units)

    def reduce_ices(self,dataV0,rd,timeint,zint,zmaxint,lonint,latint):
        #column keys
        tcol=0; zcol=1; latcol=2; loncol=3; zmaxcol=5; varcol=4
        #remove non-inf values
        vali=np.isfinite(dataV0[:,varcol])
        dataV=dataV0[vali,:]

        #initial filterering according to lat,lon,max.depth, and depth interval (if applicable)
        #time filter
        dates=np.array([rd+datetime.timedelta(0,dataV[r,tcol]) for r in range(len(dataV[:,1]))]) #reconstruct dates
        ti = (dates > timeint[0]) & (dates < timeint[1])

        #max.depth filter
        zmaxi=(dataV[:,zmaxcol]>=zmaxint[0]) | (dataV[:,zmaxcol]<=zmaxint[1])

        #depth-interval filter (if specified)
        zi = (dataV[:, zcol] >= zint[0]) & (dataV[:, zcol] <= zint[1])

        #reduce the data , such that the lat-lon filtering doesn't have to scan all the data.
        ind=ti & zmaxi & zi
        if sum(ind)==0:
            raise(Warning('no observations left after filtering for t,zmax and z'))
            return([])
        else:
            dataR=dataV[ind,:]

        #lat-lon filter
        latloni = (dataR[:, loncol] > lonint[0]) & (dataR[:, loncol] < lonint[1]) & \
                  (dataR[:, latcol] > latint[0]) & (dataR[:, latcol] < latint[1])

        dataR2=dataR[latloni,:]

        # screen outliers
        OLth=np.mean(dataR2[:,varcol])+np.std(dataR2[:,varcol])*3
        vali=np.where(dataR2[:,varcol]<=OLth)[0]
        dataR3=dataR2[vali,:]
        print (str(len(dataR2[:,0])-len(vali))+' outliers removed')

        var=dataR3[:,[tcol,zcol,loncol,latcol,zmaxcol,varcol]]
        return var

if __name__=='__main__':

    if len(sys.argv)>1:
        files=sys.argv[1].split(',')
    else:
        files = ['/home/onur/WORK/projects/GB/data/ices/lon-1-10_lat50-57/raw/2010-2014/BGC_data_block.pickle',
                 '/home/onur/WORK/projects/2013/gpmeh/sns144-GPMEH-PPZZ-P190529-fSG97dChl-sedPS/extract_skillC_sns144-GPMEH-P190529-fSG97dChl.2012-2013_zSB.nc']
        #files=[''/home/onur/WORK/projects/GB/data/ices/lon-1-10_lat50-57/raw/TS-2008-2015/TS_data_block.pickle'',
        #       '/home/onur/WORK/projects/2013/maecs/sns144-M180109-nFpBpr-Pbg2017-B180106-vsdetp4b1-AHm2-c06/extract_Mphysred_sns144-M180109-nFpBpr-Pbg2017-B180106-vsdetp4b1-AHm2-c06.12-13_zSB.nc']

    if len(files)>2:
        raise(ValueError('2 filenames must be provided'))

    if len(sys.argv)>2:
        ftypes=sys.argv[2].split(',')
    else:
        ftypes =['ICES', 'GETM-GPMEH']
        #ftypes = ['ESA-CCI', 'GETM-MAECS']

    if len(sys.argv) > 3:
        fnames = sys.argv[3].split(',')
        print('   fnames0:%s\n   fnames1:%s' % (fnames[0], fnames[1]))
    else:
        # fnames = ['ICES', 'GETM-MAECS']
        fnames = ['ICES', 'GETM-GPM.PPZZ-EH.GC']
        # fnames = ['ESA-CCI', 'GETM-MAECS']

    if len(sys.argv)>4:
        varnames=sys.argv[4].split(',')
    else:
        varnames=['Chl','NO3','NH4','DIP'] #,'DOXY']

    if len(sys.argv) > 5:
        yrs = sys.argv[5].split(',')
    else:
        yrs=[2012,2013]


    get_matchups(files,ftypes,fnames,varnames,yrs)