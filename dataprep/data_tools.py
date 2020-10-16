__author__ = 'onur'

import netCDF4
import datetime
import numpy as np
import os
from scipy import interpolate,interp

def get_colors(n,palette='GB'):
    import matplotlib.pyplot as plt
    # set the colormap (http://matplotlib.org/1.2.1/examples/pylab_examples/show_colormaps.html)
    if palette=='BW':
        colormap = plt.cm.gist_yarg
    elif palette=='GB':
        colormap = plt.cm.GnBu
    elif palette=='coolwarm':
        colormap = plt.cm.coolwarm
    elif palette=='bwr':
        colormap = plt.cm.bwr
    colorlist=[colormap(i) for i in np.linspace(0, 0.9, n)]
    return colorlist


def get_markers(n=0):
    from matplotlib.lines import Line2D
    #Line2D.markers

    # markers = []
    # for m in Line2D.markers:
    #     try:
    #         if len(m) == 1 and m != ' ':
    #             markers.append(m)
    #     except TypeError:
    #         pass

    #markers = markers + [
    markers=['o','s','*','D','p','^','v','<','>','d','x','+','x','|','_',
            r'$\lambda$',
            r'$\bowtie$',
            r'$\circlearrowleft$',
            r'$\clubsuit$',
            r'$\checkmark$']

    if (n==0) or (n> len(markers)):
        if n> len(markers):
            raise(Warning('Requested number of markers exceeds the number of registered markers. Markers will recycle'))
        return markers
    else:
        return markers[0:(n-1)]

def regrid_loop(lon,lat,var,newlon,newlat,extramask=False,proj=0,method='linear'):
    #do lateral (x-y) regridding for time and/or z
    print ('  regridding')
    values=var['val'][:]
    if len(np.shape(values))==2:
        vals=regrid(lon,lat,values,newlon[:],newlat[:],var['mv'],proj,extramask=extramask,method=method)
    if len(np.shape(values))==3:
        d0len=np.shape(values)[0]
        vals=np.ndarray(shape=(d0len,np.shape(newlon)[0],np.shape(newlon)[1]), dtype=float, order='F')
        for d0i in range(0,d0len):
            vals[d0i,:,:]=regrid(lon,lat,values[d0i,:,:],newlon[:],newlat[:],var['mv'],proj,extramask=extramask,method=method)
    if len(np.shape(values))==4:
        d0len=np.shape(values)[0]
        d1len=np.shape(values)[1]
        vals=np.ndarray(shape=(d0len,d1len,np.shape(newlon)[0],np.shape(newlon)[1]), dtype=float, order='F')
        print ('frame#',)
        for d0i in range(0,d0len):
            print (str(d0i),)
            for d1i in range(0,d1len):
                vals[d0i,d1i,:,:]=regrid(lon,lat,values[d0i,d1i,:,:],newlon[:],newlat[:],var['mv'],proj,extramask=extramask,method=method)
        print ('.')
    return(vals)

def regrid(lon,lat,values,lon2,lat2,mv_in,proj,mv_out=-9999.0,extramask=False,method='linear'):

    #convert them to 1-D arrays
    lon1=np.reshape(lon,lon.size)
    lat1=np.reshape(lat,lon.size)
    #val1=np.reshape(values,lon.size)
    val1=np.ma.masked_equal(np.reshape(values,lon.size),mv_in)

    #trim the spatially irrelevant portion
    ikeep=np.where((lon1>=lon2.min()-0.5) & (lon1<=lon2.max()+0.5) &
                   (lat1>=lat2.min()-0.5) & (lat1<=lat2.max()+0.5) &
                   np.invert(val1.mask))

    lon1=lon1[ikeep]
    lat1=lat1[ikeep]
    val1=val1[ikeep]

    if val1.size==0:
        raise ValueError('no overlap between the available & requested spatial range')

    #re-grid based on cartesian distances
    x1,y1=proj(lon1,lat1)
    coords = np.asarray(zip(x1,y1))
    if method=='NN':
        f = interpolate.NearestNDInterpolator(coords, val1)
    elif method=='linear':
        f = interpolate.LinearNDInterpolator(coords, val1,fill_value=mv_out)
    x2,y2=proj(lon2,lat2)
    res = f(x2,y2)

    #replace missing/invalid values with the specified mv_out
    maskdef=np.zeros(res.shape,dtype=bool)
    if hasattr(res,'mask'):
        maskdef=res.mask
    res[np.where((maskdef) | (np.isnan(res)) | (res==mv_in) | (extramask) )]=mv_out

    return np.ma.masked_equal(res,mv_out)

def plot_fields(lon,lat,vals,varnames,longnames,units,fname_out,proj,tv=[0],zv=[0]):
    import matplotlib.pyplot as plt
    from postplib3D import plot2Dmap

    if len(vals[0].shape)==2:
        lateralshape=vals[0][:,:].shape
    elif len(vals[0].shape)==3:
        lateralshape = vals[0][0,:, :].shape
    elif len(vals[0].shape) == 4:
        lateralshape = vals[0][0,0, :, :].shape

    if lon.shape == lateralshape:
        #re-define the lat-lon at boundaries, assuming they are originally provided at centers
        dlon=np.diff(lon,axis=1)/2.
        lonB=np.ndarray((lon.shape[0]+1,lon.shape[1]+1))
        for r in range(lon.shape[0]): #i.e., for each row (lat)
            lonB[r,:]=np.concatenate((np.array([lon[r,0]-dlon[r,0]]),lon[r,:-1]+dlon[r,:],np.array([lon[r,-1]+dlon[r,-1]])))
        lonB[r+1,:]=lonB[r,:]
        dlat=np.diff(lat,axis=0)/2.
        latB=np.ndarray((lat.shape[0]+1,lat.shape[1]+1))
        for c in range(lat.shape[1]): #i.e., for each column (lon
            latB[:,c]=np.concatenate((np.array([lat[0,c]-dlat[0,c]]),lat[:-1,c]+dlat[:,c],np.array([lat[-1,c]+dlat[-1,c]])))
        latB[:,c+1]=latB[:,c]
        x,y=proj(lonB,latB)
    else:
        x,y = proj(lon,lat)

    fh=8;fw=8;dpi=300;left=0.1;right=0.9;bottom=0.1;top=0.9; hsp=0.25; wsp=0.25

    for vari,varname in enumerate(varnames):
        if len(vals[vari].shape) == 2:
            print ('  ' + varname)
            f = plt.figure(figsize=(fw, fh), dpi=dpi)
            f.subplots_adjust(left=left, right=right, bottom=bottom, top=top, hspace=hsp, wspace=wsp)
            f.text(0.5, 0.96, longnames[vari] + ' [' + units[vari] + ']', horizontalalignment='center')

            ax = plt.subplot(1,1,1)
            vmax = np.nanmax(vals[vari])
            plot2Dmap(f, ax, [0, vmax], x[:, :], y[:, :], vals[vari][:, :], proj, titlestr=varname,
                      cbarlabel='')  # ,cmap='YlOrRd',cbarpos='right',cbarlabel='Depth [m]'

            filename = os.path.join(os.path.dirname(fname_out),
                                    os.path.basename(fname_out)[:-3] + '_' + varname + '.png')
            plt.savefig(filename, dpi=dpi)
            plt.close()
            print ('.\n' + filename)
        elif len(vals[vari].shape)==3:
            print ('  '+varname+',frame #',)
            f = plt.figure(figsize=(fw,fh), dpi=dpi)
            f.subplots_adjust(left=left,right=right,bottom=bottom,top=top, hspace=hsp, wspace=wsp)
            f.text(0.5,0.96,longnames[vari]+' ['+units[vari]+']', horizontalalignment='center')

            for ti,t in enumerate(tv):
                print (ti+1,)
                if type(t)==np.float64:
                    title='month '+str(int(t))
                else:
                    title=t.date()
                ax=plt.subplot(np.ceil(len(tv)/3.0),3,ti+1)
                vmax=np.nanmax(vals[vari])
                plot2Dmap(f,ax,[0,vmax],x[:,:],y[:,:],vals[vari][ti,:,:],proj,titlestr=title,cbarlabel='') #,cmap='YlOrRd',cbarpos='right',cbarlabel='Depth [m]'

            filename=os.path.join(os.path.dirname(fname_out),os.path.basename(fname_out)[:-3]+'_'+varname+'.png')
            plt.savefig(filename,dpi=dpi)
            plt.close()
            print ('.\n'+filename)
        elif len(vals[vari].shape)==4:
            if len(zv)!=vals[vari].shape[1]:
                raise(Exception('var.shape[1] inconsistent with len(depths)'))

            print ('  '+varname+'\nlayer# ',)
            for zi,z in enumerate(zv):
                print (str(zi+1)+', frame #',)
                f = plt.figure(figsize=(fw,fh), dpi=dpi)
                f.subplots_adjust(left=left,right=right,bottom=bottom,top=top, hspace=hsp, wspace=wsp)
                f.text(0.5,0.97,longnames[vari]+' '+units[vari], horizontalalignment='center')
                f.text(0.5,0.94,'depth='+str(z), horizontalalignment='center')

                for ti,t in enumerate(tv):
                    print (ti+1,)
                    if (type(t)==np.float64) or (type(t)==np.int64):
                        title='month '+str(int(t))
                    else:
                        title=datetime.date(t)
                    ax=plt.subplot(np.ceil(len(tv)/3.0),3,ti+1)
                    #plot2Dmap(f,ax,[0,np.amax(vals[vari])/5.],x[:,:],y[:,:],vals[vari][ti,:,:],proj,titlestr=title,cbarlabel='') #,cmap='YlOrRd',cbarpos='right',cbarlabel='Depth [m]'

                    vmax= np.nanmax(vals[vari][:,zi,:,:])
                    #std=np.nanstd(vals[vari][:,zi,:,:]) #(nanstd)
                    #vmax=vmax/std? todo
                    plot2Dmap(f,ax,[0,vmax],x[:,:],y[:,:],vals[vari][ti,zi,:,:],proj,titlestr=title,cbarlabel='') #,cmap='YlOrRd',cbarpos='right',cbarlabel='Depth [m]'
                Zsuf='_Z'+str(zi)
                filename=os.path.join(os.path.dirname(fname_out),os.path.basename(fname_out)[:-3]+'_'+varname+Zsuf+'.png')
                plt.savefig(filename,dpi=dpi)
                plt.close()
                print ('.\n'+filename)

def mapon_container_tvec(Min):
    tvecAll=[]
    for i in range(len(Min)):
        tvecAll.extend(Min[i][:,0])

    tvec=np.unique(tvecAll)
    Mout=[None]*len(Min)
    for i in range(len(Min)):
        Mout[i]=np.ones((len(tvec),1))*np.nan
        for tno,t in enumerate(Min[i][:,0]):
            Mout[i][t==tvec]=Min[i][tno,1]

    return Mout,tvec


def mapon_container_tzvec(Min):
    tvecAll=[]
    for i in range(len(Min)):
        tvecAll.extend(Min[i][:,0])
    tvec=np.unique(tvecAll)

    zvecAll=[]
    for i in range(len(Min)):
        zvecAll.extend(Min[i][:,1])
    zvec=np.unique(zvecAll)

    Mout=[None]*len(Min)
    for i in range(len(Min)):
        Mout[i]=np.ones((len(tvec),len(zvec)))*np.nan
        for rno in range(len(Min[i][:,0])):
            t=Min[i][rno,0]
            z=Min[i][rno,1]
            Mout[i][t==tvec,z==zvec] = Min[i][rno,2]

    return Mout,tvec,zvec

def csv2nc(fin,fout,coords,deftime='none',refdate=0):
    print ('processing: %s'%(fin)) #,fout

    with open(fin, 'rb') as csvfile:
        RDraw=csvfile.readlines()

    colnames=RDraw[0].split(',')
    del RDraw[0] #delete the column name row
    units=RDraw[0].split(',')
    del RDraw[0] #delete the units row

    numcolO=len(colnames)
    numrow=len(RDraw)
    print (str(numrow)+' lines,'+str(numcolO)+' columns. ')

    #read time:
    if 'Date' in colnames:
        datec=colnames.index('Date')
    else:
        raise(ValueError('no "Date" column found in the input data set'))
    if 'Time' in colnames:
        if deftime=='none':
            timeav = True
            timec=colnames.index('Time')
        else:
            timeav = False
            hours=int(deftime.split(':')[0])
            mins=int(deftime.split(':')[1])
            secs=int(deftime.split(':')[2])
    else:
        timeav=False
        hours=12
        mins=0
        secs=0
        print ('no "Time" column found in the input data set. Assuming %s:%s:%s'%(hours,mins,secs))
    if 'Depth' in colnames:
        depthc=colnames.index('Depth')
        depthav=True
    else:
        depthav=False
        print ('no "Depth" column found in the input data set: Assuming 0 (surface).')

    tvec_s=np.ndarray((numrow,1)) #seconds
    depths=np.ndarray((numrow,1))
    for  r in np.arange(numrow):
        years=int(RDraw[r].split(',')[datec].split('-')[0])
        months=int(RDraw[r].split(',')[datec].split('-')[1])
        days = int(RDraw[r].split(',')[datec].split('-')[2])
        if months>12:
            raise(ValueError('Time format (eg: %s) wrong in %s.Expected format:Y-M-D')%(RDraw[r].split(',')[datec],fin))

        if timeav:
            hours=int(RDraw[r].split(',')[timec].split(':')[0])
            mins=int(RDraw[r].split(',')[timec].split(':')[1])
            try:
                secs=int(RDraw[r].split(',')[timec].split(':')[2])
            except:
                secs=0

        if (r==0) and refdate==0:
            refdate=datetime.datetime(years,months,days,hours,mins,secs)-datetime.timedelta(days=1)
        deltadate=datetime.datetime(years,months,days,hours,mins,secs)- refdate
        tvec_s[r]=86400*deltadate.days +deltadate.seconds

        if depthav:
            depths[r]=float(RDraw[r].split(',')[depthc])
        else:
            depths[r]=0.0

    #read each variable into a time-value matrix
    V=colnames[:]
    del V[V.index('Date')]
    if timeav:
        del V[V.index('Time')]
    if depthav:
        del V[V.index('Depth')]
        valc=2
    else:
        valc = 1 #

    M = [None] *len(V)
    U = [None] * len(V)
    for varno,var in enumerate(V):
        print (var,)
        colno=colnames.index(var)

        V[varno]=var.replace('\n','') #remove problematic characters
        U[varno]=units[colno]
        if depthav:
            M[varno]=np.ndarray((numrow,3)) #seconds,depth, var
        else:
            M[varno] = np.ndarray((numrow, 2))  # seconds,depth, var

        for  r in range(numrow):
            M[varno][r,0] = tvec_s[r]
            if depthav:
                M[varno][r,1] = depths[r]

            val=RDraw[r].split(',')[colno].split('\n')[0]
            if val=='':
                v=np.nan
            elif val[0]=='<':
                v=0.0
            else:
                try: #if converting to float doesn't work, write nan.
                    v=float(val)
                except:
                    v=np.nan
            M[varno][r, valc]=v
        #get rid of empty lines
        M[varno]=M[varno][np.invert(np.isnan(M[varno][:,1])),:]

    #write in a ncdf file:
    if depthav:
        dims = {'t': 'time', 'z': 'depth', 'x': 'lon', 'y': 'lat'}
        Mv, tvec_s, zvec = mapon_container_tzvec(M)  # map all vars along a container time/z vector
    else:
        dims = {'t': 'time', 'x': 'lon', 'y': 'lat'}
        Mv,tvec_s=mapon_container_tvec(M) #map all vars along a container time/z vector
        zvec = 0

    tvec=[datetime.datetime(refdate.year,refdate.month,refdate.day,0,0,0)+datetime.timedelta(0,t_s) for t_s in tvec_s]
    lon=np.array([coords['lon']])
    lat=np.array([coords['lat']])
    #create_ncfile(fout, lon, lat, Mv, V, V, U, dims, tvec, zvec, climatology=False,
    #              refdate=refdate, missing_value=-99, notify=True)

    create_ncfile_0D(fout,tvec,-1,Mv,V,V,U,lat,lon,dims={'t':'time'},climatology=False,refdate=refdate,missing_value=-99,notify=True)

def get_fields_fromtable(data,varname_in,fp,VarnValue,varnames,units,dimlims,dims_in,dims_out,mv_in=0.0,mv_out=-9999.0):
    import pickle
    if os.path.exists(fp):
        print ('Using pickled reshaped data')
        with open(fp) as f:
            var=pickle.load(f)
    else:

        print ('Reshaping raw data:')
        Cvar=varnames.index(VarnValue)
        Ct=varnames.index(dims_in['t'])
        Cz=varnames.index(dims_in['z'])
        Clat=varnames.index(dims_in['y'])
        Clon=varnames.index(dims_in['x'])

        #reduce the processing time by reducing the data based on dimension-limits
        ind=(data[:,Ct]>dimlims['t'][0]) & (data[:,Ct]<dimlims['t'][1]) &\
            (data[:,Cz]>dimlims['z'][0]) & (data[:,Cz]<dimlims['z'][1]) &\
            (data[:,Clat]>dimlims['y'][0]) & (data[:,Clat]<dimlims['y'][1]) &\
            (data[:,Clon]>dimlims['x'][0]) & (data[:,Clon]<dimlims['x'][1])
        dataR=data[ind,:]

        var={'unit':units[Cvar],'longname':varname_in,'mv':mv_out}

        for dimi,dim_out in enumerate(dims_out):
            #print dims_in[dimi],
            var[dim_out]=np.unique(dataR[:,varnames.index(dims_in[dim_out])])

        val=np.ndarray((len(var['t']),len(var['z']),len(var['y']),len(var['x']) ))*np.nan
        i=0;th=0.1;tot=len(dataR[:,1]);
        print ('% complete:',)
        for r in range(tot):
            #print 'r:'+str(r)+' th:'+str(th*tot)
            if r>(th*tot):
                print (str(th*100),)
                th=th+0.1

            ti=np.where(var['t']==dataR[r,Ct])[0]
            zi=np.where(var['z']==dataR[r,Cz])[0]
            lati=np.where(var['y']==dataR[r,Clat])[0]
            loni=np.where(var['x']==dataR[r,Clon])[0]
            if (len(ti)==1) & (len(zi)==1) & (len(lati)==1) &(len(loni)==1) :
                val[ti,zi,lati,loni]=dataR[r,Cvar]
            else:
                print ('r: '+str(r))

        var['val']=val

        if 't' in dims_in:
            if units[Ct]=='climatological month':
                var['time_date']=[datetime.date(2000,int(mon),15) for mon in var['t']]
            else:
                utime=netcdftime.utime(varnames==dims_in['time'])
                var['time_date']=utime.num2date(data[:,Ct])

        with open(fp, 'w') as f:
            pickle.dump(var,f)
            print (' Pickled the data.')

    return(var)

def get_ncfields(fname,varname_mult,mv_in=-99,mv_out=-99):
    #replaces 'get_fields_fromfile

    import netcdftime

    nc=netCDF4.Dataset(fname)
    varnamel=varname_mult.split(' ')
    if len(varnamel)>1:
        mult = float(varnamel[0])
        varname=varnamel[1]
    else:
        mult = 1.0
        varname=varnamel[0]
    if ('+' in varname):
        varnames=varname.split('+')
        val=0
        for vari,varnamei in enumerate(varnames):
            vali=nc.variables[varnamei][:]*mult
            val=val+vali
        #val[vali==mv_in]=mv_out
        val[(val>(mv_in-1e-5)) * (val<(mv_in+1e-5)) *(val<=-98)]=np.nan
    else:
        val=nc.variables[varname][:]*mult
        #val[val==mv_in]=mv_out
        if type(val) != np.float64:
            val[(val<=-98) | (val==mv_in)]=np.nan
    try:
        unit=nc.variables[varname].units
    except:
        unit='1'
    try:
        longname=nc.variables[varname].long_name
    except:
        longname=varname

    #attempt to include 'level', if exists
    if 'h' in nc.variables.keys():
        var={'val':val,'unit':unit,'longname':longname,'mv':mv_out, 'h':nc.variables['h'][:]}
    else:
        var={'val':val,'unit':unit,'longname':longname,'mv':mv_out}

    #get dimension names and values
    dimns=nc.variables[varname].dimensions
    dims={}
    for dim in dimns:
        if dim in ['t', 'time', 'time_counter']:
            dimout='t'
        elif dim in ['z', 'depth', 'level']:
            dimout='z'
        elif dim in ['x', 'longitude', 'lons', 'lon', 'xic']:
            dimout='x'
        elif dim in ['y', 'latitude', 'lats', 'lat', 'etac']:
            dimout = 'y'

        try:
            dims[dimout] = dim
        except:
            raise(ValueError('unknown dimension type:'+dim))

        # insert the dim values, either as provided in the input file, or just as counters
        if dim not in nc.variables: #just a counter
            var[dimout] = np.arange(1, len(nc.dimensions[dim]) + 1)
        else:
            if dimout=='t': #if time, only after appropriate transformation of the data
                if nc.variables[dim].units=='climatological month':
                    var[dimout] = [datetime.date(2000, int(mon), 15) for mon in nc.variables[dim][:]]
                else:
                    utime=netcdftime.utime(nc.variables[dim].units)
                    var[dimout]=utime.num2date(list(nc.variables[dim][:]))
            else:
                var[dimout] = nc.variables[dim][:]

    nc.close()
    return(var,dims)

def get_fields_fromfile(fname,varname_mult,dims,mv_in=-99,mv_out=-99):
    import netcdftime

    nc=netCDF4.Dataset(fname)
    varnamel=varname_mult.split(' ')
    if len(varnamel)>1:
        mult = float(varnamel[0])
        varname=varnamel[1]
    else:
        mult = 1.0
        varname=varnamel[0]
    if ('+' in varname):
        varnames=varname.split('+')
        val=0
        for vari,varnamei in enumerate(varnames):
            vali=nc.variables[varnamei][:]*mult
            val=val+vali
        #val[vali==mv_in]=mv_out
        val[(val>(mv_in-1e-5)) * (val<(mv_in+1e-5)) *(val<-98)]=np.nan
    else:
        val=nc.variables[varname][:]*mult
        #val[val==mv_in]=mv_out
        val[(val<=-98) | (val==mv_in)]=np.nan

    try:
        unit=nc.variables[varname].units
    except:
        unit='1'
    try:
        longname=nc.variables[varname].long_name
    except:
        longname=varname

    var={'val':val,'unit':unit,'longname':longname,'mv':mv_out}

    #print '  dims:(',
    for dim_out,dim_in in dims.iteritems():
        if dim_out=='t':
            if nc.variables[dim_in].units=='climatological month':
                var['t'] = [datetime.date(2000, int(mon), 15) for mon in nc.variables[dim_in][:]]
            else:
                utime=netcdftime.utime(nc.variables[dim_in].units)
                var['t']=utime.num2date(nc.variables[dim_in][:])
        else:
            if dim_in in nc.variables:
                var[dim_out]=nc.variables[dim_in][:]
            else:
                var[dim_out]=np.arange(1,len(nc.dimensions[dim_in])+1)
    #print ')'
    #lat=nc.variables[dims[3]][:]
    #lon=nc.variables[dims[4]][:]

    #var={'time':time,'time_date':time_date,'lon':lon,'lat':lat,'val':val,'unit':unit,'longname':longname,'mv':-1}
    nc.close()
    return(var)

def update_dimvals(dimvals,dims,var):
    for dimout,dimin in dims.items():
        if dimout not in dimvals:
            dimvals[dimout]=var[dimout]
    return(dimvals)


def create_ncfile(fname,lon,lat,vals,names,longnames,units,dims,tvec=[0],zvec=[0],climatology=False,refdate=0,missing_value=-99,notify=False):
    if refdate==0:
        refdate=datetime.datetime(1960,1,1,0,0,0)

    nc = netCDF4.Dataset(fname,'w')


    #DIMENSIONS

    #time
    if 't' in dims:
        nc.createDimension(dims['t'],size=None)
        nctime=nc.createVariable(dims['t'],float,(dims['t']), fill_value=missing_value)
        if climatology:
            months=[t.month for t in tvec]
            nctime[:]=months
            nctime.units='climatological month'
        else:
            deltadate=[t-refdate for t in tvec]
            tvec_s=[86400*deltadate[r].days +deltadate[r].seconds for r in range(len(deltadate))]
            nctime[:]=tvec_s
            nctime.units='seconds since '+str(refdate)

    #z
    if 'z' in dims:
        nc.createDimension(dims['z'],len(zvec))
        nczax=nc.createVariable(dims['z'],float,(dims['z']), fill_value=missing_value)
        nczax[:]=zvec
        nczax.units='meters'

    #lat-lon
    if 'y' in dims:
        nc.createDimension('y',lat.shape[0])
        #nc.createDimension('y',lon.shape[0])
    if 'x' in dims:
        if len(lon.shape)==1:
            nc.createDimension('x',lon.shape[0])
        elif len(lon.shape)==2:
            nc.createDimension('x',lon.shape[1])

    if 'y' in dims:
        if len(lat.shape)==1:
            nclat=nc.createVariable(dims['y'],float,('y'), fill_value=missing_value)
        elif len(lat.shape)==2:
            #nclat=nc.createVariable(dims['y'],float,('y','x'), fill_value=missing_value)
            nclat=nc.createVariable(dims['y'],float,('y','x'), fill_value=missing_value)
        nclat[:]=lat
        nclat.units = 'degrees_north'

    if 'x' in dims:
        if len(lon.shape)==1:
            nclon=nc.createVariable(dims['x'],float,('x'), fill_value=missing_value)
        elif len(lon.shape)==2:
            #nclat=nc.createVariable(dims['y'],float,('y','x'), fill_value=missing_value)
            nclon=nc.createVariable(dims['x'],float,('y','x'), fill_value=missing_value)
        nclon[:]=lon
        nclon.units = 'degrees_east'

    #print 'filling in:'
    #for varno in range(0,len(names)):
    for varno,varname in enumerate(names):
        #varname = names[varno]

        #print names[varno],
        if 't' in dims:
            if 'z' in dims:
                if 'y' in dims:
                    #print '(t z y x)'
                    ncvar=nc.createVariable(varname,float,(dims['t'],dims['z'],'y','x'), fill_value=missing_value)
                    ncvar.coordinates=dims['t'] +' '+ dims['z'] +' '+ 'y' +' '+ 'x'
                else :
                    #print '(t z y x)'
                    ncvar=nc.createVariable(varname,float,(dims['t'],dims['z']), fill_value=missing_value)
                    ncvar.coordinates=dims['t'] +' '+ dims['z']
                #ncvar[0:len(tvec_s),:,:]=vals[varno]
            else:
                #print '(t y x)'
                if 'y' in dims:
                    ncvar=nc.createVariable(varname, float, (dims['t'], 'y', 'x'), fill_value=missing_value)
                    ncvar.coordinates= dims['t'] +' '+ 'y' +' '+ 'x'
                else:
                    ncvar=nc.createVariable(varname, float, (dims['t']), fill_value=missing_value)
                    ncvar.coordinates= dims['t']
                #ncvar[0:len(tvec_s),:,:]=vals[varno]
        else:
            if 'z' in dims:
                if 'y' in dims:
                    ncvar=nc.createVariable(varname,float,(dims['z'],'y','x'), fill_value=missing_value)
                    ncvar.coordinates=dims['z'] +' '+ 'y' +' '+ 'x'
                else:
                    ncvar=nc.createVariable(varname,float,(dims['z']), fill_value=missing_value)
                    ncvar.coordinates=dims['z']
            else:
                #print '(y x)'
                ncvar=nc.createVariable(varname,float,('y','x'), fill_value=missing_value)
                ncvar.coordinates='y' +' '+ 'x'

        ncvar[:]=vals[varno]
        ncvar.units=units[varno]
        if longnames[varno]!='':
            ncvar.longname=longnames[varno]

    nc.sync()
    nc.institution = 'Helmholtz-Zentrum Geesthacht'
    nc.contact = 'onur.kerimoglu@hzg.de'
    nc.history = str(datetime.datetime.now()) + ': first creation'
    nc.close()
    if notify:
        print ('data written:'+fname)

def create_ncfile_0D(fname,tvec,zvec,vals,names,longnames,units,lats,lons,dims,climatology=False,refdate=0,missing_value=-99,notify=False):
    if refdate==0:
        refdate=datetime.datetime(2000,1,1,0,0,0)
    nc = netCDF4.Dataset(fname,'w')

    #DIMENSIONS

    # time
    if 't' in dims:
        nc.createDimension(dims['t'], size=None)
        nctime = nc.createVariable(dims['t'], float, (dims['t']), fill_value=missing_value)
        if climatology:
            months = [t.month for t in tvec]
            nctime[:] = months
            nctime.units = 'climatological month'
        else:
            deltadate = [t - refdate for t in tvec]
            tvec_s = [86400 * deltadate[r].days + deltadate[r].seconds for r in range(len(deltadate))]
            nctime[:] = tvec_s
            nctime.units = 'seconds since ' + str(refdate)

    #z
    if 'z' in dims:
        nc.createDimension(dims['z'],len(zvec))
        nczax=nc.createVariable(dims['z'],float,(dims['z']), fill_value=missing_value)
        nczax[:]=zvec

    #print 'filling in:'
    for varno in range(0,len(names)):
        #print names[varno],
        if 't' in dims:
            if 'z' in dims:
                #print '(t z)'
                ncvar=nc.createVariable(names[varno],float,(dims['t'],dims['z']), fill_value=missing_value)
                ncvar.coordinates=dims['t'] +' '+ dims['z']
                #ncvar[0:len(tvec_s),:,:]=vals[varno]
            else:
                #print '(t)'
                ncvar=nc.createVariable(names[varno],float,(dims['t']), fill_value=missing_value)
                ncvar.coordinates=dims['t']
        else:
            if 'z' in dims:
                #print '(z, y, x)'
                ncvar=nc.createVariable(names[varno],float,(dims['z']), fill_value=missing_value)
                ncvar.coordinates=dims['z']
            else:
                raise(Exception('unknown dimension structure:',dims))

        try:
            ncvar[:]=vals[varno]
        except:
            print ('!')
        ncvar.units=units[varno]
        if longnames[varno]!='':
            ncvar.longname=longnames[varno]
        if len(lats)==1:
            ncvar.lat=lats[0]
            ncvar.lon = lons[0]
        else:
            ncvar.lat=lats[varno]
            ncvar.lon=lons[varno]

    nc.sync()
    nc.institution='Helmholtz-Zentrum Geesthacht'
    nc.contact = 'onur.kerimoglu@hzg.de'
    nc.history=str(datetime.datetime.now())+': first creation'
    nc.close()
    if notify:
        print ('data written:'+fname)

def create_nc(fname,dimvals,dimslist,vals,names,longnames,units,history='',refdate=datetime.datetime(2000,1,1,0,0,0),climatology=False,missing_value=-99,notify=False):
    #replaces 'create_ncfile' and 'create_ncfile_0D'

    nc = netCDF4.Dataset(fname,'w')

    #DIMENSIONS
    #construct an encompassing dimension dictionary
    dims={}
    for varno,varn in enumerate(names):
        dimsv=dimslist[varno]
        for dimin,dimout in dimsv.items():
            if dimin not in dims.keys():
                dims[dimin]=dimout

    #create std dims if in dict
    if 't' in dims.keys():
        nc.createDimension(dims['t'],size=None) #len(dimvals['t'])
    if 'z' in dims.keys():
        nc.createDimension(dims['z'],len(dimvals['z']))
    if 'y' in dims.keys():
        nc.createDimension(dims['y'],len(dimvals['y']))
    if 'x' in dims.keys():
        nc.createDimension(dims['x'],len(dimvals['x']))


    #VARIABLES
    for varno,varname in enumerate(names):
        #varname = names[varno]
        dims=dimslist[varno]

        #print names[varno],
        if 't' in dims:
            if 'z' in dims:
                if 'y' in dims:
                    #print '(t z y x)'
                    ncvar=nc.createVariable(varname,float,(dims['t'],dims['z'],dims['y'],dims['x']), fill_value=missing_value)
                    ncvar.coordinates=dims['t'] +' '+ dims['z'] +' '+ dims['y']+' '+ dims['x']
                else :
                    #print '(t z)'
                    ncvar=nc.createVariable(varname,float,(dims['t'],dims['z']), fill_value=missing_value)
                    ncvar.coordinates=dims['t'] +' '+ dims['z']

            else:
                if 'y' in dims:
                    # print '(t y x)'
                    ncvar=nc.createVariable(varname, float, (dims['t'], dims['y'], dims['x']), fill_value=missing_value)
                    ncvar.coordinates= dims['t'] +' '+ dims['y'] +' '+ dims['x']
                else:
                    # print '(t)'
                    try:
                        ncvar=nc.createVariable(varname, float, (dims['t']), fill_value=missing_value)
                    except:
                        print('!')
                    ncvar.coordinates= dims['t']
                #ncvar[0:len(tvec_s),:,:]=vals[varno]
        else:
            if 'z' in dims:
                if 'y' in dims:
                    # print '(z y x)'
                    ncvar=nc.createVariable(varname,float,(dims['z'],dims['y'],dims['x']), fill_value=missing_value)
                    ncvar.coordinates=dims['z'] +' '+ dims['y'] +' '+ dims['x']
                else:
                    # print '(z)'
                    ncvar=nc.createVariable(varname,float,(dims['z']), fill_value=missing_value)
                    ncvar.coordinates=dims['z']
            else:
                if ('y' in dims) and ('x' in dims):
                    ncvar=nc.createVariable(varname,float,(dims['y'],dims['x']), fill_value=missing_value)
                    ncvar.coordinates=dims['y'] +' '+ dims['x']
                elif 'y' in dims:
                    ncvar=nc.createVariable(varname,float,(dims['y']), fill_value=missing_value)
                    ncvar.coordinates=dims['y']
                elif 'x' in dims:
                    ncvar=nc.createVariable(varname,float,(dims['x']), fill_value=missing_value)
                    ncvar.coordinates=dims['x']
                else:
                    raise(ValueError('unknown dimension structure:'+dims))

        if varname in ['t','time']:
            #nctime=nc.createVariable(dims['t'],float,(dims['t']), fill_value=missing_value)
            if climatology:
                months=[t.month for t in dimvals['t']]
                ncvar[:]=months
                ncvar.units='climatological month'
            else:
                deltadate=[t-refdate for t in dimvals['t']]
                tvec_s=[86400*deltadate[r].days +deltadate[r].seconds for r in range(len(deltadate))]
                ncvar[:]=tvec_s
                ncvar.units='seconds since '+str(refdate)
        else:
            try:
                ncvar[:]=vals[varno]
            except:
                print ('!')
            ncvar.units=units[varno]
            if longnames[varno]!='':
                ncvar.longname=longnames[varno]

    nc.sync()

    #GLOBAL ATTRIBUTES
    nc.institution = 'Helmholtz-Zentrum Geesthacht'
    nc.contact = 'onur.kerimoglu@hzg.de'
    if history=='':
        nc.history = str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M")) + '- first creation'
    else:
        nc.history = str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M")) + '-'+ history

    nc.close()
    if notify:
        print ('data written:'+fname)

def get_botdepth(lon,lat,method='tree'):
    topofile = '/home/onur/WORK/projects/GB/data/topo/topo_HR2.nc'
    #topofile='/home/onur/WORK/projects/GB/data/topo/topo_GreaterNorthSea.nc'
    if os.path.exists(topofile):
        nc=netCDF4.Dataset(topofile)
        #depths=nc.variables['depth_average'][:]
        #lats = nc.variables['lat'][:]
        #lons = nc.variables['lon'][:]
        depths = nc.variables['bathymetry'][:]
        lats = nc.variables['latx'][:-1,:-1]
        lons = nc.variables['lonx'][:-1,:-1]
        nc.close()
        depth = interpval2D(lats, lons, depths, lat, lon, method)
    else:
        warnings.warn('Topo file not found in the provided location (%s). Filling nan for station depth:'%topofile)
        depth=np.nan
    return depth

def interpval2D(lats,lons,vals,lat,lon,method,domaintree=0):
    from postplib3D import getproj
    from scipy.spatial import cKDTree
    import time
    import warnings

    #meshscipy.spatial.ckdtree.cKDTree.query
    if len(lons.shape)==2:
        lonsM,latsM=lons,lats
    else:
        lonsM,latsM=np.meshgrid(lons,lats)

    #convert lat-lons to cartesian coordinates
    proj=getproj('NS')
    x1,y1 = proj(lonsM,latsM)
    x2,y2 = proj(lon, lat)

    if method=='pretree':
        if not istype(domaintree,integer):
            domaintree=domaintree
        else:
            warnings.WarningMessage('no pre-generated tree provided. Generating a new one')
            method='tree'

    if method == 'tree':
        xypairs = zip(x1.flat, y1.flat)
        # ts=time.time()
        print (' creating a tree for interpolation, this can take a while..',)
        domaintree = cKDTree(xypairs)
        # print 'tree generated in:%.1f' %(time.time() - ts)

    print (' Interpolating (%s)'%method,)
    if method in ['NN','linear']:
        coords = np.asarray(zip(np.reshape(x1,x1.size),np.reshape(y1,y1.size)))
        if method=='NN':
            f = interpolate.NearestNDInterpolator(coords, np.reshape(vals, vals.size))
        elif method=='linear':
            f = interpolate.LinearNDInterpolator(coords, np.reshape(vals,vals.size), fill_value=np.nan)
        val = f(x2, y2)
    elif method in ['tree','pretree']:
        val=interp_2d_tree_p2(vals, domaintree, x2, y2, k=4, kmin=1)
    print (' done.')

    return val

def get_2Dtree_p2(lons,lats):
    # create the query tree needed for 2-D interpolation
    from scipy.spatial import cKDTree
    if len(lons.shape)==2:
        lonsM,latsM=lons,lats
    else:
        lonsM,latsM=np.meshgrid(lons,lats)
    #lonlatpairs = list(zip(lons.flatten(), lats.flatten()))[:] # (python3 zip)
    #lonlatpairs= zip(lonsM.flatten(), latsM.flatten())
    lonlatpairs=zip(lonsM.flat, lonsM.flat)
    domaintree = cKDTree(lonlatpairs)
    return domaintree

def interp_2d_tree_p2(vals,domaintree,lon,lat,k=4,kmin=1):
    #vals: 2d field of values to be interpolated
    #domaintree: query tree generated with scipy.spatial.cKDTree
    #lon,lat:  coordinates to be used for interpolation
    #k: number of cells to use for interpolation

    dLi, gridindsLi = domaintree.query([lon, lat], k)
    wLi = 1.0 / dLi ** 2

    #remove the mask, if necessary
    if isinstance(vals,np.ma.MaskedArray):
        vals=vals.filled(np.nan)

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

def calc_osol(temp,sal,method='Weiss'):
    #Yoana's unreferenced formula
    #x = (temp + 273.16) / 100;
    #osol = np.exp(((-21.8492 * x - 173.4292) * x + 249.6339) +
    #       sal * ((-.0017 * x + 0.014259) * x - 0.033096) + 143.3484 * np.log(x))
    if method=='Weiss':
        # http://www.calcofi.org/references/347-data-algorithms.html
        A1 = -173.4292
        A2 = 249.6339
        A3 = 143.3483
        A4 = -21.8492
        B1 = -0.033096
        B2 = 0.014259
        B3 = -0.00170
        t = temp + 273.15
        osol = np.exp(A1 + A2 * (100. / t) + A3 * np.log(t / 100.) + A4 * (t / 100)
                         + sal * (B1 + B2 * (t / 100.) + B3 * (t / 100.)**2)) #ml/l
    elif method=='helcom':
        # General guidelines on quality assurance for monitoring in the Baltic Sea
        # http://www.helcom.fi/Documents/Action%20areas/Monitoring%20and%20assessment/Manuals%20and%20Guidelines/Manual%20for%20Marine%20Monitoring%20in%20the%20COMBINE%20Programme%20of%20HELCOM_PartB_AnnexB8_Appendix3.pdf
        osol_mumol_l = np.exp(-135.90205
           + (1.575701 * 10**5) / (temp + 273.15)
           - (6.642308 * 10**7) / (temp + 273.15)**2
           + (1.243800 * 10**10)/ (temp + 273.15)**3
           - (8.621949 * 10**11)/ (temp + 273.15)**4
           - sal * (0.017674 - 10.754 / (temp + 273.15) + 2140.7 / (temp + 273.15)**2))
        osol=osol_mumol_l*0.0223916 #ml/l
    else:
        raise(Exception('Unknown method (%s) for O2-solubility calculation'%method))
    return osol
