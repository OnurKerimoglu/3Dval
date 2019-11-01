#!/usr/bin/env python
#> @brief script for making 2-D plots of fabm (water and soil) variables
#> @author: kerimoglu.o@gmail.com

""" module plot_varmap_multi.py
example call from shell:

$ python plot_varmap_multi.py filename.nc var1,var2 Maverage

"""
#from pylab import *
import sys,os
sys.path.insert(1, '/h/kerimogl/local/lib/python2.7/site-packages')
import calendar
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import numpy as np
import netCDF4,netcdftime,os
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LogNorm
from getm_funcs import get_getm_bathymetry_cropped
from general_funcs import fnclean,cm2inch,getproj,discrete_cmap_tuner

knownunits={'total_chlorophyll_calculator_result':'mg/m$^{3}$','GPM_phy_Chl':'mg/m$^{3}$',
            'total_NPPR_calculator_result':'mgC/m$^2$/d','GPM_phy_NPPR':'mgC/m$^2$/d',
            'EH_abioP_o2o_pel':'mmol/m$^2$/d','EH_abioS_o2o_brm':'mmol/m$^2$/d',
            'EH_abioP_DINO3':'$\mu$M','EH_abioP_DINH4':'$\mu$M','EH_abioP_DIN':'$\mu$M',
            'EH_abioP_DIP':'$\mu$M','EH_abioP_DISi':'$\mu$M',
            'sigma_t':'kg/m$^3$','sigma0':'kg/m$^3$','rho':'kg/m$^3$', 'temp':u'\N{DEGREE SIGN}C', 'tempmean':u'\N{DEGREE SIGN}C', 'salt':'g/kg', 'saltmean':'g/kg'}
logvars=['']
primprodvars=['hzg_maecs_GPPR', 'hzg_maecs_NPPR','GPM_phy_NPPR','total_NPPR_calculator_result']

def do_2Dplotmap(fname, varnames,setup,VertMeth,TempMeth0,colmap,Nlev,mode,plottopo=True,datasource='GF'):

    if mode=='singlepanel': #single panel plots (eg, 2013 study)
        years2plot=0 #[2012,2013] #[2000,2001,2002,2003,2004,2005]
        months2plot= [] #[9,10,11]
        numcol0 = 3.0
        figh=8
        figw=8
        dpi=200
        plottopo= True
        showparmer=True
        left=0.12;right=0.76;bottom=0.1;top=0.9; hsp=0.25; wsp=0.1
    elif mode=='3panelsperrow': # 3 panels per row
        years2plot = 0  # [2012,2013] #[2000,2001,2002,2003,2004,2005]
        months2plot = [] # [4,5,6,7,8,9]
        numcol0 = 3.0
        figh = 18
        figw = 15
        dpi = 200
        plottopo= False
        showparmer = False
        left = 0.08;right = 0.8;bottom = 0.05;top = 0.88;hsp = 0.2; wsp = 0.2
    elif mode == 'halfpagewidth':  # half-width page style plot
        months2plot= [1,2,3,4,5,6,7,8,9,10,11,12]
        numcol0 = 2.0
        figw = 6.5  # 12
        figh = 12
        showparmer = False
        left=0.07;right=0.8;bottom=0.05;top=0.85; hsp=0.25; wsp=0.1


    if ('Y' in TempMeth0) and ('M' in TempMeth0):
        ystr=TempMeth0.split('-')[0]
        years2plot = [int(ystr.split('Y')[1])]

    multiYfile=True
    if os.path.basename(fname)[0:7]=='extract':
        multiYfile=False

    proj=getproj(setup=setup)
    
        # read the nc-file
    print ('opening: '+fname)
    nc=netCDF4.Dataset(fname)
    ncv=nc.variables

    #prepare a fname for the plot
    fnamedir=os.path.dirname(fname)
    plotdir=os.path.join(fnamedir,'varmaps')
    if not os.path.exists(plotdir):
        os.mkdir(plotdir)
    fnamebase = os.path.basename(fname).split('.nc')[0]
    fname=os.path.join(plotdir,fnamebase)
    
    #common variables
    H=0
    if datasource=='MERIS':
        lats=nc.variables['latitude'][:]
        lons=nc.variables['longitude'][:]
        H=nc.variables['bathymetry'][:]
        x,y=proj(lons,lats)
    elif datasource=='GF':
       ncvk=nc.variables.keys()
       if 'bathymetry' in ncvk and 'lon' in ncvk and 'lat' in ncvk:
           lons=nc.variables['lon'][:] 
           lats=nc.variables['lat'][:] 
           H=nc.variables['bathymetry'][:] 
       else:       
            topo=get_getm_bathymetry_cropped()
            lats=topo['lats']; lons=topo['lons']
            H = topo['H']
       x,y = proj(lons,lats)
   
    tv = nc.variables['time']
    utime=netcdftime.utime(tv.units)
    tvec=utime.num2date(tv[:])

    #extract days, months, years to be used for finding time indices later
    years=[tvec[ti].year for ti in range(0,len(tvec))]
    months=[tvec[ti].month for ti in range(0,len(tvec))]
    days=[tvec[ti].day for ti in range(0,len(tvec))]

    #find available years (will make one multi-panel plot per year)
    #years2plot=np.unique(years)
    if years2plot==0:
        years2plot=np.unique(years)

    print(str(len(days))+' scenes found in total, for '+ str(np.unique(years)))

    for yi,year in enumerate(years2plot):
        print ('plotting: '+ str(year),)

        #tind=np.arange(0,len(days),2) #every 3rd frame (when cut was done for every 10th time step we get ~1 plot/month)
        #tind=np.arange(1,11,2) #every 3rd frame (when cut was done for every 10th time step we get ~1 plot/month)

        # construct a time-index by finding the first day of each month
        yeari=np.where(np.in1d(years,year))
        if len(yeari[0])==0:
            raise ValueError('Specified year:'+str(year)+ ' not found in this data set')
        tvecY=tvec[yeari[0]]

        #parse the TempMeth0
        if ('Y' in TempMeth0) and ('M' in TempMeth0):
            TempMeth='YMaverage'
            mstr=TempMeth0.split('-')[1]
            months2plot = [int(mstr.split('M')[1])]
        elif '_' in TempMeth0:
            TempMeth=TempMeth0.split('_')[0]
            mstr=TempMeth0.split('_')[1]
            if mstr=='1-2-12':
                titlestr='Winter Average'
            elif mstr in ['3-4-5-6-7-8-9','4-5-6-7-8-9']:
                titlestr='Growing Season Average'
            elif mstr=='1-2-3-10-11-12':
                titlestr='Non-growing Season Average'
            elif mstr=='3-4-5':
                titlestr='Spring Average'
            elif mstr=='6-7-8':
                titlestr='Summer Average'
            else:
                titlestr='Months %s averaged'%mstr
            #months2plot=[int(mstr[i]) for i in range(len(mstr))]
            months2plot=[int(m) for m in mstr.split('-')]
        else:
            TempMeth=TempMeth0
            if months2plot == []:
                monthsY=[tvecY[ti].month for ti in range(0,len(tvecY))]
                uniquemonths=np.unique(monthsY)
                months2plot=uniquemonths

        if TempMeth in ['Yaverage','Yintegral']:
            showparmer=True
            tind=1*[None]
            if '_' not in TempMeth0: #take everything available that for the year
                tind[0]= np.where(np.in1d(years,year))[0][:]
            else:
                yri=[]
                for mi,month in enumerate(months2plot):
                    newi=np.where(np.in1d(months,month) & np.in1d(years,year))[0]
                    yri.extend(newi)
                tind[0]=yri
        elif TempMeth in ['snapshots','Maverage','YMaverage']:
            tind=len(months2plot)*[None]
            tindsuf='_'
            for mi,month in enumerate(months2plot):
                tindsuf+=str(month)
                moni=np.where(np.in1d(months,month) & np.in1d(years,year))
                #if TempMeth=='snapshots':
                #    tind[mi]=moni[0]
                #elif TempMeth in ['Maverage']:
                tind[mi]=moni[0][:]

        skipy=True
        for tindm in tind:
            if len(tindm)>0:
                skipy=False
        if skipy:
            print ('not enough data, skipping')
            continue

        tindsuf='-'+TempMeth0

        # construct a time-index by finding the prescribed days and months
        #months2plot=[6,7,8]
        #days2plot=[7,17,27]#
        #(tind,)=np.where(np.in1d(days,days2plot) & np.in1d(months,months2plot) & (years==year))

        if len(tind)<numcol0:
            numcol=float(len(tind))
        else:
            numcol=numcol0
        numrow=np.ceil(len(tind)/numcol)

        if figw==0 or figh==0:
            figw=numcol*4.; figh=numrow*3

        print ('\n Specified Methods: T:'+TempMeth + ' V:'+VertMeth)
        print (str(len(days))+' scenes found. plotting:' +str(tind))
        for varno,varname in enumerate(varnames.keys()):
            print (varname,)
            try:
                clim=varnames[varname] #[0,]
            except:
                clim=[varname+'mean']

            #read the values from nc file
            if varname == 'currs':
                if ('u' in ncv) and ('v' in ncv):
                    u = ncv['u'][:]
                    v = ncv['v'][:]
                    unitstr = ncv['u'].units
                elif ('uumean' in ncv) and ('vvmean' in ncv):
                    u = ncv['uumean'][:]
                    v = ncv['vvmean'][:]
                    unitstr = ncv['uumean'].units
                else:
                    print ('u and/or v was not found in dataset')
                longname = 'current speed'
            else:
                v,unitstr,longname=getncvar(ncv,varname,clim)

            if True: #(obsolete) i.e. not  (len(v.shape)==4) and (VertMeth in ['each','SB'])

                f = plt.figure(figsize=cm2inch((figw,figh)), dpi=dpi)
                f.subplots_adjust(left=left,right=right,bottom=bottom,top=top, hspace=hsp, wspace=wsp)

                if len(v.shape)==3:
                    #find a common clim, if not provided
                    if clim[1]==-9:
                        clim=[np.nanmin(v[:,:,:]),np.nanmax(v[:,:,:])]
                    print (' clim: '+ str(clim),)

                for mi,i in enumerate(tind):
                    ax=plt.subplot(np.ceil(len(tind)/numcol),numcol,mi+1)

                    # if len(v.shape)==3: #directly plot the variable
                    #     vI=v[i,:,:]
                    #     suffix=''
                    #     f.text(0.5,0.96,longname+' ['+unitstr+']', horizontalalignment='center')
                    #     #filename=fname.split('.nc')[0]+'_'+varname+suffix+'_'+ str(tvec[i].date()) + '.png' #pdf
                    #     #titlestr=suffix+' '+longname+'\n'+str(tvec[i].date())
                    #     titlestr=str(tvec[i].date())
                    #     scalesuf=plot2Dmap(f,ax,clim,x,y,vI,varname,proj,setup,titlestr,plottopo,H,cbarstr=unitstr)
                    #     filename=fname.split('.nc')[0]+'_'+str(year)+'_varmap'+'_'+varname+scalesuf+'.png' #pdf

                    if len(v.shape)==3:
                        suffix=''
                        if TempMeth=='snapshots':#directly plot the variable
                            #print 'scene-'+str(i),
                            if varname!='currs':
                                vI=v[i,:,:]
                                if varname=='t2':
                                    vI=vI-273.15
                            else:
                                uI=u[i,:,:]
                                vI=v[i,:,:]

                            titlestr=str(tvec[i][0].date())
                            if datasource=='MERIS':
                                month=months2plot[mi]
                                titlestr=str(calendar.month_name[month])
                                scalesuf=plot2Dmap(f,ax,clim,x,y,vI,varname,proj,setup,titlestr,plottopo,H,showparmer,unitstr,colmap,Nlev)
                            else:
                                if varname!='currs':
                                    scalesuf=plot2Dmap(f,ax,clim,x[:,:],y[:,:],vI[:,:],varname,proj,setup,titlestr,plottopo,H[:,:],showparmer,unitstr,colmap,Nlev)
                                else:
                                    scalesuf=plot2Dmap_Q(f,ax,clim,x[:,:],y[:,:],uI[:,:],vI[:,:],varname,proj,setup,titlestr,plottopo,H[:,:],showparmer,unitstr,colmap,Nlev)

                        elif TempMeth in ['Maverage','YMaverage', 'Yaverage', 'Yintegral']:

                            #print '\nM: '+str(mi)+':',
                            tindM=tind[mi]
                            vII=0; uII=0 #temporal integrand
                            for ti,i in enumerate(tindM):
                                #print 'scene-'+str(i),
                                #temporal integration
                                if varname!='currs':
                                    vI=v[i,:,:]
                                    if varname=='t2':
                                        vI=vI-273.15
                                    vII=vII+vI
                                else:
                                    uI=u[i,:,:]
                                    vI=v[i,:,:]
                                    uII=uII+uI
                                    vII=vII+vI

                            if  TempMeth in ['Maverage','YMaverage', 'Yaverage']: #transform to temporal average
                                vII=vII/len(tindM)
                                #print('divide by len(tindM)')
                                if varname=='currs':
                                    uII=uII/len(tindM)
                                if TempMeth0=='Yaverage' :
                                    titlestr='Annual Average'
                                elif TempMeth == 'YMaverage':
                                    titlestr='%s %s average'%(year,calendar.month_name[month])
                            elif TempMeth in ['Yintegral']:
                                unitstr=unitstr.replace('d','y')
                                titlestr='Annual Integral'
                                #titlestr=unitstr
                                if varname in primprodvars:
                                    vII=vII/1000*12 #[mmolC/1000molC*12gC/molC]
                                    unitstr=unitstr.replace('mmolC','gC')
                                    titlestr='GPPR ['+unitstr+'] ('+str(tvec[i].year)+')'
                                if varname=='hzg_medmac_denit_rate':
                                    titlestr='Den. Rate ['+unitstr+'] ('+str(tvec[i].year)+')'

                            if TempMeth in ['Maverage']:
                                titlestr=str(calendar.month_name[months2plot[mi]])

                            if varname!='currs':
                                scalesuf=plot2Dmap(f,ax,clim,x[:,:],y[:,:],vII[:,:],varname,proj,setup,titlestr,plottopo,H[:,:],showparmer,unitstr,colmap,Nlev)
                            else:
                                scalesuf=plot2Dmap_Q(f,ax,clim,x[:,:],y[:,:],uII[:,:],vII[:,:],varname,proj,setup,titlestr,plottopo,H[:,:],showparmer,unitstr,colmap,Nlev)
                    
                    f.text(0.5,0.96,suffix+longname+' ['+unitstr+'] ('+str(tvec[i].year)+')', horizontalalignment='center')
                    if multiYfile:
                        filename=fname.split('.nc')[0]+'_'+str(year)+'_varmap'+suffix+'_'+varname+tindsuf+'-'+setup+ scalesuf+'.png' #pdf
                    else:
                        filename=fname.split('.'+str(year)+'.nc')[0]+'_'+str(year)+'_varmap'+suffix+'_'+varname+tindsuf+'-'+setup+ scalesuf+'.png' #pdf
                filename = fnclean(filename)
                plt.savefig(filename,dpi=dpi)
                print ('\nsaved:'+filename)
                #s=show()
                plt.close()

    #close the netcdf file
    nc.close()

def getncvar(ncv, varname0,clim):

    #if varname != 'currs':
    unitstr=''
    if varname0 in knownunits.keys():
        unitstr = knownunits[varname0]
    
    varname=varname0
    if varname0 == 'EH_abioP_DIN':
        varname = 'EH_abioP_DINO3+EH_abioP_DINH4'
    if varname0 == 'EH_abioP_o2o_pel':
        varname = 'EH_abioP_o2o_bac+EH_abioP_o2o_n4n'

    if (varname in ncv) or (varname + 'mean' in ncv):
        try:
            v = ncv[varname][:]
            if unitstr == '':
                unitstr = ncv[varname].units
        except:
            v = ncv[varname + 'mean'][:]
            if unitstr == '':
                unitstr = ncv[varname+'mean'].units
        try:
            longname = ncv[varname].long_name.replace('_', ' ')
        except:
            longname = varname
    elif '/' in varname:
        vn1, vn2 = varname.split('/')
        v = ncv[vn1][:] / ncv[vn2][:]
        longname = vn1 + '/' + vn2
        if unitstr == '':
            if ncv[vn1].units == ncv[vn2].units:
                unitstr = '[-]'
            else:
                unitstr = ncv[vn1].units + '/' + ncv[vn2].units
    elif '*' in varname:
        vn1, vn2 = varname.split('*')
        v = ncv[vn1][:] * ncv[vn2][:]
        if unitstr=='':
            unitstr = ncv[vn1].units + '*' + ncv[vn2].units
        longname = ncv[vn1].long_name.replace('_', ' ') + '*' + ncv[vn2].long_name.replace('_', ' ')
    elif '+' in varname:
        vn1, vn2 = varname.split('+')
        v = ncv[vn1][:] + ncv[vn2][:]
        if unitstr == '':
            unitstr = ncv[vn1].units + '+' + ncv[vn2].units
        longname = ncv[vn1].long_name.replace('_', ' ') + '+' + ncv[vn2].long_name.replace('_', ' ')
    else:
        raise (ValueError(varname + 'was not found in dataset'))

    if unitstr == '-':
        unitstr = ''  # '[-]'
    #if '^3' in unitstr:
    #    unitstr = unitstr.replace('^3', '$^3$')
    if varname0 in primprodvars and clim[1] > 0:
        v = v * 12  # [mmolC/1000molC*12gC/molC]
        unitstr = unitstr.replace('mmolC', 'mgC')
        print('primary production conversion (mmolC -> mgC)')

    return (v,unitstr,longname)

def plot2Dmap_Q(f,ax,clim,x,y,u,v,varname,proj,setup,titlestr,plottopo,H,showparmer=False,unitstr='',colmap='viridis',Nlev=7):
    #dim(latx)=dim(lonx)=98,139
    landgr = 0.8  # color of land
    plotcont= False

    # mask the nan values
    u = np.ma.masked_array(u, mask=np.isnan(u))
    v = np.ma.masked_array(v, mask=np.isnan(v))
    speed = np.sqrt(u ** 2 + v ** 2)

    # plot
    if plotcont:
        maxval = np.max(np.max(speed))
        minval = np.min(np.min(speed))
        # print 'min:%s-max:%s'%(minval,maxval)

        if clim[0] == clim[-1]:
            clim = [np.floor(minval), np.ceil(maxval)]

        cmap, extopt, cbt, intbounds, allbounds = discrete_cmap_tuner(clim, [minval, maxval], Nlev=Nlev, colmap=colmap,
                                                                      nancol='white')
        cmap.set_bad(str(landgr))
        pcf = proj.pcolormesh(x, y, speed, cmap=cmap, vmin=intbounds[0], vmax=intbounds[-1])

    #velocity field
    r=2
    Q=proj.quiver(x[::r,::r],y[::r,::r], u[::r,::r], v[::r,::r],
                  width=0.005,units='width',scale=6)
    mpl.pyplot.quiverkey(Q,0.48,0.2, 0.2,'20 cm s$^{-1}$', labelpos='E', coordinates='figure',
                                        fontproperties={'size':10})
    #mpl.pyplot.quiverkey(Q,0.4,0.10, 0.1,'10 cm s$^{-1}$', labelpos='E', coordinates='axes',
    #                                    fontproperties={'size':9})

    if plottopo:
        proj.contour(x,y,H,colors=str(landgr),linewidths=0.5,levels=[5.0,10.,20.,30.,40.,50.,60.])

    plt.title(titlestr,size=10)

    #setup specific features
    if setup in ['NSBS','NSBSfull','SNSfull','SNS','GBight','WadSea']:
        #coastlines, etc
        #proj.drawcoastlines(color=(0.3,0.3,0.3),linewidth=0.5)
        proj.fillcontinents((landgr, landgr, landgr), lake_color=(landgr, landgr, landgr))
        if showparmer:
            proj.drawparallels(np.arange(51., 57., 1.), labels=[1, 0, 0, 0], fontsize=9)
            proj.drawmeridians(np.arange(-1., 10., 1.), labels=[0, 0, 0, 1], fontsize=9)
        else:
            proj.drawparallels(np.arange(51., 57., 1.), labels=[1, 0, 0, 0], fontsize=9,linewidth = 0)
            proj.drawmeridians(np.arange(-1., 10., 1.), labels=[0, 0, 0, 1], fontsize=9,linewidth = 0)

    # plot
    if plotcont:
        # retrieve the axes position to set the colorbar position
        pos1 = ax.get_position()
        poscbar = [pos1.x0 + pos1.width + 0.07, pos1.y0 + 0.09, 0.02, pos1.height * 0.7]
        cax = plt.axes(position=poscbar)

        #cb=plt.colorbar(pcf)
        #plt.setp(cb.ax.get_yticklabels()[::2], visible=False)
        norm = mpl.colors.BoundaryNorm(intbounds, cmap.N)
        cb = mpl.colorbar.ColorbarBase(cax, cmap=cmap,
                                       norm=norm,
                                       boundaries=allbounds,
                                       extend=extopt,
                                       ticks=cbt,
                                       # ticks= cbt[0::2],
                                       # ticks= cbt[1:-1:2],
                                       spacing='proportional')
        cb.formatter.set_scientific(True)
        cb.formatter.set_powerlimits((-4, 4))
        cb.update_ticks()
        scalesuf = '_CL_%.1f_%.1f' % (clim[0], clim[1])
        scalesuf = scalesuf.replace('.0', '')

        # unitstr
        unitstr = unitstr.replace('**', '^')
        # print unitstr
        if unitstr != '':
            cb.ax.set_title(unitstr, size=10.)
        else:
            cb.ax.set_title('[-]', size=10.)
    else:
        scalesuf=''

    return scalesuf

#plot2Dmap(f,ax,clim,x,y,v,proj,titlestr,=titlestr) #[cmap='YlOrRd',cbarpos='right',cbarlabel='Depth [m]',showparmer='True'] not there: plottopo,H,varname,setup,
def plot2Dmap(f,ax,clim,x,y,v,varname,proj,setup,titlestr,plottopo,H,showparmer=False,unitstr='',colmap='viridis',Nlev=7):
    #dim(latx)=dim(lonx)=98,139
    landgr = 0.8  # color of land

    #mask the nan values
    v=np.ma.masked_array(v,mask=np.isnan(v))

    #plot
    maxval = np.max(np.max(v))
    minval = np.min(np.min(v))
    #print 'min:%s-max:%s'%(minval,maxval)

    if clim[0] == clim[-1]:
        clim=[np.floor(minval), np.ceil(maxval)]

    if colmap in ('bwr','Blues_r','Reds'): #(clim[0]*-1)==(clim[1]):
        unitstr = '$\%\Delta$'

    if unitstr=='degC' and clim[1]<15:
        unitstr= 'K'

    if varname in logvars:
        logplot = True
        cmap = plt.get_cmap(colmap)
    else:
        logplot = False
        cmap, extopt, cbt, intbounds, allbounds = discrete_cmap_tuner(clim, [minval, maxval], Nlev=Nlev, colmap=colmap,
                                                                      nancol='white')
    cmap.set_bad(str(landgr))


    if logplot:
        #clim=[1e-1,np.amax(v)] #[np.amax(v)*.8,np.amax(v)] #
        if clim[0]<1e0:
            clim[0]=1e0
        pcf=proj.pcolormesh(x,y,v,cmap=cmap,vmin=clim[0], vmax=clim[1],norm=LogNorm(clim[0], clim[1]))
    else:
        #clim=[np.amin(v),np.amax(v)] #[np.amax(v)*.8,np.amax(v)] #
        #pcf=proj.pcolormesh(x,y,v,cmap=plt.get_cmap(colmap),vmin=clim[0], vmax=clim[1])
        pcf = proj.pcolormesh(x, y, v.squeeze(), cmap = cmap, vmin=intbounds[0], vmax=intbounds[-1])

    if plottopo:
        proj.contour(x,y,H,colors=str(landgr),linewidths=0.5,levels=[5.0,10.,20.,30.,40.,50.,60.])

    #cb=plt.colorbar(pcf)

    plt.title(titlestr,size=9)

    #setup specific features
    if setup in ['NSBS','NSBSfull','SNSfull','SNS','GBight','WadSea']:
        #coastlines, etc
        #proj.drawcoastlines(color=(0.3,0.3,0.3),linewidth=0.5) #, resolution='h')
        proj.fillcontinents((landgr,landgr,landgr),lake_color=(landgr,landgr,landgr))
        if showparmer:
            proj.drawparallels(np.arange(51.,57.,1.), labels=[1,0,0,0],fontsize=9)
            #proj.drawmeridians(np.arange(-1.,10.,1.))
            proj.drawmeridians(np.arange(-1.,10.,1.), labels=[0,0,0,1],fontsize=9)
        #ax.patch.set_facecolor((0.9,0.9,1.0))

    #retrieve the axes position to set the colorbar position
    pos1 = ax.get_position()
    poscbar = [pos1.x0 + pos1.width + 0.07, pos1.y0+0.125 , 0.02, pos1.height * 0.62]
    cax = plt.axes(position=poscbar)

    #cb = plt.colorbar(pcf, cax=cax,ticks=cbt)
    # cb.ax.tick_params(axis='both', which='major', labelsize=10, length=0)

    #if showparmer:
    #    cb=plt.colorbar(pcf,shrink=0.6)
    #else:
    #    cb=plt.colorbar(pcf,shrink=0.9)

    #norm = mpl.colors.BoundaryNorm(allbounds[1:-1], cmap.N)
    if logplot:
        cb = plt.colorbar(pcf, cax=cax, format=mpl.ticker.LogFormatter(10, labelOnlyBase=False))
        ticks = list(cb.formatter.locs[:])
        ticks.append(clim[1])
        ticks=[1,5,10,20]
        cb.set_ticks(ticks, update_ticks=True)
        scalesuf = '_CL_%d_%d_log'%(clim[0],clim[1])
    else:
        norm = mpl.colors.BoundaryNorm(intbounds, cmap.N)
        cb = mpl.colorbar.ColorbarBase(cax, cmap=cmap,
                                        norm=norm,
                                        boundaries=allbounds,
                                        extend=extopt,
                                        ticks=cbt,
                                        #ticks= cbt[0::2],
                                        #ticks= cbt[1:-1:2],
                                        spacing='proportional')
        cb.formatter.set_scientific(True)
        cb.formatter.set_powerlimits((-4, 4))
        cb.ax.tick_params(labelsize=10)
        cb.update_ticks()
        scalesuf = '_CL_%.1f_%.1f'%(clim[0],clim[1])
        scalesuf=scalesuf.replace('.0', '')

    #unitstr
    unitstr=unitstr.replace('**', '^')
    unitstr=unitstr.replace('mmol-', 'mmol')
    unitstr = unitstr.replace('mg-', 'mg')
    #print unitstr
    if unitstr != '':
        cb.ax.set_title(unitstr, size=10.)
    else:
        cb.ax.set_title('[-]', size=10.)

    #cint=(clim[1]-clim[0])/5.
    #cb.locator = MultipleLocator(cint)
    #cb.update_ticks()
    #plt.setp(cb.ax.get_yticklabels()[::2], visible=False)

    return scalesuf

if __name__=='__main__':
    varlib = {
        'kd490':[0,2],
        'MLD': [0, 50],
        'MLD/bathymetry': [0., 1.],
        'u': [-.05, .05],
        'v': [-.05, .05],
        'umean': [-.05, .05],
        'vmean': [-.05, .05],
        'currs': [0.0, 0.2],
        't2': [0, 20],
        'elev': [-1.0, 1.0],
        'temp': [12, 18.],
        'salt': [15, 33],
        'sigma_t': [-1.0, 0.0],
        'sigma0': [-1.0, 0.0],
        'EH_abioP_DIN': [0, 48],
        'EH_abioP_DINO3': [0, 50],
        'EH_abioP_DINH4': [0, 5],
        'EH_abioP_DIP': [0, 3.0],
        'EH_abioP_DISi': [0, 45.0],
        'zoo2phy':[0.4,1.0],
        'fr_diatC': [0.5, 0.9],
        'fr_meszooC': [0.1, 0.5],
        #'EH_abioP_DIP': [0, 1.5],
        #'EH_abioS_sed_nn2': [0, 400],  # yintegral
        #'EH_abioS_sed_nn2':[0,3.0], #yaverage
        'EH_abioP_O2_percSat':[50,100],
        'EH_abioP_o2o_pel':[0,50],
        'EH_abioS_o2o_brm':[0,50],
        'total_chlorophyll_calculator_result':[0,10],
        'total_NPPR_calculator_result':[0, 1000], # 1000 (mgC/m2/d) 300:Yintegral (gC/m2/y)
        'GPM_diat_C': [0., 20],
        'GPM_nf_C': [0., 20],
        'GPM_phy_C': [0., 20],
        'GPM_phy_Chl':[0.,10],
        'GPM_phy_NPPR': [0, 1000],  # 1000 (mgC/m2/d) 300:Yintegral (gC/m2/y)
        'GPM_meszoo_C':[0,5.],
        'GPM_miczoo_C': [0, 10.],
        'GPM_zoo_C': [0, 5.],
        'total_nitrogen_calculator_result': [0,40.],
        'total_phosphorus_calculator_result': [0.,1.6],
    }
    if len(sys.argv)>1:
        fname=sys.argv[1]
    else:
        fname = '/home/onur/WORK/projects/2013/gpmeh/sns144-GPMEH-G191007-Fnew3-PPZZSi-PinR-P191010-vS/extract_PZCtotVA_sns144-GPMEH-G191007-Fnew3-PPZZSi-PinR-P191010-vS.2013-mm.nc'
        #fname = '/home/onur/WORK/projects/2013/gpmeh/sns144-GPMEH-PPZZ-P190628-fSG97dChl/extract_RavgC_sns144-GPMEH-PPZZ-P190628-fSG97dChl.2012-mm.nc'
        #fname = '/home/onur/WORK/projects/2013/gpmeh/sns144-GPMEH-P190607-fSG97dChl/extract_skillCS_sns144-GPMEH-P190607-fSG97dChl.2012-mm.nc'
    if len(sys.argv)>2:
        varnames=sys.argv[2].split(',')
    else:
        #varnames = ['total_chlorophyll_calculator_result'] #, 'GPM_diat_C', 'GPM_nf_C', 'GPM_miczoo_C', 'GPM_meszoo_C']
        varnames=['zoo2phy','fr_diatC','fr_meszooC']
        #varnames = ['EH_abioP_O2_percSa']
        #varnames = ['GPM_phy_Chl','GPM_phy_C','GPM_zoo_C']
        #varnames=['EH_abioP_DIN','EH_abioP_DIP']
        #varnames=['EH_abioP_DINO3','EH_abioP_DINH4']
        #varnames=['temp']
        #varnames = ['sigma_t']
        #varnames=['currs']
        #varnames=['kd490']
        #varnames=['fsalt']
        #varnames = ['MLD/bathymetry']
        #varnames=['EH_abioP_o2o_bac', 'EH_abioP_o2o_n4n']
        #varnames=['EH_abioS_sed_nn2', 'EH_abioS_o2o_brm']
        # varnames=['Dissolved_Inorganic_Nitrogen_DIN_nutN_in_water','Dissolved_Inorganic_Phosphorus_DIP_nutP_in_water']
        # varnames =['Chl_chl_in_water']#,'practical_salinity_in_water']
        #varnames=['practical_salinity_in_water']
        #varnames=['temperature_in_water']

    if len(sys.argv)>3:
        TempMeth=sys.argv[3]
    else:
        #TempMeth='snapshots' #'Saverage','Maverage','Yintegral','snapshots'
        #TempMeth='Yaverage_1-2-3-10-11-12'
        #TempMeth = 'Yaverage_3-4-5-6-7-8-9'
        #TempMeth = 'Yaverage_1-2-3'
        #TempMeth = 'Yaverage_6-7-8'
        #TempMeth = 'Yaverage_3-4-5'
        #TempMeth='Yaverage_1-2-12'
        TempMeth='Yaverage_7'
        #TempMeth='Maverage'
        #TempMeth = 'Y2013-M7'

    if len(sys.argv) > 4:
        mode = sys.argv[4]
    else:
        mode = 'singlepanel' # '3panelsperrow'  # 'singlepanel','halfpagewidth'

    if len(sys.argv) > 5:
        VertMeth = sys.argv[5]
    else:
        VertMeth = 'surf' #'surf','avg','int','each', 'SB'
        # setup='deep_lake'

    if len(sys.argv) > 6:
        setup = sys.argv[6]
    else:
        setup = 'SNSfull'  # WadSea SNSfull GBight
        # setup='deep_lake'

    if len(sys.argv) > 7:
        colmap = sys.argv[7]
    else:
        colmap = 'viridis' #'v' 'viridis' #'jet'
        #colmap='Blues_r' #'bwr' #'Reds' #Blues_r' Reds_r

    if len(sys.argv) > 8:
        Nlev = int(sys.argv[8])
    else:
        Nlev= 6 #7 is std

    #cbar lims
    vars = {}
    if len(sys.argv) > 9:
        climstr = sys.argv[9].split(',')
        clim=[float(str) for str in climstr]
        #print('clim: ' + climstr)
        if clim[0]!=clim[1]:
            print ('using the requested colorbar lims: [%s-%s]:'%(clim[0],clim[1]))
            climfromlib=False
            for varn in varnames:
                vars[varn] = clim
        else:
            print ('keeping predefined colorbar lims')
            climfromlib=True
    else:
        print ('keeping predefined colorbar lims')
        climfromlib = True

    if climfromlib:
        for varn in varnames:
            if varn in varlib:
                vars[varn] = varlib[varn]
                #vars[varn]=[-75,5]
            else:
                vars[varn]=[0,0]
                print('clims for ' +varn + ' were not proviveded, natural limits will be used')

    do_2Dplotmap(fname,vars,setup,VertMeth,TempMeth,colmap,Nlev,mode,plottopo=True,datasource='GF')
