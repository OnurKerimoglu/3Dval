# -*- coding: utf-8 -*-
"""
Created on Mon Jun  12 12:10 2017

@author: kerimoglu.o@gmail.com

"""
import os
import numpy as np
import numpy.ma as ma
import cftime
import matplotlib as mpl
import matplotlib.pyplot as plt
from general_funcs import getproj,format_date_axis

class Style:
    def __init__(self,opt='default'):
        if opt=='TSdefault':
            self.res = 300
            #self.figwh=[12, 10] #BG2020
            self.figwh=[12, 15] #BG2020
            #self.col={'obs':'0.6','sim':['k','r','b','g']}
            self.col={'obs':'0.6','sim':['k','tomato','darkblue','deepskyblue']}
            self.line={'obs':'None','sim':['-','-','-','-']}
            self.marker={'obs':'o','sim':['None','None','None','None']}
            self.lw={'obs':1,'sim':[1,2,1,1]}

def stations_plots(plotopts,obs,sim,plotrootpath,statsets,stations,timeint,depthints):
    fnamecode= '_%s-%s' %(timeint[0].year, timeint[1].year)

    if plotopts['TS']==True:
        for statset in statsets:
            plotpath=os.path.join(plotrootpath,'3Dval_stations'+fnamecode,statset)
            stations_plots_ts(plotopts, obs, sim, plotpath, stations, timeint, depthints, fnamecode)

def stations_plots_ts(plotopts,obs,simset,plotpath,stations,timeint,depthints,fnamecode):
    print('Doing the time series plots for stations:')

    #variables to plot, definitions
    varlongnames={'temp':'Temperature', 'salt':'Salinity', 'DOs':'O2 sat.', 'DIN':'DIN', 'DIP':'DIP', 'Chl':'Chl'}
    varunits={'temp':u'\N{DEGREE SIGN}C', 'salt':'g/kg', 'DOs':'%', 'DIN':'$\mu$MN', 'NH4':'$\mu$MN', 'NO3':'$\mu$MN', 'Si':'$\mu$MSi', 'DIP':'$\mu$MP', 'Chl':'mg/m$^3$'}
    varlims_offshore={'temp':[0,20],'salt':[28,35],'NH4':[0,20],'NO3':[0,60],'DIN':[0,50],
                      'DIP':[0,2.1], 'Si':[0,50],'Chl':[0,20]}
    varticks_offshore={'temp':[0,5,10,15,20],'salt':[29,31,33,35],
                       'NH4': [0,5,10,15,20],'NO3': [0,10,20,30,40,50],'DIN': [0,10,20,30,40,50],
                       'DIP': [0,0.5,1.0,1.5,2.0],'Si': [0,10,20,30,40,50], 'Chl': [0,5,10,15,20]}
    varlims_coastal = {'temp': [-1.0, 22.], 'salt': [0, 30],'NH4':[0,50],'NO3':[0,350], 'DIN': [0, 350],
                       'DIP': [0, 3.5], 'Si':[0, 250], 'Chl': [0, 40]}
    varticks_coastal = {'temp': [0, 5, 10, 15, 20], 'salt': [0,10,20,30],
                        'NH4': [0,10,20,30,40,50],'NO3':[0,100,200,300],'DIN':[0,100,200,300],
                        'DIP':[0,1,2,3],'Si': [0,50,100,150,200,250],'Chl':[0,10,20,30,40]}
    axtune=True
    #figure parameters:
    #colnum= len(depthints.keys())
    colnum=1
    rownum=max(2,len(plotopts['varns']))
    S = Style(opt=plotopts['TSstyle'])

    #if the plotpath doesn't exist, create it
    if not os.path.exists(plotpath):
        os.makedirs(plotpath)

    #projection (for showing stations on maps)
    proj=getproj(setup='SNSfull',projpath=os.path.dirname(os.path.realpath(__file__)))
    #proj = getproj(setup='WadSea', projpath=os.path.dirname(os.path.realpath(__file__)))

    #extract the station list if not provided
    if len(stations)==0:
        stations=list(obs.keys())

    for stationno,station in enumerate(stations):
        print ('  '+station)
        if station in ['Cuxhaven','HPA-Elbe', 'Norderelbe']:
            varticks = varticks_coastal
            varlims = varlims_coastal
        else:
            varticks = varticks_offshore
            varlims = varlims_offshore

        for layerno,layer in enumerate(depthints.keys()):

            # genereate a new figure, size of which is a function of number of variables (rows) and layers (columns) to be show
            fig = prepfig(S.res, S.figwh, colnum, rownum, timeint)
            fig.subplots_adjust(hspace=.30, wspace=.2, left=0.15, right=0.7, top=0.75, bottom=0.06)
            # name of the station
            fig.text(0.76, 0.99, station + '\n$Z_{max}$=%.1f m' % obs[station]['bottom_depth'], verticalalignment='top',
                     horizontalalignment='right', size=10)
            # show the location of the station on a map in one panel
            ax = plt.axes([0.77, 0.77, 0.22, 0.22])
            markstatonmap(ax, proj, station, obs[station]['lon'], obs[station]['lat'], obs[station]['bottom_depth'])

            # if no plot is made, don't save an empty figure: to achieve this track whether any plot is made
            anyplotinfig = False
            
            for varno,varname in enumerate(plotopts['varns']): #in each panel

                #create a panel
                #ax=plt.subplot2grid((rownum,colnum),(varno,layerno))
                ax = plt.subplot2grid((rownum, colnum), (varno, 0))
                #if no series are available in the panel, legending and etc will be problamatic, so track it
                anyplotinax=False

                #if first panel, indicate layer as the title
                if varno==0: #
                    plt.title(station + ', ' + layer + '\n',size=9)

                #make the plots:
                hset = []; idset = []  #list of handles (needed for legend)

                #plot obs
                if obs[station][varname]['presence']:
                    # limit the months to include
                    months2keep = [] #[7, 8]
                    if len(months2keep)>0:
                        monsuf='_M'+'-'.join(map(str,months2keep))
                        obs[station][varname][layer] = stationdata_filter_time(obs[station][varname][layer],months2keep)
                    else:
                        monsuf=''
                    hset,idset,anyplotinax,anyplotinfig = plot_ts_panel(anyplotinax,anyplotinfig,hset,idset,'obs',ax,
                                                                        obs[station][varname][layer]['time'],
                                                                        obs[station][varname][layer]['value'],
                                                                        timeint,S,'obs')
                # plot each sim
                for simno,simname in enumerate(plotopts['sims2plot']): #enumerate(simset.keys()):
                    if simset[simname][station][varname]['presence']:
                        if 'GF-' in simname:
                            simnameleg=simname.replace('GF','sim')
                        else:
                            simnameleg=simname
                        hset,idset,anyplotinax,anyplotinfig = plot_ts_panel(anyplotinax,anyplotinfig,hset,idset,simnameleg,ax,
                                                                            simset[simname][station][varname][layer]['time'],
                                                                            simset[simname][station][varname][layer]['value'],
                                                                            timeint, S, 'sim',simno)
                        # annotate skill scores
                        if (obs[station][varname]['presence']) and (simset[simname][station][varname]['presence']):
                            skills=get_skillscores(obs[station][varname][layer],simset[simname][station][varname][layer],timeint)
                            if (simno==0) and (skills['n'] != 0):
                                if (len(plotopts['sims2plot']) - 0) == 1:
                                    y = 1.05
                                else:
                                    y = 1.08 - 0.12 * simno
                                plt.text(1.0, y, r'$B^*$:%3.2f, $\rho$:%3.2f, $n$:%d'
                                        %(np.round(skills['B*'] * 100) / 100, np.round(skills['r'] * 100) / 100, skills['n']),
                                        fontsize=9, ha='right',va='center', transform=ax.transAxes, color=S.col['sim'][simno])

                #ylabel:varname, unit
                if varname in varlongnames.keys():
                    varlongname=varlongnames[varname]
                else:
                    varlongname=varname
                plt.ylabel(varlongname+' ['+varunits[varname]+']',size=9)
                ax.get_yaxis().set_label_coords(-0.17, 0.5)

                if (axtune) and (varname in varticks.keys()):
                    yticks = varticks[varname]
                    ylims = varlims[varname]
                    ax.set_ylim([ylims[0],ylims[-1]])
                    ax.set_yticks(yticks)
                    #ax.yaxis.set_major_locator(yticks)
                    if (yticks[-1] - yticks[0]) <= 3.0:
                        ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.1))
                    elif (yticks[-1]-yticks[0])<=36:
                        ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(1.00))
                    elif (yticks[-1]-yticks[0])<=100:
                        ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(5.00))
                    elif (yticks[-1]-yticks[0])<=500:
                        ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(25.00))
                ax.tick_params(axis='y', which='minor', direction='in', labelsize=9)
                ax.tick_params(axis='y', which='major', direction='out', labelsize=9)
                ax.grid(axis='y', which='major', color='0.5', linestyle='-', linewidth=.5)

                #format date axes# convert to 1-d, assuming that all measurements are from surface
                format_date_axis(ax, timeint)
                if not varno==len(plotopts['varns'])-1:
                    ax.set_xticklabels([])

                #add legend
                if anyplotinax:
                    #ax = plt.axes([0.6, 0.75, 0.4, 0.15],visible=False) #todo: place the legend in a dedicated axis within the top margin
                    lgd = ax.legend(handles=hset, labels=idset, loc='center left',fontsize=9, numpoints=1, bbox_to_anchor=(1.02, 0.5))

            #save&close the figure
            if not anyplotinfig:
                plt.close()
            else:
                fname = os.path.join(plotpath,'TSplots%s_%s_%s%s.png' % (fnamecode, station, layer,monsuf))
                fig.savefig(fname,dpi=S.res, bbox_extra_artists=(lgd,)) #, bbox_inches='tight')
                plt.close()
                print ('figure saved:%s'%fname)
                #return
    return

def stationdata_filter_time(din,months2keep):

    #find the indices
    dates=din['time']
    months=np.array([date.month for date in dates])
    #tind=np.array([])
    tfilter=False #reject all
    for m in months2keep:
        tfilter=tfilter + (months==m)
        #tind= np.concatenate((tind,np.where(months==m)[0]),axis=0)
    tind=np.where(tfilter)[0]
    #apply the indices to all fields
    dout={}
    for f in din.keys():
        if len(din[f])==len(dates):
            dout[f]=din[f][tind]
        else:
            dout[f]=din[f]

    return dout

def get_skillscores(obs,sim,timeint):
    #reduce time
    tind = np.where((obs['time'] >= timeint[0]) * (obs['time'] <= timeint[1]))[0]
    t,o,s=match_time(obs['time'][tind], obs['value'][tind], sim['time'],sim['value'])

    # calculate statistics:
    skills={}
    skills['n']=len(o)
    from scipy.stats import pearsonr

    if len(tind)>1:
        r, p = pearsonr(s, o)
    else:
        r=np.nan
        p=np.nan

    if p < 0.001:
        pstr = '***'
    elif p < 0.01:
        pstr = '**'
    elif p < 0.05:
        pstr = '*'
    else:
        pstr = ''
    skills['r']=r
    skills['p']=p
    skills['pstr']=pstr
    # B:
    skills['B'] = s.mean() - o.mean()
    skills['B*'] = (s.mean() - o.mean()) / o.mean()
    # rmsd:
    skills['rmsd']=np.sqrt((o-s)**2)
    # MEF
    oe=(o-o.mean())**2
    se=(o-s)**2
    skills['mef']=1-se.sum()/oe.sum()

    return skills

def match_time(t1dt,v1,t2dt,v2):
    # find the common dates, at a daily resolution
    #convert the datetimes to dates
    t1d=np.array([d.date() for d in t1dt])
    t2d=np.array([d.date() for d in t2dt])
    comd = np.sort(list(set(t1d).intersection(t2d)))
    # find the indices
    lt1=list(t1d)
    lt2=list(t2d)
    i1 = [lt1.index(d) for d in comd]
    i2 = [lt2.index(d) for d in comd]
    if not all(t1d[i1] == t2d[i2]):
        raise(Exception('Error encountered while temporal-pairing the obs and sim : obs(t)!= sim(t)'))
    return(t1dt[i1],v1[i1],v2[i2])

def plot_ts_panel(anyplotinax,anyplotinfig,hset,idset,id,ax,times,values,timeint,S,seriestype,sno=-1):
    tind =np.where((times>=timeint[0]) * (times<=timeint[1]))[0]
    if len(tind)==0:
        return (hset,idset,anyplotinax,anyplotinfig)

    if seriestype=='obs':
        h, = ax.plot(times[tind], values[tind], linestyle=S.line['obs'], marker=S.marker['obs'], lw=S.lw['obs'], color=S.col['obs'], mfc=S.col['obs'], mec=S.col['obs'], markersize=3, label=id)
    else:
        if sno==-1:
            raise(Exception('Simulation # must be provided (sno)'))
        # h, = ax.plot(cftime.num2pydate(times[tind],units=days), values[tind], linestyle=S.line['sim'][sno], marker=S.marker['sim'][sno], lw=S.lw['sim'][sno], color=S.col['sim'][sno], mfc=S.col['sim'][sno], mec=S.col['sim'][sno], markersize=1,label=id)
        h, = ax.plot(times[tind], values[tind], linestyle=S.line['sim'][sno], marker=S.marker['sim'][sno], lw=S.lw['sim'][sno], color=S.col['sim'][sno], mfc=S.col['sim'][sno], mec=S.col['sim'][sno], markersize=1,label=id)
    hset.append(h)
    idset.append(id)
    anyplotinax = True
    anyplotinfig = True

    return (hset, idset, anyplotinax, anyplotinfig)

def markstatonmap(ax, proj, station, lon,lat,maxz):
    tx, ty = proj(lon, lat)
    proj.plot(tx, ty, 'r.', markersize=5, marker='d')
    #plt.text(tx, ty + 8000, ' ($z_{max}$=%s)'%maxz, size=10.0, horizontalalignment='center', verticalalignment='bottom', color='black',backgroundcolor='white')
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    # proj.drawlsmask(land_color='.8', ocean_color='w',resolution='f')
    proj.drawcoastlines(color=(0.3, 0.3, 0.3), linewidth=0.5)
    proj.fillcontinents((.8, .8, .8), lake_color=(0.6, 0.6, 1.0))

def prepfig(res,figwh,colno,rowno,timeint):
    # start a figure
    #years:
    numy=timeint[1].year-timeint[0].year+1
    if (figwh[0] == 0):
        #figwh = cm2inch(3.75 * numy * colno + 2, 4 * rowno)
        figwh = cm2inch(2.25*numy * colno + 5, 4 * rowno)
    else:
        figwh = cm2inch(figwh[0], figwh[1])
    fig = plt.figure(figsize=figwh, dpi=res)
    return fig

def cm2inch(*tupl):
    inch = 2.54
    if type(tupl[0]) == tuple:
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)
