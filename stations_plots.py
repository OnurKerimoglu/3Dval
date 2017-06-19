# -*- coding: utf-8 -*-
"""
Created on Mon Jun  12 12:10 2017

@author: kerimoglu.o@gmail.com

"""
import os
import numpy as np
import matplotlib.pyplot as plt
from general_funcs import getproj,format_date_axis

class Style:
    def __init__(self,opt='default'):
        if opt=='TSdefault':
            self.res = 150
            self.figwh=[0, 0]
            self.col={'obs':'0.3','sim':['r','b','g','k']}
            self.line={'obs':'None','sim':['-','-','-','-']}
            self.marker={'obs':'o','sim':['None','None','None','None']}
            self.lw={'obs':1,'sim':[1,1,1,1]}

def stations_plots(plotopts,obs,sim,plotrootpath,stations,timeint,depthints):
    fnamecode= '_%s-%s' %(timeint[0].year, timeint[1].year)

    if plotopts['TS']==True:
        plotpath=os.path.join(plotrootpath,'3Dval_stations'+fnamecode)
        stations_plots_ts(plotopts, obs, sim, plotpath, stations, timeint, depthints, fnamecode)


def stations_plots_ts(plotopts,obs,simset,plotpath,stations,timeint,depthints,fnamecode):
    print('Doing the time series plots for stations:')

    #variables to plot, definitions
    varlongnames={'temp':'Temperature', 'salt':'Salinity'}
    varunits={'temp':'$^\circ$C', 'salt':'psal'}

    #figure parameters:
    colnum= len(depthints.keys())
    rownum=len(plotopts['varns'])
    S = Style(opt=plotopts['TSstyle'])

    #if the plotpath doesn't exist, create it
    if not os.path.exists(plotpath):
        os.makedirs(plotpath)

    #projection (for showing stations on maps)
    proj=getproj(setup='SNSfull',projpath=os.path.dirname(os.path.realpath(__file__)))

    #extract the station list if not provided
    if len(stations)==0:
        stations=list(obs.keys())

    for stationno,station in enumerate(stations):
        print ('  '+station)

        # genereate a new figure, size of which is a function of number of variables (rows) and layers (columns) to be show
        fig = prepfig(S.res, S.figwh, colnum, rownum)
        fig.subplots_adjust(hspace=.2, wspace=.15, left=0.15, right=0.85, top=0.7, bottom=0.03)
        #name of the station
        fig.text(0.4,0.98,station+'\n$Z_{max}$=%s'%obs[station]['bottom_depth'],verticalalignment='top',horizontalalignment='left',size=10)
        # show the location of the station on a map in one panel
        ax = plt.axes([0.15, 0.75, 0.23, 0.23])
        markstatonmap(ax, proj, station, obs[station]['lon'], obs[station]['lat'], obs[station]['bottom_depth'])

        for layerno,layer in enumerate(depthints.keys()):

            # if no plot is made, don't save an empty figure, so track whether any plot is made in the figure
            anyplotinfig = False

            for varno,varname in enumerate(plotopts['varns']): #in each panel

                #create a panel
                #ax = plt.subplot(rownum, colnum, (varno+1)*colnum + varno+1)
                ax=plt.subplot2grid((rownum,colnum),(varno,layerno))
                #if no series are available in the panel, legending and etc will be problamatic, so track it
                anyplotinax=False

                #if first panel, indicate layer as the title
                if varno==0: #
                    plt.title(layer,size=9)

                #make the plots:
                hset = []; idset = []  #list of handles (needed for legend)
                #plot obs
                if obs[station][varname]['presence']:
                    hset,idset,anyplotinax,anyplotinfig = plot_ts_panel(anyplotinax,anyplotinfig,hset,idset,'obs',ax,obs[station][varname][layer]['time'],obs[station][varname][layer]['value'],timeint,S,'obs')
                # plot each sim
                for simno,simname in enumerate(plotopts['sims2plot']): #enumerate(simset.keys()):
                    if simset[simname][station][varname]['presence']:
                        hset,idset,anyplotinax,anyplotinfig = plot_ts_panel(anyplotinax,anyplotinfig,hset,idset,simname,ax,simset[simname][station][varname][layer]['time'],simset[simname][station][varname][layer]['value'],timeint, S, 'sim',simno)

                #ylabel:varname, unit
                plt.ylabel(varlongnames[varname]+' ['+varunits[varname]+']',size=9)
                ax.tick_params(axis='y', which='major', direction='out', labelsize=9)

                #format date axes
                if varno==len(plotopts['varns'])-1:
                    format_date_axis(ax, timeint)
                else:
                    ax.set_xticklabels([])

                #add legend
                if anyplotinax:
                    #ax = plt.axes([0.6, 0.75, 0.4, 0.15],visible=False) #todo: place the legend in a dedicated axis within the top margin
                    lgd = ax.legend(handles=hset, labels=idset, loc='lower right',fontsize=9, numpoints=1, bbox_to_anchor=(1.1, 0.5))

        #save&close the figure
        if not anyplotinfig:
            plt.close()
        else:
            fname = os.path.join(plotpath,'TSplots_%s_%s.png' % (fnamecode, station))
            fig.savefig(fname,dpi=S.res, bbox_extra_artists=(lgd,), bbox_inches='tight')
            plt.close()
            print ('figure saved:%s'%fname)
            #return
    return

def plot_ts_panel(anyplotinax,anyplotinfig,hset,idset,id,ax,times,values,timeint,S,seriestype,sno=-1):
    tind =np.where((times>=timeint[0]) * (times<=timeint[1]))[0]
    if len(tind)==0:
        return (hset,idset,anyplotinax,anyplotinfig)

    if seriestype=='obs':
        h, = ax.plot(times[tind], values[tind], linestyle=S.line['obs'], marker=S.marker['obs'], lw=S.lw['obs'], color=S.col['obs'], mfc=S.col['obs'], mec=S.col['obs'], markersize=1, label=id)
    else:
        if sno==-1:
            raise(Exception('Simulation # must be provided (sno)'))
        h, = ax.plot(times[tind], values[tind], linestyle=S.line['sim'][sno], marker=S.marker['sim'][sno], lw=S.lw['sim'][sno], color=S.col['sim'][sno], mfc=S.col['sim'][sno], mec=S.col['sim'][sno], markersize=1,label=id)
    hset.append(h)
    idset.append(id)
    anyplotinax = True
    anyplotinfig = True
    # calculate and annotate statistics?

    return (hset, idset, anyplotinax, anyplotinfig)

def markstatonmap(ax, proj, station, lon,lat,maxz):
    tx, ty = proj(lon, lat)
    proj.plot(tx, ty, 'k.', markersize=3, marker='d')
    #plt.text(tx, ty + 8000, ' ($z_{max}$=%s)'%maxz, size=10.0, horizontalalignment='center', verticalalignment='bottom', color='black',backgroundcolor='white')
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    # proj.drawlsmask(land_color='.8', ocean_color='w',resolution='f')
    proj.drawcoastlines(color=(0.3, 0.3, 0.3), linewidth=0.5)
    proj.fillcontinents((.8, .8, .8), lake_color=(0.6, 0.6, 1.0))

def prepfig(res,figwh,colno,rowno):
    # start a figure
    if (figwh[0] == 0):
        figwh = cm2inch(12 * colno + 2, 4 * rowno)
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