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
            self.col={'obs':'0.3','sim':['b','r','g','k']}
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
    rownum=len(plotopts['varns'])+1 #+1 is for the map with station marked on it
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
        fig.text(0.5,0.96,station,horizontalalignment='center',size=10)

        for layerno,layer in enumerate(depthints.keys()):

            anyplotinfig = False  # if no plot is made, don't save an empty figure?

            #show the location of the station on a map in one panel
            ax = plt.subplot2grid((rownum,colnum),(0,0),colspan=colnum)
            markstatonmap(ax, proj, station, obs[station]['lon'],obs[station]['lat'],obs[station]['bottom_depth'])

            for varno,varname in enumerate(plotopts['varns']): #in each panel

               #create a panel
                #ax = plt.subplot(rownum, colnum, (varno+1)*colnum + varno+1)
                ax=plt.subplot2grid((rownum,colnum),(varno+1,layerno))
                anyplotinax=False

                #if first panel, indicate layer as the title
                if varno==0: #
                    plt.title(layer,size=9)

                #make the plots:
                hset = []; idset = []  #list of handles (needed for legend)
                #plot obs
                if obs[station][varname]['presence']:
                    hset,idset,plotdone=plot_ts_panel(hset,idset,'obs',ax,obs[station][varname][layer]['time'],obs[station][varname][layer]['value'],timeint,S,'obs')
                    if plotdone:
                        anyplotinax = True; anyplotinfig = True
                # plot each sim
                for simno,simname in enumerate(simset.keys()):
                    if simset[simname][station][varname]['presence']:
                        hset,idset,plotdone= plot_ts_panel(hset,idset,simname,ax,simset[sim][station][varname][layer]['time'], obs[station][varname][layer]['value'],timeint, S, 'sims',simno)
                        if plotdone:
                            anyplotinax = True; anyplotinfig = True

                #ylabel:varname, unit
                plt.ylabel(varlongnames[varname]+' ['+varunits[varname]+']',size=9)
                ax.tick_params(axis='y', which='major', direction='out', labelsize=9)

                #format date axes
                format_date_axis(ax, timeint)

                #add legend
                if anyplotinax:
                    lgd = plt.legend(handles=hset, labels=idset, loc='center right',fontsize=9, numpoints=1, bbox_to_anchor=(1.1, 0.5))

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

def plot_ts_panel(hset,idset,id,ax,times,values,timeint,S,seriestype,sno=-1):
    tind =np.where((times>=timeint[0]) * (times<=timeint[1]))[0]
    if len(tind)>0:
        if seriestype=='obs':
            h, = ax.plot(times[tind], values[tind], linestyle=S.line['obs'], marker=S.marker['obs'], lw=S.lw['obs'], color=S.col['obs'], mfc=S.col['obs'], mec=S.col['obs'], markersize=1, label=id)
        else:
            if sno==-1:
                raise(Exception('Simulation # must be provided (sno)'))
            h, = ax.plot(times[tind], values[tind], linestyle=S.line['sim'][sno], marker=S.marker['sim'][sno], lw=S.lw['sim'][sno], color=S.col['sim'][sno], mfc=S.col['sim'][sno], mec=S.col['sim'][sno], markersize=1,label=id)
        hset.append(h)
        idset.append(id)
        plotdone=True
        # calculate and annotate statistics?
    else:
        plotdone=False

    return (hset,idset,plotdone)

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
    fig.subplots_adjust(hspace=.5, wspace=.5, left=0.15, right=0.85)  # top=0.05,bottom=0.04,
    # plt.figtext(0.5, 0.97, gentit, ha='center', va='top')
    return fig

def cm2inch(*tupl):
    inch = 2.54
    if type(tupl[0]) == tuple:
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)