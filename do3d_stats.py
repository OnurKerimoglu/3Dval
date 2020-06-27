"""
Created on 23 Nov 2016
@author: onur.kerimoglu@hzg.de
"""

import os,sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

#local modules
#sys.path.insert(1, '/home/onur/WORK/codes/python/')
from matchup_refmod import get_matchups,filter_matchups_vertloc,combine_matchedup_sets
from taylorDiagram import TaylorDiagram

def main(files,ftypes,fnames,varnames,plfpath,vertloc,yrs,seasons,split,scatter,taylor,barplotsacross,remOLmeth,demo=False):
    if demo==True:
        # get statistics
        [stdrefs,samples,x95,y95,x99,y99]=get_stats_demo()
        # color and panel arrangement
        # Colormap (see http://www.scipy.org/Cookbook/Matplotlib/Show_colormaps)
        colors = plt.matplotlib.cm.Set1(np.linspace(0, 1, len(samples['winter'])))
        rects = dict(winter=221,spring=222,summer=223,autumn=224)
        #make the Taylor plot
        fname='taylor_demo.png'
        plot_Taylor(fname,rects,colors,stdrefs,samples)
        return

    #collect statistics & do scatter plots
    print ('Calculating Statistics and Making Scatter Plots')
    print ('plotdir:' + plfpath)

    #Scatter and Taylor diagrams for each model
    MSVstats={};MSVset={}
    for mno,m in enumerate(fnames[1:]):
        print ('  Model:'+m)
        Vset,Uset=get_matchups(files,fnames,fnames,varnames,yrs)

        if barplotsacross == 'salinity':
            filesS=['/home/onur/WORK/projects/GB/data/ices/lon-1-10_lat50-57/raw/2010-2014/TS_data_block.pickle',
                    '/home/onur/WORK/projects/2013/gpmeh/sns144-GPMEH-G191216-Fnew3-PPZZSi-vS-P191223/extract_MphysC_sns144-GPMEH-G191007-Fnew3-PPZZSi-PinR-P191010-vS.2011-2014_zSB.nc']
            fnamesS=['ICES', 'GETM-GPM.PPZZ-EH.GC']
            VsetS, UsetS = get_matchups(filesS, fnamesS, fnames, ['S'], yrs)
            Vset=combine_matchedup_sets(Vset,VsetS,Uset,UsetS)

        # year interval suffix
        yrsuf = '_' + str(yrs[0]) + '-' + str(yrs[1])

        if split=='none':
            #filter for vertical location
            Vset,vlsuf = filter_matchups_vertloc(Vset,vertloc)

            # group and store into a seasons structure
            SVset=group_seasons(Vset,seasons)
            SVstats={}
            for s in seasons:
                print ('    Season:'+s)
                # remove the outliers.
                Dset,OLsuf = remove_outliers(SVset[s],remOLmeth)
                # get statistics
                SVstats[s] = get_stats(Dset)

                if barplotsacross != 'none':
                    fnameroot = files[mno + 1].split('/')[-1].replace('.nc', '') + '_vs_' + fnames[
                        0] + vlsuf + yrsuf + '_S' + s + OLsuf + '_across_' + barplotsacross
                    barplots_across(Dset, Uset, barplotsacross, plfpath, fnameroot,
                                    [ftypes[0], ftypes[mno + 1]])  # ,title=s
                if scatter:
                    fnameroot = files[mno+1].split('/')[-1].replace('.nc', '') + '_vs_' + fnames[0] + vlsuf +yrsuf +'_S'+s+OLsuf+'_scatter'
                    scatter_refmod(Dset,SVstats[s],Uset,plfpath,fnameroot,[ftypes[0],ftypes[mno+1]]) # ,title=s
                if taylor:
                    fname = os.path.join(plfpath,files[mno + 1].split('/')[-1].replace('.nc', '') + '_vs_' + fnames[0] + vlsuf + yrsuf + '_S' + s+OLsuf+'_taylor'+'_' + '_'.join(varnames)+'.png')
                    stdrefs, samples, rects, colors = arrange4Taylor(SVstats[s], 'V', s=s)
                    plot_Taylor(fname,rects,colors,stdrefs,samples,legend=False)

        elif split=='regions':
            # filter for vertical location
            Vset, vlsuf = filter_matchups_vertloc(Vset, vertloc)

            #regions=['NW','NE','SW','SE']
            regions = ['WadSea']
            #regions = ['WadSea']
            # group and store into a seasonsXregions structure
            SVRset = group_seasonsXregions(Vset, seasons,regions)
            SVRstats = {}
            for s in seasons:
                print ('    Season:' + s)
                VRstats={}
                for v in SVRset[s].keys():
                    print ('      Variable:' + v)
                    #remove the outliers
                    Dset,OLsuf=remove_outliers(SVRset[s][v],remOLmeth)
                    # get statistics
                    VRstats[v] = get_stats(Dset)

                    if taylor:
                        fname = os.path.join(plfpath,files[mno + 1].split('/')[-1].replace('.nc', '') + '_vs_' + fnames[0] + vlsuf + yrsuf + '_S' + s + '_'+v+OLsuf+ '_taylor' + '_' + '_'.join(regions) + '.png')
                        stdrefs, samples, rects, colors = arrange4Taylor(VRstats[v], 'V', s=s)
                        plot_Taylor(fname, rects, colors, stdrefs, samples, legend=False)

                    for r in VRstats[v].keys():
                        if scatter:
                            fnameroot = files[mno+1].split('/')[-1].replace('.nc', '') + '_vs_' + fnames[0] + vlsuf + yrsuf + '_S' + s + '_'+v+OLsuf+'_scatter+'+r
                            scatter_refmod({v:Dset[r]}, {v:VRstats[v][r]}, Uset, plfpath,fnameroot, [fnames[0], fnames[mno + 1]])  # ,title=s

                        if barplotsacross != 'none':
                            fnameroot = files[mno + 1].split('/')[-1].replace('.nc', '') + '_vs_' + fnames[0] + \
                                        vlsuf + yrsuf + '_S' + s + OLsuf + '_across_' + barplotsacross + '_'+r
                            barplots_across({v:Dset[r]}, Uset, barplotsacross, plfpath, fnameroot, [ftypes[0], ftypes[mno + 1]])

                SVRstats[s]=VRstats
        elif split == 'vertloc':
            vertlocs=['surf','bot']
            VZstats={}
            for vertloc in vertlocs:
                VZstats[vertloc]=filter_matchups_vertloc(Vset, Uset, vertloc)
        else:
            raise(ValueError('unknown splitting method: '+split))

                #make a multi-panel Taylor plot, each panel for a season, different variables each panel
        #fname = files[mno + 1].replace('.', '_') + '_vs_' + fnames[0] + '_S' + s + '_taylor' + '_' + '_'.join(varnames) + '.png'

        #stdrefs,samples,rects,colors=arrange4Taylor(SVstats,'SxV')
        #plot_Taylor(fname,rects,colors,stdrefs,samples)
        #print 'plotted:' + fname

        #collect for model-intercomparison
        #MSVstats[m] = SVstats
        #MSVset[m] = SVset

    #Model-intercomparison:
    #Vor each variable, make a multi-panel Taylor plot, each panel for a season, different models for each panel
    #for v in varnames
        #fname=?+'.png'
        #stdrefs,samples,rects,colors=arrange4Taylor(MSVstats,'SxM',v)
        #plot_Taylor(fname,rects,colors,stdrefs,samples)
        #print 'plotted'+fname
    print ('done.')
    return

def remove_outliers(Dsetin, method='percentile',logtransform=True):
    if method=='none':
        suf=''
        return (Dsetin,suf)

    Dsetout={}
    #for each variable
    for v in Dsetin.keys():
        if logtransform and (v in ['Chl']):
            Dsetin[v]['ref'] = np.log(Dsetin[v]['ref'])
            Dsetin[v]['model'] = np.log(Dsetin[v]['model'])
        #find the indices to keep
        if method=='percentile':
            # parameters for the percentile method
            perc = 99
            #find the treshold value based on the percentile
            HthR=np.percentile(Dsetin[v]['ref'],perc)
            HthM = np.percentile(Dsetin[v]['model'], perc)
            #find the valid data for both ref and model
            ind=(Dsetin[v]['ref']<=HthR) * (Dsetin[v]['model']<=HthM)
            nrem=len(ind)-sum(ind)
            methmsg='%dth percentile method'%perc
            # construct a suffix
            suf = '_' + str(perc) + 'perc'
        elif method=='middle':
            Mperc=90
            Lperc = (100. - Mperc) / 2.
            Hperc=100.-(100-Mperc)/2.
            # find the treshold value based on the percentile
            HthR = np.percentile(Dsetin[v]['ref'], Hperc)
            HthM = np.percentile(Dsetin[v]['model'], Hperc)
            LthR = np.percentile(Dsetin[v]['ref'], Lperc)
            LthM = np.percentile(Dsetin[v]['model'], Lperc)
            # find the valid data for both ref and model
            indR = (Dsetin[v]['ref'] <= HthR) * (Dsetin[v]['ref'] >= LthR)
            indM = (Dsetin[v]['model'] <= HthM) * (Dsetin[v]['model'] >= LthM)
            ind=indR*indM
            nrem = len(ind) - sum(ind)
            methmsg = 'middle %d method' % Mperc
            # construct a suffix
            suf = '_' + 'mid'+str(Mperc)
        else:
            raise(ValueError('unknown method:%s for outlier trimming'%method))

        #filter each field with the indices
        Dsetout[v] = {}
        for f in Dsetin[v].keys():
            if list(Dsetin[v][f])==[-1]:
                Dsetout[v][f] = [-1]
            else:
                Dsetout[v][f]=Dsetin[v][f][ind]

        if logtransform and (v in ['Chl']):
            Dsetout[v]['ref'] = np.exp(Dsetout[v]['ref'])
            Dsetout[v]['model'] = np.exp(Dsetout[v]['model'])

    print ('removed %d outliers based on %s'%(nrem, methmsg))

    return (Dsetout,suf)

def barplots_across(Vset,Uset,acrossvar,plfpath,fnameroot,realids,title=''):
    possibleenvvars=['botdepth','salinity']
    if acrossvar not in possibleenvvars:
        raise (Exception('barplots_across: unknown environmental parameter "%s". Available options are: %s.'%(acrossvar,', '.join(possibleenvvars))))

    print('barplots_across: %s:'%acrossvar)
    import pandas as pd
    import numpy as np
    import seaborn as sns

    plottype = 'box'  # box,violin

    # pal = {'2011': 'skyblue', '2012': 'darkblue', '2013': 'orange', '2014': 'red'}
    palw = {'obs': 'white', 'sim': 'gray'}
    vars2plot=Vset.keys() #i.e., everything available
    #vars2plot=['NO3','DIP','Chl']
    unitspretty = {'NO3': '$\mu$MN','NH4': '$\mu$MN', 'DIP': '$\mu$MP', 'Chl':'mgChl/m$^3$'}  # 'NO3:DIP': 'molN/molP',

    ids=['obs','sim']
    if acrossvar=='botdepth':
        nodes = ['<10', '10-20', '>20']
        xlabel = 'Bottom depth [m]'
        envmethod=''
    elif acrossvar=='salinity':
        nodes = ['<29.5', '29.5-32', '>32']
        xlabel = 'Salinity [g/kg]'
        envmethod='abs' #'rel','abs'
    else:
         raise(Exception('unknown environmental variable: %s'%acrossvar))

    for vno,v in enumerate(vars2plot):
        print('      %s'%v,end=": ")
        f = plt.figure(figsize=(3.0, 3.0), dpi=150)
        ax = f.add_axes([0.2, 0.2, 0.6, 0.6])

        # collect the data
        df0=[[],[],[]]
        minval=1e9
        maxval=0.0
        nvec=[None]*len(nodes)
        for nodeno, node in enumerate(nodes):
            for idno,id in enumerate(ids):
                envvar = 'depth'
                if id=='obs':
                    vals=Vset[v]['ref']
                elif id=='sim':
                    vals=Vset[v]['model']
                if acrossvar=='botdepth':
                    depths=Vset[v]['maxdepths']
                    if '-' in node:
                        envind=np.where((depths>=float(node.split('-')[0])) * (depths<float(node.split('-')[1])))
                    elif '<' in node:
                        envind = np.where(depths < float(node.split('<')[1]))
                    elif '>' in node:
                        envind = np.where(depths >= float(node.split('>')[1]))
                elif acrossvar == 'salinity':
                    envvar = 'salinity'
                    if envmethod == 'abs':
                        S = Vset[v]['S-ref']
                    elif envmethod =='rel':
                        if id=='obs':
                            S = Vset[v]['S-ref']
                        elif id=='sim':
                            S = Vset[v]['S-mod']
                    if '-' in node:
                        envind = np.where((S >= float(node.split('-')[0])) * (S < float(node.split('-')[1])))
                    elif '<' in node:
                        envind = np.where(S < float(node.split('<')[1]))
                    elif '>' in node:
                        envind = np.where(S >= float(node.split('>')[1]))

                vals=vals[envind]
                nvec[nodeno] = len(vals)

                minval = min(min(vals), minval)
                maxval = max(max(vals), maxval)
                df_element = [vals, [node] * len(vals), [id] * len(vals)]
                df0=np.concatenate([df0, df_element], axis=1)

        DF = pd.DataFrame(columns=['value', envvar, 'dataset'], data=df0.T)
        DF['value'] = DF['value'].astype(float)

        # do the plot
        if plottype == 'violin':
            sns.violinplot(x=envvar, y='value', data=DF, inner=None, hue='dataset', palette=palw)  # , color = ".8")
        elif plottype == 'box':
            sns.boxplot(x=envvar, y='value', data=DF, ax=ax, hue='dataset', palette=palw)
        #sns.stripplot(x='station', y='value', data=Odf, ax=ax, color='white', hue='year', palette=palw,
                      # palette= pal2,
        #              jitter=False, split=True, linewidth=1, size=8, edgecolor='black')
        #ax.title()
        plt.text(0.0, 1.02, r'$n_i$=%s; $\Sigma n_i=%s$'%(','.join(list(map(str,nvec))),sum(nvec)),
                 transform=ax.transAxes, size=9)
        ax.tick_params(axis='both', which='major', labelsize=8)
        ax.set_xlabel(xlabel,size=8)
        ax.set_ylabel(v + ' [' + unitspretty[v] + ']',size=8)
        ax.set_ylim(minval*0.98,maxval*1.02)
        #ax.get_yaxis().set_label_coords(-0.15, 0.5)

        # legend:
        # place it outside (hints from https://stackoverflow.com/questions/4700614/how-to-put-the-legend-out-of-the-plot)
        # box = ax.get_position()
        # ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        if vno == 0:
            ax.legend_.remove()
        #     handles, labels = ax.get_legend_handles_labels()
        #     lgd = ax.legend(reversed(handles[0:2]), reversed(labels[0:2]),
        #                     loc='center left',
        #                     fontsize=12,  # 'large',
        #                     handletextpad=0.5,
        #                     bbox_to_anchor=(1., 0.5))
        #     lgd.legendHandles[0]._sizes = [40]
        #     lgd.legendHandles[1]._sizes = [40]
        else:
            ax.legend_.remove()

        fname = os.path.join(plfpath, fnameroot+envmethod + '_' +plottype+'_' + v + '.png')
        plt.savefig(fname, dpi=300)
        #print('      plotted:' + fname)
        print ('ok', end =" ")
        plt.close()
    print ('')

def plot_Taylor(fname, rects, colors, stdrefs, samples, legend=True,x95=0, y95=0, x99=0, y99=0):

    # figure size according to row&column numbers
    if len(rects) < 4:
        r = 1
    elif len(rects) >= 4:
        r = 2
    c = np.ceil(len(rects) / r)
    w, h =(3.0*c,3.0*r)
    fig = plt.figure(figsize=(w, h), dpi=150)

    #fig.suptitle("Precipitations", size='x-large')

    sampleskeys=list(samples.keys())
    if legend:
        markers = ['$%d$' % (i + 1) for i in len(samples[sampleskeys[0]])]
    else:
        L=samples[sampleskeys[0]]
        markers= ['$%s$'%row[2] for row in L]
        colors=['k' for i in range(len(colors))]

    if all(np.array(list(stdrefs.values()))==1.0):
        norm=True

    for season in sampleskeys:

        #x-axis lim (std)
        Mstds = [row[0]/stdrefs[season] if not np.isnan(row[1]) else 0.0 for row in samples[season]]
        stdfact = np.amax(np.array(Mstds))*1.05

        dia = TaylorDiagram(stdrefs[season], fig=fig, rect=rects[season],
                            label='Reference',stdfact=stdfact,norm=norm)

        # lines for correlation coefficients
        for cc in np.concatenate((np.arange(0.1,0.91,0.1),[0.95,0.99])):
            smax=stdfact*stdrefs[season]
            if smax<1.0:
                smax=1.05
            dia.ax.plot([np.arccos(cc),np.arccos(cc)], [0.0,smax],color='0.5')

        # Add samples to Taylor diagram
        for i, (stddev, corrcoef, name) in enumerate(samples[season]):
            if np.isnan(corrcoef): #this means there was not enough data
                continue

            if legend:
                dia.add_sample(stddev, corrcoef,
                           #marker='$%d$' % (i + 1),
                           marker=markers[i],
                           ms=12, ls='',
                           # mfc='k', mec='k', # B&W
                           mfc=colors[i], mec=colors[i],  # Colors
                           label=name)
            else:
                dia.add_sample_wlab(stddev, corrcoef,
                               # marker='$%d$' % (i + 1),
                               label=name,
                               color=colors[i],
                               marker='o',
                               ms=6, ls='',
                               # mfc='k', mec='k', # B&W
                               mfc=colors[i], mec=colors[i],  # Colors
                               )

        # Add RMS contours, and label them
        contours = dia.add_contours(levels=5, colors='0.5')  # 5 levels
        dia.ax.clabel(contours, inline=1, fontsize=8, fmt='%.1f')

        # Tricky: ax is the polar ax (used for plots), _ax is the
        # container (used for layout)
        #title:
        #dia._ax.set_title(season,fontsize=8.) #.capitalize()

    # Add a figure legend and title. For loc option, place x,y tuple inside [ ].
    # Can also use special options here:
    # http://matplotlib.sourceforge.net/users/legend_guide.html

    if legend:
        fig.legend(dia.samplePoints,
                   [p.get_label() for p in dia.samplePoints],
                   numpoints=1, prop=dict(size='small'), loc='center')

    fig.tight_layout()

    plt.savefig(fname, dpi=150)
    #plt.show()
    print ('      plotted:' + fname)

def scatter_refmod(Vset,Vstat,Uset,plfpath,fnameroot,xylabels,title='',logtransform=True):
    # expected structure: V[var]={'dates','lats','lons','ref','sim'}
    print('scatterplots:')
    plotext='.png'
    for v in Vset.keys():
        x = Vset[v]['ref']
        y = Vset[v]['model']

        Ustr = ' [' + Uset[v] + ']'

        f = plt.figure(figsize=(3.0, 3.0), dpi=150)
        ax = f.add_axes([0.2, 0.2, 0.6, 0.6])
        cbaxes = f.add_axes([0.85, 0.2, 0.03, 0.6])
        minval=min(min(x),min(y))*0.98
        maxval=max(max(x),max(y))*1.02
        xedges, yedges = np.linspace(minval, maxval, 40), np.linspace(minval, maxval, 40)
        hist, xedges, yedges = np.histogram2d(x, y, (xedges, yedges))
        xidx = np.clip(np.digitize(x, xedges), 0, hist.shape[0] - 1)
        yidx = np.clip(np.digitize(y, yedges), 0, hist.shape[1] - 1)
        plt.sca(ax)
        cmap = truncate_colormap(plt.get_cmap('gray_r'), 0.4, 1.0) 
        #pcf=plt.scatter(x, y, c=hist[xidx, yidx],cmap=plt.get_cmap('viridis'),lw=0,s=3)
        pcf=plt.scatter(x, y, c=hist[xidx, yidx],cmap=cmap,lw=0,s=3)
        plt.plot([minval,maxval],[minval,maxval],'-k')
        ax.tick_params(axis='both', which='major', labelsize=8)
        if logtransform and (v in ['Chl']):
            ax.set_yscale('log')
            ax.set_xscale('log')
            logsuf='_log'
        else:
            logsuf=''
        ax.set_xlim([minval,maxval])
        ax.set_ylim([minval,maxval])
        plt.xlabel(xylabels[0]+' '+v+Ustr,size=8)
        plt.ylabel(xylabels[1]+' '+v+Ustr,size=8)
        ax.set_title(title)
        #plt.text(0.03, 0.92, r'$\rho$: %3.2f'%Vstat[v]['rho'], ha='left', va='center',size=8, transform=ax.transAxes)
        #plt.text(0.03, 0.85, r'$n$: %s' %Vstat[v]['n'], ha='left', va='center', size=8, transform=ax.transAxes)

        #plt.text(0.68, 0.19, r'$B*$: %3.2f' % Vstat[v]['B*'], ha='left', va='center', size=8, transform=ax.transAxes)
        #plt.text(0.68, 0.12, r'$\rho$: %3.2f' % Vstat[v]['rho'], ha='left', va='center', size=8, transform=ax.transAxes)
        #plt.text(0.68, 0.05, r'$n$: %s' %Vstat[v]['n'], ha='left', va='center', size=8, transform=ax.transAxes)
        if np.isnan(Vstat[v]['B*']):
            plt.text(0.0, 1.02,r'$\rho$: %3.2f, $n$: %d' % (Vstat[v]['B*'], Vstat[v]['rho'], Vstat[v]['n']),transform=ax.transAxes, size=9)
        else:
            plt.text(0.0,1.02, r'$B^*$: %3.2f, $\rho$: %3.2f, $n$: %d'%(Vstat[v]['B*'],Vstat[v]['rho'],Vstat[v]['n']), transform=ax.transAxes, size=9)

        cb = plt.colorbar(pcf, cax=cbaxes)  # , orientation='horizontal'
        cb.ax.tick_params(axis='both', which='major', labelsize=8)
        cb.ax.set_title('count', fontsize=8)

        fname = os.path.join(plfpath,fnameroot + '_' + v+logsuf + plotext)
        plt.savefig(fname, dpi=300)
        plt.close()
        print ('      plotted:' + fname)

def arrange4Taylor(stats, type, s='', v=''):

    if type == 'V':
        # expected structure: stats[var] = {'Rstd', 'Mstd', 'rho', 'n'}
        vars = stats.keys()

        stdrefs = {};
        samples = {};
        rects = {}
        samplesV = []
        for v in vars:
            samplesV.append([stats[v]['Mstd'] / stats[v]['Rstd'], stats[v]['rho'], v])
        samples[s] = samplesV
        stdrefs[s] = 1.0  # i.e., normalized
        rects[s] = 111  # (1,1,1)

        # colors
        # Colormap (see http://www.scipy.org/Cookbook/Matplotlib/Show_colormaps)
        colors = plt.matplotlib.cm.Set1(np.linspace(0, 1, len(vars)))

    elif type == 'SxV':
        # expected structure: stats[season][var] = {'Rstd', 'Mstd', 'rho', 'n'}
        seasons = stats.keys()
        vars = stats[seasons[0]].keys()

        # find row&column numbers
        if len(seasons) < 4:
            r = 1
        elif len(seasons) >= 4:
            r = 2
        c = np.ceil(len(seasons) / r)

        stdrefs = {};
        samples = {};
        rects = {}
        for sno, s in enumerate(seasons):
            samplesV = []
            for v in vars:
                samplesV.append([stats['Mstd'] / stats['Rstd'], stats['rho'], v])
            samples[s] = samplesV
            stdrefs[s] = 1.0  # i.e., normalized
            rects[s] = r * 100 + c * 10 + v  # todo: any py function to construct the code?

        # colors
        # Colormap (see http://www.scipy.org/Cookbook/Matplotlib/Show_colormaps)
        colors = plt.matplotlib.cm.Set1(np.linspace(0, 1, len(vars)))

    elif type == 'SxM':
        # expected structure: stats[model][season][var]={'Rstd','Mstd','rho','n'}
        models = stats.keys()
        seasons = stats[models[0]].keys()
        vars = stats[models[0]][seasons[0]].keys()
        print ('Taylor plots in ' + type + ' layout not yet implemented')

    return (stdrefs, samples, rects, colors)

def get_stats(V):
    #expected structure: V[var]={'dates','lats','lons','ref','sim'}

    Vstats={}
    for v in V.keys():
        ref=V[v]['ref']
        model=V[v]['model']
        stats={}
        stats['Rstd'] = np.std(ref)
        stats['Mstd'] = np.std(model)
        stats['n'] = len(ref)
        stats['rho']=np.corrcoef([ref,model])[0,1]
        if np.mean(ref)>-0.0001 and np.mean(ref)<0.0001: #eg, because the values are anomolies
            stats['B']=np.nan
            stats['B*']=np.nan
        else:
            stats['B'] = np.mean(model) - np.mean(ref)
            stats['B*'] = (np.mean(model) - np.mean(ref))/(np.mean(ref))
        #todo: other statistics?
        Vstats[v]=stats
    return Vstats

def group_seasonsXregions (matchups,seasons,regions):
    print ('Grouping Seasons X Regions')
    Rdict = {'NW': [[54,56.0],[0,4.0]], 'NE': [[55.0,56.0],[4.0,10.0]],
              'SW': [[51.0,54.0],[0,5.0]], 'SE': [[51.0,55.0],[5.0,10.0]],
             'W': [[51.0,56.0],[0,5.0]], 'E': [[51.0,56.0],[5.0,10.0]],
             'WadSea':[[53.2, 55.1],[6.2,9.3]], 'all': []}
    seasondefs = {'all': np.arange(1, 13), 'winter': np.array([1, 2, 3, 10, 11, 12]),
                  'growing': np.arange(4, 10), 'spring': np.arange(4, 7), 'summer': np.arange(7, 9)}
    SVRset = {}
    for s in seasons:
        months = seasondefs[s]
        VRset={}
        for v in matchups.keys():
            Rset = {}
            for r in regions:
                if r == 'all':
                    # just take it all
                    ri=np.arange(len(matchups[v]['dates']))
                else:
                    #time filter
                    datesM = [date.month for date in matchups[v]['dates']]
                    mi = np.in1d(datesM, months)
                    #regional filter
                    latint=Rdict[r][0]
                    lonint=Rdict[r][1]
                    ri=(matchups[v]['lats']>=latint[0]) * (matchups[v]['lats']<=latint[1]) *(matchups[v]['lons']>=lonint[0]) * (matchups[v]['lons']<=lonint[1])
                    mri=mi*ri
                #Rset[r] = {'dates': matchups[v]['dates'][mri], 'lats': matchups[v]['lats'][mri], 'lons': matchups[v]['lons'][mri],
                #           'ref': matchups[v]['ref'][mri], 'model': matchups[v]['model'][mri]}
                Rset[r] = {}
                for par in matchups[v].keys():
                    Rset[r][par] = matchups[v][par][mri]
            VRset[v]=Rset
        SVRset[s]=VRset
    return SVRset

def group_regions (matchups,regions):
    print ('Grouping Regions')
    Rdict = {'NW': [[53.5, 56.0], [0, 5.0]], 'NE': [[54.5, 56.0], [5.0, 10.0]],
             'SW': [[51.0, 53.5], [0, 5.0]], 'SE': [[51.0, 54.5], [5.0, 10.0]],
             'W': [[51.0, 56.0], [0, 6.0]], 'E': [[51.0, 56.0], [6.0, 10.0]],
             'WadSea':[[53.2, 55.1],[6.2,9.3]], 'all': []}
    VRset={}
    for v in matchups.keys():
        Rset = {}
        for r in regions:
            if r == 'all':
                # just take it all
                ri=np.arange(len(matchups[v]['dates']))
            else:
                latint=Rdict[r][0]
                lonint=Rdict[r][1]
                ri=(matchups[v]['lats']>=latint[0]) * (matchups[v]['lats']<=latint[1]) *(matchups[v]['lons']>=lonint[0]) * (matchups[v]['lons']<=lonint[1])
            #Rset[r] = {'dates': matchups[v]['dates'][ri], 'lats': matchups[v]['lats'][ri], 'lons': matchups[v]['lons'][ri],
            #           'ref': matchups[v]['ref'][ri], 'model': matchups[v]['model'][ri]}
            Rset[v] = {}
            for par in matchups[v].keys():
                Rset[v][par] = matchups[v][par][ri]
        VRset[v]=Rset
    return VRset

def group_seasons (matchups,seasons):
    print ('Grouping Seasons')
    seasondefs = {'all': np.arange(1, 13), 'winter':np.array([1,2,3,10,11,12]),
                  'growing':np.arange(4,10),'spring':np.arange(4,7),'summer':np.arange(7,9)}
    SVset={}
    for s in seasons:
        if s=='all':
            SVset[s]=matchups #just take it all
        else:
            Vset={}
            months=seasondefs[s]
            for v in matchups.keys():
                datesM=[date.month for date in matchups[v]['dates']]
                mi=np.in1d(datesM,months)
                #Vset[v] = {'dates': matchups[v]['dates'][mi], 'lats': matchups[v]['lats'][mi], 'lons': matchups[v]['lons'][mi],
                #                 'ref': matchups[v]['ref'][mi], 'model': matchups[v]['model'][mi]}
                Vset[v] = {}
                for par in matchups[v].keys():
                    Vset[v][par] = matchups[v][par][mi]
            SVset[s]=Vset
    return SVset

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

def get_stats_demo():
    # Reference std
    stdrefs = dict(winter=48.491,
                   spring=44.927,
                   summer=37.664,
                   autumn=41.589)

    # Sample std,rho: Be sure to check order and that correct numbers are placed!
    samples = dict(winter=[[17.831, 0.360, "CCSM CRCM"],
                           [27.062, 0.360, "CCSM MM5"]],
                   spring=[[32.174, -0.262, "CCSM CRCM"],
                           [24.042, -0.055, "CCSM MM5"]],
                   summer=[[35.863, 0.096, "CCSM CRCM"],
                           [43.771, 0.367, "CCSM MM5"]],
                   autumn=[[27.374, 0.150, "CCSM CRCM"],
                           [20.270, 0.451, "CCSM MM5"]])

    # Here set placement of the points marking 95th and 99th significance
    # levels. For more than 102 samples (degrees freedom > 100), critical
    # correlation levels are 0.195 and 0.254 for 95th and 99th
    # significance levels respectively. Set these by eyeball using the
    # standard deviation x and y axis.

    #x95 = [0.01, 0.68] # For Tair, this is for 95th level (r = 0.195)
    #y95 = [0.0, 3.45]
    #x99 = [0.01, 0.95] # For Tair, this is for 99th level (r = 0.254)
    #y99 = [0.0, 3.45]

    x95 = [0.05, 13.9] # For Prcp, this is for 95th level (r = 0.195)
    y95 = [0.0, 71.0]
    x99 = [0.05, 19.0] # For Prcp, this is for 99th level (r = 0.254)
    y99 = [0.0, 70.0]

    return (stdrefs,samples,x95,y95,x99,y99)

if __name__=='__main__':

    if len(sys.argv)>1:
        files=sys.argv[1].split(',')
        #print '   file0:%s\n   file1:%s'%(files[0],files[1])
    else:
        # ICES:
        files=['/home/onur/WORK/projects/GB/data/ices/lon-1-10_lat50-57/raw/2010-2014/BGC_data_block.pickle',
               #'/home/onur/WORK/projects/2013/gpmeh/sns144-GPMEH-PPZZ-P190529-fSG97dChl-sedPS/extract_skillC_sns144-GPMEH-P190529-fSG97dChl.2012-2013_zSB.nc'
               '/home/onur/WORK/projects/2013/gpmeh/sns144-GPMEH-G191216-Fnew3-PPZZSi-vS-P191223/extract_skillC_sns144-GPMEH-G191216-Fnew3-PPZZSi-vS-P191223.2010-2014_zSB.nc']
        #files = ['/home/onur/WORK/projects/GB/data/ices/lon-1-10_lat50-57/raw/TS-2008-2015/TS_data_block.pickle',,
        #         '/home/onur/WORK/projects/2013/gpmeh/sns144-GPMEH-Fnew/extract_MphysC_sns144-GPMEH-Fnew.2012_zSB.nc']
        # ESA-CCI:
        #files = ['/home/onur/WORK/projects/GB/data/satellite/ESA-CCI-v3.1-monthly/CCI_ALL-v3.1-MONTHLY_2008-2010_rgSNS_137x94_ymonmean_mean.nc',
        #        '/home/onur/WORK/projects/GB/maecs/3d/sns144-M161117n-P161118-bdyi3-z01mm-wAtmN/sns144-M161117n-P161118-bdyi3-z01mm-wAtmN-attchl_2008-2010_ymonmean_mean.nc']

    if len(sys.argv)>2:
        ftypes=sys.argv[2].split(',')
        for i,ftype in enumerate(ftypes):
            if 'ESA-CCI' in ftype:
                ftypes[i]='ESA-CCI'
            elif 'GETM.M-MAECS' in ftype:
                ftypes[i]='GETM-MAECS'
            elif 'GETM-GPMEH' in ftype:
                ftypes[i]='GETM-GPMEH'
    else:
        #ftypes = ['ICES', 'MOSSCO-MAECS']
        #ftypes=['ICES','GETM.M-MAECS']
        ftypes = ['ICES', 'GETM.M-GPMEH']
        #ftypes = ['ESA-CCI', 'GETM.M-MAECS']
    
    #todo: eliminate fnames (and use only ftypes)
    if len(sys.argv)>3:
        fnames=sys.argv[3].split(',')
        print ('   fnames0:%s\n   fnames1:%s' % (fnames[0], fnames[1]))
    else:
        #fnames=['ICES','MOSSCO-MAECS']
        #fnames = ['ICES', 'GETM.M-MAECS']
        fnames = ['ICES', 'GETM-GPM.PPZZ-EH.GC']
        #fnames = ['ESA-CCI', 'GETM-MAECS']

    if len(sys.argv) > 4:
        varnames = sys.argv[4].split(',')
    else:
        #varnames=['Kd']
        varnames = ['Chl'] #,'NO3','NH4','DIP']
        #varnames = ['S','T']

    if len(sys.argv) > 5:
        plfpath = sys.argv[5]
    else:
        plfpath = os.path.join(os.path.dirname(files[1]),'3Dstats')

    if not os.path.exists(plfpath):
        os.mkdir(plfpath)

    if len(sys.argv) > 6:
        vertloc = sys.argv[6]
    else:
        vertloc = 'surf' #surf,bot,all

    if len(sys.argv) > 7:
        yrsS = sys.argv[7].split(',')
        yrs=[np.int(y) for y in yrsS]   
    else:
        yrs=[2011,2014]

    if len(sys.argv) > 8:
        remOLmeth=sys.argv[8]
    else:
        remOLmeth='middle' #none,percentile,middle

    if len(sys.argv) > 9:
        seasons = sys.argv[9].split(',')
    else:
        #seasons=['growing','winter','spring','summer']
        seasons = ['all']

    if len(sys.argv) > 10:
        split = sys.argv[10]
    else:
        split='none' #none,regions,vertloc

    if len(sys.argv)>11:
        scatter=True if sys.argv[11] == '1' else False
    else:
        scatter=True

    if len(sys.argv)>12:
        taylor = True if sys.argv[12] == '1' else False
    else:
        taylor=False

    if len(sys.argv)>13:
        barplotsacross = sys.argv[13]
    else:
        barplotsacross='none' #'none',botdepth,salinity

    main(files,ftypes,fnames,varnames,plfpath,vertloc,yrs,seasons,split,scatter,taylor,barplotsacross,remOLmeth,demo=False)