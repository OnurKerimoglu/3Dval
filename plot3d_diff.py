import matplotlib as mpl
import matplotlib.pyplot as plt
import netCDF4
import os,sys
import copy
import numpy as np
import netcdftime
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# central file names and folders
dataroot= "/work/ku0646/g260105/IR/"
simroot = "sns144-GPMEH-G200124-Fnew3-PPPMZZ-vS-ICGEMO-CS-BCdcsmP-rivWS"


# central plotting switches (CHANGE_HERE):
cbarorient='horizontal'
figuresize=(16,5)
dpi=120

unitkeys = {'Chl': '$\mu g/l$', 'DIN': '$\mu MN$', 'DIP': '$\mu MP$'}
prettytitle = {'Chl': 'Summer Chl-a', 'DIP': 'Winter DIP', 'DIN': 'Winter DIN'}
#par_bounds_comp=[-100, 112.5, 12.5]
par_bounds_abs={'Chl':[-2,2.5,0.5],'DIN':[-25, 30, 5],'DIP':[-2e-2,2.5e-2,5e-3]}
par_bounds_rel={'Chl':[-15,18,3],'DIN':[-40, 50, 10],'DIP':[-2.5,3,0.5]}
#formatstr={'Chl':"%d",'DIN':"%d",'DIP':"%.1E"}

def main(scen1,scen2,varns,plfpath):
    # open data file and create variables object
    ncfile1 = f"{dataroot}{simroot}-{scen1}/extract_skillCS_{simroot}-{scen1}.2017-avgout.nc"
    ncfile2 = f"{dataroot}{simroot}-{scen2}/extract_skillCS_{simroot}-{scen2}.2017-avgout.nc"
    nc1=netCDF4.Dataset(ncfile1)
    ncv1=nc1.variables
    nc2=netCDF4.Dataset(ncfile2)
    ncv2=nc2.variables

    lons=ncv1['lon'][:]
    lats=ncv1['lat'][:]

    f = plt.figure(figsize=figuresize, dpi=dpi)

    for vari,varn in enumerate(varns):
        print('varn: ' + varn)

        ax1= f.add_subplot(2,len(varns), vari+1,projection=ccrs.PlateCarree())

        units = unitkeys[varn]
        par_bounds = par_bounds_abs[varn]
        cmap = copy.copy(mpl.cm.get_cmap("bwr"))
        extopt = 'neither'

        if ncv1[varn].ndim > 2:
            var1=np.squeeze(ncv1[varn][0,:,:])
            var2=np.squeeze(ncv2[varn][0,:,:])
        else:
            var1=np.squeeze(ncv1[varn][:,:])
            var2=np.squeeze(ncv2[varn][:,:])
        var=(var2-var1)

        bounds = np.array([i for i in np.arange(par_bounds[0], par_bounds[1], par_bounds[2])])
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

        pcf = ax1.contourf(lons, lats, var, transform=ccrs.PlateCarree(),
                          levels=bounds, norm=norm, cmap=cmap, extend=extopt)

        ax1.add_feature(cfeature.OCEAN)
        ax1.add_feature(cfeature.LAND)
        ax1.add_feature(cfeature.COASTLINE)
        ax1.add_feature(cfeature.BORDERS)

        ax1.set_extent([-0.2, 9.5, 51, 55.8], ccrs.PlateCarree())

        titlestr=prettytitle[varn] +' '+scen2+'-'+scen1
        plt.title(titlestr,size=12.0)

        #Colorbar
        pos1 = f.gca().get_position()
        if cbarorient == 'vertical':
            poscbar = [pos1.x0 + pos1.width + 0.01, pos1.y0, 0.01, pos1.height*0.9]
        elif cbarorient == 'horizontal':
            poscbar = [pos1.x0, pos1.y0-0.02, pos1.width, 0.02]
        cax = plt.axes(position=poscbar)

        if par_bounds[0]==-100:
            cbt=bounds[0::2]
        else:
            cbt=bounds[0:]
        cb= mpl.colorbar.ColorbarBase(cax,cmap=cmap,
                                 norm=norm,
                                 spacing='proportional', orientation=cbarorient,
                                 ticks=cbt, extend=extopt)

        cb.ax.tick_params(labelsize=10)
        cb.update_ticks()
        cb.ax.set_title(units, size=10.)
        cb.formatter.set_powerlimits((-2,3))

        ax1= f.add_subplot(2,len(varns), len(varns)+vari+1,projection=ccrs.PlateCarree())

        units = '%'
        par_bounds = par_bounds_rel[varn]

        var=100*(var2-var1)/var1

        bounds = np.array([i for i in np.arange(par_bounds[0], par_bounds[1], par_bounds[2])])
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

        pcf = ax1.contourf(lons, lats, var, transform=ccrs.PlateCarree(),
                          levels=bounds, norm=norm, cmap=cmap, extend=extopt)

        ax1.add_feature(cfeature.OCEAN)
        ax1.add_feature(cfeature.LAND)
        ax1.add_feature(cfeature.COASTLINE)
        ax1.add_feature(cfeature.BORDERS)

        ax1.set_extent([-0.2, 9.5, 51, 55.8], ccrs.PlateCarree())

        titlestr=prettytitle[varn] +' ('+scen2+'-'+scen1+')/'+scen1
        plt.title(titlestr,size=12.0)

        #Colorbar
        pos1 = f.gca().get_position()
        if cbarorient == 'vertical':
            poscbar = [pos1.x0 + pos1.width + 0.01, pos1.y0, 0.01, pos1.height*0.9]
        elif cbarorient == 'horizontal':
            poscbar = [pos1.x0, pos1.y0 - 0.08, pos1.width, 0.02]
        cax = plt.axes(position=poscbar)

        if par_bounds[0]==-100:
            cbt=bounds[0::2]
        else:
            cbt=bounds[0:]
        cb= mpl.colorbar.ColorbarBase(cax,cmap=cmap,
                                 norm=norm,
                                 spacing='proportional', orientation=cbarorient,
                                 ticks=cbt, extend=extopt)

        cb.ax.tick_params(labelsize=10)
        cb.update_ticks()
        cb.ax.set_title(units, size=10.)

    f.subplots_adjust(wspace=0.15, hspace=0.7)
    plfname = f"diff_%s_{scen1}-{scen2}.png" % ('-'.join(varns))
    print(plfname)
    if plfpath=='':
        plfpath =os.path.join(os.path.dirname(ncfile), 'varmap')
    os.system('mkdir -p %s' % (plfpath))
    plt.savefig(os.path.join(plfpath,plfname), dpi=dpi, bbox_inches='tight')
    plt.close(f)

if __name__=='__main__':
    # get commandline options:
    if len(sys.argv) > 1:
        #ncfile1 = sys.argv[1]
        scen1 = sys.argv[1]
    else:
        scen1 = "28"
        

    if len(sys.argv) > 2:
        scen2 = sys.argv[2]
    else:
        scen2 = "28M"
    
    if len(sys.argv) > 3:
        varns = [sys.argv[3].split(',')[0]]
    else:
        varns = ['Chl','DIN','DIP']
        #varns = ['Chl']

    if len(sys.argv) > 4:
        plfpath = os.path.join(os.path.dirname(ncfile1), sys.argv[4])
    else:
        plfpath = '/work/ku0646/g260105/IR/Harmonization/diffmaps'

    main(scen1, scen2, varns, plfpath)
