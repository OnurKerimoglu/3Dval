import matplotlib as mpl
import matplotlib.pyplot as plt
import os,sys
import copy
import numpy as np
import pickle
import netCDF4
import netcdftime
import datetime
import pandas as pd
import geopandas as gpd
from shapely import geometry
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# central plotting switches (CHANGE_HERE):
plottopo=False
cbarorient='horizontal'
figuresize=(20,5)
dpi=120
#dir_obs='/home/onur/WORK/projects/ICG-EMO/data/Assessment-area-validation/reprocessed_ICES/'
#dir_obs='/work/ku0646/g260105/IR/stations/InterReg/allNC/'
dir_obs='/work/ku0646/UBA/obsdata/reprocessed_ICES/'
unitkeys = {'Chl': '$\mu g/l$', 'DIN': '$\mu MN$', 'DIP': '$\mu MP$'}
prettytitle = {'Chl': 'Summer Chl-a', 'DIP': 'Winter DIP', 'DIN': 'Winter DIN'}
par_bounds_ref = {'Chl': [0,11,1], 'DIN':[0,55,5], 'DIP':[0,1.1,0.1]}
par_bounds_comp=[-100, 112.5, 12.5]

def main(ncfile,varns,simid,plfpath):
    # open data file and create variables object
    print(ncfile)
    nc=netCDF4.Dataset(ncfile)
    ncv=nc.variables

    lons=ncv['lon'][:]
    lats=ncv['lat'][:]
    H=ncv['bathymetry'][:]

    f = plt.figure(figsize=figuresize, dpi=dpi)

    for vari,varn in enumerate(varns):
        print('varn: ' + varn)

        ax1= f.add_subplot(1,len(varns), vari+1,projection=ccrs.PlateCarree())

        if simid=='ref':
            units = unitkeys[varn]
            par_bounds = par_bounds_ref[varn]
            cmap = copy.copy(mpl.cm.get_cmap("viridis"))
            cmap.set_over('gold')
            extopt = 'max'
        elif simid=='comp':
            units = '%'
            par_bounds = par_bounds_comp
            cmap = copy.copy(mpl.cm.get_cmap("bwr"))
            extopt = 'neither'
        
        if ncv[varn].ndim > 2:
            var=np.squeeze(ncv[varn][0,:,:])
        else:
            var=np.squeeze(ncv[varn][:,:])

        bounds = np.array([i for i in np.arange(par_bounds[0], par_bounds[1], par_bounds[2])])
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

        pcf = ax1.contourf(lons, lats, var, transform=ccrs.PlateCarree(),
                          levels=bounds, norm=norm, cmap=cmap, extend=extopt)

        if plottopo:
            proj.contour(xx,yx,H,colors='k',linewidths=0.5,levels=[5.0,10.,20.,30.,40.,50.,60.])

        if simid=='ref':
            #Plot Observations
            obs_val,obs_lat,obs_lng = get_obs(varn)
            ax1.scatter(obs_lng, obs_lat, c=obs_val, cmap=cmap,
                        norm=norm,
                        edgecolor='white', s=13, linewidth=0.2, zorder=3,
                        transform=ccrs.PlateCarree())
            sufobs='_ices'
        else:
            sufobs=''

        ax1.add_feature(cfeature.OCEAN)
        ax1.add_feature(cfeature.LAND)
        ax1.add_feature(cfeature.COASTLINE)
        ax1.add_feature(cfeature.BORDERS)

        ax1.set_extent([-0.2, 9.5, 51, 55.8], ccrs.PlateCarree())

        titlestr=prettytitle[varn] # +' ['+units+']'
        plt.title(titlestr,size=12.0)

        #Colorbar
        pos1 = f.gca().get_position()
        if cbarorient == 'vertical':
            poscbar = [pos1.x0 + pos1.width + 0.07, pos1.y0 + 0.09, 0.02, pos1.height * 0.7]
        elif cbarorient == 'horizontal':
            poscbar = [pos1.x0, pos1.y0 - 0.10, pos1.width, 0.02]
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

    plfname = '%s_%s%s.png' % (simid,'-'.join(varns),sufobs)
    if plfpath=='':
        plfpath =os.path.join(os.path.dirname(ncfile), 'varmap')
    os.system('mkdir -p %s' % (plfpath))
    plt.savefig(os.path.join(plfpath,plfname), dpi=dpi, bbox_inches='tight')
    plt.close(f)

def get_obs(varn):
    obsfile = {'DIN': 'nuts_win_ave_edit_anouk.csv', 'DIP': 'nuts_win_ave_edit_anouk.csv', 'Chl': 'chlfa_sum_ave_anouk.csv'}
    obscols = {'DIN': 'DIN', 'DIP': 'DIP', 'Chl': 'Chlorophyll_mean'}
    obscol_n = {'DIN': 'DIN_n', 'DIP': 'DIP_n', 'Chl': 'Chlorophyll_n'}
    n_var = {'DIN': 2, 'DIP': 2, 'Chl': 6}
    # import observations
    obs = import_csv(dir_obs, obsfile[varn])
    # create geodataframe, point geometry based on lat and lng

    # obs_gdf = gpd.GeoDataFrame(obs, geometry=gpd.points_from_xy(obs.Longitude_mean, obs.Latitude_mean),
    #                           crs="EPSG:4326")
    df_geometry = [geometry.Point(xy) for xy in zip(obs.Longitude_mean, obs.Latitude_mean)]
    obs_gdf  = gpd.GeoDataFrame(obs, geometry=df_geometry,crs="EPSG:4326")

    # remove points without data
    obs_var = obs_gdf[['location_name', 'Latitude_mean','Longitude_mean',
                       obscols[varn], obscol_n[varn], 'geometry']]
    obs_var = obs_var.loc[obs_var[obscol_n[varn]] > n_var[varn]].reset_index(drop=True)

    obs_vals=obs_var[obscols[varn]]

    # convert coordinates
    #obs_lng, obs_lat = map(np.array(obs_var.Longitude_mean),np.array(obs_var.Latitude_mean))
    return obs_vals,np.array(obs_var.Latitude_mean),np.array(obs_var.Longitude_mean)

def import_csv(path, file, sep = ';', **kwargs):
    ''' import data from csv file '''
    df = pd.read_csv(os.path.join(path, file), sep = sep, **kwargs)
    return df

if __name__=='__main__':
    # get commandline options:
    if len(sys.argv) > 1:
        ncfile = sys.argv[1]
    else:
        #ncfile = '/home/onur/WORK/projects/ICG-EMO/data/GPM/CS/CS_09-14.nc'
        #ncfile = '/home/onur/WORK/projects/ICG-EMO/data/GPM/CS/CS_09-14-PERCDIF-HS1_09-14-REL-CS_09-14.nc'
        #ncfile = '/home/onur/WORK/projects/ICG-EMO/data/GPM/HS1/HS1_09-14-PERCDIF-HS2_09-14-REL-CS_09-14.nc'
        ncfile = '/work/ku0646/g260105/IR/sns144-GPMEH-G200124-Fnew3-PPPMZZ-vS-ICGEMO-CS-BCdcsmP-rivWS-28M/extract_skillCS_sns144-GPMEH-G200124-Fnew3-PPPMZZ-vS-ICGEMO-CS-BCdcsmP-rivWS-28M.2017-avgout.nc'

    if len(sys.argv) > 2:
        varns = [sys.argv[2].split(',')[0]]
    else:
        varns = ['Chl','DIN','DIP']
        #varns = ['Chl']

    if len(sys.argv) > 3:
        simid = sys.argv[3]
    else:
        simid =  'ref' # 'comp'  # ref,comp

    if len(sys.argv) > 4:
        plfpath = os.path.join(os.path.dirname(ncfile), sys.argv[4])
    else:
        plfpath = ''

    main(ncfile, varns, simid, plfpath)
