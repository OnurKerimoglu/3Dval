'''
This routine creates concentration maps. It is inspired by DELTARES scripts.
created by Daniel Thewes, Oct. '22, for the InterReg Water Quality project.
'''

import copy
import os
import sys

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib as mpl
import matplotlib.pyplot as plt
import netCDF4
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1 import make_axes_locatable

from NWDM_funcs import get_obs_avg


# functions
def add_colorbar(fig, ax, data, labelname, size, **kwargs):
    ''' function to create a colorbar '''

    # add colorbar to axes
    cb = fig.colorbar(data, cax=ax, **kwargs)
    # add label to colorbar
    cb.set_label(label=labelname, fontsize=10)
    cb.ax.xaxis.set_label_position('top')
    # set ticklabel font size
    cb.ax.tick_params(labelsize=size)
    # set position and font size scientific notation
    cb.ax.xaxis.get_offset_text().set_fontsize(5)

    return cb


# def get_obs(varn):
#     obsfile = {'DIN': 'nuts_win_ave_edit_anouk.csv', 'DIP': 'nuts_win_ave_edit_anouk.csv', 'Chl': 'chlfa_sum_ave_anouk.csv'}
#     obscols = {'DIN': 'DIN', 'DIP': 'DIP', 'Chl': 'Chlorophyll_mean'}
#     obscol_n = {'DIN': 'DIN_n', 'DIP': 'DIP_n', 'Chl': 'Chlorophyll_n'}
#     n_var = {'DIN': 2, 'DIP': 2, 'Chl': 6}
#     # import observations
#     obs = import_csv(dir_obs, obsfile[varn])
#     # create geodataframe, point geometry based on lat and lng
#
#     # obs_gdf = gpd.GeoDataFrame(obs, geometry=gpd.points_from_xy(obs.Longitude_mean, obs.Latitude_mean),
#     #                           crs="EPSG:4326")
#     df_geometry = [geometry.Point(xy) for xy in zip(obs.Longitude_mean, obs.Latitude_mean)]
#     obs_gdf  = gpd.GeoDataFrame(obs, geometry=df_geometry,crs="EPSG:4326")
#
#     # remove points without data
#     obs_var = obs_gdf[['location_name', 'Latitude_mean','Longitude_mean',
#                        obscols[varn], obscol_n[varn], 'geometry']]
#     obs_var = obs_var.loc[obs_var[obscol_n[varn]] > n_var[varn]].reset_index(drop=True)
#
#     obs_vals=obs_var[obscols[varn]]
#
#     # convert coordinates
#     #obs_lng, obs_lat = map(np.array(obs_var.Longitude_mean),np.array(obs_var.Latitude_mean))
#     return obs_vals,np.array(obs_var.Latitude_mean),np.array(obs_var.Longitude_mean)

def import_csv(path, file, sep=';', **kwargs):
    ''' import data from csv file '''
    df = pd.read_csv(os.path.join(path, file), sep=sep, **kwargs)
    return df


# central file names and folders
# dataroot = "/work/ku0646/g260105/IR/"
dataroot = "/home/daniel/levante_work/IR/"

# "datafiles" needs to be with lower case "d"
# dir_obs='/work/ku0646/UBA/obsdata/reprocessed_ICES/'
# dir_obs=f'{dataroot}/Harmonization/reprocessed_ICES/'
# dir_obs = f'{dataroot}/Harmonization/NWDMdata_InterregStations_20142017/datafiles/'
dir_obs = f'{dataroot}Harmonization/NWDM/datafiles/'

# scenarios
scenout = {'2t': 'CS', '28': '2.8m', '28M': '2.8o', 'HS1': 'HS1', 'HS2': 'HS2'}
scenoutfolders = {'2t': 'CS', '28': '28', '28M': '28M', 'HS1': 'HS1', 'HS2': 'HS2'}

# central plotting switches (CHANGE_HERE):
cbarorient = 'horizontal'
figuresize = (11.69, 8.27)
dpi = 120

unitdict = {'Chl': '$\mu g/l$', 'DIN': '$\mu MN$', 'DIP': '$\mu MP$', 'salt': 'PSU'}
prettytitle = {'Chl': 'summer mean Chl-a', 'DIP': 'winter mean DIP', 'DIN': 'winter mean DIN',
               'salt': 'annual mean salinity'}
par_bounds_rel = {'Chl': [-25, 30, 5], 'DIN': [-50, 60, 10], 'DIP': [-2.5, 3, 0.5], 'salt': [15, 36, 3]}
lim_vars = {'Chl': 20, 'DIN': 100, 'DIP': 2, 'salt': 33}
cb_interval = {'Chl': 2, 'DIN': 10, 'DIP': 0.2, 'salt': 3}

# Plot extent
ext = {'SNS': {'lat': [52.5, 56.], 'lon': [4., 9.]},
       'FSK': {'lat': [52.5, 54., ], 'lon': [5., 7.5]}}


def main(scen, varns, plfpath):
    scenoutstr = scenout[scen]
    # open data file and create variables object
    if 'HS' in scen:
        simroot = f"sns144-GPMEH-G200124-Fnew3-PPPMZZ-vS-ICGEMO-{scen}-BCdcsmP-rivWS"
    else:
        simroot = f"sns144-GPMEH-G200124-Fnew3-PPPMZZ-vS-ICGEMO-CS-BCdcsmP-rivWS-{scen}"

    # ncfile1 = f"{dataroot}{simroot}/extract_skillCS_{simroot}.2017-avgout.nc"
    ncfile1 = f"{dataroot}{simroot}/extract_MphysCS_{simroot}.2017-avgout.nc"
    nc1 = netCDF4.Dataset(ncfile1)
    ncv1 = nc1.variables

    lons = ncv1['lon'][:]
    lats = ncv1['lat'][:]

    for e, coords in ext.items():
        for vari, varn in enumerate(varns):
            print('varn: ' + varn)
            f = plt.figure(figsize=figuresize, dpi=dpi)
            spec = GridSpec(1, 1, figure=f)

            ax1 = f.add_subplot(spec[0], projection=ccrs.PlateCarree())

            # create an axes
            divider = make_axes_locatable(ax1)
            # bottom of the map
            # The width of cax will be 5% of ax and the padding between cax and ax will be fixed at 0.3 inch.
            cax = divider.append_axes("bottom", size="5%", pad=0.8, axes_class=plt.Axes)

            # cmap = mpl.colors.Colormap("viridis")
            cmap = copy.copy(mpl.cm.get_cmap("viridis", lut=int(lim_vars[varn] / cb_interval[varn])))
            cmap.set_over('gold')
            # cmap.colorbar_extend='max'

            # # function to normalize values
            if varn == 'salt':
                cmin = 15
                norm = mpl.colors.BoundaryNorm(np.arange(cmin, lim_vars[varn] + cb_interval[varn], cb_interval[varn]),
                                               cmap.N + 2, extend='both')
            else:
                cmin = 0
                norm = mpl.colors.BoundaryNorm(np.arange(cmin, lim_vars[varn] + cb_interval[varn], cb_interval[varn]),
                                               cmap.N + 1, extend='max')

            # generate output variable from raw data
            if ncv1[varn].ndim > 2:
                var1 = np.squeeze(ncv1[varn][0, :, :])
            else:
                var1 = np.squeeze(ncv1[varn][:, :])
            var = var1

            # pcf = ax1.contourf(lons, lats, var, transform=ccrs.PlateCarree(),
            #                    levels=bounds, norm=norm, cmap=cmp_new, extend=extopt)
            pcf = ax1.pcolor(lons, lats, var, transform=ccrs.PlateCarree(), cmap=cmap, norm=norm)

            if scen == '2t':
                # Plot Observations
                # obs_val, obs_lat, obs_lng = get_obs(varn)
                # ax1.scatter(obs_lng, obs_lat, c=obs_val, cmap=cmap,
                #             norm=norm,
                #             edgecolor='black', s=30, linewidth=0.2, zorder=3,
                #             transform=ccrs.PlateCarree())
                if varn == 'salt':
                    gdfobs = get_obs_avg(varn, dir_obs, avgwindow='annual', avgmode='mean',
                                         years=list(range(2017, 2017 + 1)))
                    gdfobs.plot(ax=ax1, column='value', cax=cax, linewidths=1.5, edgecolor='black', cmap=cmap,
                                norm=norm)
                else:
                    gdfobs = get_obs_avg(varn, dir_obs, avgwindow='seasonal', avgmode='mean',
                                         years=list(range(2017, 2017 + 1)))
                    gdfobs.plot(ax=ax1, column='value', cax=cax, linewidths=1.5, edgecolor='black', cmap=cmap,
                                norm=norm)

            # define extent of map
            ax1.set_extent([coords['lon'][0], coords['lon'][1],
                            coords['lat'][0], coords['lat'][1]],
                           ccrs.PlateCarree())

            # add coastlines and gridlines and fill land
            ax1.coastlines(resolution='10m')
            ax1.gridlines(draw_labels=True)
            ax1.add_feature(cfeature.LAND)
            # ax1.add_feature(cfeature.OCEAN)
            # ax1.add_feature(cfeature.LAND)
            # ax1.add_feature(cfeature.COASTLINE)
            # ax1.add_feature(cfeature.BORDERS)

            ax1.set_aspect('equal')

            # enable grid
            ax1.grid(True, linewidth=.25)

            # set label
            ax1.set_xlabel('Longitude [°]', fontsize=10)
            ax1.set_ylabel('Latitude [°]', fontsize=10)

            # add colorbar
            # cbar = add_colorbar(f, cax, plt.cm.ScalarMappable(cmap=cmap, norm=norm),
            #                     unitdict[varn], 10, spacing='proportional',
            #                     orientation='horizontal', extend='max',
            #                     ticks=list(np.arange(0, int(lim_vars[varn]+cb_interval[varn]*2),
            #                                          int(2*cb_interval[varn]))))
            if varn == 'DIP':
                cbstrformat = "%1.2f"
            else:
                cbstrformat = "%d"

            if varn == 'salt':
                cbar = add_colorbar(f, cax, pcf,
                                    unitdict[varn], 10, spacing='proportional',
                                    orientation='horizontal', extend='both',
                                    format=cbstrformat,
                                    ticks=list(np.arange(cmin, lim_vars[varn] + cb_interval[varn] * 2,
                                                         2 * cb_interval[varn])))
            else:
                cbar = add_colorbar(f, cax, pcf,
                                    unitdict[varn], 10, spacing='proportional',
                                    orientation='horizontal', extend='max',
                                    format=cbstrformat,
                                    ticks=list(np.arange(0, lim_vars[varn] + cb_interval[varn] * 2,
                                                         2 * cb_interval[varn])))

            # cbar = mpl.colorbar.ColorbarBase(cax, cmap=cmap,
            #                                  spacing='proportional', orientation='horizontal',
            #                                  ticks=list(np.arange(0, int(lim_vars[varn]),
            #                                                       int(lim_vars[varn] / cb_interval[varn]))))

            titlestr = f'{prettytitle[varn]} {scenoutstr}'
            ax1.set_title(titlestr, size=11)

            # f.subplots_adjust(wspace=0.15, hspace=0.7)
            plfname = f"conc_{varn}_{scenoutstr.replace('.', '_')}_{e}.png"
            # if plfpath == '':
            #     plfpath = os.path.join(os.path.dirname(ncfile), 'varmap')
            # if not os.path.isdir(plfpath):
            #     os.system('mkdir -p %s' % (plfpath))
            plfpath = f"{dataroot}Harmonization/{scenoutfolders[scen]}/"

            plt.savefig(os.path.join(plfpath, plfname), dpi=dpi, bbox_inches='tight')
            print(f'{plfname} saved!')
            plt.close(f)


if __name__ == '__main__':
    # get commandline options:
    if len(sys.argv) > 1:
        # ncfile1 = sys.argv[1]
        scen = sys.argv[1]
    else:
        scen = "2t"
        # scen = "28"
        # scen = "28M"
        # scen = "HS1"
        # scen = "HS2"

    if len(sys.argv) > 2:
        varns = [sys.argv[2].split(',')[0]]
    else:
        # varns = ['Chl', 'DIN', 'DIP']
        varns = ['salt']

    if len(sys.argv) > 3:
        plfpath = os.path.join(os.path.dirname(ncfile1), sys.argv[3])
    else:
        # plfpath = '/work/ku0646/g260105/IR/Harmonization/diffmaps'
        plfpath = '/home/daniel/levante_work/IR/Harmonization/maps'

    main(scen, varns, plfpath)
