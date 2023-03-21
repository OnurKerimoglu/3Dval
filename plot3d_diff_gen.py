'''
This routine creates difference maps. It is inspired by DELTARES scripts.
Created by Daniel Thewes in Oct. 22, for the project InterReg Water Quality
'''

import matplotlib as mpl
from matplotlib import cm
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import netCDF4
import os, sys
import copy
import numpy as np
# import netcdftime
# import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature


# functions
def define_varlims(scen1, scen2, var):
    if 'HS' not in scen1 and 'HS' not in scen2:
        # par_bounds_rel = {'Chl': [-25, 30, 5], 'DIN': [-50, 60, 10], 'DIP': [-2.5, 3, 0.5]}
        # # par_bounds_rel = {'Chl': [-50, 60, 5], 'DIN': 2 * [-100, 120, 20], 'DIP': 2 * [-5, 6, 1]}
        lim_vars = {'Chl': 100, 'DIN': 50, 'DIP': 50}
        cb_interval_dict = {'Chl': 10, 'DIN': 5, 'DIP': 5}
        cb_interval = int(cb_interval_dict[var])
        cbticks_dict = {'Chl': 10, 'DIN': 10, 'DIP': 10}
        cbticks = cbticks_dict[var]
    elif ('HS' in scen1 and not 'HS' in scen2) or ('HS' not in scen1 and 'HS' in scen2):
        # par_bounds_rel = {'Chl': [-25, 30, 5], 'DIN': [-50, 60, 10], 'DIP': [-2.5, 3, 0.5]}
        # par_bounds_rel = {'Chl': [-50, 60, 10], 'DIN': [-100, 120, 20], 'DIP': [-5, 6, 1]}
        lim_vars = {'Chl': 50, 'DIN': 100, 'DIP': 80}
        cb_interval_dict = {'Chl': 5, 'DIN': 10, 'DIP': 10}
        cb_interval = int(cb_interval_dict[var])
        cbticks_dict = {'Chl': 10, 'DIN': 20, 'DIP': 20}
        cbticks = cbticks_dict[var]
    else:
        lim_vars = {'Chl': 100, 'DIN': 50, 'DIP': 80}
        cb_interval_dict = {'Chl': 10, 'DIN': 10, 'DIP': 10}
        cb_interval = int(cb_interval_dict[var])
        cbticks_dict = {'Chl': 10, 'DIN': 10, 'DIP': 10}
        cbticks = cbticks_dict[var]

    return lim_vars, cb_interval, cbticks


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


def cmap_to_colors(cmap):
    ''' Convert a colormap to a list of RGB tuples. '''
    # write colors to a numpy array
    colors = cmap(np.arange(0, cmap.N))
    # convert numpy array to a list of tuples
    colors = list(map(tuple, colors))

    return colors


def colors_to_cmap(name, colors, bins=256):
    ''' Convert a list of RGB tuples to a colormap. '''
    return LinearSegmentedColormap.from_list(name, colors, bins)


def get_cmap(name, bins=256):
    ''' Return a colormap of n bins.'''
    return cm.get_cmap(name, bins)


# central file names and folders
# dataroot = "/work/ku0646/g260105/IR/"
dataroot = "/home/daniel/levante_work/IR/"
# simroot = "sns144-GPMEH-G200124-Fnew3-PPPMZZ-vS-ICGEMO-CS-BCdcsmP-rivWS"

# scenarios
scenout = {'2t': 'CS', '28': '2.8m', '28M': '2.8o', 'HS1': 'HS1', 'HS2': 'HS2', 'r': 'r', 't4': 't4', '2g-CS': '2g-CS',
           '4g-CS': '4g-CS', '2g-CS-NEC': '2g-CS-NEC', '4g-CS-NEC': '4g-CS-NEC'}

# central plotting switches (CHANGE_HERE):
cbarorient = 'horizontal'
figuresize = (11.69, 8.27)
dpi = 120

units = "%"
prettytitle = {'Chl': 'summer mean Chl-a', 'DIP': 'winter mean DIP', 'DIN': 'winter mean DIN'}

# Plot extent
ext = {'SNS': {'lat': [52.5, 56.], 'lon': [4., 9.]},
       'FSK': {'lat': [52.5, 54., ], 'lon': [5., 7.5]}}


def main(scen1, scen2, varns, plfpath):
    scen1out = scenout[scen1]
    scen2out = scenout[scen2]
    # open data file and create variables object
    if 'HS' in scen1:
        simroot1 = f'sns144-GPMEH-G200124-Fnew3-PPPMZZ-vS-ICGEMO-{scen1}-BCdcsmP-rivWS'
    elif '2t' in scen1 or '28' in scen1:
        simroot1 = f'sns144-GPMEH-G200124-Fnew3-PPPMZZ-vS-ICGEMO-CS-BCdcsmP-rivWS-{scen1}'
    else:
        simroot1 = f'sns144-{scen1}'
    if 'HS' in scen2:
        simroot2 = f'sns144-GPMEH-G200124-Fnew3-PPPMZZ-vS-ICGEMO-{scen2}-BCdcsmP-rivWS'
    elif '2t' in scen1 or '28' in scen1:
        simroot2 = f'sns144-GPMEH-G200124-Fnew3-PPPMZZ-vS-ICGEMO-CS-BCdcsmP-rivWS-{scen2}'
    else:
        simroot2 = f'sns144-{scen2}'

    ncfile1 = f"{dataroot}{simroot1}/extract_skillCS_{simroot1}.2017-avgout.nc"
    ncfile2 = f"{dataroot}{simroot2}/extract_skillCS_{simroot2}.2017-avgout.nc"

    nc1 = netCDF4.Dataset(ncfile1)
    ncv1 = nc1.variables
    nc2 = netCDF4.Dataset(ncfile2)
    ncv2 = nc2.variables

    lons = ncv1['lon'][:]
    lats = ncv1['lat'][:]

    for e, coords in ext.items():
        for vari, varn in enumerate(varns):
            print('varn: ' + varn)
            lim_vars, cb_interval, cbticks = define_varlims(scen1, scen2, varn)
            f = plt.figure(figsize=figuresize, dpi=dpi)
            spec = GridSpec(1, 1, figure=f)

            ax1 = f.add_subplot(spec[0], projection=ccrs.PlateCarree())

            # create an axes
            divider = make_axes_locatable(ax1)
            # bottom of the map
            # The width of cax will be 5% of ax and the padding between cax and ax will be fixed at 0.3 inch.
            cax = divider.append_axes("bottom", size="5%", pad=0.8, axes_class=plt.Axes)

            # get colormap
            num_colors = int(lim_vars[varn] * 2 / cb_interval)
            cmp = get_cmap('RdBu_r', num_colors)
            # get list of colors
            cmp_colors = cmap_to_colors(cmp)
            # make range around zero white
            cmp_colors[int(num_colors / 2 - 1)], cmp_colors[int(num_colors / 2)] = (1.0, 1.0, 1.0, 1.0), (
                1.0, 1.0, 1.0, 1.0)
            # create new cmap
            cmp_new = colors_to_cmap('RdBu_r_edit', cmp_colors, num_colors)
            #
            # # function to normalize values
            norm = TwoSlopeNorm(vcenter=0, vmin=-1. * lim_vars[varn], vmax=lim_vars[varn])

            # extopt = 'neither'
            extopt = 'both'

            # generate output variable from raw data
            if ncv1[varn].ndim > 2:
                var1 = np.squeeze(ncv1[varn][0, :, :])
                var2 = np.squeeze(ncv2[varn][0, :, :])
            else:
                var1 = np.squeeze(ncv1[varn][:, :])
                var2 = np.squeeze(ncv2[varn][:, :])
            var = 100 * (var2 - var1) / var1

            # pcf = ax1.contourf(lons, lats, var, transform=ccrs.PlateCarree(),
            #                    levels=bounds, norm=norm, cmap=cmp_new, extend=extopt)
            pcf = ax1.pcolor(lons, lats, var, transform=ccrs.PlateCarree(),
                             norm=norm, cmap=cmp_new)

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
            cbar = add_colorbar(f, cax, plt.cm.ScalarMappable(cmap=cmp_new, norm=norm),
                                "Difference (%)", 10, spacing='proportional',
                                orientation='horizontal', extend='both',
                                ticks=list(np.arange(-lim_vars[varn], lim_vars[varn] + 1, cbticks)))
            # cbar = add_colorbar(f, cax, plt.cm.ScalarMappable(cmap=cmp_new, norm=norm),
            #                     "Difference (%)", 10, spacing='proportional',
            #                     orientation='horizontal', extend='both',
            #                     ticks=cbticks)

            titlestr = f"{prettytitle[varn]} 100*(B-A)/A, A={scen2out}, B={scen2out}"
            ax1.set_title(titlestr, size=11)

            # f.subplots_adjust(wspace=0.15, hspace=0.7)
            plfname = f"diff_{varn}_{scen1out}-{scen2out}_{e}.png"
            print(plfname)
            if plfpath == '':
                plfpath = os.path.join(os.path.dirname(ncfile1), 'figures', 'diffmaps')
            if not os.path.isdir(plfpath):
                os.system('mkdir -p %s' % (plfpath))

            plt.savefig(os.path.join(plfpath, plfname), dpi=dpi, bbox_inches='tight')
            print(f'{plfname} saved!')
            plt.close(f)


if __name__ == '__main__':
    # get commandline options:
    if len(sys.argv) > 1:
        # ncfile1 = sys.argv[1]
        scen1 = sys.argv[1]
    else:
        # scen1 = "2t"
        # scen1 = "HS1"
        # scen1 = 'r'
        # scen1 = '2g-CS'
        scen1 = '2g-CS-NEC'
    if len(sys.argv) > 2:
        scen2 = sys.argv[2]
    else:
        # scen2 = "28"
        # scen2 = "HS1"
        # scen2 = "HS2"
        # scen2 = 't4'
        # scen2 = '4g-CS'
        scen2 = '4g-CS-NEC'

    if len(sys.argv) > 3:
        varns = [sys.argv[3].split(',')[0]]
    else:
        varns = ['Chl', 'DIN', 'DIP']
        # varns = ['Chl','DIN']

    if len(sys.argv) > 4:
        plfpath = os.path.join(os.path.dirname(ncfile1), sys.argv[4])
    else:
        # plfpath = '/work/ku0646/g260105/IR/Harmonization/diffmaps'
        plfpath = ''

    main(scen1, scen2, varns, plfpath)
