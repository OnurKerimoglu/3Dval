
import copy
import os
import sys

import netCDF4
import numpy as np
import cartopy.feature as cfeature
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1 import make_axes_locatable
import cartopy.crs as ccrs

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

# central plotting switches (CHANGE_HERE):
cbarorient = 'horizontal'
figuresize = (11.69, 8.27)
dpi = 120

filename='/home/daniel/OC-CCI/CCI_CHL_l4_2017.nc'
varn='Chl'

nc = netCDF4.Dataset(filename)

var = nc.variables[varn]
lon = nc.variables['lon']
lat = nc.variables['lat']

#

unitdict = {'Chl': '$\mu g/l$', 'DIN': '$\mu MN$', 'DIP': '$\mu MP$', 'salt': 'PSU'}
prettytitle = {'Chl': 'summer mean Chl-a', 'DIP': 'winter mean DIP', 'DIN': 'winter mean DIN',
               'salt': 'annual mean salinity'}
par_bounds_rel = {'Chl': [-25, 30, 5], 'DIN': [-50, 60, 10], 'DIP': [-2.5, 3, 0.5], 'salt': [15, 36, 3]}
lim_vars = {'Chl': 25, 'DIN': 100, 'DIP': 2, 'salt': 33}
cb_interval = {'Chl': 2.5, 'DIN': 10, 'DIP': 0.2, 'salt': 3}
cmin = 0

ext = {'SNS': {'lat': [52.5, 56.], 'lon': [4., 9.]}}

for e, coords in ext.items():
    f = plt.figure(figsize=figuresize, dpi=dpi)
    spec = GridSpec(1, 1, figure=f)

    ax1 = f.add_subplot(spec[0], projection=ccrs.PlateCarree())

    # create an axes
    divider = make_axes_locatable(ax1)
    # bottom of the map
    # The width of cax will be 5% of ax and the padding between cax and ax will be fixed at 0.3 inch.
    cax = divider.append_axes("bottom", size="5%", pad=0.8, axes_class=plt.Axes)

    # cmap = mpl.colors.Colormap("viridis")
    # cmap = copy.copy(mpl.cm.get_cmap("viridis", lut=int(lim_vars[varn] / cb_interval[varn])))
    cmap = copy.copy(mpl.cm.get_cmap(mpl.colormaps["viridis"], lut=int(lim_vars[varn] / cb_interval[varn])))
    cmap.set_over('gold')

    norm = mpl.colors.BoundaryNorm(np.arange(cmin, lim_vars[varn] + cb_interval[varn], cb_interval[varn]),
                                   cmap.N + 1, extend='max')

    pcf = ax1.pcolor(lon, lat, var, transform=ccrs.PlateCarree(), cmap=cmap, norm=norm)
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
    cbstrformat = "%d"
    cbar = add_colorbar(f, cax, pcf,
            unitdict[varn], 10, spacing='proportional',
            orientation='horizontal', extend='max',
            format=cbstrformat,
            ticks=list(np.arange(0, lim_vars[varn] + cb_interval[varn] * 2,
                                 2 * cb_interval[varn])))

    titlestr = f'{prettytitle[varn]} CCI'
    ax1.set_title(titlestr, size=11)

    # f.subplots_adjust(wspace=0.15, hspace=0.7)
    plfname = f"conc_{varn}_CCI_{e}.png"
    plfpath = "/home/daniel/OC-CCI/"
    if not os.path.exists(plfpath):
        os.makedirs(plfpath)

    plt.savefig(os.path.join(plfpath, plfname), dpi=dpi, bbox_inches='tight')
print('done!')