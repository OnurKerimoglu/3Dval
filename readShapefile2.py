# Routine to read in WFD shapefile, join GETM grid with it and calculate area means
# This routine is a GETM specific reproduction of a routine doing the same thing for the DCSM model by Deltares.
# It generates two output files (per year):
# CSV-file that is directly comparable to Deltares output
# Shapefile, based on the WFD shape file, but with added columns for area means and standard deviations
# Author: Daniel Thewes, 26.7.22

import netCDF4 as nc
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
from shapely.geometry import Point, Polygon
import math
import os
import fiona
import csv

# Switches
keep_old_switch = False # True: overwrites old file, False: renames old file and creates new one

# Functions
# Function to create path if path does not exist
def makedir(path):
    if not os.path.exists(path):
        print("path doesn't exist. trying to make")
        os.makedirs(path)

# Function to calculate average and standard deviation
def calc_avg_and_std(inp,var,gdf,ind):
    x_out = np.array(gdf['x'][ind])
    y_out = np.array(gdf['y'][ind])
    a_out = np.array(gdf['area'][ind])
    if var=='EH_abioP_DINO3':
        tmp = np.array(inp.variables['EH_abioP_DINO3'])+np.array(inp.variables['EH_abioP_DINH4'])
    else:
        tmp = np.array(inp.variables[var])
    asum = sum(a_out)

    v = 0
    for i in range(len(x_out)):
        x = x_out[i]
        y = y_out[i]
        a = a_out[i]/asum
        if not np.isnan(tmp[x,y]) and tmp[x,y]>=0:
            v += tmp[x,y]*a
        else:
            print(f'lon={lon[x,y]}, lat={lat[x,y]}: NaN!')

    v2 = 0
    for i in range(len(x_out)):
        x = x_out[i]
        y = y_out[i]
        a = a_out[i] / asum
        if not np.isnan(tmp[x,y]) and tmp[x,y]>=0:
            v2 += a*(tmp[x,y]-v)**2

    return(v,math.sqrt(v2))

def write_dataframe(path, filename, df, index=False, sep=';', **kwargs):
    makedir(path)
    # Write dataframe to csv-file
    df.to_csv(os.path.join(path, filename), index=False, sep=sep, **kwargs)


# environment variables
scenario_list=["CS","28"]
yearlist=[2017]
varlistout={'DIN':'DIN','DIP':'DIP','Chl':'CHL'}

# file names
dataroot="/home/daniel/levante_work/IR/Harmonization/"
WFDname = "WFD_waterbody_all"
stationsname = "PP_Stations_waterbody"
bathynamecut = "/home/daniel/levante_work/IR/Bathymetry/topo_sns144_cut_C.nc"
bathyname = "/home/daniel/levante_work/IR/Bathymetry/topo_HR2.nc"
outp = f"{dataroot}/WFD/"

print('shapefile is: ' + WFDname)
print('bathymetry is: ' + bathyname)

# Reading files
print('reading shapefile')
sfshp = gpd.read_file(WFDname + '.shp')
idcol = "EU_CD_CW"
areaids = list(sfshp[idcol])

# Reading bathymetry
print('reading bathymetry')
bathy = nc.Dataset(bathyname)
bathycut = nc.Dataset(bathynamecut)


# get lon and lat
lonf = bathy.variables['lonx']
latf = bathy.variables['latx']
xxf = bathy.variables['xx']
yxf = bathy.variables['yx']

lonc = bathycut.variables['lon']
latc = bathycut.variables['lat']
lon = lonf[5:-1,2:-1]
lat = latf[5:-1,2:-1]
xx = np.array(xxf[5:-1,2:-1])
yx = np.array(yxf[5:-1,2:-1])

[Nxf,Nyf]=np.shape(lonf)
[Nx,Ny]=np.shape(lon)
grid_type=int(np.array(bathy.variables['grid_type'])[0])

# calculate cell areas
areaf=np.zeros((Nxf,Nyf))
for x in range(Nxf):
    for y in range(Nyf):
        if grid_type == 3:
            areaf[x, y] = abs((xxf[x, y - 1] - xxf[x - 1, y]) * (yxf[x, y] - yxf[x - 1, y - 1])
                              + (xxf[x, y] - xxf[x - 1, y - 1]) * (yxf[x - 1, y] - yxf[x, y - 1]))
        else:
            # see opt/getm/src/domain/domain.F90, subroutine metric() for directions
            print('Only grid type 3 (planar curvilinear) available right now!')

area = areaf[5:-1,2:-1]

# loop over x and y to determine to which polygon each grid point belongs
p_list = []
a_list = []
x_list =[]
y_list =[]
for x in range(Nx):
    for y in range(Ny):
        p=Point(lon[x,y],lat[x,y])
        p_list.append(p)
        a=area[x,y]
        a_list.append(a)
        x_list.append(x)
        y_list.append(y)
    # print(str(np.ceil(x/Nx*100))+"%")

crs=str(sfshp.crs)

gdfdata = {'x':x_list,'y':y_list,'area':a_list}
# create geodataframe
modshp = gpd.GeoDataFrame(data = gdfdata, geometry = p_list, crs = f"EPSG:{crs[-4:]}")
modshp = modshp.rename_axis('gridcell' ).reset_index()

# join data frames
sjoin = gpd.sjoin(modshp, sfshp, op = 'within',
                          how = 'inner').reset_index(drop=True)

# Eliminate potential duplicates, resulting from overlapping polygons
sjoin_clean = sjoin.drop_duplicates(subset=['gridcell'],
                                    keep='first')

outshp = sfshp
# list grid cell-area ids mapping based on the spatial join
gridcell_areaid = list(sjoin_clean[idcol])
valid_areas = []
for i in gridcell_areaid:
    if i not in valid_areas:
        valid_areas.append(i)

for y in yearlist:
    # create output dataframe
    df_areameans = pd.DataFrame(list(areaids), columns=['Area code'])

    for scenario in scenario_list:
        filename=f"{dataroot}{scenario}/extract_skillCS_sns144-GPMEH-G200124-Fnew3-PPPMZZ-vS-ICGEMO-CS-BCdcsmP-rivWS-{scenario}.2017-avgout.nc"

        # Reading input file
        print('reading input file')
        print('input file is: ' + filename)
        inp = nc.Dataset(filename)
        varlist = inp.variables.keys()


    # loop over time
        for var in varlist:
            if var in varlistout.keys():
                print(var)


                # pre-allocate area mean and std
                area_mean, area_std = [], []

                # iterate over areas and calculate area mean
                for area in areaids:
                    print(area)
                    # get indices of grid cells centroids that fall within area bounds
                    area_cells = [icell for icell, areaid in enumerate(gridcell_areaid) if areaid==area]
                    if len(area_cells) > 0:
                        amean,astd = calc_avg_and_std(inp,var,sjoin_clean,area_cells)
                    else:
                        amean,astd = np.nan, np.nan
                    print(amean,astd)
                    if np.isnan(amean):
                        print('Warning: '+area+' is NaN!')
                    area_mean.append(amean)
                    area_std.append(astd)

                df_areameans.insert(len(df_areameans.columns), f'{var}_{scenario}', area_mean)
                df_areameans.insert(len(df_areameans.columns), f'{var}_{scenario}_std', area_std)
                outshp.insert(len(outshp.columns), f'{var}_{scenario}', area_mean)
                outshp.insert(len(outshp.columns), f'{var}_{scenario}_std', area_std)


    write_dataframe(outp, f'WFD_assessment_areameans_{y}.csv', df_areameans)
    outshp.to_file(f'{outp}WFD_assessment_areameans_{y}.shp')

print('done!')
