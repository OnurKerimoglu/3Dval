import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import point
from shapely import wkt
import os
from matplotlib import pyplot as plt
import datetime

dataroot = '/home/daniel/levante_work/IR/Harmonization/'
obsfolder = 'NWDMdata_InterregStations_20142017/datafiles/'
obsfilelist = os.listdir(os.path.join(dataroot, obsfolder))

var = 'Chl'
avgmode = 'median'

# convfactor = {'Chl': 1.0, 'DIN': 7.14e+1, 'DIP': 3.23e+1}
vardict = {'Chl': 'Water body chlorophyll-a', 'DIN': 'Water body DIN', 'DIP': 'Water body phosphate'}
vardict2 = {'Chl': 'EPC00105', 'DIN': 'EPC00198', 'DIP': 'EPC00007'}

# gdfdata = {'geometry': [], 'value': []}
p_list = []
v_list = []
n_list = []
for fi, file in enumerate(obsfilelist):
    print(file)
    df = pd.read_csv(os.path.join(dataroot, obsfolder, file), sep=';')
    df['geom'] = df['geom'].apply(wkt.loads)
    df = df.rename(columns={'geom': 'geometry'})
    # df = gpd.read_file(os.path.join(dataroot, obsfolder, file), sep=';')
    gdf = gpd.GeoDataFrame(df, crs='epsg:4326')
    # get years and months
    yy = pd.DatetimeIndex(df['datetime']).year
    mm = pd.DatetimeIndex(df['datetime']).month

    # get only the variable we want
    tmp = (gdf['p35code'] == vardict2[var]) & ((mm>=3) & (mm<=9)) # check crs
    # dout['geometry'][fi] = gdf.geometry[0]
    # dout['value'][fi] = np.nanmean(gdf.value[tmp])
    if len(tmp)>0:
        print(gdf.location_code[0])
        p_list.append(gdf.geometry[0])
        if avgmode=='median':
            v_list.append(np.nanmedian(gdf.value[tmp]))
        elif avgmode=='mean':
            v_list.append(np.nanmean(gdf.value[tmp]))
        n_list.append(gdf.location_code[0])

gout = gpd.GeoDataFrame(data = {'value': v_list, 'area_code': n_list}, geometry=p_list, crs='epsg:4326')

print('done!')