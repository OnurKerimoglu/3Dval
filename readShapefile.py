# Routine to read in WFD shapefile and write out mask for GETM in NetCDF format
# Author: Daniel Thewes, 18.7.22

# modules
import netCDF4 as nc
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
from shapely.geometry import Point, Polygon
import os
import csv

# Functions

def write_NC_mask(outnamenc,keep_old_switch,lat,lon,WFDmask,area):
    print('writing NC outfile')
    try: outfile.close()  # just to be safe, make sure dataset is not already open.
    except: pass

    if os.path.exists(outnamenc):
        if keep_old_switch:
            os.rename(outnamenc,outnamenc[:-3]+'_old.nc')
        else:
            os.remove(outnamenc)

    outfile = nc.Dataset(outnamenc,mode = 'w', format = 'NETCDF4')

    outfile.title = 'WFD typology mask'


    x_T_dim = outfile.createDimension('x_T',137)
    y_T_dim = outfile.createDimension('y_T',94)
    x_dim = outfile.createDimension('x',137)
    y_dim = outfile.createDimension('y',94)
    # n_dim = outfile.createDimension('n',int(len(sfshp)))

    latnc = outfile.createVariable('lat', np.float32, ('y', 'x'))
    latnc.units = 'degrees north'
    latnc.long_name = 'latitude'
    lonnc = outfile.createVariable('lon', np.float32, ('y', 'x'))
    lonnc.units = 'degrees east'
    lonnc.long_name = 'longitude'
    WFDmasknc = outfile.createVariable('WFD', np.float32, ('y_T', 'x_T'))
    WFDmasknc.units = 'None'
    WFDmasknc.long_name = 'WFD typology applied to GETM SNS grid'
    areanc = outfile.createVariable('area', np.float32, ('y_T', 'x_T'))
    areanc.units = 'm2'
    areanc.long_name = 'cell surface areas'
    # WFDnamesnc = outfile.createVariable('WFDnames','c',('n',))
    # WFDnamesnc.units = 'None'
    # WFDnamesnc.long_name = 'Names of WFD areas'
    # WFDtypenc = outfile.createVariable('WFDtype','c',('n',))
    # WFDtypenc.units = 'None'
    # WFDtypenc.long_name = 'Type of WFD areas'

    latnc[:,:] = np.array(lat)
    lonnc[:,:] = np.array(lon)
    WFDmasknc[:,:] = np.array(WFDmask)
    areanc[:,:] = np.array(area)
    # WFDnamesnc[:] = WFDnames
    # WFDtypenc[:] = WFDtype

    outfile.close()  
    
def write_CSV_file(outnamecsv,sfshp,WFDnames,WFDtype,keep_old_switch):
    print('writing CSV outfile')
    
    if os.path.exists(outnamecsv):
        if keep_old_switch:
            os.rename(outnamecsv,outnamecsv[:-4]+'_old.csv')
        else:
            os.remove(outnamecsv)
            
    f = open(outnamecsv,'w', encoding='UTF8', newline='')
    writer = csv.writer(f,delimiter=';')

    header = ['Name','Type']
    writer.writerow(header)

    for i in range(len(sfshp)):
        row = [WFDnames[i],WFDtype[i]]
        # print(row)
        writer.writerow(row)

    f.close()
    
# Switches
keep_old_switch = False # True: overwrites old file, False: renames old file and creates new one

# file names
WFDname = "WFD_waterbody_all"
stationsname = "PP_Stations_waterbody"
bathynamecut = "/home/daniel/levante_work/IR/Bathymetry/topo_sns144_cut_C.nc"
bathyname = "/home/daniel/levante_work/IR/Bathymetry/topo_HR2.nc"
outnamenc = "/home/daniel/levante_work/IR/Bathymetry/WFDshapefile.nc"
outnamecsv = "/home/daniel/levante_work/IR/Bathymetry/WFDshapefile.csv"

print('shapefile is: '+WFDname)
print('bathymetry is: '+bathyname)

# Reading files
print('reading shapefile')
sfshp = gpd.read_file(WFDname+'.shp')

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
print('making mask')
WFDmask=np.zeros((Nx,Ny))
areac=np.zeros((Nx,Ny))
for x in range(Nx):
    for y in range(Ny):
        p=Point(lon[x,y],lat[x,y])
        tmp=np.array(sfshp.contains(p))
        if any(tmp == True):
            WFDmask[x,y]=int(np.where(tmp==True)[0])
    print(str(np.ceil(x/Nx*100))+"%")
WFDmask[WFDmask==0]=np.NaN
WFDnames=sfshp.EU_CD_CW.values
for i in range(len(WFDnames)): WFDnames[i] = str(WFDnames[i]).replace('\xe5','a')
WFDtype=sfshp.NEA_TYPE.values

# Writing NC outfile
write_NC_mask(outnamenc,keep_old_switch,lat,lon,WFDmask,area)
# Writing CSV outfile
write_CSV_file(outnamecsv,sfshp,WFDnames,WFDtype,keep_old_switch)

print('done!')
