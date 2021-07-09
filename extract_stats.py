import sys, os, csv
import numpy as np
from getm_funcs import get_getm_dom_vars,get_getm_dataF
from general_funcs import interpval2D,get_2Dtree,getproj
sys.path.insert(0, "dataprep")
from data_tools import create_nc

#read stations
rootvalpath='/home/onur/WORK/projects/ICG-EMO/data/validation/'
fstatlist=os.path.join(rootvalpath,'template_stations_for_validation_vs2.csv')
stats={}
with open(fstatlist) as csvfile:
    spamreader = csv.reader(csvfile, delimiter=' ', quotechar='|')
    for rowno,row in enumerate(spamreader):
        if rowno>0:
            rowl=row[0].split(',')
            stats[rowl[0]]={'lat':rowl[1],'lon':rowl[2]}

#read sim file
simf=os.path.join(rootvalpath,'extract_valICGEMO_sns144-GPMEH-G200124-Fnew3-PPZZSi-vS-P191223-ICGEMO-CS.2009-2014-mm.nc')
varns=['DIN','DIP','Chl','TN','TP']
unitlib={'DIN':'mumolN/l','DIP':'mumolP/l','Chl':'mug/l','TN':'mmolN/m^3','TP':'mmolP/m^3'}
#get the domain data, construct an interpolation tree
lons,lats,bat,ysl,xsl=get_getm_dom_vars(simf)
maxlat=lats.max();minlat=lats.min();maxlon=lons.max();minlon=lons.min()
proj = getproj(setup = 'SNSfull', projpath = os.path.dirname(os.path.realpath(__file__)))
domaintree = get_2Dtree(lons,lats,proj,Vpy=2)

#get the sim data
simdata,simtime=get_getm_dataF(simf,varns,ysl,xsl,getmv='mean',modtype='GF-PPZZ')

#extract & interpolate stations
#first variable is time
varno=0;
dimslist=[{'t': 'time'}]
vals=[simtime]
names=['time'];
units=[''];
longnames=['']
#loop over the stations and variables
#statsred={'NOORDWK70':stats['NOORDWK70'],'TERSLG235':stats['TERSLG235']}
for statno,statname in enumerate(stats.keys()):
    lat=float(stats[statname]['lat'])
    lon=float(stats[statname]['lon'])
    if lat>maxlat or lat<minlat or lon>maxlon or lon<minlon:
        print ('station coords (lat:%s, lon:%s) for %s out of the model domain'%(lat,lon,statname))
        continue
    else:
        print('%s is within domain. Extracting..'%statname)
    for varn in varns:
        print(varn)
        varno=varno+1
        data = np.zeros(len(simtime)) * np.nan
        for ti in range(len(simtime)):
            if varn in ['DIN','DIP','Chl']:
                data[ti] = interpval2D(0, 0, simdata[varn][ti, 0, :, :], lat, lon, 'pretree', proj, domaintree)
            else:
                data[ti] = interpval2D(0, 0, simdata[varn][ti, :, :], lat, lon, 'pretree', proj, domaintree)
        dimslist.append({'t': 'time'})
        vals.append(data)
        names.append(varn+'_'+statname)
        longnames.append(varn)
        units.append(unitlib[varn])

dimvals={'t':simtime}

#write nc file
fname=os.path.join(rootvalpath,'ICGEMO_validation_stations_GPM.nc')

create_nc(fname,dimvals,dimslist,vals,names,longnames,units,history='Input file:%s'%simf)