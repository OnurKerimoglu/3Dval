import netCDF4
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import point
from shapely import wkt
import io
import os
from matplotlib import pyplot as plt
import datetime
from general_funcs import get_botdepth

# get this here: python -m pip install git+https://github.com/lkschn/WQ_tools
#                python -m pip install git+https://github.com/Deltares/dfm_tools
from WQ_tools.nwdmFunctions import wfsbuild, readUrl
from WQ_tools.dwaqFunctions import get_modkey, get_modTime
#from WQ_tools.plotFunctions import plotTS_modelNWDM  # Apparent compatibility issue

# -*- coding: utf-8 -*-
"""
Created on Thu Feb 10 13:22:36 2022

dictionaries for variables and locations

@author: schne_la
"""

# How do I comment this correct?

# dictionary for the timeseries plot
# add layer to dictionary
# add minimum value as well
# add depth (surface and unknown)
varDict = {
    'Chlfa': {'mod_vars': 'Chlfa',
              'p35code': 'EPC00105',
              'lab_vars': 'CHLF-a (ug/L)',
              'out_vars': 'CHLFa',
              'conv_vars': 1.,
              'depth': 'sea level'},
    'NO3': {'mod_vars': 'NO3',
            'p35code': 'EPC00004',
            'lab_vars': 'NO\u2083 (mgN/l)',
            'out_vars': 'NO3',
            'conv_vars': 1.,
            'depth': 'sea level'},
    'NH4': {'mod_vars': 'NH4',
            'p35code': 'EPC00009',
            'lab_vars': 'NH\u2084 (mgN/l)',
            'out_vars': 'NH4',
            'conv_vars': 1.,
            'depth': 'sea level'},
    'DIN': {'mod_vars': 'DIN',
            'p35code': 'EPC00198',
            'lab_vars': 'DIN (mgN/l)',
            'out_vars': 'DIN',
            'conv_vars': 1.,
            'depth': 'sea level'},
    'PO4': {'mod_vars': 'PO4',
            'p35code': 'EPC00007',
            'lab_vars': 'PO\u2084 (mgP/l)',
            'out_vars': 'PO4',
            'conv_vars': 1.,
            'depth': 'sea level'},
    'Si': {'mod_vars': 'Si',
           'p35code': 'EPC00008',
           'lab_vars': 'Si (mgSi/L)',
           'out_vars': 'Si',
           'conv_vars': 1.,
           'depth': 'sea level'},
    'OXY': {'mod_vars': 'OXY',
            'p35code': 'EPC00002',
            'lab_vars': 'O\u2082 (mg/L)',
            'out_vars': 'O2',
            'conv_vars': 1.,
            'depth': 'sea level'},
    'IM1': {'mod_vars': 'IM1',
            'p35code': 'EXTRA004',
            'lab_vars': 'SPM (mg/L)',
            'out_vars': 'SPM',
            'conv_vars': 1.},
    'POC1': {'mod_vars': 'POC1',
             'p35code': 'EPC00157',
             'lab_vars': 'POC1 (mgC/l)',
             'out_vars': 'POC1',
             'conv_vars': 1.,
             'depth': 'sea level'},
    'PON': {'mod_vars': 'PON',
            'p35code': 'EPC00212',
            'lab_vars': 'PON (mgN/l)',
            'out_vars': 'PON',
            'conv_vars': 1.,
            'depth': 'sea level'},
    'POP1': {'mod_vars': 'POP1',
             'p35code': 'EPC00201',
             'lab_vars': 'POP1 (mgP/l)',
             'out_vars': 'POP1',
             'conv_vars': 1.,
             'depth': 'sea level'},
    'salinity': {'mod_vars': 'salinity',
                 'p35code': 'EPC00001',
                 'lab_vars': 'Salinity (g/kg)',
                 'out_vars': 'sal',
                 'conv_vars': 1.,
                 'depth': 'sea level'},
    'temperature': {'mod_vars': 'temperature',
                    'p35code': 'WATERTEMP',
                    'lab_vars': 'Temperature (celcius)',
                    'out_vars': 'temp',
                    'conv_vars': 1.,
                    'depth': 'sea level'},
    'TotN': {'mod_vars': 'TotN',
             'p35code': 'EPC00134',
             'lab_vars': 'TotN (mgN/l)',
             'out_vars': 'TotN',
             'conv_vars': 1.,
             'depth': 'sea level'},
    'TotP': {'mod_vars': 'TotP',
             'p35code': 'EPC00135',
             'lab_vars': 'TotP (mgP/l)',
             'out_vars': 'TotP',
             'conv_vars': 1.,
             'depth': 'sea level'},
    'pH': {'mod_vars': 'pH',
           'p35code': 'EPC00168',
           'lab_vars': 'pH',
           'out_vars': 'pH',
           'conv_vars': 1.,
           'depth': 'sea level'},
    'pCO2water': {'mod_vars': 'pCO2water',
                  'p35code': 'EPC00133',
                  'lab_vars': 'pCO2 (microatms)',
                  'out_vars': 'pCO2',
                  'conv_vars': 1.,
                  'depth': 'sea level'},
    'ExtVl': {'mod_vars': 'ExtVl',
              'p35code': 'EXTRA002',
              'lab_vars': 'extinction (1/m)',
              'out_vars': 'Ext',
              'conv_vars': 1.,
              'depth': 'sea level'},
    'O2_bottom': {'mod_vars': 'OXY',
                  'p35code': 'EPC00002',
                  'lab_vars': 'O\u2082 (mg/L)',
                  'out_vars': 'O2_bottom',
                  'conv_vars': 1.,
                  'depth': 'sea level'},
    'PP': {'mod_vars': 'fPPtot',
           'p35code': 'INPPPIC1',
           'lab_vars': 'PP (mgC/m2/d)',
           'out_vars': 'PP',
           'conv_vars': 1000.,
           'depth': 'sea level'}
}

varDict_TS = {
    'Chlfa': {'mod_vars': 'Chlfa',
              'p35code': 'EPC00105',
              'lab_vars': 'CHLF-a (ug/L)',
              'out_vars': 'CHLFa',
              'conv_vars': 1.,
              'depth': ['sea level', 'unknown']},
    'NO3': {'mod_vars': 'NO3',
            'p35code': 'EPC00004',
            'lab_vars': 'NO\u2083 (mgN/l)',
            'out_vars': 'NO3',
            'conv_vars': 1.,
            'depth': ['sea level', 'unknown']},
    'NH4': {'mod_vars': 'NH4',
            'p35code': 'EPC00009',
            'lab_vars': 'NH\u2084 (mgN/l)',
            'out_vars': 'NH4',
            'conv_vars': 1.,
            'depth': ['sea level', 'unknown']},
    'DIN': {'mod_vars': 'DIN',
            'p35code': 'EPC00198',
            'lab_vars': 'DIN (mgN/l)',
            'out_vars': 'DIN',
            'conv_vars': 1.,
            'depth': ['sea level', 'unknown']},
    'PO4': {'mod_vars': 'PO4',
            'p35code': 'EPC00007',
            'lab_vars': 'PO\u2084 (mgP/l)',
            'out_vars': 'PO4',
            'conv_vars': 1.,
            'depth': ['sea level', 'unknown']},
    'Si': {'mod_vars': 'Si',
           'p35code': 'EPC00008',
           'lab_vars': 'Si (mgSi/L)',
           'out_vars': 'Si',
           'conv_vars': 1.,
           'depth': ['sea level', 'unknown']},
    'IM1': {'mod_vars': 'IM1',
            'p35code': 'EXTRA004',
            'lab_vars': 'SPM (mg/L)',
            'out_vars': 'SPM',
            'conv_vars': 1.,
            'depth': ['sea level', 'unknown']},
    'salinity': {'mod_vars': 'salinity',
                 'p35code': 'EPC00001',
                 'lab_vars': 'Salinity (g/kg)',
                 'out_vars': 'sal',
                 'conv_vars': 1.,
                 'depth': ['sea level', 'unknown']},
    'temperature': {'mod_vars': 'temperature',
                    'p35code': 'WATERTEMP',
                    'lab_vars': 'Temperature (celcius)',
                    'out_vars': 'temp',
                    'conv_vars': 1.,
                    'depth': 'sea level'}
}

# dictionary translating NWDM p35 names to local variable names
p35dict = {
    'salinity': 'salt',
    'NO3': 'NO3',
    'NH4': 'NH4',
    'PO4': 'DIP',
    'DIN': 'DIN',
    'Si': 'Si',
    'Chlfa': 'Chl',
    'temperature': 'temp',
    'IM1': 'SPM',
    'ExtVl': 'zD',
    'OXY':'OXY'
}

# dictionary translating local variable names to NWDM p35 names
p35dict_reverse = {
    'salt': 'salinity',
    'NO3': 'NO3',
    'NH4': 'NH4',
    'DIP': 'PO4',
    'DIN': 'DIN',
    'Si': 'Si',
    'Chl': 'Chlfa',
    'temp': 'temperature',
    'SPM': 'IM1',
    'zD': 'ExtVl',
    'OXY':'OXY'
}

# plot locations in NWDM and model output
plot_locs = ['NOORDWK1', 'NOORDWK10', 'NOORDWK2', 'NOORDWK20', 'NOORDWK30', 'NOORDWK4', 'NOORDWK50', 'NOORDWK70',
             'ROTTMPT10', 'ROTTMPT100', 'ROTTMPT15', 'ROTTMPT20', 'ROTTMPT3', 'ROTTMPT30', 'ROTTMPT5', 'ROTTMPT50',
             'ROTTMPT70', 'TERSLG10', 'TERSLG100', 'TERSLG135', 'TERSLG175', 'TERSLG20', 'TERSLG235', 'TERSLG30',
             'TERSLG4', 'TERSLG50', 'TERSLG70', 'WALCRN1', 'WALCRN10', 'WALCRN2', 'WALCRN20', 'WALCRN30', 'WALCRN4',
             'WALCRN50', 'WALCRN70', '220006', '220017', '220052_B', '220052_S', '220057', '220065', 'Nney_W_1',
             'Nney_W_2', 'Nney_W_3', 'WeMu_W_1', 'WeMu_W2', 'Wu_KU-W1', 'ZOUTKPLZGT', 'ZUIDOLWOT', 'BOCHTVWTM',
             'BOOMKDP', 'DANTZGT', 'DOOVBWT', 'GROOTGND', 'HUIBGOT', 'MARSDND', 'VLIESM', 'Vlissingen_boei_SSVH',
             'JaBu_W_1', 'Bork_W_1', 'DOWSING', 'J1', 'J12', 'LIVBAY', 'M388', 'M394', 'M397',
             'M435', 'M55', 'M57', 'M765', 'N14_FALKENBERG', 'Point_2_SRN_Boulogne', 'P2', 'St._Peter_Eider',
             'TH1', 'W04', 'W05', 'A17', 'AMO_6', 'ANHOLTE', 'UKO1','GE_Sdl_Amrum','GE_Hrnum_Vortrapptief',
             'GE_Norderelbe','GE_Sderpiep','GE_Bsum']

# column names in in-situ database (use for measurement_view_all)
columnsNWDM = ("location_code", "date", "depth", "vertical_reference_preflabel",
               "vertical_reference_code", "p35code", "p35preflabel",
               "value", "unit_preflabel", "quality_code", "station", "geom")


def get_obs_avg(var, dir, avgwindow, avgmode, years):
    # obsfilelist = os.listdir(dir)
    obsfilelist = [f for f in os.listdir(dir) if os.path.isfile(os.path.join(dir, f)) and f[-3:] == 'csv']

    vardict2 = {'Chl': 'EPC00105', 'DIN': 'EPC00198', 'DIP': 'EPC00007', 'salt': 'EPC00001'}

    p_list = []
    v_list = []
    n_list = []
    c_list = []
    for fi, file in enumerate(obsfilelist):
        # print(file)
        # read CSV file as dataframe
        df = pd.read_csv(os.path.join(dir, file), sep=';')

        try:
            df['geom'] = df['geom'].apply(wkt.loads)
        except:
            df['geom'] = df['geom'].apply(wkt.loads)

        df = df.rename(columns={'geom': 'geometry'})
        aux = df['p35code'] == vardict2[var]

        if (not df.empty) and ((len(aux) > 0) and (any(aux))):
            startind = next(x for x in range(len(df.p35code)) if df.p35code[x] == vardict2[var])
            unitraw = df.unit_preflabel[startind]
            if (unitraw == 'Milligrams per litre') and (var == 'Chl'):
                convfactor = 1.0
            elif (unitraw == 'Milligrams per litre') and (var == 'DIN'):
                convfactor = 71.4
            elif (unitraw == 'Milligrams per litre') and (var == 'DIP'):
                convfactor = 32.3
            elif (unitraw == 'Parts per thousand') and (var == 'salt'):
                convfactor = 1.0
            else:
                convfactor = 1.0

            # Check if "datetime" in columns and if not fix
            if "datetime" not in df.columns:
                df["datetime"] = pd.to_datetime(df["date"])
            # get years and months
            yy = pd.DatetimeIndex(df['datetime']).year
            mm = pd.DatetimeIndex(df['datetime']).month

            # get only the variable we want
            if avgwindow == 'seasonal':
                if var == 'Chl':
                    tmp = (df['p35code'] == vardict2[var]) & ((mm >= 3) & (mm <= 9)) & (
                                (yy >= int(years[0])) & (yy <= (int(years[-1]))))
                elif var == 'DIN' or var == 'DIP':
                    tmp = (df['p35code'] == vardict2[var]) & ((mm <= 2) | (mm == 12)) & (
                                (yy >= int(years[0])) & (yy <= (int(years[-1]))))
                else:
                    tmp = (df['p35code'] == vardict2[var]) & ((mm <= 2) | (mm == 12)) & (
                                (yy >= int(years[0])) & (yy <= (int(years[-1]))))
            else:
                tmp = (df['p35code'] == vardict2[var]) & ((yy >= int(years[0])) & (yy <= (int(years[-1]))))

            if len(tmp) > 0:
                if any(tmp):
                    if avgmode == 'median':
                        tmp2 = df.value[tmp]
                        v_list.append(convfactor * np.nanmedian(tmp2))
                    elif avgmode == 'mean':
                        tmp2 = df.value[tmp]
                        v_list.append(convfactor * np.nanmean(tmp2))
                    n_list.append(df.location_code[0])
                    c_list.append(len(df.value[tmp]))
                    p_list.append(df.geometry[0])

    gout = gpd.GeoDataFrame(data={'value': v_list, 'area_code': n_list, 'count': c_list}, geometry=p_list,
                            crs='epsg:4326')
    outfolder = dir.split('datafiles')[0]
    gout.to_csv(f'{outfolder}{var}_{avgwindow}_{years[0]}-{years[-1]}.csv', sep=';', index=False)
    return gout


def wfsrequest(stations, timeint, vars, olf, maxz_switch=True):
    import requests
    from datetime import datetime

    # if stations to include are not specified,
    if len(stations) == 0:
        # TODO: update this some time to be automated
        stations = ['BOCHTVWTM','BOOMKDP','Bork_W_1','DANTZGT','DOOVBWT','GROOTGND','HUIBGOT','JaBu_W_1','MARSDND',
                       'Nney_W_1','Nney_W_2', 'Nney_W_3','ROTTMPT3','ROTTMPT70','TERSLG4','TERSLG10','WeMu_W_1',
                       'WeMu_W2', 'Wu_KU-W1','ZOUTKPLZGT','ZUIDOLWOT','GE_Sdl_Amrum','GE_Hrnum_Vortrapptief',
                    'GE_Norderelbe','GE_Sderpiep','GE_Bsum']
    else:
        # TODO: check the station attribute of each file, return the needed file names
        raise (Exception('sfile generation based on stations not yet implemented. leave the station list empty'))

    # access to database
    user = 'interreg'
    password = 'terreing21'

    # time
    tstart = timeint[0]
    tend = timeint[1]

    # chose which filters, variables and columns you want here:
    cql_dict = {}
    columns = ["date", "depth", "vertical_reference_preflabel",
               "vertical_reference_code",
               "value", "unit_preflabel"]

    obs = {}
    # loop over stations
    for ist, stationid in enumerate(stations):
        # extract data for specific location, variable, depth from NWDM
        url = wfsbuild(typename="NWDM:measurement_interreg",
                        cql_dict={'station': stationid},
                        outputFormat="csv", columns=columnsNWDM)

        # read url and return raw NWDM data
        obsdata_raw = readUrl(url, user, password)

        if not obsdata_raw.empty:
            #obtain coordinates
            aux2 = obsdata_raw['geom'][0].split('(')[1].split(')')[0].split(' ')

            lon = float(aux2[0])
            lat = float(aux2[1])

            if maxz_switch:
                maxz = get_botdepth(lon, lat)

                obs[stationid] = {
                    'longname': 'descriptive name of the station',
                    'lon': lon,
                    'lat': lat,
                    'bottom_depth': maxz,
                     }
            else:
                obs[stationid] = {
                    'longname': 'descriptive name of the station',
                    'lon': lon,
                    'lat': lat
                }
            # drop all depth levels except surface and/or unknown
            obsdata_raw = obsdata_raw[(obsdata_raw.vertical_reference_preflabel == 'unknown') + (obsdata_raw.vertical_reference_preflabel == 'sea level')]

            # extract data for timeperiod
            obsdata_raw['datetime'] = pd.to_datetime(obsdata_raw['date'], dayfirst=False)

            if obsdata_raw.empty:
                disp('bla!')
            obsdata = obsdata_raw.loc[(obsdata_raw['datetime'] >= tstart) &
                                        (obsdata_raw['datetime'] <= tend)].reset_index(drop=True)

            for dict_id, dict_info in varDict_TS.items():
                var = p35dict[dict_id]
                if var in vars:
                    print(dict_id)
                    obsdata_var = obsdata[obsdata.p35code == dict_info['p35code']]



                    if var == 'DIN' and obsdata_var.empty:
                        obsdata_var = combineDIN(obsdata, tstart, tend)
                        # obsdata_NO3 = obsdata[obsdata.p35code == varDict['NO3']['p35code']]
                        # obsdata_NH4 = obsdata[obsdata.p35code == varDict['NH4']['p35code']]
                        #
                        # obsdata_NO3 = obsdata_NO3.loc[(obsdata_NO3['datetime'] >= tstart) &
                        #                               (obsdata_NO3['datetime'] <= tend)].reset_index(drop=True)
                        # # Sort times
                        # obsdata_NO3 = obsdata_NO3.sort_values('datetime')
                        #
                        # obsdata_NH4 = obsdata_NH4.loc[(obsdata_NH4['datetime'] >= tstart) &
                        #                               (obsdata_NH4['datetime'] <= tend)].reset_index(drop=True)
                        # # Sort times
                        # obsdata_NH4 = obsdata_NH4.sort_values('datetime')
                        #
                        # #drop duplicates
                        # obsdata_NO3.drop_duplicates(['datetime'],inplace=True)
                        # obsdata_NH4.drop_duplicates(['datetime'],inplace=True)
                        #
                        # t1 = obsdata_NO3['datetime'].reset_index(drop=True)
                        # t2 = obsdata_NH4['datetime'].reset_index(drop=True)
                        # if len(t1) == len(t2):
                        #     try:
                        #         if (t1 == t2).all():
                        #             obsdata_var = pd.DataFrame(data={
                        #                 'datetime': obsdata_NO3['datetime'],
                        #                 'value': obsdata_NO3['value'] + obsdata_NH4['value'],
                        #                 'p35code': 'EPC00198',
                        #                 'p35preflabel': 'Water body nitrate plus ammonium',
                        #                 'unit_preflabel':'Milligrams per litre'
                        #                 })
                        #         else:
                        #             print('bla!')
                        #     except:
                        #         print('bla!')
                        #
                        #
                        # else:
                        #     aux1 = set(t1)
                        #     aux2 = set(t2)
                        #     aux11 = set.intersection(aux1, aux2)
                        #     ind_NO3 = [i for i in range(len(t1)) if t1[i] in aux11]
                        #     ind_NH4 = [i for i in range(len(t2)) if t2[i] in aux11]
                        #
                        #     obsdata_var = pd.DataFrame(data={
                        #         'datetime':obsdata_NO3['datetime'][ind_NO3],
                        #         'value':obsdata_NO3['value'][ind_NO3] + obsdata_NH4['value'][ind_NH4],
                        #         'p35code':'EPC00198',
                        #         'p35preflabel':'Water body nitrate plus ammonium',
                        #         'unit_preflabel':'Milligrams per litre'
                        #         })
                        #     print('bla!')

                    #else:

                    # Sort times
                    obsdata_var = obsdata_var.sort_values('datetime')

                    # reset index
                    time = obsdata_var['datetime'].reset_index(drop=True)

                    obs[stationid][var] = {'surface': {}, 'presence': False}
                    if not obsdata_var.empty:
                        # convert units
                        if (obsdata_var['unit_preflabel'].reset_index(drop=True)[0] == 'Micrograms per litre') and var == 'Chl':
                            conv_factor = 1.0
                        elif (obsdata_var['unit_preflabel'].reset_index(drop=True)[0] == 'Milligrams per litre') and var == 'DIN':
                            conv_factor = 1000/14 #1000mg/g*1mmolN/14gN -> mmolN/m3
                        elif (obsdata_var['unit_preflabel'].reset_index(drop=True)[0] == 'Milligrams per litre') and var == 'DIP':
                            conv_factor = 1000/31 #1000mg/g*1mmolN/14gN -> mmolN/m3
                        elif (obsdata_var['unit_preflabel'].reset_index(drop=True)[0] == 'Milligrams per litre') and var == 'Si':
                            conv_factor = 1000/28 #1000mg/g*1mmolN/14gN -> mmolN/m3
                        elif (obsdata_var['unit_preflabel'].reset_index(drop=True)[0] == 'Milligrams per litre') and var == 'OXY':
                            conv_factor = 1000/16 #1000mg/g*1mmolN/14gN -> mmolN/m3
                        else:
                            conv_factor = 1.0

                        if len(obsdata_var['datetime']) > 0:
                            mv = np.nanmean(obsdata_var['value'])
                            sv = np.nanstd(obsdata_var['value'])

                            obsdata_var.drop(obsdata_var[(obsdata_var.value <= mv - olf * sv) |
                                                         (obsdata_var.value >= mv + olf * sv)].index, inplace=True)
                            obsdata_var.dropna(subset=['value'],inplace=True)


                            obs[stationid][var]['presence'] = True
                            # try:
                            obs[stationid][var]['surface'] = {
                                'time': [obsdata_var['datetime'].reset_index(drop=True)[i].to_pydatetime()
                                         for i in range(len(obsdata_var.datetime))],
                                'value': obsdata_var['value'] * conv_factor,
                                'depth_interval': [0, 0]}
    return obs

def wfs2csv(outpath,stations, yint, vars):
    import requests
    from datetime import datetime
    #
    # if stations to include are not specified,
    if len(stations) == 0:
        # TODO: update this some time to be automated
        # stations = ['BOCHTVWTM','BOOMKDP','Bork_W_1','DANTZGT','DOOVBWT','GROOTGND','HUIBGOT','JaBu_W_1','MARSDND',
        #                'Nney_W_1','Nney_W_2', 'Nney_W_3','ROTTMPT3','ROTTMPT70','TERSLG4','TERSLG10','WeMu_W_1',
        #                'WeMu_W2', 'Wu_KU-W1','ZOUTKPLZGT','ZUIDOLWOT','GE_Sdl_Amrum','GE_Hrnum_Vortrapptief',
        #             'GE_Norderelbe','GE_Sderpiep','GE_Bsum']
        stations = plot_locs

    # access to database
    user = 'interreg'
    password = 'terreing21'

    # time
    timeint = [datetime(yint[0], 1, 1,0,0,0), datetime(yint[1], 12, 31,23,59,59)]
    tstart = timeint[0]
    tend = timeint[1]

    # chose which filters, variables and columns you want here:
    columns = ["location_code","geom","date", "depth", "vertical_reference_preflabel",
               "vertical_reference_code",
               "value", "unit_preflabel"]

    # get list of p35codes
    v_codes = []
    [v_codes.append(varDict[p35dict_reverse[var]]['p35code']) for i,var in enumerate(vars)]

    # loop over stations
    for ist, stationid in enumerate(stations):
        p_list = []
        v_list = []
        n_list = []
        u_list = []
        c_list = []
        cc_list = []
        t_list = []
        tt_list = []
        # loop over variables
        outfile = f'{outpath}{stationid}_NWDM_{yint[0]}_{yint[1]}.csv'
        for vi,var in enumerate(vars):
            print(varDict[p35dict_reverse[var]]['p35code'],var)
            # extract data for specific location, variable, depth from NWDM
            if var == 'DIN':
                url = wfsbuild(typename="NWDM:measurement_interreg",
                                cql_dict={'station': stationid},
                                outputFormat="csv", columns=columnsNWDM)

                # read url and return raw NWDM data
                obsdata = readUrl(url, user, password)
                obsdata['datetime'] = pd.to_datetime(obsdata['date'], dayfirst=False)
                obsdata_raw = combineDIN(obsdata, tstart, tend)
            else:
                url = wfsbuild(typename="NWDM:measurement_interreg",
                                cql_dict={'station': stationid, 'p35code':v_codes[vi]},
                                outputFormat="csv", columns=columnsNWDM)

                # read url and return raw NWDM data
                obsdata_raw = readUrl(url, user, password)
                obsdata_raw['datetime'] = pd.to_datetime(obsdata_raw['date'], dayfirst=False)

            if not obsdata_raw.empty:

                # drop all depth levels except surface and/or unknown
                obsdata_raw = obsdata_raw[(obsdata_raw.vertical_reference_preflabel == 'unknown') + (obsdata_raw.vertical_reference_preflabel == 'sea level')]

                # extract data for timeperiod
                obsdata_raw = obsdata_raw.loc[(obsdata_raw['datetime'] >= tstart) &
                                            (obsdata_raw['datetime'] <= tend)].reset_index(drop=True)

                # Sort times
                obsdata_raw = obsdata_raw.sort_values('datetime')
                # n_list.append(list(obsdata_raw['location_code']))
                # v_list.append(list(obsdata_raw['value']))
                # p_list.append(list(obsdata_raw['geom']))
                # u_list.append(list(obsdata_raw['unit_preflabel']))
                # c_list.append(list(obsdata_raw['p35code']))
                # t_list.append(list(obsdata_raw['date']))
                # tt_list.append(list(obsdata_raw['date']))
                # tt_list.append(list(obsdata_raw['datetime']))
                n_list.extend(list(obsdata_raw['location_code']))
                v_list.extend(list(obsdata_raw['value']))
                p_list.extend(list(obsdata_raw['geom']))
                u_list.extend(list(obsdata_raw['unit_preflabel']))
                c_list.extend(list(obsdata_raw['p35code']))
                cc_list.extend(list(obsdata_raw['p35preflabel']))
                t_list.extend(list(obsdata_raw['date']))
                tt_list.extend(list(obsdata_raw['date']))
                tt_list.extend(list(obsdata_raw['datetime']))
            else:
                print('empty dataframe!')

        obs = pd.DataFrame({'location_code': n_list, 'geom': p_list, 'date': t_list, 'p35code': c_list,
                            'p35preflabel': cc_list, 'unit_preflabel': u_list, 'station': n_list, 'value': v_list})
        obs.to_csv(outfile, sep=';', index=False)


def combineDIN(obsdata, tstart, tend):

    obsdata_NO3 = obsdata[obsdata.p35code == varDict['NO3']['p35code']]
    obsdata_NH4 = obsdata[obsdata.p35code == varDict['NH4']['p35code']]

    obsdata_NO3 = obsdata_NO3.loc[(obsdata_NO3['datetime'] >= tstart) &
                                  (obsdata_NO3['datetime'] <= tend)].reset_index(drop=True)
    # Sort times
    obsdata_NO3 = obsdata_NO3.sort_values('datetime')

    obsdata_NH4 = obsdata_NH4.loc[(obsdata_NH4['datetime'] >= tstart) &
                                  (obsdata_NH4['datetime'] <= tend)].reset_index(drop=True)
    # Sort times
    obsdata_NH4 = obsdata_NH4.sort_values('datetime')

    # drop duplicates
    obsdata_NO3.drop_duplicates(['datetime'], inplace=True)
    obsdata_NH4.drop_duplicates(['datetime'], inplace=True)

    t1 = obsdata_NO3['datetime'].reset_index(drop=True)
    t2 = obsdata_NH4['datetime'].reset_index(drop=True)
    if len(t1) == len(t2):
        if (t1 == t2).all():
            obsdata_var = pd.DataFrame(data={
                'datetime': obsdata_NO3['datetime'],
                'date': obsdata_NO3['date'],
                'value': obsdata_NO3['value'] + obsdata_NH4['value'],
                'vertical_reference_preflabel':obsdata_NO3.vertical_reference_preflabel,
                'location_code':obsdata_NO3['location_code'],
                'geom':obsdata_NO3['geom'],
                'p35code': 'EPC00198',
                'p35preflabel': 'Water body nitrate plus ammonium',
                'unit_preflabel':'Milligrams per litre'
                })
        else:
            print('bla!')


    else:
        aux1 = set(t1)
        aux2 = set(t2)
        aux11 = set.intersection(aux1, aux2)
        ind_NO3 = [i for i in range(len(t1)) if t1[i] in aux11]
        ind_NH4 = [i for i in range(len(t2)) if t2[i] in aux11]
        obsdata_NO3.reset_index(drop=True, inplace=True)
        obsdata_NH4.reset_index(drop=True, inplace=True)
        obsdata_var = pd.DataFrame(data={
            'datetime': obsdata_NO3['datetime'][ind_NO3],
            'date': obsdata_NO3['date'][ind_NO3],
            'value': obsdata_NO3['value'][ind_NO3] + obsdata_NH4['value'][ind_NH4],
            'vertical_reference_preflabel':obsdata_NO3.vertical_reference_preflabel[ind_NO3],
            'location_code':obsdata_NO3['location_code'][ind_NO3],
            'geom':obsdata_NO3['geom'][ind_NO3],
            'p35code': 'EPC00198',
            'p35preflabel': 'Water body nitrate plus ammonium',
            'unit_preflabel': 'Milligrams per litre'
        })

    return obsdata_var
