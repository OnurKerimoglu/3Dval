from NWDM_funcs import get_station_list_from_nc

dataroot = '/home/daniel/levante_work/IR/Harmonization/'

testfile = f'{dataroot}DCSM/DCSM-FM_0_5nm_waq_0000_2014-2017_CS.nc'

stations = get_station_list_from_nc(testfile)

for ist, stationid in enumerate(stations):
    if not (stationid in plot_locs):
        continue

    print(stationid)
#
#     # retrieve location of station
#     url_sts = wfsbuild(typename = "NWDM:location", cql_dict={'station' : stationid}, outputFormat = "csv")
#     obs_sts_info = readUrl(url_sts, user, password)
#     obs_sts_info = obs_sts_info.drop_duplicates(['station'])[["x", "y", "station"]]
#
#     # loop over variables in dictionary
#     for dict_id, dict_info in varDict.items():
#
#         print(dict_id)
#
#         # extract data for specific location, variable, depth from NWDM
#         url = wfsbuild(typename = "NWDM:measurement_p35_all",
#                        cql_dict={'station' : stationid,
#                                  'p35code' : dict_info['p35code'],
#                                  'vertical_reference_preflabel' : dict_info['depth']},
#                        outputFormat = "csv", columns = columnsNWDM)
#         # read url and return raw NWDM data
#         obsdata_raw = readUrl(url, user, password)
#
#         # extract data for timeperiod
#         obsdata_raw['datetime'] = pd.to_datetime(obsdata_raw['date'], dayfirst = False)
#         obsdata = obsdata_raw.loc[(obsdata_raw['datetime'] >= tstart) &
#                                   (obsdata_raw['datetime'] <= tend)].reset_index(drop=True)
#         # Sort times
#         obsdata = obsdata.sort_values('datetime')
#
#         # skip rest and clear plot if there is no in-situ data available
#         if not(len(obsdata) > 0):
#             continue
#
#         # initialize empty df for model years output
#         modData_all = pd.DataFrame()
#
#         # loop over years to extract model output
#         for iyear, year in enumerate(yearList):
#
#             hisfile = os.path.join(basedir,  f'current_{year}\DCSM-FM_0_5nm_waq_0000_his.nc')
#
#             # get list of variables and dimension from hisfile
#             varshis_pd, dimshis_pd = get_ncvardimlist(file_nc=hisfile)
#
#             # modify get_modkey to deal with other column for salinity and temperature
#             modkey = get_modkey(varshis_pd, dict_info['mod_vars'])
#
#             # Get time indexes between tstart and tend
#             mdtimes = get_modTime(hisfile)
#
#             # get model data
#             hisdata = get_ncmodeldata(file_nc=hisfile, varname=modkey, station=stationid, timestep=mdtimes, layer=-1)
#             mvals = hisdata.squeeze()
#             # convert model result unit to in-situ data unit
#             mvals = np.multiply(mvals, dict_info['conv_vars'])
#
#             # save yaers in one dataframe
#             modData_year = pd.DataFrame(data={'time': mdtimes, 'modValues': list(mvals)})
#             modData_all = modData_all.append(modData_year, ignore_index=True)
#
#         # plot timeseries for all years
#         plotTS_modelNWDM(obsdata.datetime,
#                          obsdata.value,
#                          modData_all.time,
#                          modData_all.modValues,
#                          obs_sts_info.x,
#                          obs_sts_info.y,
#                          dict_info['lab_vars'],
#                          stationid,
#                          minVal = 0,
#                          title = 'North Sea',
#                          mapExtent = [-5, 15, 50, 65]
#                          )
#         # save plot
#         plt.savefig(os.path.join(outdir, f'ts_{dict_id}_{tdate}_{stationid}.png'))
#         plt.close()
print('done!')