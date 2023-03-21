def combineDIN(obsdata,tstart,tend)

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

    #drop duplicates
    obsdata_NO3.drop_duplicates(['datetime'],inplace=True)
    obsdata_NH4.drop_duplicates(['datetime'],inplace=True)

    t1 = obsdata_NO3['datetime'].reset_index(drop=True)
    t2 = obsdata_NH4['datetime'].reset_index(drop=True)
    if len(t1) == len(t2):
            if (t1 == t2).all():
                obsdata_var = pd.DataFrame(data={
                    'datetime': obsdata_NO3['datetime'],
                    'value': obsdata_NO3['value'] + obsdata_NH4['value'],
                    'vertical_reference_preflabel':obsdata_NO3.vertical_reference_preflabel,
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

        obsdata_var = pd.DataFrame(data={
            'datetime':obsdata_NO3['datetime'][ind_NO3],
            'value':obsdata_NO3['value'][ind_NO3] + obsdata_NH4['value'][ind_NH4],
            'vertical_reference_preflabel':obsdata_NO3.vertical_reference_preflabel,
            'p35code':'EPC00198',
            'p35preflabel':'Water body nitrate plus ammonium',
            'unit_preflabel':'Milligrams per litre'
            })
        
    return obsdata_var
