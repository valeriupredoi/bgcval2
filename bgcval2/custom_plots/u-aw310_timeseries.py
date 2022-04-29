# Global : AMOC 26N, ACC@Drake,  SST
# S. Ocean mean: Air-Sea CO2 flux, surface DIC, SH Ice extent, surface nitrate, surface O2, surface salinity, SST
# S. Atlantic Ocean mean : Surface salinity, surface temp, volume integrated temperature
# Weddel and Wilkes region: SST, SSS, Mixed layer depth, volume integrated temp.

import numpy as np
from matplotlib import pyplot
from shelve import open as shopen


def movingaverage_DT(data, times, window_len=10., window_units='years'):
    window_units = window_units.lower()
    if window_units not in ['days', 'months', 'years']:
        raise ValueError("movingaverage_DT: window_units not recognised" +
                         str(window_units))

    data = np.ma.array(data)
    times = np.ma.array(times)

    if len(data) != len(times):
        raise ValueError(
            "movingaverage_DT: Data and times are different lengths.")

    #####
    # Assuming time
    if window_units in [
            'years',
    ]: window = float(window_len) / 2.
    if window_units in [
            'months',
    ]: window = float(window_len) / (2. * 12.)
    if window_units in [
            'days',
    ]: window = float(window_len) / (2. * 365.25)

    output = []  #np.ma.zeros(data.shape)
    for i, t in enumerate(times):

        tmin = t - window
        tmax = t + window
        arr = np.ma.masked_where((times < tmin) + (times > tmax), data)

        #print [i,t],[tmin,tmax],[t,data[i]], arr.mean(), 'mask:',arr.mask.sum()
        output.append(arr.mean())

    return np.array(output)


def ensemblemean(dicts, operator='mean'):
    times = {}
    for job, dic in list(dicts.items()):
        print(job, len(dic))
        for t, d in list(dic.items()):
            try:
                times[t].append(d)
            except:
                times[t] = [
                    d,
                ]
    check_lengths = {}
    for t in sorted(times):
        try:
            check_lengths[len(times[t])].append(t)
        except:
            check_lengths[len(times[t])] = [
                t,
            ]

    times_keys = sorted(times.keys())
    if operator == 'mean':
        means = [np.mean(times[t]) for t in times_keys]
        return times_keys, means
    if operator == 'min':
        for t in times_keys:
            print(operator, t, times[t], np.min(times[t]))
        means = [np.min(times[t]) for t in times_keys]
        return times_keys, means
    if operator == 'max':
        means = [np.max(times[t]) for t in times_keys]
        return times_keys, means


def getAMOCdata(j, field='AMOC'):
    prefix = '/group_workspaces/jasmin2/ukesm/BGC_data/ldemora/shelves/timeseries/' + j + '/'
    if field == 'AMOC':
        fn = prefix + j + '_AMOC_26N.shelve'
        index = ('regionless', 'layerless', 'metricless')
    if field == 'Drake':
        fn = prefix + j + '_DrakePassageTransport.shelve'
        index = ('regionless', 'layerless', 'metricless')
    if field == 'GVT':
        fn = prefix + j + '_GlobalMeanTemperature.shelve'
        index = ('regionless', 'layerless', 'metricless')
    if field == 'GMT':
        fn = prefix + j + '_Temperature.shelve'
        index = ('Global', 'Surface', 'mean')
    if field == 'AirSeaFluxCO2':
        fn = prefix + j + '_TotalAirSeaFluxCO2.shelve'
        index = ('regionless', 'layerless', 'metricless')
    if field == 'SouthernTotalIceExtent':
        fn = prefix + j + '_SouthernTotalIceExtent.shelve'
        index = ('regionless', 'layerless', 'metricless')

    # S. Ocean mean: Air-Sea CO2 flux, surface DIC, SH Ice extent, surface nitrate, surface O2, surface salinity, SST
    if field == 'SO_AirSeaFlux':
        fn = prefix + j + '_AirSeaFluxCO2.shelve'
        index = ('SouthernOcean', 'Surface', 'mean')
    if field == 'SO_DIC':
        fn = prefix + j + '_DIC.shelve'
        index = ('SouthernOcean', 'Surface', 'mean')
    if field == 'SO_Nitrate':
        fn = prefix + j + '_Nitrate.shelve'
        index = ('SouthernOcean', 'Surface', 'mean')
    if field == 'SO_Oxygen':
        fn = prefix + j + '_Oxygen.shelve'
        index = ('SouthernOcean', 'Surface', 'mean')
    if field == 'SO_Salinity':
        fn = prefix + j + '_Salinity.shelve'
        index = ('SouthernOcean', 'Surface', 'mean')
    if field == 'SO_Temperature':
        fn = prefix + j + '_Temperature.shelve'
        index = ('SouthernOcean', 'Surface', 'mean')

    # S. Atlantic Ocean mean : Surface salinity, surface temp, volume integrated temperature
    if field == 'SA_Temperature':
        fn = prefix + j + '_Temperature.shelve'
        index = ('AtlanticSOcean', 'Surface', 'mean')
    if field == 'SA_Salinity':
        fn = prefix + j + '_Salinity.shelve'
        index = ('AtlanticSOcean', 'Surface', 'mean')
    if field == 'SA_VolumeMeanTemperature':
        fn = prefix + j + '_VolumeMeanTemperature.shelve'
        index = ('AtlanticSOcean', 'layerless', 'wcvweighted')
    if field == 'SA_MLD':
        fn = prefix + j + '_MLD.shelve'
        index = ('AtlanticSOcean', 'layerless', 'mean')

    # Weddel and Wilkes region: SST, SSS, Mixed layer depth, volume integrated temp.
    if field == 'Weddel_Temperature':
        fn = prefix + j + '_Temperature.shelve'
        index = ('Weddel', 'Surface', 'mean')
    if field == 'Weddel_Salinity':
        fn = prefix + j + '_Salinity.shelve'
        index = ('Weddel', 'Surface', 'mean')
    if field == 'Weddel_VolumeMeanTemperature':
        fn = prefix + j + '_VolumeMeanTemperature.shelve'
        index = ('Weddel', 'layerless', 'wcvweighted')
    if field == 'Weddel_MLD':
        fn = prefix + j + '_MLD.shelve'
        index = ('Weddel', 'layerless', 'mean')

    if field == 'Wilkes_Temperature':
        fn = prefix + j + '_Temperature.shelve'
        index = ('Wilkes', 'Surface', 'mean')
    if field == 'Wilkes_Salinity':
        fn = prefix + j + '_Salinity.shelve'
        index = ('Wilkes', 'Surface', 'mean')
    if field == 'Wilkes_VolumeMeanTemperature':
        fn = prefix + j + '_VolumeMeanTemperature.shelve'
        index = ('Wilkes', 'layerless', 'wcvweighted')
    if field == 'Wilkes_MLD':
        fn = prefix + j + '_MLD.shelve'
        index = ('Wilkes', 'layerless', 'mean')

    shelve = shopen(fn)
    data = shelve['modeldata'][index]
    shelve.close()
    return data


def fig1(field, range_key='mine', job='u-aw310'):

    data1 = getAMOCdata(job, field=field)

    times1 = sorted(data1.keys())
    data1 = [data1[t] for t in times1]

    newd1 = movingaverage_DT(data1, times1)
    pyplot.plot(times1, newd1, 'k', lw=0.5)

    newd2 = movingaverage_DT(data1,
                             times1,
                             window_len=30.,
                             window_units='years')
    pyplot.plot(times1, newd2, 'r', lw=1.5)

    title = ''
    if field in ['AMOC', 'Drake']:
        pyplot.ylabel('Sv')

    if field in [
            'GMT',
            'GVT',
    ] or field.find('Temperature') > -1:
        pyplot.ylabel('Celsius')

    if field in ['AirSeaFluxCO2', 'SO_AirSeaFlux']:
        pyplot.ylabel('Pg/yr')

    if field.find('Salinity') > -1:
        pyplot.ylabel('PSU')

    for text in [
            'DIC',
            'Nitrate',
            'Oxygen',
    ]:
        if field.find(text) > -1:
            pyplot.ylabel('mmol m^3')

    if field == 'AMOC':
        title += 'AMOC '

    if field[:3] == 'SO_':
        title += 'Southern Ocean '

    if field[:3] == 'SA_':
        title += 'Southern Atlantic Ocean '

    for text in [
            'Weddel',
            'Wilkes',
            'DIC',
            'Salinity',
            'Nitrate',
            'Oxygen',
    ]:
        if field.find(text) > -1:
            title += ' ' + text

    if field.find('VolumeMeanTemperature') > -1:
        title += ' Volume-Mean Temperature'

    elif field.find('Temperature') > -1:
        title += ' Temperature'

    if field.find('MLD') > -1:
        title += ' Mixed Layer Depth'
        pyplot.ylabel('m')
    if field == 'Drake':
        title = 'Drake Passage Current'
        if range_key == 'colins':
            pyplot.ylim([98., 208.])  # Colins suggestion
        if range_key == 'mine':
            pyplot.ylim([117., 193.])

    if field in ['AirSeaFluxCO2', 'SO_AirSeaFlux']:
        title = ' Air Sea Flux of CO2'

    if field == 'GVT':
        title = 'Global Volume Weighted Mean Temperature'

    if field == 'SouthernTotalIceExtent':
        title = 'Southern Hermisphere Total Ice Extent'
        pyplot.ylabel('1E6 km^2')

    if field == 'GMT':
        title = 'Global Mean Surface Temperature'

    if field == 'AirSeaFluxCO2':
        title = 'Air Sea Flux of CO2'

    title += ' - 10 year moving average'
    for i in range(5):
        title = title.replace('  ', ' ')
    pyplot.title(title)
    fn = 'custom_plots/' + job + '/' + field + '.png'
    pyplot.savefig(fn, dpi=300)
    print(fn)
    pyplot.close()


# Global : AMOC 26N, ACC Drake,  SST
fig1('AMOC')
fig1('Drake', range_key='mine')
#fig1('Drake',range_key='colins')
fig1('GMT')
fig1('GVT')
fig1('AirSeaFluxCO2')

# S. Ocean mean: Air-Sea CO2 flux, surface DIC, SH Ice extent, surface nitrate, surface O2, surface salinity, SST
fig1('SouthernTotalIceExtent')
fig1('SO_AirSeaFlux')
fig1('SO_DIC')
fig1('SO_Nitrate')
fig1('SO_Oxygen')
fig1('SO_Salinity')
fig1('SO_Temperature')

# S. Atlantic Ocean mean : Surface salinity, surface temp, volume integrated temperature
fig1('SA_Salinity')
fig1('SA_Temperature')
fig1('SA_VolumeMeanTemperature')
# fig1('SA_MLD')

# Weddel and Wilkes region: SST, SSS, Mixed layer depth, volume integrated temp.
fig1('Weddel_Salinity')
fig1('Weddel_Temperature')
fig1('Weddel_VolumeMeanTemperature')
fig1('Weddel_MLD')
fig1('Wilkes_Salinity')
fig1('Wilkes_Temperature')
fig1('Wilkes_VolumeMeanTemperature')
fig1('Wilkes_MLD')
