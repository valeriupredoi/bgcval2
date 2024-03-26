#!/usr/bin/ipython

# 1. (Easy) a time series of only UKESM1 u-aw310 for the full length (~1960 to year 3???) applying either a 5 year or 10 year running mean (likely the former).

# 2. A time series (using the same time meaning as in 1...so likely 5 year running mean) that consists of the following:
#    (i) A 250 year segment of the piControl (it probably does not matter what 250 years are chosen).
#    (ii) The ensemble mean (of however many historical runs have completed...think it is now 11), 5 year running mean, of the 165 years of the (11) historical runs
#    (iii) A 5-member ensemble mean for each of the 4 Tier 1 scenarios for the 85 years of the projections.

import os
import numpy as np
from matplotlib import pyplot
from shelve import open as shOpen


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
    index = ('regionless', 'layerless', 'metricless')

    if field == 'AMOC':
        fn = '/group_workspaces/jasmin2/ukesm/BGC_data/ldemora/shelves/timeseries/' + j + '/' + j + '_AMOC_26N.shelve'
    if field == 'Drake':
        fn = '/group_workspaces/jasmin2/ukesm/BGC_data/ldemora/shelves/timeseries/' + j + '/' + j + '_DrakePassageTransport.shelve'
    if field == 'GVT':
        fn = '/group_workspaces/jasmin2/ukesm/BGC_data/ldemora/shelves/timeseries/' + j + '/' + j + '_GlobalMeanTemperature.shelve'
    if field == 'GMT':
        fn = '/group_workspaces/jasmin2/ukesm/BGC_data/ldemora/shelves/timeseries/' + j + '/' + j + '_Temperature.shelve'
        index = ('Global', 'Surface', 'mean')
    if field == 'GMT_anomaly':
        fn = '/group_workspaces/jasmin2/ukesm/BGC_data/ldemora/shelves/timeseries/' + j + '/' + j + '_Temperature.shelve'
        index = ('Global', 'Surface', 'mean')

    if field == 'AirSeaFluxCO2':
        fn = '/group_workspaces/jasmin2/ukesm/BGC_data/ldemora/shelves/timeseries/' + j + '/' + j + '_TotalAirSeaFluxCO2.shelve'

    if not os.path.exists(fn):
        print('Does not exist:', fn)
    shelve = shOpen(fn)
    data = shelve['modeldata'][index]
    shelve.close()
    #	times = sorted(data.keys())
    #	data = [data[t] for t in times]
    return data


def add_co2_data():
    co2fn = 'custom_plots/uptake_data_for_GFN_2014.txt'
    co2f = open(co2fn)
    co2datav = co2f.readlines()
    co2f.close()
    co2_time, co2uptake, co2_min, co2_max = [], [], [], []

    for line in co2datav[1:]:
        line = line.replace('\r\n', '')
        line = [float(l) for l in line.split()]
        print(line)
        co2_time.append(line[0])
        co2uptake.append(line[1])
        co2_min.append(line[2])
        co2_max.append(line[3])
    color = 'black'
    pyplot.fill_between(co2_time, co2_min, co2_max, color=color, alpha=0.15)
    pyplot.plot(co2_time,
                co2uptake,
                color,
                lw=1.2,
                ls='--',
                label='Khatiwala 2014')


def fig1():
    data1 = getAMOCdata('u-aw310')

    times1 = sorted(data1.keys())
    data1 = [data1[t] for t in times1]

    #pyplot.plot(times, amoc,'k',lw=0.3)
    newd1 = movingaverage_DT(data1, times1)
    pyplot.plot(times1, newd1, 'k', lw=1.5)
    pyplot.title('AMOC - 10 year moving average')
    pyplot.ylabel('Sv')
    pyplot.savefig('custom_plots/amoc_fig1.png', dpi=300)
    pyplot.close()


#fig1()


def fig2(field, range_key='mine'):
    piControl = [
        'u-aw310',
    ]
    historical = [
        'u-az513', 'u-az515', 'u-az524', 'u-bb075', 'u-bb277', 'u-bc179',
        'u-bc292', 'u-bc370', 'u-bc470'
    ]
    ssp126 = ['u-be509', 'u-be679', 'u-be682', 'u-be393', 'u-be397']
    ssp245 = ['u-be537', 'u-be606', 'u-be683', 'u-be394', 'u-be398']
    ssp370 = ['u-be647', 'u-be690', 'u-be684', 'u-be335', 'u-be395']
    ssp585 = ['u-be653', 'u-be693', 'u-be686', 'u-be392', 'u-be396']

    # Tier2:
    ssp119 = ['u-bh409', 'u-bh570', 'u-bh716', 'u-bh210', 'u-bh807']
    ssp434 = ['u-bh454', 'u-bh724', 'u-bh717', 'u-bh254', 'u-bh808']
    ssp534 = ['u-bh456', 'u-bh744', 'u-bh718', 'u-bh285', 'u-bh809']

    ssp_jobs = {
        'piControl': piControl,
        'historical': historical,
        'SSP 1 2.6': ssp126,
        'SSP 2 4.5': ssp245,
        'SSP 3 7.0': ssp370,
        'SSP 5 8.5': ssp585,
        'SSP 1 1.9': ssp119,
        'SSP 4 3.4': ssp434,
        'SSP 5 3.4': ssp534,
    }
    #	ssp_colours = {'piControl': 'black', 'historical': 'purple', 'SSP 1 2.6': 'blue', 'SSP 2 4.5': 'green', 'SSP 3 7.0': 'pink','SSP 5 8.5': 'red',}

    order = [
        'piControl', 'historical', 'SSP 1 1.9', 'SSP 1 2.6', 'SSP 2 4.5',
        'SSP 3 7.0', 'SSP 4 3.4', 'SSP 5 3.4', 'SSP 5 8.5'
    ]

    #        ssp_colours = {'piControl': 'black', 'historical': 'purple', 'SSP 1 2.6': (0., 52./256, 102./256), 'SSP 2 4.5': (112./256, 160./256, 205./256), 'SSP 3 7.0': (196./256, 121./256, 0.),'SSP 5 8.5': (153./256, 0., 2./256),}

    ssp_colours = {
        'piControl': 'black',
        'historical': 'purple',
        'SSP 1 2.6': '#003466',
        'SSP 2 4.5': '#70a0cd',
        'SSP 3 7.0': '#c47900',
        'SSP 5 8.5': '#990002',
        'SSP 1 1.9': 'green',
        'SSP 4 3.4': 'gold',
        'SSP 5 3.4': 'orange',
    }
    # Older colours:
    #        ssp_colours = {'piControl': 'green', 'historical': 'purple', 'SSP 1 2.6': '#003466', 'SSP 2 4.5': '#70a0cd', 'SSP 3 7.0': '#c47900','SSP 5 8.5': '#990002',		}
    if field in [
            'GMT_anomaly',
    ]:
        datas = getAMOCdata(ssp_jobs['piControl'][0], field=field)
        times = np.array(sorted(datas.keys()))
        data = [datas[t] for t in times]
        anomaly = np.mean(data)
# Fill between. pi control:
    for ssp in order:
        if ssp != 'piControl': continue
        jobs = ssp_jobs[ssp]
        runs = {}
        for job in jobs:
            datas = getAMOCdata(job, field=field)
            times = np.array(sorted(datas.keys()))
            data = [datas[t] for t in times]
            years = [
                '1960', '1995', '2020', '2050', '2120', '2165', '2200', '2250',
                '2285', '2340', '2395', '2460', '2619', '2716', '2760', '2815'
            ]

            if field in [
                    'GMT_anomaly',
            ]: data = np.array(data) - anomaly
            for start_times in years:
                jobname = job + '_' + str(start_times)
                new_times = times - float(start_times) + 1601 - 250.

                new_times = np.ma.masked_where(
                    (new_times < 1750.) + (new_times > 1850.), new_times)
                new_data = np.ma.masked_where(new_times.mask, data)
                new_times = new_times.compressed()
                new_data = new_data.compressed()

                #new_data = movingaverage_DT(new_data, new_times)
                runs[jobname] = {t: d for t, d in zip(new_times, new_data)}

        times, data_min = ensemblemean(runs, operator='min')
        times, data_max = ensemblemean(runs, operator='max')
        times, data_mean = ensemblemean(runs, operator='mean')

        data_min = movingaverage_DT(data_min, times)
        data_max = movingaverage_DT(data_max, times)
        data_mean = movingaverage_DT(data_mean, times)

        pyplot.plot(times, data_mean, ssp_colours[ssp], lw=1.5, label=ssp)
        pyplot.fill_between(times,
                            data_min,
                            data_max,
                            color=ssp_colours[ssp],
                            alpha=0.2)

# Ensemble mean lines
    for ssp in order:
        if ssp == 'piControl': continue
        jobs = ssp_jobs[ssp]
        runs = {}
        for job in jobs:
            runs[job] = getAMOCdata(job, field=field)
        times, data = ensemblemean(runs)
        data = movingaverage_DT(
            data,
            times,
        )
        if field in [
                'GMT_anomaly',
        ]: data = np.array(data) - anomaly
        pyplot.plot(times, data, ssp_colours[ssp], lw=1.5, label=ssp)

    # Fill between. not pi control:
    for ssp in order:
        if ssp == 'piControl': continue
        jobs = ssp_jobs[ssp]
        runs = {}
        for job in jobs:
            datas = getAMOCdata(job, field=field)
            times = sorted(datas.keys())
            data = [datas[t] for t in times]
            if field in [
                    'GMT_anomaly',
            ]: data = np.array(data) - anomaly
            #data = movingaverage_DT(data, times,)
            runs[job] = {t: d for t, d in zip(times, data)}

        times, data_min = ensemblemean(runs, operator='min')
        times, data_max = ensemblemean(runs, operator='max')
        #		for t, mi, ma in zip(times, data_min, data_max):
        #			print ssp, t, mi, ma
        data_min = movingaverage_DT(
            data_min,
            times,
        )
        data_max = movingaverage_DT(
            data_max,
            times,
        )
        pyplot.fill_between(times,
                            data_min,
                            data_max,
                            color=ssp_colours[ssp],
                            alpha=0.2)

    if field == 'AMOC':
        pyplot.title('AMOC - 10 year moving average')
        pyplot.ylabel('Sv')
        pyplot.legend(loc='lower left')

    if field == 'Drake':
        pyplot.title('Drake Passage Current - 10 year moving average')
        if range_key == 'colins':
            pyplot.ylim([98., 208.])  # Colins suggestion
        if range_key == 'mine':
            pyplot.ylim([117., 193.])
        pyplot.legend(loc='upper left')
        pyplot.ylabel('Sv')

    if field == 'GVT':
        pyplot.title(
            'Global Volume weighted mean Temperature - 10 year moving average')
        pyplot.ylabel('Celsius')
        pyplot.legend(loc='upper left')

    if field == 'GMT':
        pyplot.title(
            'Global Mean Surface Temperature - 10 year moving average')
        pyplot.ylabel('Celsius')
        pyplot.legend(loc='upper left')
    if field == 'GMT_anomaly':
        pyplot.title('Global Mean Surface Temperature anomaly')
        pyplot.ylabel('Celsius')
        pyplot.legend(loc='upper left')

    if field == 'AirSeaFluxCO2':
        pyplot.title('Air Sea Flux of CO2 - 10 year moving average')
        pyplot.ylabel('Pg/yr')
        add_co2_data()
        pyplot.legend(loc='upper left')

    pyplot.gca().tick_params(axis='y', left=True, right=True)
    fn = 'custom_plots/' + field + '_fig2.png'
    if range_key == 'colins':
        fn = 'custom_plots/' + field + '_fig2_colins_range.png'
    pyplot.savefig(fn, dpi=600)
    pyplot.close()

fig2('GMT_anomaly')

fig2('GMT')
fig2('AirSeaFluxCO2')
fig2('AMOC')
fig2('Drake', range_key='mine')
fig2('Drake', range_key='colins')
fig2('GVT')
