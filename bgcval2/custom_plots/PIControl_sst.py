#!/usr/bin/ipython

import os
import numpy as np
from matplotlib import pyplot
from shelve import open as shOpen
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition,
                                                  mark_inset,
                                                  zoomed_inset_axes)
"""
Plot for the CMIP6 implementation paper (section describing how we selected historical initial conditions)..

This page lists the dates of the historical ICS: https://code.metoffice.gov.uk/trac/ukcmip6/ticket/98 they span 1960 to 2815.

 
Global & Southern Ocean mean SST
6 yaer running mean (annual, runing)

LHS: 1960 to 2960 (renormalized so that 1960 = 1850)
RHS: 2550 to 2850

Mark Historical ICs

I'd suggest a plot covering 1960 to 2960 of u-aw310 (with all dates renormalized
so that 1960 = 1850). And plot both global and S. Ocean SST (with time meaning 
tested as above) and somehow include an indication of the years of the historical ICs.

Then in addition, using basically the same figure, include a blow-up of the period 
2550 to 2850 to highlight the 4 ICs for historical members r16-r19 as these 4 were 
deliberately chosen to sample the centennial timescale variability in the S.Ocean.

 
I am not sure whether we can get away with just using global SST or just S. Ocean
SST or both, so if we can look at both first and then decide, along with a 
suitable time meaning that would be great.

 
"""
job_dicts = {
    'r1i1p1f2': [
        'u-bc179',
        '2250',
    ],
    'r2i1p1f2': [
        'u-bc292',
        '2165',
    ],
    'r3i1p1f2': [
        'u-bc370',
        '2120',
    ],
    'r4i1p1f2': [
        'u-bb075',
        '1960',
    ],
    'r5i1p1f3': [
        'u-az513',
        '2020',
    ],
    'r6i1p1f3': [
        'u-az515',
        '2050',
    ],
    'r7i1p1f3': [
        'u-az524',
        '1995',
    ],
    'r8i1p1f2': [
        'u-bb277',
        '2395',
    ],
    'r9i1p1f2': [
        'u-bc470',
        '2285',
    ],
    'r10i1p1f2': [
        'u-bd288',
        '2340',
    ],
    'r11i1p1f2': [
        'u-bd416',
        '2460',
    ],
    'r12i1p1f2': [
        'u-bd483',
        '2200',
    ],
    'r13i1p1f2': [
        'u-bf935',
        '2565',
    ],
    'r14i1p1f2': [
        'u-bh100',
        '2685',
    ],
    'r15i1p1f2': [
        'u-bh101',
        '2745',
    ],
    'r16i1p1f2': [
        'u-bf647',
        '2629',
    ],
    'r17i1p1f2': [
        'u-bf656',
        '2716',
    ],
    'r18i1p1f2': [
        'u-bf703',
        '2760',
    ],
    'r19i1p1f2': [
        'u-bh162',
        '2815',
    ],
}


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


def get_data(j, field='AMOC'):
    index = ('regionless', 'layerless', 'metricless')

    if field == 'AMOC':
        fn = '/group_workspaces/jasmin2/ukesm/BGC_data/ldemora/shelves/timeseries/' + j + '/' + j + '_AMOC_26N.shelve'
    if field == 'Drake':
        fn = '/group_workspaces/jasmin2/ukesm/BGC_data/ldemora/shelves/timeseries/' + j + '/' + j + '_DrakePassageTransport.shelve'
    if field == 'GVT':
        fn = '/group_workspaces/jasmin2/ukesm/BGC_data/ldemora/shelves/timeseries/' + j + '/' + j + '_GlobalMeanTemperature.shelve'
    if field == 'SST':
        fn = '/group_workspaces/jasmin2/ukesm/BGC_data/ldemora/shelves/timeseries/' + j + '/' + j + '_Temperature.shelve'
        index = ('Global', 'Surface', 'mean')
    if field == 'SOMT':
        fn = '/group_workspaces/jasmin2/ukesm/BGC_data/ldemora/shelves/timeseries/' + j + '/' + j + '_Temperature.shelve'
        index = ('SouthernOcean', 'Surface', 'mean')
    if field == 'AirSeaFluxCO2':
        fn = '/group_workspaces/jasmin2/ukesm/BGC_data/ldemora/shelves/timeseries/' + j + '/' + j + '_TotalAirSeaFluxCO2.shelve'

    if not os.path.exists(fn):
        print('Does not exist:', fn)
    shelve = shOpen(fn)
    data = shelve['modeldata'][index]
    shelve.close()
    return data


def getClosestPoint(year, times, data):
    times = [int(t) for t in times]
    #print year, times.index(year), data[times.index(year)]
    return data[times.index(year)]


def fig1(field='SST', window_len=10):
    data1 = get_data(
        'u-aw310',
        field=field,
    )

    times1 = sorted(data1.keys())
    data1 = [data1[t] for t in times1]

    newd1 = movingaverage_DT(data1, times1, window_len=window_len)
    years = [yr[1] for cmip, yr in list(job_dicts.items())]

    pyplot.plot(times1, newd1, 'k', lw=1.5)

    for yr in years:
        value = getClosestPoint(int(yr), times1, newd1)
        pyplot.plot(float(yr), value, marker='o', ms=5, color='blue')

    pyplot.gca().set_xlim(1960., 2960.)
    if field == 'SST':
        pyplot.title('Sea Surface Temperature - ' + str(window_len) +
                     ' year moving average')
        pyplot.ylabel('Celsius')
    if field == 'SOMT':
        pyplot.title('Southern Ocean Surface Temperature - ' +
                     str(window_len) + ' year moving average')
        pyplot.ylabel('Celsius')

    pyplot.savefig(field + '_fig_' + str(window_len) + '.png', dpi=300)
    pyplot.close()


def fig2(field='SST', window_len=10):
    tdiff = -110
    data1 = get_data(
        'u-aw310',
        field=field,
    )

    times1 = np.array([int(t) for t in sorted(data1.keys())])
    data1 = [data1[t] for t in sorted(data1.keys())]

    if window_len:
        newd1 = movingaverage_DT(data1, times1, window_len=window_len)
        title_suff = ' - ' + str(window_len) + ' year moving average'
    else:
        newd1 = data1[:]
        title_suff = ''
#years = [yr[1] for cmip, yr in job_dicts.items()]
#runid = [yr[0] for cmip, yr in job_dicts.items()]

    fig, ax1 = pyplot.subplots()

    pyplot.plot(times1 + tdiff, newd1, 'k', lw=1.5)

    for run, [jobid, yr] in list(job_dicts.items()):
        print(yr, run, run[:3])
        value = getClosestPoint(int(yr), times1, newd1)
        if run[:3] in ['r16', 'r17', 'r18', 'r19']:
            color = 'red'
        else:
            color = 'blue'
        pyplot.plot(float(yr) + tdiff, value, marker='o', ms=5, color=color)

    pyplot.gca().set_xlim(1955. + tdiff, 2960. + tdiff)

    if field == 'SST':
        pyplot.gca().set_ylim(17.0, 18.1)

        pyplot.title('Sea Surface Temperature' + title_suff)
        pyplot.ylabel('Celsius')
    elif field == 'SOMT':
        pyplot.gca().set_ylim(4.2, 5.5)
        pyplot.title('Southern Ocean Surface Temperature' + title_suff)
        pyplot.ylabel('Celsius')
    elif field == 'Drake':
        ylims = pyplot.gca().get_ylim()
        ylims = [ylims[0] - 0.45 * (ylims[1] - ylims[0]), ylims[1]]
        pyplot.gca().set_ylim(ylims)
        pyplot.title('Drake Passage Current' + title_suff)
        pyplot.ylabel('Sv')
    else:
        ylims = pyplot.gca().get_ylim()
        ylims = [ylims[0] - 0.45 * (ylims[1] - ylims[0]), ylims[1]]
        pyplot.gca().set_ylim(ylims)
        pyplot.title(field + title_suff)

# second axes:
    # Create a set of inset Axes: these should fill the bounding box allocated to
    # them.
    ax2 = pyplot.axes([0, 0, 1, 1])
    # Manually set the position and relative size of the inset axes within ax1
    ip = InsetPosition(ax1, [0.1, 0.075, 0.75, 0.25])
    ax2.set_axes_locator(ip)
    # Mark the region corresponding to the inset axes on ax1 and draw lines
    # in grey linking the two axes.
    mark_inset(ax1, ax2, loc1=2, loc2=4, fc="none", ec='0.5')
    ax2.plot(times1 + tdiff, newd1, 'k', lw=1.5)
    ax2.set_xlim(2550. + tdiff, 2850. + tdiff)
    if field == 'SOMT':
        pyplot.gca().set_ylim(4.65, 5.45)

    if field == 'Drake':
        pyplot.gca().set_ylim(138., 172.)

    for run, [jobid, yr] in list(job_dicts.items()):
        value = getClosestPoint(int(yr), times1, newd1)
        if run[:3] in ['r16', 'r17', 'r18', 'r19']:
            color = 'red'
        else:
            color = 'blue'
        if int(yr) < 2550: continue
        if int(yr) > 2850: continue
        value = getClosestPoint(int(yr), times1, newd1)
        ax2.plot(float(yr) + tdiff, value, marker='o', ms=5, color=color)
    ax2.tick_params(
        axis='x',
        labelsize=7)  #direction='out', length=6, width=2, colors='r',
    ax2.tick_params(
        axis='y',
        labelsize=7)  #direction='out', length=6, width=2, colors='r',

    #grid_color='r', grid_alpha=0.5)

    fn = field + '_inset_fig_' + str(window_len) + '.png'
    print("Saving:", fn)
    pyplot.savefig(fn, dpi=300)
    pyplot.close()


def fig5(field='SST', window_len=10):
    tdiff = 0  #-110

    data1 = get_data(
        'u-aw310',
        field=field,
    )
    #data1 = {y:5 for i,y in enumerate(range(1750,3000))}
    times1 = np.array([int(t) for t in sorted(data1.keys())])
    data1 = [data1[t] for t in sorted(data1.keys())]

    if window_len:
        newd1 = movingaverage_DT(data1, times1, window_len=window_len)
        title_suff = ' - ' + str(window_len) + ' year moving average'
    else:
        newd1 = data1[:]
        title_suff = ''

    fig, ax1 = pyplot.subplots()

    pyplot.plot(times1 + tdiff, newd1, 'k', lw=1.5)

    for run, [jobid, yr] in list(job_dicts.items()):
        print(yr, run, run[:3])
        value = getClosestPoint(int(yr), times1, newd1)
        if run[:3] in ['r16', 'r17', 'r18', 'r19']:
            color = 'red'
        else:
            color = 'blue'
        pyplot.plot(float(yr) + tdiff, value, marker='o', ms=5, color=color)

    pyplot.gca().set_xlim(1955. + tdiff, 2960. + tdiff)

    if field == 'SST':
        pyplot.gca().set_ylim(17.0, 18.1)

        pyplot.title('Sea Surface Temperature' + title_suff)
        pyplot.ylabel('Celsius')
    elif field == 'SOMT':
        pyplot.gca().set_ylim(4.2, 5.5)
        pyplot.title('Southern Ocean Surface Temperature' + title_suff)
        pyplot.ylabel('Celsius')
    else:
        ylims = pyplot.gca().get_ylim()
        ylims = [ylims[0] - 0.45 * (ylims[1] - ylims[0]), ylims[1]]
        pyplot.gca().set_ylim(ylims)
        pyplot.title(field + title_suff)

    axespos = [0.12, 0.07, 0.75, 0.3]
    # Add 3rd set of axes
    ax3 = pyplot.axes([0, 0, 1, 1], label='sst')
    ax3.patch.set_alpha(0.7)
    ip3 = InsetPosition(ax1, axespos)
    ax3.set_axes_locator(ip3)
    ax3.set_ylim(17.39, 18.1)

    #data2 = {y:17.1 for i,y in enumerate(range(1750,3000))}
    data2 = get_data(
        'u-aw310',
        field='SST',
    )
    times2 = np.array([int(t) for t in sorted(data2.keys())])
    data2 = np.array([data2[t] for t in sorted(data2.keys())])
    #data2[300:] = 16.9
    color = 'tab:green'

    ax3.plot(times2 + tdiff, data2, color=color, zorder=3, lw=0.5)
    ax3.set_xlim(2550. + tdiff, 2850. + tdiff)

    ax3.set_ylabel(
        'SST',
        color=color,
    )  # fontsize = 8)
    ax3.yaxis.set_label_position("right")
    ax3.yaxis.tick_right()
    ax3.tick_params(
        axis='x',
        labelsize=7)  #direction='out', length=6, width=2, colors='r',
    ax3.tick_params(axis='y', labelcolor=color, labelsize=7, colors=color)

    # second axes:
    # Create a set of inset Axes: these should fill the bounding box allocated to
    # them.
    ax2 = pyplot.axes([0, 0, 1, 1], label='somt')
    ip = InsetPosition(ax1, axespos)
    ax2.set_axes_locator(ip)
    ax2.patch.set_alpha(0.3)
    # Mark the region corresponding to the inset axes on ax1 and draw lines
    # in grey linking the two axes.
    mark_inset(ax1, ax2, loc1=2, loc2=4, fc="none", ec='0.5', zorder=1)

    ax2.plot(times1 + tdiff, newd1, 'k', lw=1.5, zorder=5)
    ax2.set_xlim(2550. + tdiff, 2850. + tdiff)
    if field == 'SOMT':
        pyplot.gca().set_ylim(4.65, 5.45)

    for run, [jobid, yr] in list(job_dicts.items()):
        value = getClosestPoint(int(yr), times1, newd1)
        if run[:3] in ['r16', 'r17', 'r18', 'r19']:
            color = 'red'
        else:
            color = 'blue'
        if int(yr) < 2550: continue
        if int(yr) > 2850: continue
        value = getClosestPoint(int(yr), times1, newd1)
        ax2.plot(float(yr) + tdiff,
                 value,
                 marker='o',
                 ms=6,
                 color=color,
                 zorder=6)

    ax2.set_ylabel(
        'SOST',
        color='k',
    )  #fontsize = 8)
    ax2.yaxis.set_label_position("left")
    ax2.tick_params(
        axis='x',
        labelsize=7)  #direction='out', length=6, width=2, colors='r',
    ax2.tick_params(
        axis='y', labelsize=7, labelcolor='k', right=False,
        left=True)  #direction='out', length=6, width=2, colors='r',

    fn = field + '_inset_fig5_' + str(window_len) + '.png'
    print("Saving:", fn)
    pyplot.savefig(fn, dpi=300)
    pyplot.close()


#fig2(field = 'Drake', window_len=0)
fig5(field='SOMT', window_len=0)
#fig4(field = 'SOMT', window_len=0)

#for field in ['SST', 'SOMT', 'AMOC', 'Drake', 'GVT', 'AirSeaFluxCO2']:
#    for window_len in [0, 1,3,5,10]:
#        fig2(field = field, window_len=window_len)
#	fig1(field = field, window_len=window_len)
