#!/usr/bin/ipython

#
# Copyright 2017, Plymouth Marine Laboratory
#
# This file is part of the bgc-val library.
#
# bgc-val is free software: you can redistribute it and/or modify it
# under the terms of the Revised Berkeley Software Distribution (BSD) 3-clause license.

# bgc-val is distributed in the hope that it will be useful, but
# without any warranty; without even the implied warranty of merchantability
# or fitness for a particular purpose. See the revised BSD license for more details.
# You should have received a copy of the revised BSD license along with bgc-val.
# If not, see <http://opensource.org/licenses/BSD-3-Clause>.
#
# Address:
# Plymouth Marine Laboratory
# Prospect Place, The Hoe
# Plymouth, PL1 3DH, UK
#
# Email:
# ledm@pml.ac.uk
"""
.. module:: extentProfiles
   :platform: Unix
   :synopsis: A tool for producing a depth profile/Transect of contours.
.. moduleauthor:: Lee de Mora <ledm@pml.ac.uk>
"""

import numpy as np
from shelve import open as shOpen
from netCDF4 import num2date
import os
from glob import glob
import shutil
from matplotlib import pyplot, gridspec
from matplotlib.colors import LogNorm
import matplotlib.patches as mpatches
from matplotlib.ticker import Locator
import cartopy
import cartopy.crs as ccrs
import cartopy.io.shapereader as shapereader
from cartopy import img_transform, feature as cfeature

#Specific local code:
from .. import UKESMpython as ukp
from ..bgcvaltools.pftnames import getLongName
from ..bgcvaltools.dataset import dataset
from . import timeseriesTools as tst
from . import timeseriesPlots as tsp
from ..Paths import paths

try:
    defcmap = pyplot.cm.jet
    defcmapstr = 'jet'
except:
    defcmap = viridis
    defcmapstr = 'viridis'

zonalCuts = [
    'Equator',
    '10 N',
    '10 S',
]
MeridionalCuts = ['Atlantic', 'Atlantic28W', 'Pacific135W']


class MinorSymLogLocator(Locator):
    """
    Dynamically find minor tick positions based on the positions of
    major ticks for a symlog scaling.
    from: https://stackoverflow.com/questions/20470892/how-to-place-minor-ticks-on-symlog-scale
    """

    def __init__(self, linthresh):
        """
        Ticks will be placed between the major ticks.
        The placement is linear for x between -linthresh and linthresh,
        otherwise its logarithmically
        """
        self.linthresh = linthresh

    def __call__(self):
        'Return the locations of the ticks'
        majorlocs = self.axis.get_majorticklocs()

        # iterate through minor locs
        minorlocs = []

        # handle the lowest part
        for i in range(1, len(majorlocs)):
            majorstep = majorlocs[i] - majorlocs[i - 1]
            if abs(majorlocs[i - 1] + majorstep / 2) < self.linthresh:
                ndivs = 10
            else:
                ndivs = 9
            minorstep = majorstep / ndivs
            locs = np.arange(majorlocs[i - 1], majorlocs[i], minorstep)[1:]
            minorlocs.extend(locs)

        return self.raise_if_exceeds(np.array(minorlocs))

    def tick_values(self, vmin, vmax):
        raise NotImplementedError('Cannot get tick locations for a '
                                  '%s type.' % type(self))


def makeLonSafe(lon):
    while True:
        if 0. < lon <= 360.: return lon
        if lon <= 0.: lon += 360.
        if lon > 360.: lon -= 360.


def makeLonSafeArr(lon):
    if lon.ndim == 2:
        for l, lon1 in enumerate(lon):
            for ll, lon2 in enumerate(lon1):
                lon[l, ll] = makeLonSafe(lon2)
        return lon
    if lon.ndim == 1:
        for l, lon1 in enumerate(lon):
            lon[l] = makeLonSafe(lon1)
        return lon


def findClosest(arr, z, debug=False):
    """
		Calculate the index of the closest point on the array, arr and the point z.
	"""
    arr = np.array(arr)
    if arr.ndim > 1:
        print("Array should be one dimenstional!", arr.shape)
        assert 0
    d = 10000.
    best = -1
    for i, zz in enumerate(arr.squeeze()):
        #print i,z,zz,depth.shape
        d2 = abs(z - zz)
        if d2 < d:
            d = d2
            best = i
    if debug:
        print('Closest: in situ:', z, 'index:', best, 'distance:', d, ', closest model:', arr.shape, arr[
            best])
    return best


def loadKeyFromFile(
    fn,
    coords,
    nc='',
):
    if nc == '':
        nc = dataset(fn, 'r')
    dtimes = num2date(nc.variables[coords['t']][:],
                      nc.variables[coords['t']].units,
                      calendar=coords['cal'])[:]
    return str(dtimes[0].year)
    #return os.path.basename(fn).replace('u-ad371o_1y_','').replace('_ptrc_T.nc','')[:4]


def LoadZonalTransect(
    data,
    transectLat,
    lats,
    lons,
    depths,
):
    #####
    # Loads data along a horizontal transect (ie equator, 10 N....)
    data = np.ma.array(data)
    lons = makeLonSafeArr(lons)
    depths = depths[::-1] * -1.

    if lats.ndim == 1:
        loc = findClosest(lats, transectLat)
        newlons = lons
        if data.ndim == 4: dat = data[:, ::-1, loc, :].squeeze()

    if lats.ndim == 2:
        print("LoadZonalTransect:", transectLat, data.shape, lats.shape, lons.shape, depths.shape)
        dat = np.ma.zeros((data.shape[1], data.shape[3]))
        newlons = np.zeros(data.shape[3])
        #print data.shape, outDat.shape
        #assert 0
        for i in np.arange(lats.shape[1]):
            loc = findClosest(lats[:, i], transectLat, debug=False)
            #print i, loc, lats[:,i].mean(), lats[loc,i],lons[loc,i]
            dat[:, i] = data[0, ::-1, loc, i]
            newlons[i] = lons[loc, i]

    newX, newZ = np.meshgrid(newlons, depths)
    dat = np.ma.array(dat)
    dat = np.ma.masked_where((newX >= 359.2) + (newX < 0.2) + dat.mask, dat)
    dat = np.ma.masked_where(
        np.ma.array(dat).mask + (dat < 1E-10) + (dat > 1e10), dat).squeeze()

    return dat, newX, newZ


def LoadMeridionalTransect(data, transectLon, lats, lons, depths):
    #####
    # Loads data along a vertical transect (ie 0E, dateline, etc)
    data = np.ma.array(data)
    lons = makeLonSafeArr(lons)
    transectLon = makeLonSafe(transectLon)
    depths = depths[::-1] * -1.

    if lons.ndim == 1:
        loc = findClosest(lons, transectLon)
        newlats = lats
        if data.ndim == 4: dat = data[:, ::-1, :, loc].squeeze()
        print("LoadMeridionalTransect:", transectLon, data.shape, lons.shape, loc, lons[
            loc])

    if lons.ndim == 2:
        ####
        # For each line of the 2D grid, it finds the closest point in longitude.
        print("LoadMeridionalTransect:", transectLon, data.shape, lons.shape, lons.shape, depths.shape)
        dat = np.ma.zeros((data.shape[1], data.shape[2]))
        newlats = np.zeros(data.shape[2])
        print(data.shape, dat.shape)
        for i in np.arange(lons.shape[0]):
            loc = findClosest(lons[i, :], transectLon, debug=False)

            print("LoadMeridionalTransect:", i, transectLon, loc, [
                lats[i, loc],
                'N',
                lons[i, loc],
                'E',
            ])
            if transectLon - lons[i, loc] > 5.:
                dat[:, i] = np.ma.masked_all_like(depths)
                newlats[i] = lats[i, loc]
            else:
                dat[:, i] = data[0, ::-1, i, loc]
                newlats[i] = lats[i, loc]

    newX, newZ = np.meshgrid(newlats, depths)
    dat = np.ma.array(dat)
    print(dat.shape, newX.shape, newZ.shape)
    dat = np.ma.masked_where((newX >= 90.) + (newX <= -90.) + dat.mask, dat)
    dat = np.ma.masked_where(
        np.ma.array(dat).mask + (dat < 1E-10) + (dat > 1e10), dat).squeeze()
    return dat, newX, newZ


def contourplot(
    jobID,
    name,
    modelfiles,
    datafile,
    modeldetails,
    datadetails,
    modelcoords,
    datacoords,
    zmin=0.,
    zmax=400.,
    oxcutoff=80.,
    plotKey='Equator',
    title='',
    filename='',
    cbarlabel='',
):

    contours = [zmin, oxcutoff, zmax]
    plotKeyDict = {
        'Equator': 0.,
        '10 N': 10.,
        '10 S': -10.,
        'Atlantic28W': -28.,
        'Pacific135W': -135.
    }
    zonalCuts = [
        'Equator',
        '10 N',
        '10 S',
    ]
    MeridionalCuts = ['Atlantic', 'Atlantic28W', 'Pacific135W']

    #####
    # Load plot details
    pd = {}
    pd['Data'] = {
        'label': 'Data',
        'c': [
            'k',
        ],
        'lw': [
            2.,
        ],
        'ls': [
            '-',
        ]
    }
    for i, fn in enumerate(modelfiles):
        lw = 1
        key = loadKeyFromFile(
            fn,
            modelcoords,
            nc='',
        )
        color = defcmap(float(i) / float(len(modelfiles)))
        label = key
        pd[key] = {
            'label': label,
            'c': [
                color,
            ],
            'lw': [
                lw,
            ],
            'ls': [
                '-',
            ],
        }

    #####
    # Load data lats/lons, and data
    print("Loading:", datafile)
    dnc = dataset(datafile, 'r')
    dlats = dnc.variables[datacoords['lat']][:]
    dlons = dnc.variables[datacoords['lon']][:]
    ddepths = dnc.variables[datacoords['z']][:]
    do2 = ukp.extractData(dnc, datadetails)
    dnc.close()

    if plotKey in zonalCuts:
        do2, dnewX, dnewZ = LoadZonalTransect(do2, plotKeyDict[plotKey], dlats,
                                              dlons, ddepths)

    if plotKey in MeridionalCuts:
        do2, dnewX, dnewZ = LoadMeridionalTransect(do2, plotKeyDict[plotKey],
                                                   dlats, dlons, ddepths)

    #####
    # Load model lats/lons
    print("Loading:", modelfiles[0])
    nc = dataset(modelfiles[0], 'r')
    lats = nc.variables[modelcoords['lat']][:]
    lons = nc.variables[modelcoords['lon']][:]
    depths = nc.variables[modelcoords['z']][:]
    nc.close()

    #####
    # Start making the figure
    fig = pyplot.figure()
    fig.set_size_inches(10, 6)
    ax = pyplot.subplot(111)

    #####
    # Add data as a colormesh
    cmesh = pyplot.pcolormesh(dnewX,
                              dnewZ,
                              np.ma.masked_where(do2.mask, do2),
                              cmap='Blues_r',
                              vmin=zmin,
                              vmax=zmax)
    im = ax.contour(
        dnewX,
        dnewZ,
        do2,
        contours,
        colors=pd['Data']['c'],
        linewidths=pd['Data']['lw'],
        linestyles=pd['Data']['ls'],
    )

    #####
    # Add model data as a colormesh
    for fn in modelfiles:
        nc = dataset(fn, 'r')
        o2 = ukp.extractData(nc, modeldetails)
        key = loadKeyFromFile(
            fn,
            modelcoords,
            nc=nc,
        )
        nc.close()

        if plotKey in zonalCuts:
            o2, newX, newZ = LoadZonalTransect(o2, plotKeyDict[plotKey], lats,
                                               lons, depths)

        if plotKey in MeridionalCuts:
            o2, newX, newZ = LoadMeridionalTransect(o2, plotKeyDict[plotKey],
                                                    lats, lons, depths)

        im = ax.contour(
            newX,
            newZ,
            o2,
            contours,
            colors=pd[key]['c'],
            linewidths=pd[key]['lw'],
            linestyles=pd[key]['ls'],
        )

    #####
    # Add legend:
    addLegend = True
    if addLegend:
        legendSize = len(list(pd.keys())) + 1
        ncols = int(legendSize / 25) + 1
        box = ax.get_position()
        ax.set_position(
            [box.x0, box.y0, box.width * (1. - 0.1 * ncols), box.height])

        for i in sorted(pd.keys()):
            pyplot.plot(
                [],
                [],
                c=pd[i]['c'][0],
                lw=pd[i]['lw'][0],
                ls=pd[i]['ls'][0],
                label=pd[i]['label'],
            )

        legd = ax.legend(loc='center left',
                         ncol=ncols,
                         prop={'size': 10},
                         bbox_to_anchor=(1., 0.5))
        legd.draw_frame(False)
        legd.get_frame().set_alpha(0.)

    c1 = pyplot.colorbar(
        cmesh,
        orientation='horizontal',
    )  #pad=0.05,shrink=0.9 )
    if len(cbarlabel) > 0: c1.set_label(cbarlabel)

    ####
    # Sort out X axis:
    if plotKey in zonalCuts:
        pyplot.xticks([0., 60., 120., 180., 240., 300., 360.])
        ax.set_xlim([0., 360.])
        pyplot.xlabel('Longitude, degress E', )
    if plotKey in MeridionalCuts:
        pyplot.xlabel('Latitude, degress N', )
        ax.set_xlim([-90., 90.])
        pyplot.xticks([-90., -60., -30., 0., 30., 60., 90.])
    ####
    # Sort out Y axis:
    pyplot.axhline(y=-500., c='k', ls='--')
    pyplot.axhline(y=-1000., c='k', ls='--')
    ax.set_yscale('symlog', linthreshy=1000.)
    ax.set_ylim([min([dnewZ.min(), newZ.min()]), -1.])
    pyplot.yticks([-10., -100., -500., -1000., -2000., -5000.],
                  ['10', '100', '500', '1000', '2000', '5000'])
    pyplot.ylabel('Depth, m', ha='center', va='center', rotation='vertical')

    pyplot.title(title)

    #####
    #Save figure
    if filename == '':
        filename = ukp.folder(['images', jobID, 'OMZ']) + '-'.join(
            [name, plotKey, 'contour',
             str(int(oxcutoff))]) + '.png'
    print("saving", filename)
    pyplot.savefig(filename)
    pyplot.close()


def run():

    def listModelDataFiles(jobID, filekey, datafolder, annual):
        if annual:
            return sorted(
                glob(datafolder + jobID + "/" + jobID + "o_1y_*_" + filekey +
                     ".nc"))
        else:
            return sorted(
                glob(datafolder + jobID + "/" + jobID + "o_1m_*_" + filekey +
                     ".nc"))

    jobID = 'u-ai886'
    modelfiles = listModelDataFiles(jobID, 'ptrc_T', paths.ModelFolder_pref,
                                    True)
    #glob('/data/euryale7/scratch/ledm/UKESM/MEDUSA/'+jobID+'/'+jobID+'o_1y_*_ptrc_T.nc')
    datafile = paths.WOAFolder_annual + 'woa13_all_o00_01.nc'
    name = 'OMZExtent'

    medusaCoords = {
        't': 'time_counter',
        'z': 'deptht',
        'lat': 'nav_lat',
        'lon': 'nav_lon',
        'cal': '360_day',
    }  # model doesn't need time dict.
    woaCoords = {
        't': 'index_t',
        'z': 'depth',
        'lat': 'lat',
        'lon': 'lon',
        'cal': 'standard',
        'tdict': ukp.tdicts['ZeroToZero']
    }
    modeldetails = {
        'name': name,
        'vars': [
            'OXY',
        ],
        'convert': ukp.NoChange,
        'units': 'mmol O2/m^3'
    }
    datadetails = {
        'name': name,
        'vars': [
            'o_an',
        ],
        'convert': ukp.oxconvert,
        'units': 'mmol O2/m^3'
    }

    cbarlabel = 'WOA Oxygen concentration, mmol O2/m^3'

    for plotKey in [
            'Pacific135W',
            'Atlantic28W',
    ]:  #'Equator','10 N', '10 S',]:
        oxcutoffs = [
            80.,
            20.,
            50.,
        ]
        for oxcutoff in oxcutoffs:
            filename = ukp.folder(['images', jobID, 'OMZ']) + '-'.join(
                [name, plotKey, 'contour',
                 str(int(oxcutoff))]) + '.png'
            title = ' '.join([
                jobID,
                getLongName(plotKey),
                getLongName(name) + ',',
                str(oxcutoff), modeldetails['units']
            ])
            contourplot(jobID,
                        name,
                        modelfiles,
                        datafile,
                        modeldetails,
                        datadetails,
                        medusaCoords,
                        woaCoords,
                        plotKey=plotKey,
                        oxcutoff=oxcutoff,
                        filename=filename,
                        title=title,
                        cbarlabel=cbarlabel)


if __name__ == "__main__":
    run()
    print('The end.')
