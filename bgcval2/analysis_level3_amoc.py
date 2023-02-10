#!/usr/bin/ipython

#
# Copyright 2015, Plymouth Marine Laboratory
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
#
"""
.. module:: analysis_level3_amoc
   :platform: Unix
   :synopsis: A script to produce a level 3 analysis for the AMOC.

.. moduleauthor:: Lee de Mora <ledm@pml.ac.uk>

"""

#####
# Load Standard Python modules:
from sys import argv, exit
from os.path import exists
from calendar import month_name
from socket import gethostname
from netCDF4 import Dataset
from glob import glob
from scipy.interpolate import interp1d
import numpy as np
import os, sys
from getpass import getuser
from matplotlib import pyplot
from shelve import open as shOpen

#####
# Load specific local code:
from bgcval2.bgcvaltools import bv2tools as bvt
from .timeseries import timeseriesAnalysis
from .timeseries import profileAnalysis
from .timeseries import timeseriesPlots as tsp

#####
# User defined set of paths pointing towards the datasets.
from .Paths import paths as paths


def analysis_omz(jobID=''):
    annual = True

    analysisKeys = []
    analysisKeys.append('AMOC_26N')
    #analysisKeys.append('AMOC_32S')

    analysisDict = {}
    imagedir = bvt.folder(paths.imagedir + '/' + jobID + '/Level3/AMOC')
    shelvedir = bvt.folder(paths.shelvedir + '/' + jobID + '/Level3/AMOC')
    if annual: WOAFolder = paths.WOAFolder_annual
    else: WOAFolder = paths.WOAFolder

    #####
    # make a link to the time series
    for a in analysisKeys:
        level1shelveFold = paths.shelvedir + '/timeseries/' + jobID
        files = glob(level1shelveFold + '/*' + a + '*')
        for f in files:
            lnfile = shelvedir + os.path.basename(f)
            if os.path.exists(lnfile): continue
            print("linking ", f, lnfile)
            #assert 0
            os.symlink(f, lnfile)

    dataD = {}
    modeldataD = {}

    def listModelDataFiles(jobID, filekey, datafolder, annual):
        if annual:
            return sorted(
                glob(datafolder + jobID + "/" + jobID + "o_1y_*_" + filekey +
                     ".nc"))
        else:
            return sorted(
                glob(datafolder + jobID + "/" + jobID + "o_1m_*_" + filekey +
                     ".nc"))

    medusaCoords = {
        't': 'time_counter',
        'z': 'deptht',
        'lat': 'nav_lat',
        'lon': 'nav_lon',
        'cal': '360_day',
    }  # model doesn't need time dict.

    av = bvt.AutoVivification()

    if 'AMOC_26N' in analysisKeys or 'AMOC_32S' in analysisKeys:
        # Note that this will only work with the eORCAgrid.
        latslice26N = slice(227, 228)
        latslice32S = slice(137, 138)
        e3v, e1v, tmask, alttmask = {}, {}, {}, {}
        for name in ['AMOC_26N', 'AMOC_32S']:
            if name not in analysisKeys: continue

            ####
            if name == 'AMOC_26N': latslice = latslice26N
            if name == 'AMOC_32S': latslice = latslice32S

            # Load grid data
            nc = Dataset(paths.orcaGridfn, 'r')
            e3v[name] = nc.variables['e3v'][:,
                                            latslice, :]  # z level height 3D
            e1v[name] = nc.variables['e1v'][latslice, :]  #
            tmask[name] = nc.variables['tmask'][:, latslice, :]
            nc.close()

            # load basin mask
            nc = Dataset('data/basinlandmask_eORCA1.nc', 'r')
            alttmask[name] = nc.variables['tmaskatl'][
                latslice, :]  # 2D Atlantic mask
            nc.close()

        def calc_amoc32S(nc, keys):
            name = 'AMOC_32S'
            zv = np.ma.array(nc.variables['vomecrty'][...,
                                                      latslice32S, :])  # m/s
            atlmoc = np.array(np.zeros_like(zv[0, :, :, 0]))
            e2vshape = e3v[name].shape
            for la in range(e2vshape[1]):  #ji, y
                for lo in range(e2vshape[2]):  #jj , x,
                    if int(alttmask[name][la, lo]) == 0: continue
                    for z in range(e2vshape[0]):  # jk
                        if int(tmask[name][z, la, lo]) == 0: continue
                        if np.ma.is_masked(zv[0, z, la, lo]): continue
                        atlmoc[z,
                               la] = atlmoc[z, la] - e1v[name][la, lo] * e3v[
                                   name][z, la, lo] * zv[0, z, la, lo] / 1.E06

####
# Cumulative sum from the bottom up.
            for z in range(73, 1, -1):
                atlmoc[z, :] = atlmoc[z + 1, :] + atlmoc[z, :]
            return np.ma.max(atlmoc)

        def calc_amoc26N(nc, keys):
            name = 'AMOC_26N'
            zv = np.ma.array(nc.variables['vomecrty'][...,
                                                      latslice26N, :])  # m/s
            atlmoc = np.array(np.zeros_like(zv[0, :, :, 0]))
            e2vshape = e3v[name].shape
            for la in range(e2vshape[1]):  #ji, y
                for lo in range(e2vshape[2]):  #jj , x,
                    if int(alttmask[name][la, lo]) == 0: continue
                    for z in range(e2vshape[0]):  # jk
                        if int(tmask[name][z, la, lo]) == 0: continue
                        if np.ma.is_masked(zv[0, z, la, lo]): continue
                        atlmoc[z,
                               la] = atlmoc[z, la] - e1v[name][la, lo] * e3v[
                                   name][z, la, lo] * zv[0, z, la, lo] / 1.E06

####
# Cumulative sum from the bottom up.
            for z in range(73, 1, -1):
                atlmoc[z, :] = atlmoc[z + 1, :] + atlmoc[z, :]
            return np.ma.max(atlmoc)

        for name in ['AMOC_26N', 'AMOC_32S']:
            if name not in analysisKeys: continue

            av[name]['modelFiles'] = listModelDataFiles(
                jobID, 'grid_V', paths.ModelFolder_pref, annual)
            av[name]['dataFile'] = ''

            av[name]['modelcoords'] = medusaCoords
            av[name]['datacoords'] = medusaCoords

            if name == 'AMOC_26N':
                av[name]['modeldetails'] = {
                    'name': name,
                    'vars': [
                        'vomecrty',
                    ],
                    'convert': calc_amoc26N,
                    'units': 'Sv'
                }
            if name == 'AMOC_32S':
                av[name]['modeldetails'] = {
                    'name': name,
                    'vars': [
                        'vomecrty',
                    ],
                    'convert': calc_amoc32S,
                    'units': 'Sv'
                }

            av[name]['datadetails'] = {
                'name': '',
                'units': '',
            }
            av[name]['layers'] = [
                'layerless',
            ]
            av[name]['regions'] = [
                'regionless',
            ]
            av[name]['metrics'] = [
                'metricless',
            ]
            av[name]['datasource'] = ''
            av[name]['model'] = 'NEMO'
            av[name]['modelgrid'] = 'eORCA1'
            av[name]['gridFile'] = paths.orcaGridfn
            av[name]['Dimensions'] = 1

#####
# Calling timeseriesAnalysis
# This is where the above settings is passed to timeseriesAnalysis, for the actual work to begin.
# We loop over all fiels in the first layer dictionary in the autovificiation, av.
#
# Once the timeseriesAnalysis has completed, we save all the output shelves in a dictionairy.
# At the moment, this dictioary is not used, but we could for instance open the shelve to highlight specific data,
#	(ie, andy asked to produce a table showing the final year of data.

    shelves = {}
    shelves_insitu = {}
    for name in list(av.keys()):
        print(
            "------------------------------------------------------------------"
        )
        print(
            "analysis-Timeseries.py:\tBeginning to call timeseriesAnalysis for ",
            name)

        if len(av[name]['modelFiles']) == 0:
            print(
                "analysis-Timeseries.py:\tWARNING:\tmodel files are not found:",
                name, av[name]['modelFiles'])
            if strictFileCheck: assert 0

        modelfilesexists = [os.path.exists(f) for f in av[name]['modelFiles']]
        if False in modelfilesexists:
            print(
                "analysis-Timeseries.py:\tWARNING:\tnot model files do not all exist:",
                av[name]['modelFiles'])
            for f in av[name]['modelFiles']:
                if os.path.exists(f): continue
                print(f, 'does not exist')
            if strictFileCheck: assert 0

        if av[name]['dataFile'] != '':
            if not os.path.exists(av[name]['dataFile']):
                print(
                    "analysis-Timeseries.py:\tWARNING:\tdata file is not found:",
                    av[name]['dataFile'])
                if strictFileCheck: assert 0


#####
# time series and traffic lights.
        tsa = timeseriesAnalysis(
            av[name]['modelFiles'],
            av[name]['dataFile'],
            dataType=name,
            modelcoords=av[name]['modelcoords'],
            modeldetails=av[name]['modeldetails'],
            datacoords=av[name]['datacoords'],
            datadetails=av[name]['datadetails'],
            datasource=av[name]['datasource'],
            model=av[name]['model'],
            jobID=jobID,
            layers=av[name]['layers'],
            regions=av[name]['regions'],
            metrics=av[name]['metrics'],
            workingDir=shelvedir,
            imageDir=imagedir,
            grid=av[name]['modelgrid'],
            gridFile=av[name]['gridFile'],
            clean=False,
            noNewFiles=True,  # stops loading new files               
        )
        amocshelve = tsa.shelvefn

        #####
        #
        sh = shOpen(amocshelve)
        readFiles = sh['readFiles']
        modeldataD = sh['modeldata']
        sh.close()

        m = 'metricless'
        r = 'regionless'
        l = 'layerless'

        modeldataDict = modeldataD[(r, l, m)]
        times = sorted(modeldataDict.keys())
        modeldata = [modeldataDict[t] for t in times]

        substract1000time = False
        if substract1000time:
            times = [t - 1000. for t in times]

        if len(times) > 120.:
            filename = bvt.folder(imagedir) + '_'.join([
                name,
                jobID,
            ]) + '_first100yrs.png'
            tmpt = times[:120]
            tmpmd = modeldata[:120]
            xlims = [tmpt[0], tmpt[-1]]

            fig = pyplot.figure()
            ax = fig.add_subplot(111)

            #if len(tmpmd)>30:
            #	smoothing = tsp.movingaverage2(tmpmd,window_len=30,window='hanning',extrapolate='axially')
            #	pyplot.plot(tmpt,tmpmd,      c='b',ls='-',lw=0.2)
            #	pyplot.plot(tmpt,smoothing,c='b',ls='-',lw=2,label='Model')
            #else:
            pyplot.plot(
                tmpt,
                tmpmd,
                c='b',
                ls='-',
                lw=2,
                label='Model',
            )

            if np.max(tmpmd) > 20.:
                pyplot.ylim([0., np.max(tmpmd) * 1.1])
            else:
                pyplot.ylim([0., 20.])

            pyplot.xlim(xlims)

            pyplot.ylabel('AMOC, Sv')
            pyplot.xlabel('Year')

            #pyplot.axhline(y='',c='k',ls='-',lw=1,label = 'data')

            print("timeseriesPlots:\tsimpletimeseries:\tSaving:", filename)
            pyplot.savefig(filename)
            pyplot.close()

        filename = bvt.folder(imagedir) + '_'.join([
            name,
            jobID,
        ]) + '.png'
        xlims = [times[0], times[-1]]

        fig = pyplot.figure()
        ax = fig.add_subplot(111)

        if len(modeldata) > 30:
            smoothing = tsp.movingaverage2(modeldata,
                                           window_len=30,
                                           window='hanning',
                                           extrapolate='axially')
            pyplot.plot(times,
                        modeldata,
                        c='b',
                        ls='-',
                        lw=0.2,
                        label='Model1')
            pyplot.plot(times, smoothing, c='b', ls='-', lw=2, label='Model2')
        else:
            pyplot.plot(
                times,
                modeldata,
                c='b',
                ls='-',
                lw=1,
                label='Model',
            )

        if np.max(modeldata) > 20.:
            pyplot.ylim([0., np.max(modeldata) * 1.1])
        else:
            pyplot.ylim([0., 20.])

        pyplot.xlim(xlims)
        pyplot.ylabel('AMOC, Sv')
        pyplot.xlabel('Year')

        #pyplot.axhline(y='',c='k',ls='-',lw=1,label = 'data')

        print("timeseriesPlots:\tsimpletimeseries:\tSaving:", filename)
        pyplot.savefig(filename)
        pyplot.close()


def main():
    if "--help" in argv or len(argv) == 1:
        print("Running with no arguments. Exiting.")
        if "--help" in argv:
            print("Read the documentation.")
        sys.exit(0)
    try:
        jobID = argv[1]
    except:
        jobID = "u-ag543"

    if 'debug' in argv[1:]: suite = 'debug'

    analysis_omz(jobID=jobID, )  #clean=1)
    #if suite == 'all':
    #analysis_timeseries(jobID =jobID,analysisSuite='FullDepth', z_component = 'FullDepth',)#clean=1)


if __name__ == "__main__":
    main()
