#!/usr/bin/env python
########!/usr/bin/python

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
.. module:: analysis_compare
   :platform: Unix
   :synopsis: A tool that generates an intercomparison of multiple UKESM jobs time series analyses.
.. moduleauthor:: Lee de Mora <ledm@pml.ac.uk>
.. moduleauthor:: Valeriu Predoi <valeriu.predoi@ncas.ac.uk>

"""

import argparse
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')

#####
# Load Standard Python modules:
from calendar import month_name
from socket import gethostname
from netCDF4 import Dataset
from glob import glob
from scipy.interpolate import interp1d
import numpy as np
import os
import sys
import fnmatch
from getpass import getuser
from collections import defaultdict
import yaml
import random
import itertools


#####
# Load specific local code:
from . import UKESMpython as ukp
from .timeseries import timeseriesAnalysis
from .timeseries import profileAnalysis
from .timeseries import timeseriesPlots as tsp
from bgcval2.analysis_timeseries import analysis_timeseries, build_list_of_suite_keys, load_key_file
from bgcval2.download_from_mass import download_from_mass


try:
    from .bgcvaltools.pftnames import getLongName
except:
    from .pftnames import getLongName
from .bgcvaltools.mergeMonthlyFiles import mergeMonthlyFiles, meanDJF
from .netcdf_manipulation.alwaysInclude import alwaysInclude
from .bgcval2_make_report import comparehtml5Maker
#from .Paths import paths

#from .comparison.shifttimes import shifttimes as shifttimes_legacy
from .comparison.ensembles import build_ensemble
from .config.configToDict import configToDict
from .bgcvaltools.dataset import dataset
from ._runtime_config import get_run_configuration

#####
# User defined set of paths pointing towards the datasets.
from .Paths.paths import paths_setter


def titleify(ls):
    """
    Turnms a list into a title string.
    """
    return ' '.join([getLongName(i) for i in ls])


def listModelDataFiles(jobID, filekey, datafolder, annual, year=''):
    """

    """
    if year == '':
        if annual:
            return sorted(
                glob(datafolder + jobID + "/" + jobID + "o_1y_*_" + filekey +
                     ".nc"))
        else:
            return sorted(
                glob(datafolder + jobID + "/" + jobID + "o_1m_*_" + filekey +
                     ".nc"))
    else:
        if annual:
            print(datafolder + jobID + "/" + jobID + "o_1y_*" + year +
                   "????_" + filekey + ".nc")
            return sorted(
                glob(datafolder + jobID + "/" + jobID + "o_1y_*" + year +
                     "????_" + filekey + ".nc"))
        else:
            return sorted(
                glob(datafolder + jobID + "/" + jobID + "o_1m_*" + year +
                     "????_" + filekey + ".nc"))


def apply_shifttimes(mdata, jobID, shifttimes):
    """
    Replaces .comparison.shifttimes loaded as shifttimes_legacy
    which takes the model data, the jobID and the year0

    This version takes the mdata, and shifttime value - a number to add to the time axis.
    the value of the shift is provided in the yaml file.

    Outputs two lists: dates & data.
    """
    times, datas = [], []
    if not len(mdata.keys()):
        return [], []

    t0 = float(sorted(mdata.keys())[0])
    for t in sorted(mdata.keys()):
        t1 = t + float(shifttimes[jobID])
        times.append(t1)
        datas.append(mdata[t])
    return times, datas


def apply_timerange(times, datas, jobID, timeranges):
    """
    This version takes the times and dataafter apply_shifttimes ,
    and removes things outside the time range requested.
    the value of the range is provided in the yaml file.

    Outputs two lists: dates & data.
    """
    if 0 in [len(times), len(datas), ]:
       return times, datas

    print('apply_timerange', jobID, timeranges)

    timerange = timeranges.get(jobID, None)

    print('apply_timerange', timerange)    
    if timerange is None: 
       print('apply_timerange', timerange, 'is', None) 
       return times, datas

    print(jobID, timerange, 'is not None', np.min(times), np.max(times))
    
    n_times, n_datas = [], [] # to ensure they stay lists
    for ti, da in zip(times, datas):
        if ti < np.min(timerange):
            continue
        if ti > np.max(timerange):
            continue
        n_times.append(ti)
        n_datas.append(da)
        print('apply_timerange:', jobID, ti, da)

    if not len(n_times):
        print('apply_timerange: WARNING: No times made the cut?', len(times),
              'original times', [np.min(times), np.max(times)], 
              'timerange:', timerange)
        assert 0
    return n_times, n_datas


def timeseries_compare(jobs,
                       colours,
                       suites = [],
                       analysisname='',
                       shifttimes={},
                       timeranges={},   
                       jobDescriptions={},
                       lineThicknesses=defaultdict(lambda: 1),
                       linestyles=defaultdict(lambda: '-'),
                       ensembles={},
                       config_user=None):
    """
    timeseries_compare:
        Suite of tools to take pre-analyses time series model data
        then compile into single plots, then publish it to an html
        document.

    """
    ### strategy here is a simple wrapper.
    # It's a little cheat-y, as I'm copying straight from analysis_timeseries.py

    jobs = sorted(jobs)
    #jobs = sorted(colours.keys())

    for ensemble in list(ensembles.keys()):
        # ensembles names can not be the same as jobIDs
        jobs.remove(ensemble)

    # get runtime configuration
    if not config_user:
        paths_dict, config_user = get_run_configuration("defaults")
    else:
        paths_dict, config_user = get_run_configuration(config_user)

    # filter paths dict into an object that's usable below
    paths = paths_setter(paths_dict)
    if analysisname == '':
        print('ERROR: please provide an name for this analsys')
        sys.exit(0)
    else:
        imageFolder = paths.imagedir + '/TimeseriesCompare/' + analysisname

    annual = True
    strictFileCheck = False

    if not isinstance(suites, list):
        ValueError(f"Suites need to be a list, got: {suites}")
        sys.exit(1)

    analysisKeys = build_list_of_suite_keys(suites, debug=True)
    print(f'Using analysis keys {str(analysisKeys)}')

    layerList = [
        'Surface',
    ]
    metricList = [
        'mean',
    ]
    regionList = [
        'Global',
    ]

    PierceRegions = [
        'Enderby',
        'Wilkes',
        'Ross',
        'Amundsen',
        'Weddel',
    ]

    vmtregionList = [
        'Global',
        'Depth_700m',
        'Depth_2000m',
        'Depth_700-2000m',
        'ignoreInlandSeas',
        'Equator10',
        'AtlanticSOcean',
        'SouthernOcean',
        'ArcticOcean',
        'Remainder',
        'NorthernSubpolarAtlantic',
        'NorthernSubpolarPacific',
        'WeddelSea',
        'Cornwall',
    ]
    #vmtregionList = ['Global', 'ignoreInlandSeas','Equator10','AtlanticSOcean','SouthernOcean','ArcticOcean',  'Remainder','NorthernSubpolarAtlantic','NorthernSubpolarPacific','WeddelSea']
    vmtregionList.extend(PierceRegions)
    OMZRegions = [
        'EquatorialPacificOcean', 'IndianOcean', 'EquatorialAtlanticOcean'
    ]  #'Ross','Amundsen','Weddel',]
    level3 = [
        'DMS',
    ]

    #####
    # paths:
    orcaGridfn = paths.orcaGridfn  #'/group_workspaces/jasmin4/esmeval/example_data/bgc/mesh_mask_eORCA1_wrk.nc'
    if annual: WOAFolder = paths.WOAFolder_annual
    else: WOAFolder = paths.WOAFolder

    #####
    # Coordinate dictionairy
    # These are python dictionairies, one for each data source and model.
    # This is because each data provider seems to use a different set of standard names for dimensions and time.
    # The 'tdict' field is short for "time-dictionary".
    #	This is a dictionary who's indices are the values on the netcdf time dimension.
    #	The tdict indices point to a month number in python numbering (ie January = 0)
    # 	An example would be, if a netcdf uses the middle day of the month as it's time value:
    #		tdict = {15:0, 45:1 ...}

    dataD = {}
    modeldataD = {}

    for jobID in jobs:

        #####
        # Location of images directory
        # the imagedir is where the analysis images will be saved.
        imagedir = ukp.folder(paths.imagedir + '/' + jobID + '/timeseries')
        shelvedir = ukp.folder(paths.shelvedir + "/timeseries/" + jobID)

        if jobID in list(ensembles.keys()): continue
        # ensembles names can not be the same as jobIDs

        av = ukp.AutoVivification()

        # NEW STYLE keys from file:
        for key in analysisKeys:
            av[key] = load_key_file(key, paths, jobID)

        if 'GlobalExportRatio' in analysisKeys:

            def calcExportRatio(nc, keys):
                a = (nc.variables['SDT__100'][:] + nc.variables['FDT__100'][:]
                     ).sum() / (nc.variables['PRD'][:] +
                                nc.variables['PRN'][:]).sum()
                #a = np.ma.masked_where(a>1.01, a)
                return a

            name = 'ExportRatio'
            av[name]['modelFiles'] = listModelDataFiles(
                jobID, 'diad_T', paths.ModelFolder_pref, annual)

            av[name]['dataFile'] = ""
            av[name]['modelcoords'] = medusaCoords
            av[name]['datacoords'] = maredatCoords
            av[name]['modeldetails'] = {
                'name': name,
                'vars': [
                    'SDT__100',
                    'FDT__100',
                    'PRD',
                    'PRN',
                ],
                'convert': calcExportRatio,
                'units': ''
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
            av[name]['model'] = 'MEDUSA'
            av[name]['modelgrid'] = 'eORCA1'
            av[name]['gridFile'] = orcaGridfn
            av[name]['Dimensions'] = 1

        if 'CHD' in analysisKeys or 'CHN' in analysisKeys:
            for name in [
                    'CHD',
                    'CHN',
            ]:
                if name not in analysisKeys: continue

                av[name]['modelFiles'] = listModelDataFiles(
                    jobID, 'ptrc_T', paths.ModelFolder_pref, annual)
                av[name]['dataFile'] = ''

                av[name]['modelcoords'] = medusaCoords
                av[name]['datacoords'] = ''

                av[name]['modeldetails'] = {
                    'name': name,
                    'vars': [
                        name,
                    ],
                    'convert': ukp.NoChange,
                    'units': 'mg C/m^3'
                }
                av[name]['datadetails'] = {'name': '', 'units': ''}

                av[name]['layers'] = [
                    'Surface',
                    '100m',
                ]  # CCI is surface only, it's a satellite product.
                av[name]['regions'] = regionList
                av[name]['metrics'] = metricList  #['mean','median', ]

                av[name]['datasource'] = ''
                av[name]['model'] = 'MEDUSA'

                av[name]['modelgrid'] = 'eORCA1'
                av[name]['gridFile'] = paths.orcaGridfn
                av[name]['Dimensions'] = 2

        if 'DiaFrac' in analysisKeys:

            name = 'DiaFrac'

            def caldiafrac(nc, keys):
                return 100. * nc.variables[keys[0]][:].squeeze() / (
                    nc.variables[keys[0]][:].squeeze() +
                    nc.variables[keys[1]][:].squeeze())

            av[name]['modelFiles'] = listModelDataFiles(
                jobID, 'ptrc_T', paths.ModelFolder_pref, annual)
            av[name]['dataFile'] = ''

            av[name]['modelcoords'] = medusaCoords
            av[name]['datacoords'] = ''

            av[name]['modeldetails'] = {
                'name': name,
                'vars': [
                    'CHD',
                    'CHN',
                ],
                'convert': caldiafrac,
                'units': '%'
            }
            av[name]['datadetails'] = {'name': '', 'units': ''}

            av[name]['layers'] = [
                'Surface',
                '100m',
            ]  # CCI is surface only, it's a satellite product.
            av[name]['regions'] = regionList
            av[name]['metrics'] = metricList  #['mean','median', ]

            av[name]['datasource'] = ''
            av[name]['model'] = 'MEDUSA'

            av[name]['modelgrid'] = 'eORCA1'
            av[name]['gridFile'] = paths.orcaGridfn
            av[name]['Dimensions'] = 2
        if 'DTC' in analysisKeys:
            for name in [
                    'DTC',
            ]:
                if name not in analysisKeys: continue

                av[name]['modelFiles'] = listModelDataFiles(
                    jobID, 'ptrc_T', paths.ModelFolder_pref, annual)
                av[name]['dataFile'] = ''

                av[name]['modelcoords'] = medusaCoords
                av[name]['datacoords'] = ''

                av[name]['modeldetails'] = {
                    'name': name,
                    'vars': [
                        name,
                    ],
                    'convert': ukp.mul1000,
                    'units': 'umol-C/m3'
                }
                av[name]['datadetails'] = {'name': '', 'units': ''}

                av[name]['layers'] = [
                    '3000m',
                ]  #'100m',]         # CCI is surface only, it's a satellite product.
                av[name]['regions'] = regionList
                av[name]['metrics'] = metricList  #['mean','median', ]

                av[name]['datasource'] = ''
                av[name]['model'] = 'MEDUSA'

                av[name]['modelgrid'] = 'eORCA1'
                av[name]['gridFile'] = paths.orcaGridfn
                av[name]['Dimensions'] = 3

        if 'OMZMeanDepth' in analysisKeys:
            if annual:
                av['OMZMeanDepth']['modelFiles'] = sorted(
                    glob(paths.ModelFolder_pref + jobID + "/" + jobID +
                         "o_1y_*_ptrc_T.nc"))
                av['OMZMeanDepth'][
                    'dataFile'] = WOAFolder + 'woa13_all_o00_01.nc'
            else:
                print("OMZ Thickness not implemented for monthly data")
                assert 0

            nc = Dataset(paths.orcaGridfn, 'r')
            depths = np.abs(nc.variables['gdepw'][:])
            tmask = nc.variables['tmask'][:]
            nc.close()

            omzthreshold = 20.

            def modelMeanOMZdepth(nc, keys):
                o2 = nc.variables[keys[0]][:].squeeze()
                meandepth = np.ma.masked_where(
                    (o2 > omzthreshold) + o2.mask + (tmask == 0),
                    depths).mean(0)
                if meandepth.max() in [0., 0]: return np.array([
                        0.,
                ])
                return np.ma.masked_where(meandepth == 0., meandepth)

            def woaMeanOMZdepth(nc, keys):
                o2 = nc.variables[keys[0]][:].squeeze() * 44.661
                pdepths = np.zeros_like(o2)
                lons = nc.variables['lon'][:]
                lats = nc.variables['lat'][:]
                wdepths = np.abs(nc.variables['depth'][:])

                for y, lat in enumerate(lats):
                    for x, lon in enumerate(lons):
                        pdepths[:, y, x] = wdepths
                wmeanDepth = np.ma.masked_where((o2 > omzthreshold) + o2.mask,
                                                pdepths).mean(0).data
                print("woaMeanOMZdepth", wmeanDepth.min(), wmeanDepth.mean(),
                      wmeanDepth.max())
                #assert 0

                if wmeanDepth.max() in [0., 0]: return np.array([
                        1000.,
                ])
                return np.ma.masked_where(wmeanDepth == 0., wmeanDepth)

            av['OMZMeanDepth']['modelcoords'] = medusaCoords
            av['OMZMeanDepth']['datacoords'] = woaCoords

            av['OMZMeanDepth']['modeldetails'] = {
                'name': 'OMZMeanDepth',
                'vars': [
                    'OXY',
                ],
                'convert': modelMeanOMZdepth,
                'units': 'm'
            }
            av['OMZMeanDepth']['datadetails'] = {
                'name': 'OMZMeanDepth',
                'vars': [
                    'o_an',
                ],
                'convert': woaMeanOMZdepth,
                'units': 'm'
            }

            av['OMZMeanDepth']['layers'] = [
                'layerless',
            ]
            av['OMZMeanDepth']['regions'] = regionList
            av['OMZMeanDepth']['metrics'] = metricList

            av['OMZMeanDepth']['datasource'] = 'WOA'
            av['OMZMeanDepth']['model'] = 'MEDUSA'

            av['OMZMeanDepth']['modelgrid'] = 'eORCA1'
            av['OMZMeanDepth']['gridFile'] = paths.orcaGridfn
            av['OMZMeanDepth']['Dimensions'] = 2

        if 'OMZThickness' in analysisKeys or 'OMZThickness50' in analysisKeys:
            if 'OMZThickness' in analysisKeys and 'OMZThickness50' in analysisKeys:
                print("Only run one of these at a time")
                assert 0

            if annual:
                av['OMZThickness']['modelFiles'] = sorted(
                    glob(paths.ModelFolder_pref + jobID + "/" + jobID +
                         "o_1y_*_ptrc_T.nc"))
                av['OMZThickness'][
                    'dataFile'] = WOAFolder + 'woa13_all_o00_01.nc'
            else:
                print("OMZ Thickness not implemented for monthly data")
                assert 0

            nc = Dataset(paths.orcaGridfn, 'r')
            thickness = nc.variables['e3t'][:]
            tmask = nc.variables['tmask'][:]
            nc.close()

            if 'OMZThickness' in analysisKeys: omzthreshold = 20.
            if 'OMZThickness50' in analysisKeys: omzthreshold = 50.

            def modelOMZthickness(nc, keys):
                o2 = nc.variables[keys[0]][:].squeeze()
                totalthick = np.ma.masked_where(
                    (o2 > omzthreshold) + o2.mask + (tmask == 0),
                    thickness).sum(0).data
                if totalthick.max() in [0., 0]: return np.array([
                        0.,
                ])

                return np.ma.masked_where(totalthick == 0., totalthick)
                #return np.ma.masked_where((arr>omzthreshold) + (arr <0.) + arr.mask,thickness).sum(0)

            def woaOMZthickness(nc, keys):
                o2 = nc.variables[keys[0]][:].squeeze() * 44.661
                pthick = np.zeros_like(o2)
                lons = nc.variables['lon'][:]
                lats = nc.variables['lat'][:]
                zthick = np.abs(nc.variables['depth_bnds'][:, 0] -
                                nc.variables['depth_bnds'][:, 1])

                for y, lat in enumerate(lats):
                    for x, lon in enumerate(lons):
                        pthick[:, y, x] = zthick
                totalthick = np.ma.masked_where((o2 > omzthreshold) + o2.mask,
                                                pthick).sum(0).data

                if totalthick.max() in [0., 0]: return np.array([
                        0.,
                ])
                return np.ma.masked_where(totalthick == 0., totalthick)

            av['OMZThickness']['modelcoords'] = medusaCoords
            av['OMZThickness']['datacoords'] = woaCoords

            av['OMZThickness']['modeldetails'] = {
                'name': 'OMZThickness',
                'vars': [
                    'OXY',
                ],
                'convert': modelOMZthickness,
                'units': 'm'
            }
            av['OMZThickness']['datadetails'] = {
                'name': 'OMZThickness',
                'vars': [
                    'o_an',
                ],
                'convert': woaOMZthickness,
                'units': 'm'
            }

            av['OMZThickness']['layers'] = [
                'layerless',
            ]
            av['OMZThickness']['regions'] = regionList
            av['OMZThickness']['metrics'] = metricList

            av['OMZThickness']['datasource'] = 'WOA'
            av['OMZThickness']['model'] = 'MEDUSA'

            av['OMZThickness']['modelgrid'] = 'eORCA1'
            av['OMZThickness']['gridFile'] = paths.orcaGridfn
            av['OMZThickness']['Dimensions'] = 2

        if 'TotalOMZVolume' in analysisKeys or 'TotalOMZVolume50' in analysisKeys:
            if 'TotalOMZVolume' in analysisKeys and 'TotalOMZVolume50' in analysisKeys:
                print("Only run one of these at a time")
                assert 0

            if annual:
                av['TotalOMZVolume']['modelFiles'] = sorted(
                    glob(paths.ModelFolder_pref + jobID + "/" + jobID +
                         "o_1y_*_ptrc_T.nc"))
                av['TotalOMZVolume'][
                    'dataFile'] = WOAFolder + 'woa13_all_o00_01.nc'
            else:
                print("OMZ volume not implemented for monthly data")
                assert 0

            nc = Dataset(paths.orcaGridfn, 'r')
            try:
                pvol = nc.variables['pvol'][:]
                tmask = nc.variables['tmask'][:]
            except:
                tmask = nc.variables['tmask'][:]
                area = nc.variables['e2t'][:] * nc.variables['e1t'][:]
                pvol = nc.variables['e3t'][:] * area
                pvol = np.ma.masked_where(tmask == 0, pvol)
            nc.close()

            if 'TotalOMZVolume' in analysisKeys: omzthreshold = 20.
            if 'TotalOMZVolume50' in analysisKeys: omzthreshold = 50.

            def modelTotalOMZvol(nc, keys):
                arr = np.ma.array(nc.variables[keys[0]][:].squeeze())
                return np.ma.masked_where(
                    (arr > omzthreshold) + pvol.mask + arr.mask, pvol).sum()

            def woaTotalOMZvol(nc, keys):
                arr = nc.variables[keys[0]][:].squeeze() * 44.661
                #area = np.zeros_like(arr[0])
                pvol = np.zeros_like(arr)
                #np.ma.masked_wjhere(arr.mask + (arr <0.)+(arr >1E10),np.zeros_like(arr))
                lons = nc.variables['lon'][:]
                lats = nc.variables['lat'][:]
                #lonbnds = nc.variables['lon_bnds'][:]
                latbnds = nc.variables['lat_bnds'][:]
                zthick = np.abs(nc.variables['depth_bnds'][:, 0] -
                                nc.variables['depth_bnds'][:, 1])

                for y, lat in enumerate(lats):
                    area = ukp.Area([latbnds[y, 0], 0.], [latbnds[y, 1], 1.])
                    for z, thick in enumerate(zthick):
                        pvol[z, y, :] = np.ones_like(lons) * area * thick

                return np.ma.masked_where(
                    arr.mask + (arr > omzthreshold) + (arr < 0.), pvol).sum()

            av['TotalOMZVolume']['modelcoords'] = medusaCoords
            av['TotalOMZVolume']['datacoords'] = woaCoords

            av['TotalOMZVolume']['modeldetails'] = {
                'name': 'TotalOMZVolume',
                'vars': [
                    'OXY',
                ],
                'convert': modelTotalOMZvol,
                'units': 'm^3'
            }
            av['TotalOMZVolume']['datadetails'] = {
                'name': 'TotalOMZVolume',
                'vars': [
                    'o_an',
                ],
                'convert': woaTotalOMZvol,
                'units': 'm^3'
            }

            av['TotalOMZVolume']['layers'] = [
                'layerless',
            ]
            av['TotalOMZVolume']['regions'] = [
                'regionless',
            ]
            av['TotalOMZVolume']['metrics'] = [
                'metricless',
            ]

            av['TotalOMZVolume']['datasource'] = 'WOA'
            av['TotalOMZVolume']['model'] = 'MEDUSA'

            av['TotalOMZVolume']['modelgrid'] = 'eORCA1'
            av['TotalOMZVolume']['gridFile'] = paths.orcaGridfn
            av['TotalOMZVolume']['Dimensions'] = 1

        if 'VolumeMeanOxygen' in analysisKeys:
            name = 'VolumeMeanOxygen'
            av[name]['modelFiles'] = listModelDataFiles(
                jobID, 'ptrc_T', paths.ModelFolder_pref, annual)
            av[name]['dataFile'] = ''
            av[name]['modelcoords'] = medusaCoords
            av[name]['datacoords'] = woaCoords

            nc = dataset(paths.orcaGridfn, 'r')
            try:
                pvol = nc.variables['pvol'][:]
                area = nc.variables['area'][:]
                gmttmask = nc.variables['tmask'][:]
            except:
                gmttmask = nc.variables['tmask'][:]
                area = nc.variables['e2t'][:] * nc.variables['e1t'][:]
                pvol = nc.variables['e3t'][:] * area
                pvol = np.ma.masked_where(gmttmask == 0, pvol)
            nc.close()

            def sumMeanLandMask(nc, keys):
                #### works like no change, but applies a mask.
                ox = np.ma.array(nc.variables[keys[0]][:].squeeze())
                ox = np.ma.masked_where((gmttmask == 0) + (ox.mask), ox)

                try:
                    vol = np.ma.masked_where(
                        ox.mask,
                        nc('thkcello')[:].squeeze() *
                        nc('area')[:])  # preferentially use in file volume.
                except:
                    vol = np.ma.masked_where(ox.mask, pvol)

                return ((ox * vol).sum(0) / vol.sum(0))  #*(area/area.sum())

            av[name]['modeldetails'] = {
                'name': name,
                'vars': [
                    'OXY',
                ],
                'convert': sumMeanLandMask,
                'units': 'mmol O2/m^3'
            }
            av[name]['datadetails'] = {'name': '', 'units': ''}
            #av[name]['datadetails']  	= {'name': name, 'vars':['t_an',], 'convert': ukp.NoChange,'units':'degrees C'}

            oxregions = [
                'Global',
                'ignoreInlandSeas',
                'Equator10',
                'AtlanticSOcean',
                'SouthernOcean',
                'ArcticOcean',
                'Remainder',
                'NorthernSubpolarAtlantic',
                'NorthernSubpolarPacific',
            ]  #'WeddelSea']
            oxregions.extend(OMZRegions)

            av[name]['layers'] = [
                'layerless',
            ]
            av[name]['regions'] = oxregions
            av[name]['metrics'] = [
                'wcvweighted',
            ]
            av[name]['datasource'] = ''
            av[name]['model'] = 'NEMO'
            av[name]['modelgrid'] = 'eORCA1'
            av[name]['gridFile'] = paths.orcaGridfn
            av[name]['Dimensions'] = 2

        if 'ERSST' in analysisKeys:
            name = 'ERSST'
            av[name][
                'modelFiles'] = paths.ObsFolder + "/bgc/ERSST.v4/sst.mnmean.v4.nc"
            ERSSTCoords = {
                't': 'time',
                'z': '',
                'lat': 'lat',
                'lon': 'lon',
                'cal': 'standard',
                'tdict': ['ZeroToZero']
            }

            av[name]['modelcoords'] = ERSSTCoords
            #av[name]['datacoords']          = woaCoords

            av[name]['modeldetails'] = {
                'name': name,
                'vars': [
                    'sst',
                ],
                'convert': ukp.NoChange,
                'units': 'degrees C'
            }
            av[name]['datadetails'] = {'name': name, 'vars': [], 'units': ''}

            av[name]['layers'] = [
                'layerless',
            ]
            av[name]['regions'] = regionList
            av[name]['metrics'] = metricList

            av[name]['datasource'] = ''
            av[name]['model'] = 'ERSST'

            av[name]['modelgrid'] = 'ERSST_2g'
            av[name][
                'gridFile'] = paths.ObsFolder + "/bgc/ERSST.v4/sst.mnmean.v4.nc"
            av[name]['Dimensions'] = 2


        if 'VolumeMeanTemperature' in analysisKeys:
            name = 'VolumeMeanTemperature'
            av[name]['modelFiles'] = listModelDataFiles(
                jobID, 'grid_T', paths.ModelFolder_pref, annual)
            av[name]['dataFile'] = ''
            av[name]['modelcoords'] = medusaCoords
            av[name]['datacoords'] = woaCoords

            def sumMeanLandMask(nc, keys):
                assert 0

            av[name]['modeldetails'] = {
                'name': name,
                'vars': [
                    'votemper',
                ],
                'convert': sumMeanLandMask,
                'units': 'degrees C'
            }
            av[name]['datadetails'] = {'name': '', 'units': ''}
            av[name]['layers'] = [
                'layerless',
            ]
            av[name]['regions'] = vmtregionList
            av[name]['metrics'] = [
                'wcvweighted',
            ]
            av[name]['datasource'] = ''
            av[name]['model'] = 'NEMO'
            av[name]['modelgrid'] = 'eORCA1'
            av[name]['gridFile'] = paths.orcaGridfn
            av[name]['Dimensions'] = 2

        if 'scvoltot' in analysisKeys:
            name = 'scvoltot'
            av[name]['modelFiles'] = listModelDataFiles(
                jobID, 'scalar', paths.ModelFolder_pref, annual)
            av[name]['dataFile'] = ''
            av[name]['modelcoords'] = {
                'lat': False,
                'lon': False,
                'z': False,
                't': 'time_centered',
            }
            av[name]['datacoords'] = woaCoords
            av[name]['modeldetails'] = {
                'name': name,
                'vars': [
                    'scvoltot',
                ],
                'convert': ukp.NoChange,
                'units': 'm3'
            }
            av[name]['datadetails'] = {'name': '', 'units': ''}
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

        if 'soga' in analysisKeys:
            name = 'soga'
            av[name]['modelFiles'] = listModelDataFiles(
                jobID, 'scalar', paths.ModelFolder_pref, annual)
            av[name]['dataFile'] = ''
            av[name]['modelcoords'] = {
                'lat': False,
                'lon': False,
                'z': False,
                't': 'time_centered',
            }
            av[name]['datacoords'] = woaCoords
            av[name]['modeldetails'] = {
                'name': name,
                'vars': [
                    'soga',
                ],
                'convert': ukp.NoChange,
                'units': 'psu'
            }
            av[name]['datadetails'] = {'name': '', 'units': ''}
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

        if 'thetaoga' in analysisKeys:
            name = 'thetaoga'
            av[name]['modelFiles'] = listModelDataFiles(
                jobID, 'scalar', paths.ModelFolder_pref, annual)
            av[name]['dataFile'] = ''
            av[name]['modelcoords'] = {
                'lat': False,
                'lon': False,
                'z': False,
                't': 'time_centered',
            }
            av[name]['datacoords'] = woaCoords
            av[name]['modeldetails'] = {
                'name': name,
                'vars': [
                    'thetaoga',
                ],
                'convert': ukp.NoChange,
                'units': 'degrees C'
            }
            av[name]['datadetails'] = {'name': '', 'units': ''}
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

        if 'scalarHeatContent' in analysisKeys:
            name = 'scalarHeatContent'
            files = listModelDataFiles(jobID, 'scalar', paths.ModelFolder_pref,
                                       annual)

            def scalarFunction(nc, keys):
                rau0 = 1026.  #kg / m3	#		volume reference mass,
                rcp = 3991.8679571196299  #J / (K * kg)	ocean specific heat capacity
                thetaoga = nc(
                    'thetaoga'
                )[:]  #		global average seawater potential temperature
                scvoltot = nc('scvoltot')[:]  # m3		ocean volume

                return thetaoga * scvoltot * rau0 * rcp * 1e-24

            if len(files) > 0:
                av[name]['modelFiles'] = files
                av[name]['dataFile'] = ''
                av[name]['modelcoords'] = {
                    'lat': False,
                    'lon': False,
                    'z': False,
                    't': 'time_centered',
                }
                av[name]['datacoords'] = woaCoords
                av[name]['modeldetails'] = {
                    'name': name,
                    'vars': [
                        'thetaoga',
                        'scvoltot',
                    ],
                    'convert': ukp.NoChange,
                    'units': 'YottaJoules'
                }
                av[name]['datadetails'] = {'name': '', 'units': ''}
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


        if 'IcelessMeanSST' in analysisKeys:
            name = 'IcelessMeanSST'
            av[name]['modelFiles'] = listModelDataFiles(
                jobID, 'grid_T', paths.ModelFolder_pref, annual)
            av[name]['dataFile'] = ''

            av[name]['modelcoords'] = medusaCoords
            av[name]['datacoords'] = woaCoords

            nc = Dataset(paths.orcaGridfn, 'r')
            tmask = nc.variables['tmask'][:]
            area_full = nc.variables['e2t'][:] * nc.variables['e1t'][:]
            nc.close()

            def calcIcelessMeanSST(nc, keys):
                #### works like no change, but applies a mask.
                icecov = nc.variables['soicecov'][:].squeeze()
                sst = nc.variables['votemper'][:, 0, ].squeeze()
                sst = np.ma.masked_where(
                    (tmask[0] == 0) + (icecov > 0.15) + sst.mask, sst)
                area = np.ma.masked_where(sst.mask, area_full)
                val = (sst * area).sum() / (area.sum())
                print("calcIcelessMeanSST", sst.shape, area.shape, val)
                return val

            av[name]['modeldetails'] = {
                'name': name,
                'vars': [
                    'soicecov',
                    'votemper',
                ],
                'convert': calcIcelessMeanSST,
                'units': 'degrees C'
            }
            av[name]['datadetails'] = {'name': '', 'units': ''}
            #av[name]['datadetails']  	= {'name': name, 'vars':['t_an',], 'convert': ukp.NoChange,'units':'degrees C'}

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

        if 'MaxMonthlyMLD' in analysisKeys or 'MinMonthlyMLD' in analysisKeys:

            #/group_workspaces/jasmin2/ukesm/BGC_data/u-ad371/monthlyMLD/MetOffice_data_licence.325210916
            monthlyFiles = glob(paths.ModelFolder_pref + '/' + jobID +
                                '/monthlyMLD/' + jobID + 'o_1m_*_grid_T.nc')
            if len(monthlyFiles):
                maxmldfiles = mergeMonthlyFiles(monthlyFiles,
                                                outfolder='',
                                                cal=medusaCoords['cal'])

                for name in ['MaxMonthlyMLD', 'MinMonthlyMLD']:
                    if name not in analysisKeys: continue

                    def mldapplymask(nc, keys):
                        mld = np.ma.array(nc.variables[keys[0]][:]).max(0)
                        mld = np.ma.masked_where(
                            (nc.variables[keys[1]][:] == 0.) + mld.mask +
                            (mld == 1.E9), mld)
                        return mld

                    def mldmonthlymask(nc, keys):
                        mld = np.ma.array(np.ma.abs(
                            nc.variables[keys[0]][:])).max(0)
                        mld = np.ma.masked_where(mld.mask + (mld.data > 1.E10),
                                                 mld)
                        return mld

                    def mldapplymask_min(nc, keys):
                        mld = np.ma.array(nc.variables[keys[0]][:]).min(0)
                        mld = np.ma.masked_where(
                            (nc.variables[keys[1]][:] == 0.) + mld.mask +
                            (mld == 1.E9), mld)
                        return mld

                    def mldmonthlymask_min(nc, keys):
                        mld = np.ma.array(np.ma.abs(
                            nc.variables[keys[0]][:])).min(0)
                        mld = np.ma.masked_where(mld.mask + (mld.data > 1.E10),
                                                 mld)
                        return mld

                    av[name][
                        'modelFiles'] = maxmldfiles  #listModelDataFiles(jobID, 'grid_T', paths.ModelFolder_pref, annual)
                    av[name][
                        'dataFile'] = paths.MLDFolder + "mld_DT02_c1m_reg2.0.nc"  #mld_DT02_c1m_reg2.0.nc"

                    av[name]['modelcoords'] = medusaCoords
                    av[name]['datacoords'] = mldCoords

                    if name == 'MaxMonthlyMLD':
                        av[name]['modeldetails'] = {
                            'name': 'mld',
                            'vars': [
                                'somxl010',
                            ],
                            'convert': mldmonthlymask,
                            'units': 'm'
                        }
                        av[name]['datadetails'] = {
                            'name': 'mld',
                            'vars': [
                                'mld',
                                'mask',
                            ],
                            'convert': mldapplymask,
                            'units': 'm'
                        }
                    if name == 'MinMonthlyMLD':
                        av[name]['modeldetails'] = {
                            'name': 'mld',
                            'vars': [
                                'somxl010',
                            ],
                            'convert': mldmonthlymask_min,
                            'units': 'm'
                        }
                        av[name]['datadetails'] = {
                            'name': 'mld',
                            'vars': [
                                'mld',
                                'mask',
                            ],
                            'convert': mldapplymask_min,
                            'units': 'm'
                        }

                    av[name]['layers'] = [
                        'layerless',
                    ]
                    av[name]['regions'] = regionList
                    av[name]['metrics'] = metricList
                    av[name]['datasource'] = 'IFREMER'
                    av[name]['model'] = 'NEMO'
                    av[name]['modelgrid'] = 'eORCA1'
                    av[name]['gridFile'] = paths.orcaGridfn
                    av[name]['Dimensions'] = 2

        if 'FreshwaterFlux' in analysisKeys:

            #ficeberg + friver + fsitherm + pr + prsn - evs

            adds = ['ficeberg', 'friver', 'fsitherm', 'pr', 'prsn']  # - evs

            def calcFreshflux(nc, keys):
                total = -1. * nc.variables['evs'][:]
                for a in adds:
                    total += nc.variables[a][:]
                #a = (nc.variables['SDT__100'][:] +nc.variables['FDT__100'][:])/ (nc.variables['PRD'][:] +nc.variables['PRN'][:] )
                #a = np.ma.masked_where(a>1.01, a)
                return total * 1000000.

            name = 'FreshwaterFlux'
            av[name]['modelFiles'] = listModelDataFiles(
                jobID, 'grid_T', paths.ModelFolder_pref, annual)

            av[name]['dataFile'] = ""
            av[name]['modelcoords'] = medusaCoords
            av[name]['datacoords'] = maredatCoords
            av[name]['modeldetails'] = {
                'name': name,
                'vars':
                ['ficeberg', 'friver', 'fsitherm', 'pr', 'prsn', 'evs'],
                'convert': calcFreshflux,
                'units': 'mg/m2/s'
            }
            av[name]['datadetails'] = {
                'name': '',
                'units': '',
            }
            av[name]['layers'] = [
                'layerless',
            ]  #'100m','200m','Surface - 1000m','Surface - 300m',]#'depthint']

            freshregions = [
                'Global',
                'ignoreInlandSeas',
                'Equator10',
                'AtlanticSOcean',
                'SouthernOcean',
                'ArcticOcean',
                'Remainder',
                'NorthernSubpolarAtlantic',
                'NorthernSubpolarPacific',
            ]
            freshregions.extend(PierceRegions)
            av[name]['regions'] = freshregions
            av[name]['metrics'] = metricList
            av[name]['datasource'] = ''
            av[name]['model'] = 'MEDUSA'
            av[name]['modelgrid'] = 'eORCA1'
            av[name]['gridFile'] = paths.orcaGridfn
            av[name]['Dimensions'] = 2

        if 'ZonalCurrent' in analysisKeys:
            name = 'ZonalCurrent'
            av[name]['modelFiles'] = listModelDataFiles(
                jobID, 'grid_U', paths.ModelFolder_pref, annual)
            if annual:
                av[name]['dataFile'] = paths.GODASFolder + 'ucur.clim.nc'

            av[name]['modelcoords'] = medusaUCoords
            av[name]['datacoords'] = godasCoords

            av[name]['modeldetails'] = {
                'name': name,
                'vars': [
                    'vozocrtx',
                ],
                'convert': ukp.mul1000,
                'units': 'mm/s'
            }
            av[name]['datadetails'] = {
                'name': name,
                'vars': [
                    'ucur',
                ],
                'convert': ukp.NoChange,
                'units': 'mm/s'
            }

            av[name]['layers'] = layerList
            av[name]['regions'] = regionList
            av[name]['metrics'] = metricList

            av[name]['datasource'] = 'GODAS'
            av[name]['model'] = 'NEMO'

            av[name]['modelgrid'] = 'eORCA1'
            av[name]['gridFile'] = paths.orcaGridfn
            av[name]['Dimensions'] = 3

        if 'MeridionalCurrent' in analysisKeys:
            name = 'MeridionalCurrent'
            av[name]['modelFiles'] = listModelDataFiles(
                jobID, 'grid_V', paths.ModelFolder_pref, annual)
            if annual:
                av[name]['dataFile'] = paths.GODASFolder + 'vcur.clim.nc'

            av[name]['modelcoords'] = medusaVCoords
            av[name]['datacoords'] = godasCoords

            av[name]['modeldetails'] = {
                'name': name,
                'vars': [
                    'vomecrty',
                ],
                'convert': ukp.mul1000,
                'units': 'mm/s'
            }
            av[name]['datadetails'] = {
                'name': name,
                'vars': [
                    'vcur',
                ],
                'convert': ukp.NoChange,
                'units': 'mm/s'
            }

            av[name]['layers'] = layerList
            av[name]['regions'] = regionList
            av[name]['metrics'] = metricList

            av[name]['datasource'] = 'GODAS'
            av[name]['model'] = 'NEMO'

            av[name]['modelgrid'] = 'eORCA1'
            av[name]['gridFile'] = paths.orcaGridfn
            av[name]['Dimensions'] = 3

        if 'VerticalCurrent' in analysisKeys:
            name = 'VerticalCurrent'
            av[name]['modelFiles'] = listModelDataFiles(
                jobID, 'grid_W', paths.ModelFolder_pref, annual)
            if annual:
                av[name]['dataFile'] = paths.GODASFolder + 'dzdt.clim.nc'

            av[name]['modelcoords'] = medusaWCoords
            av[name]['datacoords'] = godasCoords

            av[name]['modeldetails'] = {
                'name': name,
                'vars': [
                    'vovecrtz',
                ],
                'convert': ukp.mul1000000,
                'units': 'um/s'
            }
            av[name]['datadetails'] = {
                'name': name,
                'vars': [
                    'dzdt',
                ],
                'convert': ukp.NoChange,
                'units': 'um/s'
            }

            av[name]['layers'] = layerList
            av[name]['regions'] = regionList
            av[name]['metrics'] = metricList

            av[name]['datasource'] = 'GODAS'
            av[name]['model'] = 'NEMO'

            av[name]['modelgrid'] = 'eORCA1'
            av[name]['gridFile'] = paths.orcaGridfn
            av[name]['Dimensions'] = 3

        if 'DMS' in analysisKeys:
            name = 'DMS'
            av[name]['modelFiles'] = listModelDataFiles(
                jobID, 'diad_T', paths.ModelFolder_pref, annual)  #[::30]
            if annual:
                av[name]['dataFile'] = paths.DMSDir + 'DMSclim_mean.nc'
            else:
                av[name]['dataFile'] = ''

            av[name]['modelcoords'] = medusaCoords
            av[name]['datacoords'] = dmsCoords

            av[name]['modeldetails'] = {
                'name': name,
                'vars': [
                    'DMS_ARAN',
                ],
                'convert': ukp.mul1000000,
                'units': 'umol/m3'
            }
            av[name]['datadetails'] = {
                'name': name,
                'vars': [
                    'DMS',
                ],
                'convert': ukp.NoChange,
                'units': 'umol/m3'
            }

            av[name]['layers'] = [
                'layerless',
            ]
            av[name]['regions'] = regionList
            av[name]['metrics'] = metricList

            av[name]['datasource'] = 'Lana'
            av[name]['model'] = 'MEDUSA'

            av[name]['modelgrid'] = 'eORCA1'
            av[name]['gridFile'] = paths.orcaGridfn
            av[name]['Dimensions'] = 2


        #####
        # North Atlantic Salinity
        #sowaflup = "Net Upward Water Flux" ;
        #sohefldo = "Net Downward Heat Flux" ;
        #sofmflup = "Water flux due to freezing/melting" ;
        #sosfldow = "Downward salt flux" ;

        if 'HeatFlux' in analysisKeys:
            name = 'HeatFlux'
            av[name]['modelFiles'] = listModelDataFiles(
                jobID, 'ptrc_T', paths.ModelFolder_pref, annual)
            av[name]['dataFile'] = ''
            av[name]['modelcoords'] = medusaCoords
            av[name]['datacoords'] = takahashiCoords
            av[name]['modeldetails'] = {
                'name': 'HeatFlux',
                'vars': [
                    'hfds',
                ],
                'convert': ukp.NoChange,
                'units': 'W/m2'
            }
            av[name]['datadetails'] = {'name': '', 'units': ''}
            av[name]['layers'] = [
                'layerless',
            ]
            av[name]['regions'] = regionList
            av[name]['metrics'] = metricList
            av[name]['datasource'] = ''
            av[name]['model'] = 'MEDUSA'
            av[name]['modelgrid'] = 'eORCA1'
            av[name]['gridFile'] = paths.orcaGridfn
            av[name]['Dimensions'] = 2

        if 'TotalHeatFlux' in analysisKeys:
            name = 'TotalHeatFlux'
            nc = dataset(paths.orcaGridfn, 'r')
            try:
                ncarea = nc.variables['area'][:]
                surfmask = nc.variables['tmask'][:]
            except:
                surfmask = nc.variables['tmask'][0]
                ncarea = nc.variables['e2t'][:] * nc.variables['e1t'][:]
            nc.close()

            def areatotal(nc, keys):
                if area in list(nc.variables.keys()):
                    area = nc.variables['area'][:]
                else:
                    area = ncarea
                flux = np.ma.array(nc.variables[keys[0]][:].squeeze()) * ncarea
                flux = np.ma.masked_where((surfmask == 0) + (flux.mask), flux)
                return flux.sum() * 1e-12

            av[name]['modelcoords'] = medusaCoords
            av[name]['modelFiles'] = listModelDataFiles(
                jobID, 'diad_T', paths.ModelFolder_pref, annual)
            av[name]['modeldetails'] = {
                'name': name,
                'vars': [
                    'hfds',
                ],
                'convert': areatotal,
                'units': 'TW'
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
            av[name]['model'] = 'NEMO'
            av[name]['modelgrid'] = 'eORCA1'
            av[name]['gridFile'] = paths.orcaGridfn
            av[name]['Dimensions'] = 2
            av[name]['datacoords'] = {'name': '', 'units': ''}
            av[name]['datadetails'] = {'name': '', 'units': ''}
            av[name]['dataFile'] = ''
            av[name]['datasource'] = ''

        naskeys = [
            'sowaflup', 'sohefldo', 'sofmflup', 'sosfldow', 'soicecov',
            'sossheig'
        ]
        if len(set(naskeys).intersection(set(analysisKeys))):
            for name in naskeys:
                if name not in analysisKeys: continue

                #nc = Dataset(paths.orcaGridfn,'r')
                #area = nc.variables['e2t'][:] * nc.variables['e1t'][:]
                #tmask = nc.variables['tmask'][0,:,:]
                #lat = nc.variables['nav_lat'][:,:]
                #nc.close()

                nas_files = listModelDataFiles(jobID, 'grid_T',
                                               paths.ModelFolder_pref, annual)
                try:
                    nc = Dataset(nas_files[0], 'r')
                except:
                    print("nc does not exist:", name)
                    continue

                if name not in list(nc.variables.keys()):
                    print("analysis_timeseries.py:\tWARNING: ", name,
                          "is not in the model file.")
                    continue
                av[name]['modelFiles'] = nas_files
                av[name]['dataFile'] = ''

                av[name]['modelcoords'] = medusaCoords
                av[name]['datacoords'] = medusaCoords

                nasUnits = {
                    'sowaflup': "kg/m2/s",
                    'sohefldo': "W/m2",
                    'sofmflup': "kg/m2/s",
                    'sosfldow': "PSU/m2/s",
                    'soicecov': '',
                    'sossheig': 'm',
                }

                av[name]['modeldetails'] = {
                    'name': name[:],
                    'vars': [
                        name[:],
                    ],
                    'convert': ukp.NoChange,
                    'units': nasUnits[name][:]
                }

                av[name]['regions'] = [
                    'NordicSea',
                    'LabradorSea',
                    'NorwegianSea',
                    'Global',
                ]

                av[name]['datadetails'] = {
                    'name': '',
                    'units': '',
                }
                av[name]['layers'] = [
                    'layerless',
                ]
                av[name]['metrics'] = metricList
                av[name]['datasource'] = ''
                av[name]['model'] = 'NEMO'
                av[name]['modelgrid'] = 'eORCA1'
                av[name]['gridFile'] = paths.orcaGridfn
                av[name]['Dimensions'] = 2

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
                    av[name]['modelFiles'], jobID)
                if strictFileCheck: assert 0

            modelfilesexists = [
                os.path.exists(f) for f in av[name]['modelFiles']
            ]
            if False in modelfilesexists:
                print(
                    "analysis-Timeseries.py:\tWARNING:\tnot model files do not all exist:",
                    av[name]['modelFiles'])
                if strictFileCheck: assert 0

            if 'dataFile' in av[name]:
                if not os.path.exists(av[name]['dataFile']):
                    print(
                        "analysis-Timeseries.py:\tWARNING:\tdata file is not found:",
                        av[name]['dataFile'])
                    if strictFileCheck: assert 0

            #####
            # time series and traffic lights.
            tsa = timeseriesAnalysis(
                av[name]['modelFiles'],
                av[name].get('dataFile', None),
                dataType=name,
                modelcoords=av[name]['modelcoords'],
                modeldetails=av[name]['modeldetails'],
                datacoords=av[name].get('datacoords', None),
                datadetails=av[name].get('datadetails', None),
                datasource=av[name].get('datasource', None),
                model=av[name].get('model', None),
                jobID=jobID,
                layers=av[name]['layers'],
                regions=av[name]['regions'],
                metrics=av[name]['metrics'],
                workingDir=shelvedir,
                imageDir=imagedir,
                grid=av[name]['modelgrid'],
                gridFile=av[name]['gridFile'],
                clean=False,
                noNewFiles=True,
            )
            #dataD[(jobID,name )] = tsa.dataD
            modeldataD[(jobID, name)] = tsa.modeldataD

    #####
    # Data now loaded, making plots next:
    for k in list(modeldataD.keys()):
        print("Model Data D:", k)

    ####
    for name in av.keys():
#   for name in [
#           'Temperature', 'Salinity', 'MLD', 'FreshwaterFlux',
#           'AirSeaFluxCO2', 'AirSeaFlux', 'Chlorophyll', 'Nitrate',
#           'Alkalinity', 'pH'
#   ]:
#       if name not in list(av.keys()): continue
        regions = av[name]['regions']
        layers = av[name]['layers']
        metrics = av[name]['metrics']
        for region, layer, metric  in itertools.product(regions, layers, metrics):
            timesD = {}
            arrD = {}
            for jobID in jobs:
                try:
                    mdata = modeldataD[(jobID, name)][(region, layer, metric)]
                except:
                    continue
                title = titleify([region, layer, metric, name])

                    #timesD[jobID] 	= sorted(mdata.keys())
                    #arrD[jobID]	= [mdata[t] for t in timesD[jobID]]
                
                times, datas = apply_shifttimes(mdata, jobID, shifttimes)
                print('post apply_shifttimes:', len(times), len(datas))
                times, datas = apply_timerange(times, datas, jobID, timeranges)
                timesD[jobID] = times  #mdata.keys())
                arrD[jobID] = datas  #t] for t in timesD[jobID]]
                print(jobID, region, layer, metric, len(times), len(datas))
            timesD, arrD = build_ensemble(timesD, arrD, ensembles)

            if len(list(arrD.keys())) == 0: 
                continue
            units = av[name]['modeldetails']['units']

            ts = 'Together'
            for ls in ['DataOnly', ]: #  'movingav30years']
                tsp.multitimeseries(
                    timesD,  # model times (in floats)
                    arrD,  # model time series
                    data=-999,  # in situ data distribution
                    title=title,
                    filename=ukp.folder(imageFolder) +
                        '_'.join([name, region, layer, ts, ls + '.png']),
                    units=units,
                    plotStyle=ts,
                    lineStyle=ls,
                    colours=colours,
                    thicknesses=lineThicknesses,
                    linestyles=linestyles,
                )

    # Generate a list of comparison images:
    method_images = 'oswalk'
    AllImages = []
    if method_images == 'glob':
        AllImages = glob(imageFolder, recursive=True)
        print('AllImages:','glob', AllImages)
    elif method_images == 'oswalk':
        for root, dirnames, filenames in os.walk(imageFolder):
            for filename in fnmatch.filter(filenames, '*.png'):
                AllImages.append(os.path.join(root, filename))
                print('AllImages:','fors', root, dirnames, filenames, filename)

    if ensembles != {}:
        jobs = list(ensembles.keys())

    # Senmd everything to the comparison maker:
    comparehtml5Maker(
        jobIDs=jobs,
        reportdir=ukp.folder('CompareReports2/' + analysisname),
        files=AllImages,
        clean=False,
        doZip=False,
        jobDescriptions=jobDescriptions,
        jobColours=colours,
        paths=paths,
    )


def flatten(lats, lons, dataA, dataB):
    m = np.ma.array(dataA).mask
    m += np.ma.array(dataB).mask
    m += np.ma.masked_invalid(dataA / dataB).mask

    return  np.ma.masked_where(m, lats).compressed(),\
     np.ma.masked_where(m, lons).compressed(),\
     np.ma.masked_where(m, dataA).compressed(),\
     np.ma.masked_where(m, dataB).compressed()


def load_comparison_yml(master_compare_yml_fn):
    """
    Load the config yaml.
    Takes a file path string
    Returns:
        Details dict.
    """
    with open(master_compare_yml_fn, 'r') as openfile:
        input_yml_dict = yaml.safe_load(openfile)

    if not input_yml_dict or not isinstance(input_yml_dict, dict):
        print(f"Configuration file {master_compare_yml_fn} "
              "is either empty or corrupt, please check its contents")
        sys.exit(1)

    details = {}
    details['name'] = input_yml_dict.get('name', False)
    details['jobs'] = input_yml_dict.get('jobs', False)

    if not details['name']:
        print('Please provide a name for your analysis. In your yaml, this is:')
        print('name: MyAnalysisName')
        sys.exit(0)
    if not details['jobs']:
        print('Please provide at least one JobID for your analysis. In your yaml, this is:')
        print('jobs: ')
        print('    u-ab123:')
        print('        description: "Job descrition"')
        print('        colour: "red"')
        print('        thickness: 0.7')
        print("        linestyle: '-'")
        print('        shifttime: 0.')
        print('        timerange: [1950, 2000]')
        sys.exit(0)

    details['do_analysis_timeseries'] = input_yml_dict.get('do_analysis_timeseries', False)
    details['do_mass_download'] = input_yml_dict.get('do_mass_download', False)
    details['master_suites'] = input_yml_dict.get('master_suites', [])

    # auto download, can differ for each job.
    auto_download = input_yml_dict.get('auto_download', True)
    auto_download_dict = {jobID: auto_download for jobID in details['jobs'].keys()}

    default_thickness = 0.7
    default_linestyle = 'solid'
    default_suite = 'kmf'

    thicknesses = {}
    linestyles = {}
    colours = {}
    suites = {}
    descriptions = {}
    shifttimes = {} # number of years to shift time axis.
    timeranges = {}

    for jobID, job_dict in details['jobs'].items():
        if job_dict.get('colour', False):
            colours[jobID] = job_dict['colour']
        else:
            colours[jobID] = ''.join(['#', "%06x" % random.randint(0, 0xFFFFFF)])
            print('WARNING: No colour provided, setting to random hex colour:', colours[jobID])

        descriptions[jobID] = job_dict.get('description', '')
        thicknesses[jobID] = job_dict.get('thickness', default_thickness)
        linestyles[jobID] = job_dict.get('linestyle', default_linestyle)
        shifttimes[jobID] = float(job_dict.get('shifttime', 0.))
        timeranges[jobID] = job_dict.get('timerange', None)
        suites[jobID] = job_dict.get('suite', default_suite)
        auto_download_dict[jobID] = job_dict.get('auto_download', auto_download_dict[jobID]) 

    details['colours'] = colours
    details['descriptions'] = descriptions
    details['thicknesses'] = thicknesses
    details['linestyles'] = linestyles
    details['shifttimes'] = shifttimes
    details['timeranges'] = timeranges
    details['suites'] = suites
    details['auto_download'] = auto_download_dict
#    print(details)
#    assert 0
    return details


def load_yml_and_run(compare_yml, config_user):
    """
    Loads the comparison yaml file and run compare_yml.

    """
    # Below here is analysis
    details = load_comparison_yml(compare_yml)

    jobs = details['jobs']
    analysis_name = details['name']
    do_analysis_timeseries = details['do_analysis_timeseries']
    do_mass_download = details['do_mass_download']
    master_suites = details['master_suites']

    colours = details['colours']
    thicknesses = details['thicknesses']
    linestyles = details['linestyles']
    descriptions = details['descriptions']
    shifttimes = details['shifttimes']
    timeranges = details['timeranges']
    suites = details['suites']
    auto_download = details['auto_download']

    print('---------------------')
    print('timeseries_compare:',  analysis_name)
    print('job ids:', jobs.keys())
    for jobID in jobs:
        print(jobID, 'description:',descriptions[jobID])
        print(jobID, 'colour:',colours[jobID])
        print(jobID, 'line thickness & style:',thicknesses[jobID], linestyles[jobID])
        print(jobID, 'Shift time by', shifttimes[jobID])
        print(jobID, 'Time range (None means all):', timeranges.get(jobID, None))
        print(jobID, 'suite:', suites[jobID])
        print(jobID, 'auto_download', auto_download[jobID])

    for jobID in jobs:
        # even if you don't want to download, we run this
        # as it clears up the path and ensures recently downloed data is
        # correctly symlinked.
        download_from_mass(jobID, doMoo=do_mass_download, auto_download=auto_download[jobID], config_user=config_user)

    if do_analysis_timeseries:
        for jobID in jobs:
            analysis_timeseries(
                jobID=jobID,
                suites=suites[jobID],
                config_user=config_user
            )

    # Master suite leys:
    if not master_suites:
        master_suites=['physics', 'bgc']  # Defaults

    # make sure its a list:
    if isinstance(master_suites, list) :
        master_suites = [m.lower() for m in master_suites]
    if isinstance(master_suites, str):
        master_suites = master_suites.lower()
        for split_char in [' ', ',', ':']:
            master_suites = master_suites.replace(split_char, ';')
        master_suites = master_suites.split(';')

    timeseries_compare(
        jobs,
        colours = colours,
        suites = master_suites,
        shifttimes=shifttimes,
        timeranges=timeranges,
        jobDescriptions=descriptions,
        analysisname=analysis_name,
        lineThicknesses=thicknesses,
        linestyles=linestyles,
        config_user=config_user
    )


def get_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-y',
                        '--compare_yml',
                        nargs='+',
                        type=str,
                        help='One or more Comparison Analysis configuration file, for examples see bgcval2 input_yml directory.',
                        required=True,
                        )

    parser.add_argument('-c',
                        '--config-file',
                        default=os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                             'default-bgcval2-config.yml'),
                        help='User configuration file (for paths).',
                        required=False)

    args = parser.parse_args()

    return args


def main():
    """Run the main routine."""
    args = get_args()

    # This has a sensible default value.
    config_user=args.config_file

    # This shouldn't fail as it's a required argument.
    compare_ymls = args.compare_yml

    for compare_yml in compare_ymls:
        print(f"analysis_timeseries: Comparison config file {compare_yml}")

        if not os.path.isfile(compare_yml):
            print(f"analysis_timeseries: Could not find comparison config file {compare_yml}")
            sys.exit(1)

        load_yml_and_run(compare_yml, config_user)

    print("Finished... ")


if __name__ == "__main__":
    main()
