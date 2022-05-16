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
   :synopsis: A script to produce an intercomparison of multiple runs the time series analyses.
.. moduleauthor:: Lee de Mora <ledm@pml.ac.uk>

"""

import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')

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
import os, sys, fnmatch
from getpass import getuser
from collections import defaultdict
import yaml

#####
# Load specific local code:
from . import UKESMpython as ukp
from .timeseries import timeseriesAnalysis
from .timeseries import profileAnalysis
from .timeseries import timeseriesPlots as tsp
try:
    from .bgcvaltools.pftnames import getLongName
except:
    from .pftnames import getLongName
from .bgcvaltools.mergeMonthlyFiles import mergeMonthlyFiles, meanDJF
from .netcdf_manipulation.alwaysInclude import alwaysInclude
from .makeReport import comparehtml5Maker
#'from .Paths import paths

from .comparison.shifttimes import shifttimes
from .comparison.ensembles import build_ensemble
from .config.configToDict import configToDict
from .bgcvaltools.dataset import dataset
#from ._runtime_config import _read_yaml
from ._runtime_config import get_run_configuration



#####
# User defined set of paths pointing towards the datasets.
from .Paths.paths import paths_setter


def titleify(ls):
    return ' '.join([getLongName(i) for i in ls])


def listModelDataFiles(jobID, filekey, datafolder, annual, year=''):
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
            print((datafolder + jobID + "/" + jobID + "o_1y_*" + year +
                   "????_" + filekey + ".nc"))
            return sorted(
                glob(datafolder + jobID + "/" + jobID + "o_1y_*" + year +
                     "????_" + filekey + ".nc"))
        else:
            return sorted(
                glob(datafolder + jobID + "/" + jobID + "o_1m_*" + year +
                     "????_" + filekey + ".nc"))


def timeseries_compare(colours,
                       physics=True,
                       bio=False,
                       debug=False,
                       year0=defaultdict(lambda: 0),
                       analysisname='',
                       jobDescriptions=defaultdict(lambda: ''),
                       lineThicknesses=defaultdict(lambda: 1),
                       linestyles=defaultdict(lambda: '-'),
                       ensembles={},#
                       config_user=''):
    ### strategy here is a simple wrapper.
    # It's a little cheat-y, as I'm copying straight from analysis_timeseries.py


    # get runtime configuration
    if config_user:
        paths_dict, config_user = get_run_configuration(config_user)
    else:
        paths_dict, config_user = get_run_configuration("defaults")

    # filter paths dict into an object that's usable below
    paths = paths_setter(paths_dict)

    jobs = sorted(colours.keys())
    for ensemble in list(ensembles.keys()):
        # ensembles names can not be the same as jobIDs
        jobs.remove(ensemble)
    
    if analysisname == '':
        imageFolder = paths.imagedir + '/TimeseriesCompare/'
        if len(jobs) == 1: imageFolder += jobs[0]
        elif len(jobs) == 2: imageFolder += jobs[0] + '_and_' + jobs[1]
        else: imageFolder += str(len(jobs)) + 'Runs_' + jobs[0]
    else:
        imageFolder = paths.imagedir + '/TimeseriesCompare/' + analysisname
    if debug:
        imageFolder = imageFolder + '_debug'
        analysisname = analysisname + '_debug'

    annual = True
    strictFileCheck = False

    analysisKeys = []

    if physics:
        analysisKeys.append('DrakePassageTransport')  # DrakePassageTransport
        analysisKeys.append('TotalIceArea')  # TotalIceArea
        analysisKeys.append('NorthernTotalIceArea')  # North TotalIceArea
        analysisKeys.append('SouthernTotalIceArea')  # South TotalIceArea
        analysisKeys.append('TotalIceExtent')  # work in progress
        analysisKeys.append('NorthernTotalIceExtent')  # work in progress
        analysisKeys.append('SouthernTotalIceExtent')  # work in progress
        analysisKeys.append('WeddelIceExent')  # work in progress
        #analysisKeys.append('NorthernMIZArea')
        #analysisKeys.append('SouthernMIZArea')
        #analysisKeys.append('TotalMIZArea')
        analysisKeys.append('NorthernMIZfraction')
        analysisKeys.append('SouthernMIZfraction')
        analysisKeys.append('TotalMIZfraction')

        analysisKeys.append('AMOC_26N')
        analysisKeys.append('AMOC_32S')
        analysisKeys.append('ADRC_26N')  # AMOC 26N
        analysisKeys.append('Temperature')  # WOA Temperature
        analysisKeys.append('Salinity')  # WOA Salinity
        analysisKeys.append('MLD')  # MLD
        analysisKeys.append('MaxMonthlyMLD')  # MLD Monthly max
        analysisKeys.append('MinMonthlyMLD')  # MLD Monthly min

        analysisKeys.append('ZonalCurrent')  # Zonal Veloctity
        analysisKeys.append('MeridionalCurrent')  # Meridional Veloctity
        analysisKeys.append('VerticalCurrent')  # Vertical Veloctity
        analysisKeys.append('GlobalMeanTemperature')  # Global Mean Temperature
        analysisKeys.append('VolumeMeanTemperature')  # Global Mean Temperature
        analysisKeys.append('GlobalMeanSalinity')  # Global Mean Salinity
        analysisKeys.append(
            'IcelessMeanSST')  # Global Mean Surface Temperature with no ice

        analysisKeys.append('sowaflup')  # Net Upward Water Flux
        analysisKeys.append('sohefldo')  # Net downward Water Flux
        analysisKeys.append('sofmflup')  # Water flux due to freezing/melting
        analysisKeys.append('sosfldow')  # Downward salt flux
        analysisKeys.append('soicecov')  # Ice fraction
        analysisKeys.append('sossheig')  # SSH
        analysisKeys.append('FreshwaterFlux')  # Fresh water flux
        analysisKeys.append('HeatFlux')
        analysisKeys.append('TotalHeatFlux')
        analysisKeys.append('scvoltot')
        analysisKeys.append('soga')
        analysisKeys.append('thetaoga')
        analysisKeys.append('scalarHeatContent')

    if bio:
        analysisKeys.append('TotalAirSeaFluxCO2')  # work in progress
        analysisKeys.append('NoCaspianAirSeaFluxCO2')  # work in progress

        analysisKeys.append('AirSeaFlux')  # work in progress
        analysisKeys.append('IntPP_OSU')  # OSU Integrated primpary production
        analysisKeys.append('GlobalExportRatio')

        analysisKeys.append('N')  # WOA Nitrate
        analysisKeys.append('Si')  # WOA Siliate
        analysisKeys.append('O2')  # WOA Oxygen
        analysisKeys.append('Iron')
        analysisKeys.append('Alk')
        analysisKeys.append('DIC')
        analysisKeys.append('pH')

        analysisKeys.append('CHD')
        analysisKeys.append('CHN')
        analysisKeys.append('CHL')
        analysisKeys.append('DiaFrac')
        analysisKeys.append('DMS')
        analysisKeys.append('Dust')  # Dust
        analysisKeys.append('TotalDust')  # Total Dust

        #                analysisKeys.append('DMS_ARAN')

        analysisKeys.append('DTC')

        analysisKeys.append(
            'TotalOMZVolume')  # Total Oxygen Minimum zone Volume
        analysisKeys.append('OMZThickness')  # Oxygen Minimum Zone Thickness
        analysisKeys.append('OMZMeanDepth')  # Oxygen Minimum Zone mean depth
        analysisKeys.append('VolumeMeanOxygen')  # Volume weighted mean Oxygen

    if debug:
        ####
        # Supercedes other flags.
        analysisKeys = []
        #		analysisKeys.append('ERSST')
        analysisKeys.append('DrakePassageTransport')  # DrakePassageTransport
        #                analysisKeys.append('VolumeMeanOxygen')         # Volume weighted mean Oxygen
        analysisKeys.append('AMOC_26N')
        #                analysisKeys.append('NoCaspianAirSeaFluxCO2')   # work in progress
        analysisKeys.append('VolumeMeanTemperature')  # Global Mean Temperature
        #		analysisKeys.append('CHD')
        #		analysisKeys.append('CHN')
        #		analysisKeys.append('DiaFrac')
        #                analysisKeys.append('AMOC_26N')
        #                analysisKeys.append('MLD')
        #                analysisKeys.append('Temperature')             # WOA Temperature
        #                analysisKeys.append('Salinity')                # WOA Salinity
        #                analysisKeys.append('DMS')
        #                analysisKeys.append('N')                        # WOA Nitrate
        #                analysisKeys.append('Si')                       # WOA Siliate

        #                analysisKeys.append('TotalAirSeaFluxCO2')          # work in progress
        #                analysisKeys.append('AirSeaFlux')          # work in progress
        #                analysisKeys.append('AirSeaFluxCO2')          # work in progress

        #                analysisKeys.append('scvoltot')
        #                analysisKeys.append('soga')
        #                analysisKeys.append('thetaoga')
        #                analysisKeys.append('scalarHeatContent')

        #analysisKeys.append('ADRC_26N')                # AMOC 26N
        #                analysisKeys.append('VerticalCurrent')          # Vertical Veloctity
        #                analysisKeys.append('sossheig')                 # SSH
        #                analysisKeys.append('NorthernTotalIceArea')     # North TotalIceArea
        #                analysisKeys.append('SouthernTotalIceArea')     # South TotalIceArea
        #                analysisKeys.append('TotalIceArea')     	#  TotalIceArea
        #                analysisKeys.append('NorthernTotalIceExtent')   # work in progress
        #                analysisKeys.append('SouthernTotalIceExtent')   # work in progress
        #                analysisKeys.append('WeddelIceExent')   # work in progress
        #                analysisKeys.append('TotalIceExtent')           # work in progress
        #                analysisKeys.append('NorthernMIZArea')
        #                analysisKeys.append('SouthernMIZArea')
        #                analysisKeys.append('TotalMIZArea')

        #                        analysisKeys.append('NorthernMIZArea')
        #                        analysisKeys.append('SouthernMIZArea')
        #                        analysisKeys.append('TotalMIZArea')
        #                analysisKeys.append('NorthernMIZfraction')
        #                analysisKeys.append('SouthernMIZfraction')
        #                analysisKeys.append('TotalMIZfraction')

        #		analysisKeys.append('FreshwaterFlux')		# Fresh water flux
        analysisKeys.append('GlobalMeanTemperature')
#                analysisKeys.append('GlobalMeanSalinity')

#                analysisKeys.append('HeatFlux')
#                analysisKeys.append('TotalHeatFlux')

#               	analysisKeys.append('quickSST')    		# Area Weighted Mean Surface Temperature
#       	  	analysisKeys.append('TotalOMZVolume')           # Total Oxygen Minimum zone Volume
#       	 	analysisKeys.append('OMZThickness')             # Oxygen Minimum Zone Thickness
#        	analysisKeys.append('OMZMeanDepth')             # Oxygen Minimum Zone mean depth
#		analysisKeys.append('O2')                       # WOA Oxygen
#       	if bio ==False:return
#       	if physics == True:return

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

    medusaCoords = {
        't': 'index_t',
        'z': 'deptht',
        'lat': 'nav_lat',
        'lon': 'nav_lon',
        'cal': '360_day',
    }  # model doesn't need time dict.
    medusaUCoords = {
        't': 'index_t',
        'z': 'depthu',
        'lat': 'nav_lat',
        'lon': 'nav_lon',
        'cal': '360_day',
    }  # model doesn't need time dict.
    medusaVCoords = {
        't': 'index_t',
        'z': 'depthv',
        'lat': 'nav_lat',
        'lon': 'nav_lon',
        'cal': '360_day',
    }  # model doesn't need time dict.
    medusaWCoords = {
        't': 'index_t',
        'z': 'depthw',
        'lat': 'nav_lat',
        'lon': 'nav_lon',
        'cal': '360_day',
    }  # model doesn't need time dict.

    icCoords = {
        't': 'time_counter',
        'z': 'nav_lev',
        'lat': 'nav_lat',
        'lon': 'nav_lon',
        'cal': '360_day',
    }  # model doesn't need time dict.
    maredatCoords = {
        't': 'index_t',
        'z': 'DEPTH',
        'lat': 'LATITUDE',
        'lon': 'LONGITUDE',
        'cal': 'standard',
        'tdict': ukp.tdicts['ZeroToZero']
    }
    takahashiCoords = {
        't': 'index_t',
        'z': 'index_z',
        'lat': 'LAT',
        'lon': 'LON',
        'cal': 'standard',
        'tdict': ukp.tdicts['ZeroToZero']
    }
    woaCoords = {
        't': 'index_t',
        'z': 'depth',
        'lat': 'lat',
        'lon': 'lon',
        'cal': 'standard',
        'tdict': ukp.tdicts['ZeroToZero']
    }
    glodapCoords = {
        't': 'index_t',
        'z': 'depth',
        'lat': 'latitude',
        'lon': 'longitude',
        'cal': 'standard',
        'tdict': []
    }
    glodapv2Coords = {
        't': 'time',
        'z': 'Pressure',
        'lat': 'lat',
        'lon': 'lon',
        'cal': '',
        'tdict': {
            0: 0,
        }
    }
    mldCoords = {
        't': 'index_t',
        'z': 'index_z',
        'lat': 'lat',
        'lon': 'lon',
        'cal': 'standard',
        'tdict': ukp.tdicts['ZeroToZero']
    }
    cciCoords = {
        't': 'index_t',
        'z': 'index_z',
        'lat': 'lat',
        'lon': 'lon',
        'cal': 'standard',
        'tdict': ['ZeroToZero']
    }
    dmsCoords = {
        't': 'time',
        'z': 'depth',
        'lat': 'Latitude',
        'lon': 'Longitude',
        'cal': 'standard',
        'tdict': ukp.tdicts['ZeroToZero']
    }
    godasCoords = {
        't': 'index_t',
        'z': 'level',
        'lat': 'lat',
        'lon': 'lon',
        'cal': 'standard',
        'tdict': ['ZeroToZero']
    }

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
        if 'DrakePassageTransport' in analysisKeys:
            name = 'DrakePassageTransport'
            ####
            # Note that this will only work with the eORCA1grid.

            # coordinates of Drake Passage
            LON = 219
            LAT0 = 79
            LAT1 = 109

            nc = Dataset(orcaGridfn, 'r')
            e2u = nc.variables['e2u'][LAT0:LAT1, LON]
            umask = nc.variables['umask'][:, LAT0:LAT1, LON]
            nc.close()

            def drake(nc, keys):
                e3u = nc.variables['e3u'][0, :, LAT0:LAT1, LON]
                velo = nc.variables['vozocrtx'][0, :, LAT0:LAT1, LON]
                return np.sum(velo * e3u * e2u * umask) * 1.e-6

            av[name]['modelFiles'] = listModelDataFiles(
                jobID, 'grid_U', paths.ModelFolder_pref, annual)
            av[name]['dataFile'] = ''

            av[name]['modelcoords'] = medusaCoords
            av[name]['datacoords'] = medusaCoords

            av[name]['modeldetails'] = {
                'name': name,
                'vars': [
                    'e3u',
                    'vozocrtx',
                ],
                'convert': drake,
                'units': 'Sv'
            }

            av[name]['regions'] = [
                'regionless',
            ]
            av[name]['datadetails'] = {
                'name': '',
                'units': '',
            }
            av[name]['layers'] = [
                'layerless',
            ]
            av[name]['metrics'] = [
                'metricless',
            ]
            av[name]['datasource'] = ''
            av[name]['model'] = 'NEMO'
            av[name]['modelgrid'] = 'eORCA1'
            av[name]['gridFile'] = orcaGridfn
            av[name]['Dimensions'] = 1

        icekeys = [
            'NorthernTotalIceArea', 'SouthernTotalIceArea',
            'WeddelTotalIceArea', 'TotalIceArea', 'NorthernTotalIceExtent',
            'WeddelIceExent', 'SouthernTotalIceExtent', 'TotalIceExtent',
            'NorthernMIZArea', 'SouthernMIZArea', 'TotalMIZArea',
            'NorthernMIZfraction', 'SouthernMIZfraction', 'TotalMIZfraction'
        ]
        if len(set(icekeys).intersection(set(analysisKeys))):
            for name in icekeys:
                if name not in analysisKeys: continue

                nc = dataset(paths.orcaGridfn, 'r')
                area = nc.variables['e2t'][:] * nc.variables['e1t'][:]
                tmask = nc.variables['tmask'][0, :, :]
                lat = nc.variables['nav_lat'][:, :]
                lon = nc.variables['nav_lon'][:, :]
                nc.close()

                def calcTotalIceArea(nc, keys):  #Global
                    arr = nc.variables[keys[0]][:].squeeze() * area
                    return np.ma.masked_where(tmask == 0, arr).sum() / 1E12

                def calcTotalIceAreaN(nc, keys):  # North
                    arr = nc.variables[keys[0]][:].squeeze() * area
                    return np.ma.masked_where(
                        (tmask == 0) + (lat < 0.), arr).sum() / 1E12

                def calcTotalIceAreaS(nc, keys):  # South
                    arr = nc.variables[keys[0]][:].squeeze() * area
                    return np.ma.masked_where(
                        (tmask == 0) + (lat > 0.), arr).sum() / 1E12

                def calcTotalIceAreaWS(nc, keys):
                    arr = nc.variables[keys[0]][:].squeeze() * area
                    return np.ma.masked_where(
                        (tmask == 0) + weddelmask, arr).sum() / 1E12

                def calcMIZArea(nc, keys):  #Global
                    arr = nc.variables[keys[0]][:].squeeze()
                    return np.ma.masked_where(
                        tmask == 0 + (arr < 0.15) +
                        (arr > 0.80), arr * area).sum() / 1E12

                def calcMIZAreaN(nc, keys):  # North
                    arr = nc.variables[keys[0]][:].squeeze()
                    return np.ma.masked_where(
                        (tmask == 0) + (lat < 0.) + (arr < 0.15) +
                        (arr > 0.80), arr * area).sum() / 1E12

                def calcMIZAreaS(nc, keys):  # South
                    arr = nc.variables[keys[0]][:].squeeze()
                    return np.ma.masked_where(
                        (tmask == 0) + (lat > 0.) + (arr < 0.15) +
                        (arr > 0.80), arr * area).sum() / 1E12

                def calcMIZfraction(nc, keys):  #Global
                    arr = nc.variables[keys[0]][:].squeeze()
                    arr = np.ma.masked_where(tmask == 0 + (arr == 0.), arr)
                    value = np.ma.masked_where(
                        arr.mask + (arr < 0.15) + (arr > 0.80), area).sum()
                    denom = np.ma.masked_where(arr.mask, area).sum()
                    return value / denom

                def calcMIZfractionN(nc, keys):  # North
                    arr = nc.variables[keys[0]][:].squeeze()
                    arr = np.ma.masked_where(
                        tmask == 0 + (arr == 0.) + (lat < 0.), arr)
                    value = np.ma.masked_where(
                        arr.mask + (arr < 0.15) + (arr > 0.80), area).sum()
                    denom = np.ma.masked_where(arr.mask, area).sum()
                    return value / denom

                def calcMIZfractionS(nc, keys):  # South
                    arr = nc.variables[keys[0]][:].squeeze()
                    arr = np.ma.masked_where(
                        tmask == 0 + (arr == 0.) + (lat > 0.), arr)
                    value = np.ma.masked_where(
                        arr.mask + (arr < 0.15) + (arr > 0.80), area).sum()
                    denom = np.ma.masked_where(arr.mask, area).sum()
                    return value / denom

                def calcTotalIceExtent(nc, keys):  #Global
                    return np.ma.masked_where(
                        (tmask == 0) +
                        (nc.variables[keys[0]][:].squeeze() < 0.15),
                        area).sum() / 1E12

                def calcTotalIceExtentN(nc, keys):  # North
                    return np.ma.masked_where(
                        (tmask == 0) +
                        (nc.variables[keys[0]][:].squeeze() < 0.15) +
                        (lat < 0.), area).sum() / 1E12

                def calcTotalIceExtentS(nc, keys):  # South
                    return np.ma.masked_where(
                        (tmask == 0) +
                        (nc.variables[keys[0]][:].squeeze() < 0.15) +
                        (lat > 0.), area).sum() / 1E12

                weddelmask = (lat < -80.) + (lat > -65.) + (lon < -60.) + (
                    lon > -20.)

                def calcTotalIceExtentWS(nc, keys):  # South
                    return np.ma.masked_where(
                        (tmask == 0) +
                        (nc.variables[keys[0]][:].squeeze() < 0.15) +
                        weddelmask, area).sum() / 1E12

                if jobID == 'u-as462monthly':
                    av[name]['modelFiles'] = sorted(
                        glob(
                            '/group_workspaces/jasmin2/ukesm/BGC_data/u-as462/monthly/*.nc'
                        ))
                elif jobID == 'u-ar977monthly':
                    av[name]['modelFiles'] = sorted(
                        glob(
                            '/group_workspaces/jasmin2/ukesm/BGC_data/u-ar977/monthly/*.nc'
                        ))
                else:
                    av[name]['modelFiles'] = listModelDataFiles(
                        jobID, 'grid_T', paths.ModelFolder_pref, annual)
                av[name]['dataFile'] = ''

                av[name]['modelcoords'] = medusaCoords
                av[name]['datacoords'] = medusaCoords

                if name in [
                        'NorthernTotalIceArea',
                ]:
                    av[name]['modeldetails'] = {
                        'name': name,
                        'vars': [
                            'soicecov',
                        ],
                        'convert': calcTotalIceAreaN,
                        'units': '1E6 km^2'
                    }
        #	av[name]['regions'] 		=  ['NorthHemisphere',]

                if name in [
                        'SouthernTotalIceArea',
                ]:
                    av[name]['modeldetails'] = {
                        'name': name,
                        'vars': [
                            'soicecov',
                        ],
                        'convert': calcTotalIceAreaS,
                        'units': '1E6 km^2'
                    }

                if name in [
                        'WeddelTotalIceArea',
                ]:
                    av[name]['modeldetails'] = {
                        'name': name,
                        'vars': [
                            'soicecov',
                        ],
                        'convert': calcTotalIceAreaWS,
                        'units': '1E6 km^2'
                    }

                if name in [
                        'TotalIceArea',
                ]:
                    av[name]['modeldetails'] = {
                        'name': name,
                        'vars': [
                            'soicecov',
                        ],
                        'convert': calcTotalIceArea,
                        'units': '1E6 km^2'
                    }
        #	av[name]['regions'] 		=  ['Global',]

                if name in [
                        'NorthernTotalIceExtent',
                ]:
                    av[name]['modeldetails'] = {
                        'name': name,
                        'vars': [
                            'soicecov',
                        ],
                        'convert': calcTotalIceExtentN,
                        'units': '1E6 km^2'
                    }
        #	av[name]['regions'] 		=  ['NorthHemisphere',]

                if name in [
                        'SouthernTotalIceExtent',
                ]:
                    av[name]['modeldetails'] = {
                        'name': name,
                        'vars': [
                            'soicecov',
                        ],
                        'convert': calcTotalIceExtentS,
                        'units': '1E6 km^2'
                    }
        #	av[name]['regions'] 		=  ['SouthHemisphere',]

                if name in [
                        'WeddelIceExent',
                ]:
                    av[name]['modeldetails'] = {
                        'name': name,
                        'vars': [
                            'soicecov',
                        ],
                        'convert': calcTotalIceExtentWS,
                        'units': '1E6 km^2'
                    }
        #	av[name]['regions'] 		=  ['SouthHemisphere',]

                if name in [
                        'TotalIceExtent',
                ]:
                    av[name]['modeldetails'] = {
                        'name': name,
                        'vars': [
                            'soicecov',
                        ],
                        'convert': calcTotalIceExtent,
                        'units': '1E6 km^2'
                    }
        #	av[name]['regions'] 		=  ['Global',]

                if name in [
                        'NorthernMIZArea',
                ]:
                    av[name]['modeldetails'] = {
                        'name': name,
                        'vars': [
                            'soicecov',
                        ],
                        'convert': calcMIZAreaN,
                        'units': '1E6 km^2'
                    }

                if name in [
                        'SouthernMIZArea',
                ]:
                    av[name]['modeldetails'] = {
                        'name': name,
                        'vars': [
                            'soicecov',
                        ],
                        'convert': calcMIZAreaS,
                        'units': '1E6 km^2'
                    }

                if name in [
                        'TotalMIZArea',
                ]:
                    av[name]['modeldetails'] = {
                        'name': name,
                        'vars': [
                            'soicecov',
                        ],
                        'convert': calcMIZArea,
                        'units': '1E6 km^2'
                    }

                if name in [
                        'NorthernMIZfraction',
                ]:
                    av[name]['modeldetails'] = {
                        'name': name,
                        'vars': [
                            'soicecov',
                        ],
                        'convert': calcMIZfractionN,
                        'units': ''
                    }

                if name in [
                        'SouthernMIZfraction',
                ]:
                    av[name]['modeldetails'] = {
                        'name': name,
                        'vars': [
                            'soicecov',
                        ],
                        'convert': calcMIZfractionS,
                        'units': ''
                    }

                if name in [
                        'TotalMIZfraction',
                ]:
                    av[name]['modeldetails'] = {
                        'name': name,
                        'vars': [
                            'soicecov',
                        ],
                        'convert': calcMIZfraction,
                        'units': ''
                    }

                av[name]['regions'] = [
                    'regionless',
                ]

                av[name]['datadetails'] = {
                    'name': '',
                    'units': '',
                }
                #av[name]['layers'] 		=  ['Surface',]
                av[name]['layers'] = [
                    'layerless',
                ]
                av[name]['metrics'] = [
                    'metricless',
                ]
                av[name]['datasource'] = ''
                av[name]['model'] = 'CICE'
                av[name]['modelgrid'] = 'eORCA1'
                av[name]['gridFile'] = paths.orcaGridfn
                av[name]['Dimensions'] = 1

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

        #####
        # Total
        if 'IntPP_OSU' in analysisKeys:
            noOSU = True
            nc = Dataset(orcaGridfn, 'r')
            area = nc.variables['e1t'][:] * nc.variables['e2t'][:]
            nc.close()

            def medusadepthInt(nc, keys):
                #	 mmolN/m2/d        [mg C /m2/d]   [mgC/m2/yr] [gC/m2/yr]     Gt/m2/yr
                factor = 1. * 6.625 * 12.011 * 365. / 1000. / 1E15
                arr = (nc.variables[keys[0]][:] +
                       nc.variables[keys[1]][:]).squeeze() * factor
                if arr.ndim == 3:
                    for i in np.arange(arr.shape[0]):
                        arr[i] = arr[i] * area
                elif arr.ndim == 2:
                    arr = arr * area
                else:
                    assert 0
                return arr.sum()

            name = 'TotalIntegratedPrimaryProduction'
            if annual:
                av[name]['modelFiles'] = listModelDataFiles(
                    jobID, 'diad_T', paths.ModelFolder_pref, annual)
                if noOSU: av[name]['dataFile'] = ''
                else:
                    av[name][
                        'dataFile'] = paths.OSUDir + "/standard_VGPM.SeaWIFS.global.average.nc"

            av[name]['modelcoords'] = medusaCoords
            av[name]['datacoords'] = glodapCoords

            av[name]['modeldetails'] = {
                'name': 'IntPP',
                'vars': ['PRN', 'PRD'],
                'convert': medusadepthInt,
                'units': 'Gt/yr'
            }
            if noOSU:
                av[name]['datadetails'] = {'name': '', 'units': ''}

            else:
                nc = Dataset(av[name]['dataFile'], 'r')
                lats = nc.variables['latitude'][:]
                osuareas = np.zeros((1080, 2160))
                osuarea = (111100. /
                           6.)**2.  # area of a pixel at equator. in m2
                for a in np.arange(1080):
                    osuareas[a] = np.ones(
                        (2160, )) * osuarea * np.cos(np.deg2rad(lats[a]))

                def osuconvert(nc, keys):
                    arr = nc.variables[keys[0]][:, :, :]
                    tlen = arr.shape[0]
                    arr = arr.sum(0) / tlen * 365. / 1000. / 1E15
                    if arr.ndim == 3:
                        for i in np.arange(arr.shape[0]):
                            arr[i] = arr[i] * osuarea
                    elif arr.ndim == 2:
                        arr = arr * osuarea
                    else:
                        assert 0
                    return arr.sum()

                av[name]['datadetails'] = {
                    'name': 'IntPP',
                    'vars': [
                        'NPP',
                    ],
                    'convert': osuconvert,
                    'units': 'Gt/yr'
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
            if noOSU: av[name]['datasource'] = ''
            else: av[name]['datasource'] = 'OSU'
            av[name]['model'] = 'MEDUSA'
            av[name]['modelgrid'] = 'eORCA1'
            av[name]['gridFile'] = orcaGridfn
            av[name]['Dimensions'] = 1

        if 'TotalAirSeaFluxCO2' in analysisKeys:
            name = 'TotalAirSeaFluxCO2'
            nc = Dataset(orcaGridfn, 'r')
            area = nc.variables['e1t'][:] * nc.variables['e2t'][:]
            nc.close()

            def eOrcaTotal(nc, keys):
                factor = 365.25 * 12. / 1000. / 1.E15
                arr = nc.variables['CO2FLUX'][:].squeeze(
                ) * factor  # mmolC/m2/d
                if arr.ndim == 3:
                    for i in np.arange(arr.shape[0]):
                        arr[i] = arr[i] * area
                elif arr.ndim == 2:
                    arr = arr * area
                else:
                    assert 0
                return arr.sum()

            def takaTotal(nc, keys):
                arr = nc.variables['TFLUXSW06'][:].squeeze(
                )  # 10^12 g Carbon year^-1
                arr = -1.E12 * arr / 1.E15  #/ 365.				#g Carbon/day
                #area = nc.variables['AREA_MKM2'][:].squeeze() *1E12	# 10^6 km^2
                #fluxperarea = arr/area
                return arr.sum()
                # area 10^6 km^2
                # flux:  10^15 g Carbon month^-1. (GT)/m2/month

            av[name]['modelFiles'] = listModelDataFiles(
                jobID, 'diad_T', paths.ModelFolder_pref, annual)
            if annual:
                av[name][
                    'dataFile'] = paths.TakahashiFolder + 'takahashi_2009_Anual_sumflux_2006c_noHead.nc'
            else:
                av[name][
                    'dataFile'] = paths.TakahashiFolder + 'takahashi2009_month_flux_pCO2_2006c_noHead.nc'
                print("Air Sea Flux CO2 monthly not implemented")
                assert 0
                #av[name]['dataFile'] 		=  paths.TakahashiFolder+'takahashi2009_month_flux_pCO2_2006c_noHead.nc'

            av[name]['modelcoords'] = medusaCoords
            av[name]['datacoords'] = takahashiCoords
            av[name]['modeldetails'] = {
                'name': 'AirSeaFluxCO2',
                'vars': [
                    'CO2FLUX',
                ],
                'convert': eOrcaTotal,
                'units': 'Pg C/yr'
            }
            av[name]['datadetails'] = {
                'name': 'AirSeaFluxCO2',
                'vars': ['TFLUXSW06', 'AREA_MKM2'],
                'convert': takaTotal,
                'units': 'Pg C/yr'
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
            av[name]['Dimensions'] = 2

            noTaka = True
            if noTaka:
                av[name]['datadetails'] = {'name': '', 'units': ''}
                av[name]['dataFile'] = ''
                av[name]['datasource'] = ''

        if 'NoCaspianAirSeaFluxCO2' in analysisKeys:
            name = 'NoCaspianAirSeaFluxCO2'
            nc = Dataset(paths.orcaGridfn, 'r')
            area = nc.variables['e1t'][:] * nc.variables['e2t'][:]
            nc.close()

            def eOrcaTotal(nc, keys):
                factor = 365.25 * 12. / 1000. / 1.E15
                arr = nc.variables['CO2FLUX'][:].squeeze(
                ) * factor  # mmolC/m2/d
                if arr.ndim == 3:
                    for i in np.arange(arr.shape[0]):
                        arr[i] = arr[i] * area
                elif arr.ndim == 2:
                    arr = arr * area
                else:
                    assert 0
                return arr

            av[name]['modelFiles'] = listModelDataFiles(
                jobID, 'diad_T', paths.ModelFolder_pref, annual)

            av[name]['modelcoords'] = medusaCoords

            av[name]['modeldetails'] = {
                'name': 'AirSeaFluxCO2',
                'vars': [
                    'CO2FLUX',
                ],
                'convert': eOrcaTotal,
                'units': 'Pg C/yr'
            }

            av[name]['layers'] = [
                'layerless',
            ]
            av[name]['regions'] = [
                'ignoreCaspian',
            ]
            av[name]['metrics'] = [
                'sum',
            ]
            av[name]['model'] = 'MEDUSA'
            av[name]['modelgrid'] = 'eORCA1'
            av[name]['gridFile'] = paths.orcaGridfn
            av[name]['Dimensions'] = 2

            av[name]['datacoords'] = {'name': '', 'units': ''}
            av[name]['datadetails'] = {'name': '', 'units': ''}
            av[name]['dataFile'] = ''
            av[name]['datasource'] = ''

        if 'AirSeaFlux' in analysisKeys:

            #nc = dataset(paths.orcaGridfn,'r')
            #area = nc.variables['e1t'][:]*nc.variables['e2t'][:]
            #nc.close()

            def eOrcaTotal(nc, keys):
                factor = 12. / 1000.
                arr = nc.variables['CO2FLUX'][:].squeeze()  # mmolC/m2/d
                #if arr.ndim ==3:
                #       for i in np.arange(arr.shape[0]):
                #               arr[i] = arr[i]*area
                #elif arr.ndim ==2: arr = arr*area
                #else: assert 0
                return arr * factor

            def takaTotal(nc, keys):
                arr = nc.variables['TFLUXSW06'][:].squeeze(
                )  # 10^12 g Carbon year^-1
                arr = -1.E12 * arr / 365.  #g Carbon/day
                factor = -1.E12 / (365.)  # convert to #/ 1.E12
                area = nc.variables['AREA_MKM2'][:].squeeze(
                ) * 1E12  # 10^6 km^2
                fluxperarea = arr / area
                #arr = arr*area #* 1.E24        # converts area into m^2
                #print arr.sum(), arr.sum()*factor
                return fluxperarea
                # area 10^6 km^2
                # flux:  10^15 g Carbon month^-1. (GT)/m2/month

            name = 'AirSeaFluxCO2'

            av[name]['modelFiles'] = listModelDataFiles(
                jobID, 'diad_T', paths.ModelFolder_pref, annual)
            if annual:
                av[name][
                    'dataFile'] = paths.TakahashiFolder + 'takahashi_2009_Anual_sumflux_2006c_noHead.nc'
            else:
                av[name][
                    'dataFile'] = paths.TakahashiFolder + 'takahashi2009_month_flux_pCO2_2006c_noHead.nc'
                print("Air Sea Flux CO2 monthly not implemented")
                assert 0
                #av[name]['dataFile']           =  paths.TakahashiFolder+'takahashi2009_month_flux_pCO2_2006c_noHead.nc'

            av[name]['modelcoords'] = medusaCoords
            av[name]['datacoords'] = takahashiCoords
            av[name]['modeldetails'] = {
                'name': 'AirSeaFluxCO2',
                'vars': [
                    'CO2FLUX',
                ],
                'convert': eOrcaTotal,
                'units': 'g C/m2/day'
            }
            av[name]['datadetails'] = {
                'name': 'AirSeaFluxCO2',
                'vars': ['TFLUXSW06', 'AREA_MKM2'],
                'convert': takaTotal,
                'units': 'g C/m2/day'
            }
            av[name]['layers'] = [
                'Surface',
            ]
            av[name]['regions'] = regionList
            av[name]['metrics'] = metricList
            av[name]['datasource'] = ''
            av[name]['model'] = 'MEDUSA'
            av[name]['modelgrid'] = 'eORCA1'
            av[name]['gridFile'] = paths.orcaGridfn
            av[name]['Dimensions'] = 2

        if 'pH' in analysisKeys:

            def convertmeqm3TOumolkg(nc, keys):
                return nc.variables[keys[0]][:] * 1.027

            name = 'pH'
            if annual:
                av[name]['modelFiles'] = listModelDataFiles(
                    jobID, 'diad_T', paths.ModelFolder_pref, annual)
                #av[name]['dataFile'] 	=  paths.GlodapDir+'p.nc'
            else:
                print("pH data not available for monthly Analysis")
                assert 0

            av[name]['modelcoords'] = medusaCoords
            av[name]['datacoords'] = {'name': '', 'units': ''}
            av[name]['datadetails'] = {'name': '', 'units': ''}
            av[name]['dataFile'] = ''
            av[name]['datasource'] = ''

            av[name]['modeldetails'] = {
                'name': name,
                'vars': [
                    'PH3',
                ],
                'convert': ukp.NoChange,
                'units': 'pH',
            }
            #av[name]['datadetails']  	= {'name': name, 'vars':['Alk',], 'convert': convertmeqm3TOumolkg,'units':'meq/m^3',}

            av[name]['layers'] = [
                'Surface',
            ]  #'100m','300m','1000m',]
            #	av[name]['regions'] 		= regionList
            #av[name]['layers'] 		= layerList
            if debug:
                av[name]['regions'] = [
                    'Global',
                ]
            else:
                av[name]['regions'] = regionList
            av[name]['metrics'] = metricList

            av[name]['model'] = 'MEDUSA'

            av[name]['modelgrid'] = 'eORCA1'
            av[name]['gridFile'] = paths.orcaGridfn
            av[name]['Dimensions'] = 3

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

        if 'CHL' in analysisKeys:
            name = 'Chlorophyll'
            av[name]['modelFiles'] = listModelDataFiles(
                jobID, 'ptrc_T', paths.ModelFolder_pref, annual)
            av[name]['dataFile'] = ''

            av[name]['modelcoords'] = medusaCoords
            av[name]['datacoords'] = maredatCoords

            av[name]['modeldetails'] = {
                'name': name,
                'vars': ['CHN', 'CHD'],
                'convert': ukp.sums,
                'units': 'mg C/m^3'
            }
            av[name]['datadetails'] = {'name': '', 'units': ''}

            av[name]['layers'] = [
                'Surface',
                '100m',
                '200m',
            ]
            av[name]['regions'] = regionList
            av[name]['metrics'] = metricList

            av[name]['datasource'] = ''
            av[name]['model'] = 'MEDUSA'

            av[name]['modelgrid'] = 'eORCA1'
            av[name]['gridFile'] = paths.orcaGridfn
            av[name]['Dimensions'] = 3

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

        if 'Iron' in analysisKeys:
            name = 'Iron'
            av[name]['modelFiles'] = listModelDataFiles(
                jobID, 'ptrc_T', paths.ModelFolder_pref, annual)

            av[name][
                'dataFile'] = paths.icFold + "/UKESM_fields_1860_eORCA1_small.nc"
            av[name]['modelcoords'] = medusaCoords
            av[name]['datacoords'] = icCoords
            av[name]['modeldetails'] = {
                'name': name,
                'vars': [
                    'FER',
                ],
                'convert': ukp.mul1000,
                'units': 'umolFe/m3'
            }
            av[name]['datadetails'] = {
                'name': name,
                'vars': [
                    'FER',
                ],
                'convert': ukp.mul1000,
                'units': 'umolFe/m3'
            }
            av[name]['layers'] = layerList
            av[name]['regions'] = regionList
            av[name]['metrics'] = metricList
            av[name]['datasource'] = 'InititialCondition'
            av[name]['model'] = 'MEDUSA'
            av[name]['modelgrid'] = 'eORCA1'
            av[name]['gridFile'] = orcaGridfn
            av[name]['Dimensions'] = 3

        if 'N' in analysisKeys:
            name = 'Nitrate'
            av[name]['modelFiles'] = listModelDataFiles(
                jobID, 'ptrc_T', paths.ModelFolder_pref, annual)
            if annual:
                av[name]['dataFile'] = WOAFolder + '/woa13_all_n00_01.nc'
            else:
                av[name]['dataFile'] = WOAFolder + '/nitrate_monthly_1deg.nc'

            av[name]['modelcoords'] = medusaCoords
            av[name]['datacoords'] = woaCoords

            av[name]['modeldetails'] = {
                'name': name,
                'vars': [
                    'DIN',
                ],
                'convert': ukp.NoChange,
                'units': 'mmol N/m^3'
            }
            av[name]['datadetails'] = {
                'name': name,
                'vars': [
                    'n_an',
                ],
                'convert': ukp.NoChange,
                'units': 'mmol N/m^3'
            }

            av[name]['layers'] = layerList
            av[name]['regions'] = regionList

            #av[name]['layers'] 		= ['Surface','300m',]#'1000m',]#'Surface - 300m',]'100m',
            #av[name]['regions'] 		= regionList#['Global',]#'NorthAtlanticOcean','SouthAtlanticOcean',]#'NorthAtlantic']
            av[name]['metrics'] = metricList  #['mean','median', ]

            av[name]['datasource'] = 'WOA'
            av[name]['model'] = 'MEDUSA'

            av[name]['modelgrid'] = 'eORCA1'
            av[name]['gridFile'] = orcaGridfn
            av[name]['Dimensions'] = 3

        if 'Si' in analysisKeys:
            name = 'Silicate'
            av[name]['modelFiles'] = listModelDataFiles(
                jobID, 'ptrc_T', paths.ModelFolder_pref, annual)
            if annual:
                av[name]['dataFile'] = WOAFolder + 'woa13_all_i00_01.nc'
            else:
                av[name]['dataFile'] = WOAFolder + 'wsilicate_monthly_1deg.nc'
            av[name]['modelcoords'] = medusaCoords
            av[name]['datacoords'] = woaCoords

            av[name]['modeldetails'] = {
                'name': name,
                'vars': [
                    'SIL',
                ],
                'convert': ukp.NoChange,
                'units': 'mmol Si/m^3'
            }
            av[name]['datadetails'] = {
                'name': name,
                'vars': [
                    'i_an',
                ],
                'convert': ukp.NoChange,
                'units': 'mmol Si/m^3'
            }

            av[name]['layers'] = layerList
            av[name]['regions'] = regionList
            av[name]['metrics'] = metricList

            av[name]['datasource'] = 'WOA'
            av[name]['model'] = 'MEDUSA'

            av[name]['modelgrid'] = 'eORCA1'
            av[name]['gridFile'] = orcaGridfn
            av[name]['Dimensions'] = 3

        if 'O2' in analysisKeys:
            name = 'Oxygen'
            if annual:
                av[name]['modelFiles'] = listModelDataFiles(
                    jobID, 'ptrc_T', paths.ModelFolder_pref, annual)
                av[name]['dataFile'] = WOAFolder + 'woa13_all_o00_01.nc'

            av[name]['modelcoords'] = medusaCoords
            av[name]['datacoords'] = woaCoords

            av[name]['modeldetails'] = {
                'name': name,
                'vars': [
                    'OXY',
                ],
                'convert': ukp.NoChange,
                'units': 'mmol O2/m^3'
            }
            av[name]['datadetails'] = {
                'name': name,
                'vars': [
                    'o_an',
                ],
                'convert': ukp.oxconvert,
                'units': 'mmol O2/m^3'
            }

            av[name]['layers'] = layerList
            av[name]['regions'] = regionList
            av[name]['metrics'] = metricList

            av[name]['datasource'] = 'WOA'
            av[name]['model'] = 'MEDUSA'

            av[name]['modelgrid'] = 'eORCA1'
            av[name]['gridFile'] = orcaGridfn
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

        if 'DIC' in analysisKeys:

            def convertkgToM3(nc, keys):
                return nc.variables[keys[0]][:] * 1.027

            name = 'DIC'

            av[name]['modelFiles'] = listModelDataFiles(
                jobID, 'ptrcG_T', paths.ModelFolder_pref, annual)
            av[name][
                'dataFile'] = paths.GLODAPv2Dir + 'GLODAPv2.tco2.historic.nc'

            av[name]['modelcoords'] = medusaCoords
            av[name]['datacoords'] = glodapv2Coords

            av[name]['modeldetails'] = {
                'name': 'DIC',
                'vars': [
                    'DIC',
                ],
                'convert': ukp.NoChange,
                'units': 'mmol C/m^3'
            }
            av[name]['datadetails'] = {
                'name': 'DIC',
                'vars': [
                    'tco2',
                ],
                'convert': ukp.convertkgToM3,
                'units': 'mmol C/m^3'
            }

            av[name]['layers'] = layerList
            av[name]['regions'] = regionList
            av[name]['metrics'] = metricList

            av[name]['datasource'] = 'GLODAP'
            av[name]['model'] = 'MEDUSA'

            av[name]['modelgrid'] = 'eORCA1'
            av[name]['gridFile'] = orcaGridfn
            av[name]['Dimensions'] = 3

        if 'Alk' in analysisKeys:

            def convertmeqm3TOumolkg(nc, keys):
                return nc.variables[keys[0]][:] * 1.027

            name = 'Alkalinity'
            if annual:
                av[name]['modelFiles'] = listModelDataFiles(
                    jobID, 'ptrc_T', paths.ModelFolder_pref, annual)
                av[name]['dataFile'] = paths.GlodapDir + 'Alk.nc'
            else:
                print("Alkalinity data not available for monthly Analysis")
                assert 0

            av[name]['modelcoords'] = medusaCoords
            av[name]['datacoords'] = glodapCoords

            av[name]['modeldetails'] = {
                'name': name,
                'vars': [
                    'ALK',
                ],
                'convert': ukp.NoChange,
                'units': 'meq/m^3',
            }
            av[name]['datadetails'] = {
                'name': name,
                'vars': [
                    'Alk',
                ],
                'convert': convertmeqm3TOumolkg,
                'units': 'meq/m^3',
            }

            #	av[name]['layers'] 		=  ['Surface','100m','300m','1000m',]
            #	av[name]['regions'] 		= regionList
            av[name]['layers'] = layerList
            av[name]['regions'] = regionList
            av[name]['metrics'] = metricList

            av[name]['datasource'] = 'GLODAP'
            av[name]['model'] = 'MEDUSA'

            av[name]['modelgrid'] = 'eORCA1'
            av[name]['gridFile'] = orcaGridfn
            av[name]['Dimensions'] = 3

        if 'Temperature' in analysisKeys:
            name = 'Temperature'
            av[name]['modelFiles'] = listModelDataFiles(
                jobID, 'grid_T', paths.ModelFolder_pref, annual)
            if annual:
                #av[name]['modelFiles']  	= sorted(glob(paths.ModelFolder_pref+jobID+"/"+jobID+"o_1y_*_grid_T.nc"))
                av[name]['dataFile'] = WOAFolder + 'woa13_decav_t00_01v2.nc'
            else:
                #av[name]['modelFiles']  	= sorted(glob(paths.ModelFolder_pref+jobID+"/"+jobID+"o_1m_*_grid_T.nc"))
                av[name][
                    'dataFile'] = WOAFolder + 'temperature_monthly_1deg.nc'
            av[name]['modelcoords'] = medusaCoords
            av[name]['datacoords'] = woaCoords

            av[name]['modeldetails'] = {
                'name': name,
                'vars': [
                    'votemper',
                ],
                'convert': ukp.NoChange,
                'units': 'degrees C'
            }
            av[name]['datadetails'] = {
                'name': name,
                'vars': [
                    't_an',
                ],
                'convert': ukp.NoChange,
                'units': 'degrees C'
            }

            av[name]['layers'] = layerList
            av[name]['regions'] = vmtregionList
            av[name]['metrics'] = metricList

            av[name]['datasource'] = 'WOA'
            av[name]['model'] = 'NEMO'

            av[name]['modelgrid'] = 'eORCA1'
            av[name]['gridFile'] = orcaGridfn
            av[name]['Dimensions'] = 3

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

        if 'GlobalMeanTemperature' in analysisKeys:
            name = 'GlobalMeanTemperature'
            av[name]['modelFiles'] = listModelDataFiles(
                jobID, 'grid_T', paths.ModelFolder_pref, annual)
            av[name]['dataFile'] = ''

            av[name]['modelcoords'] = medusaCoords
            av[name]['datacoords'] = woaCoords

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

            def sumMeanLandMask(nc, keys):
                temperature = np.ma.masked_where(
                    tmask == 0, nc.variables[keys[0]][:].squeeze())
                return (temperature * pvol).sum() / (pvol.sum())

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

        if 'GlobalMeanSalinity' in analysisKeys:
            name = 'GlobalMeanSalinity'
            av[name]['modelFiles'] = listModelDataFiles(
                jobID, 'grid_T', paths.ModelFolder_pref, annual)
            av[name]['dataFile'] = ''

            av[name]['modelcoords'] = medusaCoords
            av[name]['datacoords'] = woaCoords

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

            def sumMeanLandMask(nc, keys):
                temperature = np.ma.masked_where(
                    tmask == 0, nc.variables[keys[0]][:].squeeze())
                return (temperature * pvol).sum() / (pvol.sum())

            av[name]['modeldetails'] = {
                'name': name,
                'vars': [
                    'vosaline',
                ],
                'convert': sumMeanLandMask,
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

        if 'quickSST' in analysisKeys:
            name = 'quickSST'
            av[name]['modelFiles'] = listModelDataFiles(
                jobID, 'grid_T', paths.ModelFolder_pref, annual)

            nc = Dataset(paths.orcaGridfn, 'r')
            ssttmask = nc.variables['tmask'][0]
            area = nc.variables['e2t'][:] * nc.variables['e1t'][:]
            area = np.ma.masked_where(ssttmask == 0, area)
            nc.close()

            def meanLandMask(nc, keys):
                #### works like no change, but applies a mask.
                #print "meanLandMask:",ssttmask.shape,nc.variables[keys[0]][0,0].shape
                temperature = np.ma.masked_where(
                    ssttmask == 0, nc.variables[keys[0]][0, 0].squeeze())
                print("meanLandMask:", nc.variables['time_counter'][:],
                      temperature.mean(),
                      (temperature * area).sum() / (area.sum()))
                return (temperature * area).sum() / (area.sum())

            if annual:
                #av[name]['modelFiles']  	= sorted(glob(paths.ModelFolder_pref+jobID+"/"+jobID+"o_1y_*_grid_T.nc"))
                av[name]['dataFile'] = ''  #WOAFolder+'woa13_decav_t00_01v2.nc'
            else:
                #av[name]['modelFiles']  	= sorted(glob(paths.ModelFolder_pref+jobID+"/"+jobID+"o_1m_*_grid_T.nc"))
                av[name][
                    'dataFile'] = ''  #WOAFolder+'temperature_monthly_1deg.nc'

            av[name]['modelcoords'] = medusaCoords
            av[name]['datacoords'] = woaCoords

            av[name]['modeldetails'] = {
                'name': name,
                'vars': [
                    'votemper',
                ],
                'convert': meanLandMask,
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

        if 'Salinity' in analysisKeys:
            name = 'Salinity'
            av[name]['modelFiles'] = listModelDataFiles(
                jobID, 'grid_T', paths.ModelFolder_pref, annual)
            if annual:
                #av[name]['modelFiles']  	= sorted(glob(paths.ModelFolder_pref+jobID+"/"+jobID+"o_1y_*_grid_T.nc"))
                av[name]['dataFile'] = WOAFolder + 'woa13_decav_s00_01v2.nc'
            else:
                #av[name]['modelFiles']  	= sorted(glob(paths.ModelFolder_pref+jobID+"/"+jobID+"o_1m_*_grid_T.nc"))
                av[name]['dataFile'] = WOAFolder + 'salinity_monthly_1deg.nc'

            av[name]['modelcoords'] = medusaCoords
            av[name]['datacoords'] = woaCoords

            av[name]['modeldetails'] = {
                'name': name,
                'vars': [
                    'vosaline',
                ],
                'convert': ukp.NoChange,
                'units': 'PSU'
            }
            av[name]['datadetails'] = {
                'name': name,
                'vars': [
                    's_an',
                ],
                'convert': ukp.NoChange,
                'units': 'PSU'
            }

            av[name]['layers'] = layerList
            av[name]['regions'] = vmtregionList
            av[name]['metrics'] = metricList

            av[name]['datasource'] = 'WOA'
            av[name]['model'] = 'NEMO'

            av[name]['modelgrid'] = 'eORCA1'
            av[name]['gridFile'] = orcaGridfn
            av[name]['Dimensions'] = 3

        if 'MLD' in analysisKeys:

            def mldapplymask(nc, keys):
                mld = np.ma.array(nc.variables[keys[0]][:])
                return np.ma.masked_where((nc.variables[keys[1]][:] == 0.) +
                                          mld.mask + (mld == 1.E9), mld)

            name = 'MLD'
            av[name]['modelFiles'] = listModelDataFiles(
                jobID, 'grid_T', paths.ModelFolder_pref, annual)
            av[name][
                'dataFile'] = paths.MLDFolder + "mld_DT02_c1m_reg2.0-annual.nc"  #mld_DT02_c1m_reg2.0.nc"

            av[name]['modelcoords'] = medusaCoords
            av[name]['datacoords'] = mldCoords

            av[name]['modeldetails'] = {
                'name': 'mld',
                'vars': [
                    'somxl010',
                ],
                'convert': ukp.NoChange,
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

            av[name]['layers'] = [
                'layerless',
            ]  #'Surface - 1000m','Surface - 300m',]#'depthint']

            av[name]['regions'] = vmtregionList
            av[name]['metrics'] = metricList

            av[name]['datasource'] = 'IFREMER'
            av[name]['model'] = 'NEMO'

            av[name]['modelgrid'] = 'eORCA1'
            av[name]['gridFile'] = paths.orcaGridfn
            av[name]['Dimensions'] = 2

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

        if 'AMOC_26N' in analysisKeys or 'AMOC_32S' in analysisKeys or 'ADRC_26N' in analysisKeys:
            # Note that this will only work with the eORCAgrid.
            latslice26N = slice(227, 228)
            latslice32S = slice(137, 138)
            e3v, e1v, tmask, alttmask = {}, {}, {}, {}
            for name in ['AMOC_26N', 'AMOC_32S', 'ADRC_26N']:
                if name not in analysisKeys: continue

                ####
                if name in ['AMOC_26N', 'ADRC_26N']: latslice = latslice26N
                if name == 'AMOC_32S': latslice = latslice32S

                # Load grid data
                nc = Dataset(paths.orcaGridfn, 'r')
                e3v[name] = nc.variables[
                    'e3v'][:, latslice, :]  # z level height 3D
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
                zv = np.ma.array(
                    nc.variables['vomecrty'][..., latslice32S, :])  # m/s
                atlmoc = np.array(np.zeros_like(zv[0, :, :, 0]))
                e2vshape = e3v[name].shape
                for la in range(e2vshape[1]):  #ji, y
                    for lo in range(e2vshape[2]):  #jj , x,
                        if int(alttmask[name][la, lo]) == 0: continue
                        for z in range(e2vshape[0]):  # jk
                            if int(tmask[name][z, la, lo]) == 0: continue
                            if np.ma.is_masked(zv[0, z, la, lo]): continue
                            atlmoc[z, la] = atlmoc[z, la] - e1v[name][
                                la, lo] * e3v[name][z, la, lo] * zv[0, z, la,
                                                                    lo] / 1.E06

    ####
    # Cumulative sum from the bottom up.
                for z in range(73, 1, -1):
                    atlmoc[z, :] = atlmoc[z + 1, :] + atlmoc[z, :]
                return np.ma.max(atlmoc)

            def amoc26N_array(nc, keys, amocname='AMOC_26N'):
                zv = np.ma.array(
                    nc.variables['vomecrty'][..., latslice26N, :])  # m/s
                atlmoc = np.array(np.zeros_like(zv[0, :, :, 0]))
                e2vshape = e3v[amocname].shape
                for la in range(e2vshape[1]):  #ji, y
                    for lo in range(e2vshape[2]):  #jj , x,
                        if int(alttmask[amocname][la, lo]) == 0: continue
                        for z in range(e2vshape[0]):  # jk
                            if int(tmask[amocname][z, la, lo]) == 0: continue
                            if np.ma.is_masked(zv[0, z, la, lo]): continue
                            atlmoc[z, la] = atlmoc[
                                z, la] - e1v[amocname][la, lo] * e3v[amocname][
                                    z, la, lo] * zv[0, z, la, lo] / 1.E06

    ####
    # Cumulative sum from the bottom up.
                for z in range(73, 1, -1):
                    atlmoc[z, :] = atlmoc[z + 1, :] + atlmoc[z, :]
                #return np.ma.max(atlmoc)
                return atlmoc

            def calc_amoc26N(nc, keys):
                return np.ma.max(amoc26N_array(nc, keys, amocname='AMOC_26N'))

            def calc_min_amoc26N(nc, keys):
                return np.ma.min(amoc26N_array(nc, keys, amocname='ADRC_26N'))

            for name in ['AMOC_26N', 'AMOC_32S', 'ADRC_26N']:
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
                if name == 'ADRC_26N':
                    av[name]['modeldetails'] = {
                        'name': name,
                        'vars': [
                            'vomecrty',
                        ],
                        'convert': calc_min_amoc26N,
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

        if 'TotalDust' in analysisKeys:
            name = 'TotalDust'
            av[name]['modelFiles'] = listModelDataFiles(
                jobID, 'diad_T', paths.ModelFolder_pref, annual)[:]
            av[name]['dataFile'] = paths.Dustdir + 'mahowald.orca100_annual.nc'

            av[name]['modelcoords'] = medusaCoords
            av[name]['datacoords'] = medusaCoords

            nc = Dataset(paths.orcaGridfn, 'r')
            masked_area = nc.variables['e2t'][:] * nc.variables[
                'e1t'][:] * nc.variables['tmask'][0]
            nc.close()

            def datadustsum(nc, keys):
                #factors are:
                # 0.035: iron as a fraction of total dust
                # 1e6: convert from kmol -> mmol
                # 1e-12: convert from mol to Gmol
                # 0.00532: solubility factor for iron
                # 55.845: atmoic mass of iron (g>mol conversion)
                # (24.*60.*60.*365.25): per second to per year

                dust = nc.variables[keys[0]][:]
                dust[:, :, 234:296, 295:348] = 0.
                dust[:, :, 234:248, 285:295] = 0.
                dust[:, :, 228:256, 290:304] = 0.
                return (masked_area *
                        dust).sum() * 0.035 * 1.e6 * 1.e-12 * 0.00532 * (
                            24. * 60. * 60. * 365.25) / 55.845

            def modeldustsum(nc, keys):
                dust = nc.variables[keys[0]][:]
                dust[:, 234:296, 295:348] = 0.
                dust[:, 234:248, 285:295] = 0.
                dust[:, 228:256, 290:304] = 0.
                return (masked_area * dust).sum() * 1.E-12 * 365.25

            av[name]['modeldetails'] = {
                'name': name,
                'vars': [
                    'AEOLIAN',
                ],
                'convert': modeldustsum,
                'units': 'Gmol Fe/yr'
            }
            av[name]['datadetails'] = {
                'name': name,
                'vars': [
                    'dust_ann',
                ],
                'convert': datadustsum,
                'units': 'Gmol Fe/yr'
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

            av[name]['datasource'] = 'Mahowald'
            av[name]['model'] = 'MEDUSA'

            av[name]['modelgrid'] = 'eORCA1'
            av[name]['gridFile'] = paths.orcaGridfn
            av[name]['Dimensions'] = 1

        if 'Dust' in analysisKeys:
            name = 'Dust'
            av[name]['modelFiles'] = listModelDataFiles(
                jobID, 'diad_T', paths.ModelFolder_pref, annual)[:]
            av[name]['dataFile'] = paths.Dustdir + 'mahowald.orca100_annual.nc'

            av[name]['modelcoords'] = medusaCoords
            av[name]['datacoords'] = medusaCoords
            av[name]['modeldetails'] = {
                'name': name,
                'vars': [
                    'AEOLIAN',
                ],
                'convert': ukp.NoChange,
                'units': 'mmol Fe/m2/d'
            }

            def mahodatadust(nc, keys):
                #factors are:
                # 0.035: iron as a fraction of total dust
                # 1e6: convert from kmol -> mmol
                # 0.00532: solubility factor or iron
                # 55.845: atmoic mass of iron (g>mol conversion)
                # (24.*60.*60.): per second to per day
                dust = nc.variables[keys[0]][:]
                dust[0, 0, 194:256, 295:348] = 0.
                dust[0, 0, 194:208, 285:295] = 0.
                dust[0, 0, 188:216, 290:304] = 0.
                return dust * 0.035 * 1.e6 * 0.00532 * (24. * 60. *
                                                        60.) / 55.845

            av[name]['datadetails'] = {
                'name': name,
                'vars': [
                    'dust_ann',
                ],
                'convert': mahodatadust,
                'units': 'mmol Fe/m2/d'
            }

            av[name]['layers'] = [
                'layerless',
            ]
            av[name]['regions'] = regionList
            av[name]['metrics'] = metricList

            av[name]['datasource'] = 'Mahowald'
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
                noNewFiles=True,
            )
            #dataD[(jobID,name )] = tsa.dataD
            modeldataD[(jobID, name)] = tsa.modeldataD

    #####
    # Data now loaded, making plots next:
    for k in list(modeldataD.keys()):
        print("Model Data D:", k)

    ####
    for name in [
            'Temperature', 'Salinity', 'MLD', 'FreshwaterFlux',
            'AirSeaFluxCO2', 'AirSeaFlux', 'Chlorophyll', 'Nitrate',
            'Alkalinity', 'pH'
    ]:
        if name not in list(av.keys()): continue
        for region in vmtregionList:
            for layer in [
                    'Surface',
                    '500m',
                    '1000m',
                    'layerless',
            ]:
                timesD = {}
                arrD = {}

                for jobID in jobs:
                    try:
                        mdata = modeldataD[(jobID, name)][(region, layer,
                                                           'mean')]
                    except:
                        continue
                    title = titleify([region, layer, 'Mean', name])

                    #timesD[jobID] 	= sorted(mdata.keys())
                    #arrD[jobID]	= [mdata[t] for t in timesD[jobID]]
                    times, datas = shifttimes(mdata, jobID, year0=year0)
                    timesD[jobID] = times  #mdata.keys())
                    arrD[jobID] = datas  #t] for t in timesD[jobID]]

                    if jobID == 'u-aj588':
                        arrD[jobID] = np.ma.masked_where(
                            arrD[jobID] == 0., arrD[jobID])

    #times,datas = shifttimes(mdata, jobID,year0=year0)
                timesD, arrD = build_ensemble(timesD, arrD, ensembles)

                if len(list(arrD.keys())) == 0: continue
                units = av[name]['modeldetails']['units']

                for ts in [
                        'Together',
                ]:
                    if name == 'FreshwaterFlux': ls = 'movingav30years'
                    else: ls = 'DataOnly'
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


#		assert 0

####
# Standard surface:
    for name in list(av.keys()):
        timesD = {}
        arrD = {}
        units = av[name]['modeldetails']['units']
        title = ''
        for jobID in jobs:
            if name in [
                    'Iron', 'Nitrate', 'Silicate', 'Oxygen', 'Temperature',
                    'Salinity', 'O2', 'Alkalinity', 'DIC', 'CHD', 'CHN',
                    'DiaFrac', 'CHL', 'Chlorophyll', 'pH', 'ZonalCurrent',
                    'MeridionalCurrent', 'VerticalCurrent'
            ]:
                mdata = modeldataD[(jobID, name)][('Global', 'Surface',
                                                   'mean')]
                title = ' '.join(
                    ['Global', 'Surface', 'Mean',
                     getLongName(name)])
            elif name in [
                    'OMZThickness',
                    'OMZMeanDepth',
                    'DMS',
                    'DMS_ARAN',
                    'Dust',
                    'MaxMonthlyMLD',
                    'MinMonthlyMLD',
                    'MLD',
            ]:
                try:
                    mdata = modeldataD[(jobID, name)][('Global', 'layerless',
                                                       'mean')]
                    title = ' '.join(['Global', getLongName(name)])
                except:
                    continue
            elif name in [
                    'DTC',
            ]:
                mdata = modeldataD[(jobID, name)][('Global', '3000m', 'mean')]
                title = ' '.join(
                    ['Global', '3000m', 'Mean',
                     getLongName(name)])
            elif name in [
                    'AirSeaFlux',
                    'AirSeaFluxCO2',
                    'Alk',
                    'Alkalinity',
            ]:
                try:
                    mdata = modeldataD[(jobID, name)][('ignoreInlandSeas',
                                                       'Surface', 'mean')]
                    title = ' '.join([
                        'ignoreInlandSeas', 'Surface', 'Mean',
                        getLongName(name)
                    ])
                except:
                    continue
            elif name in [
                    'NoCaspianAirSeaFluxCO2',
            ]:
                try:
                    mdata = modeldataD[(jobID, name)][('ignoreCaspian',
                                                       'layerless', 'sum')]
                    title = ' '.join(
                        ['Total', getLongName(name), '(No Caspian)'])
                except:
                    continue
            elif name in [
                    'VolumeMeanTemperature',
                    'VolumeMeanSalinity',
                    'VolumeMeanOxygen',
            ]:
                try:
                    mdata = modeldataD[(jobID, name)][('Global', 'layerless',
                                                       'wcvweighted')]
                    title = ' '.join([
                        'Global',
                        getLongName(name),
                    ])
                except:
                    continue

            elif name in [
                    'sowaflup', 'sohefldo', 'sofmflup', 'sosfldow', 'sossheig',
                    'soicecov', 'FreshwaterFlux', 'HeatFlux'
            ]:

                try:
                    mdata = modeldataD[(jobID, name)][('Global', 'layerless',
                                                       'mean')]
                except:
                    continue
                title = ' '.join(['Global mean', getLongName(name)])

            #####
            # Special hack for these guys.
            #nasregionList	= ['NordicSea', 'LabradorSea', 'NorwegianSea'	]
            #mdata = modeldataD[(jobID,name )][('regionless', 'layerless', 'mean')]
            #title = getLongName(name)
            else:
                try:
                    mdata = modeldataD[(jobID,
                                        name)][('regionless', 'layerless',
                                                'metricless')]
                except:
                    continue

                title = getLongName(name)

            times, datas = shifttimes(mdata, jobID, year0=year0)
            timesD[jobID] = times
            arrD[jobID] = datas
        timesD, arrD = build_ensemble(timesD, arrD, ensembles)
        #####
        # To account for changing units.
        if name in ['TotalAirSeaFluxCO2', 'TotalAirSeaFlux']:
            for j in list(arrD.keys()):
                if j in [
                        'u-ad980', 'u-af123', 'u-af725', 'u-ae742', 'u-af139',
                        'u-af578', 'u-af728'
                ]:
                    arrD[j] = np.ma.array(arrD[j]) * 5.09369e-7

        if name in [
                'DMS',
        ]:
            for j in list(arrD.keys()):
                if j in [
                        'u-ag914',
                ]:
                    for i, (t, d) in enumerate(zip(timesD[jobID],
                                                   arrD[jobID])):
                        if float(t) < 1600.: continue
                        arrD[j][i] = d / 1000000.

        for ts in [
                'Together',
        ]:  #'Separate']:
            for ls in [
                    'DataOnly',
            ]:  #'','Both',]:
                if ls == '' and name not in level3: continue

                tsp.multitimeseries(
                    timesD,  # model times (in floats)
                    arrD,  # model time series
                    data=-999,  # in situ data distribution
                    title=title,
                    filename=ukp.folder(imageFolder) + name + '_' + ts + '_' +
                    ls + '.png',
                    units=units,
                    plotStyle=ts,
                    lineStyle=ls,
                    colours=colours,
                    thicknesses=lineThicknesses,
                    linestyles=linestyles,
                )

                if name in [
                        'NoCaspianAirSeaFluxCO2',
                        'AMOC_26N',
                        'ADRC_26N',
                        'AMOC_32S',
                        'TotalAirSeaFluxCO2',
                ]:
                    targetdict = {
                        'NoCaspianAirSeaFluxCO2': 0.,
                        'AMOC_26N': 17.,
                        'ADRC_26N': -999,
                        'AMOC_32S': -999,
                        'TotalAirSeaFluxCO2': 0.,
                    }
                    for ls in ['movingav30years', 'movingav100years']:
                        tsp.multitimeseries(
                            timesD,  # model times (in floats)
                            arrD,  # model time series
                            data=targetdict[name],  # in situ gata distribution
                            dataname='Target',
                            title=title,
                            filename=ukp.folder(imageFolder) + name + '_' +
                            ts + '_' + ls + '.png',
                            units=units,
                            plotStyle=ts,
                            lineStyle=ls,
                            colours=colours,
                            thicknesses=lineThicknesses,
                            linestyles=linestyles,
                        )

    ####
    # Oxygen at Depth:
    regionList = [
        'Global',
        'ignoreInlandSeas',
        'SouthernOcean',
        'ArcticOcean',
        'AtlanticSOcean',
        'Equator10',
        'Remainder',
        'NorthernSubpolarAtlantic',
        'NorthernSubpolarPacific',
    ]
    for name in [
            'Oxygen',
            'O2',
    ]:
        if name not in list(av.keys()): continue
        for region in regionList:
            for layer in ['Surface', '500m', '1000m']:
                timesD = {}
                arrD = {}

                for jobID in jobs:
                    try:
                        mdata = modeldataD[(jobID, name)][(region, layer,
                                                           'mean')]
                    except:
                        continue
                    title = ' '.join(
                        [region, layer, 'Mean',
                         getLongName(name)])
                    times, datas = shifttimes(mdata, jobID, year0=year0)
                    timesD[jobID] = times
                    arrD[jobID] = datas

                timesD, arrD = build_ensemble(timesD, arrD, ensembles)
                units = av[name]['modeldetails']['units']

                for ts in [
                        'Together',
                ]:  #'Separate']:
                    for ls in [
                            'DataOnly',
                    ]:  #'','Both',]:
                        tsp.multitimeseries(
                            timesD,  # model times (in floats)
                            arrD,  # model time series
                            data=-999,  # in situ data distribution
                            title=title,
                            filename=ukp.folder(imageFolder + '/Oxygen') +
                            '_'.join([name, region, layer, ts, ls + '.png']),
                            units=units,
                            plotStyle=ts,
                            lineStyle=ls,
                            colours=colours,
                            thicknesses=lineThicknesses,
                            linestyles=linestyles,
                        )
    wcvwRegions = vmtregionList[:]
    wcvwRegions.extend(OMZRegions)
    for name in [
            'VolumeMeanTemperature',
            'VolumeMeanSalinity',
            'VolumeMeanOxygen',
    ]:
        if name not in list(av.keys()): continue
        for region in wcvwRegions:
            for layer in [
                    'layerless',
            ]:
                timesD = {}
                arrD = {}

                for jobID in jobs:
                    try:
                        mdata = modeldataD[(jobID, name)][(region, layer,
                                                           'wcvweighted')]
                    except:
                        continue
                    title = ' '.join([region, layer, getLongName(name)])
                    times, datas = shifttimes(mdata, jobID, year0=year0)
                    timesD[jobID] = times
                    arrD[jobID] = datas
                if len(list(arrD.keys())) == 0: continue
                timesD, arrD = build_ensemble(timesD, arrD, ensembles)
                units = av[name]['modeldetails']['units']
                for ts in [
                        'Together',
                ]:  #'Separate']:
                    for ls in [
                            'DataOnly',
                    ]:  #'','Both',]:
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

    for name in [
            'DiaFrac',
            'CHD',
            'CHN',
            'CHL',
            'N',
            'Si',
            'Iron',
            'Alk',
            'DIC',
            'Chlorophyll',
            'pH',
            'DMS',
            'Nitrate',
            'Salinity',
            'Silicate',
            'MaxMonthlyMLD',
            'MinMonthlyMLD',
            'Dust',
    ]:
        if name not in list(av.keys()): continue

        for region in regionList:
            for layer in [
                    'Surface',
                    '100m',
                    '200m',
                    'layerless',
            ]:

                timesD = {}
                arrD = {}

                for jobID in jobs:
                    try:
                        mdata = modeldataD[(jobID, name)][(region, layer,
                                                           'mean')]
                    except:
                        continue
                    title = ' '.join(
                        [region, layer, 'Mean',
                         getLongName(name)])
                    times, datas = shifttimes(mdata, jobID, year0=year0)
                    timesD[jobID] = times
                    arrD[jobID] = datas

                timesD, arrD = build_ensemble(timesD, arrD, ensembles)
                units = av[name]['modeldetails']['units']

                for ts in [
                        'Together',
                ]:  #'Separate']:
                    for ls in [
                            'DataOnly',
                    ]:  #'','Both',]:
                        tsp.multitimeseries(
                            timesD,  # model times (in floats)
                            arrD,  # model time series
                            data=-999,  # in situ data distribution
                            title=title,
                            filename=ukp.folder(imageFolder + '/BGC') +
                            '_'.join([name, region, layer, ts, ls + '.png']),
                            units=units,
                            plotStyle=ts,
                            lineStyle=ls,
                            colours=colours,
                            thicknesses=lineThicknesses,
                            linestyles=linestyles,
                        )

    for name in [
            'HeatFlux',
            'AirSeaFlux',
    ]:
        if name not in list(av.keys()): continue
        for region in vmtregionList:
            for layer in [
                    'layerless',
            ]:
                timesD = {}
                arrD = {}
                for jobID in jobs:
                    try:
                        mdata = modeldataD[(jobID, name)][(region, layer,
                                                           'mean')]
                    except:
                        continue
                    title = ' '.join(
                        [region, layer, 'Mean',
                         getLongName(name)])
                    times, datas = shifttimes(mdata, jobID, year0=year0)
                    timesD[jobID] = times
                    arrD[jobID] = datas
                timesD, arrD = build_ensemble(timesD, arrD, ensembles)
                units = av[name]['modeldetails']['units']
                for ts in [
                        'Together',
                ]:  #'Separate']:
                    for ls in [
                            'DataOnly',
                    ]:  #'','Both',]:
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

    for name in [
            'scalarHeatContent',
    ]:
        if name not in list(av.keys()): continue
        region = 'regionless'
        layer = 'layerless'
        metric = 'metricless'
        timesD = {}
        arrD = {}
        for jobID in jobs:
            try:
                mdata = modeldataD[(jobID, name)][(region, layer, metric)]
            except:
                continue
            title = ' '.join(['Year to year change in ', getLongName(name)])
            times, datas = shifttimes(mdata, jobID, year0=year0)
            if len(times) < 3: continue
            dtimes, ddatas = [], []

            for i, t in enumerate(times[:-1]):
                tdiff = times[i + 1] - t  # time difference in years
                dtimes.append((times[i + 1] + t) / 2.)  # midpoint
                ddatas.append((datas[i + 1] - datas[i]) / tdiff)
            timesD[jobID] = dtimes
            arrD[jobID] = ddatas
        timesD, arrD = build_ensemble(timesD, arrD, ensembles)

        units = r'$\Delta$ ' + av[name]['modeldetails'][
            'units'] + ' y' + r'$^{-1}$'
        for ts in [
                'Together',
        ]:  #'Separate']:
            for ls in ['DataOnly', 'movingav30years']:  #'','Both',]:
                tsp.multitimeseries(
                    timesD,  # model times (in floats)
                    arrD,  # model time series
                    data=0.,  # in situ data distribution
                    dataname='Zero',
                    title=title,
                    filename=ukp.folder(imageFolder) + '_'.join(
                        ['Change_in', name, region, layer, ts, ls + '.png']),
                    units=units,
                    plotStyle=ts,
                    lineStyle=ls,
                    colours=colours,
                    thicknesses=lineThicknesses,
                    linestyles=linestyles,
                )

    for name in [
            'MaxMonthlyMLD',
            'MinMonthlyMLD',
    ]:
        if name not in list(av.keys()): continue
        for region in vmtregionList:
            for layer in [
                    'layerless',
            ]:

                timesD = {}
                arrD = {}

                for jobID in jobs:
                    try:
                        mdata = modeldataD[(jobID, name)][(region, layer,
                                                           'mean')]
                    except:
                        continue
                    title = ' '.join(
                        [region, layer, 'Mean',
                         getLongName(name)])
                    times, datas = shifttimes(mdata, jobID, year0=year0)
                    timesD[jobID] = times
                    arrD[jobID] = datas

                units = av[name]['modeldetails']['units']
                timesD, arrD = build_ensemble(timesD, arrD, ensembles)
                for ts in [
                        'Together',
                ]:  #'Separate']:
                    for ls in [
                            'DataOnly',
                    ]:  #'','Both',]:
                        tsp.multitimeseries(
                            timesD,  # model times (in floats)
                            arrD,  # model time series
                            data=-999,  # in situ data distribution
                            title=title,
                            filename=ukp.folder(imageFolder + '/MLD') +
                            '_'.join([name, region, layer, ts, ls + '.png']),
                            units=units,
                            plotStyle=ts,
                            lineStyle=ls,
                            colours=colours,
                            thicknesses=lineThicknesses,
                            linestyles=linestyles,
                        )

    for name in [
            'DMS',
    ]:
        continue
        if name not in list(av.keys()): continue
        for region in regionList:
            for layer in [
                    'Surface',
                    '100m',
                    '200m',
            ]:

                timesD = {}
                arrD = {}

                for jobID in jobs:
                    try:
                        mdata = modeldataD[(jobID, name)][(region, layer,
                                                           'mean')]
                    except:
                        continue
                    title = ' '.join(
                        [region, layer, 'Mean',
                         getLongName(name)])

                    timesD[jobID] = sorted(mdata.keys())
                    arrD[jobID] = [mdata[t] for t in timesD[jobID]]
                units = av[name]['modeldetails']['units']
                for ts in [
                        'Together',
                ]:  #'Separate']:
                    for ls in [
                            'DataOnly',
                    ]:  #'','Both',]:
                        tsp.multitimeseries(
                            timesD,  # model times (in floats)
                            arrD,  # model time series
                            data=-999,  # in situ data distribution
                            title=title,
                            filename=ukp.folder(imageFolder + '/DMS') +
                            '_'.join([name, region, layer, ts, ls + '.png']),
                            units=units,
                            plotStyle=ts,
                            lineStyle=ls,
                            colours=colours,
                            thicknesses=lineThicknesses,
                            linestyles=linestyles,
                        )
    try:
        AllImages = glob(imageFolder, recursive=True)
    except:
        AllImages = []
        for root, dirnames, filenames in os.walk(imageFolder):
            for filename in fnmatch.filter(filenames, '*.png'):
                AllImages.append(os.path.join(root, filename))

    if ensembles != {}:
        jobs = list(ensembles.keys())
    comparehtml5Maker(
        jobIDs=jobs,
        reportdir=ukp.folder('CompareReports2/' + analysisname),
        files=AllImages,
        clean=False,
        doZip=False,
        jobDescriptions=jobDescriptions,
        jobColours=colours,
    )


def flatten(lats, lons, dataA, dataB):
    m = np.ma.array(dataA).mask
    m += np.ma.array(dataB).mask
    m += np.ma.masked_invalid(dataA / dataB).mask

    return  np.ma.masked_where(m, lats).compressed(),\
     np.ma.masked_where(m, lons).compressed(),\
     np.ma.masked_where(m, dataA).compressed(),\
     np.ma.masked_where(m, dataB).compressed()


def CompareTwoRuns(jobIDA,
                   jobIDB,
                   physics=True,
                   bio=True,
                   yearA='',
                   yearB='',
                   debug=True):
    #
    #spatial maps of specific fields

    filetype = []
    if physics:
        filetype.extend([
            'grid_T',
            'grid_U',
            'grid_V',
            'grid_W',
        ])
    if bio:
        filetype.extend([
            'diad_T',
            'ptrc_T',
        ])
    filesA = {}
    filesB = {}

    imageFolder = 'images/TimeseriesCompare/'
    imageFolder += jobIDA + '-' + yearA + '_and_' + jobIDB + '-' + yearB

    for ft in filetype:
        filesA[ft] = listModelDataFiles(jobIDA,
                                        ft,
                                        paths.ModelFolder_pref,
                                        True,
                                        year=yearA)[0]
        filesB[ft] = listModelDataFiles(jobIDB,
                                        ft,
                                        paths.ModelFolder_pref,
                                        True,
                                        year=yearB)[0]
        print("files A:", filesA[ft])
        print("files B:", filesB[ft])

        ncA = Dataset(filesA[ft], 'r')
        ncB = Dataset(filesB[ft], 'r')
        keys = ukp.intersection(list(ncA.variables.keys()),
                                list(ncB.variables.keys()))

        lats = ncA.variables['nav_lat'][:]
        lons = ncA.variables['nav_lon'][:]

        for key in keys:
            filename = ukp.folder(imageFolder + '/' +
                                  ft) + ft + '-' + key + '.png'
            if os.path.exists(filename): continue

            dataA = 0.
            dataB = 0.
            if key in alwaysInclude: continue
            if key in ['bounds_lon', 'bounds_lat']: continue
            if ncA.variables[key].ndim == 4:
                dataA = ncA.variables[key][0, 0, :, :]
                dataB = ncB.variables[key][0, 0, :, :]

            elif ncA.variables[key].ndim == 3:
                dataA = ncA.variables[key][0, :, :]
                dataB = ncB.variables[key][0, :, :]
            elif ncA.variables[key].ndim == 2:
                continue
                dataA = ncA.variables[key][:, :]
                dataB = ncB.variables[key][:, :]
            else:
                print("can't plot:", key, ncA.variables[key].ndim)
                continue
            try:
                title = ncA.variables[key].long_name
            except:
                title = getLongName(key)

            dmin = min([dataA.min(), dataB.min()])
            dmax = max([dataA.max(), dataB.max()])

            print(key, lats.shape, lons.shape, dataA.shape, dataB.shape)
            la, lo, data, datb = flatten(lats, lons, dataA, dataB)
            print(key, la.shape, lo.shape, data.shape, datb.shape)

            if 0 in [len(la), len(lo), len(data), len(datb)]: continue
            #		filename = ukp.folder(imageFolder+'/'+ft)+ft+'-'+key+'.png'
            #		if os.path.exists(filename):continue
            ukp.robinPlotQuad(
                lo,
                la,
                data,
                datb,
                filename,
                titles=[jobIDA + ' ' + yearA, jobIDB + ' ' + yearB],
                title=title,
                lon0=0.,
                marble=False,
                drawCbar=True,
                cbarlabel='',
                doLog=False,
                vmin=dmin,
                vmax=dmax,
                maptype='Cartopy',
                scatter=False)


def main():
# For now, there are no arguments needed. 
    if "--help" in argv or len(argv) == 1:
        print("Running with no arguments. Exiting.")
        if "--help" in argv:
            print("Read the documentation.")
        exit(0)

    yml_fn = argv[1]
    if not os.path.exists(yml_fn):
         print("ERROR:\tInput yml file", yml_fn, "does not exist. Exiting.")
         exit(0)

    if len(argv)>2: 
         config_user = argv[2]
    else: 
        config_user = ''

    """
    Sample input file:

    analysis_name: hotsauce
    keys:
        level1
    jobIDs : 
       u-ab123: 
        description: job number 1
        colour: green
        linestyle: -
        linewidth: 1.2
        year_diff: none
       u-ab122: 
        description: job 1 with extra hot sauce
        colour: red
        linestyle: -
        linewidth: 1.2
        year_diff: +75
       u-ab121: 
        description:  job 1 hold the hot sauce
        colour: blue
        linestyle: -
        linewidth: 1.2
        year_diff: -1850
    """

    with open(yml_fn, 'r') as file:
        yml = yaml.safe_load(file)
    #input_yml = _read_yaml(yml_fn)
    name = yml.get('analysis_name', None)
    keys = yml.get('keys', 'level1')
    jobIDs= yml.get('jobIDs', None).keys()
    if None in [name, keys, jobIDs]:
        print('ERROR:\t problem with input yaml file,', fn, ':', )
        print('name:', name)
        print('keys:', keys)
        print('jobIDs', jobIDs)
        exit(0)

    if isinstance(keys, str): keys = [keys, ]

    physics = False
    bio = False
    debug = False

    if 'physics' in keys:
        physics = True
    if 'bio'  in keys:
        bio = True
    if 'debug' in keys:
        debug = True
    if 'level1' in keys:
        physics = True
        bio = True

    randomcolours = [
                'red', 'blue', 'purple', 'green', 'gold', 'sienna', 'orange',
                'black', 'navy', 'goldenrod', 'dodgerblue',
            ]
    colours = {}
    for i, jobid in enumerate(sorted(yml['jobIDs'].keys())):
        colours[jobid] = yml['jobIDs'][jobid].get('colour', randomcolours[i])

    descriptions = {jobid:di.get('description', '') for jobid,di in yml['jobIDs'].items()}
    linewidths = {jobid:di.get('linewidth', 1.0) for jobid,di in yml['jobIDs'].items()}
    linestyles = {jobid:di.get('linestyle', '-') for jobid,di in yml['jobIDs'].items()}
    time_dict = {jobid:di.get('year_diff', 0.) for jobid,di in yml['jobIDs'].items()}

    for job in yml['jobIDs'].keys():
        print('plotting', job,  yml['jobIDs'][job])
    timeseries_compare(
        colours,
        analysisname=name,
        physics=physics,
        bio=bio,
        debug=debug,
        year0=time_dict, 
        jobDescriptions=descriptions,
        lineThicknesses=linewidths,
        linestyles=linestyles,
        config_user=config_user, 
        )


    exit(0)

    if 1:
        if suite == 'all':
            phys = 1
            bio = 1
            debug = 0
        if suite == 'physics':
            phys = 1
            bio = 0
            debug = 0
        if suite == 'bio':
            phys = 0
            bio = 1
            debug = 0
        if suite == 'debug':
            phys = 0
            bio = 0
            debug = 1
        try:
            colours = {i: standards[i] for i in jobIDs}
        except:
            colours = {}
            randomcolours = [
                'red', 'blue', 'purple', 'green', 'gold', 'sienna', 'orange',
                'black'
            ]
            for i, job in enumerate(jobIDs):
                colours[job] = randomcolours[i]
        name = '_'.join(jobIDs)
        timeseries_compare(colours,
                           physics=phys,
                           bio=bio,
                           year0=False,
                           debug=debug,
                           analysisname=name,
                           jobDescriptions=jobDescriptions)
        print("Successful command line comparison")

    return 


    #standards = configToDict('config/jobIDcolours.ini')
    thicknesses = defaultdict(lambda: 0.75)
    thicknesses['u-ar783'] = 2.2
    thicknesses['u-at793'] = 2.2
    #        thicknesses['u-at760'] = 2.2
    #        thicknesses['u-at628'] = 2.2
    #        thicknesses['u-at629'] = 2.2
    #/        thicknesses['u-at572'] = 2.2
    #        thicknesses['u-au027'] = 2.2

    hjthicknesses = defaultdict(lambda: 1.75)
    hjthicknesses['u-at793'] = 1.
    hjthicknesses['u-at760'] = 1.

    #jobDescriptions = configToDict('config/jobIDdescriptions.ini')
    #live_jobs = configToDict('RemoteScripts/jobids_config.ini')

    try:
        args = argv[1:]
        jobIDs = []
        suite = 'all'
        for job in args:
            if job == 'debug': suite = 'debug'
            elif job == 'physics': suite = 'physics'
            elif job in ['bgc', 'bio']: suite = 'bgc'
            else: jobIDs.append(job)

    except:
        jobsIDs = []
        suite = ''

    if len(jobIDs):
        print("analysis_compare.py:", jobIDs, suite)
        if suite == 'all':
            phys = 1
            bio = 1
            debug = 0
        if suite == 'physics':
            phys = 1
            bio = 0
            debug = 0
        if suite == 'bio':
            phys = 0
            bio = 1
            debug = 0
        if suite == 'debug':
            phys = 0
            bio = 0
            debug = 1
        try:
            colours = {i: standards[i] for i in jobIDs}
        except:
            colours = {}
            randomcolours = [
                'red', 'blue', 'purple', 'green', 'gold', 'sienna', 'orange',
                'black'
            ]
            for i, job in enumerate(jobIDs):
                colours[job] = randomcolours[i]
        name = '_'.join(jobIDs)
        timeseries_compare(colours,
                           physics=phys,
                           bio=bio,
                           year0=False,
                           debug=debug,
                           analysisname=name,
                           jobDescriptions=jobDescriptions)
        print("Successful command line comparison")
        exit
    else:

        UKESM11_fast = True
        if UKESM11_fast:

            customColours = {
                'u-by230': 'black',  #standard UKESM1.1
                'u-ck416': 'green',  #CN-fast standard radiation
                'u-cg799' : 'red', #CN-fast rad=3h, two_fsd=1.7
                # 'u-cg843' : 'blue', #CN-fast, rad=3h, two_fsd=1.65
            }
            cnthicknesses = defaultdict(lambda: 0.7)
            linestyles = defaultdict(lambda: '-')
            time_dict = defaultdict(lambda: None)
            time_dict['u-by230'] = 50.
            time_dict['u-ck416'] = -50.


            timeseries_compare(
                customColours,
                physics=1,
                bio=1,
                debug=0,
                year0=time_dict, #'UKESM11_Fast_piControl_2',
                #obDescriptions=defaultdict(lambda: {} ),
                analysisname='UKESM11_Fast_piControl_2',
                lineThicknesses=cnthicknesses,
                linestyles=linestyles,
            )

        """
        Below here is the old system.
        """
        return
        assert 0

        UKESM11_scenarios = False
        if UKESM11_scenarios:
            customColours = {
                'u-cb261':
                'green',  #                 :  UKESM1.1 ScenarioMIP ssp126 1
                'u-cb584':
                'blue',  #                  : UKESM1.1 ScenarioMIP ssp126 2
                'u-cb586':
                'purple',  #                  : UKESM1.1 ScenarioMIP ssp126 3
                'u-cb180':
                'orange',  #                  : UKESM1.1 ScenarioMIP ssp370 1
                'u-cb581':
                'red',  #                  : UKESM1.1 ScenarioMIP ssp370 2
                'u-cb585':
                'gold',  #                  : UKESM1.1 ScenarioMIP  ssp370 3
            }
            cnthicknesses = defaultdict(lambda: 0.7)
            linestyles = defaultdict(lambda: '-')
            timeseries_compare(
                customColours,
                physics=1,
                bio=1,
                debug=0,
                year0=False,
                jobDescriptions=jobDescriptions,
                analysisname='UKESM11_scenarios',
                lineThicknesses=cnthicknesses,
                linestyles=linestyles,
            )

        UKESM11_pic = True
        if UKESM11_pic:
            customColours = {
                'u-bx188': 'black',
                'u-cb375': 'red',
                'u-bz866': 'blue',
                'u-cb737': 'orange',
            }
            cnthicknesses = defaultdict(lambda: 0.7)
            linestyles = defaultdict(lambda: '-')
            timeseries_compare(
                customColours,
                physics=1,
                bio=1,
                debug=0,
                year0=False,
                jobDescriptions=jobDescriptions,
                analysisname='UKESM11_pic_compare',
                lineThicknesses=cnthicknesses,
                linestyles=linestyles,
            )
        UKESM_fast = False
        if UKESM_fast:
            customColours = {
                'u-aw310': 'black',
                'u-by242': 'red',
            }
            cnthicknesses = defaultdict(lambda: 0.7)
            linestyles = defaultdict(lambda: '-')
            timeseries_compare(
                customColours,
                physics=1,
                bio=1,
                debug=0,
                year0=False,
                jobDescriptions=jobDescriptions,
                analysisname='UKESM_Fast_piControl',
                lineThicknesses=cnthicknesses,
                linestyles=linestyles,
            )

            customColours = {
                'u-bc179': 'black',
                'u-bc292': 'black',
                'u-bc370': 'black',
                'u-bb075': 'black',
                'u-az513': 'black',
                'u-az515': 'black',
                'u-bz321': 'red',
                'u-bz360': 'red',
                'u-bz361': 'red',
                'u-bz362': 'red',
                'u-bz363': 'red',
                'u-bz364': 'red',
            }
            cnthicknesses = defaultdict(lambda: 0.7)
            linestyles = defaultdict(lambda: '-')
            timeseries_compare(
                customColours,
                physics=1,
                bio=1,
                debug=0,
                year0=False,
                jobDescriptions=jobDescriptions,
                analysisname='UKESM_Fast_historical',
                lineThicknesses=cnthicknesses,
                linestyles=linestyles,
            )

            customColours = { # u-bb446 vs. u-bz385
                    'u-bb446': 'black',
                    'u-bz385': 'red',
                    }
            cnthicknesses = defaultdict(lambda: 0.7)
            linestyles = defaultdict(lambda: '-')
            timeseries_compare(
                customColours,
                physics=1,
                bio=1,
                debug=0,
                year0=False,
                jobDescriptions=jobDescriptions,
                analysisname='UKESM_Fast_4xCO2',
                lineThicknesses=cnthicknesses,
                linestyles=linestyles,
            )

            customColours = { # u-bb448 vs. u-bz393
                    'u-bb448': 'black',
                    'u-bz393': 'red',
                    }
            cnthicknesses = defaultdict(lambda: 0.7)
            linestyles = defaultdict(lambda: '-')
            timeseries_compare(
                customColours,
                physics=1,
                bio=1,
                debug=0,
                year0=False,
                jobDescriptions=jobDescriptions,
                analysisname='UKESM_Fast_1pcCO2',
                lineThicknesses=cnthicknesses,
                linestyles=linestyles,
            )

#                UKESM_fast = True
#                if UKESM_fast:
#                        customColours = {
#                                'u-bz499': 'black',
#                                'u-bz705': 'red',
#                                }
#                        cnthicknesses = defaultdict(lambda: 0.7)
#                        linestyles = defaultdict(lambda: '-')
#                        timeseries_compare(
#                                customColours,
#                                physics=1,
#                                bio=1,
#                                debug=0,
#                                year0=False,
#                                jobDescriptions=jobDescriptions,
#                                analysisname='UKESM_Fast',
#                                lineThicknesses= cnthicknesses,
#                                linestyles = linestyles,)
#
#

        UKESM11_historical2 = False
        if UKESM11_historical2:
            customColours = {
                'u-by230': 'black',
                'u-by791': 'green',
                'u-bz502': 'red',
                'u-bz897': 'blue',
                'u-ca306': 'purple',
                'u-ca811': 'orange',
                'u-cb799': 'goldenrod',  # ukesm1.1 histo #6
                'u-bc179': 'black',  # R1 historical
                'u-bc292': 'black',  # R2 hist
                'u-bc370': 'black',  # R3 histroccal
                'u-bd483': 'black',  # R? hist
                'u-bc470': 'black',  # R? hist
                'u-bd288': 'black',  # R? hist
                # 'u-bb075': 'black', # R4 hist
                # 'u-bb277': 'black', # R5 hist
            }

            cnthicknesses = defaultdict(lambda: 0.7)
            cnthicknesses['u-by791'] = 2.
            cnthicknesses['u-bz502'] = 2.
            cnthicknesses['u-bz897'] = 2.
            cnthicknesses['u-ca306'] = 2.
            cnthicknesses['u-ca811'] = 2.
            cnthicknesses['u-cb799'] = 2.
            cnthicknesses['u-by230'] = 1.5

            linestyles = defaultdict(lambda: '-')
            linestyles['u-by230'] = ':'
            timeseries_compare(
                customColours,
                physics=1,
                bio=1,
                debug=0,
                year0='UKESM11_historical1',
                jobDescriptions=jobDescriptions,
                analysisname='UKESM11_historical2',
                lineThicknesses=cnthicknesses,
                linestyles=linestyles,
            )
        UKESM11_historical = False
        if UKESM11_historical:
            customColours = {
                'u-by230': 'black',
                'u-by791': 'green',
                'u-bz502': 'red',
                'u-bz897': 'blue',
                'u-ca306': 'purple',
                'u-ca811': 'orange',
                'u-ca730': 'goldenrod',  # ukesm1.1 histo #6
                'u-bc179': 'black',  # R1 historical
                'u-bc292': 'black',  # R2 hist
                'u-bc370': 'black',  # R3 histroccal
                'u-bb075': 'black',  # R4 hist
                'u-bb277': 'black',  # R5 hist
            }
            cnthicknesses = defaultdict(lambda: 0.7)
            cnthicknesses['u-by791'] = 2.
            cnthicknesses['u-bz502'] = 2.
            cnthicknesses['u-bz897'] = 2.
            cnthicknesses['u-ca306'] = 2.
            cnthicknesses['u-ca811'] = 2.
            cnthicknesses['u-by230'] = 1.5

            linestyles = defaultdict(lambda: '-')
            linestyles['u-by230'] = ':'
            timeseries_compare(
                customColours,
                physics=1,
                bio=1,
                debug=0,
                year0='UKESM11_historical1',
                jobDescriptions=jobDescriptions,
                analysisname='UKESM11_historical1',
                lineThicknesses=cnthicknesses,
                linestyles=linestyles,
            )

        UKESM_vs_HadGEM3 = True
        if UKESM_vs_HadGEM3:
            customColours = {
                'u-aw310': 'black',  # UKESM1 piControl
                'u-ar766': 'green',  # HadGEM3 piControl:
                'u-bb446': 'blue',  # UKESM1 abrupt-4xCO2
                'u-bg555': 'orange',  # HadGEM3 abrupt-4xCO2
            }
            cnthicknesses = defaultdict(lambda: 1.1)
            linestyles = defaultdict(lambda: '-')
            linestyles['u-aw310'] = ':'
            linestyles['u-ar766'] = ':'
            timeseries_compare(
                customColours,
                physics=1,
                bio=1,
                debug=0,
                year0=None,
                jobDescriptions=jobDescriptions,
                analysisname='UKESM_vs_HadGEM3',
                lineThicknesses=cnthicknesses,
                linestyles=linestyles,
            )

        UKESM_fast2 = True
        if UKESM_fast2:
            picontrol = {}
            Co2_4x = {}
            year0 = {}
            picontrol['N96L85'] = 'u-aw310'
            #picontrol['N48L85'] = 'u-bk171'
            #picontrol['N48L38'] = 'u-bs733'
            Co2_4x['N96L85'] = 'u-bb446'
            Co2_4x['N48L85'] = 'u-bk093'
            Co2_4x['N48L38'] = 'u-bt089'
            year0['N96L85'] = 'AlignToDECK2100'
            year0['N48L85'] = None
            year0['N48L38'] = 'AlignToDECK2300'
            others = {
                'N48L38': {
                    'u-bt670': 'blue',
                    'u-bu847': 'green',
                    'u-bv334': 'purple',
                    'u-bv270': 'orange'
                },
            }
            for configuration in list(picontrol.keys()):
                customColours = {
                    picontrol[configuration]: 'black',
                    Co2_4x[configuration]: 'red',
                }
                if configuration in list(others.keys()):
                    customColours.update(others[configuration])

                cnthicknesses = defaultdict(lambda: 1.1)
                linestyles = defaultdict(lambda: '-')
                linestyles[picontrol[configuration]] = ':'
                timeseries_compare(
                    customColours,
                    physics=1,
                    bio=0,
                    debug=0,
                    year0=year0[configuration],
                    jobDescriptions=jobDescriptions,
                    analysisname='UKESM_fast_' + configuration,
                    lineThicknesses=cnthicknesses,
                    linestyles=linestyles,
                )

        UKESM_fast = False
        if UKESM_fast:
            customColours = {
                'u-aw310': 'black',
                'u-bs733': 'red',
            }

            cnthicknesses = defaultdict(lambda: 1.1)
            linestyles = defaultdict(lambda: '-')
            linestyles['u-aw310'] = ':'
            timeseries_compare(
                customColours,
                physics=1,
                bio=0,
                debug=0,
                year0='AlignToDECK2100',
                jobDescriptions=jobDescriptions,
                analysisname='UKESM_fast_final',
                lineThicknesses=cnthicknesses,
                linestyles=linestyles,
            )

        UKESM11_picontrol = False
        if UKESM11_picontrol:
            customColours = {
                #'u-bs463': 'green',
                'u-aw310': 'black',
                'u-bs522': 'blue',
                'u-bs704': 'red',
                'u-bt233': 'green',
                'u-bt320': 'purple',
                'u-bt931': 'orange',
                'u-bu504': 'brown',
            }
            descripts = {
                'u-aw310': 'UKESM1.0',
                # 'u-bs463' : 'UKESM1.0.1',
                'u-bs522': 'UKESM1.0.1 + snow mods',
                'u-bs704': 'UKESM1.0.2',
                'u-bt233': '',
                'u-bt320': '',
                'u-bt931':
                'continuation of u-bt233 with a slight tweak on TOA radiation',
                'u-bu504': '',
            }
            cnthicknesses = defaultdict(lambda: 1.1)
            linestyles = defaultdict(lambda: '-')
            linestyles['u-aw310'] = ':'
            timeseries_compare(
                customColours,
                physics=1,
                bio=1,
                debug=0,
                year0='2300-2700',
                jobDescriptions=descripts,
                analysisname='UKESM11_picontrol',
                lineThicknesses=cnthicknesses,
                linestyles=linestyles,
            )

        UKESM11_picontrol2 = False
        if UKESM11_picontrol2:
            customColours = {
                'u-bt233': 'green',
                'u-bu504': 'red',
                'u-bu737': 'purple',
                'u-bu794': 'gold',
                'u-bv334': 'blue',
                #'u-bv270': 'orange',
            }
            descripts = {
                'u-bt233': 'UKESM1.1 pi control',
                'u-bu504':
                'UKESM1.1  - based on u-bt233 + UKCA & JULES temporary logicals activated',
                'u-bu737':
                'UKESM1.1  - based on u-bt233 + sea ice bugfix #661 included',
                'u-bu794':
                'UKESM1.1  - based on u-bt233 + CDNC correction included',
                'u-bv334': 'UKESM 1.0.4 - two_fsd  = 1.56, ICs from aw310',
                #'u-bv270': 'two_fsd = 1.56 ICs from bt670',
            }

            cnthicknesses = defaultdict(lambda: 1.1)
            linestyles = defaultdict(lambda: '-')
            timeseries_compare(
                customColours,
                physics=1,
                bio=1,
                debug=0,
                year0=None,
                jobDescriptions=descripts,
                analysisname='UKESM11_picontrol_2',
                lineThicknesses=cnthicknesses,
                linestyles=linestyles,
            )

        UKESM11_picontrol3 = True
        if UKESM11_picontrol3:
            customColours = {
                'u-aw310': 'black',
                # 'u-bv334': 'red',
                'u-bv936': 'green',
                'u-bw462': 'blue',
                'u-bw837': 'orange',
                'u-bx188': 'darkorange',
                'u-bx499': 'magenta',
                'u-by230': 'purple',
                'u-ca730': 'red',  # ukesm1.1 histo #6
            }
            descripts = {
                'u-aw310': 'UKESM1',
                # 'u-bv334': 'UKESM1.0.4',
                'u-bv936': 'UKESM1.0.5.0',
                'u-bw462': 'UKESM1.0.5.1',
                'u-bw837': '',
                'u-bx188': '',
                'u-bx499':
                'Variant of bx188 with scaling on dust to iron for MEDUSA, parameter xfe_sol * 0.96',
                'u-by230': 'our "chilled" UKESM1.1 configuration',
                'u-ca730': '',  # ukesm1.1 histo #6
            }

            cnthicknesses = defaultdict(lambda: 1.1)
            linestyles = defaultdict(lambda: '-')
            timeseries_compare(
                customColours,
                physics=1,
                bio=1,
                debug=0,
                year0='2300-2900',
                jobDescriptions=descripts,
                analysisname='UKESM11_picontrol_3',
                lineThicknesses=cnthicknesses,
                linestyles=linestyles,
            )

        N96L85 = False
        if N96L85:
            customColours = {
                'u-bk712': 'green',
                'u-bk575': 'red',
                'u-br173': 'orange',
                'u-bm585': 'mediumpurple',
                'u-bk741': 'blue',
                'u-bm448': 'deeppink',
            }
            cnthicknesses = defaultdict(lambda: 1.1)
            linestyles = defaultdict(lambda: '-')
            linestyles['u-aw310'] = ':'
            timeseries_compare(
                customColours,
                physics=1,
                bio=1,
                debug=0,
                year0='1950-2020',
                jobDescriptions=jobDescriptions,
                analysisname=
                'N96L85_mod_snow_albedo_701_parents_without_snow_mod',
                lineThicknesses=cnthicknesses,
                linestyles=linestyles,
            )

        N96L85 = False
        if N96L85:
            customColours = {
                'u-br106': 'blue',
                'u-br173': 'orange',
                'u-br894': 'purple',
                'u-br896': 'red',
            }
            cnthicknesses = defaultdict(lambda: 1.1)
            linestyles = defaultdict(lambda: '-')
            linestyles['u-aw310'] = ':'
            timeseries_compare(
                customColours,
                physics=1,
                bio=1,
                debug=0,
                year0='1950-2020',
                jobDescriptions=jobDescriptions,
                analysisname='N96L85_mod_snow_albedo_701_Difference_mods',
                lineThicknesses=cnthicknesses,
                linestyles=linestyles,
            )

        N96L85 = False
        if N96L85:
            customColours = {
                'u-aw310': 'black',
                'u-bk575': 'red',
                'u-br173': 'orange',
                'u-bk712': 'green',
                'u-br106': 'blue',
                'u-bm585': 'mediumpurple',
                'u-bm585': 'dodgerblue',
                'u-br894': 'purple',
                'u-bm448': 'deeppink',
                'u-br896': 'goldenrod',
            }
            cnthicknesses = defaultdict(lambda: 1.1)
            linestyles = defaultdict(lambda: '-')
            linestyles['u-aw310'] = ':'
            timeseries_compare(
                customColours,
                physics=1,
                bio=1,
                debug=0,
                year0='AlignToDECK2200',
                jobDescriptions=jobDescriptions,
                analysisname='N96L85_mod_snow_albedo_701',
                lineThicknesses=cnthicknesses,
                linestyles=linestyles,
            )

        extensions = False
        if extensions:
            r4_colours = {
                'u-bb075': 'black',
                'u-be393': 'blue',
                'u-be394': 'purple',
                'u-be335': 'pink',
                'u-be392': 'red',
                'u-bh210': 'green',
                'u-bh254': 'gold',
                'u-bh285': 'orange',
            }

            timeseries_compare(r4_colours,
                               physics=1,
                               bio=1,
                               debug=0,
                               year0=False,
                               jobDescriptions=jobDescriptions,
                               analysisname='UKESM_scenarios_extensions',
                               lineThicknesses=hjthicknesses)

        UKESM_fast = False
        if UKESM_fast:
            customColours = {
                'u-aw310': 'black',
                'u-bp531': 'red',
                'u-bp705': 'purple',
            }

            cnthicknesses = defaultdict(lambda: 1.1)
            linestyles = defaultdict(lambda: '-')
            linestyles['u-aw310'] = ':'
            timeseries_compare(
                customColours,
                physics=1,
                bio=0,
                debug=0,
                year0='AlignToDECK2100',
                jobDescriptions=jobDescriptions,
                analysisname='UKESM_fast_semiofficial',
                lineThicknesses=cnthicknesses,
                linestyles=linestyles,
            )

        UKESM_fast2 = False
        if UKESM_fast2:
            customColours = {
                'u-aw310': 'black',
                #'u-bp531': 'red',
                #'u-bp705': 'purple',
                'u-bq621': 'green',
            }
            cnthicknesses = defaultdict(lambda: 1.1)
            linestyles = defaultdict(lambda: '-')
            linestyles['u-aw310'] = ':'
            timeseries_compare(
                customColours,
                physics=1,
                bio=0,
                debug=0,
                year0='AlignToDECK2020',
                jobDescriptions=jobDescriptions,
                analysisname='UKESM_fast_piControl',
                lineThicknesses=cnthicknesses,
                linestyles=linestyles,
            )

        UKESM_fast3 = False
        if UKESM_fast3:
            customColours = {
                'u-aw310': 'black',
                #'u-bp531': 'red',
                #'u-bp705': 'purple',
                'u-bq621': 'green',
                'u-bq834': 'red',
                'u-bb446': 'teal',
            }
            cnthicknesses = defaultdict(lambda: 1.1)
            linestyles = defaultdict(lambda: '-')
            linestyles['u-aw310'] = ':'
            timeseries_compare(
                customColours,
                physics=1,
                bio=0,
                debug=0,
                year0='AlignToDECK2020',
                jobDescriptions=jobDescriptions,
                analysisname='UKESM_fast_Abrupt4x',
                lineThicknesses=cnthicknesses,
                linestyles=linestyles,
            )

        UKESM_fast4 = True
        if UKESM_fast4:
            customColours = {
                'u-aw310': 'black',
                'u-bx082': 'red',
                'u-bw717': 'purple',
                #                          'u-bq621': 'green',
                #                          'u-bq834': 'red',
                #                          'u-bb446': 'teal',
            }
            cnthicknesses = defaultdict(lambda: 1.1)
            linestyles = defaultdict(lambda: '-')
            linestyles['u-aw310'] = ':'
            timeseries_compare(
                customColours,
                physics=1,
                bio=0,
                debug=0,
                year0='fast4',
                jobDescriptions=jobDescriptions,
                analysisname='UKESM_fast_4',
                lineThicknesses=cnthicknesses,
                linestyles=linestyles,
            )

        GeoengineeringFudge = False
        if GeoengineeringFudge:
            #GeoengineeringFudge
            #u-bj141			: control, emissions-driven esm-ssp585
            #u-bo540			: experiment 1, default geoengineering, stops at 2050
            #u-bp104			: experiment 2, modified geoengineering, no stop
            #jobs = ['u-bj141', 'u-bo540', 'u-bp104', ]
            customColours = {
                'u-bj141': 'black',
                'u-bo540': 'red',
                'u-bp104': 'blue',
            }
            cnthicknesses = defaultdict(lambda: 1.1)
            linestyles = defaultdict(lambda: '-')
            linestyles['u-bj141'] = ':'
            timeseries_compare(
                customColours,
                physics=1,
                bio=1,
                debug=0,
                year0=False,  #'AlignToDECK2020',
                jobDescriptions=jobDescriptions,
                analysisname='GeoengineeringFudge',
                lineThicknesses=cnthicknesses,
                linestyles=linestyles,
            )

        N48ORCA1 = False
        if N48ORCA1:
            customColours = {
                'u-aw310': 'black',
                #u-bc179': 'red',
                'u-bp179': 'green',
                'u-bk171': 'blue',
                #'u-bo523': 'green',
                'u-bp360': 'orange',
            }
            cnthicknesses = defaultdict(lambda: 1.1)
            linestyles = defaultdict(lambda: '-')
            linestyles['u-aw310'] = ':'
            timeseries_compare(
                customColours,
                physics=1,
                bio=1,
                debug=0,
                year0='AlignToDECK2100',
                jobDescriptions=jobDescriptions,
                analysisname='Proto_UKESM1_fast',
                lineThicknesses=cnthicknesses,
                linestyles=linestyles,
            )

        N48ORCA1 = False
        if N48ORCA1:
            customColours = {
                'u-aw310': 'black',
                'u-bc179': 'black',
                'u-bm926': 'lightcoral',
                'u-bo523': 'seagreen',
            }
            cnthicknesses = defaultdict(lambda: 1.1)
            linestyles = defaultdict(lambda: '-')
            linestyles['u-aw310'] = ':'
            linestyles['u-bm926'] = ':'
            timeseries_compare(
                customColours,
                physics=1,
                bio=1,
                debug=0,
                year0='AlignToDECK2020',
                jobDescriptions=jobDescriptions,
                analysisname='N48ORCA1',
                lineThicknesses=cnthicknesses,
                linestyles=linestyles,
            )

        N96L85 = False
        if N96L85:
            jobs = [
                'u-aw310', 'u-bc179', 'u-bn824', 'u-bk437', 'u-bk575',
                'u-bk712', 'u-bk713', 'u-bk741'
            ]
            customColours = {
                'u-aw310': 'black',
                'u-bc179': 'black',
                'u-bn824': 'navy',
                'u-bk437': 'navy',
                'u-bk575': 'red',
                'u-bk712': 'green',
                'u-bk713': 'purple',
                'u-bk741': 'orange',
            }
            cnthicknesses = defaultdict(lambda: 1.1)
            linestyles = defaultdict(lambda: '-')
            linestyles['u-aw310'] = ':'
            linestyles['u-bn824'] = ':'
            linestyles['u-bk437'] = ':'
            timeseries_compare(
                customColours,
                physics=1,
                bio=1,
                debug=0,
                year0='AlignToDECK2200',
                jobDescriptions=jobDescriptions,
                analysisname='N96L85',
                lineThicknesses=cnthicknesses,
                linestyles=linestyles,
            )

        return

        ocean_only = 0
        if ocean_only:
            jobs = [
                'u-bc179',
                'u-bc292',
                'u-bc370',
                'u-bg720',
                'u-bg742',
                'u-bg743',
            ]
            customColours = {
                'u-bc179': 'lightcoral',
                'u-bc292': 'seagreen',
                'u-bc370': 'dodgerblue',
                'u-bg720': 'red',
                'u-bg742': 'green',
                'u-bg743': 'blue',
            }
            descr = {
                'u-bc179':
                jobDescriptions['u-bc179'],
                'u-bc292':
                jobDescriptions['u-bc292'],
                'u-bc370':
                jobDescriptions['u-bc370'],
                'u-bg720':
                'AerchemMIP experiments branched from UKESM1 historical (u-bc179) with CFC and HCFC emissions fixed at 1950 values.',
                'u-bg742':
                'AerchemMIP experiments branched from UKESM1 historical (u-bc292) with CFC and HCFC emissions fixed at 1950 values.',
                'u-bg743':
                'AerchemMIP experiments branched from UKESM1 historical (u-bc370) with CFC and HCFC emissions fixed at 1950 values.',
            }
            cnthicknesses = defaultdict(lambda: 1.1)
            linestyles = defaultdict(lambda: '-')
            timeseries_compare(
                customColours,
                physics=1,
                bio=1,
                debug=0,
                year0=None,
                jobDescriptions=descr,
                analysisname='AerchemMIP_with_CFC_HCFC_fixed_1950',
                lineThicknesses=cnthicknesses,
                linestyles=linestyles,
            )

        faster_UKESM = True
        if faster_UKESM:
            jobs = [
                'u-aw310',
                'u-bk171',
                'u-bj890',
                'u-bl696',
            ]
            customColours = {
                'u-aw310': 'black',
                'u-bk171': 'red',
                'u-bj890': 'blue',
                'u-bl696': 'green',
            }
            cnthicknesses = defaultdict(lambda: 1.1)
            linestyles = defaultdict(lambda: '-')
            timeseries_compare(
                customColours,
                physics=1,
                bio=1,
                debug=0,
                year0='Historical',
                jobDescriptions=jobDescriptions,
                analysisname='Progressively_faster_UKESM',
                lineThicknesses=cnthicknesses,
                linestyles=linestyles,
            )

        a = 1
        if a:
            newemissionscolours = {
                'u-av651': 'black',
                'u-aw310': 'black',
                'u-az021': 'darkgreen',
                'u-az417': 'rebeccapurple',
                'u-az418': 'dodgerblue',
                'u-az513': 'darkorange',  # UKESM1 fifth Historical run (2020)
                'u-az515':
                'deeppink',  #      UKESM1 sixth Historical run (2050)
                'u-az524': 'goldenrod',  # UKESM1 seventh Historical run (1995)
                'u-bb075': 'green',
                'u-bb277': 'dodgerblue',
                'u-bc179':
                'pink',  #                  :standard ; UKESM1 Historical run (2250)
                'u-bc292':
                'brown',  #                :standard ; UKESM1 Historical run(2165)
                'u-bc370':
                'slateblue',  #                :standard ; UKESM1 Historical run(2120)
                'u-bc470':
                'gold',  #                :standard ; UKESM1 Historical run(2120)
                'u-bd288': 'orchid',  # UKESM1 Historical run (2340)
                'u-bd416': 'navy',  # UKESM1 Historical run (2460)
                'u-bd483': 'olive',  # UKESM1 Historical run (2200)
                #					'u-bf935': 'fushia', # UKESM1 Historical run '2565', ],
                'u-bh100': 'yellow',  # UKESM1 Historical run '2685', ],
                'u-bh101': 'darkgold',  # UKESM1 Historical run '2745', ],
                'u-bf647': 'lightgreen',  # UKESM1 Historical run '2629', ],
                'u-bf656': 'cyan',  # UKESM1 Historical run '2716', ],
                'u-bf703': 'crimson',  # UKESM1 Historical run '2760', ],
                'u-bh162': 'lavender',  # UKESM1 Historical run '2815', ],}
            }
            linestyles = defaultdict(lambda: '-')
            linestyles['u-av651'] = '--'
            linestyles['u-aw310'] = ':'
            cnthicknesses = defaultdict(lambda: 1.1)
            cnthicknesses['u-aw310'] = 1.7
            timeseries_compare(
                newemissionscolours,
                physics=1,
                bio=1,
                debug=0,
                year0='ControlAligned',
                jobDescriptions=jobDescriptions,
                analysisname='UKESM1_ControlAligned',
                lineThicknesses=cnthicknesses,
                linestyles=linestyles,
            )

        pi_aerosols = False
        if pi_aerosols:
            jobs = [
                'u-bg705', 'u-bg723', 'u-bg724', 'u-bc179', 'u-bc292',
                'u-bc370'
            ]
            customColours = {
                'u-bg705': 'orange',
                'u-bg723': 'orange',
                'u-bg724': 'orange',
                'u-bc179': 'black',
                'u-bc292': 'black',
                'u-bc370': 'black',
            }
            cnthicknesses = defaultdict(lambda: 1.1)
            linestyles = defaultdict(lambda: '-')
            linestyles['u-bc179'] = '--'
            linestyles['u-bc292'] = '--'
            linestyles['u-bc370'] = '--'
            timeseries_compare(
                customColours,  #{i:standards[i] for i in jobs},
                physics=1,
                bio=1,
                debug=0,
                year0=False,
                jobDescriptions=jobDescriptions,
                analysisname='UKESM1_pi_aerosols',
                lineThicknesses=cnthicknesses,
                linestyles=linestyles,
            )

        scenarios = 0
        if scenarios:
            r1_colours = {
                'u-bc179': 'black',
                'u-be509': 'blue',
                'u-be537': 'purple',
                'u-be647': 'pink',
                'u-be653': 'red',
                'u-bh409': 'green',
                'u-bh454': 'gold',
                'u-bh456': 'orange',
            }
            r2_colours = {
                'u-bc292': 'black',
                'u-be679': 'blue',
                'u-be606': 'purple',
                'u-be690': 'pink',
                'u-be693': 'red',
                'u-bh570': 'green',
                'u-bh712': 'gold',
                'u-bh744': 'orange',
            }
            r3_colours = {
                'u-bc370': 'black',
                'u-be682': 'blue',
                'u-be683': 'purple',
                'u-be684': 'pink',
                'u-be686': 'red',
                'u-bh716': 'green',
                'u-bh717': 'gold',
                'u-bh718': 'orange',
            }

            r4_colours = {
                'u-bb075': 'black',
                'u-be393': 'blue',
                'u-be394': 'purple',
                'u-be335': 'pink',
                'u-be392': 'red',
                'u-bh210': 'green',
                'u-bh254': 'gold',
                'u-bh285': 'orange',
            }
            r8_colours = {
                'u-bb277': 'black',
                'u-be397': 'blue',
                'u-be398': 'purple',
                'u-be395': 'pink',
                'u-be396': 'red',
                'u-bh807': 'green',
                'u-bh808': 'gold',
                'u-bh809': 'orange',
            }
            run_names = {
                'Run1': r1_colours,
                'Run2': r2_colours,
                'Run3': r3_colours,
                'Run4': r4_colours,
                'Run8': r8_colours
            }
            #                        run_dicts = [: r1_colours, r2_colours, r3_colours, r4_colours, r5_colours]
            #                        run_names = ['Run1', 'Run2', 'Run3', 'Run4', 'Run5']
            completed_runs = []  #'Run1', 'Run2', 'Run3', 'Run4', 'Run8']

            for run_name, run_colours, in list(run_names.items()):

                if run_name not in completed_runs:
                    timeseries_compare(run_colours,
                                       physics=1,
                                       bio=1,
                                       debug=0,
                                       year0=False,
                                       jobDescriptions=jobDescriptions,
                                       analysisname='UKESM_scenarios_' +
                                       run_name,
                                       lineThicknesses=hjthicknesses)

        scenarios_figures = 0
        if scenarios_figures:
            ssp126 = {
                'r1': 'u-be509',
                'r2': 'u-be679',
                'r3': 'u-be682',
                'r4': 'u-be393',
                'r8': 'u-be397'
            }
            ssp245 = {
                'r1': 'u-be537',
                'r2': 'u-be606',
                'r3': 'u-be683',
                'r4': 'u-be394',
                'r8': 'u-be398'
            }
            ssp370 = {
                'r1': 'u-be647',
                'r2': 'u-be690',
                'r3': 'u-be684',
                'r4': 'u-be335',
                'r8': 'u-be395'
            }
            ssp585 = {
                'r1': 'u-be653',
                'r2': 'u-be693',
                'r3': 'u-be686',
                'r4': 'u-be392',
                'r8': 'u-be396'
            }

            # Tier2:
            ssp119 = {
                'r1': 'u-bh409',
                'r2': 'u-bh570',
                'r3': 'u-bh716',
                'r4': 'u-bh210',
                'r8': 'u-bh807'
            }
            ssp434 = {
                'r1': 'u-bh454',
                'r2': 'u-bh724',
                'r3': 'u-bh717',
                'r4': 'u-bh254',
                'r8': 'u-bh808'
            }
            ssp534 = {
                'r1': 'u-bh456',
                'r2': 'u-bh744',
                'r3': 'u-bh718',
                'r4': 'u-bh285',
                'r8': 'u-bh809'
            }

            run_colours = {
                'r1': 'red',
                'r2': 'orange',
                'r3': 'goldenrod',
                'r4': 'green',
                'r8': 'blue',
            }

            #scenarios = {'ssp126': ssp126,  'ssp245': ssp245, 'ssp370': ssp370, 'ssp585': ssp585, 'ssp119':ssp119, 'ssp434':ssp434, 'ssp534':ssp534}
            scenarios = {'ssp119': ssp119, 'ssp434': ssp434, 'ssp534': ssp534}
            for scean, sc_list in list(scenarios.items()):
                sc_colours = {
                    jobID: run_colours[run_i]
                    for run_i, jobID in list(sc_list.items())
                }
                sc_jobDescriptions = {
                    jobID: ' '.join(['UKESM', scean,
                                     run_i.upper()])
                    for run_i, jobID in list(sc_list.items())
                }
                timeseries_compare(sc_colours,
                                   physics=1,
                                   bio=1,
                                   debug=0,
                                   year0=False,
                                   jobDescriptions=sc_jobDescriptions,
                                   analysisname='UKESM_scenario_' + scean,
                                   lineThicknesses=hjthicknesses)

#			runs = run_dicts[0].copy()
#			for r_dict in run_dicts:
#	                        runs = dict(runs, **r_dict)
#		        linestyles = defaultdict(lambda: '-')
#		        for k in runs.keys():
#		        	if k in r4_colours: linestyles[k] = '-'
#		        	if k in r5_colours: linestyles[k] = ':'
#                        timeseries_compare(
#                                 runs,
#                                 physics=1,
#                                 bio=1,
#                                 debug=0,
#                                 year0=False,
#                                 jobDescriptions=jobDescriptions,
#                                 analysisname='UKESM_scenarios_r4_r5',
#                                 lineThicknesses= hjthicknesses,
#                                 linestyles = linestyles,)

        a = False
        if a:
            cncolours = {
                'u-aw310': 'black',
                'u-az508': 'purple',
                'u-bb458': 'red',
            }
            timeseries_compare(cncolours,
                               physics=1,
                               bio=1,
                               debug=0,
                               year0='control_2100',
                               jobDescriptions=jobDescriptions,
                               analysisname='UKESM_control_2100',
                               lineThicknesses=hjthicknesses)

#####
# plot an individual job against the pi control.
        customColours = {
            'u-bb075':
            'teal',  # UKESM1 Historical run (1960) with new SO2 emissions height
            'u-az524':
            'green',  # UKESM1 Historical run (1995) with new SO2 emissions height
            'u-az513':
            'blue',  # UKESM1 Historical run (2020) with new SO2 emissions height
            'u-az515':
            'purple',  # UKESM1 Historical run (2050) with new SO2 emissions height
            'u-bb277':
            'orange',  # UKESM1 Historical run (2395) with new SO2 emissions height
            'u-bb446':
            'pink',  # UKESM1 4xCO2 run (1960) with new SO2 emissions height
            'u-bb448':
            'gold',  # UKESM1 1%CO2 run (1960) with new SO2 emissions height
            'u-bc179':
            'dodgerblue',  #                  :standard ; UKESM1 Historical run (2250)
            'u-bc292':
            'brown',  #                :standard ; UKESM1 Historical run(2165)
            'u-bc370': 'slateblue',  # UKESM1 Historical run (2120)
            'u-bc470': 'gold',  # UKESM1 Historical run (2285)
            'u-bd288': 'orchid',  # UKESM1 Historical run (2340)
            'u-bd416': 'navy',  # UKESM1 Historical run (2460)
            'u-bd483': 'olive',  # UKESM1 Historical run (2460)
            'u-bf647': 'red',  # ;UKESM1 Historical run (2619)
            'u-bf656': 'black',  # ;UKESM1 Historical run (2716)
            'u-bf703': 'darkgoldenrod',  # ;UKESM1 Historical run (2760)
            # 'u-bf705': 'yellow',	# ;UKESM1 Historical run (2815) = dead
            'u-bh162': 'yellow',  # ;UKESM1 Historical run (2815)
        }
        start_year = {
            'u-bb075': 1960,
            'u-az524': 1995,
            'u-az513': 2020,
            'u-az515': 2050,
            'u-bb277': 2395,
            'u-bb446': 1960,
            'u-bb448': 1960,
            'u-bc179': 2250,
            'u-bc292': 2165,
            'u-bc370': 2120,
            'u-bc470': 2285,
            'u-bd288': 2340,
            'u-bd416': 2460,
            'u-bd483': 2200,
            'u-bf647': 2619,
            'u-bf656': 2716,
            'u-bf703': 2760,
            #'u-bf705': 2815,
            'u-bh162': 2815,
        }

        cr_name = {  # Compare report name
            'u-bb075':
            'hist',  # UKESM1 Historical run (1960) with new SO2 emissions height
            'u-az524':
            'hist',  # UKESM1 Historical run (1995) with new SO2 emissions height
            'u-az513':
            'hist',  # UKESM1 Historical run (2020) with new SO2 emissions height
            'u-az515':
            'hist',  # UKESM1 Historical run (2050) with new SO2 emissions height
            'u-bb277':
            'hist',  # UKESM1 Historical run (2395) with new SO2 emissions height
            'u-bc179': 'hist',  # UKESM1 Historical run (2250)
            'u-bc292': 'hist',  # UKESM1 Historical run (2165)
            'u-bc370': 'hist',  # UKESM1 Historical run (2120)
            'u-bc470': 'hist',  # UKESM1 Historical run (2285)
            'u-bd288': 'hist',  # UKESM1 Historical run (2340)
            'u-bd416': 'hist',  # UKESM1 Historical run (2460)
            'u-bd483': 'hist',  # UKESM1 Historical run (2200)
            'u-bf647': 'hist',  # ;UKESM1 Historical run (2619)
            'u-bf656': 'hist',  # ;UKESM1 Historical run (2716)
            'u-bf703': 'hist',  # ;UKESM1 Historical run (2760)
            #'u-bf705': 'hist',	# ;UKESM1 Historical run (2815) - dead
            'u-bh162': 'hist',  # ;UKESM1 Historical run (2815)
            'u-bb446':
            '4xCO2',  # UKESM1 4xCO2 run (1960) with new SO2 emissions height
            'u-bb448':
            '1pcCO2',  # UKESM1 1%CO2 run (1960) with new SO2 emissions height
        }
        oldemssions = {}
        oldemssionscolours = {}

        for jobID, yr in list(start_year.items()):
            if jobID not in list(live_jobs.keys()): continue
            colourpair = {jobID: customColours[jobID], 'u-aw310': 'black'}
            lineThicknesses = {jobID: 1.7, 'u-aw310': 1.7}
            linestyles = {jobID: '-', 'u-aw310': '-'}
            if jobID in list(oldemssions.keys()):
                newrun = oldemssions[jobID]
                colourpair[newrun] = oldemssionscolours[newrun]
                lineThicknesses[newrun] = 2.4
                linestyles[newrun] = '--'
            timeseries_compare(
                colourpair,
                physics=1,
                bio=1,
                debug=0,
                year0='hist_vs_pi_' + str(yr),
                analysisname='UKESM1_pi_vs_' + cr_name[jobID] + '_' + str(yr),
                jobDescriptions=jobDescriptions,
                lineThicknesses=lineThicknesses,
                linestyles=linestyles,
            )
#hist_new_emissions = {j:customColours[j] for (j,h) in cr_name.items() if h == 'hist'}
        a = 0
        if a:
            timeseries_compare(
                {
                    j: customColours[j]
                    for (j, h) in list(cr_name.items()) if h == 'hist'
                },
                physics=1,
                bio=1,
                debug=0,
                year0='new_emissions',  #'from2240',#False, #'4800-5100',
                jobDescriptions=jobDescriptions,
                analysisname='UKESM_hist_new_emissions',
                lineThicknesses=hjthicknesses)

        a = False
        if a:
            jobs = [
                'u-aw310',
                'u-bb446',
            ]
            linestyles = defaultdict(lambda: '-')
            customColours = {
                'u-aw310': 'black',
                'u-bb446': 'green',
            }
            cnthicknesses = defaultdict(lambda: 1.1)
            timeseries_compare(
                customColours,  #{i:standards[i] for i in jobs},
                physics=1,
                bio=1,
                debug=0,
                year0='AlignControl',
                jobDescriptions=jobDescriptions,
                analysisname='UKESM1_4xco2_new_emssions',
                lineThicknesses=cnthicknesses,
                linestyles=linestyles,
            )
        a = False
        if a:
            jobs = [
                'u-aw310',
                'u-bb448',
            ]
            linestyles = defaultdict(lambda: '-')
            customColours = {
                'u-aw310': 'black',
                'u-bb448': 'red',
            }
            cnthicknesses = defaultdict(lambda: 1.1)
            timeseries_compare(
                customColours,  #{i:standards[i] for i in jobs},
                physics=1,
                bio=1,
                debug=0,
                year0='AlignControl',
                jobDescriptions=jobDescriptions,
                analysisname='UKESM1_1pcCO2_new_emssions',
                lineThicknesses=cnthicknesses,
                linestyles=linestyles,
            )

        a = 0
        if a:
            cncolours = {
                'u-av651': 'black',
                'u-aw310': 'green',
                'u-ay124': 'purple',
                'u-ay694': 'violet',
            }
            timeseries_compare(cncolours,
                               physics=1,
                               bio=1,
                               debug=0,
                               year0='UKESM_CN_control',
                               jobDescriptions=jobDescriptions,
                               analysisname='UKESM_CN_control',
                               lineThicknesses=hjthicknesses)


#		a = 1
#		if a:
#		        jobs = ['u-az513','u-az515','u-az524','u-az508', 'u-az021', 'u-az417', 'u-az418', 'u-bb075', 'u-bb277']
#		        timeseries_compare({
#		                 i:standards[i] for i in jobs},
#		                 physics=1,
#		                 bio=1,
#		                 debug=0,
#		                 year0='new_emissions', #'from2240',#False, #'4800-5100',
#		                 jobDescriptions=jobDescriptions,
#		                 analysisname='UKESM_hist_new_emissions',
#		                 lineThicknesses= hjthicknesses)

        a = 0
        if a:
            newemissionscolours = {
                'u-av651': 'black',
                'u-aw310': 'black',
                'u-az021': 'darkgreen',
                'u-az417': 'rebeccapurple',
                'u-az418': 'dodgerblue',
                'u-az513': 'darkorange',  # UKESM1 fifth Historical run (2020)
                'u-az515':
                'deeppink',  #      UKESM1 sixth Historical run (2050)
                'u-az524': 'goldenrod',  # UKESM1 seventh Historical run (1995)
                'u-bb075': 'green',
                'u-bb277': 'dodgerblue',
                'u-bc179':
                'pink',  #                  :standard ; UKESM1 Historical run (2250)
                'u-bc292':
                'brown',  #                :standard ; UKESM1 Historical run(2165)
                'u-bc370':
                'slateblue',  #                :standard ; UKESM1 Historical run(2120)
                'u-bc470':
                'gold',  #                :standard ; UKESM1 Historical run(2120)
                'u-bd288': 'orchid',  # UKESM1 Historical run (2340)
                'u-bd416': 'navy',  # UKESM1 Historical run (2460)
                'u-bd483': 'olive',  # UKESM1 Historical run (2200)
                'u-bf935': 'fushia',  # UKESM1 Historical run '2565', ],
                'u-bh100': 'yellow',  # UKESM1 Historical run '2685', ],
                'u-bh101': 'darkgold',  # UKESM1 Historical run '2745', ],
                'u-bf647': 'lightgreen',  # UKESM1 Historical run '2629', ],
                'u-bf656': 'cyan',  # UKESM1 Historical run '2716', ],
                'u-bf703': 'crimson',  # UKESM1 Historical run '2760', ],
                'u-bh162': 'orchid',  # UKESM1 Historical run '2815', ],}
            }
            linestyles = defaultdict(lambda: '-')
            linestyles['u-av651'] = '--'
            linestyles['u-aw310'] = ':'
            cnthicknesses = defaultdict(lambda: 1.1)
            cnthicknesses['u-aw310'] = 1.7
            timeseries_compare(
                newemissionscolours,
                physics=1,
                bio=1,
                debug=0,
                year0='ControlAligned',
                jobDescriptions=jobDescriptions,
                analysisname='UKESM1_newEmissions_ControlAligned',
                lineThicknesses=cnthicknesses,
                linestyles=linestyles,
            )

        a = False
        if a:
            ensembles = {}
            ensembles['PI Control 3'] = [
                'u-aw310',
            ]  #'u-av651',
            ensembles['PI Control 6'] = [
                'u-aw310',
            ]  #'u-av651',
            ensembles['PI Control 7'] = [
                'u-aw310',
            ]  #'u-av651',
            ensembles['New Emissions'] = [
                'u-az021', 'u-az417', 'u-az418', 'u-az513', 'u-az515',
                'u-az524'
            ]
            ensembles['Old Emissions'] = [
                'u-aw331', 'u-ax195', 'u-ax589', 'u-ax718', 'u-ay078',
                'u-ay167', 'u-ay491'
            ]
            jobs = []
            for ensemble in list(ensembles.keys()):
                jobs.extend(ensembles[ensemble])

            customColours = {i: standards[i] for i in jobs}
            customColours['PI Control 3'] = 'black'
            customColours['PI Control 6'] = 'black'
            customColours['PI Control 7'] = 'black'
            customColours['New Emissions'] = 'red'
            customColours['Old Emissions'] = 'blue'

            jobDescriptions[
                'PI Control 3'] = 'Pre industrial control 1850-1920 with same starting points as 3 new emissions runs'
            jobDescriptions[
                'PI Control 6'] = 'Pre industrial control 1920-present with same starting points as all 6 new emissions runs'
            jobDescriptions[
                'PI Control 7'] = 'Pre industrial control 1850-present with same starting points as all 7 initial historical runs.'
            jobDescriptions[
                'New Emissions'] = 'New SO2 emissions historial ensemble'
            jobDescriptions[
                'Old Emissions'] = 'Initial CMIP6 historical ensemble'

            linestyles = defaultdict(lambda: '-')
            linestyles['PI Control 3'] = ':'
            linestyles['PI Control 6'] = ':'
            linestyles['PI Control 7'] = '-'

            cnthicknesses = defaultdict(lambda: 1.1)
            #cnthicknesses['PI Control 3'] = 1.7
            #cnthicknesses['PI Control 6'] = 1.7

            timeseries_compare(customColours,
                               physics=1,
                               bio=1,
                               debug=1,
                               year0='EnsembleAlign',
                               jobDescriptions=jobDescriptions,
                               analysisname='UKESM1_CMIP6_ensembles_3',
                               lineThicknesses=cnthicknesses,
                               linestyles=linestyles,
                               ensembles=ensembles)

        a = False
        if a:
            jobs = [
                'u-aw310',
                'u-ar766',
                'u-av651',
                'u-aq853',
                'u-av450',
            ]
            linestyles = defaultdict(lambda: '-')
            linestyles['u-av651'] = '--'
            linestyles['u-aq853'] = '--'
            customColours = {
                'u-av651': 'black',
                'u-ar766': 'red',
                'u-aq853': 'red',
                'u-av450': 'blue',
                'u-aw310': 'black',
            }
            cnthicknesses = defaultdict(lambda: 1.1)
            cnthicknesses['u-aw310'] = 1.7
            timeseries_compare(
                customColours,  #{i:standards[i] for i in jobs},
                physics=1,
                bio=1,
                debug=0,
                year0='AlignToDECK',
                jobDescriptions=jobDescriptions,
                analysisname='UKESM1_piControl',
                lineThicknesses=cnthicknesses,
                linestyles=linestyles,
            )

        a = False
        if a:
            linestyles = defaultdict(lambda: '-')
            customColours = {
                'u-aw310': 'black',
                'u-aw448': 'blue',
                'u-ax202': 'magenta',
                'u-ax663': 'red',
                'u-ax725': 'purple',
                'u-bb448':
                'orange',  # UKESM1 1%CO2 run (1960) with new SO2 emissions height
            }
            linestyles['u-bb448'] = ':'
            cnthicknesses = defaultdict(lambda: 1.1)
            timeseries_compare(
                customColours,  #{i:standards[i] for i in jobs},
                physics=1,
                bio=1,
                debug=0,
                year0='AlignToDECK2050',
                jobDescriptions=jobDescriptions,
                analysisname='UKESM1_1pco2_old_and_new_emissions',
                lineThicknesses=cnthicknesses,
                linestyles=linestyles,
            )

            linestyles = defaultdict(lambda: '-')
            customColours = {
                'u-aw310': 'black',
                'u-aw447': 'green',
                'u-bb446':
                'orange',  # UKESM1 4xCO2 run (1960) with new SO2 emissions height
            }
            linestyles['u-bb446'] = ':'
            cnthicknesses = defaultdict(lambda: 1.1)
            timeseries_compare(
                customColours,  #{i:standards[i] for i in jobs},
                physics=1,
                bio=1,
                debug=0,
                year0='AlignToDECK2400',
                jobDescriptions=jobDescriptions,
                analysisname='UKESM1_4xco2_old_and_new_emissions',
                lineThicknesses=cnthicknesses,
                linestyles=linestyles,
            )

        print("Finished... ")
        return

        emissions_driven = True
        if emissions_driven:
            jobs = ['u-bh519', 'u-az513', 'u-bf703']
            customColours = {
                'u-bh519': 'blue',
                'u-bf703': 'brown',
                'u-az513': 'green',
            }
            cnthicknesses = defaultdict(lambda: 1.1)
            linestyles = defaultdict(lambda: '-')
            timeseries_compare(
                customColours,
                physics=1,
                bio=1,
                debug=0,
                year0='Historical',
                jobDescriptions=jobDescriptions,
                analysisname='Emissions_driven',
                lineThicknesses=cnthicknesses,
                linestyles=linestyles,
            )

        ocean_only = True
        if ocean_only:
            jobs = [
                'u-bc370', 'u-bl710', 'u-bl711', 'u-bl712', 'u-bl713',
                'u-bl714', 'u-bl715'
            ]
            customColours = {
                'u-bc370': 'black',
                'u-bl710': 'red',
                'u-bl711': 'blue',
                'u-bl712': 'purple',
                'u-bl713': 'green',
                'u-bl714': 'orange',
                'u-bl715': 'goldenrod',
            }
            descr = {
                'u-bc370': 'coupled',
                'u-bl710': 'T, S relaxation',
                'u-bl711': 'T, S relaxation (time-interp)',
                'u-bl712': 'T-relax only',
                'u-bl713': 'S-relax only',
                'u-bl714': 'no relaxation',
                'u-bl715': 'no relaxation (no FWB)',
            }
            cnthicknesses = defaultdict(lambda: 1.1)
            linestyles = defaultdict(lambda: '-')
            timeseries_compare(
                customColours,
                physics=1,
                bio=1,
                debug=0,
                year0=None,
                jobDescriptions=descr,
                analysisname='ocean_only',
                lineThicknesses=cnthicknesses,
                linestyles=linestyles,
            )

        ocean_only_TS = True
        if ocean_only_TS:
            jobs = [
                'u-bc370',
                'u-bl710',
                'u-bm176',
                'u-bm177',
                'u-bm178',
                'u-bm179',
            ]
            customColours = {
                'u-bc370': 'black',
                'u-bl710': 'red',
                'u-bm176': 'blue',
                'u-bm177': 'purple',
                'u-bm178': 'green',
                'u-bm179': 'orange',
            }
            descr = {
                'u-bc370': 'coupled',
                'u-bl710': 'T, S relaxation',
                'u-bm176': 'T, S relaxation; geographically uniform',
                'u-bm177': 'as bm176, but relaxation scaled x1/3',
                'u-bm178': 'as bm176, but relaxation scaled x3',
                'u-bm179': 'as bm176, but no S relaxation under ice',
            }

            cnthicknesses = defaultdict(lambda: 1.1)
            linestyles = defaultdict(lambda: '-')
            timeseries_compare(
                customColours,
                physics=1,
                bio=1,
                debug=0,
                year0=None,
                jobDescriptions=descr,
                analysisname='ocean_only_TS_relaxation',
                lineThicknesses=cnthicknesses,
                linestyles=linestyles,
            )

        ocean_only_S = True
        if ocean_only_S:
            jobs = [
                'u-bc370',
                'u-bl710',
                'u-bm180',
                'u-bm182',
                'u-bm183',
                'u-bm184',
            ]
            customColours = {
                'u-bc370': 'black',
                'u-bl710': 'red',
                'u-bm180': 'blue',
                'u-bm182': 'purple',
                'u-bm183': 'green',
                'u-bm184': 'orange',
            }
            descr = {
                'u-bc370': 'coupled',
                'u-bl710': 'T, S relaxation',
                'u-bm180': 'S relaxation; geographically uniform',
                'u-bm182': 'as bm180, but relaxation scaled x1/3',
                'u-bm183':
                'as bm180, but relaxation scaled x3 [this run ended fatally, but it still managed about a decade]',
                'u-bm184': 'as bm180, but no S relaxation under ice',
            }

            cnthicknesses = defaultdict(lambda: 1.1)
            linestyles = defaultdict(lambda: '-')
            timeseries_compare(
                customColours,
                physics=1,
                bio=1,
                debug=0,
                year0=None,
                jobDescriptions=descr,
                analysisname='ocean_only_S_relaxation',
                lineThicknesses=cnthicknesses,
                linestyles=linestyles,
            )

        low_res_atmosphere = True
        if low_res_atmosphere:
            jobs = ['u-bi700', 'u-aw310', 'u-bj002', 'u-bi831']  #'u-bi481',
            customColours = {
                #'u-bi481': 'blue',
                'u-bi700': 'red',
                'u-aw310': 'black',
                'u-bj002': 'orange',
                'u-bi831': 'purple',
            }
            cnthicknesses = defaultdict(lambda: 1.1)
            linestyles = defaultdict(lambda: '-')
            linestyles['u-aw310'] = '--'
            timeseries_compare(
                customColours,
                physics=1,
                bio=1,
                debug=0,
                year0='Historical',
                jobDescriptions=jobDescriptions,
                analysisname='N48_eORCA1',
                lineThicknesses=cnthicknesses,
                linestyles=linestyles,
            )

        low_res_atmosphere_timePlots = True
        if low_res_atmosphere_timePlots:
            jobs = [
                'u-bi700',
                'u-bi913',
                'u-bi914',
            ]
            customColours = {
                #'u-bi481': 'blue',
                'u-bi700': 'red',
                'u-bi913': 'orange',
                'u-bi914': 'purple',
            }
            cnthicknesses = defaultdict(lambda: 1.1)
            linestyles = defaultdict(lambda: '-')
            timeseries_compare(
                customColours,
                physics=1,
                bio=1,
                debug=0,
                year0='Historical',
                jobDescriptions=jobDescriptions,
                analysisname='N48_eORCA1_timesteps',
                lineThicknesses=cnthicknesses,
                linestyles=linestyles,
            )

        linestyles = defaultdict(lambda: '-')
        linestyles['u-av651'] = '--'
        linestyles['u-aq853'] = '--'
        customColours = {
            'u-av651': 'black',
            'u-as371': 'red',
            'u-aq853': 'red',
            'u-aw331': 'black',
            'u-ax195': 'green',
            'u-ax589': 'blue',
            'u-ax718': 'purple',
            'u-ay078': 'orange',
            'u-ay167': 'pink',
            'u-ay491': 'gold',
            'u-az942': 'darkgoldenrod',
        }
        cnthicknesses = defaultdict(lambda: 1.1)
        cnthicknesses['u-aw331'] = 1.7
        timeseries_compare(
            customColours,
            physics=1,
            bio=1,
            debug=0,
            year0='HistoricalDECK2050',
            jobDescriptions=jobDescriptions,
            analysisname='UKESM1_historical',
            lineThicknesses=cnthicknesses,
            linestyles=linestyles,
        )

        customColours = {
            'u-av651': 'black',
            'u-aw310': 'black',
            'u-aw331': 'teal',
            'u-ax195': 'green',
            'u-ax589': 'blue',
            'u-ax718': 'purple',
            'u-ay078': 'orange',
            'u-ay167': 'pink',
            'u-ay491': 'gold',
            'u-az942': 'darkgoldenrod',
        }
        linestyles = defaultdict(lambda: '-')
        linestyles['u-av651'] = '--'
        linestyles['u-aw310'] = ':'
        cnthicknesses = defaultdict(lambda: 1.1)
        cnthicknesses['u-aw310'] = 1.7
        timeseries_compare(
            customColours,
            physics=1,
            bio=1,
            debug=0,
            year0='ControlAligned',
            jobDescriptions=jobDescriptions,
            analysisname='UKESM1_historical_ControlAligned',
            lineThicknesses=cnthicknesses,
            linestyles=linestyles,
        )
        timeseries_compare(
            customColours,  #{i:standards[i] for i in jobs},
            physics=1,
            bio=1,
            debug=0,
            year0='AlignToDECK2050',
            jobDescriptions=jobDescriptions,
            analysisname='UKESM1_historical_pi',
            lineThicknesses=cnthicknesses,
            linestyles=linestyles,
        )

        #CompareTwoRuns('u-ba811','u-aw448',yearA='1980',yearB='1980')
        jobs = [
            'u-ba811',
            'u-aw448',
        ]  #-az524','u-az508', 'u-az021', 'u-az417', 'u-a8']
        timeseries_compare(
            {i: standards[i]
             for i in jobs},
            physics=1,
            bio=1,
            debug=0,
            year0=False,  #issions', #'from2240',#False, #'4800-5100',
            jobDescriptions=jobDescriptions,
            analysisname='UKESM_CO2_mmr_bug',
            lineThicknesses=hjthicknesses)

        #                jobs = ['u-bc016', 'u-bc057', 'u-bc058', 'u-bc060']
        #                timeseries_compare({
        #                        i:standards[i] for i in jobs},
        #                        physics=1,
        #                        bio=1,
        #                        debug=0,
        #                        year0=False,
        #                        jobDescriptions=jobDescriptions,
        #                        analysisname='CRESCENDO_OO_test_8',
        #                        lineThicknesses= hjthicknesses,
        #                        )

        #                jobs = ['u-au984','u-av079', 'u-ax628', 'u-ax629', ]
        #                timeseries_compare({
        #                        i:standards[i] for i in jobs},
        #                        physics=1,
        #                        bio=1,
        #                        debug=0,
        #                        year0=False,
        #                        jobDescriptions=jobDescriptions,
        #                        analysisname='CRESCENDO_OO_test_4',
        #                        lineThicknesses= hjthicknesses,
        #                        )

        #                jobs = ['u-av079', 'u-aw721', 'u-ay123', 'u-ay108']
        #                timeseries_compare({
        #                        i:standards[i] for i in jobs},
        #                        physics=1,
        #                        bio=1,
        #                        debug=0,
        #                        year0=False,
        #                        jobDescriptions=jobDescriptions,
        #                        analysisname='CRESCENDO_OO_test_5',
        #                        lineThicknesses= hjthicknesses,
        #                        )

        #                jobs = ['u-aw721', 'u-ax134']
        #		linestyles={}
        #                linestyles['u-aw721'] = '-'
        #                linestyles['u-ax134'] = '--'
        #                timeseries_compare({
        #                        i:standards[i] for i in jobs},
        #                        physics=1,
        #                        bio=1,
        #                        debug=0,
        #                        year0=False,
        #                        jobDescriptions=jobDescriptions,
        #                        analysisname='CRESCENDO_OO_test_6',
        #                        lineThicknesses= hjthicknesses,
        #			linestyles = linestyles,
        #                        )

        jobs = ['u-aw310', 'u-aw448', 'u-ax202', 'u-ax663', 'u-ax725']
        linestyles = defaultdict(lambda: '-')
        customColours = {
            'u-aw310': 'black',
            'u-aw448': 'blue',
            'u-ax202': 'magenta',
            'u-ax663': 'red',
            'u-ax725': 'purple'
        }
        cnthicknesses = defaultdict(lambda: 1.1)
        timeseries_compare(
            customColours,  #{i:standards[i] for i in jobs},
            physics=1,
            bio=1,
            debug=0,
            year0='AlignToDECK2050',
            jobDescriptions=jobDescriptions,
            analysisname='UKESM1_1pco2',
            lineThicknesses=cnthicknesses,
            linestyles=linestyles,
        )

        #                jobs = ['u-av937','u-aw310','u-aw072', 'u-aw331']
        #                timeseries_compare({
        #                         i:standards[i] for i in jobs},
        #                         physics=1,
        #                         bio=1,
        #                         debug=0,
        #                         year0=False,
        #                         jobDescriptions=jobDescriptions,
        #                         analysisname='UKESM1_DECK',
        #                         lineThicknesses= hjthicknesses)

        jobs = ['u-av651', 'u-aq853', 'u-as371', 'u-aw331', 'u-ax231']
        linestyles = defaultdict(lambda: '-')
        linestyles['u-av651'] = '--'
        linestyles['u-aq853'] = '--'
        customColours = {
            'u-av651': 'black',
            'u-aq853': 'red',
            'u-as371': 'red',
            'u-aw331': 'black',
            'u-ax231': 'green',
        }
        cnthicknesses = defaultdict(lambda: 1.75)
        cnthicknesses['u-av651'] = 1.1
        cnthicknesses['u-aq853'] = 1.1

        timeseries_compare(
            customColours,  #{i:standards[i] for i in jobs},
            physics=1,
            bio=1,
            debug=0,
            year0='AlignToDECK',
            jobDescriptions=jobDescriptions,
            analysisname='UKESM1_CN_DECK',
            lineThicknesses=cnthicknesses,
            linestyles=linestyles,
        )

        jobs = [
            'u-aw310',
            'u-aw447',
        ]
        linestyles = defaultdict(lambda: '-')
        customColours = {
            'u-aw310': 'black',
            'u-aw447': 'green',
        }
        cnthicknesses = defaultdict(lambda: 1.1)
        timeseries_compare(
            customColours,  #{i:standards[i] for i in jobs},
            physics=1,
            bio=1,
            debug=0,
            year0='AlignToDECK2200',
            jobDescriptions=jobDescriptions,
            analysisname='UKESM1_4xco2',
            lineThicknesses=cnthicknesses,
            linestyles=linestyles,
        )

        linestyles = defaultdict(lambda: '-')
        linestyles['u-aw310'] = '--'
        linestyles['u-ax941'] = '--'
        linestyles['u-av450'] = '--'
        customColours = {
            'u-aw310': 'black',
            'u-aw447': 'black',
            'u-ax672': 'blue',
            'u-av450': 'blue',
            'u-ax941': 'red',
            'u-ax945': 'red',
        }
        cnthicknesses = defaultdict(lambda: 1.1)
        timeseries_compare(
            customColours,  #{i:standards[i] for i in jobs},
            physics=1,
            bio=1,
            debug=0,
            year0='AlignToDECK2050',
            jobDescriptions=jobDescriptions,
            analysisname='UKESM1_4xco2_2',
            lineThicknesses=cnthicknesses,
            linestyles=linestyles,
        )

        jobs = ['u-am927i', 'u-aq853', 'u-ar783', 'u-au835', 'u-av450']
        colours = {i: standards[i] for i in jobs}
        thicknesses3 = defaultdict(lambda: 0.75)
        thicknesses3['u-ar766'] = 1.5
        thicknesses3['u-ar783'] = 2.
        thicknesses3['u-au835'] = 2.
        thicknesses3['u-at628'] = 2.
        thicknesses3['u-at760'] = 2.
        thicknesses3['u-at572'] = 2.
        thicknesses3['u-au027'] = 2.
        thicknesses3['u-au835'] = 2
        thicknesses3['u-au756'] = 2
        thicknesses3['u-au828'] = 2
        thicknesses3['u-av450'] = 2.6
        timeseries_compare(
            {i: standards[i]
             for i in jobs},
            physics=1,
            bio=1,
            debug=0,
            year0='ransom2',
            jobDescriptions=jobDescriptions,
            analysisname=
            'HCCP_C2.3',  # Called ransom because Colin requested this in exchange for help with my CMR.
            lineThicknesses=thicknesses3)

        jobs = [
            'u-aw310', 'u-ar766', 'u-av651', 'u-aq853', 'u-ar783', 'u-au835',
            'u-av472', 'u-aw700'
        ]
        linestyles = defaultdict(lambda: '-')
        linestyles['u-av651'] = '-.'
        linestyles['u-aq853'] = '--'
        linestyles['u-av472'] = '--'
        linestyles['u-ar783'] = '--'
        linestyles['u-au835'] = ':'

        customColours = {
            'u-av651': 'black',
            'u-ar766': 'red',
            'u-aq853': 'red',
            'u-aw310': 'black',
            'u-ar783': 'black',
            'u-au835': 'black',
            'u-av472': 'black',
            'u-aw700': 'deepskyblue',
        }
        cnthicknesses = defaultdict(lambda: 1.1)
        cnthicknesses['u-aw310'] = 1.7
        cnthicknesses['u-av472'] = 1.7

        timeseries_compare(
            customColours,  #{i:standards[i] for i in jobs},
            physics=1,
            bio=1,
            debug=0,
            year0='AlignToDECK1600-2250',
            jobDescriptions=jobDescriptions,
            analysisname='UKESM1_piControl_1600-1930',
            lineThicknesses=cnthicknesses,
            linestyles=linestyles,
        )

        print("Finished... ")
        return

        jobs = ['u-an869', 'u-ak900', 'u-ar538', 'u-ar977']
        colours = {i: standards[i] for i in jobs}
        timeseries_compare({i: standards[i]
                            for i in jobs},
                           physics=1,
                           bio=1,
                           debug=0,
                           year0='FullSpinUp',
                           jobDescriptions=jobDescriptions,
                           analysisname='OriginalOceanOnlySpinUp')

        jobs = [
            'u-an869', 'u-ak900', 'u-ar538', 'u-ar783', 'u-au835', 'u-av450',
            'u-aj588'
        ]
        colours = {i: standards[i] for i in jobs}
        timeseries_compare({i: standards[i]
                            for i in jobs},
                           physics=1,
                           bio=1,
                           debug=0,
                           year0='UKESMv1SpinUp',
                           jobDescriptions=jobDescriptions,
                           analysisname='UKESMv1_fullSpinUp')

        print("Finished... ")
        return

if __name__ == "__main__": 
    main()
