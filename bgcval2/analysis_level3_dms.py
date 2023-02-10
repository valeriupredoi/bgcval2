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
.. module:: analysis_level3_dms
   :platform: Unix
   :synopsis: A script to produce a level 3 analysis for DMS.

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

#####
# Load specific local code:
from bgcval2.bgcvaltools import bv2tools as bvt
from .timeseries import timeseriesAnalysis
from .timeseries import profileAnalysis
from .timeseries import timeseriesPlots as tsp

#####
# User defined set of paths pointing towards the datasets.
from .Paths import paths as paths

medusaCoords = {
    't': 'time_counter',
    'z': 'deptht',
    'lat': 'nav_lat',
    'lon': 'nav_lon',
    'cal': '360_day',
}  # model doesn't need time dict.
dmsCoords = {
    't': 'time',
    'z': 'depth',
    'lat': 'Latitude',
    'lon': 'Longitude',
    'cal': 'standard',
    'tdict': bvt.tdicts['ZeroToZero']
}


def analysis_dms(jobID=''):
    annual = True

    analysisDict = {}
    imagedir = bvt.folder(paths.imagedir + '/' + jobID + '/Level3/DMS')
    shelvedir = bvt.folder(paths.shelvedir + '/' + jobID + '/Level3/DMS')

    regionList = [
        'Global',
        'ignoreInlandSeas',
        'SouthernOcean',
        'ArcticOcean',
        'Equator10',
        'Remainder',
        'NorthernSubpolarAtlantic',
        'NorthernSubpolarPacific',
    ]
    metricList = [
        'mean',
    ]
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

    #####
    # A time series analysis for the DMS fields.

    for name in ['DMS_ARAN', 'DMS_ANDR', 'DMS_SIMO', 'DMS_HALL']:

        dmsfiles = listModelDataFiles(jobID, 'diad_T', paths.ModelFolder_pref,
                                      annual)
        if name == 'DMS_ARAN':
            analysisDict['modelFiles'] = dmsfiles
        else:
            analysisDict['modelFiles'] = bvt.listFiles(dmsfiles,
                                                       want=100,
                                                       listType='backloaded',
                                                       first=30,
                                                       last=10)

        if annual:
            analysisDict['dataFile'] = paths.DMSDir + 'DMSclim_mean.nc'
        else:
            analysisDict['dataFile'] = ''

        analysisDict['modelcoords'] = medusaCoords
        analysisDict['datacoords'] = dmsCoords

        analysisDict['modeldetails'] = {
            'name': name,
            'vars': [
                'DMS_ARAN',
            ],
            'convert': bvt.mul1000000,
            'units': 'umol/m3'
        }
        analysisDict['datadetails'] = {
            'name': name,
            'vars': [
                'DMS',
            ],
            'convert': bvt.NoChange,
            'units': 'umol/m3'
        }

        analysisDict['layers'] = [
            'layerless',
        ]
        analysisDict['regions'] = regionList
        analysisDict['metrics'] = metricList

        analysisDict['datasource'] = 'Lana'
        analysisDict['model'] = 'MEDUSA'

        analysisDict['modelgrid'] = 'eORCA1'
        analysisDict['gridFile'] = paths.orcaGridfn
        analysisDict['Dimensions'] = 2

        tsa = timeseriesAnalysis(
            analysisDict['modelFiles'],
            analysisDict['dataFile'],
            dataType=name,
            modelcoords=analysisDict['modelcoords'],
            modeldetails=analysisDict['modeldetails'],
            datacoords=analysisDict['datacoords'],
            datadetails=analysisDict['datadetails'],
            datasource=analysisDict['datasource'],
            model=analysisDict['model'],
            jobID=jobID,
            layers=analysisDict['layers'],
            regions=analysisDict['regions'],
            metrics=analysisDict['metrics'],
            workingDir=shelvedir,
            imageDir=imagedir,
            grid=analysisDict['modelgrid'],
            gridFile=analysisDict['gridFile'],
            clean=False,
        )

        dataD[name] = tsa.dataD
        modeldataD[name] = tsa.modeldataD

    #####
    # Prepare a time series comparison of the four DMS types.

    timesD = {}
    arrD = {}

    region = 'Global'
    for name in list(dataD.keys()):
        try:
            mdata = modeldataD[name][(region, 'layerless', 'mean')]
        except:
            continue

        timesD[name] = sorted(mdata.keys())
        arrD[name] = [mdata[t] for t in timesD[name]]

    colours = {
        'DMS_ARAN': 'red',
        'DMS_SIMO': 'orange',
        'DMS_ANDR': 'blue',
        'DMS_HALL': 'purple',
    }
    title = 'DMS ' + region

    for ls in ['Both', 'smoothed', 'movingaverage']:
        tsp.multitimeseries(
            timesD,  # model times (in floats)
            arrD,  # model time series
            data=-999,  # in situ data distribution
            title=title,
            filename=bvt.folder(imagedir) + 'DMS_' + region + '_' + ls +
            '.png',
            units='',
            plotStyle='Together',
            lineStyle=ls,
            colours=colours,
        )


def main():
    if "--help" in argv or len(argv) == 1:
        print("Running with no arguments. Exiting.")
        if "--help" in argv:
            print("Read the documentation.")
        sys.exit(0)
    try:
        jobID = argv[1]
    except:
        jobID = "u-ab749"

    if 'debug' in argv[1:]: suite = 'debug'

    analysis_dms(jobID=jobID, )  #clean=1)
    #if suite == 'all':
    #analysis_timeseries(jobID =jobID,analysisSuite='FullDepth', z_component = 'FullDepth',)#clean=1)


if __name__ == "__main__":
    main()
