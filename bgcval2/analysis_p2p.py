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
.. module:: analysis_p2p
   :platform: Unix
   :synopsis: A script to produce point to point analysis of model vs data.
.. moduleauthor:: Lee de Mora <ledm@pml.ac.uk>

"""

#####
#Standard Python modules:
from os.path import exists
from calendar import month_name
from socket import gethostname
from getpass import getuser
from glob import glob
from netCDF4 import Dataset
import numpy as np
import sys
import os
import argparse

#####
#Specific local code:
from . import UKESMpython as ukp
from .p2p import makePatternStatsPlots, testsuite_p2p
from .p2p.summaryTargets import summaryTargets
from .p2p.patternAnalyses import InterAnnualPatterns, BGCvsPhysics
from .bgcvaltools.pftnames import months
from .p2p.shelveToDictionary import shelveToDictionary

#####
# User defined set of paths pointing towards the datasets.
from .Paths import paths as paths

#####
# code plan:
#	This is a the script that calls testsuite_p2p now.
#	Now all code is run though that testsuite.
#	the idea being that each analysis produces a new one of these analysis tools.
#

# 	from

#####
# suites:
p2pKeys = [
    'T',
    'S',  #'MLD',
    'Chl_pig',
    'Chl_CCI',
    'N',
    'Si',
    'O2',
    'Alk',
    'DIC',
    'AirSeaFlux',
    'IntPP_OSU',
    'Diatoms',
    'Microzoo',
    'Mesozoo',
]

p2pKeys_annual = [
    'T',
    'S',  #'MLD',
    'Chl_CCI',
    'N',
    'Si',
    'O2',
    'Alk',
    'DIC',
    'AirSeaFlux',
    'IntPP_OSU',
]

p2pKeys_level2 = [
    'Chl_CCI',
    'N',
    'Si',
    'O2',
    'Alk',
    'DIC',
    'AirSeaFlux',
    'IntPP_OSU',
    'Dust',
    'T',
    'S',  #'MLD',
    #'ZonalCurrent','MeridionalCurrent','VerticalCurrent'
]

p2pKeys_physics = [
    'T',
    'S',  #'MLD',
    'ZonalCurrent',
    'MeridionalCurrent',
    'VerticalCurrent'
]

p2pDict = {i: n for i, n in enumerate(p2pKeys)}
p2pDict_annual = {i: n for i, n in enumerate(p2pKeys_annual)}
p2pDict_level2 = {i: n for i, n in enumerate(p2pKeys_level2)}
p2pDict_physics = {i: n for i, n in enumerate(p2pKeys_physics)}


def analysis_p2p(
    models=[
        'NEMO',
        'MEDUSA',
    ],
    jobID='u-ad980',
    years=['1077'],  #'2075','2076',
    modelGrid='eORCA1',
    annual=True,
    noPlots=False,
    analysisSuite='default',
    noTargets=True,
):
    """

	"""
    #####
    # Switches:
    # These are some booleans that allow us to choose which analysis to run.
    # This lets up give a list of keys one at a time, or in parrallel.
    if type(analysisSuite) == type(['Its', 'A', 'list!']):
        analysisKeys = analysisSuite
        print("analysisSuite is a list", analysisSuite, analysisKeys)
    #####
    # Switches:
    # These are some preset switches to run in series.
    if type(analysisSuite) == type('Its_A_string'):

        analysisKeys = []
        if analysisSuite.lower() in [
                'all',
                'default',
        ]:
            analysisKeys.append('Chl_CCI')  # CCI Chlorophyll
            analysisKeys.append(
                'Chl_pig')  # Chlorophyll from pigments (MAREDAT)
            analysisKeys.append(
                'Diatoms')  # Chlorophyll from pigments (MAREDAT)
            analysisKeys.append(
                'Microzoo')  # Chlorophyll from pigments (MAREDAT)
            analysisKeys.append(
                'Mesozoo')  # Chlorophyll from pigments (MAREDAT)
            analysisKeys.append('N')  # WOA Nitrate
            analysisKeys.append('Si')  # WOA Siliate
            analysisKeys.append('O2')  # WOA Oxygen
            analysisKeys.append('Fe')  # Iron
            analysisKeys.append('Alk')  # Glodap Alkalinity
            analysisKeys.append('DIC')  # Globap tCO2
            analysisKeys.append('AirSeaFlux')  # work in progress
            analysisKeys.append(
                'IntPP_OSU')  # OSU Integrated primpary production
            #####
            # Physics switches:
            analysisKeys.append('T')  # WOA Temperature
            analysisKeys.append('S')  # WOA Salinity
            analysisKeys.append(
                'MLD')  # iFERMER Mixed Layer Depth - work in prgress

        if analysisSuite.lower() in [
                'annual',
        ]:
            analysisKeys.extend(p2pKeys_annual)

        if analysisSuite.lower() in [
                'level2',
        ]:
            analysisKeys.extend(p2pKeys_level2)

        if analysisSuite.lower() in [
                'physics',
        ]:
            analysisKeys.extend(p2pKeys_physics)

        if analysisSuite.lower() in [
                'debug',
        ]:
            #analysisKeys.append('ZonalCurrent')             # Zonal Veloctity
            #analysisKeys.append('MeridionalCurrent')        # Meridional Veloctity
            #analysisKeys.append('VerticalCurrent')          # Vertical Veloctity
            #analysisKeys.append('Dust')          # Vertical Veloctity
            #analysisKeys.append('AirSeaFlux')          # work in progress
            analysisKeys.append('Chl_CCI')  # CCI Chlorophyll

        print("analysisSuite is a string", analysisSuite, analysisKeys)

    #####
    # Location of data files.
    if gethostname().find('pmpc') > -1:
        print("analysis-p2p.py:\tBeing run at PML on ", gethostname())

        if annual:
            WOAFolder = paths.WOAFolder_annual
            ModelFolder_pref = paths.ModelFolder_pref
        else:
            WOAFolder = paths.WOAFolder
            ModelFolder_pref = paths.ModelFolder_pref + "/" + jobID + "_postProc/"

    #	if annual:
    #		#####
    # No need to stitch together multiple months into one file:
    #WOAFolder 		= "/data/euryale7/scratch/ledm/WOA/annual/"
    #ModelFolder_pref	= "/data/euryale7/scratch/ledm/UKESM/MEDUSA/"
    #	else:
    #WOAFolder 	= "/data/euryale7/scratch/ledm/WOA/"
    #ModelFolder_pref	= "/data/euryale7/scratch/ledm/UKESM/MEDUSA/"+jobID+"_postProc/"
    #MAREDATFolder 	= "/data/euryale7/scratch/ledm/MAREDAT/MAREDAT/"
    #GEOTRACESFolder = "/data/euryale7/scratch/ledm/GEOTRACES/GEOTRACES_PostProccessed/"
    #TakahashiFolder = "/data/euryale7/scratch/ledm/Takahashi2009_pCO2/"
    #MLDFolder	= "/data/euryale7/scratch/ledm/IFREMER-MLD/"

    # New eORCA1 grid
    #orcaGridfn 	= '/data/euryale7/scratch/ledm/UKESM/MEDUSA/mesh_mask_eORCA1_wrk.nc'

        imgDir = paths.imagedir

        if annual: WOAFolder = paths.WOAFolder_annual
        else: WOAFolder = paths.WOAFolder

    if gethostname().find('ceda.ac.uk') > -1:
        print("analysis-p2p.py:\tBeing run at CEDA on ", gethostname())

        ObsFolder = paths.ObsFolder  #"/group_workspaces/jasmin/esmeval/example_data/bgc/"
        modelFolder = paths.ModelFolder_pref  #"/group_workspaces/jasmin2/ukesm/BGC_data/"
        #####
        # Location of model files.
        ModelFolder_pref = ukp.folder(modelFolder)

        #####
        # Location of data files.
        if annual: WOAFolder = ukp.folder(ObsFolder + "WOA/annual")
        else: WOAFolder = ukp.folder(ObsFolder + "WOA/")

        MAREDATFolder = ObsFolder + "/MAREDAT/MAREDAT/"
        GEOTRACESFolder = ObsFolder + "/GEOTRACES/GEOTRACES_PostProccessed/"
        TakahashiFolder = ObsFolder + "/Takahashi2009_pCO2/"
        MLDFolder = ObsFolder + "/IFREMER-MLD/"
        iMarNetFolder = ObsFolder + "/LestersReportData/"
        GlodapDir = ObsFolder + "/GLODAP/"
        GLODAPv2Dir = ObsFolder + "/GLODAPv2/GLODAPv2_Mapped_Climatologies/"
        OSUDir = ObsFolder + "OSU/"
        CCIDir = ObsFolder + "CCI/"

        # Directory for output files:
        imgDir = paths.imagedir

        # eORCA1 grid
        #orcaGridfn 	= '/group_workspaces/jasmin/esmeval/example_data/bgc/mesh_mask_eORCA1_wrk.nc'

        # MONSOON
    if gethostname().find('monsoon') > -1:
        print("Please set up paths.py")
        assert 0
#                print "analysis-timeseries.py:\tBeing run at the Met Office on ",gethostname()
#               machinelocation = 'MONSOON'
#
#               ObsFolder       = "/projects/ukesm/ldmora/BGC-data/"
#              ModelFolder       = "/projects/ukesm/ldmora/UKESM"
#             #####
#            # Location of model files.
#           ModelFolder_pref       = ukp.folder(ModelFolder)
#
#               #####
#              # Location of data files.
#                if annual:      WOAFolder       = ukp.folder(ObsFolder+"WOA/annual")
#             else:           WOAFolder       = ukp.folder(ObsFolder+"WOA/")
#
#               MAREDATFolder   = ObsFolder+"/MAREDAT/MAREDAT/"
#              GEOTRACESFolder = ObsFolder+"/GEOTRACES/GEOTRACES_PostProccessed/"
#             TakahashiFolder = ObsFolder+"/Takahashi2009_pCO2/"
#            MLDFolder       = ObsFolder+"/IFREMER-MLD/"
#           iMarNetFolder   = ObsFolder+"/LestersReportData/"
#          GlodapDir       = ObsFolder+"/GLODAP/"
#         GLODAPv2Dir     = ObsFolder+"/GLODAPv2/GLODAPv2_Mapped_Climatologies/"
#        OSUDir          = ObsFolder+"OSU/"
#       CCIDir          = ObsFolder+"CCI/"
#      if jobID in ["xkrus",]:
#             # Old school ORCA1 grid
#            orcaGridfn      =ModelFolder+'/mesh_mask_ORCA1_75.nc'
#   else:
#          # New eORCA1 grid
#         orcaGridfn      = ModelFolder+'/mesh_mask_eORCA1_wrk.nc'
# paths.p2p_ppDir         = "/projects/ukesm/"+getuser()+"/UKESM_postprocessed"
#imgDir          = ukp.folder('images')

    def listModelDataFiles(jobID, filekey, datafolder, annual, yr):
        print("listing model data files:", jobID, filekey, datafolder, annual)
        if annual:
            keystr = datafolder + jobID + "/" + jobID + "o_1y_*1201[-_]" + yr + '????_' + filekey + ".nc"
            print("listModelDataFiles:", keystr)
            return sorted(glob(keystr))[0]
        else:
            return sorted(
                glob(datafolder + jobID + "/" + jobID + "o_1m_*" + yr +
                     "????_" + filekey + ".nc"))[-1]

#####
# Because we can never be sure someone won't randomly rename the
# time dimension without saying anything.
# if jobID in ['u-am515','u-am927','u-am064','u-an326',]:

    print(jobID, 'grid_T', paths.ModelFolder_pref, annual)
    tmpModelFiles = listModelDataFiles(jobID, 'grid_T', paths.ModelFolder_pref,
                                       annual, '*')

    try:
        tmpModelFiles = listModelDataFiles(jobID, 'grid_T',
                                           paths.ModelFolder_pref, annual, '*')
    except:
        print(
            "No grid_T Model files available to figure out what naming convention is used."
        )
        tmpModelFiles = []
    ukesmkeys = {}
    if len(tmpModelFiles):
        print('test opening:', tmpModelFiles)
        nctmp = Dataset(tmpModelFiles, 'r')
        nctmpkeys = list(nctmp.variables.keys())
        nctmp.close()
        if 'votemper' in nctmpkeys:
            ukesmkeys = {}
            ukesmkeys['time'] = 'time_counter'
            ukesmkeys['temp3d'] = 'votemper'
            ukesmkeys['sst'] = ''
            ukesmkeys['sal3d'] = 'vosaline'
            ukesmkeys['sss'] = ''
            ukesmkeys['v3d'] = 'vomecrty'
            ukesmkeys['u3d'] = 'vozocrtx'
            ukesmkeys['e3u'] = 'e3u'
            ukesmkeys['w3d'] = 'vovecrtz'
        else:
            ukesmkeys['time'] = 'time_centered'
            ukesmkeys['temp3d'] = 'thetao'
            ukesmkeys['sst'] = 'tos'
            ukesmkeys['sal3d'] = 'so'
            ukesmkeys['sss'] = 'sos'
            ukesmkeys['v3d'] = 'vo'
            ukesmkeys['u3d'] = 'uo'
            ukesmkeys['e3u'] = 'thkcello'
            ukesmkeys['w3d'] = 'wo'
    else:
        print("No grid_T files Found")
    print('ukesmkeys[sal3d]:', ukesmkeys['sal3d'])

    #####
    # Set which spatial and temporal limitations to plot.
    transects = [
        'AtlanticTransect',
        'PacificTransect',
    ]
    justAll = [
        'All',
    ]  # All is not a slice, it has no cut on location, time, or depth.
    AllStandard = ['All', 'Standard', 'ignoreInlandSeas']
    HighLatWinter = [
        'All',
        'HighLatWinter',
    ]
    tsRegions = [
        'Global', 'Equator10', 'Remainder', 'ArcticOcean',
        'NorthernSubpolarAtlantic', 'NorthernSubpolarPacific',
        'ignoreInlandSeas', 'SouthernOcean', 'AtlanticSOcean'
    ]

    depthLevels = [
        'Surface',
        '500m',
        '1000m',
        'Transect',
        'PTransect',
        'SOTransect',
        'ArcTransect',
        'AntTransect',
        'CanRusTransect',
    ]

    medusaCoords = {
        't': 'index_t',
        'z': 'deptht',
        'lat': 'nav_lat',
        'lon': 'nav_lon',
        'cal': '360_day',
        'tdict': ukp.tdicts['ZeroToZero']
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

    maredatCoords = {
        't': 'index_t',
        'z': 'DEPTH',
        'lat': 'LATITUDE',
        'lon': 'LONGITUDE',
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
    cciCoords = {
        't': 'index_t',
        'z': 'index_z',
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
        'tdict': ukp.tdicts['ZeroToZero']
    }
    osuCoords = {
        't': 'index_t',
        'z': 'index_z',
        'lat': 'latitude',
        'lon': 'longitude',
        'cal': 'standard',
        'tdict': ukp.tdicts['ZeroToZero']
    }
    glodapv2Coords = {
        't': 'index_t',
        'z': 'Pressure',
        'lat': 'lat',
        'lon': 'lon',
        'cal': '',
        'tdict': {
            0: 0,
        }
    }
    takahashiCoords = {
        't': 'index_t',
        'z': 'index_z',
        'lat': 'LAT',
        'lon': 'LON',
        'cal': 'standard',
        'tdict': ukp.tdicts['ZeroToZero']
    }
    godasCoords = {
        't': 'index_t',
        'z': 'level',
        'lat': 'lat',
        'lon': 'lon',
        'cal': 'standard',
        'tdict': ukp.tdicts['ZeroToZero']
    }

    shelvesAV = []
    for year in years:
        #####
        # Location of model files.
        if annual:
            ModelFolder = ModelFolder_pref + jobID + "/"
        else:
            ModelFolder = ModelFolder_pref + year + '/'

        #####
        # AutoVivification is a form of nested dictionary.
        # We use AutoVivification here to determine which files to analyse and which fields in those files.
        # depthLevel is added, because some WOA files are huges and my desktop can not run the p2p analysis of that data.
        av = ukp.AutoVivification()
        #if 'Chl_pig' in analysisKeys:
        #	name = 'Chlorophyll_pig'
        #	av[name]['Data']['File'] 		= paths.MAREDATFolder+"MarEDat20121001Pigments.nc"
        #	if modelGrid == 'ORCA1':		av[name]['MEDUSA']['File'] 	= ModelFolder+jobID+'_' + year+"_CHL.nc"
        #	if modelGrid == 'ORCA025':		av[name]['MEDUSA']['File']	= ModelFolder+"xjwki_1979_CH.nc"
        #
        #	av[name]['Data']['coords'] 		= maredatCoords
        #	av[name]['MEDUSA']['coords']		= medusaCoords
        #
        #	av[name]['MEDUSA']['details']		= {'name': 'CHL', 'vars':['CHL',], 'convert': ukp.NoChange,'units':'mg C/m^3'}
        #	av[name]['Data']['details']		= {'name': 'Chlorophylla', 'vars':['Chlorophylla',], 'convert': ukp.div1000,'units':'ug/L'}
        #	av[name]['Data']['source'] 		= 'MAREDAT'
        #	av[name]['MEDUSA']['source']		= 'MEDUSA'
        #	av[name]['depthLevels'] 		= ['',]
        #	av[name]['MEDUSA']['grid']		= modelGrid
        #	av[name]['plottingSlices'] 		= tsRegions

        if 'Chl_CCI' in analysisKeys:
            name = 'Chlorophyll_cci'
            if annual:
                av[name]['Data'][
                    'File'] = paths.CCIDir + "ESACCI-OC-L3S-OC_PRODUCTS-CLIMATOLOGY-16Y_MONTHLY_1degree_GEO_PML_OC4v6_QAA-annual-fv2.0.nc"
                av[name]['MEDUSA']['File'] = listModelDataFiles(
                    jobID, 'ptrc_T', paths.ModelFolder_pref, annual, year)
#			else:
#				av[name]['Data']['File'] 	= paths.CCIDir+'ESACCI-OC-L3S-OC_PRODUCTS-CLIMATOLOGY-16Y_MONTHLY_1degree_GEO_PML_OC4v6_QAA-all-fv2.0.nc'
#				av[name]['MEDUSA']['File'] 	= ModelFolder+jobID+'_' + year+"_CHL.nc"

            av[name]['MEDUSA']['grid'] = modelGrid
            av[name]['depthLevels'] = [
                '',
            ]
            if annual: av[name]['plottingSlices'] = tsRegions
            else: av[name]['plottingSlices'] = HighLatWinter

            av[name]['Data']['coords'] = cciCoords
            av[name]['MEDUSA']['coords'] = medusaCoords

            av[name]['Data']['source'] = 'CCI'
            av[name]['MEDUSA']['source'] = 'MEDUSA'

            av[name]['MEDUSA']['details'] = {
                'name': name,
                'vars': ['CHN', 'CHD'],
                'convert': ukp.sums,
                'units': 'mg C/m^3'
            }
            av[name]['Data']['details'] = {
                'name': name,
                'vars': [
                    'chlor_a',
                ],
                'convert': ukp.NoChange,
                'units': 'mg C/m^3'
            }

        if 'Diatoms' in analysisKeys:
            name = 'Diatoms'
            if annual:
                print("No diatoms iron file", end=' ')
                assert 0
            av[name]['Data'][
                'File'] = paths.MAREDATFolder + "MarEDat20120716Diatoms.nc"
            av[name]['MEDUSA'][
                'File'] = ModelFolder + jobID + '_' + year + "_PHD.nc"

            av[name]['depthLevels'] = [
                '',
            ]
            av[name]['MEDUSA']['grid'] = modelGrid
            av[name]['plottingSlices'] = AllStandard

            av[name]['Data']['coords'] = maredatCoords
            av[name]['MEDUSA']['coords'] = medusaCoords

            av[name]['MEDUSA']['details'] = {
                'name': name,
                'vars': [
                    'PHD',
                ],
                'convert': ukp.N2Biomass,
                'units': 'mg C/m^3'
            }
            av[name]['Data']['details'] = {
                'name': name,
                'vars': [
                    'BIOMASS',
                ],
                'convert': ukp.NoChange,
                'units': 'mg C/m^3'
            }

            av[name]['Data']['source'] = 'MAREDAT'
            av[name]['MEDUSA']['source'] = 'MEDUSA'

        if 'Microzoo' in analysisKeys:
            name = 'Microzoo'
            if annual:
                print("No microzoo iron file", end=' ')
                assert 0
            av[name]['Data'][
                'File'] = paths.MAREDATFolder + "MarEDat20120424Microzooplankton.nc"
            av[name]['MEDUSA'][
                'File'] = ModelFolder + jobID + '_' + year + "_ZMI.nc"

            av[name]['MEDUSA']['grid'] = modelGrid
            av[name]['depthLevels'] = [
                '',
            ]
            av[name]['plottingSlices'] = AllStandard

            av[name]['Data']['coords'] = maredatCoords
            av[name]['MEDUSA']['coords'] = medusaCoords

            av[name]['MEDUSA']['details'] = {
                'name': name,
                'vars': [
                    'ZMI',
                ],
                'convert': ukp.N2Biomass,
                'units': 'mg C/m^3'
            }
            av[name]['Data']['details'] = {
                'name': name,
                'vars': [
                    'BIOMASS',
                ],
                'convert': ukp.NoChange,
                'units': 'mg C/m^3'
            }

            av[name]['Data']['source'] = 'MAREDAT'
            av[name]['MEDUSA']['source'] = 'MEDUSA'

        if 'Mesozoo' in analysisKeys:
            name = 'Mesozoo'
            if annual:
                print("No mesozoo iron file", end=' ')
                assert 0

            av[name]['Data'][
                'File'] = paths.MAREDATFolder + "MarEDat20120705Mesozooplankton.nc"
            av[name]['MEDUSA'][
                'File'] = ModelFolder + jobID + '_' + year + "_ZME.nc"

            av[name]['MEDUSA']['grid'] = modelGrid
            av[name]['depthLevels'] = [
                '',
            ]
            av[name]['plottingSlices'] = AllStandard

            av[name]['Data']['coords'] = maredatCoords
            av[name]['MEDUSA']['coords'] = medusaCoords

            av[name]['MEDUSA']['details'] = {
                'name': name,
                'vars': [
                    'ZME',
                ],
                'convert': ukp.N2Biomass,
                'units': 'mg C/m^3'
            }
            av[name]['Data']['details'] = {
                'name': name,
                'vars': [
                    'BIOMASS',
                ],
                'convert': ukp.NoChange,
                'units': 'mg C/m^3'
            }

            av[name]['Data']['source'] = 'MAREDAT'
            av[name]['MEDUSA']['source'] = 'MEDUSA'

        if 'N' in analysisKeys:
            name = 'Nitrate'
            if annual:
                av[name]['Data']['File'] = WOAFolder + 'woa13_all_n00_01.nc'
                av[name]['MEDUSA']['File'] = listModelDataFiles(
                    jobID, 'ptrc_T', paths.ModelFolder_pref, annual, year)
            else:
                av[name]['Data'][
                    'File'] = WOAFolder + 'nitrate_monthly_1deg.nc'
                if modelGrid == 'ORCA1':
                    av[name]['MEDUSA'][
                        'File'] = ModelFolder + jobID + '_' + year + "_DIN.nc"
                if modelGrid == 'ORCA025':
                    av[name]['MEDUSA'][
                        'File'] = ModelFolder + jobID + '_' + year + "_DIN.nc"

            av[name]['MEDUSA']['grid'] = modelGrid
            av[name]['depthLevels'] = depthLevels
            if annual: av[name]['plottingSlices'] = tsRegions
            else: av[name]['plottingSlices'] = HighLatWinter

            av[name]['Data']['coords'] = woaCoords
            av[name]['MEDUSA']['coords'] = medusaCoords

            av[name]['Data']['source'] = 'WOA'
            av[name]['MEDUSA']['source'] = 'MEDUSA'

            av[name]['MEDUSA']['details'] = {
                'name': name,
                'vars': [
                    'DIN',
                ],
                'convert': ukp.NoChange,
            }
            av[name]['Data']['details'] = {
                'name': name,
                'vars': [
                    'n_an',
                ],
                'convert': ukp.NoChange,
            }  # no units?

        if 'Si' in analysisKeys:
            name = 'Silicate'
            if annual:
                av[name]['Data']['File'] = WOAFolder + 'woa13_all_i00_01.nc'
                av[name]['MEDUSA']['File'] = listModelDataFiles(
                    jobID, 'ptrc_T', paths.ModelFolder_pref, annual, year)
            else:
                av[name]['Data'][
                    'File'] = WOAFolder + 'silicate_monthly_1deg.nc'
                av[name]['MEDUSA'][
                    'File'] = ModelFolder + jobID + '_' + year + "_SIL.nc"

            av[name]['MEDUSA']['grid'] = modelGrid
            av[name]['depthLevels'] = depthLevels
            if annual: av[name]['plottingSlices'] = tsRegions
            else: av[name]['plottingSlices'] = HighLatWinter

            av[name]['Data']['coords'] = woaCoords
            av[name]['MEDUSA']['coords'] = medusaCoords

            av[name]['Data']['source'] = 'WOA'
            av[name]['MEDUSA']['source'] = 'MEDUSA'

            av[name]['MEDUSA']['details'] = {
                'name': name,
                'vars': [
                    'SIL',
                ],
                'convert': ukp.NoChange,
            }
            av[name]['Data']['details'] = {
                'name': name,
                'vars': [
                    'i_an',
                ],
                'convert': ukp.NoChange,
            }  # no units?

        if 'Fe' in analysisKeys:
            name = 'Iron'
            if annual:
                print("No annual iron file", end=' ')
                assert 0

            av[name]['Data'][
                'File'] = paths.GEOTRACESFolder + "Iron_GEOTRACES_IDP2014_Discrete_Sample_Data_ascii.nc"
            av[name]['MEDUSA']['File'] = listModelDataFiles(
                jobID, 'ptrc_T', paths.ModelFolder_pref, annual, year)

            av[name]['depthLevels'] = [
                '',
            ]
            av[name]['MEDUSA']['grid'] = modelGrid
            av[name]['plottingSlices'] = justAll

            av[name]['Data']['coords'] = {
                't': 'MONTH',
                'z': 'DEPTH',
                'lat': 'Latitude',
                'lon': 'Longitude',
                'cal': 'standard',
                'tdict': ukp.tdicts['OneToZero']
            }
            av[name]['MEDUSA']['coords'] = medusaCoords

            av[name]['Data']['source'] = 'GEOTRACES'
            av[name]['MEDUSA']['source'] = 'MEDUSA'

            av[name]['MEDUSA']['details'] = {
                'name': name,
                'vars': [
                    'FER',
                ],
                'convert': ukp.mul1000,
                'units': 'umol F/m^3'
            }
            av[name]['Data']['details'] = {
                'name': name,
                'vars': [
                    'Fe_D_CONC_BOTTLE',
                ],
                'convert': ukp.NoChange,
            }  # no units?

        if 'O2' in analysisKeys:
            name = 'Oxygen'
            if annual:
                av[name]['MEDUSA']['File'] = listModelDataFiles(
                    jobID, 'ptrc_T', paths.ModelFolder_pref, annual, year)
                av[name]['Data']['File'] = WOAFolder + 'woa13_all_o00_01.nc'
            else:
                av[name]['Data']['File'] = WOAFolder + 'oxygen-woa13.nc'
                av[name]['MEDUSA'][
                    'File'] = ModelFolder + jobID + "_" + year + "_OXY.nc"

            av[name]['MEDUSA']['grid'] = modelGrid
            av[name]['depthLevels'] = depthLevels
            if annual: av[name]['plottingSlices'] = tsRegions
            else: av[name]['plottingSlices'] = HighLatWinter

            av[name]['Data']['coords'] = woaCoords
            av[name]['MEDUSA']['coords'] = medusaCoords

            av[name]['Data']['source'] = 'WOA'
            av[name]['MEDUSA']['source'] = 'MEDUSA'

            av[name]['MEDUSA']['details'] = {
                'name': name,
                'vars': [
                    'OXY',
                ],
                'convert': ukp.NoChange,
            }
            av[name]['Data']['details'] = {
                'name': name,
                'vars': [
                    'o_an',
                ],
                'convert': ukp.oxconvert,
                'units': 'mmol/m^3'
            }

        if 'Alk' in analysisKeys:
            name = 'Alkalinity'

            def convertmeqm3TOumolkg(nc, keys):
                return nc.variables[keys[0]][:] * 1.027

            if annual:
                av[name]['MEDUSA']['File'] = listModelDataFiles(
                    jobID, 'ptrc_T', paths.ModelFolder_pref, annual, year)
                av[name]['Data']['File'] = paths.GlodapDir + 'Alk.nc'
            else:
                print("Alkalinity data not available for monthly Analysis")
                assert 0

            av[name]['MEDUSA']['grid'] = modelGrid
            av[name]['depthLevels'] = [
                'Surface',
                '500m',
                '1000m',
                'Transect',
                'PTransect',
                'SOTransect',
                'AntTransect',
            ]

            alkregions = [
                'Global',
                'Equator10',
                'Remainder',
                'NorthernSubpolarAtlantic',
                'NorthernSubpolarPacific',
                'ignoreInlandSeas',
                'SouthernOcean',
                'AtlanticSOcean',
            ]
            #### very little arctic alkalinty
            if annual: av[name]['plottingSlices'] = alkregions
            else: av[name]['plottingSlices'] = HighLatWinter

            av[name]['Data']['coords'] = glodapCoords
            av[name]['MEDUSA']['coords'] = medusaCoords

            av[name]['Data']['source'] = 'GLODAP'
            av[name]['MEDUSA']['source'] = 'MEDUSA'

            av[name]['MEDUSA']['details'] = {
                'name': name,
                'vars': [
                    'ALK',
                ],
                'convert': ukp.NoChange,
                'units': 'meq/m^3',
            }
            av[name]['Data']['details'] = {
                'name': name,
                'vars': [
                    'Alk',
                ],
                'convert': convertmeqm3TOumolkg,
                'units': 'meq/m^3',
            }

        if 'DIC' in analysisKeys:
            name = 'DIC'

            if annual:
                av[name]['MEDUSA']['File'] = listModelDataFiles(
                    jobID, 'ptrc_T', paths.ModelFolder_pref, annual, year)
                av[name]['Data'][
                    'File'] = paths.GLODAPv2Dir + 'GLODAPv2.tco2.historic.nc'
            else:
                print("DIC data not available for monthly Analysis")
                assert 0

            av[name]['MEDUSA']['grid'] = modelGrid
            av[name]['depthLevels'] = depthLevels

            if annual: av[name]['plottingSlices'] = tsRegions
            else: av[name]['plottingSlices'] = HighLatWinter

            av[name]['Data']['coords'] = glodapv2Coords
            av[name]['MEDUSA']['coords'] = medusaCoords

            av[name]['Data']['source'] = 'GLODAPv2'
            av[name]['MEDUSA']['source'] = 'MEDUSA'

            av[name]['MEDUSA']['details'] = {
                'name': name,
                'vars': [
                    'DIC',
                ],
                'convert': ukp.NoChange,
                'units': 'mmol C/m^3'
            }
            av[name]['Data']['details'] = {
                'name': name,
                'vars': [
                    'tco2',
                ],
                'convert': ukp.convertkgToM3,
                'units': 'mmol C/m^3'
            }

        if 'IntPP_OSU' in analysisKeys:
            name = 'IntegratedPrimaryProduction_OSU'

            #####
            # Files:
            if annual:
                av[name]['MEDUSA']['File'] = listModelDataFiles(
                    jobID, 'diad_T', paths.ModelFolder_pref, annual, year)
                av[name]['Data'][
                    'File'] = paths.OSUDir + "/standard_VGPM.SeaWIFS.global.average.nc"
            else:
                print(
                    "IntegratedPrimaryProduction (OSU) data not available for monthly Analysis"
                )
                assert 0

            #####
            # Calculating depth in PP in medusa
            nc = Dataset(paths.orcaGridfn, 'r')
            area = nc.variables['e1t'][:] * nc.variables['e2t'][:]
            nc.close()

            def medusadepthInt(nc, keys):
                #	 mmolN/m2/d        [mg C /m2/d]   [mgC/m2/yr] [gC/m2/yr]     Gt/m2/yr
                factor = 1. * 6.625 * 12.011  #* 365.	      / 1000.   /     1E15
                arr = (nc.variables[keys[0]][:] +
                       nc.variables[keys[1]][:]) * factor

                #if arr.ndim ==3:
                #	for i in np.arange(arr.shape[0]):
                #		arr[i] = arr[i]*area
                #elif arr.ndim ==2: arr = arr*area
                #elif arr.ndim==1:
                #	index_x = nc.variables['index_x'][:]
                #	index_y = nc.variables['index_y'][:]
                #	for i,a in enumerate(arr):
                #		arr[i] = a * area[index_y[i],index_x[i]]
                #else: assert 0
                return arr

            #####
            # converting data to same units.
            nc = Dataset(av[name]['Data']['File'], 'r')
            lats = nc.variables['latitude'][:]
            osuareas = np.zeros((1080, 2160))
            osuarea = (111100. / 6.)**2.  # area of a pixel at equator. in m2
            for a in np.arange(1080):
                osuareas[a] = np.ones(
                    (2160, )) * osuarea * np.cos(np.deg2rad(lats[a]))

            def osuconvert(nc, keys):
                # Already in
                arr = nc.variables[keys[0]][:]
                #tlen = 1 # arr.shape[0]
                #arr  = arr/tlen * 365.	/ 1000. /     1E15
                #if arr.ndim ==3:
                #	for i in np.arange(arr.shape[0]):
                #		arr[i] = arr[i]*osuareas
                #elif arr.ndim ==2: arr = arr*osuareas
                #elif arr.ndim ==1:
                #	index_x = nc.variables['index_x'][:]
                #	index_y = nc.variables['index_y'][:]
                #	for i,a in enumerate(arr):
                #		#print i,a,[index_y[i],index_x[i]]
                #		arr[i] = a * osuareas[index_y[i],index_x[i]]
                #else:
                #	assert 0
                return arr

            av[name]['MEDUSA']['coords'] = medusaCoords
            av[name]['Data']['coords'] = osuCoords

            av[name]['MEDUSA']['grid'] = modelGrid
            av[name]['depthLevels'] = [
                '',
            ]

            if annual: av[name]['plottingSlices'] = tsRegions
            else: av[name]['plottingSlices'] = HighLatWinter

            av[name]['Data']['source'] = 'OSU'
            av[name]['MEDUSA']['source'] = 'MEDUSA'

            av[name]['MEDUSA']['details'] = {
                'name': name,
                'vars': ['PRN', 'PRD'],
                'convert': medusadepthInt,
                'units': 'mgC/m^2/day'
            }
            av[name]['Data']['details'] = {
                'name': name,
                'vars': [
                    'NPP',
                ],
                'convert': osuconvert,
                'units': 'mgC/m^2/day'
            }

        if 'AirSeaFlux' in analysisKeys:

            name = 'AirSeaFluxCO2'
            if annual:
                av[name]['MEDUSA']['File'] = listModelDataFiles(
                    jobID, 'diad_T', paths.ModelFolder_pref, annual, year)
                av[name]['Data'][
                    'File'] = paths.TakahashiFolder + 'takahashi_2009_Anual_sumflux_2006c_noHead.nc'
            else:
                av[name]['Data'][
                    'File'] = paths.TakahashiFolder + 'takahashi2009_month_flux_pCO2_2006c_noHead.nc'

                print("Air Sea Flux CO2 monthly not implemented")
                assert 0

            def eOrcaTotal(nc, keys):
                factor = 12. / 1000.  #/ 1.E12
                arr = nc.variables['CO2FLUX'][:].squeeze()  # mmolC/m2/d
                return arr * factor

            def takaTotal(nc, keys):
                arr = nc.variables['TFLUXSW06'][:].squeeze(
                )  # 10^12 g Carbon year^-1
                arr = -1.E12 * arr / 365.  #g Carbon/day
                area = nc.variables['AREA_MKM2'][:].squeeze(
                ) * 1E12  # 10^6 km^2
                fluxperarea = arr / area
                return fluxperarea

            av[name]['MEDUSA']['coords'] = medusaCoords
            av[name]['Data']['coords'] = takahashiCoords

            av[name]['MEDUSA']['grid'] = modelGrid
            av[name]['depthLevels'] = [
                '',
            ]

            if annual: av[name]['plottingSlices'] = tsRegions
            else: av[name]['plottingSlices'] = HighLatWinter

            av[name]['Data']['source'] = 'Takahashi2009'
            av[name]['MEDUSA']['source'] = 'MEDUSA'

            av[name]['MEDUSA']['details'] = {
                'name': name,
                'vars': [
                    'CO2FLUX',
                ],
                'convert': eOrcaTotal,
                'units': 'g C/m2/yr'
            }
            av[name]['Data']['details'] = {
                'name': name,
                'vars': ['TFLUXSW06', 'AREA_MKM2'],
                'convert': takaTotal,
                'units': 'g C/m2/yr'
            }

#		#if 'PCO2:
#			av['pCO2']['Data']['File'] 	= paths.TakahashiFolder + "takahashi2009_month_flux_pCO2_2006c_noHead.nc"
#			av['pCO2']['MEDUSA']['File'] 	= ModelFolder+"medusa_bio_"+year+".nc"
#			av['pCO2']['Data']['Vars'] 	= ['PCO2_SW',] 		#l+'_mn',
#			av['pCO2']['MEDUSA']['Vars'] 	= ['OCN_PCO2',]
#			av['pCO2']['depthLevels'] 	= ['',]
#			av['pCO2']['MEDUSA']['grid']	= modelGrid
#			#av['pCO2']['plottingSlices'] 	= []

        if 'S' in analysisKeys:
            name = 'Salinity'
            if annual:
                av[name]['NEMO']['File'] = listModelDataFiles(
                    jobID, 'grid_T', paths.ModelFolder_pref, annual, year)
                av[name]['Data'][
                    'File'] = WOAFolder + 'woa13_decav_s00_01v2.nc'
            else:
                av[name]['Data'][
                    'File'] = WOAFolder + 'salinity_monthly_1deg.nc'
                av[name]['NEMO'][
                    'File'] = ModelFolder + jobID + "_" + year + '_SAL.nc'

            av[name]['NEMO']['grid'] = modelGrid
            av[name]['depthLevels'] = depthLevels
            av[name]['plottingSlices'] = tsRegions

            av[name]['Data']['coords'] = woaCoords
            av[name]['NEMO']['coords'] = medusaCoords
            av[name]['Data']['source'] = 'WOA'
            av[name]['NEMO']['source'] = 'NEMO'

            av[name]['NEMO']['details'] = {
                'name': name,
                'vars': [
                    ukesmkeys['sal3d'],
                ],
                'convert': ukp.NoChange,
            }
            av[name]['Data']['details'] = {
                'name': name,
                'vars': [
                    's_an',
                ],
                'convert': ukp.NoChange,
            }  # no units?

        if 'T' in analysisKeys:
            name = 'Temperature'
            if annual:
                av[name]['NEMO']['File'] = listModelDataFiles(
                    jobID, 'grid_T', paths.ModelFolder_pref, annual, year)
                av[name]['Data'][
                    'File'] = WOAFolder + 'woa13_decav_t00_01v2.nc'
            else:
                av[name]['Data'][
                    'File'] = WOAFolder + 'temperature_monthly_1deg.nc'
                av[name]['NEMO'][
                    'File'] = ModelFolder + jobID + "_" + year + '_TEMP.nc'

            av[name]['NEMO']['grid'] = modelGrid
            av[name][
                'depthLevels'] = depthLevels  #['Surface','Transect','PTransect','SOTransect','ArcTransect','1000m',]
            av[name]['plottingSlices'] = tsRegions

            av[name]['Data']['coords'] = woaCoords
            av[name]['NEMO']['coords'] = medusaCoords
            av[name]['Data']['source'] = 'WOA'
            av[name]['NEMO']['source'] = 'NEMO'

            av[name]['NEMO']['details'] = {
                'name': name,
                'vars': [
                    ukesmkeys['temp3d'],
                ],
                'convert': ukp.NoChange,
            }
            av[name]['Data']['details'] = {
                'name': name,
                'vars': [
                    't_an',
                ],
                'convert': ukp.NoChange,
            }  # no units?

        if 'ZonalCurrent' in analysisKeys:
            name = 'ZonalCurrent'
            if annual:
                av[name]['NEMO']['File'] = listModelDataFiles(
                    jobID, 'grid_U', paths.ModelFolder_pref, annual, year)
                av[name]['Data']['File'] = paths.GODASFolder + 'ucur.clim.nc'
            else:
                assert 0
#				av[name]['Data']['File'] 	= WOAFolder+'temperature_monthly_1deg.nc'
#				av[name]['NEMO']['File'] 	= ModelFolder+jobID+"_"+year+'_TEMP.nc'

            av[name]['NEMO']['grid'] = modelGrid
            av[name][
                'depthLevels'] = depthLevels  #['Surface','Transect','PTransect','SOTransect','ArcTransect','1000m',]
            av[name]['plottingSlices'] = tsRegions

            av[name]['Data']['coords'] = godasCoords
            av[name]['NEMO']['coords'] = medusaUCoords
            av[name]['Data']['source'] = 'GODAS'
            av[name]['NEMO']['source'] = 'NEMO'

            av[name]['NEMO']['details'] = {
                'name': name,
                'vars': [
                    ukesmkeys['u3d'],
                ],
                'convert': ukp.mul1000,
                'units': 'mm/s'
            }
            av[name]['Data']['details'] = {
                'name': name,
                'vars': [
                    'ucur',
                ],
                'convert': ukp.NoChange,
                'units': 'mm/s'
            }

        if 'MeridionalCurrent' in analysisKeys:
            name = 'MeridionalCurrent'
            if annual:
                av[name]['NEMO']['File'] = listModelDataFiles(
                    jobID, 'grid_V', paths.ModelFolder_pref, annual, year)
                av[name]['Data']['File'] = paths.GODASFolder + 'vcur.clim.nc'
            else:
                assert 0
#				av[name]['Data']['File'] 	= WOAFolder+'temperature_monthly_1deg.nc'
#				av[name]['NEMO']['File'] 	= ModelFolder+jobID+"_"+year+'_TEMP.nc'

            av[name]['NEMO']['grid'] = modelGrid
            av[name][
                'depthLevels'] = depthLevels  #['Surface','Transect','PTransect','SOTransect','ArcTransect','1000m',]
            av[name]['plottingSlices'] = tsRegions

            av[name]['Data']['coords'] = godasCoords
            av[name]['NEMO']['coords'] = medusaVCoords
            av[name]['Data']['source'] = 'GODAS'
            av[name]['NEMO']['source'] = 'NEMO'

            av[name]['NEMO']['details'] = {
                'name': name,
                'vars': [
                    ukesmkeys['v3d'],
                ],
                'convert': ukp.mul1000,
                'units': 'mm/s'
            }
            av[name]['Data']['details'] = {
                'name': name,
                'vars': [
                    'vcur',
                ],
                'convert': ukp.NoChange,
                'units': 'mm/s'
            }

        if 'VerticalCurrent' in analysisKeys:
            name = 'VerticalCurrent'
            if annual:
                av[name]['NEMO']['File'] = listModelDataFiles(
                    jobID, 'grid_W', paths.ModelFolder_pref, annual, year)
                av[name]['Data']['File'] = paths.GODASFolder + 'dzdt.clim.nc'
            else:
                assert 0


#				av[name]['Data']['File'] 	= WOAFolder+'temperature_monthly_1deg.nc'
#				av[name]['NEMO']['File'] 	= ModelFolder+jobID+"_"+year+'_TEMP.nc'

            av[name]['NEMO']['grid'] = modelGrid
            av[name][
                'depthLevels'] = depthLevels  #['Surface','Transect','PTransect','SOTransect','ArcTransect','1000m',]
            av[name]['plottingSlices'] = tsRegions

            av[name]['Data']['coords'] = godasCoords
            av[name]['NEMO']['coords'] = medusaWCoords
            av[name]['Data']['source'] = 'GODAS'
            av[name]['NEMO']['source'] = 'NEMO'

            av[name]['NEMO']['details'] = {
                'name': name,
                'vars': [
                    ukesmkeys['w3d'],
                ],
                'convert': ukp.mul1000000,
                'units': 'um/s'
            }
            av[name]['Data']['details'] = {
                'name': name,
                'vars': [
                    'dzdt',
                ],
                'convert': ukp.NoChange,
                'units': 'um/s'
            }

        if 'MLD' in analysisKeys:
            name = 'MLD'
            if annual:
                av[name]['NEMO']['File'] = listModelDataFiles(
                    jobID, 'grid_T', paths.ModelFolder_pref, annual, year)
                av[name]['Data'][
                    'File'] = paths.MLDFolder + "mld_DT02_c1m_reg2.0-annual.nc"
            else:
                av[name]['Data'][
                    'File'] = paths.MLDFolder + "mld_DT02_c1m_reg2.0.nc"
                av[name]['NEMO'][
                    'File'] = ModelFolder + jobID + "_" + year + '_MLD.nc'

            av[name]['NEMO']['grid'] = modelGrid
            av[name]['depthLevels'] = [
                '',
            ]
            av[name]['plottingSlices'] = tsRegions

            av[name]['Data']['coords'] = {
                't': 'index_t',
                'z': 'index_z',
                'lat': 'lat',
                'lon': 'lon',
                'cal': 'standard',
                'tdict': ukp.tdicts['ZeroToZero']
            }
            av[name]['NEMO']['coords'] = medusaCoords
            av[name]['Data']['source'] = 'IFREMER'
            av[name]['NEMO']['source'] = 'NEMO'

            av[name]['NEMO']['details'] = {
                'name': name,
                'vars': [
                    'somxl010',
                ],
                'convert': ukp.NoChange,
                'units': 'm'
            }
            av[name]['Data']['details'] = {
                'name': name,
                'vars': [
                    'mld',
                    'mask',
                ],
                'convert': ukp.applymask,
                'units': 'm'
            }  # no units?

        if 'Dust' in analysisKeys:
            name = 'Dust'
            av[name]['MEDUSA']['File'] = listModelDataFiles(
                jobID, 'diad_T', paths.ModelFolder_pref, annual, year)
            av[name]['Data'][
                'File'] = paths.Dustdir + 'mahowald.orca100_annual.nc'

            av[name]['MEDUSA']['coords'] = medusaCoords
            av[name]['Data']['coords'] = medusaCoords

            av[name]['MEDUSA']['details'] = {
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
                #		dust[:,:,194:256,295:348] = 0.
                #		dust[:,:,194:208,285:295] = 0.
                #		dust[:,:,188:216,290:304] = 0.
                return dust * 0.035 * 1.e6 * 0.00532 * (24. * 60. *
                                                        60.) / 55.845

            def modeldustsum(nc, keys):
                dust = nc.variables[keys[0]][:]
                #              dust[:,234:296,295:348] = 0.
                #             dust[:,234:248,285:295] = 0.
                #            dust[:,228:256,290:304] = 0.
                return dust * 1.E-12 * 365.25

            if annual:
                av[name]['Data']['details'] = {
                    'name': name,
                    'vars': [
                        'dust_ann',
                    ],
                    'convert': mahodatadust,
                    'units': 'mmol Fe/m2/d'
                }
            else:
                av[name]['Data']['details'] = {
                    'name': name,
                    'vars': [
                        'dust',
                    ],
                    'convert': mahodatadust,
                    'units': 'mmol Fe/m2/d'
                }

            av[name]['Data']['source'] = 'Mahowald'
            av[name]['MEDUSA']['source'] = 'MEDUSA'

            av[name]['MEDUSA']['grid'] = modelGrid
            av[name]['depthLevels'] = [
                '',
            ]
            av[name]['plottingSlices'] = tsRegions

        #if 'AOU' in analysisKeys:
        #	name = 'AOU'
        #	if annual:
        #		av[name]['NEMO']['File'] 	= listModelDataFiles(jobID, 'grid_T', paths.ModelFolder_pref, annual,year)
        #		av[name]['Data']['File'] 	= paths.MLDFolder+"mld_DT02_c1m_reg2.0-annual.nc"
        #	else:
        #		av[name]['Data']['File'] 	= paths.MLDFolder+"mld_DT02_c1m_reg2.0.nc"
        #		av[name]['NEMO']['File'] 	= ModelFolder+jobID+"_"+year+'_MLD.nc'

        #	av[name]['NEMO']['grid'] 		= modelGrid
        #	av[name]['depthLevels'] 		= ['',]
        #	av[name]['plottingSlices'] 		= tsRegions

        #	av[name]['Data']['coords'] 	= {'t':'index_t', 'z':'index_z','lat':'lat','lon':'lon','cal': 'standard','tdict':ukp.tdicts['ZeroToZero']}
        #	av[name]['NEMO']['coords']	= medusaCoords
        #	av[name]['Data']['source'] 	= 'IFREMER'
        #	av[name]['NEMO']['source']	= 'NEMO'

        #	av[name]['NEMO']['details']	= {'name': name, 'vars':['somxl010',], 'convert': ukp.NoChange,'units':'m'}
        #	av[name]['Data']['details']	= {'name': name, 'vars':['mld','mask',], 'convert': ukp.applymask,'units':'m'}	# no units?

        for model in models:
            workingDir = ukp.folder(paths.p2p_ppDir + '/' + model + '-' +
                                    jobID + '-' + year)
            imageFolder = ukp.folder(imgDir + '/' + jobID)

            shelvesAV.extend(
                testsuite_p2p(
                    model=model,
                    jobID=jobID,
                    year=year,
                    av=av,
                    plottingSlices=
                    [],  # set this so that testsuite_p2p reads the slice list from the av.
                    workingDir=workingDir,
                    imageFolder=imageFolder,
                    noPlots=
                    noPlots,  # turns off plot making to save space and compute time.
                    gridFile=paths.orcaGridfn,  # enforces custom gridfile.
                    annual=annual,
                    noTargets=noTargets,
                ))
        if len(list(av.keys())) == 1: return
        if noTargets: return
        ######
        # Summary Target diagrams:
        imageFold = ukp.folder(imageFolder + '/Targets/' + year + '/Summary')
        summaryTargets(shelvesAV, imageFold, year)


def single_p2p(jobID, key, year):
    print("single_p2p:", jobID, key, year)
    try:
        analysis_p2p(
            models=[
                'NEMO',
                'MEDUSA',
            ],
            jobID=jobID,
            years=[
                year,
            ],  #'2075','2076',
            modelGrid='eORCA1',
            annual=True,
            noPlots=False,
            analysisSuite=[
                key,
            ],
        )
    except:
        print("Failed single_p2p", (jobID, key, year))
        print("Error: %s" % sys.exc_info()[0])


def get_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    accepted_keys = ['debug', 'physics', 'leve2', 'annual']

    parser.add_argument('-j',
                        '--jobID',
                        nargs='+',
                        type=str,
                        default=None,
                        help='One or more UKESM Job IDs (automatically generated by cycl/rose suite).',
                        required=True,
                        )
    parser.add_argument('-k',
                        '--keys',
                        default=['level2', ],
                        nargs='+',
                        type=str,
                        help=''.join(['Runtime keys - each key links to a pre-determined list of variables to analyse. ',
                                      'Keys are: ', ', '.join( accepted_keys)]),
                        required=False,
                        )
    parser.add_argument('-y',
                        '--years',
                        default=['best',],
                        nargs='+',
                        type=str,
                        help=''.join(['Years to analysis. Default is "best". ']),
                        required=False,
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
    from ._version import __version__
    print(f'BGCVal2: {__version__}')
    args = get_args()

    jobIDs = args.jobID
    keys = args.keys
    years = args.years

    accepted_keys = ['debug', 'physics', 'leve2', 'annual']
    for jobID, year, suite in itertools.product(jobIDs, years, keys):
        if year == 'best':
            best_year = False
            for divby in [100, 50, 25, 10, 5, 1]:
                print("analysis_p2p:\t find best year", divby, year,best_year )
                if best_year:
                    continue
                best_year = findLastFinishedYear(jobID,
                                            dividby=divby,) #  numberfiles=numberfiles)
            if best_year == False:
                continue
            else:
                year = best_year

        print('analysis_p2p:', jobID, year, suite)
        single_p2p(jobID, suite, year)
        # analysis_p2p(
        #     models=[
        #         'NEMO',
        #         'MEDUSA',
        #     ],
        #     jobID=jobID,
        #     years=[
        #         year,
        #     ],
        #     modelGrid='eORCA1',
        #     annual=True,
        #     noPlots=False,
        #     analysisSuite=suite,
        # )
    print("Finished p2p... ")


if __name__ == "__main__":
    main()
