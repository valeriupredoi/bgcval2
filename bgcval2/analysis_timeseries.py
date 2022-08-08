#!/usr/bin/env python
#####!/usr/bin/python

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
.. module:: analysis_timeseries
   :platform: Unix
   :synopsis: A script to produce analysis for time series.

.. moduleauthor:: Lee de Mora <ledm@pml.ac.uk>

"""
import argparse
import matplotlib as mpl

mpl.use('Agg')

#####
# Load Standard Python modules:
from sys import exit
from os.path import exists
from calendar import month_name
from socket import gethostname
from glob import glob
from scipy.interpolate import interp1d
import numpy as np
import os, sys
from getpass import getuser

#####
# Load specific local code:
from . import UKESMpython as ukp
from .timeseries import timeseriesAnalysis
from .timeseries import profileAnalysis
from .timeseries import timeseriesTools as tst

from .bgcvaltools.mergeMonthlyFiles import mergeMonthlyFiles, meanDJF
from .bgcvaltools.AOU import AOU
from .bgcvaltools.dataset import dataset
from ._runtime_config import get_run_configuration

#####
# User defined set of paths pointing towards the datasets.
from .Paths.paths import paths_setter


#####
# Biogeochemistry keys
bgcKeys = []
if True:
    bgcKeys.append('N')  # WOA Nitrate
    bgcKeys.append('Si')  # WOA Siliate
    bgcKeys.append('O2')  # WOA Oxygen
    bgcKeys.append('Alk')  # Glodap Alkalinity
    bgcKeys.append('DIC')  # Globap tCO2
    bgcKeys.append('AirSeaFluxCO2')  # Air Sea Flux
    bgcKeys.append('TotalAirSeaFluxCO2')  # Total global air sea flux
    bgcKeys.append('IntPP_OSU')  # OSU Integrated primpary production
    bgcKeys.append('PP_OSU')  # OSU Integrated primpary production
    #	bgcKeys.append('LocalExportRatio')         # Export ratio (no data)
    #	bgcKeys.append('GlobalExportRatio')        # Export ratio (no data)
    bgcKeys.append('TotalOMZVolume')  # Total Oxygen Minimum zone Volume
    #	bgcKeys.append('OMZThickness')             # Oxygen Minimum Zone Thickness
    #	bgcKeys.append('OMZMeanDepth')             # Oxygen Minimum Zone mean depth
    bgcKeys.append(
        'VolumeMeanOxygen')  # Volune Mean regional oxygen concentation
    #	bgcKeys.append('AOU')                      # Apparent Oxygen Usage
    bgcKeys.append('Iron')  # Iron
    bgcKeys.append('Dust')  # Dust
    #	bgcKeys.append('TotalDust')                # Total Dust
    #	bgcKeys.append('DiaFrac')                  # Diatom Fraction
    #        bgcKeys.append('DTC')                      # Detrital carbon
    bgcKeys.append('CHL')  # Total Chlorophyll
    #        bgcKeys.append('DMS_ARAN')                      # Total Chlorophyll
    bgcKeys.append('pH')

bgcKeysDict = {i: n for i, n in enumerate(bgcKeys)}

#####
# Physical keys
physKeys = []
if True:

    physKeys.append('Temperature')  # WOA Temperature
    #physKeys.append('VolWeightedT')			# Volume weighted WOA Temperature
    physKeys.append('GlobalMeanTemperature')  # Global Mean Temperature
    physKeys.append('GlobalMeanTemperature_700')
    physKeys.append('GlobalMeanTemperature_2000')

    physKeys.append('VolumeMeanTemperature')  # Global Mean Temperature
    physKeys.append('GlobalMeanSalinity')  # Global Mean Salinity
    #	physKeys.append('IcelessMeanSST')    		# Global Mean Surface Temperature with no ice
    physKeys.append('Salinity')  # WOA Salinity
    physKeys.append('MLD')  # iFERMER Mixed Layer Depth

    #	physKeys.append('TotalIceArea')			# work in progress
    #	physKeys.append('NorthernTotalIceArea')		# work in progress
    #	physKeys.append('SouthernTotalIceArea')		# work in progress
    #	physKeys.append('WeddelTotalIceArea')
    physKeys.append('TotalIceExtent')  # work in progress
    physKeys.append('NorthernTotalIceExtent')  # work in progress
    physKeys.append('SouthernTotalIceExtent')  # work in progress
    #        physKeys.append('WeddelIceExent')       # work in progress
    #physKeys.append('NorthernMIZArea')
    #physKeys.append('SouthernMIZArea')
    #physKeys.append('TotalMIZArea')
    #        physKeys.append('NorthernMIZfraction')
    #        physKeys.append('SouthernMIZfraction')
    #        physKeys.append('TotalMIZfraction')

    physKeys.append('DrakePassageTransport')  # DrakePassageTransport
    #	physKeys.append('AMOC_32S')                 	# AMOC 32S
    physKeys.append('AMOC_26N')  # AMOC 26N
    #        physKeys.append('AMOC_26N_nomexico')            # AMOC 26N
    #	physKeys.append('ADRC_26N')                 	# ADRC 26N
    #	physKeys.append('ZonalCurrent')             	# Zonal Veloctity
    #	physKeys.append('MeridionalCurrent')        	# Meridional Veloctity
    #	physKeys.append('VerticalCurrent')          	# Vertical Veloctity

    #physKeys.append('FreshwaterFlux')  # Freshwater flux
    physKeys.append('sowaflup')  # Net Upward Water Flux
    #	physKeys.append('soicecov')			# Ice fraction

    #####
    # unused:

    #        physKeys.append('MaxMonthlyMLD')               # MLD Monthly max
    #        physKeys.append('MinMonthlyMLD')               # MLD Monthly min
    #        physKeys.append('HeatFlux')
    physKeys.append('TotalHeatFlux')
    physKeys.append('scvoltot')
    physKeys.append('soga')
    physKeys.append('thetaoga')

#       physKeys.append('WindStress')                   # Wind Stress
#       physKeys.append('sohefldo')                     # Net downward Water Flux
#       physKeys.append('sofmflup')                     # Water flux due to freezing/melting
#       physKeys.append('sosfldow')                     # Downward salt flux
#       physKeys.append('sossheig')                 # Sea surface height

physKeysDict = {i: n for i, n in enumerate(physKeys)}

fastKeys = []
if True:
    fastKeys.append('N')  # WOA Nitrate
    fastKeys.append('Si')  # WOA Siliate
    fastKeys.append('O2')  # WOA Oxygen
    fastKeys.append('Alk')  # Glodap Alkalinity
    fastKeys.append('DIC')  # Globap tCO2
    fastKeys.append('AirSeaFluxCO2')  # Air Sea Flux
    fastKeys.append('TotalAirSeaFluxCO2')  # Total global air sea flux
    fastKeys.append('IntPP_OSU')  # OSU Integrated primpary production
    fastKeys.append('PP_OSU')  # OSU Integrated primpary production
    # fastKeys.append('TotalOMZVolume')           # Total Oxygen Minimum zone Volume
    fastKeys.append(
        'VolumeMeanOxygen')  # Volune Mean regional oxygen concentation
    fastKeys.append('Iron')  # Iron
    fastKeys.append('Dust')  # Dust
    fastKeys.append('CHL')  # Total Chlorophyll
    fastKeys.append('pH')

    fastKeys.append('Temperature')  # WOA Temperature
    fastKeys.append('GlobalMeanTemperature')  # Global Mean Temperature
    fastKeys.append('GlobalMeanTemperature_700')
    fastKeys.append('GlobalMeanTemperature_2000')

    fastKeys.append('VolumeMeanTemperature')  # Global Mean Temperature
    fastKeys.append('GlobalMeanSalinity')  # Global Mean Salinity
    fastKeys.append('Salinity')  # WOA Salinity
    fastKeys.append('MLD')  # iFERMER Mixed Layer Depth

    fastKeys.append('TotalIceExtent')  # work in progress
    fastKeys.append('NorthernTotalIceExtent')  # work in progress
    fastKeys.append('SouthernTotalIceExtent')  # work in progress

    fastKeys.append('DrakePassageTransport')  # DrakePassageTransport
    fastKeys.append('AMOC_26N')  # AMOC 26N

    fastKeys.append('FreshwaterFlux')  # Freshwater flux
    fastKeys.append('sowaflup')  # Net Upward Water Flux

    fastKeys.append('TotalHeatFlux')
    fastKeys.append('scvoltot')
    fastKeys.append('soga')
    fastKeys.append('thetaoga')

#####
# Level 1 keys
level1Keys = []
level1Keys.extend(physKeys)
level1Keys.extend(bgcKeys)
level1KeysDict = {i: n for i, n in enumerate(level1Keys)}

#####
# The important keys
keymetricsfirstKeys = [
    'TotalAirSeaFluxCO2',
    #                'NoCaspianAirSeaFluxCO2',
    #		'IntPP_OSU',
    #		'GlobalExportRatio',
    #                'TotalIceExtent',
    #                'NorthernTotalIceExtent',
    #                'SouthernTotalIceExtent',
    'DrakePassageTransport',
    'AMOC_26N',
    'GlobalMeanTemperature',
    #'GlobalMeanSalinity',
]
keymetricsfirstDict = {i: n for i, n in enumerate(keymetricsfirstKeys)}


def listModelDataFiles(jobID, filekey, datafolder, annual):
    print("listing model data files:\njobID:\t", jobID, '\nfile key:\t',
         filekey, '\ndata folder:\t', datafolder, '\nannual flag:\t', annual)
    if annual:
        print("listing model data files:",
               datafolder + jobID + "/" + jobID + "o_1y_*_" + filekey + ".nc")
        datafolder = os.path.join(datafolder, jobID)
        model_files = sorted(
            glob(datafolder + "/" + jobID + "o_1y_*_" + filekey +
                 ".nc"))
    else:
        print("listing model data files:",
               datafolder + jobID + "/" + jobID + "o_1m_*_" + filekey + ".nc")
        datafolder = os.path.join(datafolder, jobID)
        model_file = sorted(
            glob(datafolder + "/" + jobID + "o_1m_*_" + filekey +
                 ".nc"))

    return model_files


def analysis_timeseries(
    jobID="u-ab671",
    analysisSuite='all',
    regions='all',
    clean=0,
    annual=True,
    strictFileCheck=True,
    config_user=None
):
    """
	The role of this code is to produce time series analysis.
	The jobID is the monsoon/UM job id and it looks for files with a specific format

	The clean flag allows you to start the analysis without loading previous data.

	The annual flag means that we look at annual (True) or monthly (False) data.

	The strictFileCheck switch checks that the data/model netcdf files exist.
	It fails if the switch is on and the files no not exist.

	analysisSuite chooses a set of fields to look at.

	regions selects a list of regions, default is 'all', which is the list supplied by Andy Yool.

	:param jobID: the jobID
	:param clean: deletes old images if true
	:param annual: Flag for monthly or annual model data.
	:param strictFileCheck: CStrickt check for model and data files. Asserts if no files are found.
	:param analysisSuite: Which data to analyse, ie level1, physics only, debug, etc
	:param regions:

	"""

    print('-----------------------')
    print('Starting analysis_timeseries')
    print('jobID:', jobID)
    print('analysisSuite:',analysisSuite)
    print('regions:', regions)
    print('clean:',  clean, 'annual:',annual, 'strictFileCheck:', strictFileCheck)
    print('config_user:', config_user)   

    # get runtime configuration
    if config_user:
        paths_dict, config_user = get_run_configuration(config_user)
    else:
        paths_dict, config_user = get_run_configuration("defaults")

    # filter paths dict into an object that's usable below
    paths = paths_setter(paths_dict)

    #print "analysis_p2p:",	jobID,clean, annual,strictFileCheck,analysisSuite,regions
    #assert 0

    #####
    # Switches:
    # These are some booleans that allow us to choose which analysis to run.
    # This lets up give a list of keys one at a time, or in parrallel.
    if type(analysisSuite) == type(['Its', 'A', 'list!']):
        analysisKeys = analysisSuite

    #####
    # Switches:
    # These are some preset switches to run in series.
    if type(analysisSuite) == type('Its_A_string'):
        analysisKeys = []

        if analysisSuite.lower() in [
                'keymetricsfirst', 'kmf',
        ]:
            analysisKeys.extend(keymetricsfirstKeys)

        if analysisSuite.lower() in [
                'level1',
        ]:
            analysisKeys.extend(level1Keys)

        if analysisSuite.lower() in [
                'fast',
        ]:
            analysisKeys.extend(fastKeys)

        if analysisSuite.lower() in [
                'bgc',
        ]:
            analysisKeys.extend(bgcKeys)

        if analysisSuite.lower() in [
                'physics',
        ]:
            analysisKeys.extend(physKeys)

        if analysisSuite.lower() in [
                'level3',
        ]:
            analysisKeys.append('DMS_ARAN')  # DMS Aranami Tsunogai

        if analysisSuite.lower() in [
                'spinup',
        ]:
            analysisKeys.append('O2')  # WOA Oxygen
            analysisKeys.append('DIC')  # work in progress
            analysisKeys.append('Alk')  # Glodap Alkalinity
            analysisKeys.append('Iron')  # work in progress
            analysisKeys.append('N')  # WOA Nitrate
            analysisKeys.append('Si')  # WOA Nitrate
            analysisKeys.append('Temperature')  #             # WOA Temperature
            analysisKeys.append('Salinity')  #             # WOA Salinity
        if analysisSuite.lower() in [
                'salinity',
        ]:
            analysisKeys.append('Salinity')  #             # WOA Salinity

        if analysisSuite.lower() in [
                'debug',
        ]:
            analysisKeys.append('AMOC_26N')                # AMOC 26N

            #analysisKeys.append('AirSeaFlux')		# work in progress
            #analysisKeys.append('TotalAirSeaFluxCO2')	# work in progress
            #analysisKeys.append('NoCaspianAirSeaFluxCO2')	# work in progress
            #analysisKeys.append('TotalOMZVolume')		# work in progress
            #analysisKeys.append('TotalOMZVolume50')	# work in progress
            #analysisKeys.append('OMZMeanDepth')		# work in progress
            #analysisKeys.append('OMZThickness')            # Oxygen Minimum Zone Thickness
            #analysisKeys.append('TotalOMZVolume')		# work in progress
            #analysisKeys.append('O2')  # WOA Oxygen
            #analysisKeys.append('AOU')                      # Apparent Oxygen Usage
            #analysisKeys.append('WindStress')               # Wind Stress
            #analysisKeys.append('Dust')                    # Dust
            #analysisKeys.append('TotalDust')               # Total Dust
            #analysisKeys.append('TotalDust_nomask')
            #analysisKeys.append('DIC')  # work in progress
            #analysisKeys.append('DrakePassageTransport')	# DrakePassageTransport
            #analysisKeys.append('TotalIceArea')		# work in progress
            #analysisKeys.append('CHN')
            #analysisKeys.append('CHD')
            #analysisKeys.append('CHL')
            #analysisKeys.append('pH')
            #analysisKeys.append('Alk')  # Glodap Alkalinity

            #if jobID in ['u-am004','u-am005']:
            #        analysisKeys.append('DMS_ANDR')                 # DMS Anderson
            #else:   analysisKeys.append('DMS_ARAN')                 # DMS Aranami Tsunogai

            #analysisKeys.append('DiaFrac')			# work in progress
            #analysisKeys.append('Iron')			# work in progress
            #analysisKeys.append('DTC')                 # work in progress

            #analysisKeys.append('Iron')  # work in progress
            #analysisKeys.append('N')  # WOA Nitrate
            #analysisKeys.append('Si')  # WOA Nitrate
            #analysisKeys.append('IntPP_OSU')               # OSU Integrated primpary production
            #analysisKeys.append('Chl_CCI')
        #analysisKeys.append('CHL_MAM')
        #analysisKeys.append('CHL_JJA')
        #analysisKeys.append('CHL_SON')
        #analysisKeys.append('CHL_DJF')
        #analysisKeys.append('GC_CHL_MAM')
        #analysisKeys.append('GC_CHL_JJA')
        #analysisKeys.append('GC_CHL_SON')
        #analysisKeys.append('GC_CHL_DJF')
        #####
        # Physics switches:
        # analysisKeys.append('Temperature')          #             # WOA Temperature
#                        analysisKeys.append('HeatFlux')
#                        analysisKeys.append('TotalHeatFlux')

#                        analysisKeys.append('scvoltot')
#                        analysisKeys.append('soga')
#                        analysisKeys.append('thetaoga')
#                        analysisKeys.append('scalarHeatContent')

#analysisKeys.append('VolumeMeanTemperature')#
#                        analysisKeys.append('GlobalMeanTemperature_700')
#                        analysisKeys.append('GlobalMeanTemperature_2000')
#			analysisKeys.append('WeddelIceExent')
#analysisKeys.append('Salinity')                        # WOA Salinity
#                        analysisKeys.append('MLD')                      # MLD
#analysisKeys.append('MaxMonthlyMLD')            # MLD
#analysisKeys.append('MinMonthlyMLD')

#analysisKeys.append('NorthernTotalIceArea')    # work in progress
#analysisKeys.append('SouthernTotalIceArea')    # work in progress
#analysisKeys.append('WeddelTotalIceArea')
#                        analysisKeys.append('NorthernMIZArea')
#                        analysisKeys.append('SouthernMIZArea')
#                        analysisKeys.append('TotalMIZArea')
#analysisKeys.append('NorthernMIZfraction')
#analysisKeys.append('SouthernMIZfraction')
#analysisKeys.append('TotalMIZfraction')

#analysisKeys.append('TotalIceArea')            # work in progress
#analysisKeys.append('TotalIceExtent')		# work in progress
#analysisKeys.append('NorthernTotalIceExtent')	# work in progress
#analysisKeys.append('SouthernTotalIceExtent')	# work in progress
#analysisKeys.append('AMOC_32S')                # AMOC 32S
#analysisKeys.append('AMOC_26N')                # AMOC 26N
#analysisKeys.append('AMOC_26N_nomexico')
#analysisKeys.append('ADRC_26N')                # AMOC 26N

#                       analysisKeys.append('ERSST')    		# Global Surface Mean Temperature
#analysisKeys.append('VolumeMeanOxygen')

#                        analysisKeys.append('GlobalMeanTemperature')    # Global Mean Temperature#
#			analysisKeys.append('GlobalMeanSalinity')    	# Global Mean Salinity
#analysisKeys.append('IcelessMeanSST')    	# Global Mean Surface Temperature with no ice
#analysisKeys.append('quickSST')    		# Area Weighted Mean Surface Temperature

#analysisKeys.append('ZonalCurrent')             # Zonal Veloctity
#analysisKeys.append('MeridionalCurrent')        # Meridional Veloctity
#analysisKeys.append('VerticalCurrent')          # Vertical Veloctity

#analysisKeys.append('sowaflup')			# Net Upward Water Flux
#analysisKeys.append('sohefldo')			# Net downward Water Flux
#			analysisKeys.append('sofmflup')			# Water flux due to freezing/melting
#			analysisKeys.append('sosfldow')			# Downward salt flux
#			analysisKeys.append('soicecov')			# Ice fraction
#                       analysisKeys.append('sossheig')                 # Sea surface height
#analysisKeys.append('FreshwaterFlux')		# Fresh water flux
#analysisKeys.append('max_soshfldo')		# Max short wave radiation.

#####
# Physics switches:
#    if jobID in [
#            'u-aj588',
#            'u-ak900',
#            'u-ar538',
#            'u-an869',
#            'u-ar977',
#    ]:
#        try:
#            analysisKeys.remove('FreshwaterFlux')
#        except:
#            pass

    #####
    # Some lists of region.
    # This are pre-made lists of regions that can be investigated.
    # Note that each analysis below can be given its own set of regions.
    layerList = [
        'Surface',
    ]  #'500m','1000m',]
    metricList = [
        'mean',
    ]  #'median', '10pc','20pc','30pc','40pc','50pc','60pc','70pc','80pc','90pc','min','max']

    if regions == 'all':
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
    if regions == 'short':
        regionList = [
            'Global',
            'SouthernHemisphere',
            'NorthernHemisphere',
        ]

    if analysisSuite.lower() == 'debug':
        regionList = ['Global', 'ArcticOcean']

    if analysisSuite.lower() in [
            'spinup',
            'salinity',
    ]:
        regionList = [
            'Global',
        ]
        metricList = [
            'mean',
        ]
        layerList = ['500m', '1000m', '2000m', '4000m']

    # Regions from Pierce 1995 - https://doi.org/10.1175/1520-0485(1995)025<2046:CROHAF>2.0.CO;2
    PierceRegions = [
        'Enderby',
        'Wilkes',
        'Ross',
        'Amundsen',
        'Weddel',
    ]
    OMZRegions = [
        'EquatorialPacificOcean', 'IndianOcean', 'EquatorialAtlanticOcean'
    ]  #'Ross','Amundsen','Weddel',]

    #if analysisSuite.lower() in ['debug',]:
    #        regionList      = ['Global', 'ArcticOcean',]

    #####
    # The z_component custom command:
    # This flag sets a list of layers and metrics.
    # It's not advised to run all the metrics and all the layers, as it'll slow down the analysis.
    # if z_component in ['SurfaceOnly',]:

    #	if z_component in ['FullDepth',]:
    #		layerList = [0,2,5,10,15,20,25,30,35,40,45,50,55,60,70,]
    #		metricList = ['mean','median',]

    #####
    # Location of images directory
    # the imagedir is where the analysis images will be saved.

    #####
    # Location of shelves folder
    # The shelve directory is where the intermediate processing files are saved in python's shelve format.
    # This allows us to put away a python open to be re-opened later.
    # This means that we can interupt the analysis without loosing lots of data and processing time,
    # or we can append new simulation years to the end of the analysis without starting from scratch each time.
    #shelvedir 	= ukp.folder('shelves/timeseries/'+jobID)

    #####
    # Location of data files.
    # The first thing that this function does is to check which machine it is being run.
    # This is we can run the same code on multiple machines withouht having to make many copies of this file.
    # So far, this has been run on the following machines:
    #	PML
    #	JASMIN
    #	Charybdis (Julien's machine at NOCS)
    #
    # Feel free to add other macihines onto this list, if need be.
    machinelocation = ''


    shelvedir = ukp.folder(paths.shelvedir + "/timeseries/" + jobID)

    #####
    # PML
    hostname = gethostname()

    if hostname.find('pmpc') > -1:
        print("analysis-timeseries.py:\tBeing run at PML on ", gethostname())

        imagedir = ukp.folder(paths.imagedir + '/' + jobID + '/timeseries')

        if annual: WOAFolder = paths.WOAFolder_annual
        else: WOAFolder = paths.WOAFolder

        #shelvedir 	= ukp.folder(paths.shelvedir+'/'+jobID+'/timeseries/'+jobID)
        #shelvedir = ukp.folder(paths.shelvedir + "/timeseries/" + jobID)
    #####
    # JASMIN
    if hostname.find('ceda.ac.uk') > -1 or hostname.find(
            'jasmin') > -1 or hostname.find('jc.rl.ac.uk') > -1:
        print("analysis-timeseries.py:\tBeing run at CEDA on ", hostname)
        #machinelocation = 'JASMIN'

        #try:	shelvedir 	= ukp.folder("/group_workspaces/jasmin2/ukesm/BGC_data/"+getuser()+"/shelves/timeseries/"+jobID)
        #except: shelvedir       =            "/group_workspaces/jasmin2/ukesm/BGC_data/"+getuser()+"/shelves/timeseries/"+jobID
        #try:
        #    shelvedir = ukp.folder("/gws/nopw/j04/ukesm/BGC_data/" +
        #                           getuser() + "/shelves/timeseries/" + jobID)
        #except:
        #    shelvedir = "/gws/nopw/j04/ukesm/BGC_data/" + getuser(
        #    ) + "/shelves/timeseries/" + jobID

        if annual: WOAFolder = paths.WOAFolder_annual
        else: WOAFolder = paths.WOAFolder

        try:
            imagedir = ukp.folder(paths.imagedir + '/' + jobID + '/timeseries')
        except:
            imagedir = paths.imagedir + '/' + jobID + '/timeseries'

    if hostname.find('monsoon') > -1:
        print("Please set up paths.py")
        assert 0

#print "analysis-timeseries.py:\tBeing run at the Met Office on ",gethostname()
#machinelocation = 'MONSOON'

#ObsFolder       = "/projects/ukesm/ldmora/BGC-data/"
#ModelFolder       = "/projects/ukesm/ldmora/UKESM"
#####
# Location of model files.
#MEDUSAFolder_pref       = ukp.folder(ModelFolder)

#####
# Location of data files.
#if annual:      WOAFolder       = ukp.folder(ObsFolder+"WOA/annual")
#else:           WOAFolder       = ukp.folder(ObsFolder+"WOA/")

#MAREDATFolder   = ObsFolder+"/MAREDAT/MAREDAT/"
#GEOTRACESFolder = ObsFolder+"/GEOTRACES/GEOTRACES_PostProccessed/"
#TakahashiFolder = ObsFolder+"/Takahashi2009_pCO2/"
#MLDFolder       = ObsFolder+"/IFREMER-MLD/"
#iMarNetFolder   = ObsFolder+"/LestersReportData/"
#GlodapDir       = ObsFolder+"/GLODAP/"
#GLODAPv2Dir     = ObsFolder+"/GLODAPv2/GLODAPv2_Mapped_Climatologies/"
#OSUDir          = ObsFolder+"OSU/"
#CCIDir          = ObsFolder+"CCI/"
#orcaGridfn      = ModelFolder+'/mesh_mask_eORCA1_wrk.nc'

#####
# Unable to find location of files/data.
    if not paths.machinelocation:
        print(
            "analysis-timeseries.py:\tFATAL:\tWas unable to determine location of host: ",
            gethostname())
        print("Please set up paths.py, based on Paths/paths_template.py")
        #FIXME this is just for local testing by V
        if gethostname() == "valeriu-PORTEGE-Z30-C":
            assert True
        else:
            assert False

#####
# Because we can never be sure someone won't randomly rename the
# time dimension without saying anything.
# if jobID in ['u-am515','u-am927','u-am064','u-an326',]:
    try:
        tmpModelFiles = listModelDataFiles(jobID, 'grid_T',
                                           paths.ModelFolder_pref, annual)
    except:
        print(
            "No grid_T Model files available to figure out what naming convention is used."
        )
        tmpModelFiles = []
    ukesmkeys = {}
    if len(tmpModelFiles):
        nctmp = dataset(tmpModelFiles[0], 'r')
        nctmpkeys = list(nctmp.variables.keys())
        nctmp.close()
        if 'votemper' in nctmpkeys:
            ukesmkeys['time'] = 'time_counter'
            ukesmkeys['temp3d'] = 'votemper'
            ukesmkeys['sst'] = ''
            ukesmkeys['sal3d'] = 'vosaline'
            ukesmkeys['sss'] = ''
            ukesmkeys['v3d'] = 'vomecrty'
            ukesmkeys['u3d'] = 'vozocrtx'
            ukesmkeys['e3u'] = 'e3u'
            ukesmkeys['w3d'] = 'vovecrtz'
            ukesmkeys['MLD'] = 'ssomxl010'
        elif 'thetao_con' in nctmpkeys: # Added keys for UKESM2.
            ukesmkeys['time'] = 'time_counter'
            ukesmkeys['temp3d'] = 'thetao_con'
            ukesmkeys['sst'] = 'tos_con'
            ukesmkeys['sal3d'] = 'so_abs'
            ukesmkeys['sss'] = 'sos_abs'
            ukesmkeys['v3d'] = 'vo'
            ukesmkeys['u3d'] = 'uo'
            ukesmkeys['e3u'] = 'thkcello'
            ukesmkeys['w3d'] = 'wo'
            ukesmkeys['MLD'] = 'somxzint1'
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
            ukesmkeys['MLD'] = 'ssomxl010'

#	else:
#                        ukesmkeys['time']       = 'time_centered'
#                        ukesmkeys['temp3d']     = 'thetao'
#                        ukesmkeys['sst']        = 'tos'
#                        ukesmkeys['sal3d']     = 'so'
#                        ukesmkeys['sss']        = 'sos'
#                        ukesmkeys['v3d']     = 'vo'
#                        ukesmkeys['u3d']     = 'uo'
#                        ukesmkeys['e3u']    = 'thkcello'
#                        ukesmkeys['w3d']     = 'wo'

#                        ukesmkeys['time'] = 'time_counter'
#                        ukesmkeys['temp3d']     = 'votemper'
#                        ukesmkeys['sst']        = ''
#                        ukesmkeys['sal3d']     = 'vosaline'
#                        ukesmkeys['sss']        = ''
#                        ukesmkeys['v3d']     = 'vomecrty'
#                        ukesmkeys['u3d']     = 'vozocrtx'
#                        ukesmkeys['e3u']    = 'e3u'
#                        ukesmkeys['w3d']     = 'vovecrtz'

#	if jobID > 'u-am514' and jobID not in ['u-an619','u-an629','u-an631','u-an869', 'u-an908', 'u-an911','u-an989',]:
#		# There are other changes here too.
#                #####
#                # Because we can never be sure someone won't randomly rename the
#                # time dimension without saying anything.
#		ukesmkeys={}
#                ukesmkeys['time'] 	= 'time_centered'
#		ukesmkeys['temp3d'] 	= 'thetao'
#                ukesmkeys['sst'] 	= 'tos'
#                ukesmkeys['sal3d']     = 'so'
#                ukesmkeys['sss']        = 'sos'
#                ukesmkeys['v3d']     = 'vo'
#                ukesmkeys['u3d']     = 'uo'
#                ukesmkeys['e3u']    = 'thkcello'
#                ukesmkeys['w3d']     = 'wo'
#
#        else:
#                ukesmkeys={}
#                ukesmkeys['time'] = 'time_counter'
#                ukesmkeys['temp3d']     = 'votemper'
#                ukesmkeys['sst']        = ''
#                ukesmkeys['sal3d']     = 'vosaline'
#                ukesmkeys['sss']        = ''
#                ukesmkeys['v3d']     = 'vomecrty'
#                ukesmkeys['u3d']     = 'vozocrtx'
#                ukesmkeys['e3u']    = 'e3u'
#                ukesmkeys['w3d']     = 'vovecrtz'

#####
# Coordinate dictionairy
# These are python dictionairies, one for each data source and model.
# This is because each data provider seems to use a different set of standard names for dimensions and time.
# The 'tdict' field is short for "time-dictionary".
#	This is a dictionary who's indices are the values on the netcdf time dimension.
#	The tdict indices point to a month number in python numbering (ie January = 0)
# 	An example would be, if a netcdf uses the middle day of the month as it's time value:
#		tdict = {15:0, 45:1 ...}

    timekey = ukesmkeys['time']
    medusaCoords = {
        't': timekey,
        'z': 'deptht',
        'lat': 'nav_lat',
        'lon': 'nav_lon',
        'cal': '360_day',
    }  # model doesn't need time dict.
    medusaUCoords = {
        't': timekey,
        'z': 'depthu',
        'lat': 'nav_lat',
        'lon': 'nav_lon',
        'cal': '360_day',
    }  # model doesn't need time dict.
    medusaVCoords = {
        't': timekey,
        'z': 'depthv',
        'lat': 'nav_lat',
        'lon': 'nav_lon',
        'cal': '360_day',
    }  # model doesn't need time dict.
    medusaWCoords = {
        't': timekey,
        'z': 'depthw',
        'lat': 'nav_lat',
        'lon': 'nav_lon',
        'cal': '360_day',
    }  # model doesn't need time dict.

    icCoords = {
        't': timekey,
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
    osuCoords = {
        't': 'index_t',
        'z': '',
        'lat': 'latitude',
        'lon': 'longitude',
        'cal': 'standard',
        'tdict': []
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
    dmsCoords = {
        't': 'time',
        'z': 'depth',
        'lat': 'Latitude',
        'lon': 'Longitude',
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
    mogcCoords = {
        't': 'index_t',
        'z': 'index_z',
        'lat': 'latitude',
        'lon': 'longitude',
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

    #	def listModelDataFiles(jobID, filekey, datafolder, annual):
    #		print "listing model data files:\njobID:\t",jobID, '\nfile key:\t',filekey,'\ndata folder:\t', datafolder, '\nannual flag:\t',annual
    #		if annual:
    #			print "listing model data files:",datafolder+jobID+"/"+jobID+"o_1y_*_"+filekey+".nc"
    #			return sorted(glob(datafolder+jobID+"/"+jobID+"o_1y_*_"+filekey+".nc"))
    #		else:
    #                        print "listing model data files:",datafolder+jobID+"/"+jobID+"o_1m_*_"+filekey+".nc"
    #			return sorted(glob(datafolder+jobID+"/"+jobID+"o_1m_*_"+filekey+".nc"))

    masknc = dataset(paths.orcaGridfn, 'r')
    tlandmask = masknc.variables['tmask'][:]
    masknc.close()

    def applyLandMask(nc, keys):
        #### works like no change, but applies a mask.
        return np.ma.masked_where(tlandmask == 0,
                                  nc.variables[keys[0]][:].squeeze())

    def applySurfaceMask(nc, keys):
        #### works like no change, but applies a mask.
        return np.ma.masked_where(tlandmask[0, :, :] == 0,
                                  nc.variables[keys[0]][:].squeeze())

    def applyLandMask1e3(nc, keys):
        return applyLandMask(nc, keys) * 1000.
#####
# The analysis settings:
# Below here is a list of analysis settings.
# The settings are passed to timeseriesAnalysis using a nested dictionary (called an autovivification, here).
#
# These analysis were switched on or off at the start of the function.
# Each analysis requires:
#	model files
#	data files
#	model and data coordinate dictionaries, (defines above)
#	model and data details (a set of instructions of what to analyse:
#		name: 		field name
#		vars:		variable names in the netcdf
#		convert: 	a function to manipuate the data (ie change units, or add two fields together.
#				There are some standard ones in UKESMPython.py, but you can write your own here.
#		units: 		the units after the convert function has been applied.
#		layers:		which vertical layers to look at (ie, surface, 100m etc...)
#		regions:	which regions to look at. Can be speficied here, or use a pre-defined list (from above)
#		metrics:	what metric to look at:  mean, median or sum
#		model and data source: 	the name of source of the model/data (for plotting)
#		model grid: 	the model grid, usually eORCA1
#		the model grid file: 	the file path for the model mesh file (contains cell area/volume/masks, etc)
#
#	Note that the analysis can be run with just the model, it doesn't require a data file.
#	If so, just set to data file to an empty string:
#		av[name]['dataFile']  = ''

    av = ukp.AutoVivification()
    if 'Chl_pig' in analysisKeys:
        name = 'Chlorophyll_pig'
        av[name]['modelFiles'] = sorted(
            glob(paths.ModelFolder_pref + jobID + "/" + jobID +
                 "o_1y_*_ptrc_T.nc"))
        av[name][
            'dataFile'] = paths.MAREDATFolder + "MarEDat20121001Pigments.nc"

        av[name]['modelcoords'] = medusaCoords
        av[name]['datacoords'] = maredatCoords

        av[name]['modeldetails'] = {
            'name': name,
            'vars': ['CHN', 'CHD'],
            'convert': ukp.sums,
            'units': 'mg C/m^3'
        }
        av[name]['datadetails'] = {
            'name': name,
            'vars': [
                'Chlorophylla',
            ],
            'convert': ukp.div1000,
            'units': 'ug/L'
        }

        av[name]['layers'] = layerList
        av[name]['regions'] = regionList
        av[name]['metrics'] = metricList

        av[name]['datasource'] = 'MAREDAT'
        av[name]['model'] = 'MEDUSA'

        av[name]['modelgrid'] = 'eORCA1'
        av[name]['gridFile'] = paths.orcaGridfn
        av[name]['Dimensions'] = 3

    if 'CHL' in analysisKeys:
        name = 'Chlorophyll'
        av[name]['modelFiles'] = listModelDataFiles(jobID, 'ptrc_T',
                                                    paths.ModelFolder_pref,
                                                    annual)
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

    if 'Chl_CCI' in analysisKeys:
        name = 'Chlorophyll_cci'
        #####
        # Not that this is the 1 degree resolution dataset, but higher resolution data are also available.

        av[name]['modelFiles'] = listModelDataFiles(jobID, 'ptrc_T',
                                                    paths.ModelFolder_pref,
                                                    annual)
        if annual:
            av[name][
                'dataFile'] = paths.CCIDir + "ESACCI-OC-L3S-OC_PRODUCTS-CLIMATOLOGY-16Y_MONTHLY_1degree_GEO_PML_OC4v6_QAA-annual-fv2.0.nc"
            print(paths.ModelFolder_pref + "/" + jobID + "o_1y_*_ptrc_T.nc")
        else:
            av[name][
                'dataFile'] = paths.CCIDir + 'ESACCI-OC-L3S-OC_PRODUCTS-CLIMATOLOGY-16Y_MONTHLY_1degree_GEO_PML_OC4v6_QAA-all-fv2.0.nc'

        av[name]['modelcoords'] = medusaCoords
        av[name]['datacoords'] = cciCoords

        av[name]['modeldetails'] = {
            'name': name,
            'vars': ['CHN', 'CHD'],
            'convert': ukp.sums,
            'units': 'mg C/m^3'
        }
        av[name]['datadetails'] = {
            'name': name,
            'vars': [
                'chlor_a',
            ],
            'convert': ukp.NoChange,
            'units': 'mg C/m^3'
        }

        av[name]['layers'] = [
            'Surface',
        ]  # CCI is surface only, it's a satellite product.
        av[name]['regions'] = regionList
        av[name]['metrics'] = metricList  #['mean','median', ]

        av[name]['datasource'] = 'CCI'
        av[name]['model'] = 'MEDUSA'

        av[name]['modelgrid'] = 'eORCA1'
        av[name]['gridFile'] = paths.orcaGridfn
        av[name]['Dimensions'] = 2

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
            ]  #'100m',] 	# CCI is surface only, it's a satellite product.
            av[name]['regions'] = regionList
            av[name]['metrics'] = metricList  #['mean','median', ]

            av[name]['datasource'] = ''
            av[name]['model'] = 'MEDUSA'

            av[name]['modelgrid'] = 'eORCA1'
            av[name]['gridFile'] = paths.orcaGridfn
            av[name]['Dimensions'] = 3

    for name in [
            'CHL_JJA', 'CHL_DJF', 'CHL_SON', 'CHL_MAM', 'GC_CHL_JJA',
            'GC_CHL_DJF', 'GC_CHL_SON', 'GC_CHL_MAM'
    ]:
        if name not in analysisKeys: continue

        if name in ['CHL_MAM', 'GC_CHL_MAM']:
            monthlyFiles = glob(paths.ModelFolder_pref + '/' + jobID +
                                '/monthlyCHL/' + jobID +
                                'o_1m_????0[345]*_ptrc_T.nc')
        if name in [
                'CHL_JJA',
                'GC_CHL_JJA',
        ]:
            monthlyFiles = glob(paths.ModelFolder_pref + '/' + jobID +
                                '/monthlyCHL/' + jobID +
                                'o_1m_????0[678]*_ptrc_T.nc')
        if name in [
                'CHL_SON',
                'GC_CHL_SON',
        ]:
            monthlyFiles = []
            for month in ['09', '10', '11']:
                monthlyFiles.extend(
                    glob(paths.ModelFolder_pref + '/' + jobID +
                         '/monthlyCHL/' + jobID + 'o_1m_????' + month +
                         '*_ptrc_T.nc'))
        if name in ['CHL_DJF', 'GC_CHL_DJF']:
            monthlyFiles = []
            for month in ['12', '01', '02']:
                monthlyFiles.extend(
                    glob(paths.ModelFolder_pref + '/' + jobID +
                         '/monthlyCHL/' + jobID + 'o_1m_????' + month +
                         '*_ptrc_T.nc'))

        if len(monthlyFiles):
            if name in [
                    'CHL_JJA', 'CHL_SON', 'CHL_MAM', 'GC_CHL_JJA',
                    'GC_CHL_SON', 'GC_CHL_MAM'
            ]:
                chlfiles = mergeMonthlyFiles(monthlyFiles,
                                             outfolder='',
                                             cal=medusaCoords['cal'],
                                             timeAverage=True,
                                             expectedNumberOfFiles=3)
            if name in [
                    'CHL_DJF',
                    'GC_CHL_DJF',
            ]:
                chlfiles = meanDJF(monthlyFiles,
                                   outfolder='',
                                   cal=medusaCoords['cal'],
                                   timeAverage=True)

            def CHL_MODEL(nc, keys):
                chl = nc.variables[keys[0]][:, 0] + nc.variables[keys[1]][:, 0]
                chl = np.ma.array(chl)
                chl = np.ma.masked_where(chl.mask + (chl > 10E10), chl)
                return chl

            def CHLMAM_Data(nc, keys):
                chl = nc.variables[keys[0]][2:5].mean(0).squeeze()
                chl = np.ma.array(chl)
                chl = np.ma.masked_where(chl.mask + (chl > 10E10), chl)
                return chl

            def CHLJJA_Data(nc, keys):
                chl = nc.variables[keys[0]][5:8].mean(0).squeeze()
                chl = np.ma.array(chl)
                chl = np.ma.masked_where(chl.mask + (chl > 10E10), chl)
                return chl

            def CHLSON_Data(nc, keys):
                chl = nc.variables[keys[0]][8:11].mean(0)
                chl = np.ma.masked_where(chl.mask + (chl > 10E10), chl)
                return chl

            def CHLDJF_Data(nc, keys):
                chl = []
                chl.append(nc.variables[keys[0]][0])
                chl.append(nc.variables[keys[0]][1])
                chl.append(nc.variables[keys[0]][11])
                chl = np.ma.array(chl).mean(0)
                chl = np.ma.masked_where(chl.mask + (chl > 10E10), chl)
                return chl

            av[name]['modelFiles'] = chlfiles
            av[name]['modelcoords'] = medusaCoords

            if name[:2] == 'GC':
                av[name][
                    'dataFile'] = paths.ObsFolder + 'MO-GlobColour/qrclim_globcolour_masked.sea.nc'
                av[name]['datacoords'] = mogcCoords
                av[name]['datasource'] = 'MO-GlobColour'
                chldatakey = 'chl'
            else:
                av[name][
                    'dataFile'] = paths.CCIDir + 'ESACCI-OC-L3S-OC_PRODUCTS-CLIMATOLOGY-16Y_MONTHLY_1degree_GEO_PML_OC4v6_QAA-all-fv2.0.nc'
                av[name]['datacoords'] = cciCoords
                av[name]['datasource'] = 'CCI'
                chldatakey = 'chlor_a'

            av[name]['modeldetails'] = {
                'name': name,
                'vars': ['CHN', 'CHD'],
                'convert': CHL_MODEL,
                'units': 'mg C/m^3'
            }

            if name in ['CHL_MAM', 'GC_CHL_MAM']:
                av[name]['datadetails'] = {
                    'name': name,
                    'vars': [
                        chldatakey,
                    ],
                    'convert': CHLMAM_Data,
                    'units': 'mg C/m^3'
                }
                av[name]['regions'] = [
                    'Equator10',
                    'Remainder',
                    'NorthernSubpolarAtlantic',
                    'NorthernSubpolarPacific',
                    'Global',
                    'SouthernOcean',
                ]  #'CCI_MAM']
            if name in [
                    'CHL_JJA',
                    'GC_CHL_JJA',
            ]:
                av[name]['datadetails'] = {
                    'name': name,
                    'vars': [
                        chldatakey,
                    ],
                    'convert': CHLJJA_Data,
                    'units': 'mg C/m^3'
                }
                av[name]['regions'] = [
                    'Equator10',
                    'Remainder',
                    'NorthernSubpolarAtlantic',
                    'NorthernSubpolarPacific',
                    'CCI_JJA',
                    'Global',
                    'SouthernOcean',
                ]

            if name in [
                    'CHL_SON',
                    'GC_CHL_SON',
            ]:
                av[name]['datadetails'] = {
                    'name': name,
                    'vars': [
                        chldatakey,
                    ],
                    'convert': CHLSON_Data,
                    'units': 'mg C/m^3'
                }
                av[name]['regions'] = [
                    'Equator10',
                    'Remainder',
                    'NorthernSubpolarAtlantic',
                    'NorthernSubpolarPacific',
                    'Global',
                    'CCI_SON',
                    'SouthernOcean',
                ]

            if name in ['CHL_DJF', 'GC_CHL_DJF']:
                av[name]['datadetails'] = {
                    'name': name,
                    'vars': [
                        chldatakey,
                    ],
                    'convert': CHLDJF_Data,
                    'units': 'mg C/m^3'
                }
                av[name]['regions'] = [
                    'Equator10',
                    'Remainder',
                    'NorthernSubpolarAtlantic',
                    'NorthernSubpolarPacific',
                    'CCI_DJF',
                    'Global',
                    'SouthernOcean',
                ]

            av[name]['layers'] = [
                'Surface',
            ]

            av[name]['metrics'] = metricList

            av[name]['model'] = 'NEMO'

            av[name]['modelgrid'] = 'eORCA1'
            av[name]['gridFile'] = paths.orcaGridfn
            av[name]['Dimensions'] = 2

        else:
            print("No monthly CHL files found")

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

    if 'DiaFrac' in analysisKeys:

        name = 'DiaFrac'

        def caldiafrac(nc, keys):
            chd = applyLandMask(nc, [
                keys[0],
            ]).squeeze()
            chn = applyLandMask(nc, [
                keys[1],
            ]).squeeze()
            return 100. * chd / (chd + chn)

        av[name]['modelFiles'] = listModelDataFiles(jobID, 'ptrc_T',
                                                    paths.ModelFolder_pref,
                                                    annual)
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
        av[name]['Dimensions'] = 3

    if 'N' in analysisKeys:
        name = 'Nitrate'
        av[name]['modelFiles'] = listModelDataFiles(jobID, 'ptrc_T',
                                                    paths.ModelFolder_pref,
                                                    annual)
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
        av[name]['gridFile'] = paths.orcaGridfn
        av[name]['Dimensions'] = 3

    if 'Si' in analysisKeys:
        name = 'Silicate'
        av[name]['modelFiles'] = listModelDataFiles(jobID, 'ptrc_T',
                                                    paths.ModelFolder_pref,
                                                    annual)
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
        av[name]['gridFile'] = paths.orcaGridfn
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

        if analysisSuite.lower() == 'spinup':
            av[name]['layers'] = layerList
        else:
            av[name]['layers'] = [
                'Surface',
                '500m',
                '1000m',
            ]  #layerList
        av[name]['regions'] = regionList
        av[name]['metrics'] = metricList

        av[name]['datasource'] = 'WOA'
        av[name]['model'] = 'MEDUSA'

        av[name]['modelgrid'] = 'eORCA1'
        av[name]['gridFile'] = paths.orcaGridfn
        av[name]['Dimensions'] = 3

    if 'OMZMeanDepth' in analysisKeys:
        if annual:
            av['OMZMeanDepth']['modelFiles'] = sorted(
                glob(paths.ModelFolder_pref + jobID + "/" + jobID +
                     "o_1y_*_ptrc_T.nc"))
            av['OMZMeanDepth']['dataFile'] = WOAFolder + 'woa13_all_o00_01.nc'
        else:
            print("OMZ Thickness not implemented for monthly data")
            assert 0

        nc = dataset(paths.orcaGridfn, 'r')
        depths = np.abs(nc.variables['gdepw'][:])
        tmask = nc.variables['tmask'][:]
        nc.close()

        omzthreshold = 20.

        def modelMeanOMZdepth(nc, keys):
            o2 = np.ma.array(nc.variables[keys[0]][:].squeeze())
            meandepth = np.ma.masked_where(
                (o2 > omzthreshold) + o2.mask + (tmask == 0), depths).mean(0)
            if meandepth.max() in [0., 0]: return np.array([
                    0.,
            ])
            return np.ma.masked_where(meandepth == 0., meandepth)

        def woaMeanOMZdepth(nc, keys):
            o2 = np.ma.array(nc.variables[keys[0]][:].squeeze() * 44.661)
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
            av['OMZThickness']['dataFile'] = WOAFolder + 'woa13_all_o00_01.nc'
        else:
            print("OMZ Thickness not implemented for monthly data")
            assert 0

        nc = dataset(paths.orcaGridfn, 'r')
        thickness = nc.variables['e3t'][:]
        tmask = nc.variables['tmask'][:]
        nc.close()

        if 'OMZThickness' in analysisKeys: omzthreshold = 20.
        if 'OMZThickness50' in analysisKeys: omzthreshold = 50.

        def modelOMZthickness(nc, keys):
            o2 = np.ma.array(nc.variables[keys[0]][:].squeeze())
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
            tomzv_dir = os.path.join(paths.ModelFolder_pref, jobID)
            av['TotalOMZVolume']['modelFiles'] = sorted(
                glob(tomzv_dir + "/" + jobID +
                     "o_1y_*_ptrc_T.nc"))
            av['TotalOMZVolume'][
                'dataFile'] = WOAFolder + 'woa13_all_o00_01.nc'
        else:
            print("OMZ volume not implemented for monthly data")
            assert 0

        nc = dataset(paths.orcaGridfn, 'r')
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
            arr = np.ma.array(nc.variables[keys[0]][:].squeeze() * 44.661)
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
        av[name]['modelFiles'] = listModelDataFiles(jobID, 'ptrc_T',
                                                    paths.ModelFolder_pref,
                                                    annual)
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

    if 'AOU' in analysisKeys:
        name = 'AOU'
        if annual:
            av[name]['modelFiles'] = listModelDataFiles(
                jobID, 'ptrc_T', paths.ModelFolder_pref, annual)
            av[name]['dataFile'] = ''  #WOAFolder+'woa13_all_o00_01.nc'

        def modelAOU(nc, keys):
            o2 = nc.variables[keys[0]][:]

            ncpath = nc.filename
            print(ncpath)

            newpath = ncpath.replace('ptrc_T', 'grid_T')

            print("modelAOU:", ncpath, newpath)
            nc2 = dataset(newpath, 'r')
            print("Loaded", newpath)
            temp = nc2.variables[ukesmkeys['temp3d']][:]
            sal = nc2.variables[ukesmkeys['sal3d']][:]
            return AOU(temp, sal, o2)

        av[name]['modelcoords'] = medusaCoords
        av[name]['datacoords'] = woaCoords

        av[name]['modeldetails'] = {
            'name': name,
            'vars': ['OXY', ukesmkeys['temp3d'], ukesmkeys['sal3d']],
            'convert': modelAOU,
            'units': 'mmol O2/m^3'
        }
        av[name]['datadetails'] = {
            'name': '',
            'units': '',
        }

        av[name]['layers'] = layerList
        av[name]['regions'] = regionList
        av[name]['metrics'] = metricList

        av[name]['datasource'] = ''
        av[name]['model'] = 'MEDUSA'

        av[name]['modelgrid'] = 'eORCA1'
        av[name]['gridFile'] = paths.orcaGridfn
        av[name]['Dimensions'] = 3

    if 'DIC' in analysisKeys:

        def convertkgToM3(nc, keys):
            return nc.variables[keys[0]][:] * 1.027

        name = 'DIC'

        av[name]['modelFiles'] = listModelDataFiles(jobID, 'ptrc_T',
                                                    paths.ModelFolder_pref,
                                                    annual)
        av[name]['dataFile'] = paths.GLODAPv2Dir + 'GLODAPv2.tco2.historic.nc'

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
        av[name]['gridFile'] = paths.orcaGridfn
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
        av[name]['gridFile'] = paths.orcaGridfn
        av[name]['Dimensions'] = 3

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
        if analysisSuite.lower() == 'fast':
            av[name]['Dimensions'] = 2
            av[name]['modeldetails'] = {
                'name': name,
                'vars': [
                    'OCN_PH',
                ],
                'convert': ukp.NoChange,
                'units': 'pH',
            }
        else:
            av[name]['Dimensions'] = 3
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
        av[name]['regions'] = regionList
        av[name]['metrics'] = metricList

        av[name]['model'] = 'MEDUSA'

        av[name]['modelgrid'] = 'eORCA1'
        av[name]['gridFile'] = paths.orcaGridfn
        av[name]['Dimensions'] = 3

    if 'AirSeaFluxCO2' in analysisKeys:

        #nc = dataset(paths.orcaGridfn,'r')
        #area = nc.variables['e1t'][:]*nc.variables['e2t'][:]
        #nc.close()

        def eOrcaTotal(nc, keys):
            factor = 12. / 1000.
            arr = nc.variables['CO2FLUX'][:].squeeze()  # mmolC/m2/d
            #if arr.ndim ==3:
            #	for i in np.arange(arr.shape[0]):
            #		arr[i] = arr[i]*area
            #elif arr.ndim ==2: arr = arr*area
            #else: assert 0
            return arr * factor

        def takaTotal(nc, keys):
            arr = nc.variables['TFLUXSW06'][:].squeeze(
            )  # 10^12 g Carbon year^-1
            arr = -1.E12 * arr / 365.  #g Carbon/day
            factor = -1.E12 / (365.)  # convert to #/ 1.E12
            area = nc.variables['AREA_MKM2'][:].squeeze() * 1E12  # 10^6 km^2
            fluxperarea = arr / area
            #arr = arr*area #* 1.E24 	# converts area into m^2
            #print arr.sum(), arr.sum()*factor
            return fluxperarea
            # area 10^6 km^2
            # flux:  10^15 g Carbon month^-1. (GT)/m2/month

        name = 'AirSeaFluxCO2'

        av[name]['modelFiles'] = listModelDataFiles(jobID, 'diad_T',
                                                    paths.ModelFolder_pref,
                                                    annual)
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

    if 'TotalAirSeaFluxCO2' in analysisKeys:
        name = 'TotalAirSeaFluxCO2'
        nc = dataset(paths.orcaGridfn, 'r')
        area = nc.variables['e1t'][:] * nc.variables['e2t'][:]
        nc.close()

        def eOrcaTotal(nc, keys):
            factor = 365.25 * 12. / 1000. / 1.E15
            arr = nc.variables['CO2FLUX'][:].squeeze() * factor  # mmolC/m2/d
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

        av[name]['modelFiles'] = listModelDataFiles(jobID, 'diad_T',
                                                    paths.ModelFolder_pref,
                                                    annual)
        av[name]['modelcoords'] = medusaCoords
        #if jobID == 'u-am515':	av[name]['modeldetails']        = {'name': 'AirSeaFluxCO2', 'vars':['CO2FLUX',], 'convert': uam515_TotalCO2,'units':'Pg C/yr'}
        #else:
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
            'regionless',
        ]
        av[name]['metrics'] = [
            'metricless',
        ]
        av[name]['model'] = 'MEDUSA'
        av[name]['modelgrid'] = 'eORCA1'
        av[name]['gridFile'] = paths.orcaGridfn
        av[name]['Dimensions'] = 2

        noTaka = True
        if noTaka:
            av[name]['datacoords'] = {'name': '', 'units': ''}
            av[name]['datadetails'] = {'name': '', 'units': ''}
            av[name]['dataFile'] = ''
            av[name]['datasource'] = ''
        else:
            if annual:
                av[name][
                    'dataFile'] = paths.TakahashiFolder + 'takahashi_2009_Anual_sumflux_2006c_noHead.nc'
            else:
                av[name][
                    'dataFile'] = paths.TakahashiFolder + 'takahashi2009_month_flux_pCO2_2006c_noHead.nc'
                print("Air Sea Flux CO2 monthly not implemented")
                assert 0
                #av[name]['dataFile'] 		=  paths.TakahashiFolder+'takahashi2009_month_flux_pCO2_2006c_noHead.nc'
            av[name]['datadetails'] = {
                'name': 'AirSeaFluxCO2',
                'vars': ['TFLUXSW06', 'AREA_MKM2'],
                'convert': takaTotal,
                'units': 'Pg C/yr'
            }
            av[name]['datasource'] = 'Takahashi'
            av[name]['datacoords'] = takahashiCoords

    if 'NoCaspianAirSeaFluxCO2' in analysisKeys:
        name = 'NoCaspianAirSeaFluxCO2'
        nc = dataset(paths.orcaGridfn, 'r')
        area = nc.variables['e1t'][:] * nc.variables['e2t'][:]
        nc.close()

        def eOrcaTotal(nc, keys):
            factor = 365.25 * 12. / 1000. / 1.E15
            arr = nc.variables['CO2FLUX'][:].squeeze() * factor  # mmolC/m2/d
            if arr.ndim == 3:
                for i in np.arange(arr.shape[0]):
                    arr[i] = arr[i] * area
            elif arr.ndim == 2:
                arr = arr * area
            else:
                assert 0
            return arr

        av[name]['modelFiles'] = listModelDataFiles(jobID, 'diad_T',
                                                    paths.ModelFolder_pref,
                                                    annual)

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

    if 'IntPP_iMarNet' in analysisKeys:
        name = 'IntegratedPrimaryProduction_1x1'

        def medusadepthInt(nc, keys):
            return (nc.variables[keys[0]][:] +
                    nc.variables[keys[1]][:]) * 6.625 * 12.011 / 1000.

        av[name]['modelFiles'] = listModelDataFiles(jobID, 'diad_T',
                                                    paths.ModelFolder_pref,
                                                    annual)
        #av[name]['modelFiles']  	= sorted(glob(paths.ModelFolder_pref+jobID+"/"+jobID+"o_1y_*_diad_T.nc"))
        av[name]['dataFile'] = paths.iMarNetFolder + "/PPint_1deg.nc"

        av[name]['modelcoords'] = medusaCoords
        av[name]['datacoords'] = maredatCoords

        av[name]['modeldetails'] = {
            'name': 'IntPP',
            'vars': ['PRN', 'PRD'],
            'convert': medusadepthInt,
            'units': 'mg C/m^3'
        }
        #av[name]['datadetails']  	= {'name': 'IntPP', 'vars':['Chlorophylla',], 'convert': ukp.div1000,'units':'ug/L'}
        av[name]['datadetails'] = {
            'name': 'IntPP',
            'vars': [
                'PPint',
            ],
            'convert': ukp.div1000,
            'units': '[ug/L/d'
        }

        av[name]['layers'] = [
            'Surface',
        ]  #'100m','200m','Surface - 1000m','Surface - 300m',]#'depthint']
        av[name]['regions'] = regionList
        av[name]['metrics'] = metricList

        av[name]['datasource'] = 'MAREDAT'
        av[name]['model'] = 'MEDUSA'

        av[name]['modelgrid'] = 'eORCA1'
        av[name]['gridFile'] = paths.orcaGridfn
        av[name]['Dimensions'] = 2

    if 'PP_OSU' in analysisKeys:
        nc = dataset(paths.orcaGridfn, 'r')
        area = nc.variables['e1t'][:] * nc.variables['e2t'][:]
        nc.close()

        def medusadepthInt(nc, keys):
            #	 mmolN/m2/d        [mg C /m2/d]   [mgC/m2/yr] [gC/m2/yr]     Gt/m2/yr
            factor = 1. * 6.625 * 12.011  #* 365.	      / 1000.   /     1E15
            arr = (nc.variables[keys[0]][:] +
                   nc.variables[keys[1]][:]).squeeze() * factor
            #if arr.ndim ==3:
            #	for i in np.arange(arr.shape[0]):
            #		arr[i] = arr[i]*area
            #elif arr.ndim ==2: arr = arr*area
            #else: assert 0
            return arr

        name = 'IntegratedPrimaryProduction_OSU'
        if annual:
            av[name]['modelFiles'] = listModelDataFiles(
                jobID, 'diad_T', paths.ModelFolder_pref, annual)
            av[name][
                'dataFile'] = paths.OSUDir + "/standard_VGPM.SeaWIFS.global.average.nc"
#		else:
#			print ""

        av[name]['modelcoords'] = medusaCoords
        av[name]['datacoords'] = osuCoords

        nc = dataset(av[name]['dataFile'], 'r')
        lats = nc.variables['latitude'][:]
        osuareas = np.zeros((1080, 2160))
        osuarea = (111100. / 6.)**2.  # area of a pixel at equator. in m2
        for a in np.arange(1080):
            osuareas[a] = np.ones(
                (2160, )) * osuarea * np.cos(np.deg2rad(lats[a]))

        def osuconvert(nc, keys):
            arr = nc.variables[keys[0]][:, :, :]
            #tlen = arr.shape[0]

            #arr  = arr.sum(0)/tlen * 365.	/ 1000. /     1E15
            #if arr.ndim ==3:
            #	for i in np.arange(arr.shape[0]):
            #		arr[i] = arr[i]*osuarea
            #elif arr.ndim ==2: arr = arr*osuarea
            #else: assert 0
            return arr

        av[name]['modeldetails'] = {
            'name': name,
            'vars': ['PRN', 'PRD'],
            'convert': medusadepthInt,
            'units': 'mgC/m^2/day'
        }
        av[name]['datadetails'] = {
            'name': name,
            'vars': [
                'NPP',
            ],
            'convert': osuconvert,
            'units': 'mgC/m^2/day'
        }

        av[name]['layers'] = [
            'Surface',
        ]
        av[name]['regions'] = regionList
        av[name]['metrics'] = metricList
        av[name]['datasource'] = 'OSU'
        av[name]['model'] = 'MEDUSA'
        av[name]['modelgrid'] = 'eORCA1'
        av[name]['gridFile'] = paths.orcaGridfn
        av[name]['Dimensions'] = 2

    #####
    # Total
    if 'IntPP_OSU' in analysisKeys:
        noOSU = True
        nc = dataset(paths.orcaGridfn, 'r')
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
            'name': 'IntPP_OSU',
            'vars': ['PRN', 'PRD'],
            'convert': medusadepthInt,
            'units': 'Gt/yr'
        }
        if noOSU:
            av[name]['datadetails'] = {'name': '', 'units': ''}

        else:
            nc = dataset(av[name]['dataFile'], 'r')
            lats = nc.variables['latitude'][:]
            osuareas = np.zeros((1080, 2160))
            osuarea = (111100. / 6.)**2.  # area of a pixel at equator. in m2
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
                'name': 'IntPP_OSU',
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
        av[name]['gridFile'] = paths.orcaGridfn
        av[name]['Dimensions'] = 1

    if 'GlobalExportRatio' in analysisKeys:

        def calcExportRatio(nc, keys):
            a = (nc.variables['SDT__100'][:] +
                 nc.variables['FDT__100'][:]).sum() / (
                     nc.variables['PRD'][:] + nc.variables['PRN'][:]).sum()
            #a = np.ma.masked_where(a>1.01, a)
            return a

        name = 'ExportRatio'
        av[name]['modelFiles'] = listModelDataFiles(jobID, 'diad_T',
                                                    paths.ModelFolder_pref,
                                                    annual)

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
        av[name]['gridFile'] = paths.orcaGridfn
        av[name]['Dimensions'] = 1

    if 'LocalExportRatio' in analysisKeys:

        def calcExportRatio(nc, keys):
            a = (nc.variables['SDT__100'][:] + nc.variables['FDT__100'][:]) / (
                nc.variables['PRD'][:] + nc.variables['PRN'][:])
            a = np.ma.masked_where(a > 1.01, a)
            return a

        name = 'LocalExportRatio'
        av[name]['modelFiles'] = listModelDataFiles(jobID, 'diad_T',
                                                    paths.ModelFolder_pref,
                                                    annual)

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
        ]  #'100m','200m','Surface - 1000m','Surface - 300m',]#'depthint']
        av[name]['regions'] = regionList
        av[name]['metrics'] = metricList
        av[name]['datasource'] = ''
        av[name]['model'] = 'MEDUSA'
        av[name]['modelgrid'] = 'eORCA1'
        av[name]['gridFile'] = paths.orcaGridfn
        av[name]['Dimensions'] = 2

    if 'Iron' in analysisKeys:

        name = 'Iron'
        av[name]['modelFiles'] = listModelDataFiles(jobID, 'ptrc_T',
                                                    paths.ModelFolder_pref,
                                                    annual)

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
        av[name]['gridFile'] = paths.orcaGridfn
        av[name]['Dimensions'] = 3

    if 'GlobalMeanTemperature' in analysisKeys:
        name = 'GlobalMeanTemperature'
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

        #av[name]['modelFiles']  = listModelDataFiles(jobID, 'grid_T', paths.ModelFolder_pref, annual)
        av[name]['dataFile'] = ''
        av[name]['modelcoords'] = medusaCoords
        av[name]['datacoords'] = woaCoords

        nc = dataset(paths.orcaGridfn, 'r')
        try:
            pvol = nc.variables['pvol'][:]
            gmttmask = nc.variables['tmask'][:]
        except:
            gmttmask = nc.variables['tmask'][:]
            area = nc.variables['e2t'][:] * nc.variables['e1t'][:]
            pvol = nc.variables['e3t'][:] * area
            pvol = np.ma.masked_where(gmttmask == 0, pvol)
        nc.close()

        def sumMeanLandMask(nc, keys):
            #### works like no change, but applies a mask.
            temp = np.ma.array(nc.variables[keys[0]][:].squeeze())
            temp = np.ma.masked_where((gmttmask == 0) + (temp.mask), temp)
            #try:
            vol = np.ma.masked_where(
                temp.mask,
                nc('thkcello')[:].squeeze() *
                nc('area')[:])  # preferentially use in file volume.
            #except: vol = np.ma.masked_where(temp.mask, pvol)
            return (temp * vol).sum() / (vol.sum())

        av[name]['modeldetails'] = {
            'name': name,
            'vars': [
                ukesmkeys['temp3d'],
            ],
            'convert': sumMeanLandMask,
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

    if 'GlobalMeanTemperature_700' in analysisKeys:
        name = 'GlobalMeanTemperature_700'
        av[name]['modelFiles'] = listModelDataFiles(jobID, 'grid_T',
                                                    paths.ModelFolder_pref,
                                                    annual)

        #av[name]['modelFiles']  = listModelDataFiles(jobID, 'grid_T', paths.ModelFolder_pref, annual)
        av[name]['dataFile'] = ''
        av[name]['modelcoords'] = medusaCoords
        av[name]['datacoords'] = woaCoords

        nc = dataset(paths.orcaGridfn, 'r')
        try:
            pvol = nc.variables['pvol'][:]
            gmttmask = nc.variables['tmask'][:]
        except:
            gmttmask = nc.variables['tmask'][:]
            area = nc.variables['e2t'][:] * nc.variables['e1t'][:]
            pvol = nc.variables['e3t'][:] * area
            pvol = np.ma.masked_where(gmttmask == 0, pvol)
        nc.close()

        def sumMeanLandMask(nc, keys):
            #### works like no change, but applies a mask.
            temp = np.ma.array(nc.variables[keys[0]][:].squeeze())
            temp = np.ma.masked_where((gmttmask == 0) + (temp.mask), temp)
            temp.mask[43:, :, :] = True
            #epth = nc.variables[medusaCoords['z']]
            vol = np.ma.masked_where(
                temp.mask,
                nc('thkcello')[:].squeeze() *
                nc('area')[:])  # preferentially use in file volume.
            return (temp * vol).sum() / (vol.sum())

        av[name]['modeldetails'] = {
            'name': name,
            'vars': [
                ukesmkeys['temp3d'],
            ],
            'convert': sumMeanLandMask,
            'units': 'degrees C'
        }
        av[name]['datadetails'] = {'name': '', 'units': ''}
        #av[name]['datadetails']        = {'name': name, 'vars':['t_an',], 'convert': ukp.NoChange,'units':'degrees C'}

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

    if 'GlobalMeanTemperature_2000' in analysisKeys:
        name = 'GlobalMeanTemperature_2000'
        av[name]['modelFiles'] = listModelDataFiles(jobID, 'grid_T',
                                                    paths.ModelFolder_pref,
                                                    annual)

        #av[name]['modelFiles']  = listModelDataFiles(jobID, 'grid_T', paths.ModelFolder_pref, annual)
        av[name]['dataFile'] = ''
        av[name]['modelcoords'] = medusaCoords
        av[name]['datacoords'] = woaCoords

        nc = dataset(paths.orcaGridfn, 'r')
        try:
            pvol = nc.variables['pvol'][:]
            gmttmask = nc.variables['tmask'][:]
        except:
            gmttmask = nc.variables['tmask'][:]
            area = nc.variables['e2t'][:] * nc.variables['e1t'][:]
            pvol = nc.variables['e3t'][:] * area
            pvol = np.ma.masked_where(gmttmask == 0, pvol)
        nc.close()

        def sumMeanLandMask(nc, keys):
            #### works like no change, but applies a mask.
            temp = np.ma.array(nc.variables[keys[0]][:].squeeze())
            temp = np.ma.masked_where((gmttmask == 0) + (temp.mask), temp)
            temp.mask[54:, :, :] = True
            #epth = nc.variables[medusaCoords['z']]
            vol = np.ma.masked_where(
                temp.mask,
                nc('thkcello')[:].squeeze() *
                nc('area')[:])  # preferentially use in file volume.
            return (temp * vol).sum() / (vol.sum())

        av[name]['modeldetails'] = {
            'name': name,
            'vars': [
                ukesmkeys['temp3d'],
            ],
            'convert': sumMeanLandMask,
            'units': 'degrees C'
        }
        av[name]['datadetails'] = {'name': '', 'units': ''}
        #av[name]['datadetails']        = {'name': name, 'vars':['t_an',], 'convert': ukp.NoChange,'units':'degrees C'}

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
        av[name]['modelFiles'] = listModelDataFiles(jobID, 'grid_T',
                                                    paths.ModelFolder_pref,
                                                    annual)
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
            temp = np.ma.array(nc.variables[keys[0]][:].squeeze())
            temp = np.ma.masked_where((gmttmask == 0) + (temp.mask), temp)

            #try:
            vol = np.ma.masked_where(
                temp.mask,
                nc('thkcello')[:].squeeze() *
                nc('area')[:])  # preferentially use in file volume.
            #except: vol = np.ma.masked_where(temp.mask, pvol)

            #vol = np.ma.masked_where(temp.mask, pvol)
            return ((temp * vol).sum(0) / vol.sum(0))  #*(area/area.sum())
            #return (((temp*vol).sum(0)/(vol.sum(0))) * (vol.sum(0)/vol.sum()))#.sum()

        av[name]['modeldetails'] = {
            'name': name,
            'vars': [
                ukesmkeys['temp3d'],
            ],
            'convert': sumMeanLandMask,
            'units': 'degrees C'
        }
        av[name]['datadetails'] = {'name': '', 'units': ''}
        #av[name]['datadetails']  	= {'name': name, 'vars':['t_an',], 'convert': ukp.NoChange,'units':'degrees C'}

        vmtregions = [
            'Global', 'ignoreInlandSeas', 'Equator10', 'AtlanticSOcean',
            'SouthernOcean', 'ArcticOcean', 'Remainder',
            'NorthernSubpolarAtlantic', 'NorthernSubpolarPacific', 'WeddelSea'
        ]
        vmtregions.extend(PierceRegions)

        av[name]['layers'] = [
            'layerless',
        ]
        av[name]['regions'] = vmtregions
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
        files = listModelDataFiles(jobID, 'scalar', paths.ModelFolder_pref,
                                   annual)
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
        files = listModelDataFiles(jobID, 'scalar', paths.ModelFolder_pref,
                                   annual)
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
        files = listModelDataFiles(jobID, 'scalar', paths.ModelFolder_pref,
                                   annual)
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
            #                	print nc('time_centered'), thetaoga ,scvoltot , rau0 , rcp, thetaoga * scvoltot * rau0 * rcp
            #			assert 0
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
                'convert': scalarFunction,
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

    if 'HeatFlux' in analysisKeys:
        name = 'HeatFlux'
        av[name]['modelFiles'] = listModelDataFiles(jobID, 'grid_T',
                                                    paths.ModelFolder_pref,
                                                    annual)
        av[name]['dataFile'] = ''
        av[name]['modelcoords'] = medusaCoords
        #av[name]['datacoords']          = takahashiCoords
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
            if 'area' in list(nc.variables.keys()):
                area = nc.variables['area'][:]
            else:
                area = ncarea
            flux = np.ma.array(nc.variables[keys[0]][:].squeeze()) * ncarea
            flux = np.ma.masked_where((surfmask == 0) + (flux.mask), flux)
            return flux.sum() * 1e-12

        av[name]['modelcoords'] = medusaCoords
        av[name]['modelFiles'] = listModelDataFiles(jobID, 'grid_T',
                                                    paths.ModelFolder_pref,
                                                    annual)
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

    if 'ERSST' in analysisKeys:
        name = 'ERSST'
        av[name]['modelFiles'] = [
            "/group_workspaces/jasmin4/esmeval/example_data/bgc/ERSST.v4/sst.mnmean.v4.nc",
        ]
        av[name]['dataFile'] = ''

        ERSSTCoords = {
            't': 'time',
            'z': '',
            'lat': 'lat',
            'lon': 'lon',
            'cal': 'standard',
            'tdict': ['ZeroToZero']
        }

        av[name]['modelcoords'] = ERSSTCoords
        av[name]['datacoords'] = woaCoords

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
            'gridFile'] = '/group_workspaces/jasmin4/esmeval/example_data/bgc/ERSST.v4/ERSST_sst_grid.nc'
        av[name]['Dimensions'] = 2

    if 'GlobalMeanSalinity' in analysisKeys:
        name = 'GlobalMeanSalinity'
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

        #av[name]['modelFiles']  = listModelDataFiles(jobID, 'grid_T', paths.ModelFolder_pref, annual)
        av[name]['dataFile'] = ''

        av[name]['modelcoords'] = medusaCoords
        av[name]['datacoords'] = woaCoords

        nc = dataset(paths.orcaGridfn, 'r')
        try:
            pvol = nc.variables['pvol'][:]
            gmstmask = nc.variables['tmask'][:]
        except:
            gmstmask = nc.variables['tmask'][:]
            area = nc.variables['e2t'][:] * nc.variables['e1t'][:]
            pvol = nc.variables['e3t'][:] * area
            pvol = np.ma.masked_where(gmstmask == 0, pvol)
        nc.close()

        def sumMeanLandMask(nc, keys):
            #### works like no change, but applies a mask.
            sal = np.ma.array(nc.variables[keys[0]][:].squeeze())
            sal = np.ma.masked_where((gmstmask == 0) + (sal.mask), sal)
            try:
                vol = np.ma.masked_where(
                    sal.mask,
                    nc('thkcello')[:].squeeze() *
                    nc('area')[:])  # preferentially use in file volume.
            except:
                vol = np.ma.masked_where(sal.mask, pvol)
            return (sal * vol).sum() / (vol.sum())

        av[name]['modeldetails'] = {
            'name': name,
            'vars': [
                ukesmkeys['sal3d'],
            ],
            'convert': sumMeanLandMask,
            'units': 'PSU'
        }
        av[name]['datadetails'] = {'name': '', 'units': ''}
        #av[name]['datadetails']  	globalVolumeMean= {'name': name, 'vars':['t_an',], 'convert': ukp.NoChange,'units':'degrees C'}

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
        av[name]['modelFiles'] = listModelDataFiles(jobID, 'grid_T',
                                                    paths.ModelFolder_pref,
                                                    annual)
        av[name]['dataFile'] = ''

        av[name]['modelcoords'] = medusaCoords
        av[name]['datacoords'] = woaCoords

        nc = dataset(paths.orcaGridfn, 'r')
        icetmask = nc.variables['tmask'][:]
        area_full = nc.variables['e2t'][:] * nc.variables['e1t'][:]
        nc.close()

        def calcIcelessMeanSST(nc, keys):
            #### works like no change, but applies a mask.
            icecov = nc.variables['soicecov'][:].squeeze()
            sst = np.ma.array(nc.variables[ukesmkeys['temp3d']][:,
                                                                0, ].squeeze())
            sst = np.ma.masked_where(
                (icetmask[0] == 0) + (icecov > 0.15) + sst.mask, sst)
            area = np.ma.masked_where(sst.mask, area_full)
            val = (sst * area).sum() / (area.sum())
            print("calcIcelessMeanSST", sst.shape, area.shape, val)
            return val

        av[name]['modeldetails'] = {
            'name': name,
            'vars': [
                'soicecov',
                ukesmkeys['temp3d'],
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

        if 'quickSST' in analysisKeys:
            name = 'quickSST'
            av[name]['modelFiles'] = listModelDataFiles(
                jobID, 'grid_T', paths.ModelFolder_pref, annual)

            nc = dataset(paths.orcaGridfn, 'r')
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
                    ukesmkeys['temp3d'],
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
            av[name]['gridFile'] = paths.orcaGridfn
            av[name]['Dimensions'] = 1

    if 'Temperature' in analysisKeys:
        name = 'Temperature'

        av[name]['modelFiles'] = listModelDataFiles(jobID, 'grid_T',
                                                    paths.ModelFolder_pref,
                                                    annual)
        if annual:
            #av[name]['modelFiles']  	= sorted(glob(paths.ModelFolder_pref+jobID+"/"+jobID+"o_1y_*_grid_T.nc"))
            av[name]['dataFile'] = WOAFolder + 'woa13_decav_t00_01v2.nc'
        else:
            #av[name]['modelFiles']  	= sorted(glob(paths.ModelFolder_pref+jobID+"/"+jobID+"o_1m_*_grid_T.nc"))
            av[name]['dataFile'] = WOAFolder + 'temperature_monthly_1deg.nc'
        av[name]['modelcoords'] = medusaCoords
        av[name]['datacoords'] = woaCoords

        av[name]['modeldetails'] = {
            'name': name,
            'vars': [
                ukesmkeys['temp3d'],
            ],
            'convert': applyLandMask,
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

        tregions = [
            'Global', 'ignoreInlandSeas', 'Equator10', 'AtlanticSOcean',
            'SouthernOcean', 'ArcticOcean', 'Remainder',
            'NorthernSubpolarAtlantic', 'NorthernSubpolarPacific', 'WeddelSea',
            'Cornwall'
        ]
        tregions.extend(PierceRegions)
        #tregions = ['Global',]
        if analysisSuite.lower() == 'spinup':
            av[name]['layers'] = layerList
            av[name]['regions'] = regionList
        else:
            av[name]['layers'] = layerList
            av[name]['regions'] = tregions
        av[name]['metrics'] = metricList

        #try:
        #	if analysisSuite.lower() in ['debug',]:
        #               av[name]['layers']              =  ['Surface',]
        #	        av[name]['regions']             =  ['Global',]
        #except:pass

        av[name]['datasource'] = 'WOA'
        av[name]['model'] = 'NEMO'

        av[name]['modelgrid'] = 'eORCA1'
        av[name]['gridFile'] = paths.orcaGridfn
        av[name]['Dimensions'] = 3

    if 'Salinity' in analysisKeys:
        name = 'Salinity'
        av[name]['modelFiles'] = listModelDataFiles(jobID, 'grid_T',
                                                    paths.ModelFolder_pref,
                                                    annual)
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
                ukesmkeys['sal3d'],
            ],
            'convert': applyLandMask,
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

        salregions = [
            'Global', 'ignoreInlandSeas', 'Equator10', 'AtlanticSOcean',
            'SouthernOcean', 'ArcticOcean', 'Remainder',
            'NorthernSubpolarAtlantic', 'NorthernSubpolarPacific', 'WeddelSea'
        ]
        salregions.extend(PierceRegions)
        if analysisSuite.lower() in ['spinup', 'salinity']:
            av[name]['layers'] = layerList
            av[name]['regions'] = regionList
        else:
            av[name]['layers'] = layerList
            av[name]['regions'] = salregions
        av[name]['metrics'] = metricList

        av[name]['datasource'] = 'WOA'
        av[name]['model'] = 'NEMO'

        av[name]['modelgrid'] = 'eORCA1'
        av[name]['gridFile'] = paths.orcaGridfn
        av[name]['Dimensions'] = 3

    if 'ZonalCurrent' in analysisKeys:
        name = 'ZonalCurrent'
        av[name]['modelFiles'] = listModelDataFiles(jobID, 'grid_U',
                                                    paths.ModelFolder_pref,
                                                    annual)
        if annual:
            av[name]['dataFile'] = paths.GODASFolder + 'ucur.clim.nc'

        av[name]['modelcoords'] = medusaUCoords
        av[name]['datacoords'] = godasCoords
        av[name]['modeldetails'] = {
            'name': name,
            'vars': [
                ukesmkeys['u3d'],
            ],
            'convert': applyLandMask1e3,
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
        av[name]['gridFile'] = os.path.join(paths.bgcval2_repo, 'bgcval2/data/eORCA1_gridU_mesh.nc')
        av[name]['Dimensions'] = 3

    if 'MeridionalCurrent' in analysisKeys:
        name = 'MeridionalCurrent'
        av[name]['modelFiles'] = listModelDataFiles(jobID, 'grid_V',
                                                    paths.ModelFolder_pref,
                                                    annual)
        if annual:
            av[name]['dataFile'] = paths.GODASFolder + 'vcur.clim.nc'

        av[name]['modelcoords'] = medusaVCoords
        av[name]['datacoords'] = godasCoords

        av[name]['modeldetails'] = {
            'name': name,
            'vars': [
                ukesmkeys['v3d'],
            ],
            'convert': applyLandMask1e3,
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
        av[name]['gridFile'] = os.path.join(paths.bgcval2_repo, 'bgcval2/data/eORCA1_gridV_mesh.nc')
        av[name]['Dimensions'] = 3

    if 'VerticalCurrent' in analysisKeys:
        name = 'VerticalCurrent'
        av[name]['modelFiles'] = listModelDataFiles(jobID, 'grid_W',
                                                    paths.ModelFolder_pref,
                                                    annual)
        if annual:
            av[name]['dataFile'] = paths.GODASFolder + 'dzdt.clim.nc'

        av[name]['modelcoords'] = medusaWCoords
        av[name]['datacoords'] = godasCoords

        def applyLandMask1e6(nc, keys):
            return applyLandMask(nc, keys) * 1000000.

        av[name]['modeldetails'] = {
            'name': name,
            'vars': [
                ukesmkeys['w3d'],
            ],
            'convert': applyLandMask1e6,
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

        vregions = regionList
        #                vregions.extend(['NordicSea', 'LabradorSea', 'NorwegianSea'])

        av[name]['layers'] = layerList
        av[name]['regions'] = vregions
        av[name]['metrics'] = metricList

        av[name]['datasource'] = 'GODAS'
        av[name]['model'] = 'NEMO'

        av[name]['modelgrid'] = 'eORCA1'
        av[name]['gridFile'] = os.path.join(paths.bgcval2_repo, 'bgcval2/data/eORCA1_gridW_mesh.nc')
        av[name]['Dimensions'] = 3

    if 'WindStress' in analysisKeys:
        name = 'WindStress'
        av[name]['modelFiles'] = listModelDataFiles(jobID, 'grid_U',
                                                    paths.ModelFolder_pref,
                                                    annual)
        av[name]['dataFile'] = ''

        #paths.GODASFolder+'ucur.clim.nc'

        def calcwind(nc, keys):
            taux = applySurfaceMask(nc, [
                'sozotaux',
            ])

            ncpath = nc.filename
            newpath = ncpath.replace('grid_U', 'grid_V')

            nc2 = dataset(newpath, 'r')
            print("Loaded", newpath)
            tauy = applySurfaceMask(nc2, [
                'sometauy',
            ])

            return np.ma.sqrt(taux * taux + tauy * tauy)

        av[name]['modelcoords'] = medusaUCoords
        av[name]['datacoords'] = godasCoords

        av[name]['modeldetails'] = {
            'name': name,
            'vars': ['sozotaux', 'sometauy'],
            'convert': calcwind,
            'units': 'N/m2'
        }
        av[name]['datadetails'] = {'name': '', 'units': ''}

        av[name]['layers'] = [
            'layerless',
        ]
        av[name]['regions'] = regionList
        av[name]['metrics'] = metricList

        av[name]['datasource'] = ''
        av[name]['model'] = 'NEMO'

        av[name]['modelgrid'] = 'eORCA1'
        av[name]['gridFile'] = os.path.join(paths.bgcval2_repo, 'bgcval2./data/eORCA1_gridU_mesh.nc')
        av[name]['Dimensions'] = 2

    #####
    # North Atlantic Salinity
    #sowaflup = "Net Upward Water Flux" ;
    #sohefldo = "Net Downward Heat Flux" ;
    #sofmflup = "Water flux due to freezing/melting" ;
    #sosfldow = "Downward salt flux" ;

    naskeys = [
        'sowaflup',
        'sohefldo',
        'sofmflup',
        'sosfldow',
        'soicecov',
        'sossheig',
    ]
    if len(set(naskeys).intersection(set(analysisKeys))):
        for name in naskeys:
            if name not in analysisKeys: continue

            #nc = dataset(paths.orcaGridfn,'r')
            #area = nc.variables['e2t'][:] * nc.variables['e1t'][:]
            #tmask = nc.variables['tmask'][0,:,:]
            #lat = nc.variables['nav_lat'][:,:]
            #nc.close()

            nas_files = listModelDataFiles(jobID, 'grid_T',
                                           paths.ModelFolder_pref, annual)
            nc = dataset(nas_files[0], 'r')
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
                'convert': applySurfaceMask,
                'units': nasUnits[name][:]
            }

            #av[name]['regions'] 		=  ['NordicSea', 'LabradorSea', 'NorwegianSea','Global', ]
            av[name]['regions'] = [
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

    naskeys = [
        'max_soshfldo',
    ]  #'sohefldo','sofmflup','sosfldow','soicecov']
    if len(set(naskeys).intersection(set(analysisKeys))):
        for name in naskeys:
            if name not in analysisKeys: continue

            cutname = name[:].replace('max_', '')
            #nc = dataset(paths.orcaGridfn,'r')
            #area = nc.variables['e2t'][:] * nc.variables['e1t'][:]
            #tmask = nc.variables['tmask'][0,:,:]
            #lat = nc.variables['nav_lat'][:,:]
            #nc.close()

            nas_files = listModelDataFiles(jobID, 'grid_T',
                                           paths.ModelFolder_pref, annual)
            nc = dataset(nas_files[0], 'r')
            if cutname not in list(nc.variables.keys()):
                print("analysis_timeseries.py:\tWARNING: ", cutname,
                       "is not in the model file.")
                continue
            av[name]['modelFiles'] = nas_files
            av[name]['dataFile'] = ''

            av[name]['modelcoords'] = medusaCoords
            av[name]['datacoords'] = medusaCoords

            maxUnits = {
                'max_soshfldo': "W/m2",
            }

            def getMax(nc, keys):
                return applySurfaceMask(nc, keys).max()

            av[name]['modeldetails'] = {
                'name': name[:],
                'vars': [
                    cutname,
                ],
                'convert': getMax,
                'units': maxUnits[name][:]
            }

            av[name]['regions'] = [
                'Global',
            ]  # 'LabradorSea', 'NorwegianSea', ]

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
            av[name]['gridFile'] = paths.orcaGridfn
            av[name]['Dimensions'] = 1

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
        av[name]['modelFiles'] = listModelDataFiles(jobID, 'grid_T',
                                                    paths.ModelFolder_pref,
                                                    annual)

        av[name]['dataFile'] = ""
        av[name]['modelcoords'] = medusaCoords
        av[name]['datacoords'] = maredatCoords
        av[name]['modeldetails'] = {
            'name': name,
            'vars': ['ficeberg', 'friver', 'fsitherm', 'pr', 'prsn', 'evs'],
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

    if 'MLD' in analysisKeys:

        def mldapplymask(nc, keys):
            mld = np.ma.array(nc.variables[keys[0]][:])
            return np.ma.masked_where(
                (nc.variables[keys[1]][:] == 0.) + mld.mask + (mld == 1.E9),
                mld)
#eturn np.ma.masked_where((np.tile(nc.variables[keys[1]][:],(12,1,1))==0.)+mld.mask+(mld==1.E9),mld)

#nc = dataset(paths.orcaGridfn,'r')
#depth = nc.variables['nav_lev'][:]#
#nc.close()
#def calcMLD(nc,keys):
#	temp = nc.variables[keys[0]][:,:,:,:]
#	f_out = interp1d(depth[7:9],temp[7:9], axis=1)
#	tcrit = 0.2
#	t10m =  f_out(10.)
#	t10m = np.ma.masked_where(t10m>1E20, t10m) - tcrit
#	# linear regression to extrapolate below this level to find the first?
#	f_out = interp1d(temp, depth, axis=1)
#	#t_out = f_out(newdepth)

        name = 'MLD'
        av[name]['modelFiles'] = listModelDataFiles(jobID, 'grid_T',
                                                    paths.ModelFolder_pref,
                                                    annual)
        #av[name]['modelFiles']  	= sorted(glob(paths.ModelFolder_pref+jobID+"/"+jobID+"o_1y_*_grid_T.nc"))
        av[name][
            'dataFile'] = paths.MLDFolder + "mld_DT02_c1m_reg2.0-annual.nc"  #mld_DT02_c1m_reg2.0.nc"
        #MLD_DT02 = depth where (T = T_10m +/- 0.2 degC)

        av[name]['modelcoords'] = medusaCoords
        av[name]['datacoords'] = mldCoords

        av[name]['modeldetails'] = {
            'name': 'mld',
            'vars': [ ukesmkeys['MLD'],],
            'convert': applySurfaceMask,
            'units': 'm'
        }
        #av[name]['modeldetails'] 	= {'name': 'mld', 'vars':[ukesmkeys['temp3d'],],   'convert': calcMLD,'units':'m'}
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
        ]  #'Surface - 1000m','Surface - 300m',]#'depthint'

        mldregions = [
            'Global', 'ignoreInlandSeas', 'Equator10', 'AtlanticSOcean',
            'SouthernOcean', 'ArcticOcean', 'Remainder',
            'NorthernSubpolarAtlantic', 'NorthernSubpolarPacific', 'WeddelSea'
        ]
        mldregions.extend(PierceRegions)

        av[name]['regions'] = mldregions
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
                    mld = np.ma.masked_where((nc.variables[keys[1]][:] == 0.) +
                                             mld.mask + (mld == 1.E9), mld)
                    return mld

                def mldmonthlymask(nc, keys):
                    mld = np.ma.array(np.ma.abs(
                        nc.variables[keys[0]][:])).max(0)
                    mld = np.ma.masked_where(mld.mask + (mld.data > 1.E10),
                                             mld)
                    return mld

                def mldapplymask_min(nc, keys):
                    mld = np.ma.array(nc.variables[keys[0]][:]).min(0)
                    mld = np.ma.masked_where((nc.variables[keys[1]][:] == 0.) +
                                             mld.mask + (mld == 1.E9), mld)
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

        else:
            print("No monthly MLD files found")

    icekeys = [
        'NorthernTotalIceArea', 'SouthernTotalIceArea', 'WeddelTotalIceArea',
        'TotalIceArea', 'NorthernTotalIceExtent', 'WeddelIceExent',
        'SouthernTotalIceExtent', 'TotalIceExtent', 'NorthernMIZArea',
        'SouthernMIZArea', 'TotalMIZArea', 'NorthernMIZfraction',
        'SouthernMIZfraction', 'TotalMIZfraction'
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
                ice = nc.variables[keys[0]][:].squeeze()
                new_area = nc.variables['area'][:].squeeze()
                miz_area = np.ma.masked_where(
                    (ice < 0.15) + (ice > 0.8) + ice.mask, new_area)
                total_area = np.ma.masked_where((ice < 0.15) + ice.mask,
                                                new_area)
                return miz_area.sum() / total_area.sum()

            def calcMIZfractionN(nc, keys):  # North
                ice = nc.variables[keys[0]][:].squeeze()
                new_area = nc.variables['area'][:].squeeze()
                miz_area = np.ma.masked_where(
                    (lat < 0.) + (ice < 0.15) + (ice > 0.8) + ice.mask,
                    new_area)
                total_area = np.ma.masked_where(
                    (lat < 0.) + (ice < 0.15) + ice.mask, new_area)
                return miz_area.sum() / total_area.sum()

            def calcMIZfractionS(nc, keys):  # South
                ice = nc.variables[keys[0]][:].squeeze()
                new_area = nc.variables['area'][:].squeeze()
                miz_area = np.ma.masked_where(
                    (lat > 0.) + (ice < 0.15) + (ice > 0.8) + ice.mask,
                    new_area)
                total_area = np.ma.masked_where(
                    (lat > 0.) + (ice < 0.15) + ice.mask, new_area)
                return miz_area.sum() / total_area.sum()

            def calcTotalIceExtent(nc, keys):  #Global
                return np.ma.masked_where(
                    (tmask == 0) + (nc.variables[keys[0]][:].squeeze() < 0.15),
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

            weddelmask = (lat < -80.) + (lat > -65.) + (lon < -60.) + (lon >
                                                                       -20.)

            def calcTotalIceExtentWS(nc, keys):  # South
                return np.ma.masked_where(
                    (tmask == 0) +
                    (nc.variables[keys[0]][:].squeeze() < 0.15) + weddelmask,
                    area).sum() / 1E12

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

    if 'DrakePassageTransport' in analysisKeys:
        name = 'DrakePassageTransport'
        ####
        # Note that this will only work with the eORCAgrid.

        # coordinates of Drake Passage
        LON = 219
        LAT0 = 79
        LAT1 = 109

        nc = dataset(paths.orcaGridfn, 'r')
        e2u = nc.variables['e2u'][LAT0:LAT1, LON]
        umask = nc.variables['umask'][:, LAT0:LAT1, LON]
        nc.close()

        def drake(nc, keys):
            e3u = nc.variables[ukesmkeys['e3u']][0, :, LAT0:LAT1, LON]
            velo = nc.variables[ukesmkeys['u3d']][0, :, LAT0:LAT1, LON]
            return np.sum(velo * e3u * e2u * umask) * 1.e-6

        av[name]['modelFiles'] = listModelDataFiles(jobID, 'grid_U',
                                                    paths.ModelFolder_pref,
                                                    annual)
        av[name]['dataFile'] = ''

        av[name]['modelcoords'] = medusaCoords
        av[name]['datacoords'] = medusaCoords

        av[name]['modeldetails'] = {
            'name': name,
            'vars': [
                ukesmkeys['e3u'],
                ukesmkeys['u3d'],
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
        av[name]['gridFile'] = paths.orcaGridfn
        av[name]['Dimensions'] = 1

    if 'AMOC_26N' in analysisKeys or 'AMOC_32S' in analysisKeys or 'ADRC_26N' in analysisKeys or 'AMOC_26N_nomexico' in analysisKeys:
        # Note that this will only work with the eORCAgrid.
        latslice26N = slice(227, 228)
        latslice26Nnm = slice(228, 229)
        latslice32S = slice(137, 138)
        e3v, e1v, tmask, alttmask = {}, {}, {}, {}
        for name in ['AMOC_26N', 'AMOC_32S', 'ADRC_26N', 'AMOC_26N_nomexico']:
            if name not in analysisKeys: continue

            ####
            if name in ['AMOC_26N', 'ADRC_26N']: latslice = latslice26N
            if name == 'AMOC_32S': latslice = latslice32S
            if name in [
                    'AMOC_26N_nomexico',
            ]:
                latslice = latslice26Nnm
            # Load grid data
            nc = dataset(paths.orcaGridfn, 'r')
            e3v[name] = nc.variables['e3v'][:,
                                            latslice, :]  # z level height 3D
            e1v[name] = nc.variables['e1v'][latslice, :]  #
            tmask[name] = nc.variables['tmask'][:, latslice, :]
            nc.close()

            # load basin mask
            nc = dataset(os.path.join(paths.bgcval2_repo, 'bgcval2/data/basinlandmask_eORCA1.nc'), 'r')
            alttmask[name] = nc.variables['tmaskatl'][
                latslice, :]  # 2D Atlantic mask
            if name == [
                    'AMOC_26N_nomexico',
            ]:
                alttmask[name][228, 180:208] = 0.
            nc.close()

        def calc_amoc32S(nc, keys):
            name = 'AMOC_32S'
            zv = np.ma.array(
                nc.variables[ukesmkeys['v3d']][..., latslice32S, :])  # m/s
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

        def amoc26N_array(nc, keys, amocname='AMOC_26N'):
            zv = np.ma.array(
                nc.variables[ukesmkeys['v3d']][..., latslice26N, :])  # m/s
            atlmoc = np.array(np.zeros_like(zv[0, :, :, 0]))
            e2vshape = e3v[amocname].shape
            for la in range(e2vshape[1]):  #ji, y
                for lo in range(e2vshape[2]):  #jj , x,
                    if int(alttmask[amocname][la, lo]) == 0: continue
                    for z in range(e2vshape[0]):  # jk
                        if int(tmask[amocname][z, la, lo]) == 0: continue
                        if np.ma.is_masked(zv[0, z, la, lo]): continue
                        atlmoc[z, la] = atlmoc[z, la] - e1v[amocname][
                            la, lo] * e3v[amocname][z, la, lo] * zv[0, z, la,
                                                                    lo] / 1.E06

####
# Cumulative sum from the bottom up.
            for z in range(73, 1, -1):
                atlmoc[z, :] = atlmoc[z + 1, :] + atlmoc[z, :]
            #return np.ma.max(atlmoc)
            return atlmoc

        def calc_amoc26N(nc, keys):
            return np.ma.max(amoc26N_array(nc, keys, amocname='AMOC_26N'))

        def calc_amoc26Nnm(nc, keys):
            return np.ma.max(
                amoc26N_array(nc, keys, amocname='AMOC_26N_nomexico'))

        def calc_min_amoc26N(nc, keys):
            return np.ma.min(amoc26N_array(nc, keys, amocname='ADRC_26N'))

        for name in ['AMOC_26N', 'AMOC_32S', 'ADRC_26N', 'AMOC_26N_nomexico']:
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
                        ukesmkeys['v3d'],
                    ],
                    'convert': calc_amoc26N,
                    'units': 'Sv'
                }
            if name in [
                    'AMOC_26N_nomexico',
            ]:
                av[name]['modeldetails'] = {
                    'name': name,
                    'vars': [
                        ukesmkeys['v3d'],
                    ],
                    'convert': calc_amoc26Nnm,
                    'units': 'Sv'
                }
            if name == 'ADRC_26N':
                av[name]['modeldetails'] = {
                    'name': name,
                    'vars': [
                        ukesmkeys['v3d'],
                    ],
                    'convert': calc_min_amoc26N,
                    'units': 'Sv'
                }
            if name == 'AMOC_32S':
                av[name]['modeldetails'] = {
                    'name': name,
                    'vars': [
                        ukesmkeys['v3d'],
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

    if 'DMS_ARAN' in analysisKeys:
        name = 'DMS'
        av[name]['modelFiles'] = listModelDataFiles(jobID, 'diad_T',
                                                    paths.ModelFolder_pref,
                                                    annual)[:]
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
            'convert': ukp.NoChange,
            'units': 'nmol/L'
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

    if 'DMS_ANDR' in analysisKeys:
        name = 'DMS'
        av[name]['modelFiles'] = listModelDataFiles(jobID, 'diad_T',
                                                    paths.ModelFolder_pref,
                                                    annual)[:]
        if annual:
            av[name]['dataFile'] = paths.DMSDir + 'DMSclim_mean.nc'
        else:
            av[name]['dataFile'] = ''

        av[name]['modelcoords'] = medusaCoords
        av[name]['datacoords'] = dmsCoords

        av[name]['modeldetails'] = {
            'name': name,
            'vars': [
                'DMS_ANDR',
            ],
            'convert': ukp.NoChange,
            'units': 'nmol/L'
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

    if 'Dust' in analysisKeys:
        name = 'Dust'
        av[name]['modelFiles'] = listModelDataFiles(jobID, 'diad_T',
                                                    paths.ModelFolder_pref,
                                                    annual)[:]
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
            return dust * 0.035 * 1.e6 * 0.00532 * (24. * 60. * 60.) / 55.845

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

    if 'TotalDust' in analysisKeys:
        name = 'TotalDust'
        av[name]['modelFiles'] = listModelDataFiles(jobID, 'diad_T',
                                                    paths.ModelFolder_pref,
                                                    annual)[:]
        av[name]['dataFile'] = paths.Dustdir + 'mahowald.orca100_annual.nc'

        av[name]['modelcoords'] = medusaCoords
        av[name]['datacoords'] = medusaCoords

        nc = dataset(paths.orcaGridfn, 'r')
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

    if 'TotalDust_nomask' in analysisKeys:
        name = 'TotalDust_nomask'
        av[name]['modelFiles'] = listModelDataFiles(jobID, 'diad_T',
                                                    paths.ModelFolder_pref,
                                                    annual)[:]
        av[name]['dataFile'] = paths.Dustdir + 'mahowald.orca100_annual.nc'

        av[name]['modelcoords'] = medusaCoords
        av[name]['datacoords'] = medusaCoords

        nc = dataset(paths.orcaGridfn, 'r')
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
            #dust[194:256,295:348] = 0.
            #dust[194:208,285:295] = 0.
            #dust[188:216,290:304] = 0.
            return (masked_area *
                    dust).sum() * 0.035 * 1.e6 * 1.e-12 * 0.00532 * (
                        24. * 60. * 60. * 365.25) / 55.845

        def modeldustsum(nc, keys):
            dust = nc.variables[keys[0]][:]
            #dust[:,194:256,295:348] = 0.
            #dust[:,194:208,285:295] = 0.
            #dust[:,188:216,290:304] = 0.
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
            clean=clean,
        )

        #####
        # Profile plots
        if av[name]['Dimensions'] == 3 and name not in ['Iron', 'Fe']:
            continue
            profa = profileAnalysis(
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
                layers=list(
                    np.arange(102)
                ),  # 102 because that is the number of layers in WOA Oxygen
                regions=av[name]['regions'],
                metrics=[
                    'mean',
                ],
                workingDir=shelvedir,
                imageDir=imagedir,
                grid=av[name]['modelgrid'],
                gridFile=av[name]['gridFile'],
                clean=False,
            )
            #shelves[name] = profa.shelvefn
            #shelves_insitu[name] = profa.shelvefn_insitu

        #shelves[name] = tsa.shelvefn
        #shelves_insitu[name] = tsa.shelvefn_insitu


def singleTimeSeriesProfile(jobID, key):

    FullDepths = [
        'T',
        'S',
        'Chl_pig',
        'N',
        'Si',
        'O2',
        'Alk',
        'DIC',
        'Iron',
    ]
    if key in FullDepths:
        analysis_timeseries(
            jobID=jobID,
            analysisSuite=[
                key,
            ],
        )


def singleTimeSeries(
    jobID,
    key,
):
    #	try:
    analysis_timeseries(jobID=jobID,
                        analysisSuite=[
                            key,
                        ],
                        strictFileCheck=False)  #clean=1)


#	except:
#		print "Failed singleTimeSeries",(jobID,key)
#		print "Error: %s" % sys.exc_info()[0]


def get_args():
    """Parse command line arguments."""
    accepted_keys = ['kmf', 'physics','bgc', 'debug', 'spinup', 'salinity', 'level1', 'level3', ]
    
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-c',
                        '--config-file',
                        default=os.path.join(os.getcwd(),
                                             'bgcval2-config-user.yml'),
                        help='User configuration file',
                        required=False)
    parser.add_argument('-j',
                        '--jobID',
                        default=None,
                        help='Job ID (automatically generated by cycl/rose suite for UKESM jobs).',
                        required=True)
    parser.add_argument('-k',
                        '--keys',
                        default=['kmf', 'level1',],
                        nargs='+', 
                        type=str,
                        help=''.join(['Runtime keys - each key links to a pre-determined list of variables to analyse. ',
                                      'Keys are: ', ', '.join( accepted_keys)]),
                        required=False)

    args = parser.parse_args()

    return args


def main():
    from ._version import __version__
    print(f'BGCVal2: {__version__}')
    args = get_args()
    jobID = args.job_id
    keys = args.keys
    print('Running analysis_imerseries.\tjobID:', jobID, '\tkeys:', keys)
    assert 0

    accepted_keys = ['kmf', 'physics','bgc', 'debug', 'spinup', 'salinity', 'level1', 'level3', ]
    good_keys = True
    for key in keys:
        if key not in accepted_keys:
            print('Key Argument [',key,'] nor recognised. Accepted keys are:', accepted_keys)
            good_keys= False
    if good_keys: 
        sys.exit(0)

    """
    if 'debug' in keys:
        suite = 'debug'
        #elif 'all' in argv[1:]:	suite = 'all'
    elif 'spinup' in keys:
        suite = 'spinup'
    elif 'salinity' in keys:
        suite = 'salinity'
    elif 'level1' in keys:
        suite = 'level1'
    elif 'fast' in keys:
        suite = 'fast'
    elif 'level3' in keys:
        suite = 'level3'
    elif 'physics' in keys:
        suite = 'physics'
    elif 'bgc' in keys:
        suite = 'bgc'
    elif 'kmf' in keys or 'keymetricsfirst' in keys:
        suite = 'keymetricsfirst'
    else:
        suite = 'level1'
    """ 
    if args.config_file:
        config_user = os.path.join(os.getcwd(), args.config_file)
        print(f"analysis_timeseries: Using user config file {config_user}")
    else:
        config_user = os.path.join(os.getcwd(), "bgcval2-config-user.yml")
        print(f"analysis_timeseries: Using user default file {config_user}")
    if not os.path.isfile(config_user):
        print(f"analysis_timeseries: Could not find configuration file {config_user}."
              "Will proceed with defaults.")
        config_user = None

    for suite in keys:
        analysis_timeseries(
            jobID=jobID,
            analysisSuite=suite,
            config_user=config_user
        )  #clean=1)
    #if suite == 'all':
    #analysis_timeseries(jobID =jobID,analysisSuite='FullDepth', z_component = 'FullDepth',)#clean=1)


if __name__ == "__main__":
    main()
