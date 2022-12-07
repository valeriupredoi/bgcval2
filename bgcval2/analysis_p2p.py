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
import matplotlib as mpl
mpl.use('Agg')
import argparse
import os
from glob import glob
from netCDF4 import Dataset
import numpy as np
import sys

#####
#Specific local code:
from bgcval2._runtime_config import get_run_configuration
from bgcval2 import UKESMpython as ukp
#from bgcval2.p2p import makePatternStatsPlots, testsuite_p2p
from bgcval2.p2p.testsuite_p2p import testsuite_p2p
from bgcval2.p2p.summaryTargets import summaryTargets
#from bgcval2.p2p.patternAnalyses import InterAnnualPatterns, BGCvsPhysics
from bgcval2.bgcvaltools.pftnames import months
#from bgcval2.p2p.shelveToDictionary import shelveToDictionary
from bgcval2.analysis_timeseries import build_list_of_suite_keys, load_key_file

#####
# User defined set of paths pointing towards the datasets.
from bgcval2.Paths.paths import paths_setter


def analysis_p2p(
    jobID='u-ad980',
    years=['1997'],  #'2075','2076',
    modelGrid='eORCA1',
    annual=True,
    noPlots=False,
    suites='default',
    noTargets=True,
    config_user=None
):
    """

	"""
    # get runtime configuration
    if config_user:
        paths_dict, config_user = get_run_configuration(config_user)
    else:
        paths_dict, config_user = get_run_configuration("defaults")

    # filter paths dict into an object that's usable below
    paths = paths_setter(paths_dict)

    #####
    # Switches:
    # These are some booleans that allow us to choose which analysis to run.
    # This lets up give a list of keys one at a time, or in parrallel.
    #if type(suites) == type(['Its', 'A', 'list!']):
    if isinstance(suites, str):
        suites = suites.replace(',', ' ').replace('\'', '').replace('"', '')
        suites = suites.split(' ')

    analysisKeys = build_list_of_suite_keys(suites, debug=True)
    print('analysisKeys', analysisKeys)

    ModelFolder_pref = paths.ModelFolder_pref

    # NEW STYLE keys from file:
    av = ukp.AutoVivification()
    for key in analysisKeys:
        print('Loading:', key)
        av[key] = load_key_file(key, paths, jobID)

    print(analysisKeys, av)
    # def listModelDataFiles(jobID, filekey, datafolder, annual, yr):
    #     print("listing model data files:", jobID, filekey, datafolder, annual)
    #     if annual:
    #         keystr = datafolder + jobID + "/" + jobID + "o_1y_*1201[-_]" + yr + '????_' + filekey + ".nc"
    #         print("listModelDataFiles:", keystr)
    #         return sorted(glob(keystr))[0]
    #     else:
    #         return sorted(
    #             glob(datafolder + jobID + "/" + jobID + "o_1m_*" + yr +
    #                  "????_" + filekey + ".nc"))[-1]

#####
# Because we can never be sure someone won't randomly rename the
# time dimension without saying anything.
# if jobID in ['u-am515','u-am927','u-am064','u-an326',]:

    # print(jobID, 'grid_T', paths.ModelFolder_pref, annual)
    # tmpModelFiles = listModelDataFiles(jobID, 'grid_T', paths.ModelFolder_pref,
    #                                    annual, '*')
    #
    # try:
    #     tmpModelFiles = listModelDataFiles(jobID, 'grid_T',
    #                                        paths.ModelFolder_pref, annual, '*')
    # except:
    #     print(
    #         "No grid_T Model files available to figure out what naming convention is used."
    #     )
    #     tmpModelFiles = []
    # ukesmkeys = {}
    # if len(tmpModelFiles):
    #     print('test opening:', tmpModelFiles)
    #     nctmp = Dataset(tmpModelFiles, 'r')
    #     nctmpkeys = list(nctmp.variables.keys())
    #     nctmp.close()
    #     if 'votemper' in nctmpkeys:
    #         ukesmkeys = {}
    #         ukesmkeys['time'] = 'time_counter'
    #         ukesmkeys['temp3d'] = 'votemper'
    #         ukesmkeys['sst'] = ''
    #         ukesmkeys['sal3d'] = 'vosaline'
    #         ukesmkeys['sss'] = ''
    #         ukesmkeys['v3d'] = 'vomecrty'
    #         ukesmkeys['u3d'] = 'vozocrtx'
    #         ukesmkeys['e3u'] = 'e3u'
    #         ukesmkeys['w3d'] = 'vovecrtz'
    #     else:
    #         ukesmkeys['time'] = 'time_centered'
    #         ukesmkeys['temp3d'] = 'thetao'
    #         ukesmkeys['sst'] = 'tos'
    #         ukesmkeys['sal3d'] = 'so'
    #         ukesmkeys['sss'] = 'sos'
    #         ukesmkeys['v3d'] = 'vo'
    #         ukesmkeys['u3d'] = 'uo'
    #         ukesmkeys['e3u'] = 'thkcello'
    #         ukesmkeys['w3d'] = 'wo'
    # else:
    #     print("No grid_T files Found")
    # print('ukesmkeys[sal3d]:', ukesmkeys['sal3d'])

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

    # medusaCoords = {
    #     't': 'index_t',
    #     'z': 'deptht',
    #     'lat': 'nav_lat',
    #     'lon': 'nav_lon',
    #     'cal': '360_day',
    #     'tdict': ukp.tdicts['ZeroToZero']
    # }  # model doesn't need time dict.
    # medusaUCoords = {
    #     't': 'index_t',
    #     'z': 'depthu',
    #     'lat': 'nav_lat',
    #     'lon': 'nav_lon',
    #     'cal': '360_day',
    # }  # model doesn't need time dict.
    # medusaVCoords = {
    #     't': 'index_t',
    #     'z': 'depthv',
    #     'lat': 'nav_lat',
    #     'lon': 'nav_lon',
    #     'cal': '360_day',
    # }  # model doesn't need time dict.
    # medusaWCoords = {
    #     't': 'index_t',
    #     'z': 'depthw',
    #     'lat': 'nav_lat',
    #     'lon': 'nav_lon',
    #     'cal': '360_day',
    # }  # model doesn't need time dict.
    #
    # maredatCoords = {
    #     't': 'index_t',
    #     'z': 'DEPTH',
    #     'lat': 'LATITUDE',
    #     'lon': 'LONGITUDE',
    #     'cal': 'standard',
    #     'tdict': ukp.tdicts['ZeroToZero']
    # }
    # woaCoords = {
    #     't': 'index_t',
    #     'z': 'depth',
    #     'lat': 'lat',
    #     'lon': 'lon',
    #     'cal': 'standard',
    #     'tdict': ukp.tdicts['ZeroToZero']
    # }
    # cciCoords = {
    #     't': 'index_t',
    #     'z': 'index_z',
    #     'lat': 'lat',
    #     'lon': 'lon',
    #     'cal': 'standard',
    #     'tdict': ukp.tdicts['ZeroToZero']
    # }
    # glodapCoords = {
    #     't': 'index_t',
    #     'z': 'depth',
    #     'lat': 'latitude',
    #     'lon': 'longitude',
    #     'cal': 'standard',
    #     'tdict': ukp.tdicts['ZeroToZero']
    # }
    # osuCoords = {
    #     't': 'index_t',
    #     'z': 'index_z',
    #     'lat': 'latitude',
    #     'lon': 'longitude',
    #     'cal': 'standard',
    #     'tdict': ukp.tdicts['ZeroToZero']
    # }
    # glodapv2Coords = {
    #     't': 'index_t',
    #     'z': 'Pressure',
    #     'lat': 'lat',
    #     'lon': 'lon',
    #     'cal': '',
    #     'tdict': {
    #         0: 0,
    #     }
    # }
    # takahashiCoords = {
    #     't': 'index_t',
    #     'z': 'index_z',
    #     'lat': 'LAT',
    #     'lon': 'LON',
    #     'cal': 'standard',
    #     'tdict': ukp.tdicts['ZeroToZero']
    # }
    # godasCoords = {
    #     't': 'index_t',
    #     'z': 'level',
    #     'lat': 'lat',
    #     'lon': 'lon',
    #     'cal': 'standard',
    #     'tdict': ukp.tdicts['ZeroToZero']
    # }

    shelvesAV = []
    for year in years:
        #####
        # Location of model files.
        if annual:
            ModelFolder = ModelFolder_pref + jobID + "/"
        else:
            ModelFolder = ModelFolder_pref + str(int(year)) + '/'
        
        workingDir = ''.join([paths.p2p_ppDir, '/', jobID, '-', str(int(year))])
        imageFolder = ''.join([paths.imagedir, '/', jobID])

        print('P2P working dir:', workingDir)
        print('P2P image dir:', imageFolder)

        # Attempt to make the directories:
        workingDir = ukp.folder(workingDir)
        imageFolder = ukp.folder(imageFolder)

        shelvesAV.extend(
            testsuite_p2p(
                model='NEMO',
                jobID=jobID,
                year=year,
                av=av,
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


def get_args():
    """Parse command line arguments. """
    #accepted_keys = ['kmf', 'physics','bgc', 'debug', 'spinup', 'salinity', 'fast', 'level1', 'level3', 'nowmaps']

    accepted_keys = [os.path.splitext(os.path.basename(fn))[0] for fn in glob('key_lists/*.yml')]
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-j',
                        '--jobID',
                        nargs='+',
                        type=str,
                        default=None,
                        help='One or more UKESM Job IDs (automatically generated by the cylc/rose suite).',
                        required=True,
                        )
    parser.add_argument('-y',
                        '--years',
                        nargs='+',
                        type=int,
                        default=None,
                        help='One or more years to focus in on.',
                        required=True,
                        )

    parser.add_argument('-k',
                        '--keys',
                        default=['kmf', 'level1',],
                        nargs='+',
                        type=str,
                        help=''.join(['Runtime keys - each key links to a pre-determined list of variables to analyse. ',
                                      'Keys are: ', ', '.join( accepted_keys)]),
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


def run():
    from ._version import __version__
    print(f'BGCVal2: {__version__}')

    args = get_args()
    jobIDs = args.jobID
    keys = args.keys
    years = args.years

    # Note that these may not all work, as p2p will be different
    accepted_keys = [os.path.splitext(os.path.basename(fn))[0] for fn in glob('key_lists/*.yml')]

    good_keys = True
    for key in keys:
        if key not in accepted_keys:
            print('Key Argument [',key,'] nor recognised. Accepted keys are:', accepted_keys)
            good_keys= False
    if not good_keys:
        sys.exit(1)

    # user confiug file:
    if os.path.isfile(args.config_file):
        config_user = args.config_file
        print(f"analysis_timeseries: Using user config file {config_user}")
    else:
        print(f"analysis_timeseries: Could not find configuration file {config_user}."
              "Will proceed with defaults.")
        config_user = None

    print('Running analysis_p2p.\tjobID:', jobIDs, '\tkeys:', keys)


    for jobID in jobIDs:
        analysis_p2p(
            jobID=jobID,
            years=years,
            modelGrid='eORCA1',
            annual=True,
            noPlots=False,
            suites=keys,
            config_user=config_user,
            )
