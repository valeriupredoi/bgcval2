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
from bgcval2.p2p.testsuite_p2p import testsuite_p2p
from bgcval2.p2p.summaryTargets import summaryTargets
from bgcval2.bgcvaltools.pftnames import months
from bgcval2.analysis_timeseries import build_list_of_suite_keys, load_key_file

#####
# User defined set of paths pointing towards the datasets.
from bgcval2.Paths.paths import paths_setter


def analysis_p2p(
    jobID='u-ad980',
    years=['1997'],  
    modelGrid='eORCA1',
    annual=True,
    noPlots=False,
    suites='default',
    noTargets=False, 
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

    shelvesAV = []
    #####
    # Location of model files.
    ModelFolder = ModelFolder_pref + jobID + "/"
    workingDir = ''.join([paths.p2p_ppDir, '/', jobID, ])
    imageFolder = ''.join([paths.imagedir, '/', jobID])

    print('P2P working dir:', workingDir)
    print('P2P image dir:', imageFolder)

    # Attempt to make the directories:
    workingDir = ukp.folder(workingDir)
    imageFolder = ukp.folder(imageFolder)
    
    print('p2p: av keys:', av.keys())

    shelvesAV.extend(
        testsuite_p2p(
            model='NEMO',
            jobID=jobID,
            years=years,
            av=av,
            workingDir=workingDir,
            imageFolder=imageFolder,
            noPlots=
            noPlots,  # turns off plot making to save space and compute time.
            gridFile=paths.orcaGridfn,  # enforces custom gridfile.
            annual=annual,
            noTargets=noTargets,
        ))
    if noTargets: return

    ######
    # Summary Target diagrams:
    imageFold = ukp.folder(imageFolder + '/Targets/Summary')
    summaryTargets(shelvesAV, imageFold, years)


def get_args():
    """Parse command line arguments. """
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
