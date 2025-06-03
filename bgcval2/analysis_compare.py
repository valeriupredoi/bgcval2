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
from .bgcvaltools import bv2tools as bvt
from .timeseries import timeseriesAnalysis
from .timeseries import profileAnalysis
from .timeseries import timeseriesPlots as tsp
from .timeseries import timeseriesTools as tst

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
       print('apply_timerange: timerange is', None) 
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
        print('apply_timerange: WARNING: No times made the cut?', jobID, len(times),
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
                       labels = {},
                       ensembles={},
                       config_user=None,
                       dpi=None,
                       savepdf=False,
                       savejson=False,
    ):
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
        csvFolder = paths.imagedir + '/TimeseriesCompare_CSV/' + analysisname

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
        imagedir = bvt.folder(paths.imagedir + '/' + jobID + '/timeseries')
        shelvedir = bvt.folder(paths.shelvedir + "/timeseries/" + jobID)

        if jobID in list(ensembles.keys()): continue
        # ensembles names can not be the same as jobIDs

        av = bvt.AutoVivification()

        # NEW STYLE keys from file:
        for key in analysisKeys:
            av[key] = load_key_file(key, paths, jobID)

        for name in list(av.keys()):
            print(
                "------------------------------------------------------------------"
            )
            print(
                "analysis-Timeseries.py:\tBeginning to call timeseriesAnalysis for ",
                name)

            if not av[name]['modelFiles']:
                file_err = "analysis-Timeseries.py:\tWARNING:\tmodel files are not " \
                           f"found: {av[name]['modelFiles']} for {jobID}"
                print(file_err)
                if strictFileCheck:
                    raise FileNotFoundError(file_err)

            modelfilesexists = [
                os.path.exists(f) for f in av[name]['modelFiles']
            ]
            if False in modelfilesexists:
                print(
                    "analysis-Timeseries.py:\tWARNING:\tnot model files do not all exist:",
                    av[name]['modelFiles'])
                if strictFileCheck:
                    raise FileError('Model Files are not found jobID: %s, name: %s', jobID, name)

            if 'dataFile' in av[name] and not os.path.exists(av[name]['dataFile']):
                print(
                    "analysis-Timeseries.py:\tWARNING:\tdata file is not found:",
                    av[name]['dataFile'])
                if strictFileCheck:
                    raise FileError('Data Files are not found jobID: %s, name: %s', jobID, name)

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
        regions = av[name]['regions']
        layers = av[name]['layers']
        metrics = av[name]['metrics']
        smoothings = av[name]['smoothings']
        for region, layer, metric  in itertools.product(regions, layers, metrics):
            timesD = {}
            arrD = {}
            for jobID in jobs:
                try:
                    mdata = modeldataD[(jobID, name)][(region, layer, metric)]
                except:
                    continue
                title = titleify([region, layer, metric, name])

                times, datas = apply_shifttimes(mdata, jobID, shifttimes)
                print('post apply_shifttimes:', jobID, len(times), len(datas))
                times, datas = apply_timerange(times, datas, jobID, timeranges)
                timesD[jobID] = times 
                arrD[jobID] = datas  
                print(jobID, region, layer, metric, len(times), len(datas))
            timesD, arrD = build_ensemble(timesD, arrD, ensembles)

            if len(list(arrD.keys())) == 0: 
                continue
            units = av[name]['modeldetails']['units']

            ts = 'Together'

            if savejson:
                tst.save_json(
                        timesD,
                        arrD, 
                        analysisname,
                        colours=colours,
                        linestyles=linestyles,
                        thicknesses=lineThicknesses,
                        labels=labels,
                        name=name,
                        units=units,
                        region=region,
                        layer=layer,
                        metric=metric,
                        title=title,
                        ts=ts,
                        csvFolder=csvFolder)
                    
            for smoothing in smoothings:
            #or ls in ['DataOnly', ]:
                tsp.multitimeseries(
                    timesD,  # model times (in floats)
                    arrD,  # model time series
                    data=-999,  # in situ data distribution
                    title=title,
                    filename=bvt.folder(imageFolder) +
                        '_'.join([name, region, layer, metric, ts, smoothing + '.png']),
                    units=units,
                    plotStyle=ts,
                    smoothing=smoothing,
                    colours=colours,
                    thicknesses=lineThicknesses,
                    linestyles=linestyles,
                    labels=labels,
                    dpi=dpi,
                    savepdf=savepdf,
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

    # Send everything to the comparison maker:
    comparehtml5Maker(
        jobIDs=jobs,
        reportdir=bvt.folder('CompareReports2/' + analysisname),
        files=AllImages,
        clean=False,
        jobDescriptions=jobDescriptions,
        jobColours=colours,
        paths=paths,
        analysisKeys=analysisKeys,
    )
    print('End of timeseries_compare')


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
    details['strictFileCheck'] = input_yml_dict.get('strict_file_check', True)

    # Image output settings:
    # dpi: pixels per inch (image resolution)
    # savepdf: also save the image as a pdf. 
    # savejson: Save the data that appears in the image.
    details['dpi'] = input_yml_dict.get('dpi', None)
    details['savepdf'] = input_yml_dict.get('savepdf', False)
    details['savejson'] = input_yml_dict.get('savejson', False)



    if details['dpi']: # None is valid!
        try: 
            int(details['dpi']) 
        except: 
            raise ValueError(''.join(["Loading yml error: `dpi` needs to be an integer. Current value:",
                                      str(details['dpi'])]))

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
    labels = {}
   
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
        labels[jobID] = job_dict.get('label', jobID)
       
    details['colours'] = colours
    details['descriptions'] = descriptions
    details['thicknesses'] = thicknesses
    details['linestyles'] = linestyles
    details['shifttimes'] = shifttimes
    details['timeranges'] = timeranges
    details['labels'] = labels
    details['suites'] = suites
    details['auto_download'] = auto_download_dict
    return details


def load_yml_and_run(compare_yml, config_user, skip_timeseries):
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

    if skip_timeseries is None:
        pass
    else:
        do_analysis_timeseries = not skip_timeseries

    colours = details['colours']
    thicknesses = details['thicknesses']
    linestyles = details['linestyles']
    descriptions = details['descriptions']
    shifttimes = details['shifttimes']
    timeranges = details['timeranges']
    labels = details['labels']
    suites = details['suites']
    auto_download = details['auto_download']
    strictFileCheck = details.get('strictFileCheck', True)
    dpi = details.get('dpi', None)
    savepdf = details.get('savepdf', False)
    savejson = details.get('savejson', False)

    print('---------------------')
    print('timeseries_compare:',  analysis_name)
    print('job ids:', jobs.keys())
    for jobID in jobs:
        print(jobID, 'description:',descriptions[jobID])
        print(jobID, 'colour:',colours[jobID])
        print(jobID, 'line thickness & style:',thicknesses[jobID], linestyles[jobID])
        print(jobID, 'Shift time by', shifttimes[jobID])
        print(jobID, 'Label: ', labels[jobID])
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
                config_user=config_user,
                strictFileCheck=strictFileCheck,
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
        labels=labels,
        config_user=config_user,
        dpi=dpi,
        savepdf=savepdf,
        savejson=savejson,

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

    parser.add_argument('--skip_timeseries',
                        '-s',
                        default=None, 
                        help='When True: skip the new timeseries analyses and make the html report. Overwrites the do_analysis_timeseries flag in input_yml.',
                        action=argparse.BooleanOptionalAction,
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
        print(f"analysis_compare: Comparison config file {compare_yml}")

        if not os.path.isfile(compare_yml):
            print(f"analysis_compare: Could not find comparison config file {compare_yml}")
            sys.exit(1)
        skip_timeseries = args.skip_timeseries  
        load_yml_and_run(compare_yml, config_user, skip_timeseries)

    print("Finished... ")


if __name__ == "__main__":
    main()
