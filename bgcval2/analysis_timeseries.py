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
import getpass
import itertools
import yaml
import importlib
import pathlib

#####
# Load specific local code:
from bgcval2 import UKESMpython as ukp
from bgcval2.timeseries import timeseriesAnalysis
from bgcval2.timeseries import profileAnalysis
from bgcval2.timeseries import timeseriesTools as tst

from bgcval2.bgcvaltools.mergeMonthlyFiles import mergeMonthlyFiles, meanDJF
from bgcval2.bgcvaltools.AOU import AOU
from bgcval2.bgcvaltools.dataset import dataset
from bgcval2._runtime_config import get_run_configuration
from bgcval2.functions.standard_functions import std_functions

#####
# User defined set of paths pointing towards the datasets.
from bgcval2.Paths.paths import paths_setter


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
        model_files = sorted(
            glob(datafolder + "/" + jobID + "o_1m_*_" + filekey +
                 ".nc"))

    return model_files


def list_input_files(files_path, key_dict, paths):
    """
    Generate a list of files from a path, which may include
    several $PATH values.
    """
    #####
    # Replace some values for $FLAGS in the path:
    flags = ['USERNAME','basedir_model', 'basedir_obs','PATHS_GRIDFILE', 'PATHS_BGCVAL2']
    flag_values = [getpass.getuser(), paths.ModelFolder_pref, paths.ObsFolder, paths.orcaGridfn, paths.bgcval2_repo]

    for flag in ['jobID', 'model', 'years','year', 'scenario', 'name']:
        if key_dict.get(flag, False):
            flags.append(flag.upper())
            flag_values.append(key_dict[flag])

    for flag, flag_value in zip(flags, flag_values):
        files_path = findReplaceFlag(files_path, flag, flag_value)

    basedir = os.path.dirname(files_path)
    if not glob(basedir):
        raise OSError(f"Base {basedir} is not a valid directory.")
    input_files = sorted(glob(files_path))
    if not input_files:
        raise FileNotFoundError(f"Data dir {basedir} does not contain and file "
                                f"that matches pattern {files_path}, {flag}")
    return input_files


def build_list_of_suite_keys(suites, debug=True):
    """
    Generate a list of keys from a list of suites.

    """
    print('analysis_timeseries: Calling build_list_of_suite_keys to build list of keys')
    paths_dir = os.path.dirname(os.path.realpath(__file__))
    key_lists_dir = os.path.join(os.path.dirname(paths_dir), 'key_lists')

    print(f'analysis_timeseries: Directory where keys are stored: {key_lists_dir}')
    analysis_keys = {}
    for suite in suites:
        if not suite or suite in [' ', ',']:
            continue

        # look for a list in keys_list directory:
        suite_yml = os.path.join(key_lists_dir, ''.join([suite.lower(),'.yml']))
        if debug:
            print('build_list_of_suite_keys:\tlooking for suite yaml:', suite_yml)

        if not os.path.exists(suite_yml):
            print(f'analysis_timeseries: build_list_of_suite_keys:\tERROR: suite yaml: {suite_yml} file does not exist')
            sys.exit(1)

        # Open yml file:
        with open(suite_yml, 'r') as openfile:
            suite_dict = yaml.safe_load(openfile)

        if not suite_dict or not isinstance(suite_dict, dict):
            print(f"Configuration file {suite_yml} "
                  "is either empty or corrupt, please check its contents")
            sys.exit(1)

        keys_dict = suite_dict.get('keys', {})

        for key, keybool in keys_dict.items():
            if debug and key in analysis_keys:
                print(f'build_list_of_suite_keys:\tKey {key} exists in multiple suites:')

            if key in analysis_keys and keybool != analysis_keys[key]:
                print(f'build_list_of_suite_keys:\tERROR: conflict in input yamls: {key}, {keybool} != {analysis_keys[key]}')
                sys.exit(1)

            if keybool:
                if debug:
                    print('build_list_of_suite_keys:\tAdding key:', key)

                analysis_keys[key] = keybool
    analysis_keys = [key for key in analysis_keys.keys()]
    return analysis_keys


def get_kwargs_from_dict(convert_dict, avoids = ['path', 'function']):
    """
    Get the KWARGS fromn a dict
    """
    kwargs = {}
    for key, value in convert_dict.items():
        if key in avoids:
            continue
        kwargs[key] = value
    return kwargs


def load_function(convert):
    """
    Using the named function in the key yaml, load that function and return it and it kwargs.
    """
#    if convert.find(' ')>-1: 
#        converts = convert.split(' ')
#        return load_function[converts[0]](load_function[converts[1]])
    # function is listed as a standard function string
    if isinstance(convert, str) and convert in std_functions:
        print( "Standard Function Found:", convert)
        return std_functions[convert], {}

    functionname = convert.get('function', False)
    if not functionname:
        raise KeyError(f'Function name not provided in convert dict: {convert}')

    # function is listed as a standard function string with kwargs:
    if functionname in std_functions:
        print( "Standard Function Found:", functionname)
        kwargs = get_kwargs_from_dict(convert)
        return std_functions[functionname], kwargs

    func_path = convert.get('path', False)
    if not func_path:
        raise KeyError(f'No path provided in {str(convert)} dictionary')

    if func_path[:7] == 'bgcval2':
        # path is relative to bgcval2 repo.
        repo_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
        func_path = os.path.join(repo_dir, func_path)

    if not os.path.exists(func_path):
        raise FileNotFoundError(f'Path {func_path} to custom function file not found.')

    # load function from Python file in path
    modulename = os.path.splitext(os.path.basename(func_path))[0]
    loader = importlib.machinery.SourceFileLoader(modulename, func_path)
    module = loader.load_module()
    func = getattr(module, functionname)

    kwargs = get_kwargs_from_dict(convert)

    print('Successfully loaded function:', functionname, 'from', func_path, 'with kwargs:', kwargs)
    return func, kwargs


def findReplaceFlag(filepath, flag, new_value):
    """
    Looking for $FLAG in the filepath.
    If found, we replace $FLAG with new_value
    """
    lookingFor = ''.join(['$', flag.upper()])
    if filepath.find(lookingFor) == -1:
        print('NOT Changing FLAG:', lookingFor, 'to', new_value, 'in', filepath)

        return filepath
    print('Changing FLAG:', lookingFor, 'to', new_value, 'in', filepath)
    filepath = filepath.replace(lookingFor, new_value)
    return filepath


def parse_list_from_string(list1):
    """
    This tool takes a string and returns it as a list.
    """
    if isinstance(list1, list):
       return list1
    list1 = list1.replace('\t', ' ').replace(',', ' ').replace('  ', ' ')
    list1 = list1.replace('\'', '').replace('\"', '').replace(';', '')
    while list1.count('  ')>0:
        list1 = list1.replace('  ', ' ')
    #if len(list1)==0: return []
    return list1.split(' ')


def load_key_file(key, paths, jobID):
    """
    load_key_file:
        takes an input file from key_files directory
        and outputs a dictionary.
    """
    output_dict = {}
    print('load_key_file: checking {key}', key)
    key_yml_path = os.path.join(os.path.join(paths.bgcval2_repo, 'key_files', ''.join([key.lower(),'.yml'])))
    if os.path.exists(key_yml_path):
        print('key_yml_path exists:', key_yml_path)
    else:
        raise FileNotFoundError(f'{key_yml_path} does not exist')

    # Open yml file:
    with open(key_yml_path, 'r') as openfile:
        key_dict = yaml.safe_load(openfile)

    if not key_dict or not isinstance(key_dict, dict):
        print(f"Key Yaml file {key_yml_path} "
              "is either empty or corrupt, please check its contents")
        sys.exit(1)

    # add jobID to dict:
    key_dict['jobID'] = jobID

    # Load basic fields:
    for field in ['name', 'units', 'dimensions', 'modelgrid']:
        output_dict[field] = key_dict.get(field, None)
        print('Adding', field,':', output_dict[field])

    # Lists:
    for field in ['layers', 'regions', ]:
        output_dict[field] = key_dict.get(field, '')
        output_dict[field] = parse_list_from_string(output_dict[field])
        print('Adding', field,':', output_dict[field])

    # Metrics
    if key_dict['dimensions'] in  [1,]:
        metricList = ['metricless',]
    else:
        metricList = ['mean', ] #'median', '10pc','20pc','30pc','40pc','50pc','60pc','70pc','80pc','90pc','min','max']
    output_dict['metrics'] = key_dict.get('metrics', metricList)
    output_dict['metrics'] = parse_list_from_string(output_dict['metrics'])

    # Load Grid:
    gridFile = key_dict.get('gridFile', paths.orcaGridfn)
    output_dict['gridFile'] = list_input_files(gridFile, key_dict, paths)[0]

    # load model or data specific parts:
    for model_or_data in ['model', 'data']:
        md_vars = key_dict.get(''.join([model_or_data, '_vars']), False)
        if model_or_data == 'data' and md_vars is False:
            # Some analyses don't have observational data.
            continue

        if model_or_data == 'model' and md_vars is False:
            raise KeyError(f'What are you trying to analyse: model_vars {md_vars} is empty in yml_file: {key_yml_path}')

        md_vars = parse_list_from_string(md_vars)
        
        md_convert = key_dict.get(''.join([model_or_data, '_convert']), 'NoChange')
        convert_func, kwargs = load_function(md_convert)
        # source:
        model_or_data_name = key_dict.get(model_or_data, model_or_data)
        source = key_dict.get(''.join([model_or_data, '_source']), model_or_data)
        output_dict[model_or_data] = model_or_data_name

        output_dict[''.join([model_or_data, '_source'])] = source
        output_dict[''.join([model_or_data, 'details'])] = {
            'name': key_dict['name'],
            'vars': md_vars ,
            'convert': convert_func,
            'units': key_dict['units'],
            }
        for kwarg, kwarg_value in kwargs.items():
            if isinstance(kwarg_value, str) and kwarg.lower().find('file')>-1:
                output_dict[''.join([model_or_data,'details'])][kwarg] = list_input_files(kwarg_value, key_dict, paths)
            else:
                output_dict[''.join([model_or_data,'details'])][kwarg] = kwarg_value

        if model_or_data == 'model':
            file_path = key_dict[''.join([model_or_data, 'Files'])]
            mdfile = list_input_files(file_path, key_dict, paths)
            output_dict[''.join([model_or_data, 'Files'])] = mdfile
        else:
            # get data file path
            file_path = key_dict[''.join([model_or_data, 'File'])]
            mdfile = list_input_files(file_path, key_dict, paths)
            if isinstance(mdfile, list) and len(mdfile) == 1:
                mdfile = mdfile[0]

        output_dict[''.join([model_or_data, 'Files'])] = mdfile

        coords = ukp.load_coords_from_netcdf(mdfile)
        output_dict[''.join([model_or_data, 'coords'])] = {
            'tdict': ukp.tdicts[key_dict.get('tdict', 'ZeroToZero')],
            }
        for coord, value in coords.items():
            # Coordinate names are guessed, but can be over-written in the yaml.
            coord_in_yml = ''.join([model_or_data, '_', coord])
            if  coord_in_yml in key_dict:
                value = key_dict[coord_in_yml]
                print('loading coord coord_in_yml:',coord_in_yml, value)
            output_dict[''.join([model_or_data, 'coords'])][coord] = value
    return output_dict


def analysis_timeseries(
    jobID="u-ab671",
    suites=['all', ],
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

	suites chooses a set of fields to look at.

	regions selects a list of regions, default is 'all', which is the list supplied by Andy Yool.

	:param jobID: the jobID
	:param clean: deletes old images if true
	:param annual: Flag for monthly or annual model data.
	:param strictFileCheck: CStrickt check for model and data files. Asserts if no files are found.
	:param suites: Which data to analyse, ie level1, physics only, debug, etc
	:param regions:

	"""

    print('-----------------------')
    print('Starting analysis_timeseries')
    print('jobID:', jobID)
    print('suites:', suites)
    print('regions:', regions)
    print(f'clean: {clean},  annual: {annual}, strictFileCheck: {strictFileCheck}')
    print('config_user:', config_user)

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

    if regions in ['debug', 'Global', 'spinup']:
        regionList = ['Global', ]

        metricList = [
            'mean',
        ]

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

    shelvedir = ukp.folder([paths.shelvedir, "timeseries", jobID])
    imagedir = ukp.folder([paths.imagedir, jobID, 'timeseries'])

    if annual: WOAFolder = paths.WOAFolder_annual
    else: WOAFolder = paths.WOAFolder

    #####
    # PML
    hostname = gethostname()

    if hostname.find('pmpc') > -1:
        print("analysis_timeseries:\tBeing run at PML on ", gethostname())

    #####
    # JASMIN
    if hostname.find('ceda.ac.uk') > -1 or hostname.find(
            'jasmin') > -1 or hostname.find('jc.rl.ac.uk') > -1:
        print("analysis_timeseries:\tBeing run at CEDA on ", hostname)

    #####
    # Unable to find location of files/data.
    if not paths.machinelocation:
        print(
            "analysis_timeseries:\tFATAL:\tWas unable to determine location of host: ",
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


#####
# Coordinate dictionary
# These are python dictionairies, one for each data source and model.
# This is because each data provider seems to use a different set of standard names for dimensions and time.
# The 'tdict' field is short for "time-dictionary".
#	This is a dictionary who's indices are the values on the netcdf time dimension.
#	The tdict indices point to a month number in python numbering (ie January = 0)
# 	An example would be, if a netcdf uses the middle day of the month as it's time value:
#		tdict = {15:0, 45:1 ...}

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

    # NEW STYLE keys from file:
    av = ukp.AutoVivification()
    for key in analysisKeys:
        av[key] = load_key_file(key, paths, jobID)

    #####
    # Calling timeseriesAnalysis
    # This is where the above settings is passed to timeseriesAnalysis, for the actual work to begin.
    # We loop over all fiels in the first layer dictionary in the autovificiation, av.
    #
    # Once the timeseriesAnalysis has completed, we save all the output shelves in a dictionary.
    # At the moment, this dictioary is not used, but we could for instance open the shelve to highlight specific data,
    #	(ie, andy asked to produce a table showing the final year of data.

    shelves = {}
    shelves_insitu = {}
    for name in list(av.keys()):
        print("---------------------------------------------------------------")
        print("analysis_timeseries:\tBeginning timeseriesAnalysis:", name)

        if 'modelFiles' not in av[name]:
            print(
                "analysis_timeseries:\tWARNING:\tmodel files are not found:",
                name, av[name]['modelFiles'])
            if strictFileCheck:
                raise FileNotFoundError(f'Model files are not provided for {name}')

        modelfilesexists = [os.path.exists(f) for f in av[name]['modelFiles']]
        if False in modelfilesexists:
            print(
                "analysis_timeseries:\tWARNING:\tModel files do not all exist:",
                av[name]['modelFiles'])
            missing_files = []
            for f in av[name]['modelFiles']:
                if os.path.exists(f):
                    continue
                print(f, 'does not exist')
                missing_files.append(f)
            if strictFileCheck:
                raise FileNotFoundError(f'Model files are provided but not found for {name} : {missing_files}')

        if 'dataFile' in av[name]:
            print(name, 'dataFile', av[name]['dataFile'])
            if not os.path.exists(av[name]['dataFile']):
                print(
                    "analysis_timeseries:\tWARNING:\tdata file is not found:",
                    av[name]['dataFile'])
                if strictFileCheck:
                    raise FileNotFoundError(f'Obs data files are provided but not found for {name} : {av[name]["dataFile"]}')

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
            clean=clean,
        )

        #####
        # Profile plots
        if av[name]['dimensions'] == 3:
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
            suites=[
                key,
            ],
        )


def singleTimeSeries(
    jobID,
    key,
):
    #	try:
    analysis_timeseries(jobID=jobID,
                        suites=[
                            key,
                        ],
                        strictFileCheck=False)  #clean=1)


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


def main():
    """
    Main function that does all the heavy lifting.
    """
    args = get_args()
    jobIDs = args.jobID
    keys = args.keys
    print('Running analysis_imerseries.\tjobID:', jobIDs, '\tkeys:', keys)
    accepted_keys = [os.path.splitext(os.path.basename(fn))[0] for fn in glob('key_lists/*.yml')]

    good_keys = True
    for key in keys:
        if key not in accepted_keys:
            print('Key Argument [',key,'] nor recognised. Accepted keys are:', accepted_keys)
            good_keys= False
    if not good_keys:
        sys.exit(1)

    if os.path.isfile(args.config_file):
        config_user = args.config_file
        print(f"analysis_timeseries: Using user config file {config_user}")
    else:
        print(f"analysis_timeseries: Could not find configuration file {config_user}."
              "Will proceed with defaults.")
        config_user = None

    for jobID in jobIDs:
        analysis_timeseries(
            jobID=jobID,
            suites=keys,
            config_user=config_user
        )


if __name__ == "__main__":
    from ._version import __version__
    print(f'BGCVal2: {__version__}')
    main()
