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
.. module:: analysis_level3_sameYear
   :platform: Unix
   :synopsis: A script to produce a level 3 analysis comparing the same year in two jobs.

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
import re

#####
# Load specific local code:
from bgcval2.bgcvaltools import UKESMpython as ukp
from .timeseries import timeseriesAnalysis
from .timeseries import profileAnalysis
from .timeseries import timeseriesPlots as tsp

#####
# User defined set of paths pointing towards the datasets.
from .Paths import paths as paths


def listModelDataFiles(jobID, filekey, datafolder, annual):
    if annual:
        return sorted(
            glob(datafolder + jobID + "/" + jobID + "o_1y_*_" + filekey +
                 ".nc"))
    else:
        return sorted(
            glob(datafolder + jobID + "/" + jobID + "o_1m_*_" + filekey +
                 ".nc"))


def overlapyears(files1, files2):
    """ 
		Finds out which years overlap.
	"""
    filepairs = {}
    for f1 in sorted(files1):
        f1b = os.path.basename(f1)
        f1y = re.findall(r'\d+', f1b)[2]
        for f2 in files2:
            f2b = os.path.basename(f2)
            f2y = re.findall(r'\d+', f2b)[2]

            if f2y == f1y: filepairs[f2y] = [f1, f2]
    return filepairs


def maskAndCompress(la, lo, d1, d2):
    m = np.ma.array(d1).mask
    m += np.ma.array(d2).mask

    return np.ma.masked_where(m, la).compressed(), np.ma.masked_where(
        m, lo).compressed(), np.ma.masked_where(
            m, d1).compressed(), np.ma.masked_where(m, d2).compressed()


def analysis_sy(
    jobID1='u-af983',
    jobID2='u-ah531',
):
    annual = True
    """
	# Build this list with:	
	for k in nc.variables.keys():
		ndim = nc.variables[k].ndim
		if ndim<=2:continue	
		print "analysisKeys.append('"+k+"')"
	"""
    analysisKeys = []
    #	analysisKeys.append('SST')
    #	analysisKeys.append('SSS')
    #        analysisKeys.append('Ice')
    analysisKeys.append('votemper')
    analysisKeys.append('votemper2')
    analysisKeys.append('vosaline')
    analysisKeys.append('vosaline2')
    analysisKeys.append('ttrd_ldf')
    analysisKeys.append('ttrd_zdf')
    analysisKeys.append('strd_ldf')
    analysisKeys.append('strd_zdf')
    analysisKeys.append('sossheig')
    analysisKeys.append('zossq')
    analysisKeys.append('sowaflup')
    analysisKeys.append('sosafldo')
    analysisKeys.append('sohefldo')
    analysisKeys.append('soshfldo')
    analysisKeys.append('somixhgt')
    analysisKeys.append('somxl010')
    analysisKeys.append('sokaraml')
    analysisKeys.append('somlddbm')
    analysisKeys.append('soicecov')
    analysisKeys.append('sowindsp')
    analysisKeys.append('sohflisf')
    analysisKeys.append('sowflisf')
    analysisKeys.append('berg_total_melt')
    analysisKeys.append('berg_total_heat_flux')
    analysisKeys.append('sorunoff')

    analysisDict = {}
    imagedir = ukp.folder(paths.imagedir + '/' + jobID1 + '-' + jobID2 +
                          '/Level3/')
    #shelvedir 	= ukp.folder(paths.shelvedir+'/'+jobID+'/Level3/'+jobID1+'-'+jobID2)

    dataD = {}
    modeldataD = {}

    files1 = listModelDataFiles(jobID1, 'grid_T', paths.ModelFolder_pref,
                                annual)
    files2 = listModelDataFiles(jobID2, 'grid_T', paths.ModelFolder_pref,
                                annual)

    filepairs = overlapyears(files1, files2)
    """
	# Build this list with:
	for k in nc.variables.keys():
		ndim = nc.variables[k].ndim
		if ndim<=2:continue
		print "plotDetails['"+k+"'] = {'name':'"+k+"', 'key':'"+k+"', 'ndim':",
		try: print nc.variables[k].ndim,", 'longname': '",
		except: print "'longname': '",
		try: print nc.variables[k].long_name,"'}"
		except: print "'}"

	"""
    plotDetails = {}
    #	plotDetails['SST'] = {'name'='SST', 'key' = 'votemper', 'longname' = 'Sea Surface Temperature'}
    #	plotDetails['SSS'] = {'name'='SSS', 'key' = 'vosaline', 'longname' = 'Sea Surface Salinity'}
    #        plotDetails['Ice'] = {'name':'Ice', 'key':'soicecov', 'ndim':4, 'longname':'Ice fraction'}

    plotDetails['votemper'] = {
        'name': 'votemper',
        'key': 'votemper',
        'ndim': 4,
        'longname': ' temperature '
    }
    plotDetails['votemper2'] = {
        'name': 'votemper2',
        'key': 'votemper2',
        'ndim': 4,
        'longname': ' temperature '
    }
    plotDetails['vosaline'] = {
        'name': 'vosaline',
        'key': 'vosaline',
        'ndim': 4,
        'longname': ' salinity '
    }
    plotDetails['vosaline2'] = {
        'name': 'vosaline2',
        'key': 'vosaline2',
        'ndim': 4,
        'longname': ' salinity '
    }
    plotDetails['ttrd_ldf'] = {
        'name': 'ttrd_ldf',
        'key': 'ttrd_ldf',
        'ndim': 4,
        'longname': ' temperature-trend: lateral  diffusion '
    }
    plotDetails['ttrd_zdf'] = {
        'name': 'ttrd_zdf',
        'key': 'ttrd_zdf',
        'ndim': 4,
        'longname': ' temperature-trend: vertical diffusion '
    }
    plotDetails['strd_ldf'] = {
        'name': 'strd_ldf',
        'key': 'strd_ldf',
        'ndim': 4,
        'longname': ' salinity   -trend: lateral  diffusion '
    }
    plotDetails['strd_zdf'] = {
        'name': 'strd_zdf',
        'key': 'strd_zdf',
        'ndim': 4,
        'longname': ' salinity   -trend: vertical diffusion '
    }
    plotDetails['sossheig'] = {
        'name': 'sossheig',
        'key': 'sossheig',
        'ndim': 3,
        'longname': ' sea surface height '
    }
    plotDetails['zossq'] = {
        'name': 'zossq',
        'key': 'zossq',
        'ndim': 3,
        'longname': ' square of sea surface height '
    }
    plotDetails['sowaflup'] = {
        'name': 'sowaflup',
        'key': 'sowaflup',
        'ndim': 3,
        'longname': ' Net Upward Water Flux '
    }
    plotDetails['sosafldo'] = {
        'name': 'sosafldo',
        'key': 'sosafldo',
        'ndim': 3,
        'longname': ' Downward salt flux '
    }
    plotDetails['sohefldo'] = {
        'name': 'sohefldo',
        'key': 'sohefldo',
        'ndim': 3,
        'longname': ' Net Downward Heat Flux '
    }
    plotDetails['soshfldo'] = {
        'name': 'soshfldo',
        'key': 'soshfldo',
        'ndim': 3,
        'longname': ' Shortwave Radiation '
    }
    plotDetails['somixhgt'] = {
        'name': 'somixhgt',
        'key': 'somixhgt',
        'ndim': 3,
        'longname': ' Turbocline depth (Kz = 5e-4) '
    }
    plotDetails['somxl010'] = {
        'name': 'somxl010',
        'key': 'somxl010',
        'ndim': 3,
        'longname': ' Mixed Layer Depth (dsigma = 0.01 wrt 10m) '
    }
    plotDetails['sokaraml'] = {
        'name': 'sokaraml',
        'key': 'sokaraml',
        'ndim': 3,
        'longname': ' Mixed Layer Depth interpolated '
    }
    plotDetails['somlddbm'] = {
        'name': 'somlddbm',
        'key': 'somlddbm',
        'ndim': 3,
        'longname': ' Mixed Layer Depth interpolated '
    }
    plotDetails['soicecov'] = {
        'name': 'soicecov',
        'key': 'soicecov',
        'ndim': 3,
        'longname': ' Ice fraction '
    }
    plotDetails['sowindsp'] = {
        'name': 'sowindsp',
        'key': 'sowindsp',
        'ndim': 3,
        'longname': ' wind speed module '
    }
    plotDetails['sohflisf'] = {
        'name': 'sohflisf',
        'key': 'sohflisf',
        'ndim': 3,
        'longname': ' Ice Shelf Heat Flux '
    }
    plotDetails['sowflisf'] = {
        'name': 'sowflisf',
        'key': 'sowflisf',
        'ndim': 3,
        'longname': ' Ice shelf melting '
    }
    plotDetails['berg_total_melt'] = {
        'name': 'berg_total_melt',
        'key': 'berg_total_melt',
        'ndim': 3,
        'longname': ' icb melt rate 2 of icebergs '
    }
    plotDetails['berg_total_heat_flux'] = {
        'name': 'berg_total_heat_flux',
        'key': 'berg_total_heat_flux',
        'ndim': 3,
        'longname': ' icb latent heat of melting of icebergs '
    }
    plotDetails['sorunoff'] = {
        'name': 'sorunoff',
        'key': 'sorunoff',
        'ndim': 3,
        'longname': ' River Runoffs '
    }

    for ystr, [fp1, fp2] in list(filepairs.items()):
        print(ystr, [fp1, fp2])
        nc1 = Dataset(fp1, 'r')
        nc2 = Dataset(fp2, 'r')

        lons_cc = nc1.variables['nav_lon'][:]
        lats_cc = nc1.variables['nav_lat'][:]

        for n in analysisKeys:

            filename = ukp.folder(
                imagedir + ystr) + plotDetails[n]['name'] + '_' + ystr + '.png'
            if plotDetails[n]['ndim'] == 4:
                data1 = nc1.variables[plotDetails[n]['key']][0, 0]
                data2 = nc2.variables[plotDetails[n]['key']][0, 0]
            if plotDetails[n]['ndim'] == 3:
                data1 = nc1.variables[plotDetails[n]['key']][0]
                data2 = nc2.variables[plotDetails[n]['key']][0]

            lons, lats, data1, data2 = maskAndCompress(lons_cc, lats_cc, data1,
                                                       data2)

            try:
                ukp.robinPlotQuad(
                    lons,
                    lats,
                    data1,
                    data2,
                    filename,
                    titles=[jobID1, jobID2],
                    title=plotDetails[n]['longname'] + ' ' + ystr[:4] + '-' +
                    ystr[4:6] + '-' + ystr[6:],
                    vmin='',
                    vmax='',
                )  #maptype='Basemap')
            except:
                print("didn't work:", filename)


def main():
    if "--help" in argv or len(argv) == 1:
        print("Running with no arguments. Exiting.")
        if "--help" in argv:
            print("Read the documentation.")
        sys.exit(0)
    try:
        jobID1 = argv[1]
        jobID2 = argv[2]
    except:
        jobID = "u-ad371"
        jobID = "u-ad371"

    if 'debug' in argv[1:]: suite = 'debug'
    else: suite = 'normal'

    analysis_sy(
        jobID1=jobID1,
        jobID2=jobID2,
    )


if __name__ == "__main__":
    main()
