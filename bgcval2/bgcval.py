#!/usr/bin/env python

#
# Copyright 2014, Plymouth Marine Laboratory
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
#####
#
"""
.. module:: theWholePackage
   :platform: Unix
   :synopsis: A script the run the entire suite of analyses and produce an html report.
.. moduleauthor:: Lee de Mora <ledm@pml.ac.uk>
"""

import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import sys
from sys import argv, exit
from multiprocessing import Pool

from .bgcvaltools.downloadFromMass import downloadMass, findLastFinishedYear
from .analysis_timeseries import analysis_timeseries, singleTimeSeries, singleTimeSeriesProfile
from .analysis_timeseries import level1KeysDict, physKeysDict, bgcKeysDict, keymetricsfirstDict
from .analysis_p2p import analysis_p2p, p2pDict_level2, p2pDict_physics, single_p2p
from .makeReport import html5Maker
from .UKESMpython import folder


def timeseriesParrallelL1(index):
    print("timeseriesParrallel", index, jobID, 'START')
    key = level1KeysDict[index]
    singleTimeSeries(
        jobID,
        key,
    )
    print("timeseriesParrallel", index, jobID, 'SUCESS', key)


def timeseriesParrallelKMF(index):
    print("timeseriesParrallel", index, jobID, 'START')
    key = keymetricsfirstDict[index]
    singleTimeSeries(
        jobID,
        key,
    )
    print("timeseriesParrallel", index, jobID, 'SUCESS', key)


def timeseriesParrallelPhys(index):
    key = physKeysDict[index]
    print("timeseriesParrallelPhys", index, jobID, 'START', key, index)
    try:
        singleTimeSeries(
            jobID,
            key,
        )
    except:
        print("timeseriesParrallelPhys failed for", index, jobID, key)
        assert 0
    print("timeseriesParrallelPhys", index, jobID, 'SUCESS', key)


def timeseriesParrallelBGC(index):
    key = bgcKeysDict[index]
    print("timeseriesParrallelBGC", index, jobID, 'START', key, index)
    try:
        singleTimeSeries(
            jobID,
            key,
        )
    except:
        print("timeseriesParrallelBGC failed for", index, jobID, key)
        assert 0
    print("timeseriesParrallelBGC", index, jobID, 'SUCESS', key)


def p2pParrallel(index):
    print("p2pParrallel", index, jobID, 'START')
    key = p2pDict_level2[index]
    single_p2p(jobID, key, year)
    print("p2pParrallel", index, jobID, 'SUCESS', key)


def p2pParrallel_phys(index):
    print("p2pParrallel_phys", index, jobID, 'START')
    key = p2pDict_physics[index]
    single_p2p(jobID, key, year)
    print("p2pParrallel_phys", index, jobID, 'SUCESS', key)


def theWholePackage(jobID, year=False, suite='level1'):
    """
	theWholePackage function. This function runs the whole default package analysis suite
	and outputs it to an html report.
	
	:param jobID: The jobID of the run.
	:param year: The year of the run, False will search for a year, '*' will not run point to point.	    
	:param suite: The type of analysis. Options are 'Physics',
	
	This function calls:
		the html5Maker
		analysis_timeseries
		analysis_p2p
		
	"""
    #if year in [False,  '*']:
    #        year = findLastFinishedYear(jobID,dividby=25)
    #elif type(year) in [type(1000), type(1000.)]:
    #	year = str(year)

    print("########\nThe Whole Package:\tStarting job", jobID, year)
    parrallel = False
    cores = 1

    if suite == 'physics':
        physicsOnly = True
        parrallel = False
    else:
        physicsOnly = False
#        print "########\nThe Whole Package:\tmaking Summary report"
#       html5Maker(jobID =jobID,
#                 reportdir=folder('reports/'+jobID),
#                year = year,
#               clean=True,
#              physicsOnly=physicsOnly
#             )
#	return

    print("########\nThe Whole Package:\tStarting Time series (surface only)",
          jobID)
    if parrallel:

        #####
        # Running key metrics first.
        if suite == 'level1':
            remaining = sorted(keymetricsfirstDict.keys())[:]
            print('Running key metrics first:', jobID, suite, remaining)
            #assert 0
            p = Pool(cores)
            p.map(timeseriesParrallelKMF, remaining)
            p.close()

        #if suite =='all':	remaining = sorted(timeseriesDict.keys())[:]
        if suite == 'level1': remaining = sorted(level1KeysDict.keys())[:]
        if suite == 'physics': remaining = sorted(physKeysDict.keys())[:]
        if suite == 'bgc': remaining = sorted(bgcKeysDict.keys())[:]

        p = Pool(cores)
        if suite == 'level1': p.map(timeseriesParrallelL1, remaining)
        if suite == 'physics': p.map(timeseriesParrallelPhys, remaining)
        if suite == 'bgc': p.map(timeseriesParrallelBGC, remaining)
        p.close()
    else:
        analysis_timeseries(
            jobID=jobID,
            analysisSuite=suite,
        )  #z_component = 'SurfaceOnly',)


# Turning off P2P analysis:
    year = False
    if year not in ['*', False]:
        if suite == 'physics': pass
        else: suite = 'level2'

        print(
            "########\nThe Whole Package:\tRunning point to point analysis of",
            jobID, "on", year, 'in suite:', suite)

        if parrallel:
            p1 = Pool(cores)
            if suite == 'physics':
                remaining = sorted(p2pDict_physics.keys())[:]
                p1.map(p2pParrallel_phys, remaining)
            else:
                remaining = sorted(p2pDict_level2.keys())[:]
                p1.map(p2pParrallel, remaining)
            p1.close()

            #####
            # And once over to make the summary target diagrams.
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
                noPlots=True,
                analysisSuite=suite,
            )

        else:
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
                analysisSuite=suite,
            )

    else:
        print(
            "########\nThe Whole Package:\tNot Running point to point analysis of",
            jobID, " because year is:", year)

    print("########\nThe Whole Package:\tmaking Final Summary report")

    if suite == 'physics': physicsOnly = True
    else: physicsOnly = False
    html5Maker(jobID=jobID,
               reportdir=folder('reports/' + jobID),
               year=year,
               clean=True,
               physicsOnly=physicsOnly)

def run():
    from ._version import __version__
    print(f'BGCVal2: {__version__}')
    if "--help" in argv or len(argv) == 1:
        print("Running with no arguments. Exiting.")
        if "--help" in argv:
            print("Read the documentation.")
        sys.exit(0)
    try:
        jobID = argv[1]
    except:
        print("Please provide a job ID")
        exit()
    #if 'ReportOnly' in argv[:]:ReportOnly=True
    #else:	ReportOnly = False

    if 'physics' in argv[:]:
        physicsOnly = True
        numberfiles = 4
    else:
        physicsOnly = False
        numberfiles = 6

    year = False
    for ar in argv:
        try:
            ar = int(ar)
        except:
            continue
        year = str(ar)
    for divby in [100, 50, 25, 10, 5, 1]:
        print("main", divby, year)
        if year: continue
        year = findLastFinishedYear(jobID,
                                    dividby=divby,
                                    numberfiles=numberfiles)
    print("########\nThe Whole Package:\tmain:", jobID, year)

    #if not ReportOnly:
    #year = False
    print(
        "Warning:\tturning off p2p analysis to save time and speed up time series."
    )
    if physicsOnly: theWholePackage(jobID, year=year, suite='physics')
    else: theWholePackage(jobID, year=year)

#  if year == False: year = '*'
#	print "########\nThe Whole Package:\tmaking Summary report", jobID,year
#	html5Maker(jobID =jobID,
#		   reportdir=folder('reports/'+jobID),
#		   year = year,
#		   clean=True,
#		   physicsOnly=physicsOnly,
#		   )
