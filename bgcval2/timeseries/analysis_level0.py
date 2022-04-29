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
"""
.. module:: analysis_level0
   :platform: Unix
   :synopsis: A tool to produce a value for the html5 level 0 table.
.. moduleauthor:: Lee de Mora <ledm@pml.ac.uk>

"""

from shelve import open as shopen
from sys import argv, exit
import os
import numpy as np
from ..bgcvaltools.pftnames import getLongName
from ..Paths.paths import shelvedir

def printableName(field, region, layer, metric):
    #####
    #Produce a printable name
    name = field
    for v in [region, layer, metric]:
        if v in ['regionless', 'layerless', 'metricless']: continue
        name += ' ' + getLongName(v)
    return name


def analysis_level0(
    jobID='',
    field="AMOC_26N",
    region='regionless',
    layer='layerless',
    metric='metricless',
    debug=False,
):
    """
	Analysis_level0 loads the result of the shelves in  
	
	The path to the shelves is determined by the path.py symbolic link.
	
	:param jobID: A job run ID string
	:param field: A specific field to analyse.
	:param region:
	:param layer:  The depth layer, ie Surface, 1000m, etc
	:param metric:	Mean, max, median or min. The metric 
	"""
    #####
    #Produce a printable name
    name = printableName(field, region, layer, metric)

    #####
    # Start and load shelve
    if debug: print('analysis_level0:', jobID, field, region, layer, metric)
    shelvefn = shelvedir + "/timeseries/" + jobID + "/" + jobID + "_" + field + ".shelve"

    if debug: print('analysis_level0:', shelvefn, os.path.exists(shelvefn))
    if not os.path.exists(shelvefn):
        print("This shelve fn doesn't exists", shelvefn)
        return name, False, False
    try:
        s = shopen(shelvefn)
        modeldata = s['modeldata']
        readFiles = s['readFiles']
        s.close()
    except:
        print("This shelve fn doesn't work properly", shelvefn)
        return name, False, False

    #####
    # Load data
    try:
        rlmData = modeldata[(region, layer, metric)]
    except:
        return name, False, False
    #if debug:print 'analysis_level0:', rlmData

    #####
    # Load times and sorted data array.
    times = sorted(rlmData.keys())
    tdata = [rlmData[t] for t in times]
    #if debug:print 'analysis_level0:', "times:", times
    #if debug:print 'analysis_level0:', "data:", tdata

    #####
    # Deetermine the time range and take the mean of the data.
    finalData = {}
    if len(times) < 30:
        mean = np.ma.mean(tdata)
        timeRange = [times[0], times[-1]]
        if debug:
            print('analysis_level0:', "times:", times)
            print('analysis_level0:', "data:", tdata)

    else:
        mean = np.ma.mean(tdata[-30:])
        timeRange = [times[-30], times[-1]]
        if debug:
            print('analysis_level0: (last 30)', "times:", times[-30:])
            print('analysis_level0: (last 30)', "data:", tdata[-30:])

    timestr = '-'.join([str(int(y)) for y in timeRange])

    ####
    # Finish up and return fields: name, data, timerange
    if debug: print('analysis_level0:', name, float(mean), timestr)
    return name, float(mean), timestr


def analysis_level0_insitu(
    jobID='',
    field="Nitrate",
    region='regionless',
    layer='layerless',
    metric='metricless',
    debug=True,
):
    #####
    # Start and load shelve
    if debug:
        print('analysis_level0_insitu:', jobID, field, region, layer, metric)
    shelvefn = shelvedir + "/timeseries/" + jobID + "/" + jobID + "_" + field + "_insitu.shelve"

    if debug:
        print('analysis_level0_insitu:', shelvefn, os.path.exists(shelvefn))
    if not os.path.exists(shelvefn):
        print("This shelve fn doesn't exists", shelvefn)
        return False
    try:
        s = shopen(shelvefn)
        dataD = s['dataD']
        s.close()
    except:
        print("This shelve fn doesn't work properly", shelvefn)
        try:
            s = shopen(shelvefn)
            print("s.keys:", list(s.keys()))
            s.close()
        except:
            print("Can not open shelve")
        return False

    #####
    # Load data
    try:
        rlData = dataD[(region, layer)]
    except:
        print("This region/layer not in the in situ data: ", (
            jobID, field), ':', (region, layer))
        return False

    if metric == 'mean': rlmData = rlData.mean()
    if metric == 'median': rlmData = np.median(rlData)
    if metric == 'min': rlmData = np.min(rlData)
    if metric == 'max': rlmData = np.max(rlData)
    if metric == 'metricless': rlmData = rlData[0]

    if debug: print('analysis_level0_insitu:', rlmData)

    ####
    # Finish up and return fields: data
    return rlmData


def main():
    try:
        jobID = argv[1]
    except:
        jobID = "u-ad371"

    analysis_level0_insitu(jobID=jobID, )
    #analysis_level0(jobID =jobID,)


if __name__ == "__main__":
    main()
