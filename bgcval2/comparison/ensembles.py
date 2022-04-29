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
import numpy as np


def ensembleMean_PI(
        times,
        data,
        startingYears=[1850., 1880., 1960., 1922., 2020., 2050., 1995.]):
    #####
    # This is for the ensemble mean Colin requested.
    # Basically it is a mean of the same job with 7 versions.
    # However, this is run after shiftimes.
    #times, data, are equal lenth arrays/lists

    outputTimeRange = np.arange(1840, 2030)
    outTimes = {t: [] for t in outputTimeRange}
    outData = {t: [] for t in outputTimeRange}

    for startyear in startingYears:
        for outyr in outputTimeRange:
            for t, d in zip(
                    times,
                    data,
            ):
                yr = int(t)
                # trange =  yr - startyear + 1850.
                # if trange < outyr < trange + 1.:
                shifted = t - startyear + 1850.
                if outyr < shifted < outyr + 1.:
                    outTimes[outyr].append(shifted)
                    outData[outyr].append(d)
            #print 'FOUND ONE', t, yr, startyear, shifted

    for o, arr in list(outTimes.items()):
        #print 'times', o ,arr, np.ma.mean(arr), outData[o], np.ma.mean(outData[o])
        if len(arr) == 0: continue
        if len(outData[o]) == 0: continue
        outTimes[o] = np.ma.mean(arr)
        outData[o] = np.ma.mean(outData[o])

    finalTimes = []
    finalData = []
    for outyr in sorted(outputTimeRange):
        if outTimes[outyr] and outData[outyr]:
            finalTimes.append(outTimes[outyr])
            finalData.append(outData[outyr])
    #print len(finalTimes), len(finalData)
    if len(finalTimes) != len(finalData): assert 0
    #print finalTimes, finalData
    #assert 0
    return finalTimes, finalData


def build_ensemble(timesD, arrD, ensembles={}):

    #####
    # Return no change if no ensembles requested
    if ensembles == {}: return timesD, arrD

    #####
    # Determine time range
    timeRange = {}
    for jobID in sorted(timesD.keys()):
        try:
            minT = np.min(timesD[jobID])
            maxT = np.max(timesD[jobID])
        except:
            continue
        timeRange[minT] = True
        timeRange[maxT] = True

    if len(list(timeRange.keys())) == 0:
        return timesD, arrD

    timeRange = [
        np.min(list(timeRange.keys())),
        np.max(list(timeRange.keys()))
    ]
    allyears = np.arange(int(timeRange[0]), int(timeRange[1]) + 1)

    #####
    # Make empty dicts:
    newTimes = {}
    newArr = {}
    outTimes = {}
    outArr = {}
    for name, ensemble in list(ensembles.items()):
        newTimes[name] = {yr: [] for yr in allyears}
        newArr[name] = {yr: [] for yr in allyears}
        outTimes[name] = []
        outArr[name] = []

    #####
    # Load the data
    for name, ensemble in list(ensembles.items()):
        for jobID in sorted(timesD.keys()):
            if jobID not in ensemble: continue
            for t, d in zip(timesD[jobID], arrD[jobID]):
                yr = int(t)
                newTimes[name][yr].append(t)
                newArr[name][yr].append(d)

    #####
    # Take the mean of the ensemble
    for name, ensemble in list(ensembles.items()):
        for yr in sorted(allyears):
            if len(newTimes[name][yr]) == 0: continue
            if len(newArr[name][yr]) == 0: continue

            outTimes[name].append(np.ma.mean(newTimes[name][yr]))
            outArr[name].append(np.ma.mean(newArr[name][yr]))

    #####
    # Add additional years into PI control.
    if 'PI Control' in list(
            ensembles.keys()) and 'u-aw310' in ensembles['PI Control']:
        jobID = 'u-aw310'
        name = 'PI Control'
        outTimes[name], outArr[name] = ensembleMean_PI(outTimes[name],
                                                       outArr[name])
    if 'PI Control 7' in list(
            ensembles.keys()) and 'u-aw310' in ensembles['PI Control 7']:
        jobID = 'u-aw310'
        name = 'PI Control 7'
        outTimes[name], outArr[name] = ensembleMean_PI(outTimes[name],
                                                       outArr[name])

    if 'PI Control 3' in list(
            ensembles.keys()) and 'u-aw310' in ensembles['PI Control 3']:
        jobID = 'u-aw310'
        name = 'PI Control 3'

        times, arr = ensembleMean_PI(outTimes[name],
                                     outArr[name],
                                     startingYears=[2020., 2050., 1995.])
        times = np.array(times)
        arr = np.array(arr)
        times = np.ma.masked_where((times < 1850.) + (times > 1920.), times)
        arr = np.ma.masked_where((times < 1850.) + (times > 1920.), arr)
        outTimes[name], outArr[name] = times.compressed(), arr.compressed()

    if 'PI Control 6' in list(
            ensembles.keys()) and 'u-aw310' in ensembles['PI Control 6']:
        jobID = 'u-aw310'
        name = 'PI Control 6'
        times, arr = ensembleMean_PI(
            outTimes[name],
            outArr[name],
            startingYears=[1850., 1922., 1960., 2020., 2050., 1995.])
        times = np.array(times)
        arr = np.array(arr)
        times = np.ma.masked_where((times < 1920.), times)
        arr = np.ma.masked_where((times < 1920.), arr)
        outTimes[name], outArr[name] = times.compressed(), arr.compressed()

    return outTimes, outArr
