#!/usr/bin/ipython
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
"""
.. module:: mergeMonthlyFiles
   :platform: Unix
   :synopsis: A tool for stiching together multiple months of model data into one annual file. 
.. moduleauthor:: Lee de Mora <ledm@pml.ac.uk>

"""
import os
from re import findall

from .. import UKESMpython as ukp
from ..netcdf_manipulation.mergeNC import mergeNC


def getYearFromFile(fn):
    """ 
	Takes a file anem, and looks for 8 consequetive numbers, then removes those that are months, and returns the year.
	"""
    a = findall(r'\d\d\d\d\d\d\d\d', fn)
    mns = [ukp.mnStr(m) for m in range(1, 13)]
    for i in a:
        if i[-2:] in ['28', '29', '30', '31', '01'] and i[-4:-2] in mns:
            yr = i[:4]
            return yr

    return False


def getDatesFromFilename(fn):
    """ 
	Takes a file anem, and looks for 8 consequetive numbers, then removes those that are months, and returns the year.
	"""
    a = findall(r'\d\d\d\d\d\d\d\d', fn)
    dates = []
    for i in a:
        i = int(i)
        if i < 10000101: continue
        if i > 99999999: continue
        dates.append(i)
    return dates


def getMonthsFromFilename(fn):
    """ 
	Takes a file anem, and looks for 8 consequetive numbers, then removes those that are months, and returns the year.
	"""
    a = findall(r'\d\d\d\d\d\d\d\d', fn)
    #dates = []
    year = getYearFromFile(fn)
    mns = [ukp.mnStr(m) for m in range(1, 13)]
    months = []
    for i in a:
        mn = i[4:6]
        if mn not in mns: continue
        months.append(mn)
    return months


def getAnnualFilename(files, outfolder, year):
    files = sorted(files)
    #####
    # determinig filename on the way out
    if outfolder == '':
        outfolder = ukp.folder(os.path.dirname(files[0]) + '/Annual')

    mintime = ''
    maxtime = ''

    #####
    # Fiund out the range of file names
    for f in files:
        f = os.path.basename(f)
        startstop = findall(r'\d\d\d\d\d\d\d\d', f)
        if mintime == '': mintime = min(startstop)
        if maxtime == '': maxtime = max(startstop)
        if min(startstop) < mintime: mintime = min(startstop)
        if max(startstop) > maxtime: maxtime = max(startstop)

    #####
    # Create a filename that reflects the new times.
    basefile = 'Annual-' + os.path.basename(files[0])
    startstop = findall(r'\d\d\d\d\d\d\d\d', basefile)
    if len(startstop) == 2:
        basefile = basefile.replace(startstop[0], mintime)
        basefile = basefile.replace(startstop[1], maxtime)

    filenameOut = outfolder + basefile
    return filenameOut


def mergeMonthlyFiles(files,
                      outfolder='',
                      cal='360_day',
                      timeAverage=False,
                      expectedNumberOfFiles=12):
    #####
    # This assuemd that the files have already been split up using the moo filter tool
    # done in the the bgcvalTools/downloadFromMass.py

    filesOut = []
    years = {}

    #####
    # Load file
    for fn in sorted(files):
        yr = getYearFromFile(fn)
        try:
            years[yr].append(fn)
        except:
            years[yr] = [
                fn,
            ]

    #####
    #
    for yr in sorted(years.keys()):
        yearFiles = sorted(years[yr])
        if len(yearFiles) != expectedNumberOfFiles:
            print("Not enough files in ", yr, len(years[yr]), 'expecting:',
                  expectedNumberOfFiles)
            continue

        filenameOut = getAnnualFilename(yearFiles, outfolder, yr)

        if ukp.shouldIMakeFile(yearFiles, filenameOut):
            m = mergeNC(years[yr],
                        filenameOut, [],
                        timeAverage=timeAverage,
                        debug=True,
                        calendar=cal)

        filesOut.append(filenameOut)
    return filesOut


def meanDJF(
    files,
    outfolder='',
    cal='360_day',
    timeAverage=False,
):
    #####
    # This assuemd that the files have already been split up using the moo filter tool
    # done in the the bgcvalTools/downloadFromMass.py
    filesOut = []
    years = {}

    #####
    # Load file
    for fn in sorted(files):
        dates = getDatesFromFilename(fn)
        year = getYearFromFile(fn)
        months = getMonthsFromFilename(fn)
        if months in [
            ['12', '01'],
            ['01', '12'],
        ]:
            yrstr = str(int(year) + 1) + '_DJF'
        else:
            yrstr = str(int(year)) + '_DJF'

        try:
            years[yrstr].append(fn)
        except:
            years[yrstr] = [
                fn,
            ]

    #####
    for yr in sorted(years.keys()):
        yearFiles = sorted(years[yr])
        if len(yearFiles) != 3:
            print("Not enough files in ", yr, len(years[yr]))
            continue

        filenameOut = getAnnualFilename(yearFiles, outfolder, yr)
        if ukp.shouldIMakeFile(yearFiles, filenameOut):
            m = mergeNC(years[yr],
                        filenameOut, [],
                        timeAverage=timeAverage,
                        debug=True,
                        calendar=cal)

        filesOut.append(filenameOut)
    return filesOut
