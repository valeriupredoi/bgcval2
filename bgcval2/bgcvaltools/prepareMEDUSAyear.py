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
.. module:: prepareMEDUSAyear
   :platform: Unix
   :synopsis: A tool for stiching together multiple months of MEDUSA data into one annual file. 
.. moduleauthor:: Lee de Mora <ledm@pml.ac.uk>

"""

from glob import glob
from os.path import basename, exists
from sys import argv
import numpy as np
from .. import UKESMpython as ukp
from netCDF4 import Dataset
from ..netcdf_manipulation.changeNC import changeNC, AutoVivification
from ..netcdf_manipulation.mergeNC import mergeNC
from ..netcdf_manipulation.pruneNC import pruneNC
"""	The goal of this code is to have a simple way to make climatology dataset.
"""


def setup(jobID, key, runType, foldIn, foldOut):

    try:
        float(key)
        yearKey = True
        print("Key is a specific year:", key)
        years = [
            key,
        ]
    except:
        assert False
        #yearKey=False
        #years = [str(y) for y in np.arange(1997,2008)]

    baseline = [
        'deptht',
        'nav_lat',
        'nav_lon',
        'time_counter',
    ]

    if runType in [
            'DIN',
            'FER',
            'SIL',
            'OXY',
            'PHD',
            'ZMI',
            'ZME',
    ]:
        keys = [
            runType,
        ]
        L = '_ptrc_T'

    if runType == 'CHL':
        keys = [
            'CHD',
            'CHN',
        ]
        L = '_ptrc_T'

    if runType == 'SAL':
        keys = [
            'vosaline',
        ]
        L = 'grid_T'

    if runType == 'TEMP':
        keys = [
            'votemper',
        ]
        L = 'grid_T'

    if runType == 'MLD':
        keys = [
            'somxl010',
        ]
        L = 'grid_T'

    if runType == 'U':
        keys = [
            'vozocrtx',
        ]
        L = 'U'

    if runType == 'V':
        keys = [
            'vomecrty',
        ]
        L = 'V'

    if runType == 'W':
        keys = [
            'vovecrtz',
        ]
        L = 'W'

    #months = sorted(['0121', '0821','0321','0921', '0421','1021', '0521','1121', '0621','1221', '0721','1221'])
    cal = '365_day'

    #if jobID[:4]=='xjez' and key in ['2001', ]:
    #	print jobID, key
    #	months = sorted([ '20010301', '20010501', '20010701', '20010901', '20011101', '20010130', '20010330', '20010530', '20010730', '20010930', '20011130','20020101',])
    #	cal = '360_day'

    mergedFiles = []

    #fns = foldIn+'/'+jobID+'*_1m_'+key+'*'+L+'*.nc'
    fns = foldIn + '/' + jobID + '*_' + key + '*' + L + '.nc'
    filesIn = sorted(glob(fns))

    print("filesIn:", fns, filesIn)
    filenameOut = ukp.folder(foldOut +
                             key) + jobID + '_' + key + '_' + runType + '.nc'
    run(jobID, key, keys, runType, filesIn, foldOut)


def run(jobID, key, keys, runType, filesIn, filenameOut):

    #filenameOut 	= ukp.folder(foldOut+key)+jobID+'_'+key+'_'+runType+'.nc'
    #ukp.folder(foldOut+key+'-annual')+jobID+'_'+key+'-annual'+'_'+runType+'.nc'
    filenameAnnual = filenameOut.replace(key, key + '-annual')

    if exists(filenameOut) and exists(filenameAnnual):
        print("Already exist:", filenameOut, '\n\tand:', filenameAnnual)
        return

    for fn in filesIn:

        prunedfn = ukp.folder('/tmp/outNetCDF/tmp-Clims') + basename(
            fn)[:-3] + '_' + key + '_' + runType + '.nc'
        print(fn, '--->', prunedfn)

        if not exists(prunedfn):
            m = pruneNC(fn, prunedfn, keys, debug=True)  #,calendar=cal)
        if runType == 'CHL':
            nc = Dataset(prunedfn, 'r')

            fileOut = prunedfn.replace('.nc', '_CHL.nc')
            if not exists(fileOut):
                av = AutoVivification()
                av['CHN']['name'] = 'False'
                av['CHD']['name'] = 'CHL'
                av['CHD']['units'] = '[mg Chl/m3]'
                av['CHD']['long_name'] = 'Total Chlorophyll'
                av['CHD']['newDims'] = ('time_counter', 'deptht', 'y', 'x')
                av['CHD']['newData'] = nc.variables['CHD'][:] + nc.variables[
                    'CHN'][:]
                nc.close()
                c = changeNC(prunedfn, fileOut, av, debug=True)
            print(fileOut)
            #print prunedfn
        else:
            fileOut = prunedfn

        mergedFiles.append(fileOut)

    if runType == 'CHL': keys = [
            'CHL',
    ]

    if ukp.shouldIMakeFile(mergedFiles, filenameOut):
        m = mergeNC(mergedFiles,
                    filenameOut,
                    keys,
                    timeAverage=False,
                    debug=True,
                    calendar=cal)

    if ukp.shouldIMakeFile(mergedFiles, filenameAnnual):
        m = mergeNC(mergedFiles,
                    filenameAnnual,
                    keys,
                    timeAverage=True,
                    debug=True,
                    calendar=cal)


def main():
    runTypes = [
        'PHD',
        'ZMI',
        'ZME',
        'SAL',
        'TEMP',
        'MLD',
        'OXY',
        'DIN',
        'CHL',
        'SIL',
        'FER',
        'W',
        'U',
        'V',
    ]

    #'ERSEMNuts','ERSEMphytoBm','ERSEMphytoChl','ERSEMzoo', 'ERSEMMisc','ERSEMbac']
    #'SalTempWind','ERSEMFull','ERSEMphyto','Detritus', ]#'SalTempWind', ]# ]#]#]

    try:
        jobID = argv[1]
        key = argv[2]
        print("Using command line arguments:", jobID, key)
    except:
        jobID = 'xjwki'
        key = '1979'
        print("Not using command line arguments")
        jobs = [
            'xjwki',
        ]
        for j in jobs:
            for r in runTypes:
                #foldIn = '/data/euryale7/scratch/ledm/iMarNet/'+j+'/MEANS/'
                foldIn = '/data/euryale7/scratch/ledm/UKESM/MEDUSA-ORCA025/' + j
                foldOut = ukp.folder(
                    '/data/euryale7/scratch/ledm/UKESM/MEDUSA/' + j +
                    '_postProc/')
                run(j, key, r, foldIn, foldOut)
        return

    for r in runTypes:
        #foldIn = '/data/euryale7/scratch/ledm/iMarNet/'+jobID+'/MEANS/'
        foldIn = '/data/euryale7/scratch/ledm/UKESM/MEDUSA/' + jobID
        foldOut = ukp.folder('/data/euryale7/scratch/ledm/UKESM/MEDUSA/' +
                             jobID + '_postProc/')
        setup(jobID, key, r, foldIn, foldOut)


if __name__ == "__main__":
    main()
