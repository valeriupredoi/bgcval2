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
# This code removes a specific mask from the shelve.
# Untested!
"""
.. module:: removeMaskFromshelve
   :platform: Unix
   :synopsis: Removes data from a shelve file. 
.. moduleauthor:: Lee de Mora <ledm@pml.ac.uk>

"""

from shelve import open as shopen
from glob import glob


def removeFromShelves(fn, removeRegions):
    print('removing:', removeRegions, 'from', fn)
    sh = shopen(fn)

    modeldata = sh['modeldata']

    for key in list(modeldata.keys()):
        try:
            (r, l, m) = key
        except:
            continue
        if r in removeRegions:
            print('modeldata[', (r, l, m), '] will be deleted')
            del modeldata[(r, l, m)]

    sh['modeldata'] = modeldata
    sh.close()


removeRegions = [
    'Remainder',
]  #'ignoreInlandSeas',

for fn in glob('shelves/timeseries/u-ab749/u-ab749_*'):
    if fn.find('insitu') > -1: continue
    removeFromShelves(fn, removeRegions)
