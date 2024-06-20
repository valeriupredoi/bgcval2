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

from shelve import open as shOpen
from glob import glob

from bgcval2._runtime_config import get_run_configuration
from bgcval2.Paths.paths import paths_setter
from bgcval2.bgcvaltools import bv2tools as bvt


def removeFromShelves(fn, removeRegions):
    print('removing:', removeRegions, 'from', fn)
    with shOpen(fn) as sh:
        modeldata = sh['modeldata']
        for key in list(modeldata.keys()):
            try:
                (r, l, m) = key
            except ValueError:
                continue
            if r in removeRegions:
                print('modeldata[', (r, l, m), '] will be deleted')
                if (r, l, m) in modeldata:
                    del modeldata[(r, l, m)]

        sh['modeldata'] = modeldata
    # sh.close()



def main(config_user=None):

    # get runtime configuration
    if config_user:
        paths_dict, config_user = get_run_configuration(config_user)
    else:
        paths_dict, config_user = get_run_configuration("defaults")

    # filter paths dict into an object that's usable below
    paths = paths_setter(paths_dict)
        
    removeRegions = [
    #    'Remainder',
        'STNA', 'SubtropicNorthAtlantic',
        'Pitcairn',
    ]  #'ignoreInlandSeas',


    jobIDs = [
        'u-cs495', 'u-cs568', 'u-cy623', 'u-da914', 'u-da916', 
        'u-da917', 'u-cy690', 'u-cy691', 'u-cy692', 'u-cy693', 
        'u-cz152', 'u-cz014', 'u-cx209', 'u-cw988', 'u-cw989', 
        'u-cw990', 'u-cz826', 'u-cy837', 'u-cy838', 'u-cz374', 
        'u-cz375', 'u-cz376', 'u-cz377', 'u-cz378', 'u-cz834', 
        'u-cz855', 'u-cz859', 'u-db587', 'u-db723', 'u-db731', 
        'u-da087', 'u-da266', 'u-db597', 'u-db733', 'u-dc324', 
        'u-cz944', 'u-da800', 'u-da697', 'u-da892', 'u-db223', 
        'u-df453', 'u-dc251', 'u-dc051', 'u-dc052', 'u-dc248', 
        'u-dc249', 'u-db956', 'u-dc565', 'u-dd210', 'u-dc032', 
        'u-df028', 'u-dc123', 'u-dc130', 'u-df025', 'u-df027', 
        'u-df021', 'u-df023', 'u-de943', 'u-de962', 'u-de963', 
        'u-dc163', 'u-df028', ]
    for jobID in jobIDs:
        shelvedir = bvt.folder([paths.shelvedir, "timeseries", jobID])

        for fn in glob(shelvedir+'/*shelve.dir'):
            if fn.find('insitu') > -1: continue
            fn = fn[:-4]
            removeFromShelves(fn, removeRegions)

if __name__ == "__main__":
    main()
