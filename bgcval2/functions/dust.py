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
.. module:: calcExportRatio 
   :platform: Unix
   :synopsis: This function calculates the export ratio for the MEDUSA model in the eORCA grid.

.. moduleauthor:: Lee de Mora <ledm@pml.ac.uk>

"""

import numpy as np
from bgcval2.bgcvaltools.dataset import dataset
from bgcval2.functions.get_kwarg_file import get_kwarg_file


masked_area = 0
loadedArea = False


def loadDataMask(gridfn):
    global loadedArea
    global masked_area
    nc = dataset(gridfn,'r')
    masked_area = nc.variables['e2t'][:] * nc.variables['e1t'][:]*nc.variables['tmask'][0]
    nc.close()
    loadedArea = True


def modeldustsum(nc, keys, **kwargs):
    #factors are:
    # 1e-12: convert from mol to Gmol

    areafile = get_kwarg_file(kwargs, 'areafile')

    if not loadedArea:
        loadDataMask(areafile)

    dust = nc.variables[keys[0]][:]
    # fix broken parts
    dust[:,234:296,295:348] = 0.
    dust[:,234:248,285:295] = 0.
    dust[:,228:256,290:304] = 0.
    return (masked_area*dust).sum() *1.E-12 *365.25
