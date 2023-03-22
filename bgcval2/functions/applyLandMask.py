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
.. module:: applyLandMask
   :platform: Unix
   :synopsis: This function applies a land mask. Note that this will not work for p2p code.

.. moduleauthor:: Lee de Mora <ledm@pml.ac.uk>

"""

import numpy as np
from bgcval2.bgcvaltools.dataset import dataset
from bgcval2.functions.get_kwarg_file import get_kwarg_file
from bgcval2.functions.standard_functions import choose_best_var

tmask = {}


def loadDataMask(areafile, maskname):
    global tmask
    nc = dataset(areafile, 'r')        
    tmask[(areafile, maskname)] = np.array(nc.variables[maskname][:].squeeze(), dtype=bool)
    nc.close()
    return tmask[(areafile, maskname)]


def applyLandMask(nc, keys, **kwargs):
    """
    Useful for when the mask is 3D, but the field is only 2D,
    so you want to apply the surface layer of the mask.
    """
    areafile = get_kwarg_file(kwargs, 'areafile')
    maskname = kwargs.get('maskname', 'tmask')
    mask = tmask.get((areafile, maskname), loadDataMask(areafile, maskname))
        
    #for key in keys:
    #    if key not in nc.variables.keys():
    #        print(f'key {key} not in file {nc.filename}')
    #        print('Available keys:options:', nc.variables.keys())
    #        raise KeyError(f'applyLandMask: key {key} not in file {nc.filename}.')
    
    print(f"loading variable {keys[0]} from file {nc.filename}")
    arr = choose_best_var(nc, keys).squeeze()
    new_mask = np.ma.masked_invalid(arr).mask + arr.mask
    if arr.ndim == 2 and mask.ndim == 3:
        new_mask += mask[0]
    else:
        new_mask += mask
    return np.ma.masked_where(new_mask, arr)


