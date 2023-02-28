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
.. module:: TotalInvtAlk
   :platform: Unix
   :synopsis: This function calculates the total quantity of dissolved alkalinity for the MEDUSA model on the eORCA grid.

.. moduleauthor:: Andrew Yool <axy@noc.ac.uk>

"""

import numpy as np
from bgcval2.bgcvaltools.dataset import dataset
from bgcval2.functions.get_kwarg_file import get_kwarg_file

global loaded_area_and_mask
global tmask
global area

area = 0
loaded_area_and_mask = False


def load_area_and_mask(areafile):
    """Load area file and its mask."""
    if isinstance(areafile, list) and len(areafile)==1:
        areafile = areafile[0]
    else:
        raise ValueError(f"{areafile} must be a length=1 list containing a single file descriptor.")
    nc = dataset(areafile, 'r')
    area = nc.variables['e2t'][:] * nc.variables['e1t'][:]
    nc.close()
    loaded_area_and_mask = True
    return area


def calc_total_invtalk(nc, keys, **kwargs):
    """
    This function calculates the total quantity of dissolved alkalinity for the MEDUSA model on the eORCA grid.
    """
    areafile = get_kwarg_file(kwargs, 'areafile')

    if not loaded_area_and_mask:
        area = load_area_and_mask(areafile)

    factor =  1. / 1.E18   # meq / m2 -> Peq
    try:   
        arr = nc.variables[keys[0]][:].squeeze() * factor    # meq / m2 -> Peq
    except: 
        raise ValueError(f"TotalInvtAlk: Not able to load {keys[0]}) from {nc.filename}")
    
    if arr.ndim ==2: 
        arr = arr * area
    else: 
        raise ValueError(f"TotalInvtAlk: {keys[0]} from {nc.filename} has an unexpected number of dimensions: {arr.ndim}")

    return arr.sum()
