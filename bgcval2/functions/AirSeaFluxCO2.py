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
.. module:: TotalAirSeaFluxCO2
   :platform: Unix
   :synopsis: This function calculated the total Air Sea FLux for the MEDUSA model in the eORCA grid.

.. moduleauthor:: Lee de Mora <ledm@pml.ac.uk>

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
    nc = dataset(areafile, 'r')
    area = nc.variables['e2t'][:] * nc.variables['e1t'][:]
    nc.close()
    loaded_area_and_mask = True
    return area


def calc_total_airseafluxco2(nc, keys, **kwargs):
    """
    This function calculated the total Air Sea FLux for the MEDUSA model in the eORCA grid.
    """
    areafile = get_kwarg_file(kwargs, 'areafile')

    if not loaded_area_and_mask:
        area = load_area_and_mask(areafile)

    factor =  365.25 * 12. / 1000. / 1.E15
    try:   
        arr = nc.variables[keys[0]][:].squeeze() * factor    # mmolC/m2/d
    except: 
        raise ValueError(f"TotalAirSeaFluxCO2: Not able to load {keys[0]}) from {nc.filename}")
    
    if arr.ndim ==2: 
        arr = arr * area
    else: 
        raise ValueError(f"TotalAirSeaFluxCO2: {keys[0]} from {nc.filename} has an unexpected number of dimensions: {arr.ndim}")

    return arr.sum()


def calc_ma_total_airseafluxco2(nc, keys, **kwargs):
    """
    This function calculated the total Air Sea FLux for the ERSEM model in the eORCA025 grid.
    Mission Atlantic code
    """
    factor =  365.25 * 12. / 1000. / 1.E15 # convert mmolC/m2/d to Pg/yr
    if  keys[0] not in nc.variables: return 0.
    arr = nc.variables[keys[0]][:].squeeze() * nc.variables['area'][:] *factor    # mmolC/m2/d
    return arr.sum()
    
    
