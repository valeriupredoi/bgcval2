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
.. module:: ice
   :platform: Unix
   :synopsis: This function calculates the ice fields.

.. moduleauthor:: Lee de Mora <ledm@pml.ac.uk>

"""
import os
import numpy as np
from bgcval2.bgcvaltools.dataset import dataset
from bgcval2.functions.get_kwarg_file import get_kwarg_file

global loaded_area_and_mask
global tmask
global area
global lat


tmask = 0
area = 0
lat = 0
loaded_area_and_mask = False


def load_area_and_mask(gridfn, maskname,):
    if isinstance(gridfn, list) and len(gridfn)==1:
        gridfn = gridfn[0]
    nc = dataset(gridfn, 'r')
    tmask = nc.variables[maskname][0]
    area = nc.variables['e2t'][:] * nc.variables['e1t'][:]
    lat = nc.variables['nav_lat'][:]
    nc.close()
    loaded_area_and_mask = True
    return area, tmask, lat


def calculate_ice_extent(nc, keys, **kwargs):  
    """
    Calculate the total ice extent.
    This is the total area of all model cells with more than 15% ice coverage. 
    
    kwargs:
        areafile: path to area file
        maskname: name of the mask in the input file. (default is tmask)
        minIce: fraction of ice coverage, default is 15%.
        region: global, North hemisphere or Southern Hemisphere. 
    """
    areafile = get_kwarg_file(kwargs, 'areafile')

    maskname = kwargs.get('maskname', 'tmask')
    if not loaded_area_and_mask:
        area, tmask, lat = load_area_and_mask(areafile, maskname)

    minIce = kwargs.get('minIce', 0.15)
    region = kwargs.get('region', 'Global')
    region = region.lower().replace('-', '').replace(' ', '').replace('_', '')

    # create mask:
    new_mask = (tmask==0) + (nc.variables[keys[0]][:].squeeze()<minIce)
    if region in ['n', 'north', 'northern', 'northhemisphere', 'northernhemisphere']:
        new_mask += (lat<0.)
    
    if region in ['s', 'south', 'southern', 'southhemisphere', 'southernhemisphere']:
        new_mask += (lat>0.)

    # calculate sum of unmasked area in km2. 
    return np.ma.masked_where(new_mask, area).sum() / 1E12
