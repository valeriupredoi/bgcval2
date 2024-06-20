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
    return area.copy(), tmask, lat.copy()


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
    if 'area' not in nc.variables and not loaded_area_and_mask:
        area, tmask, lat = load_area_and_mask(areafile, maskname)
    else:
        area = nc.variables['area'][:]
        lat = nc.variables['nav_lat'][:]
        tmask = nc.variables[keys[0]][:].squeeze().mask

    print('tmask', tmask.min(), tmask.max(), tmask.sum())
    minIce = kwargs.get('minIce', 0.15)
    region = kwargs.get('region', 'Global')
    region = region.lower().replace('-', '').replace(' ', '').replace('_', '')

    # create mask:
    masked_area = np.ma.masked_where(tmask + (nc.variables[keys[0]][:].squeeze() < minIce), area)

    if region in ['n', 'north', 'northern', 'northhemisphere', 'northernhemisphere']:
        masked_area = np.ma.masked_where(masked_area.mask + (lat<0.), masked_area) 
    
    if region in ['s', 'south', 'southern', 'southhemisphere', 'southernhemisphere']:
        masked_area = np.ma.masked_where(masked_area.mask + (lat>0.), masked_area)

    return masked_area.sum()/1.E12


def calculate_ice_area(nc, keys, **kwargs):
    """
    Calculate the total ice area.
    This is the total area of all sea ice. 

    kwargs:
        areafile: path to area file
        maskname: name of the mask in the input file. (default is tmask)
        #minIce: fraction of ice coverage, default is 15%.
        region: global, North hemisphere or Southern Hemisphere.
    """
    areafile = get_kwarg_file(kwargs, 'areafile')

    maskname = kwargs.get('maskname', 'tmask')
    if 'area' not in nc.variables and not loaded_area_and_mask:
        area, tmask, lat = load_area_and_mask(areafile, maskname)
    else:
        area = nc.variables['area'][:]
        lat = nc.variables['nav_lat'][:]
        tmask = nc.variables[keys[0]][:].squeeze().mask
    ice_fraction = nc.variables[keys[0]][:].squeeze()
 
    if ice_fraction.max() > 1. or ice_fraction.min() < 0. :
        raise ValueError('Ice coverage fraction incorrect range. Are you using the correct variable?')

    region = kwargs.get('region', 'Global')
    region = region.lower().replace('-', '').replace(' ', '').replace('_', '')

    # create mask of ice values:
    masked_area = np.ma.masked_where(tmask, ice_fraction * area)

    if region in ['n', 'north', 'northern', 'northhemisphere', 'northernhemisphere']:
        masked_area = np.ma.masked_where(masked_area.mask + (lat<0.), masked_area)

    if region in ['s', 'south', 'southern', 'southhemisphere', 'southernhemisphere']:
        masked_area = np.ma.masked_where(masked_area.mask + (lat>0.), masked_area)

    return masked_area.sum()/1.E12


