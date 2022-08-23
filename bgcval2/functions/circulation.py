#
# Copyright 2017, Plymouth Marine Laboratory
#
# This file is part of the bgc-val library.
#
# bgc-val is free software: you can redistribute it and/or modify it
# under the terms of the Revised Berkeley Software Distribution (BSD) 3-clause license. 

# bgc-val is distributed in the hope that it will be useful, but
# without any warranty; without even the implied warranty of merchantability
# or fitness for a particular purpose. See the revised BSD license for more details.
# You should have received a copy of the revised BSD license aeORCA1_LONg with bgc-val.
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
.. module:: Circulation tools
   :platform: Unix
   :synopsis: Calculates the various current in the eORCA1 grid.
.. moduleauthor:: Lee de Mora <ledm@pml.ac.uk>
"""
import os
import numpy as np
from bgcval2.bgcvaltools.dataset import dataset
from bgcval2.UKESMpython import maenumerate
from bgcval2.functions.get_kwarg_file import get_kwarg_file

# coordinates of Drake Passage in eORCA1

eORCA1_LON=219
eORCA1_LAT0=79
eORCA1_LAT1=109

eORCA1_latslice26N = slice(227,228)
eORCA1_latslice26Nnm = slice(228,229)
eORCA1_latslice32S = slice(137,138)

umask_drake = 0
e2u_drake = 0
e3v_AMOC26N = 0
e1v_AMOC26N = 0
tmask_AMOC26N = 0
alttmask_AMOC26N = 0
loadedArea = False
loadedAltMask = False


def loadDataMask(gridfn, maskname, grid):
    """
    Load files and some key fields:
    """
    global umask_drake
    global e2u_drake
    global e3v_AMOC26N
    global e1v_AMOC26N
    global tmask_AMOC26N
    global loadedArea
    global loadedAltMask

    if grid == 'eORCA1':
        LON = eORCA1_LON
        LAT0 = eORCA1_LAT0
        LAT1 = eORCA1_LAT1
        latslice26Nnm = eORCA1_latslice26Nnm

    nc = dataset(gridfn, 'r')        
    e2u_drake = nc.variables['e2u'][LAT0:LAT1, LON]
    umask_drake = nc.variables['umask'][:, LAT0:LAT1, LON]
    e3v_AMOC26N = nc.variables['e3v'][:, latslice26Nnm, :]   # z level height 3D
    print('e3v_AMOC26N:', e3v_AMOC26N)

    e1v_AMOC26N = nc.variables['e1v'][latslice26Nnm, :]     #
    tmask_AMOC26N = nc.variables['tmask'][:, latslice26Nnm, :]    
    nc.close()
    loadedArea = True


def loadAtlanticMask(altmaskfile, maskname='tmaskatl', grid = 'eORCA1'):
    """
    Load the atlantic ocean mask.
    """
    global alttmask_AMOC26N
    if grid == 'eORCA1':
        latslice26Nnm = eORCA1_latslice26Nnm

    nc = dataset(altmaskfile, 'r')        
    alttmask_AMOC26N = nc.variables[maskname][latslice26Nnm, :]
    nc.close()
    loadedAltMask = True


def find_keys_in_nc(nc, keys):
    """
    Find some keys in a netcdf.
    """
    intersection_keys = list(set(keys) & set(nc.variables.keys()))
    if not len(intersection_keys):
        raise KeyError(f'Unable to find keys {keys} in {nc.filename}')
    if len(intersection_keys)>1:
        print('WARNING: found more than one key {intersection_keys} in {nc.filename}')
    return intersection_keys



def drakePassage(nc, keys, **kwargs):
    """
    This function calculates drake passage for eORCA1. 
    
    nc: a netcdf openned as a dataset.
    keys: a list of keys to use in this function.
    
    """
    areafile = get_kwarg_file(kwargs, 'areafile')
    maskname = kwargs.get('maskname', 'tmask')
    grid = kwargs.get('grid', 'eORCA1')

    if not loadedArea:
        loadDataMask(areafile, maskname, grid)

    if grid == 'eORCA1':
        LON = eORCA1_LON
        LAT0 = eORCA1_LAT0
        LAT1 = eORCA1_LAT1

    all_e3u_keys = ['thkcello', 'e3u']
    e3u_keys = find_keys_in_nc(nc, all_e3u_keys)
    e3u = nc.variables[e3u_keys[0]][0, :, LAT0:LAT1, LON]

    all_velo_keys = ['uo', 'u3d']
    velo_keys = find_keys_in_nc(nc, all_velo_keys)
    velo = nc.variables[velo_keys[0]][0, :, LAT0:LAT1, LON]    

    drake = np.sum(velo * e3u * e2u_drake * umask_drake) * 1.e-6
    return drake


def TwentySixNorth(nc,keys,**kwargs):
    """
    This function loads the AMOC/ADRC array that is used for eORCA
    
    nc: a netcdf openned as a dataset.
    keys: a list of keys to use in this function.
    
    """
    areafile = get_kwarg_file(kwargs, 'areafile')
    maskname = kwargs.get('maskname', 'tmask')
    grid = kwargs.get('grid', 'eORCA1')

    if not loadedArea:
        loadDataMask(areafile, maskname, grid)

    altmaskfile = get_kwarg_file(kwargs, 'altmaskfile', default = 'bgcval2/data/basinlandmask_eORCA1.nc')
     
    if not loadedAltMask: 
        loadAtlanticMask(altmaskfile, maskname='tmaskatl', grid=grid)

    if grid == 'eORCA1':
        latslice26Nnm = eORCA1_latslice26Nnm
    
    zv = np.ma.array(nc.variables[keys[0]][..., latslice26Nnm, :]) # m/s
    atlmoc = np.array(np.zeros_like(zv[0, :, :, 0]))

    for (z, la, lo), _ in np.ndenumerate(e3v_AMOC26N):
        if not alttmask_AMOC26N[la, lo]:
            continue
        if not tmask_AMOC26N[z, la, lo]:
            continue
        if np.ma.is_masked(zv[0, z, la, lo]):
                continue
        atlmoc[z, la] = atlmoc[z, la] - e1v_AMOC26N[la, lo] * e3v_AMOC26N[z, la, lo] * zv[0, z, la, lo] / 1.E06

    for z in range(e3v_AMOC26N.shape[0] -2, 1, -1): # add from the bottom up
        atlmoc[z, :] = atlmoc[z+1, :] + atlmoc[z, :]
    return atlmoc
       
 
def AMOC26N(nc, keys, **kwargs):
    atlmoc = TwentySixNorth(nc, keys, **kwargs)
    return atlmoc.max()

