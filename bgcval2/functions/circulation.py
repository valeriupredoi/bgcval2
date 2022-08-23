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

LON=219
LAT0=79
LAT1=109

latslice26N = slice(227,228)
latslice26Nnm = slice(228,229)
latslice32S = slice(137,138)

umask_drake = 0
e2u_drake = 0
e3v_AMOC26N = 0
e1v_AMOC26N = 0
tmask_AMOC26N = 0
alttmask_AMOC26N = 0
loadedArea = False
loadedAltMask = False


def loadDataMask(gridfn, maskname):
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

    nc = dataset(gridfn, 'r')        
    e2u_drake = nc.variables['e2u'][LAT0:LAT1, LON]
    umask_drake = nc.variables['umask'][:, LAT0:LAT1, LON]
    e3v_AMOC26N = nc.variables['e3v'][:, latslice26Nnm, :]   # z level height 3D
    print('e3v_AMOC26N:', e3v_AMOC26N)

    e1v_AMOC26N = nc.variables['e1v'][latslice26Nnm, :]     #
    tmask_AMOC26N = nc.variables['tmask'][:, latslice26Nnm, :]    
    nc.close()
    loadedArea = True


def loadAtlanticMask(altmaskfile, maskname='tmaskatl'):
    global alttmask_AMOC26N

    nc = dataset(altmaskfile, 'r')        
    alttmask_AMOC26N = nc.variables[maskname][latslice26Nnm, :]
    nc.close()
    loadedAltMask = True


def drakePassage(nc, keys, **kwargs):
    """
    This function. 
    
    nc: a netcdf openned as a dataset.
    keys: a list of keys to use in this function.
    
    """
    areafile = get_kwarg_file(kwargs, 'areafile')

    maskname = kwargs.get('maskname', 'tmask')
            
    if not loadedArea:
        loadDataMask(areafile, maskname)

    try:
        e3u = nc.variables['thkcello'][0, :, LAT0:LAT1, LON]    
    except:
        e3u = nc.variables['e3u'     ][0, :, LAT0:LAT1, LON]

    try:
        velo = nc.variables['uo' ][0, :, LAT0:LAT1, LON]    
    except:
        velo = nc.variables['u3d'][0, :, LAT0:LAT1, LON]

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

    if not loadedArea:
        loadDataMask(areafile, maskname)

    altmaskfile = get_kwarg_file(kwargs, 'altmaskfile',  default = 'bgcval2/data/basinlandmask_eORCA1.nc')
     
    if not loadedAltMask: 
        loadAtlanticMask(altmaskfile, maskname='tmaskatl')
    
    zv = np.ma.array(nc.variables[keys[0]][..., latslice26Nnm, :]) # m/s
    atlmoc = np.array(np.zeros_like(zv[0, :, :, 0]))
    print('e3v_AMOC26N:', e3v_AMOC26N)
    e2vshape = e3v_AMOC26N.shape
    xsectArea = (e1v_AMOC26N * e3v_AMOC26N)
    print('TwentySixNorth:', e1v_AMOC26N.shape, e3v_AMOC26N.shape, xsectArea.shape, xsectArea.sum())

    TotalXsection = 0
    for la in range(e2vshape[1]):           #ji, y
        for lo in range(e2vshape[2]):         #jj , x,
            if int(alttmask_AMOC26N[la, lo]) == 0:
                continue
            for z in range(e2vshape[0]):        # jk
                if int(tmask_AMOC26N[z, la, lo]) == 0:
                    continue
                if np.ma.is_masked(zv[0, z, la, lo]):
                    continue
                atlmoc[z, la] = atlmoc[z, la] - e1v_AMOC26N[la, lo]*e3v_AMOC26N[z, la, lo]*zv[0, z, la, lo]/1.E06
                TotalXsection += e1v_AMOC26N[la, lo]*e3v_AMOC26N[z, la, lo]

    print("TotalXsection:", TotalXsection)
    # Cumulative sum from the bottom up.
    for z in range(73, 1, -1):
        atlmoc[z, :] = atlmoc[z+1, :] + atlmoc[z, :]
    return atlmoc
       
 
def AMOC26N(nc, keys, **kwargs):
    atlmoc = TwentySixNorth(nc, keys, **kwargs)
    return atlmoc.max()

