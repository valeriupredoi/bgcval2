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


tmask = {}


def loadDataMask(gridfn, maskname):
    global tmask
    nc = dataset(gridfn, 'r')        
    tmask[(gridfn, maskname)] = nc.variables[maskname][:].squeeze()
    nc.close()


def applyLandMask(nc, keys, **kwargs):
    try:
        gridfn = kwargs['areafile']
    except:
        raise AssertionError("applyLandMask:\t Needs an `areafile` kwarg to apply mask")    
    
    try:
        maskname = kwargs['maskname']
    except:
        maskname = 'tmask'

    try:
        mask = tmask[(gridfn, maskname)]
    except:
        loadDataMask(gridfn, maskname)
        mask = tmask[(gridfn, maskname)]
    
    try:
        shape = mask.shape
    except:
        raise AssertionError("applyLandMask.py:\t Model mask not loaded correctly")

    arr = np.ma.array(nc.variables[keys[0]][:]).squeeze()
    m   = np.ma.masked_where(mask + arr.mask, arr)
    m  += np.ma.masked_invalid(arr).mask
    arr = np.ma.masked_where(m>1, arr)
    print("applyLandMask:\t ", arr.mean(), arr.min(), arr.max(), arr.shape, m.mean())
    assert 0
    return arr


def applyLandMask2D(nc, keys, **kwargs):
    """
    Useful for when the mask is 3D, but the field is only 2D,
    so you want to apply the surface layer of the mask.
    """
    try:
        gridfn = kwargs['areafile']
    except:
        raise AssertionError("applyLandMask:\t Needs an `areafile` kwarg to apply mask")    
    
    try:
        maskname = kwargs['maskname']
    except:
        maskname = 'tmask'

    try:
        mask = tmask[(gridfn, maskname)][0]
    except:
        loadDataMask(gridfn, maskname)
        mask = tmask[(gridfn, maskname)][0]
            
    try:
        shape = mask.shape
    except:
        raise AssertionError("applyLandMask.py:\t Model mask not loaded correctly")

    arr = np.ma.array(nc.variables[keys[0]][:]).squeeze()
    m  = np.ma.masked_where(mask + arr.mask, arr)
    m += np.ma.masked_invalid(arr).mask
    return np.ma.masked_where(m, arr)


def applyLandMask_maskInFile(nc, keys, **kwargs):
    """
    Useful for when the mask is already in the same netcdf.
    """
    
    try:
        maskname = kwargs['maskname']
    except:
        maskname = 'tmask'

    try:
        mask = tmask[(nc.filename, maskname)]
    except:
        loadDataMask(nc.filename, maskname)
        mask = tmask[(nc.filename, maskname)]
    
    print("applyLandMask_maskInFile", nc.filename)
    
    try:
        shape = mask.shape
    except:
        raise AssertionError("applyLandMask_maskInFile.py:\t Model mask not loaded correctly")

    arr = np.ma.array(nc.variables[keys[0]][:]).squeeze()
    m = np.ma.masked_where(mask + arr.mask, arr)
    m += np.ma.masked_invalid(arr).mask
    return np.ma.masked_where(m, arr)
    


