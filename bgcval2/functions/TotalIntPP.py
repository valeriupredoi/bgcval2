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
.. module:: TotalIntPP
   :platform: Unix
   :synopsis: This function calculated the total primary production for the MEDUSA model in the eORCA grid.
.. moduleauthor:: Lee de Mora <ledm@pml.ac.uk>
"""

import numpy as np
from bgcval2.bgcvaltools.dataset import dataset
from bgcval2.functions.get_kwarg_file import get_kwarg_file
from bgcval2.functions.standard_functions import find_best_var

global loadedArea
global model_area

model_area = 0
loadedArea = False


def loadDataMask(gridfn):
    nc = dataset(gridfn,'r')
    try:
        model_area = nc.variables['area'][:]
    except:
         model_area = nc.variables['e1t'][:]*nc.variables['e2t'][:]
    nc.close()
    loadedArea = True
    return model_area


def TotalIntPP(nc, keys, **kwargs):
    """
    This function calculated the total primary production for the MEDUSA model in the eORCA grid.
    """
    areafile = get_kwarg_file(kwargs, 'areafile')

    if not loadedArea:
        model_area = loadDataMask(areafile)

    #    mmolN/m2/d        [mg C /m2/d]   [mgC/m2/yr] [gC/m2/yr]     Gt/m2/yr # MEDUSA UNITS
    factor = 1.     * 6.625 * 12.011 * 365.       / 1000.   /     1E15
    arr = (nc.variables[keys[0]][:]+ nc.variables[keys[1]][:]).squeeze()*factor
    if arr.ndim == 3:
        for i in np.arange(arr.shape[0]):
            arr[i] = arr[i]*model_area
    elif arr.ndim == 2:
        arr = arr*model_area
    else:
        raise ValueError(f'TotalIntPP: arr shape not recognised, should be either 2 or 3, got value {arr.shape} in {nc.filename}')
    return arr.sum()


def MA_TotalIntPP(nc, keys, **kwargs):
    """
    This function calculated the total primary production for the ERSEM  model.
    """
    # mg C/m^3/d (supposedly - but possibly /m2 ?)
    factor = 365.25 / 1000. / 1E15
    area_key = find_best_var(nc, ['area', 'area_grid_T'])
    area = nc.variables[area_key][:]
    
    # This will only work for NEMO models:
    thkfn = nc.filename.replace('diag_T', 'grid_T').replace('ptrc_T', 'grid_T')
    print('MA_TotalIntPP: opening:', thkfn)
    thknc = dataset(thkfn, 'r')
    thick = thknc.variables['thkcello'][:]
    thknc.close()
    
    arr = thick * nc.variables[keys[0]][:]*factor
    arr = arr.sum(axis=1)
    arr = arr * area[None, ...]
    return arr.sum()

