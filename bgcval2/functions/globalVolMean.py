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
.. module:: globalVolMean
   :platform: Unix
   :synopsis: This function calculates the volume weighted mean of a field.

.. moduleauthor:: Lee de Mora <ledm@pml.ac.uk>

"""

import numpy as np
from bgcval2.functions.standard_functions import choose_best_var 
from bgcval2.functions.tools import load_area

def calc_vol(nc):
    """
    Calculate volume from the (grid-T) netcdf file.
    """
    area = load_area(area)
    thkcello = nc.variables['thkcello'][:]
    if area.ndim == 4:
        area = area[None, None,:, :]
    if area.ndim == 3:
        area = area[None, :, :]
    return thkcello * area


def globalVolumeMean(nc, keys, **kwargs):
    """
    Calculate the global volume mean.
    """
    # first check kwargs
    if not kwargs:
        raise FileNotFoundError("No keyword args supplied to function.")

    # To add a constant value to the data (usually Kelvin to Celcius)
    addvalue = kwargs.get('addvalue', 0.)

    # Multiply the data by some factor (ie to change units)
    multiplyBy = kwargs.get('multiplyBy', 1.)

    temp = choose_best_var(nc, keys)
    temp = np.ma.masked_where((temp == 0.) + temp.mask, temp)
    volume = calc_vol(nc)
    volume = np.ma.masked_where(volume.mask + temp.mask, volume)
    temp = temp * multiplyBy + addvalue

    if temp.shape == volume.shape:
        vol = np.ma.masked_where(temp.mask, volume)
        return (temp*volume).sum()/(volume.sum())

    raise ValueError(f"Data or mask has unexpected shape {temp.shape}, {volume.shape}")


