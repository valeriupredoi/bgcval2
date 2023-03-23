# Copyright 2023, Plymouth Marine Laboratory
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
.. module:: tools
   :platform: Unix
   :synopsis: This uses some tools for the functions library.

.. moduleauthor:: Lee de Mora <ledm@pml.ac.uk>

"""



import numpy as np
from bgcval2.functions.standard_functions import choose_best_var


def load_area(nc):
    """
    Generic tool for loading area:
    """
    area_keys = ['area', 'area_grid_T', 'area_grid_W', 'area_grid_V', 'area_grid_U']
    if set(area_keys).intersection(set(nc.variables.keys())): 
        area = choose_best_var(nc, area_keys)
    elif set(['e1t', 'e2t']).intersection(set(nc.variables.keys())):
        area = nc.variables['e1t'][:]*nc.variables['e2t'][:]
    else:
        raise ValueError('Unable to load or calculate area from this file.')
    return area

