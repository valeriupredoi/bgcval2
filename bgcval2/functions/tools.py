

import numpy as np
from bgcval2.functions.standard_functions import choose_best_var


def load_area(nc):
    """
    Generic tool for loading area:
    """
    area_keys = ['area', 'area_grid_T', 'area_grid_W', 'area_grid_V', 'area_grid_U']
    if set(area_keys).intersection(set(nc.variables.keys()): 
         return choose_best_var(nc, area_keys)
    else:
        area = nc.variables['e1t'][:]*nc.variables['e2t'][:]


