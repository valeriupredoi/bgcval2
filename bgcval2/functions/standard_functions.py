#
# Copyright 2014, Plymouth Marine Laboratory
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
.. module:: standard_functions
   :platform: Unix
   :synopsis: This module is a dictionairy of functions that can be applied to the data.
       The expectation is that users will either add to this list.
.. moduleauthor:: Lee de Mora <ledm@pml.ac.uk>

"""
import numpy as np
from bgcval2.bgcvaltools.dataset import dataset


#####
# This is the algorithm that applied the evaluation functions below.
def extractData(nc, details, key=['',], debug=False):
    """
    This extracts the data the the netcdf based on the instructions from details dictionairy.
    If you want to do something funky to the data before plotting it, just create a new convert function in this directory.
    Exampple details dict:
        {'name': 'Chlorophylla', 'vars':['Chlorophylla',], 'convert': std_functions['div1000'],'units':'ug/L'}
    """
    if isinstance(details, dict):
          if debug: print(f"extractData: details is a dict of keys {str(details.keys())}")
    else:
         raise TypeError("std_functions:\textractData:\t Details type is not a dictionary!")

    #####
    # Finding key word arguments (kwargs):
    # We define these as extra items in the details dict, and they can literally be anything
    # that the convert function needs.
    kwargs={}
    for k, value in details.items():
        if k in ['name','vars','convert', 'units']:
            continue
        kwargs[k] = details[k]

    #####
    # Apply pre-specified convert function.
    xd = details['convert'](nc,details['vars'], **kwargs)

    #####
    # Return masked array of the loaded data.
    return np.ma.array(xd)


####
# Some functions for maniulating data:
def NoChange(nc,keys):
    """
    Loads keys[0] from the netcdf, but applies no change.
    """
    return nc.variables[keys[0]]


def N2Biomass(nc,keys):
    """
    Loads keys[0] from the netcdf, but multiplies by 79.572 (to convert Nitrogen into biomass).
    """
    return nc.variables[keys[0]] * 79.573


def KtoC(nc,keys):
    """
    Loads keys[0] from the netcdf, and converts from Kelvin to Celcius.
    """
    return nc.variables[keys[0]][:] - 273.15


def mul1000(nc,keys):
    """
    Loads keys[0] from the netcdf, but multiplies by 1000.
    """
    return nc.variables[keys[0]][:]* 1000.


def mul1000000(nc,keys):
    """
    Loads keys[0] from the netcdf, but multiplies by 1000000.
    """
    return nc.variables[keys[0]][:]* 1000000.


def div1000(nc,keys):
    """
    Loads keys[0] from the netcdf, then divides by 1000.
    """
    return nc.variables[keys[0]][:]/ 1000.


def div1e6(nc,keys):
    """
    Loads keys[0] from the netcdf, but divides by 1.e6.
    """
    return nc.variables[keys[0]][:]/ 1.e6


def applymask(nc,keys):
    """
    Loads keys[0] from the netcdf, but applies a mask.
    """
    return np.ma.masked_where(nc.variables[keys[1]][:]==0.,nc.variables[keys[0]][:])


def sums(nc,keys):
    """
    Loads Key[0] from the netcdf, then sums the other keys.
    """
    a = nc.variables[keys[0]][:]
    for k in keys[1:]:
        a += nc.variables[k][:]
    return a


def oxconvert(nc,keys):
    """
    Loads keys[0] from the netcdf, but multiplies by 44.771 (to convert oxygen units ).
    """
    return nc.variables[keys[0]][:] *44.661


def convertkgToM3(nc,keys):
    """
    Loads keys[0] from the netcdf, but multiplies by 1.027 (to convert from density kg to volume).
    """
    return nc.variables[keys[0]][:]* 1.027


#####
# kwargs functions:
def multiplyBy(nc,keys, **kwargs):
    """
    Loads keys[0] from the netcdf, but multiplies by the field in kwargs , "factor".
    """
    if 'factor' not in kwargs:
        raise KeyError(f"std_functions:\tmultiplyBy:\t Did not get key word argument, 'factor' in kwargs {str(kwargs)}")
    return nc.variables[keys[0]][:]* float(kwargs['factor'])


def addValue(nc,keys, **kwargs):
    """
    Loads keys[0] from the netcdf, but adds by the field in kwargs key "value".
    """
    if 'value' not in list(kwargs.keys()):
        raise KeyError(f"std_functions:\taddValue:\t Did not get key word argument, 'value' in kwargs {str(kwargs)}")
    return nc.variables[keys[0]][:] + float(kwargs['value'])




#####
std_functions = {}
std_functions[''] = ''
std_functions['NoChange'] = NoChange
std_functions['N2Biomass'] = N2Biomass
std_functions['KtoC'] = KtoC
std_functions['mul1000'] = mul1000
std_functions['mul1000000'] = mul1000000
std_functions['div1000'] = div1000
std_functions['div1e6'] = div1e6
std_functions['applymask'] = applymask
std_functions['sums'] = sums
std_functions['sum'] = sums
std_functions['oxconvert'] = oxconvert
std_functions['convertkgToM3'] = convertkgToM3
std_functions['multiplyBy'] = multiplyBy
std_functions['addValue'] = addValue

#####
# Add lower case, upper, Title, etc...
for key in list(std_functions.keys()):
    func     = std_functions[key]
    std_functions[key.lower()] = func
    std_functions[key.upper()] = func
    std_functions[key.title()] = func
    if len(key)>1:
        std_functions[key[0].upper()+key[1:]] = func
