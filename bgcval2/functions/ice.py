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

import numpy as np
from bgcval2.bgcvaltools.dataset import dataset


tmask     = 0
area    = 0
lat    = 0
loadedArea = False


def loadDataMask(gridfn, maskname,):
    global loadedArea
    global tmask
    global area
    global lat
    if isinstance(gridfn, list) and len(gridfn)==1:
        gridfn = gridfn[0]
    nc = dataset(gridfn, 'r')
    tmask = nc.variables[maskname][0]
    area = nc.variables['e2t'][:] * nc.variables['e1t'][:]
    lat = nc.variables['nav_lat'][:]
    nc.close()
    loadedArea = True


def calcTotalIceArea(nc,keys, **kwargs):    #Global
    if 'areafile' not in kwargs.keys():
        raise AssertionError("calcTotalIceArea:\t Needs an `areafile` kwarg to run calculation.")

    try:
        maskname = kwargs['maskname']
    except:
        maskname = 'tmask'

    if not loadedArea:
        loadDataMask(kwargs['areafile'], maskname)

    arr = nc.variables[keys[0]][:].squeeze() * area
    return np.ma.masked_where(tmask==0, arr).sum() / 1E12


def calcTotalIceAreaN(nc,keys, **kwargs): # North
    if 'areafile' not in kwargs.keys():
        raise AssertionError("calcTotalIceAreaN:\t Needs an `areafile` kwarg to run calculation.")
    try:
        maskname = kwargs['maskname']
    except:
        maskname = 'tmask'

    if not loadedArea:
        loadDataMask(kwargs['areafile'], maskname)

    arr = nc.variables[keys[0]][:].squeeze() * area
    return np.ma.masked_where((tmask==0) + (lat<0.), arr).sum() / 1E12


def calcTotalIceAreaS(nc, keys, **kwargs): # South
    if 'areafile' not in kwargs.keys():
        raise AssertionError("calcTotalIceAreaS:\t Needs an `areafile` kwarg to run calculation.")
    try:
        maskname = kwargs['maskname']
    except:
        maskname = 'tmask'
    if not loadedArea:
        loadDataMask(kwargs['areafile'], maskname)

    arr = nc.variables[keys[0]][:].squeeze() * area
    return np.ma.masked_where((tmask==0)+(lat>0.), arr).sum() / 1E12


def calcTotalIceExtent(nc,keys, **kwargs):    #Global
    if 'areafile' not in kwargs.keys():
        raise AssertionError("calcTotalIceExtent:\t Needs an `areafile` kwarg to run calculation.")
    try:
        maskname = kwargs['maskname']
    except:
        maskname = 'tmask'
    if not loadedArea:
        loadDataMask(kwargs['areafile'], maskname)

    try:
        minIce = float(kwargs['minIce'])
    except:
        minIce = 0.15

    if not loadedArea:
        loadDataMask(kwargs['areafile'],maskname)
    
    return np.ma.masked_where((tmask==0) + (nc.variables[keys[0]][:].squeeze()<minIce), area).sum() / 1E12


def calcTotalIceExtentN(nc, keys, **kwargs): # North
    if 'areafile' not in kwargs.keys():
        raise AssertionError("calcTotalIceExtentN:\t Needs an `areafile` kwarg to run calculation.")
    try:
        maskname = kwargs['maskname']
    except:
        maskname = 'tmask'
    if not loadedArea:
        loadDataMask(kwargs['areafile'], maskname)

    try:
        minIce = float(kwargs['minIce'])
    except:
        minIce = 0.15

    return np.ma.masked_where((tmask==0) + (nc.variables[keys[0]][:].squeeze()<minIce) + (lat<0.), area).sum() / 1E12


def calcTotalIceExtentS(nc, keys, **kwargs): # South

    if 'areafile' not in kwargs.keys():
        raise AssertionError("calcTotalIceExtentS:\t Needs an `areafile` kwarg to run calculation.")
    try:
        maskname = kwargs['maskname']
    except:
        maskname = 'tmask'
    if not loadedArea:
        loadDataMask(kwargs['areafile'], maskname)

    try:
        minIce = float(kwargs['minIce'])
    except:
        minIce = 0.15

    return np.ma.masked_where((tmask==0) + (nc.variables[keys[0]][:].squeeze()<minIce) + (lat>0.), area).sum()/1E12


icedetails = {}
def loadArea(gridfn):
    global icedetails
    nc = dataset(gridfn, 'r')
    lat    = nc.variables['lat'][:]
    try:
        area   = nc.variables['area'][:]
    except:
        area = nc.variables['e2t'][:] * nc.variables['e1t'][:]
    if lat.ndim ==1:
        lon = nc.variables['lon'][:]
        if area.shape[0]==lat.shape[0]:
            lats, lond = np.meshgrid(lon, lat)
        if area.shape[1]==lat.shape[0]:
            lats, lond = np.meshgrid(lat, lon)
    else:
        lats = lat
    nc.close()
    icedetails[(gridfn, 'area')] = area
    icedetails[(gridfn, 'lat')] = lats


def cmipTotalIceExtent(nc, keys, **kwargs):    #Global
    if 'gridfile' not in kwargs.keys():
        raise AssertionError("cmipTotalIceExtent:\t Needs an `gridfile` kwarg to run calculation.")

    gridfn = kwargs['gridfile']
    try:
        area = icedetails[(gridfn, 'area')]
        lat  = icedetails[(gridfn, 'lat')]
    except:
        loadArea(gridfn)
        area = icedetails[(gridfn, 'area')]
        lat  = icedetails[(gridfn, 'lat')]

    try:
        minIce = float(kwargs['minIce'])
    except:
        minIce = 0.15

    try:
        hemisphere = kwargs['hemisphere']
    except:
        hemisphere = 'both'

    sic = np.ma.array(nc.variables[keys[0]][:].squeeze())
    sic = np.ma.masked_where((sic < minIce) + sic.mask, sic)

    northern = ['north', 'northern', 'northernhemisphere', 'northhemisphere']
    southern = ['south', 'southern', 'southernhemisphere', 'southhemisphere']

    if sic.ndim ==3:
        out = []
        for i in np.arange(sic.shape[0]):
            si = sic[i]*area
            if hemisphere.lower() in northern:
                si = np.ma.masked_where(lat<0., si)
            if hemisphere.lower() in southern:
                si = np.ma.masked_where(lat>0., si)

            out.append(si.sum() * 1.E-12)

        return np.ma.array(out)
    elif sic.ndim ==2:
        if hemisphere.lower() in northern:
            sic = np.ma.masked_where(lat<0., sic)
        if hemisphere.lower() in southern:
            sic = np.ma.masked_where(lat>0., sic)
        return (sic*area).sum() * 1.E-12
