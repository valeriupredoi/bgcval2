#!/usr/bin/ipython

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

import numpy as np
from netCDF4 import num2date
import os
from datetime import datetime, timedelta
#from pyproj import Proj

#Specific local code:
from .. import UKESMpython as ukp
from ..netcdf_manipulation import convertToOneDNC
from ..bgcvaltools.dataset import dataset
from ..bgcvaltools.makeMask import makeMask
from ..functions.standard_functions import extractData as std_extractData
"""
.. module:: timeseriesTools
   :platform: Unix
   :synopsis: A swiss army knife set of tools for the time series analysis.
.. moduleauthor:: Lee de Mora <ledm@pml.ac.uk>
"""


def getTimes(nc, coords):
    """
	Loads the times from the netcdf.
	"""
    if type(nc) == type('filename'):
        nc = dataset(nc, 'r')
    try:
        cal = nc.variables[coords['t']].calendar
    except:
        cal = calendar = coords['cal']

    dtimes = num2date(nc.variables[coords['t']][:],
                      nc.variables[coords['t']].units,
                      calendar=cal)
    #nc.variables[coords['t']].calendar)[:]
    #dtimes = num2date(nc.variables[coords['t']][:], nc.variables[coords['t']].units,calendar=coords['cal'])[:]
    print(dtimes)
    try:
        ts = np.array([float(dt.year) + dt.dayofyr / 365. for dt in dtimes])
        return ts
    except:
        pass
    ts = []
    for dt in dtimes:
        t = float(dt.year)
        tdel = dt - datetime(dt.year, 1, 1, 0, 0)
        t += tdel.days / 365.
        ts.append(t)
        print(dt, t)
    return np.array(ts)


def loadData(nc, details):
    """
	Loads the times from the netcdf.
	"""

    if type(nc) == type('filename'):
        nc = dataset(nc, 'r')
    return std_extractData(nc, details)[:]


def ApplyDepthSlice(arr, k):
    """
    Extracts a flat layer, k either
    """
    if arr.ndim < 3:
        return arr
    return arr[..., k, :, :]


def ApplyDepthrange(arr, k1, k2):
    if arr.ndim == 4: return arr[:, k1:k2, :, :]
    if arr.ndim == 3: return arr[k1:k2, :, :]
    if arr.ndim == 2: return arr


def getHorizontalSlice(nc, coords, details, layer, data=''):
    if type(nc) == type('filename'):
        nc = dataset(nc, 'r')
    #####
    # In the case that there is no depth field provided, or the depth field is not the netcdf.
    # We just attempt to extract the data.
    # This is useful
    if coords['z'] == '' or coords['z'] not in list(nc.variables.keys()):
        print("getHorizontalSlice:\tNo depth field in", details['name'])
        if isinstance(data, str):
            data = std_extractData(nc, details)
        return data

    ####
    #
    if len(nc.variables[coords['z']][:]) == 1 and layer in [
            'Surface',
    ]:
        print("getHorizontalSlice:\tNo depth field only 1 value",
              details['name'])
        if data == '':
            data = std_extractData(nc, details)
        return ApplyDepthSlice(data, 0)

    if layer in [
            'layerless',
    ]:
        print("getHorizontalSlice:\tNo layer data requested", layer)
        if isinstance(data, str):
            data = std_extractData(nc, details)
        return data

    #####
    # This is when there is only one dimension in the data/model file.
    if len(list(
            nc.dimensions.keys())) == 1 and layer in ['Surface', 'layerless']:
        print("getHorizontalSlice:\tOne D file", details['name'])
        if data == '':
            data = std_extractData(nc, details)
        data = np.ma.masked_where(nc.variables[coords['z']][:] > 0, data)
        return data
        #return ApplyDepthSlice(data, 0)

    if layer in [
            'Surface',
            '100m',
            '200m',
            '300m',
            '500m',
            '1000m',
            '2000m',
            '3000m',
            '4000m',
    ]:
        if layer == 'Surface': z = 0.
        if layer == '100m': z = 100.
        if layer == '200m': z = 200.
        if layer == '300m': z = 300.
        if layer == '500m': z = 500.
        if layer == '1000m': z = 1000.
        if layer == '2000m': z = 2000.
        if layer == '3000m': z = 3000.
        if layer == '4000m': z = 4000.
        print(z)
        k = ukp.getORCAdepth(z, nc.variables[coords['z']][:], debug=False)
        if isinstance(data, str):
            data = std_extractData(nc, details)
        print("getHorizontalSlice:\tSpecific depth field requested",
              details['name'], layer, [k], nc.variables[coords['z']][k],
              data.shape)
        return ApplyDepthSlice(data, k)

    elif layer in ['Surface - 1000m', 'Surface - 300m']:
        if layer == 'Surface - 300m': z = 300.
        if layer == 'Surface - 1000m': z = 1000.
        k_surf = ukp.getORCAdepth(0.,
                                  nc.variables[coords['z']][:],
                                  debug=False)
        k_low = ukp.getORCAdepth(z, nc.variables[coords['z']][:], debug=False)
        print("getHorizontalSlice:\t", layer, "surface:", k_surf, '-->', k_low)
        if data == '':
            return ApplyDepthSlice(std_extractData(nc, details),
                                   k_surf) - ApplyDepthSlice(
                                       std_extractData(nc, details), k_low)
        return ApplyDepthSlice(data, k_surf) - ApplyDepthSlice(data, k_low)

    elif layer in [
            'Surface to 100m', 'Surface to 300m', 'Surface to 700m',
            'Surface to 2000m'
    ]:
        if layer == 'Surface to 100m': z = 100.
        if layer == 'Surface to 300m': z = 300.
        if layer == 'Surface to 500m': z = 500.
        if layer == 'Surface to 700m': z = 700.
        if layer == 'Surface to 2000m': z = 2000.

        k_surf = ukp.getORCAdepth(0.,
                                  nc.variables[coords['z']][:],
                                  debug=False)
        k_low = ukp.getORCAdepth(z, nc.variables[coords['z']][:], debug=False)
        print("getHorizontalSlice:\t", layer, "surface:", k_surf, '-->', k_low)
        if len(data) == 0:
            return ApplyDepthrange(std_extractData(nc, details), k_surf, k_low)
        return ApplyDepthrange(data, k_surf, k_low)
    elif layer == 'depthint':
        print(
            "getHorizontalSlice\t:ERROR:\tDepth in should be done in the extractData phase by passing a function in the details dictionary."
        )
        assert 0

    if type(layer) in [
            type(0),
            np.int64,
            np.int,
    ]:
        k = layer
        try:
            z = nc.variables[coords['z']][k]
        except:
            return []
        if data == '':
            data = std_extractData(nc, details)
        print("getHorizontalSlice:\tSpecific depth level requested",
              details['name'], layer, nc.variables[coords['z']][k], data.shape)
        return ApplyDepthSlice(data, k)

    if layer in nc.variables[coords['z']][:]:
        z = layer
        k = ukp.getORCAdepth(z, nc.variables[coords['z']][:], debug=False)
        if data == '':
            data = std_extractData(nc, details)
        print("getHorizontalSlice:\tSpecific depth requested", details['name'],
              layer, nc.variables[coords['z']][k], data.shape)
        return ApplyDepthSlice(data, k)

    print("getHorizontalSlice\t:ERROR:\tunrecoginised layer instructions",
          layer, coords, type(layer))
    assert 0


class DataLoader:

    def __init__(self,
                 fn,
                 nc,
                 coords,
                 details,
                 regions=[
                     'Global',
                 ],
                 layers=[
                     'Surface',
                 ],
                 data=''):
        self.fn = fn
        if type(nc) == type('filename'):
            nc = dataset(fn, 'r')
        self.nc = nc
        self.coords = coords
        self.details = details
        self.regions = regions
        self.layers = layers
        self.name = self.details['name']
        if data == '': data = std_extractData(nc, self.details)
        self.Fulldata = data
        self.__lay__ = -999.
        self._makeTimeDict_()
        self.run()

    def run(self):
        self.load = {}
        try:
            depths = {
                i: z
                for i, z in enumerate(self.nc.variables[self.coords['z']][:])
            }
        except:
            depths = {}
        #print "self.nc.variables[self.coords[z]][:]", self.nc.variables[self.coords['z']][:]
        #print "self.coords[z]", self.coords['z']
        #print "depths",depths
        #print "layers",self.layers

        #  	assert 0
        lays = self.layers[:]
        lays.reverse()

        for l in lays:  #self.layers:
            try:
                layer = int(l)
            except:
                layer = l

            print(l, layer, type(layer),
                  type(layer) in [
                      type(1),
                      type(1.),
                  ], layer not in list(depths.keys()))
            if type(layer) in [
                    type(1),
                    type(1.),
            ] and layer not in list(depths.keys()):
                print("DataLoader:\tLayer outside depths. Setting:", layer,
                      'to a masked value.')
                for region in self.regions:
                    self.maskedload(region, layer)
                continue

            for region in self.regions:
                arr, arr_t, arr_z, arr_lat, arr_lon = self.createDataArray(
                    region, layer)
                if len(arr) == 0:
                    self.maskedload(region, layer)
                    continue

                if arr.mask.all():
                    self.maskedload(region, layer)
                    continue

                self.load[(region, layer)] = arr
                self.load[(region, layer, 't')] = arr_t
                self.load[(region, layer, 'z')] = arr_z
                self.load[(region, layer, 'lat')] = arr_lat
                self.load[(region, layer, 'lon')] = arr_lon

                print("DataLoader:\tLoaded", self.name, 'in', end=' ')
                print('{:<24} layer: {:<8}'.format(region, layer), end=' ')
                print('\tdata length:',
                      len(self.load[(region, layer)]),
                      end=' ')
                print('\tmean:',
                      self.load[(region, layer)].mean(),
                      'of',
                      len(self.load[(region, layer)]),
                      end=' ')
                print('\trange:', [
                    self.load[(region, layer)].min(),
                    self.load[(region, layer)].max()
                ])

    def _makeTimeDict_(self, ):
        """ Make a dictionairy linking the time index with the float time.
	"""
        try:
            ts = getTimes(self.nc, self.coords)
        except:
            print(
                "DataLoaded:\t_makeTimeDict_:\tUnable to load time array, probably due to time zero in file."
            )
            return
        self.timedict = {i: t for i, t in enumerate(ts)}
        self.timedict_ti = {t: i for i, t in enumerate(ts)}

    def maskedload(self, region, layer):
        """ Quick in line tool to set a layer/region to masked.
  	"""
        maskedValue = np.ma.array([
            -999.,
        ], mask=[
            True,
        ])
        self.load[(region, layer)] = maskedValue
        self.load[(region, layer, 't')] = maskedValue
        self.load[(region, layer, 'z')] = maskedValue
        self.load[(region, layer, 'lat')] = maskedValue
        self.load[(region, layer, 'lon')] = maskedValue
        print("DataLoader:\tLoaded empty", self.name, 'in', end=' ')
        print('{:<24} layer: {:<8}'.format(region, layer), end=' ')
        print('\tdata length:', len(self.load[(region, layer)]))  #,
#print #, '\tmean:',np.ma.mean(self.load[(region,layer)])

    def __getlayerDat__(self, layer):
        """ Minimise quick load and save to minimise disk-reading time.
  	"""

        if self.__lay__ == layer:
            return self.__layDat__
        else:
            self.__layDat__ = np.ma.array(
                getHorizontalSlice(self.nc,
                                   self.coords,
                                   self.details,
                                   layer,
                                   data=self.Fulldata))
            print("DataLoader:\tgetlayerDat:", self.name, layer)
            self.__lay__ = layer
            return self.__layDat__

    def createDataArray(self, region, layer):
        """	
  		This creates a set of 1D arrays of the dat and 4D coordinates for the required region.
  		The leg work is done in makeMask.py
  	"""

        #print 'DataLoader:\tcreateDataArray:\t',self.details['name'],region,layer

        self.createOneDDataArray(layer)

        m = makeMask(
            self.details['name'],
            region,
            self.oneDData['arr_t'],
            self.oneDData['arr_z'],
            self.oneDData['arr_lat'],
            self.oneDData['arr_lon'],
            self.oneDData['arr'],
        )

        return  np.ma.masked_where(m,self.oneDData['arr']),\
         np.ma.masked_where(m,self.oneDData['arr_t']),\
         np.ma.masked_where(m,self.oneDData['arr_z']),\
         np.ma.masked_where(m,self.oneDData['arr_lat']),\
         np.ma.masked_where(m,self.oneDData['arr_lon'])

    def createOneDDataArray(self, layer):
        """ 	This is a relatively simple routine that takes a layer and makes a series of 1D arrays containing points.
  		These output 1D arrays are then passed to UKESMpython.py's makemasks toolkit.
  	"""

        #####
        # load lat, lon and data.
        if self.coords['lat'] == self.coords['lon'] == False:
            ####
            # scalar fields
            dat = self.__getlayerDat__(layer)

            lat = np.zeros_like(dat)
            lon = np.zeros_like(dat)
            dims = self.nc.variables[self.details['vars'][0]].dimensions

        else:
            if self.coords['lat'] not in self.nc.variables or self.coords['lon'] not in self.nc.variables:
                raise KeyError(f"ERROR: coordinates provided do not match coordinates in file: {self.coords['lat']}, {self.coords['lon']}")
            lat = self.nc.variables[self.coords['lat']][:]
            lon = ukp.makeLonSafeArr(self.nc.variables[self.coords['lon']]
                                     [:])  # makes sure it's between +/-180

            dims = self.nc.variables[self.details['vars'][0]].dimensions
            dat = self.__getlayerDat__(layer)
        if dat.ndim == 2: dat = dat[None, :, :]

        try:
            l = len(dat)
        except:
            dat = np.ma.array([
                dat,
            ])
            l = len(dat)
            print("createOneDDataArray: \tWarning:\tdata was a single float:",
                  self.name, dat, l, dat.shape, dat.ndim)
        if l == 0:
            a = np.ma.array([
                -999.,
            ], mask=[
                True,
            ])
            return a, a, a, a, a


        #####
        # Create Temporary Output Arrays.
        arr = []
        arr_lat = []
        arr_lon = []
        arr_t = []
        arr_z = []

        latnames = [
            'lat',
            'latitude',
            'latbnd',
            'nav_lat',
            'y',
            'lat',
        ]
        lonnames = [
            'lon',
            'longitude',
            'lonbnd',
            'nav_lon',
            'x',
            'lon',
        ]

        ####
        # Different data has differnt shapes, and order or dimensions, this takes those differences into account.
        if dat.ndim > 2:
            if dims[-2].lower() in latnames and dims[-1].lower() in lonnames:

                #print 'createDataArray',self.details['name'],layer, "Sensible dimsions order:",dims
                for index, v in ukp.maenumerate(dat):
                    try:
                        (t, z, y, x) = index
                    except:
                        (t, y, x) = index
                        z = 0

                    try:
                        la = lat[y, x]
                    except:
                        la = lat[y]
                    try:
                        lo = lon[y, x]
                    except:
                        lo = lon[x]
                    arr.append(v)
                    arr_t.append(t)
                    arr_z.append(z)
                    arr_lat.append(la)
                    arr_lon.append(lo)

            elif dims[-2].lower() in lonnames and dims[-1].lower() in latnames:

                for index, v in ukp.maenumerate(dat):
                    try:
                        (t, z, x, y) = index
                    except:
                        (t, x, y) = index
                        z = 0

                    try:
                        la = lat[y, x]
                    except:
                        la = lat[y]
                    try:
                        lo = lon[y, x]
                    except:
                        lo = lon[x]

                    arr.append(v)
                    arr_t.append(t)
                    arr_z.append(z)
                    arr_lat.append(la)
                    arr_lon.append(lo)

            else:
                print("Unknown dimensions order", dims)
                assert False

        elif dat.ndim == 1:
            if dims[0] == 'index':
                print('createDataArray', self.name, layer, "1 D data:", dims,
                      dat.shape)
                for i, v in enumerate(dat):
                    la = lat[i]
                    lo = lon[i]

                    arr.append(v)
                    arr_t.append(0)
                    arr_z.append(0)
                    arr_lat.append(la)
                    arr_lon.append(lo)
            elif len(dat) == 1:
                print('createDataArray', self.name, layer,
                      "single point data:", dims, dat.shape)
                arr = dat
                arr_t = [
                    0,
                ]
                arr_z = [
                    0.,
                ]
                arr_lat = [
                    0.,
                ]
                arr_lon = [
                    0.,
                ]

            else:
                print("Unknown dimensions order", dims)
                assert False
        else:
            print("Unknown dimensions order", dims)
            assert False

        arr = np.ma.masked_invalid(np.ma.array(arr))
        mask = np.ma.masked_where((arr > 1E20) + arr.mask, arr).mask

        self.oneDData = {}
        self.oneDData['arr_lat'] = np.ma.masked_where(mask,
                                                      arr_lat).compressed()
        self.oneDData['arr_lon'] = np.ma.masked_where(mask,
                                                      arr_lon).compressed()
        self.oneDData['arr_z'] = np.ma.masked_where(mask, arr_z).compressed()
        self.oneDData['arr_t'] = np.ma.masked_where(mask, arr_t).compressed()
        self.oneDData['arr'] = np.ma.masked_where(mask, arr).compressed()


def makeArea(fn, coordsdict):
    nc = dataset(fn, 'r')
    lats = nc.variables[coordsdict['lat']][:]
    lons = nc.variables[coordsdict['lon']][:]
    #depths = nc.variables[coordsdict['z']][:]
    nc.close()
    if lats.ndim == 1:
        #lat2d,lon2d = np.meshgrid(lats,lons)
        area = np.zeros((len(lats), len(lons)))
        meanLatDiff = np.abs(lats[:-1] - lats[1:]).mean()
        meanLonDiff = np.abs(lons[:-1] - lons[1:]).mean()
        for a in np.arange(len(lats)):
            #area[a] = np.ones(len(lats)*calculateArea(lats[a]-meanLatDiff,lats[a]+meanLatDiff,-meanLonDiff,meanLonDiff))
            print(a, area.shape, len(lats), len(lons))
            area[a, :] = np.ones(len(lons)) * ukp.Area(
                [lats[a] - meanLatDiff / 2., -meanLonDiff / 2.],
                [lats[a] + meanLatDiff / 2., meanLonDiff / 2.])
        return area
    elif lats.ndim == 2:
        print(
            "timeseriesTools.py:\tWARNING: Setting area to flat for uneven grid! "
        )
        return np.ones_like(lats)

    else:
        print("timeseriesTools.py:\tNot implemeted makeArea for this grid. ",
              lats.ndim, coordsdict)
        assert 0, 'timeseriesTools.py:\tNot implemeted makeArea for this grid. ' + str(
            lats.ndim)


#def calculateArea(lat0,lat1,lon0,lon1):
#		co = {"type": "Polygon", "coordinates": [
#		    [(lon0, lat0), #('lon', 'lat')
#		     (lon0, lat1),
#		     (lon1, lat1),
#		     (lon1, lat0),
#		     (lon0, lat0)]]}
#		clon, clat = zip(*co['coordinates'][0])
#
#		pa = Proj("+proj=aea +lat_1="+str(lat0)+" +lat_2="+str(lat1)+"+lon_0="+str(lon0)+" +lon_1="+str(lon0))
#		x, y = pa(clon, clat)
#		cop = {"type": "Polygon", "coordinates": [zip(x, y)]}
#		area = shape(cop).area
#		return area
#
#
#def calculateVol(lat0,lat1,lon0,lon1,d0,d1):
#		a = calculateArea(lat0,lat1,lon0,lon1)
#		return area*abs(d1-d0)


def calcCuSum(times, arr):
    newt, cusum = [], []

    c = 0.
    for i, ti in enumerate(times):
        if i == 0: continue
        t0 = times[i - 1]
        t = t0 + (ti - t0) / 2.
        c += arr[i] - arr[i - 1]
        newt.append(t)
        cusum.append(c)

    return newt, cusum
