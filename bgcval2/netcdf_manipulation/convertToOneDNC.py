#
# Copyright 2015, Plymouth Marine Laboratory
#
# Address:
# Plymouth Marine Laboratory
# Prospect Place, The Hoe
# Plymouth, PL1 3DH, UK
#
# Email:
# ledm@pml.ac.uk
#
# This file is part of the netcdf_manip library.
#
# netcdf_manip is free software: you can redistribute it and/or modify it
# under the terms of the Revised Berkeley Software Distribution (BSD) 3-clause license.
#
# netcdf_manip is distributed in the hope that it will be useful, but
# without any warranty; without even the implied warranty of merchantability
# or fitness for a particular purpose. See the revised BSD license for more details.
# You should have received a copy of the revised BSD license along with netcdf_manip.
# If not, see <http://opensource.org/licenses/BSD-3-Clause>.
#

try:
    from netCDF4 import Dataset, default_fillvals
except:
    from netCDF4 import Dataset, _default_fillvals
from datetime import date
from getpass import getuser
from os.path import exists
import numpy as np
from numpy import array, int64
from numpy.ma import array as marray, nonzero, masked_where, compressed, zeros
from .pruneNC import todaystr
from operator import itemgetter
from .alwaysInclude import alwaysInclude, depthNames, timeNames, latnames, lonnames
"""	This routine takes a netcdf input and creates a new one with only 1 dimension. """

# list of variables to save, assuming some conventions
#alwaysInclude = ['time', 'lat','lon', 'latbnd', 'lonbnd','LONGITUDE','LATITUDE','DEPTH','TIME','nav_lat','nav_lon','nav_lev','time_counter','deptht','depth', 'latitude', 'longitude','month','mask']
#'crs',]'lat_bnds','lon_bnds',


def getCoordsToKeep_new(nc, variables, newMask='', debug=False):
    """This routine takes an ncdfView object, a list of variables,
	   and tests which coordinates should be saved.
	   A Mask can be applied instead.
	"""
    CoordsToKeep = {}
    for var in variables:
        if var in alwaysInclude: continue
        arr = nc.variables[var][:]
        if len(newMask):
            out = []
            if newMask.shape != arr.shape:
                if debug: 'getCoordsToKeep:\t', var, 'Wrong shape'
                continue
            try:  # 1D arrays
                for i, m in enumerate(newMask):
                    if not m: continue
                    out.append(arr[i])
                arr = array(out).squeeze()
            except:
                # multi dimensional arrays
                arr = masked_where(newMask, array(arr))
        #nearlyZero = 1.E-6
        nzdims = len(arr)
        for i, a in enumerate(arr[0]):
            coords = tuple([arr[j][i] for j in range(nzdims)])
            try:
                if i in CoordsToKeep[coords]: pass
            except:
                try:
                    CoordsToKeep[coords].append(i)
                except:
                    CoordsToKeep[coords] = [
                        i,
                    ]
        if debug:
            print("getCoordsToKeep:\t", var, "\tndims:", nzdims, len(nz[0]),
                  "\tNumber of Coords:", len(list(CoordsToKeep.keys())))
    return CoordsToKeep, variables


def getCoordsToKeep(nc, variables, newMask='', debug=False):
    """This routine takes an ncdfView object, a list of variables,
	   and tests which coordinates should be saved.
	   A Mask can be applied instead.
	"""
    CoordsToKeep = {}
    if debug: print('getCoordsToKeep:', variables)
    if not len(variables):
        raise ValueError('getCoordsToKeep: not variables provided')
    for var in variables:
        if var in alwaysInclude: continue
        arr = nc.variables[var][:]
        print('checking out var:', var, arr.shape, arr.mean())
        if len(newMask):
            out = []
            if newMask.shape != arr.shape:
                raise ValueError('getCoordsToKeep:\t', var, 'Wrong shape')
            arr = masked_where(newMask, arr)
            print('getCoordsToKeep: now masked:', arr.shape, arr.mean())

        nz = nonzero(arr)  
        if not len(nz):
            print('Var failed:', var)
            variables.remove(var)
            continue
        nzdims = len(nz)
        if not nzdims:
            raise ValueError('No non-zero values found!', nzdims, nz)
        print('nz:', nz)
        for i, a in enumerate(nz[0]):
            if var in ['OBSERVATIONS', ]:
                coords = tuple([nz[j][i] for j in range(nzdims)])
                print("NEW OBSERVATION LOCATION:",i,coords)
                try:
                    if i in CoordsToKeep[coords]: pass
                except:
                    try:
                        CoordsToKeep[coords].append(i)
                    except:
                        CoordsToKeep[coords] = [
                            i,
                        ]
            else:
                print('getCoordsToKeep',i, a, nzdims)
                coords = tuple([nz[j][i] for j in range(nzdims)])
                print(coords)
                try:
                    if i in CoordsToKeep[coords]: pass
                except:
                    try:
                        CoordsToKeep[coords].append(i)
                    except:
                        CoordsToKeep[coords] = [
                            i,
                        ]
        if debug:
            print("getCoordsToKeep:\t", var, "\tndims:", nzdims, len(nz[0]),
                  "\tNumber of Coords:", len(list(CoordsToKeep.keys())))
    return CoordsToKeep, variables


class convertToOneDNC:

    def __init__(self,
                 filenameIn,
                 filenameOut,
                 newMask='',
                 variables=[],
                 debug=False,
                 dictToKeep=''):
        self.fni = filenameIn
        self.fno = filenameOut
        self.vars = variables
        self.newMask = newMask
        self.debug = debug
        self.dictToKeep = dictToKeep
        self.run()

    def run(self):
        nci = Dataset(self.fni, 'r')
        if not self.vars:
            assert 0 
        if self.vars == 'all':
            self.vars = list(nci.variables.keys())  # save all

        #create dataset and header.
        if self.debug:
            print('convertToOneDNC:\tINFO:\tCreating a new dataset:\t',
                  self.fno)
        nco = Dataset(self.fno, 'w')
        for a in nci.ncattrs():
            if self.debug:
                print(
                    'convertToOneDNC:\tINFO:\tcopying attribute: \t\"' + a +
                    '\":\t', nci.getncattr(a))
            nco.setncattr(a, nci.getncattr(a))

        appendToDesc = 'Reprocessed on ' + todaystr() + ' by ' + getuser(
        ) + ' using convertToOneDNC.py'
        try:
            nco.Notes = nci.Notes + '\n\t\t' + appendToDesc
        except:
            nco.Notes = appendToDesc

        save = list(set(nci.variables.keys()).intersection(set(alwaysInclude)))
        save = sorted(list(set(sorted(save + self.vars))))
        save = find_vars_in_nc(nci, save)

        # test to find out which coordinates should be saved.
        if not self.dictToKeep:
            CoordsToKeep, save = getCoordsToKeep(nci,
                                                 save,
                                                 newMask=self.newMask,
                                                 debug=self.debug)
        else:
            CoordsToKeep = self.dictToKeep

        print("convertToOneDNC:\tinfo:\tCoordsToKeep:",
              save)  
        # create dimensions:
        nco.createDimension('index', None)

        # create Variables:
        nco.createVariable('index', int64, [
            'index',
        ], zlib=True, complevel=5)  
        nco.createVariable('index_t',
                           int64, [
                               'index',
                           ],
                           zlib=True,
                           complevel=5)
        nco.createVariable('index_z',
                           int64, [
                               'index',
                           ],
                           zlib=True,
                           complevel=5)
        nco.createVariable('index_y',
                           int64, [
                               'index',
                           ],
                           zlib=True,
                           complevel=5) 
        nco.createVariable('index_x',
                           int64, [
                               'index',
                           ],
                           zlib=True,
                           complevel=5)
        for var in save:
            if var in ['index', 'index_t', 'index_z', 'index_y', 'index_x']:
                continue
            nco.createVariable(var,
                               nci.variables[var].dtype, [
                                   'index',
                               ],
                               zlib=True,
                               complevel=5)

        # Long Names:
        nco.variables['index'].long_name = 'index'
        nco.variables['index_t'].long_name = 'index - time'
        nco.variables['index_z'].long_name = 'index - depth'
        nco.variables['index_y'].long_name = 'index - latitude'
        nco.variables['index_x'].long_name = 'index - longitude'

        for var in save:
            try:
                long_name = nci.variables[var].long_name
            except:
                if self.debug:
                    print('convertToOneDNC:\tWARNING:\tNo long_name for ', var)
                long_name = var

            nco.variables[var].long_name = long_name
            if self.debug:
                print('convertToOneDNC:\tINFO:\tAdding long_name for ', var,
                      '\t', long_name)

        # Units:
        nco.variables['index'].units = ''
        nco.variables['index_t'].units = ''
        nco.variables['index_z'].units = ''
        nco.variables['index_y'].units = ''
        nco.variables['index_x'].units = ''

        for var in save:
            try:
                nco.variables[var].units = nci.variables[var].units
                if self.debug:
                    print('convertToOneDNC:\tINFO:\tAdding units for ', var,
                          '\t', nci.variables[var].units)
            except:
                print('convertToOneDNC:\tWARNING:\tNo units for ', var)

        # Fill Values:
        def itemsgetter(a):
            return a[1][0]

        if not len(CoordsToKeep.keys()):
            raise ValueError('did not find enough coordinates worth keeping:', CoordsToKeep)

        sorted_Coords = sorted(iter(CoordsToKeep.items()), key=itemsgetter)

        print("convertToOneDNC:\tINFO:\tsorted_Coords:", sorted_Coords)
        if not len(sorted_Coords):
            raise ValueError('Sorting Coordinates failed')
        data = {}
        if self.debug:
            print('convertToOneDNC:\tINFO:\tCopying index  ...',
                  len(sorted_Coords))

        nco.variables['index'][:] = array([a[1][0] for a in sorted_Coords])
        nco.sync()

        # 4D coordinates:
        if len(sorted_Coords[0][0]) == 4:
            print("4D:", sorted_Coords[0][0])
            if self.debug:
                print('convertToOneDNC:\tINFO:\tCopying index t ...')
            nco.variables['index_t'][:] = array(
                [a[0][0] for a in sorted_Coords])
            nco.sync()
            if self.debug:
                print('convertToOneDNC:\tINFO:\tCopying index z ...')
            nco.variables['index_z'][:] = array(
                [a[0][1] for a in sorted_Coords])
            nco.sync()
            if self.debug:
                print('convertToOneDNC:\tINFO:\tCopying index y ...')
            nco.variables['index_y'][:] = array(
                [a[0][2] for a in sorted_Coords])
            nco.sync()
            if self.debug:
                print('convertToOneDNC:\tINFO:\tCopying index x ...')
            nco.variables['index_x'][:] = array(
                [a[0][3] for a in sorted_Coords])
            nco.sync()

        if len(sorted_Coords[0][0]) == 3:
            print("3D:", sorted_Coords[0][0])
            if self.debug:
                print('convertToOneDNC:\tINFO:\tCopying index t ...')
            nco.variables['index_t'][:] = array(
                [a[0][0] for a in sorted_Coords])
            nco.sync()

            nco.variables['index_z'][:] = zeros(len(sorted_Coords))
            if self.debug:
                print('convertToOneDNC:\tINFO:\tCopying index y ...')
            nco.variables['index_y'][:] = array(
                [a[0][1] for a in sorted_Coords])
            nco.sync()
            if self.debug:
                print('convertToOneDNC:\tINFO:\tCopying index x ...')
            for a in sorted_Coords:
            	print('sorted_coords, a', a)

            nco.sync()
        errorScript = "If it is failing here, then you need to check that your dimensions are names correctly in netcdf_manip/alwaysInclude.py:timeNames, depthNames, lonnames and latnames"

        for var in save:
            if self.debug:
                print('convertToOneDNC:\tINFO:\tCopying ', var, ' ...')
            arr = nci.variables[var][:]


            if len(arr) == 0:
                print('convertToOneDNC:\tWarning:\tIt looks like the netcdf ',
                      self.fni, 'does not contain the variable', var)

            outarr = []
            if arr.ndim == 1 and len(sorted_Coords[0][0]) == 4:
                if var.lower() in [
                        'time', 'time_counter', 't', 'month', 'time_centered', 'index_t',
                ]:
                    d = 0

                elif var.lower() in timeNames: d = 0
                elif var.lower() in depthNames: d = 1
                elif var.lower() in latnames: d = 2
                elif var.lower() in lonnames: d = 3
                else:
                    raise ValueError('Dimension not recognised:  %s' % var)
                    assert 0

                print(var,
                      'convertToOneDNC:\tINFO:\tndim: (1-4)',
                      arr.ndim, arr.shape,
                      var,
                      sorted_Coords[0][0],
                      d
                      )
                for c in sorted_Coords:
                    outarr.append(arr[c[0][d]])

                try:
                    print(var, d)
                except:
                    var, "not found"

            elif arr.ndim == 1 and len(sorted_Coords[0][0]) == 1:
                d = 0
                print(var, 'convertToOneDNC:\tINFO:\tndim: (1-1)', arr.ndim,
                      var, sorted_Coords[0][0], d)
                for c in sorted_Coords:
                    try:
                        outarr.append(arr[c[0][d]])
                    except:
                        raise AssertionError(errorScript)
                try:
                    print(var, d)
                except:
                    var, "not found"

            elif arr.ndim == 1 and len(sorted_Coords[0][0]) == 3:
                if var.lower() in [
                        'time',
                        'time_counter',
                        't',
                        'month',
                        'time_centered',
                        'index_t',
                ]:
                    d = 0
                if var.lower() in timeNames: d = 0
                if var.lower() in latnames: d = 1
                if var.lower() in lonnames: d = 2

                try:
                    print(var, d)
                except:
                    print(var, "not found")

                print(var, 'convertToOneDNC:\tINFO:\tndim: (1-3)', arr.ndim,
                      var, sorted_Coords[0][0], d)
                for c in sorted_Coords:
                        print('sorted_coords:', d, c, var, c[0], len(arr), arr.shape, arr.ndim)
                        outarr.append(arr[c[0][d]])

            elif arr.ndim == 2:
                if var.lower() in ['nav_lat', 'nav_lon']:
                    d = (2, 3)
                elif var.lower() in [
                        'deptht',
                        'depthu',
                        'depthv',
                        'latitude',
                        'longitude',
                        'depth',
                        'index_z',
                ]:
                    d = (0, 1)
                elif var.lower() in ['mask'] and len(sorted_Coords[0][0]) == 3:
                    d = (1, 2)  # because of MLD datasets.
                else:
                    d = (0, 1)

                print(var, 'convertToOneDNC:\tINFO:\tndim: (2)', arr.ndim, var,
                      sorted_Coords[0][0], d)
                for c in sorted_Coords:  
                    try:
                        outarr.append(arr[(c[0][d[0]], c[0][d[1]])])
                    except:
                        raise AssertionError(errorScript)
            elif arr.ndim == 3:
                print(var, 'convertToOneDNC:\tINFO:\tndim: (3)', arr.ndim, var,
                      sorted_Coords[0][0])
                for c in sorted_Coords:
                    if len(c[0]) == 4:
                        outarr.append(arr[(c[0][0], c[0][2], c[0][3])])
                    if len(c[0]) == 3:
                        outarr.append(arr[(c[0][0], c[0][1], c[0][2])])
            elif arr.ndim == 4:
                print(var, 'convertToOneDNC:\tINFO:\tndim: (4)', arr.ndim, var,
                      sorted_Coords[0][0])
                for c in sorted_Coords:
                    outarr.append(arr[(c[0][0], c[0][1], c[0][2], c[0][3])])
            else:
                print("How many dimensions?", arr.ndim,
                      len(sorted_Coords[0][0]))
                assert False
            outarrnm = np.array(outarr)
            print(var, [outarrnm.min(), outarrnm.max()], '(no mask)')
            outarr = marray(outarr)
            print(var, [outarr.min(), outarr.max()], '(masked)')
            if self.debug:
                print('convertToOneDNC:\tINFO:\tSaving var:', var, arr.shape,
                      '->', outarr.shape, 'coords:', len(sorted_Coords),
                      [outarr.min(), outarr.max()])
            nco.variables[var][:] = outarr
            nco.sync()
        # Close netcdfs:
        nco.close()
        nci.close()
        if self.debug:
            print('convertToOneDNC:\tINFO:\tsuccessfully created:\t', self.fno)


def find_vars_in_nc(nc, keys):
    """
    Takes the list of keys and chooses the first one that exists in the input file.
    Useful if fields change for no reason.
    """
    outlist=[]
    for key in keys:
        if key not in nc.variables.keys():
            continue
        outlist.append(key)
    if outlist: 
        print('find_vars_in_nc: found these keys:',outlist)
        return outlist
    raise KeyError(f'find_vars_in_nc: unable to find any variable in {keys} in {nc.filename}')
