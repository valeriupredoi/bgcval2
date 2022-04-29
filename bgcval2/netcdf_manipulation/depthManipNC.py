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

# this file takes a netcdf in, a netcdf out and a list of variables to save.
# it automatically saves the dimensions, and the header.
# it also reduces the depth dimension to length 1. either via depth integration, or removing all but the surface layer or the deepst layer for each variable independently.

#	the input file containts a 5 dimensional object zbnd has
#	the 2nd dimension of that object is depth
#	the input variables

# This class takes an input netcdf filename, an output netcdf filename and a list of variables to keep. Usually ['var', 'lat','lon',time','depth'].
# and a list of -1,0,or 1 that is the same length as vars, where
#options: #	-1 : deepest layer only
#		 0 : surface,
#		 1 : depth integrated (sum). Remember to keep the mixed layer depth too.
#		-15: Bottom 15 m meter
# 		other: a slice in depth: ie s = slice(0,5) for top 4 layers.

# It creates a new netcdf that contains only the variables that you gave it.

#timemean = True will reduce all elements with a time dimension to size 1 and take the mean.

#from ncdfView import ncdfView
try:
    from netCDF4 import Dataset, default_fillvals
except:
    from netCDF4 import Dataset, _default_fillvals
from datetime import date
from getpass import getuser
from os.path import exists
from numpy import zeros, array, repeat, where
from numpy.ma import array as marray, masked_where
from collections import defaultdict


class depthManipNC:

    def __init__(self,
                 filenameIn,
                 filenameOut,
                 variables,
                 depthFlags='',
                 timemean=False,
                 debug=False,
                 alwaysInclude=['time', 'lat', 'lon', 'latbnd', 'lonbnd']):
        self.fni = filenameIn
        self.fno = filenameOut
        self.vars = variables
        self.depthFlags = depthFlags
        self.timemean = timemean
        self.debug = debug
        self.alwaysInclude = alwaysInclude

        def returnEmptyStr():
            return ''

        self.depthStrings = defaultdict(returnEmptyStr, (
            ['-2', '(Deepest layer)'],
            ['-1', '(Deepest layer)'],
            ['0', '(Surface)'],
            ['1', '(Depth Integrated)'],
            ['-15', '(Deepest layer)'],
        ))
        self.run()

    def run(self):
        if not self.vars:
            print('depthManipNC:\tERROR:\tvariables to save are no use:',
                  self.vars)
            return
        if not exists(self.fni):
            print('depthManipNC:\tERROR:\tinputfile name does not exists:',
                  self.fni)
            return

        if self.depthFlags == '':
            print(
                'depthManipNC:\tWARNING:\tNo depth flags given, assuming surface values only.'
            )
            self.depthFlags = zeros(len(self.vars), dtype=int)

        if len(self.vars) != len(self.depthFlags):
            print('depthManipNC:\tERROR:\tVariables do not match depth flags:',
                  len(self.vars), '!=', len(self.depthFlags))
            return
        self.varflag = {}
        for var, flag in zip(self.vars, self.depthFlags):
            self.varflag[var] = flag

        if self.debug:
            print('depthManipNC:\tINFO:\topening dataset:\t', self.fni)
        nci = Dataset(self.fni, 'r')  #Quiet =True)
        #if self.depthFlags and 'zbnd' not in nci.variables.keys():
        #	print 'depthManipNC:\tERROR:\tdepthFlags is ',self.depthFlags,'but inputfile name does contain \'zbnd\''
        #	return

        #check that there are some overlap between input vars and nci:
        for v in self.vars:
            if v in list(nci.variables.keys()): continue
            print('depthManipNC:\tERROR:\tvariable,', v, ', not found in ',
                  self.fni)
            return

        #create dataset and header.
        if self.debug:
            print('depthManipNC:\tINFO:\tCreating a new dataset:\t', self.fno)
        nco = Dataset(self.fno, 'w')
        for a in nci.ncattrs():
            if self.debug:
                print(
                    'depthManipNC:\tINFO:\tcopying attribute: \t\"' + a +
                    '\":\t', nci.getncattr(a))
            nco.setncattr(a, nci.getncattr(a))
        appendToDesc = 'Reprocessed on ' + todaystr() + ' by ' + getuser(
        ) + ' using depthManipNC.py'
        try:
            nco.Notes = nci.Notes + '\n\t\t' + appendToDesc
        except:
            nco.Notes = appendToDesc

        # list of variables to save, assuming some conventions
        save = list(
            set(nci.variables.keys()).intersection(set(self.alwaysInclude)))
        save = list(set(sorted(save + self.vars)))

        # create dimensions:
        for d in list(nci.dimensions.keys()):
            if d in [
                    'time',
            ]:
                nco.createDimension(d, None)
            elif d in [
                    'depth',
                    'z',
            ]:
                nco.createDimension(d, 1)
            else:
                nco.createDimension(d, len(nci.dimensions[d]))

        # create Variables:
        for var in save:
            nco.createVariable(var,
                               nci.variables[var].dtype,
                               nci.variables[var].dimensions,
                               zlib=True,
                               complevel=5)

        # Long Names:
        for var in save:
            varln = ''
            long_name = ''
            try:
                long_name = nci.variables[var].long_name
            except:
                if self.debug:
                    print('depthManipNC:\tWarning:\tNo long_name for ', var)
            if var in self.vars:
                long_name += ' ' + self.depthStrings[str(self.varflag[var])]
            if self.timemean: long_name.replace('Daily', '')
            if self.timemean: long_name.replace('Monthly', '')
            nco.variables[var].long_name = long_name
            if self.debug:
                print('depthManipNC:\tInfo:\tAdding long_name:', var,
                      long_name)
        # Units:
        for var in save:
            units = ''
            try:
                units = nci.variables[var].units
            except:
                if self.debug:
                    print('depthManipNC:\tWarning:\tNo units for ', var)
            if var in self.vars:
                if self.varflag[var] == 1: units = units.replace('m^3', 'm^2')

            nco.variables[var].units = units
            if self.debug:
                print('depthManipNC:\tInfo:\tAdding units:', var, units)

        if 'zbnd' in list(nci.variables.keys()):
            self.zbnd = nci.variables['zbnd'][:]
        if 'bathymetry' in list(nci.variables.keys()):
            self.bathy = nci.variables['bathymetry'][:]

        # Fill Values:
        for var in save:
            if var not in self.vars:  #no change
                arr = nci.variables[var][:]
            else:
                flag = self.varflag[var]
                if self.debug:
                    print('depthManipNC:\tInfo:\tFilling var:', var, 'flag:',
                          flag)
                if flag == 1:
                    arr = (nci.variables[var][:] * abs(
                        (self.zbnd[:, :, :, :, 1] - self.zbnd[:, :, :, :, 0]))
                           ).sum(1)
                    arr = arr[:, None, :, :]
                elif flag in [-2, -1, 0]:
                    arr = nci.variables[var][:, flag, :, :]
                elif flag in [
                        -15,
                ]:
                    arr = self.bottomLayer(nci, var)
                else:
                    arr = nci.variables[var][:, flag, :, :].mean(1)
                    arr = arr[:, None, :, :]

                #while len(arr.shape) < len(nci.variables[var].dimensions): arr = marray(arr[None,:])

            timeav = False
            if self.timemean and len(
                    intersection(['time', 't'],
                                 nci.variables[var].dimensions)):
                timeav = True
            if timeav:
                if self.debug:
                    print('depthManipNC:\tInfo:\tSaving time averaged var:',
                          var)
                arr = marray([
                    arr.mean(0),
                ])
                while len(arr.shape) < len(nci.variables[var].dimensions):
                    arr = marray(arr[None, :])

            if self.debug:
                print('depthManipNC:\tInfo:\tSaving var:', var, arr.shape,
                      '\tdims:', nci.variables[var].dimensions)
            nco.variables[var][:] = arr

        # Close netcdfs:
        nco.close()
        nci.close()
        if self.debug:
            print('depthManipNC:\tINFO:\tsuccessfully created:\t', self.fno)
        return

    def bottomLayer(self, nci, var):
        if self.debug: print('depthManipNC:\tbottomLayer:\tbegin')
        Shape = nci.variables[var][:].shape
        print(var, Shape)
        bathy = -self.bathy + 10.  #(-5985m -> -3m)
        bathy = bathy[None, None, :]
        bathy = repeat(bathy, Shape[0], axis=0)
        bathy = repeat(bathy, Shape[1], axis=1)
        print('bathy', bathy.shape)
        print('zbnd', self.zbnd.shape)
        #=where((zbnd[:,:,:,:,0]>=-10.)*(zbnd[:,:,:,:,1]<=-10.),MLDV,-1e30)
        arr = where((self.zbnd[:, :, :, :, 0] >= bathy) *
                    (self.zbnd[:, :, :, :, 1] <= bathy),
                    nci(var)[:], -1e30)
        arr = arr.max(1)
        arr = arr[:, None, :]
        arr = masked_where(abs(arr) > 1.e19, arr)
        print(var, arr.shape, arr.max(), arr.min())
        return arr


def todaystr():
    # returns string: DD/MM/YYYY on todays date
    return str(date.today().day) + '/' + str(date.today().month) + '/' + str(
        date.today().year)


def intersection(list1, list2):
    """return the overlap between two lists"""
    return list(set(list1).intersection(set(list2)))
