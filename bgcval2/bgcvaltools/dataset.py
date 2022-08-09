#!/usr/bin/ipython -i

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
# ledm@pml.ac.uk / momm@pml.ac.uk

import numpy as np
import netCDF4
from sys import argv

#####
# I wish this class wasn't neccesairy.
# We had to do this because netCDF4.Dataset caused a stack trace when trying to load the field "filepath()".
# This class should work the same way as netCDF4.Dataset, except that this object contains
# a netCDF4.Dataset instead of inheriting from it. As such, there are two layers of datasets.
# This is based on Momme's ncdfView.py class


class dataset:

    def __init__(self, filename, readflag='r', Quiet=True, skip_file=False):
        self.__filename__ = filename
        self.netcdfPath = filename
        self.filename = filename
        self.skip_file = skip_file
        try:
            self.dataset = netCDF4.Dataset(filename, 'r')
            self.dataset = netCDF4.Dataset(filename, 'r', format='NETCDF4')
        except OSError as oserr:
            print(oserr)
            print(f"dataset:\tUnable to open {filename}")
            if oserr.errno == -101:
                print(f"File {filename} appears to be corrupted")
                if self.skip_file:
                    print("Skipping it.")
                    pass
                else:
                    print("This file can not be skipped, exiting, " \
                          "please check the file!")
                    raise
            

        #####
        # link to various fields, so that the user experience is similar.
        self.variables = self.dataset.variables
        self.dimensions = self.dataset.dimensions
        try:
            self.title = self.dataset.title
        except:
            self.title = ''
        try:
            self.name = self.dataset.name
        except:
            self.name = ''
        try:
            self.Conventions = self.dataset.Conventions
        except:
            self.Conventions = ''
        try:
            self.TimeStamp = self.dataset.TimeStamp
        except:
            self.TimeStamp = ''
        try:
            self.description = self.dataset.description
        except:
            self.description = ''
        try:
            self.ncattrs = self.dataset.ncattrs
        except:
            self.ncattrs = ''

        if not Quiet: print(self.__unicode__())

    def close(self):
        try:
            self.dataset.close()
        except:
            pass

    def __call__(self, varStr, Squeeze=False, Object=False):
        if Object:
            return self.variables[varStr]
        else:
            if Squeeze:
                a = self.variables[varStr][:].squeeze()
                return a
            else:
                a = self.variables[varStr][:]
                return a

    def __unicode__(self):
        infoStr='-----------------\n'+\
                   'netCDF Object:\n'+\
                   '-----------------'
        for key in self.ncattrs():
            infoStr += '\n\n' + key.encode("utf-8") + ':\t'
            try:
                infoStr += str(getattr(self, key))
            except:
                pass
        dimList = list(self.dimensions.items())
        dimList = [(key, dim, len(dim)) for key, dim in dimList]
        dimList.sort(cmp=lambda x, y: cmp(x[0], y[0]))
        for key, dim, size in dimList:
            infoStr += '\n\t' + key
            if dim.isunlimited():
                infoStr += '\tUNLIMITED => ' + str(size)
            else:
                infoStr += '\t' + str(size)
        infoStr += '\n\n' + 'Variables:\n'
        varList = list(self.variables.items())
        varList.sort(cmp=lambda x, y: cmp(x[0], y[0]))
        for key, var in varList:
            infoStr += '\n\t' + key + ':'
            for k in var.ncattrs():
                if k == 'long_name':
                    infoStr += '\t' + str(getattr(var, k))
                elif k == 'units':
                    infoStr += '\t' + '[' + str(getattr(var, k) + ']')
            infoStr += '\n\t\t' + str(var.dimensions) + '=' + str(var.shape)
            infoStr += '\t' + str(var.dtype)
        return infoStr


if __name__ == "__main__":
    fn = None
    nc = dataset(argv[1], Quiet=False)
