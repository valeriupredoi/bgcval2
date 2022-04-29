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

#from ncdfView import ncdfView
try:
    from netCDF4 import Dataset, default_fillvals
except:
    from netCDF4 import Dataset, _default_fillvals

from datetime import date
from getpass import getuser
from os.path import exists
from .alwaysInclude import alwaysInclude


class changeNC:
    """
  This class takes a netcdf filename in, a netcdf filename out and an AutoVivification dictionary for changes.
  it copies the netcdf, but makes the changes to metadata requrested in a dictionairy.

  For a given variable, units, long_name and variable name can be changed here.
  a lambda function can be assigned to a variable for simple conversions.

  Further details:
  	from changeNC import changeNC, AutoVivification
	#To add or change a netcdf attribute:
	av = AutoVivification()
	av['att']['Description'] = 'New description'

	to remove an attribute:
	av['att']['Description'] = ''

	to change a dimension name:
	av['dim']['oldDimensionName']= 'newDimensionName'

	to change a dimension's size:
	av['dim']['oldDimensionName']['newSize']= 15

	to change a variable, i.e. change old var,'t', to new var:'time'.
	av['oldVarName']['name']	='newVarName'
	av['oldVarName']['units']	='newVarUnits'
	av['oldVarName']['long_name']	='newVarLongName'
	av['oldVarName']['newDims']	=('x','y',z') # these need to exist aready 
	
	to add a new variable:
	av['newVar'][var]['name']
	av['newVar'][var]['newDims']	
	av['newVar'][var]['dtype']
	av['newVar'][var]['long_name']
	av['newVar'][var]['units']
	av['newVar'][var]['newData']

	to remove a variable from the file, set it's name to "False", "Remove", or "Delete"
	av['t']['name']='False'

	A new data array can be set with: 
	av['oldVarName']['newData'] = NewDataArray

	
	Alternatively, the values in the array can be manipulated with a lambda function: 
	av['t']['convert']=lambda t:t/2. # divide time by two.

	or a predefined function:
	def addOne(arr): return arr+1
	av['t']['convert'] =addOne

	The debug flag prints more when set to True.
	
	The datasetFormat flag allows the netcdf to be written in a different format,
	following the conventions of netCDF4.Dataset.	
  """

    def __init__(self,
                 filenameIn,
                 filenameOut,
                 av,
                 debug=True,
                 datasetFormat='NETCDF4'):
        self.fni = filenameIn
        self.fno = filenameOut
        self.av = av
        self.debug = debug
        self.datasetFormat = datasetFormat
        self.run()

    def run(self):
        # starting checks:
        if not self.av:
            print('changeNC:\tERROR:\tchanges to make are no use:', self.av)
            return
        if not exists(self.fni):
            print('changeNC:\tERROR:\tinputfile name does not exist:',
                  self.fni)
            return

        if self.debug: print('changeNC:\tINFO:\tOpening dataset:\t', self.fni)
        nci = Dataset(self.fni, 'r')  #Quiet =True)

        # create dataset and netcdf attributes.
        if self.debug:
            print('changeNC:\tINFO:\tCreating a new dataset:\t', self.fno)
        nco = Dataset(self.fno, 'w', format=self.datasetFormat)
        attributes = nci.ncattrs()
        attributes.extend(list(self.av['att']))
        for att in list(set(attributes)):
            attribute = ''
            try:
                if att in self.av['att']: attribute = self.av['att'][att]
                else: attribute = nci.getncattr(att)
                if self.debug:
                    print('changeNC:\tINFO:\tadding attribute: \t\"', att,
                          '\":\t', attribute)
            except:
                if self.debug:
                    print(
                        'changeNC:\twarning:\tThat attribute probably isn\'t using ASCII characters!'
                    )
            nco.setncattr(att, attribute)

        # create new dimensions:
        dimensions = list(nci.dimensions.keys())
        if len(list(self.av['dim'].keys())) > 0:
            for dim in list(self.av['dim'].keys()):
                if dim not in dimensions:
                    if self.debug:
                        print('changeNC:\tINFO:\t New dimension: ', dim,
                              'not in ', dimensions)
                    dimSize = self.av['dim'][dim]['dimSize']
                    nco.createDimension(dim, dimSize)
                    if self.debug:
                        print('changeNC:\tINFO:\tadding New dimension: ', dim,
                              '\t(', dimSize, ')')

        # manipulate old dimensions:
        for dim in dimensions:
            newDim = dim
            dimSize = len(nci.dimensions[dim])
            if nci.dimensions[dim].isunlimited(): dimSize = None
            if self.av['dim'][dim]['name']:
                newDim = self.av['dim'][dim]['name']
            if self.av['dim'][dim]['newSize']:
                dimSize = self.av['dim'][dim]['newSize']
            nco.createDimension(newDim, dimSize)
            if self.debug:
                print('changeNC:\tINFO:\tadding dimension: ', dim, '-->',
                      newDim, '\t(', dimSize, ')')

        # list of variables to save
        keys = list(nci.variables.keys())

        try:
            newVars = list(self.av['newVar'].keys())
        except:
            newVars = []

        #for var in keys:
        #	if self.av[var]['name']:newname = self.av[var]['name']
        #	if newname.lower() in ['false', 'none','remove', 'delete', 0]:
        #		if self.debug: print 'changeNC:\tINFO:\tremoving variable: ',var
        #		continue

        # create Variables:
        for var in newVars:
            dimensions = self.av['newVar'][var]['newDims']
            vartype = self.av['newVar'][var]['dtype']
            if self.debug:
                print('changeNC:\tINFO:\tadding new variable: ', var, '\t(',
                      dimensions, ')')
            nco.createVariable(var,
                               vartype,
                               tuple(dimensions),
                               zlib=True,
                               complevel=5)

        for var in keys:
            if var in newVars: continue
            newname = var
            if self.av[var]['name']: newname = self.av[var]['name']
            if newname.lower() in ['false', 'none', 'remove', 'delete', 0]:
                if self.debug:
                    print('changeNC:\tINFO:\tremoving variable: ', var)
                continue
            dimensions = list(nci.variables[var].dimensions)

            if self.av[var]['dtype'] == {}:
                vartype = nci.variables[var].dtype
            else:
                vartype = self.av[var]['dtype']

            for d, dim in enumerate(dimensions):
                if self.av['dim'][dim]['name']:
                    dimensions[d] = self.av['dim'][dim]['name']

            if self.av[var]['newDims']: dimensions = self.av[var]['newDims']

            nco.createVariable(newname,
                               vartype,
                               tuple(dimensions),
                               zlib=True,
                               complevel=5)
            if self.debug:
                print('changeNC:\tINFO:\tadding variable: ', var, '-->',
                      newname, '\t(', dimensions, ')', '\tvartype:', vartype)

        # Long Names:
        for var in newVars:
            nco.variables[var].long_name = self.av['newVar'][var]['long_name']
            if self.debug:
                print('changeNC:\tINFO:\tadding new long_name: ', var, '\t(',
                      self.av['newVar'][var]['long_name'], ')')

        for var in keys:
            if var in newVars: continue
            long_name = ''
            newname = var
            if self.av[var]['name']: newname = self.av[var]['name']
            if newname.lower() in ['false', 'none', 'remove', 'delete', 0]:
                continue
            if self.av[var]['long_name']: long_name = self.av[var]['long_name']
            else:
                try:
                    long_name = nci.variables[var].long_name
                except:
                    print('changeNC:\tWarning:\tNo long_name for ', var)
            if long_name: nco.variables[newname].long_name = long_name
            if self.debug:
                print('changeNC:\tINFO:\tadding long_name: ', var, '-->',
                      newname, '\t(', long_name, ')')

        # Units:
        for var in newVars:
            nco.variables[var].units = self.av['newVar'][var]['units']
            if self.debug:
                print('changeNC:\tINFO:\tadding units: ', var, '\t(',
                      self.av['newVar'][var]['units'], ')')
        for var in keys:
            if var in newVars: continue
            units = ''
            newname = var
            if self.av[var]['name']: newname = self.av[var]['name']
            if newname.lower() in ['false', 'none', 'remove', 'delete', 0]:
                continue
            if self.av[var]['units']: units = self.av[var]['units']
            else:
                try:
                    units = nci.variables[var].units
                except:
                    print('changeNC:\tWarning:\tNo units for ', var)
            if units: nco.variables[newname].units = units
            if self.debug:
                print('changeNC:\tINFO:\tadding units: ', var, '-->', newname,
                      '\t(', units, ')')

        # Fill Values:
        if self.debug:
            print('changeNC:\tINFO:\tAbout to start Filling with new data:',
                  newVars)
        for var in newVars:
            if self.debug:
                print('changeNC:\tINFO:\tFilling ', var, ' ...',
                      self.av['newVar'][var]['newData'].shape, 'with new data',
                      nco.variables[var][:].shape)
            nco.variables[var][:] = self.av['newVar'][var]['newData']
        for var in keys:
            if var in newVars: continue
            newname = var
            func = lambda x: x
            if self.av[var]['name']: newname = self.av[var]['name']
            if newname.lower() in ['false', 'none', 'remove', 'delete', 0]:
                continue
            if self.av[var]['convert']: func = self.av[var]['convert']
            if len(self.av[var]['newData']):
                if self.debug:
                    print('changeNC:\tINFO:\tAdded newData into ', var)
                arr = self.av[var]['newData']
            else:
                arr = nci.variables[var][:]
            if self.debug:
                print('changeNC:\tINFO:\tCopying ',
                      var,
                      ' ...',
                      newname,
                      arr.shape,
                      end=' ')
            nco.variables[newname][:] = func(arr)
            #if self.debug: print '->', nco.variables[newname][:].shape
        # Close netcdfs:
        nco.close()
        nci.close()
        print('changeNC:\tINFO:\tsuccessfully created:\t', self.fno)
        return


def todaystr():
    # returns string: DD/MM/YYYY on todays date
    return str(date.today().day) + '/' + str(date.today().month) + '/' + str(
        date.today().year)


class AutoVivification(dict):
    """Implementation of perl's autovivification feature."""

    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value


# class to auto add objects to dictionaries
# ie:
# a = AutoVivification()
# a['a']['b']='ab'
# a['a']['c']='ac'
# print a['a'],a['a']['b'],a['a']['c']
