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
from netCDF4 import Dataset, num2date, date2num
try:
    from netCDF4 import default_fillvals
except:
    from netCDF4 import _default_fillvals as default_fillvals
from datetime import date
from getpass import getuser
from os.path import exists
from numpy.ma import array, masked_all
from numpy import append, mean, int32, int16
import numpy as np
from glob import glob
from .alwaysInclude import alwaysInclude as alwaysIncludList, timeNames

# a new comment to test github


class mergeNC:
    """ This class takes a list of input netcdf filenames, an output netcdf filename and a list of variables.
   it automatically saves the dimensions, and the header.
   It creates a new netcdf that contains the variables that you gave it from all the files in the input file list in chronological.
  """

    def __init__(self,
                 filesIn,
                 filenameOut,
                 variables,
                 timeAverage=False,
                 debug=False,
                 calendar='standard',
                 fullCheck=False):
        self.fnsi = filesIn
        self.fno = filenameOut
        self.vars = variables
        self.cal = calendar
        self.timeAverage = timeAverage
        self.fullCheck = fullCheck
        self.debug = debug
        self.run()

    def run(self):
        if type(self.fnsi) == type('abc'):
            self.fnsi = glob(self.fnsi)

        if not exists(self.fnsi[0]):
            print('mergeNC:\tERROR:\tinputfile name does not exists:',
                  self.fnsi[0])
            return
        if self.debug:
            print('mergeNC:\tINFO:\topening dataset:\t', self.fnsi[0])
        nci = Dataset(self.fnsi[0], 'r')  #Quiet =True)

        if self.timeAverage:
            print('mergeNC:\tWARNING:\ttimeAverage is not yet debugged. '
                  )  # are no use:', self.vars
            #return
        if not self.vars:
            if self.debug:
                print(
                    'mergeNC:\tINFO:\tvariables to save are empty, saving all.'
                )
            self.vars = list(nci.variables.keys())

        if self.vars == 'all':
            if self.debug:
                print(
                    'mergeNC:\tINFO:\tvariables to save:  \'all\' requested. ')
            self.vars = list(nci.variables.keys())

        if self.cal != 'standard':
            if self.debug:
                print('mergeNC:\tINFO:\tUsing non-standard calendar:',
                      self.cal)

        #check that there are some overlap between input vars and nci:
        for v in self.vars:
            if v in list(nci.variables.keys()): continue
            print('mergeNC:\tERROR:\tvariable,', v, ', not found in ',
                  self.fnsi[0])
            return

        #create dataset and header.
        if self.debug:
            print('mergeNC:\tINFO:\tCreating a new dataset:\t', self.fno)
        nco = Dataset(self.fno, 'w')
        for a in nci.ncattrs():
            try:
                if self.debug:
                    print(
                        'mergeNC:\tINFO:\tcopying attribute: \t\"' + str(a) +
                        '\":\t', str(nci.getncattr(a)))
                nco.setncattr(a, nci.getncattr(a))
            except:
                if self.debug:
                    print(
                        'changeNC:\twarning:\tThat attribute probably isn\'t using ASCII characters!'
                    )
        appendToDesc = 'Reprocessed on ' + todaystr() + ' by ' + getuser(
        ) + ' using mergeNC.py'
        try:
            nco.Notes = nci.Notes + '\n\t\t' + appendToDesc
        except:
            nco.Notes = appendToDesc

        # list of variables to save, assuming some conventions
        #	alwaysInclude = ['time', 'lat','lon', 'latbnd', 'lonbnd', 'latitude', 'Latitude', 'longitude', 'Longitude',
        #			 't','nav_lat','nav_lon', 'time_counter','time_centered',
        #			 'deptht','depth','depthu','depthv', 'depthw','z','month','bathymetry','Depth','deptht_bounds',
        #			  'lat_bnds',  'lon_bnds', 'depth_bnds',
        #			  'ensemble',]
        #	alwaysInclude = intersection(nci.variables.keys(),alwaysInclude)
        alwaysInclude = intersection(list(nci.variables.keys()),
                                     alwaysIncludList)
        save = list(set(sorted(alwaysInclude + self.vars)))
        time = intersection([
            'time',
            't',
            'time_counter',
            'month',
            'time_centered',
        ], alwaysInclude)
        tvars = []
        if len(time) == 1:
            tvar = time[0]
            tvars = [
                tvar,
            ]
        elif not len(time):
            tvar = 'time'
            tvars = [
                tvar,
            ]
        else:
            tvar = time[0]
            tvars = time

        missing = {}
        if self.fullCheck:
            if self.debug:
                print(
                    'mergeNC:\tINFO:\tPerforming full check for missing entries'
                )
            for t, fni in enumerate(self.fnsi):
                #if self.debug: print 'mergeNC:\tINFO:\tOpening ', fni, ' ...', t
                nci = Dataset(fni, 'r')
                keys = list(nci.variables.keys())
                for s in save:
                    if s in alwaysInclude: continue
                    if s not in keys:
                        print('mergeNC:\tWARNING:\tFull check: ', s,
                              ' is missing from ', fni)
                        try:
                            missing[s].append(fni)
                        except:
                            missing[s] = [
                                fni,
                            ]
                nci.close()

            for s in list(missing.keys()):
                #remove key:
                #print 'mergeNC:\tWARNING:\tFull check:\tremoving',s,' from ',save
                #save.remove(s)

                #remove missing files:
                for fn in missing[s]:
                    print('mergeNC:\tWARNING:\tFull check:\tremoving', fni,
                          ' from files')
                    try:
                        self.fnsi.remove(fn)
                    except:
                        print('mergeNC:\tWARNING:\tFull check:\t', fni,
                              ' already removed from files')

        # create dimensions:
        nci = Dataset(self.fnsi[0], 'r')  #Quiet =True)
        for d in list(nci.dimensions.keys()):
            if nci.dimensions[d].isunlimited() or d.lower() in [
                    'time', 'time_counter', time
            ]:
                dimSize = None
            else:
                dimSize = len(nci.dimensions[d])
            nco.createDimension(d, dimSize)
            if self.debug:
                print('mergeNC:\tINFO:\tCreating Dimension:', d, dimSize)

        # create Variables:
        for var in save:
            dt = nci.variables[var].dtype

            if dt in [
                    int32([
                        5,
                    ]).dtype,
            ]: dfkey = 'i8'
            elif dt in [
                    int16([
                        5,
                    ]).dtype,
            ]: dfkey = 'i4'
            else: dfkey = 'f8'

            if self.debug:
                print('mergeNC:\tINFO:\tCreating Variable:',
                      var,
                      dt,
                      nci.variables[var].dimensions,
                      end=' ')
                print("zlib=True,complevel=5,fill_value=",
                      default_fillvals[dfkey], dfkey)

            nco.createVariable(var,
                               dt,
                               nci.variables[var].dimensions,
                               zlib=True,
                               complevel=5,
                               fill_value=default_fillvals[dfkey])
        #try:	nco.createVariable(var, dt, nci.variables[var].dimensions,zlib=True,complevel=5,fill_value=default_fillvals[dfkey])
        #except: nco.createVariable(var, dt, nci.variables[var].dimensions,zlib=True,complevel=5,fill_value=default_fillvals[dfkey])

        # Long Names:
        for var in save:
            try:
                nco.variables[var].long_name = nci.variables[var].long_name
            except:
                if self.debug:
                    print('mergeNC:\tWarning:\tNo long_name for ', var)

        # Units:
        for var in save:
            #if var in time and self.timeAverage: nco.variables[var].units='Month'
            #else:
            try:
                nco.variables[var].units = nci.variables[var].units
            except:
                if self.debug: print('mergeNC:\tWarning:\tNo units for ', var)

        # Fill Values:
        for var in alwaysInclude:
            #if var in time:continue
            if var in tvars:
                continue  # there may be more than one time variable: ie time and month.
            if self.debug:
                print('mergeNC:\tINFO:\tCopying ', var, ' ...',
                      nci.variables[var][:].shape)
            try:
                nco.variables[var][:] = nci.variables[var][:].data
            except:
                nco.variables[var][:] = nci.variables[var][:]
        nci.close()

        a = {}
        for t in tvars:
            a[t] = []
        for var in save:
            if var in alwaysInclude: continue
            a[var] = []

        for t, fni in enumerate(self.fnsi):
            if self.debug: print('mergeNC:\tINFO:\tOpening ', fni, ' ...', t)
            nci = Dataset(fni, 'r')

            #times:
            for t in tvars:
                try:
                    tval = num2date(nci.variables[t][:],
                                    nci.variables[t].units,
                                    calendar=self.cal)
                    a[t].extend(
                        date2num(tval,
                                 nco.variables[t].units,
                                 calendar=self.cal))
                except:
                    a[t].extend(nci.variables[t][:])

            if self.debug:
                print('mergeNC:\tINFO:\tTIME:', t, tvar, array(a[tvar]).shape)

            # not time:
            for var in list(a.keys()):
                #if var in time:continue
                if var in tvars:
                    continue  # there may be more than one time variable: ie time and month.
                if var in list(nci.variables.keys()):
                    arr = nci.variables[var][:]
                else:
                    print('mergeNC:\tWARNING:', fni, ' is missing variable:',
                          var, nco.variables[var][0, :].shape)
                    arr = masked_all(nco.variables[var][0, :].shape)
                if not self.timeAverage:
                    if not len(a[var]): a[var] = arr
                    else: a[var] = append(a[var], arr, axis=0)
                else:
                    if not len(a[var]): a[var] = arr
                    else: a[var] = append(a[var], arr, axis=0)

                if self.debug:
                    print('mergeNC:\tINFO\tvar:', t, var, 'len:', len(a[var]),
                          arr.shape,
                          array(a[var]).shape)

            nci.close()

        if self.timeAverage:
            for var in list(a.keys()):
                if self.debug: print("mergeNC:\tINFO\tTime Average:", var)
                if var in timeNames:
                    nco.variables[tvar][:] = [
                        array(a[var]).mean(),
                    ]
                else:
                    #nco.variables[var][:] = array(a[var])[None,:]/float(len(self.fnsi)) # assumes one month per file.
                    if self.debug:
                        print("mergeNC:\tINFO\tTime Average: shape shift: ",
                              var)
                    if self.debug:
                        print("mergeNC:\tINFO\tTime Average: shape shift: ",
                              array(a[var]).shape, '-->',
                              array(a[var]).mean(0)[None, :].shape)
                    try:

                        timeAverageArr = np.ma.array(a[var])
                        timeAverageArr = np.ma.masked_where(
                            timeAverageArr.mask + (timeAverageArr > 9.969e+36),
                            timeAverageArr)
                        nco.variables[var][:] = timeAverageArr.mean(0)[None, :]
                    except:
                        assert 0
                        if tvar not in nco.variables[var].dimensions:
                            print('mergeNC:\tWARNING:', var,
                                  ' has no time dimension:',
                                  nco.variables[var].dimensions)
                            nco.variables[var][:] = array(a[var])
                    if self.debug:
                        print("mergeNC:\tINFO\tTime Average:", var,
                              nco.variables[var][:].shape,
                              array(a[var]).shape,
                              nco.variables[var].dimensions)
                    if self.debug:
                        print("mergeNC:\tINFO\tTime Average: min-max range",
                              var, nco.variables[var][:].min(), '-->',
                              nco.variables[var][:].max())

        else:  # No time averaging.
            for var in list(a.keys()):
                if self.debug:
                    print('mergeNC:\tINFO:\tsaving ', var, ' ...',
                          nco.variables[var][:].shape,
                          array(a[var]).shape,
                          nco.variables[var].dimensions)  #, a[var][0]
                nco.variables[var][:] = array(a[var])

        # Close output netcdfs:
        nco.close()
        if self.debug:
            print('mergeNC:\tINFO:\tsuccessfully created:\t', self.fno)
        return


def todaystr():
    # returns string: DD/MM/YYYY on todays date
    return str(date.today().day) + '/' + str(date.today().month) + '/' + str(
        date.today().year)


def intersection(list1, list2):
    return list(set(list1).intersection(set(list2)))
