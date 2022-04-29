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
"""
.. module:: Coordinates Lists
   :platform: Unix
   :synopsis: A set of lists for the code to understand the differences 
   	      between dimensions, with lists of lats, lon, depths, etc.
.. moduleauthor:: Lee de Mora <ledm@pml.ac.uk>

"""


def completeList(alist):
    """
	This takes a list of strings, and extends the list
	with lots of permutations of uppercase/lowercase/title formattings.
	"""
    alist = list(alist)
    alist.extend([l + 's' for l in alist])  # add an s on the end, just in case
    list2 = alist[:]
    list2.extend([t.lower() for t in alist])
    list2.extend([t.upper() for t in alist])
    list2.extend([t.title() for t in alist])
    d = {l: 1 for l in list2}  # cast to a dict to remove duplicates
    return sorted(d.keys())


depthNames = [
    'depth',
    'DEPTH',
    'Depth',
    'deptht',
    'depthu',
    'depthv',
    'depthw',
    'lev',
    'nav_lev',
    'level',
    'Pressure',
    'pressure',
    'PRESSURE',
    'index_z',
    'level',
    'z',
]

depthNames = completeList(depthNames)

timeNames = [
    'time',
    'date',
    'month',
    'index_t',
    'timePlot',
    'time_centered',
    'time_counter',
]
timeNames = completeList(timeNames)

latnames = ['lat', 'latitude', 'latbnd', 'nav_lat', 'y', 'lat', 'rlat', 'j']
latnames = completeList(latnames)

lonnames = [
    'lon',
    'longitude',
    'lonbnd',
    'nav_lon',
    'x',
    'lon',
    'rlon',
    'i',
]
lonnames = completeList(lonnames)

alwaysInclude = [
    'time',
    'lat',
    'lon',
    'Time',
    'Lat',
    'Lon',
    'TIME',
    'LAT',
    'LON',
    'latbnd',
    'lonbnd',
    'time_centered',
    'LONGITUDE',
    'LATITUDE',
    'DEPTH',
    'TIME',
    'MONTH',
    'Longitude',
    'Latitude',
    'Depth',
    'Time',
    'Month',
    'nav_lat',
    'nav_lon',
    'nav_lev',
    'deptht',
    'depthu',
    'depthv',
    'depthw',
    'lev',  #'deptht_bounds',
    'Pressure',
    'pressure',
    'PRESSURE',
    'time_counter',
    'depth',
    'latitude',
    'longitude',
    'month',
    'mask',
    'index',
    'index_x',
    'index_y',
    'index_z',
    'index_t',
    'level',
]
alwaysInclude.extend(timeNames)
alwaysInclude.extend(lonnames)
alwaysInclude.extend(latnames)
alwaysInclude.extend(depthNames)
alwaysInclude = completeList(alwaysInclude)
