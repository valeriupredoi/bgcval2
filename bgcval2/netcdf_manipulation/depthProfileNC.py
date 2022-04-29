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




# this file takes a netcdf in, a netcdf out and a list of variables to save, and a position to save it.
# it automatically saves the dimensions, and the header.

# based on pruneNC, but also takes a profile location to output.
# if the values are integers, it assumes them to be a grid location.
# if the values are doubles, it assumes them to be lon-lat coords.

# 

# this could be extended into a 2D slice maker, but no inclination at the moment

# This class takes an input netcdf filename, an output netcdf filename and a list of variables to keep. Usually ['var', 'lat','lon',time','depth'].
# It creates a new netcdf that contains only the variables that you gave it. 
# It can also perform depth integration, if needed. (Assuming that the 2nd column is depth, and that there's a 5D 'zbnd' variable in the input file.)

#from ncdfView import ncdfView
try:	from netCDF4 import Dataset, default_fillvals
except: from netCDF4 import Dataset, _default_fillvals
from datetime import date
from getpass import getuser
from os.path import exists
from numpy.ma import array as marray
try:
	from numpyXtns import point_inside_polygon
assert: print "numpyXtns: not available."

class depthProfileNC:
  def __init__(self, filenameIn, filenameOut, variables, profilelocation=[0,0],  debug=False):
  	#timemean=False,
	self.fni=filenameIn
	self.fno=filenameOut
	self.vars=variables
	self.loc=profilelocation
	#self.timemean=timemean		#take mean of entire time series. Ideal for turning a daily file into a monthly file.
	self.debug=debug
	self.run()

  def run(self):	
	if not self.vars:
		print 'depthProfileNC:\tERROR:\tvariables to save are no use:', self.vars
		return
	if not exists(self.fni):
		print 'depthProfileNC:\tERROR:\tinputfile name does not exists:', self.fni
		return
				
	nci = Dataset(self.fni,'r')#Quiet =True)
	
	#check for location in the map.
	# if the values are integers, it assumes them to be a grid location.
	# if the values are doubles, it assumes them to be lon-lat coords.
	if isinstance(self.loc[0], int) and isinstance(self.loc[1], int):
		self.slice = self.loc
	else:
		if len(nci.variables['latbnd'].shape) == 2:
			#ORCA grid
			la,lo = getNemoIndex(self.loc[0],self.loc[1], nci.variables['latbnd'], nci.variables['lonbnd'])
		else:	
			#1D latbnd,lonbnd 
			la = getIndex(self.loc[0], nci.variables['latbnd'])
			lo = getIndex(self.loc[1], nci.variables['lonbnd'])		
		self.slice = [la,lo]
	if self.debug: print 'Saving the location: ',self.slice

	# special key word to save all the variables
	if self.vars == 'all':
		self.vars = nci.variables.keys()
		
	#check that there are some overlap between input vars and nci:
	for v in self.vars:
		if v in nci.variables.keys():continue
		print 'depthProfileNC:\tERROR:\tvariable,' ,v,', not found in ',self.fni
		nci.close()
		return
	
	
	
	#create dataset and header.
	if self.debug: print 'depthProfileNC:\tINFO:\tCreating a new dataset:\t', self.fno
	nco = Dataset(self.fno,'w')
	for a in nci.ncattrs():
		if self.debug: print 'depthProfileNC:\tINFO:\tcopying attribute: \t\"'+a+'\":\t', nci.getncattr(a)
		nco.setncattr(a,nci.getncattr(a))	
	appendToDesc= 'Reprocessed on '+todaystr()+' by '+getuser()+' using depthProfileNC.py'
	try: nco.Notes = nci.Notes + '\n\t\t'+appendToDesc
	except: nco.Notes = appendToDesc
	
	# list of variables to save, assuming some conventions
	alwaysInclude = ['time', 'lat','lon', 'latbnd', 'lonbnd', 'longitude', 'nav_lat','nav_lat', 'latitude','time_counter', 'deptht',]
	timeNames = ['time', 'time_counter', 't']
	latNames = ['lat', 'latbnd','nav_lat','x','latitude',]
	lonNames = ['lon', 'lonbnd','nav_lon','y','longitude',]
	save =   list(set(nci.variables.keys()).intersection(set(alwaysInclude) ) ) 
	save = list(set(sorted(save + self.vars)))
	
	# create dimensions:
	for d in nci.dimensions.keys():
	  if d in timeNames:  nco.createDimension(d, None)
	  elif d in latNames: nco.createDimension(d, 1)
	  elif d in lonNames: nco.createDimension(d, 1)	  
	  else:		      nco.createDimension(d, len(nci.dimensions[d]))

	# create Variables:
	for var in save:  nco.createVariable(var, nci.variables[var].dtype, nci.variables[var].dimensions,zlib=True,complevel=5)

	# Long Names:
	for var in save: 
		try:  	long_name=nci.variables[var].long_name
		except:	
		  if self.debug: print 'depthProfileNC:\tWarning:\tNo long_name for ', var
		  long_name = var
		  
		#if self.timemean: long_name.replace('Daily', 'Monthly')	
		nco.variables[var].long_name=long_name
		if self.debug: print 'depthProfileNC:\t Adding long_name for ', var, long_name
		  
	# Units:
	for var in save: 
		try:  	nco.variables[var].units=nci.variables[var].units
		except: print 'depthProfileNC:\tWarning:\tNo units for ', var	
		
	# Fill Values:
	for var in save:
		if self.debug: print 'depthProfileNC:\tINFO:\tCopying ', var, ' ...'
		shape  = nci.variables[var].shape
		if len(shape) == 4:
			arr = nci.variables[var][:,:,self.slice[0],self.slice[1]]
			arr = arr[:,:,None,None] #add extra empty dimensions:
		elif len(nci.variables[var].shape) == 3:
			try:
				arr = nci.variables[var][:,self.slice[0],self.slice[1]]
				arr = arr[:,None,None] #add extra empty dimensions:
			except: 
				print 'depthProfileNC:\tERROR:\tCoping failed due to unusual shape', var,shape, self.slice
				return
		elif  len(nci.variables[var].shape) ==2: # or var in latNames or var in lonNames:
			arr = nci.variables[var][self.slice[0],self.slice[1]]
			arr = arr[None,None]
		else:
			arr = nci.variables[var][:]
		
		#if self.timemean and len(intersection(['time','t'], nci.variables[var].dimensions)):
		#	if self.debug: print 'depthProfileNC:\tInfo:\tSaving time averaged var:',var
		#	arr = marray([arr.mean(0),])
		#	while len(arr.shape) < len(nci.variables[var].dimensions): arr = marray(arr[None,:])
			
		if self.debug: 
			print 'depthProfileNC:\tInfo:\tSaving var:',var, arr.shape, 
			print '[nco is expecting:',nco.variables[var][:].shape,']\tdims:', nci.variables[var].dimensions
			# this is not a bug if its a time.
			#if nco.variables[var][:].shape != arr.shape:
		#	print 'depthProfileNC:\tERROR: \toutput array does not match requested array shape'
		#	print 'depthProfileNC:\tERROR:\tvar:',var, arr.shape, 
		#	print 'but output netcdf is expecting:',nco.variables[var][:].shape,'\tdims:', nci.variables[var].dimensions	, nco.variables[var].dimensions
		#	#assert False
			
		nco.variables[var][:] =arr

	# Close netcdfs:
	nco.close()
	nci.close()
	if self.debug: print 'depthProfileNC:\tINFO:\tsuccessfully created:\t', self.fno
	return


def getIndex(x, arr):
  	#takes a depth and an array of box boundaries, returns depth index 
  	#if -1, coordinates are outside grid.
  	x = float(x)
	if ((arr[:,0]-x)*(arr[:,1]-x)).min() >0:return -1
	return ((arr[:,0]-x)*(arr[:,1]-x)).argmin()
	
def getNemoIndex(lat,lon, latbnd, lonbnd):
	# takes a lat and long coordinate, an returns the position of that coordinate in the NemoERSEM (ORCA1) grid.

	if lat < latbnd.min(): return -1,-1
	if lat > latbnd.max(): return -1,-1
	if lon < lonbnd.min(): return -1,-1
	if lon > lonbnd.max(): return -1,-1
		
	a,b = latbnd.shape[0], latbnd.shape[1]
	for la in range(a):
	  if not (latbnd[la].min()<=lat<latbnd[la].max()):continue 
	  for lo in range(b):
	    if la ==0 and lo == 107: continue # this pixel is broken
	    lab = latbnd[la,lo]
	    lob = lonbnd[la,lo]
	    #for i,j in enumerate(lob): lob[i]= safeLon(j)
	    if lob[0]>lob[1] or lob[2]<lob[3]:lob = lob+[-360.,0.,0.,-360.]
	    if not (lob.min()<=lon<lob.max()):continue
	    if lat in lab or lon in lob: #Lat,long matches box exactly:
	       print la,lo,'\ton latbnd,lonbnd edge:', lat, lab, '\t:', lon, lob
	       return la,lo
	    if point_inside_polygon(lat,lon,zip(lab, lob)):

	    	# what about if its on the edge?
	       	#print 'getNemoIndex:point ',la,lo,'is inside polygon:',lat,lon ,zip(lab, lob)
	        return la,lo
	
	print 'getNemoIndex: \tWARNING:\t didnt find point ', lat,lon 
	return -1,-1	
		
	
def todaystr():
	# returns string: DD/MM/YYYY on todays date
	return str(date.today().day)+'/'+str(date.today().month)+'/'+str(date.today().year)
	
def intersection(list1,list2):
	"""return the overlap between two lists"""
	return list(set(list1).intersection(set(list2) ) )
