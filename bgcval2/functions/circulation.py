#
# Copyright 2017, Plymouth Marine Laboratory
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
.. module:: Circulation tools
   :platform: Unix
   :synopsis: Calculates the various current in the eORCA1 grid.
.. moduleauthor:: Lee de Mora <ledm@pml.ac.uk>
"""

import numpy as np
from bgcval2.bgcvaltools.dataset import dataset
from bgcval2.UKESMpython import maenumerate
# coordinates of Drake Passage in eORCA1

LON=219
LAT0=79
LAT1=109

latslice26N = slice(227,228)
latslice26Nnm = slice(228,229)
latslice32S = slice(137,138)

umask_drake     = 0
e2u_drake    = 0
e3v_AMOC26N     = 0
e1v_AMOC26N    = 0
tmask_AMOC26N    = 0
alttmask_AMOC26N= 0
loadedArea     = False
loadedAltMask     = False


def loadDataMask(gridfn,maskname,):
    global umask_drake
    global e2u_drake
    global loadedArea
    global e3v_AMOC26N
    global e1v_AMOC26N
    global tmask_AMOC26N    
    if isinstance(gridfn, list) and len(gridfn)==1:
        gridfn = gridfn[0]

    nc = dataset(gridfn,'r')        
    e2u_drake = nc.variables['e2u'][LAT0:LAT1,LON]
    umask_drake = nc.variables['umask'][:,LAT0:LAT1,LON]
    e3v_AMOC26N = nc.variables['e3v'][:,latslice26Nnm,:]   # z level height 3D
    e1v_AMOC26N = nc.variables['e1v'][latslice26Nnm,:]     #
    tmask_AMOC26N = nc.variables['tmask'][:,latslice26Nnm,:]    
    nc.close()
    loadedArea = True


def loadAtlanticMask(altmaskfile,maskname='tmaskatl',):
    global alttmask_AMOC26N
    global loadedAltMask
    nc = dataset(altmaskfile,'r')        
    alttmask_AMOC26N = nc.variables[maskname][latslice26Nnm,:]   # 2D Atlantic mask
 #       alttmask_AMOC26N[228,180:208]=0.
    nc.close()
    loadedAltMask = True


def drakePassage(nc,keys,**kwargs):
    """
    This function. 
    
    nc: a netcdf openned as a dataset.
    keys: a list of keys to use in this function.
    
    """
    
    if 'areafile' not in list(kwargs.keys()):
        raise AssertionError("drakePassage:\t Needs an `areafile` kwarg to run calculation.")    

    try:    maskname = kwargs['maskname']
    except:    maskname = 'tmask'
            
    if not loadedArea: loadDataMask(kwargs['areafile'],maskname)
    try:    e3u = nc.variables['thkcello'][0,:,LAT0:LAT1,LON]    
    except:    e3u = nc.variables['e3u'     ][0,:,LAT0:LAT1,LON]

    try:    velo = nc.variables['uo' ][0,:,LAT0:LAT1,LON]    
    except:    velo = nc.variables['u3d'][0,:,LAT0:LAT1,LON]

    drake = np.sum(velo*e3u*e2u_drake*umask_drake)*1.e-6
    
    return drake



drakedetails     = {}
def loadDrakeDetails(gridfn):
    global drakedetails
    print("cmip5DrakePassage:\topening grid file", gridfn)    
    nc = dataset(gridfn,'r')
    if 'Drake_A' in list(nc.variables.keys()): 
        drakeKey = 'Drake_A'
    elif 'Drake' in list(nc.variables.keys()): 
        drakeKey = 'Drake'
            
    Drake_A = nc.variables[drakeKey][:]        
    tmask     = nc.variables['tmask'][:]
    drakedetails[(gridfn,'Drake_A')] = Drake_A
    drakedetails[(gridfn,'tmask')]     = tmask
    nc.close()
       
 
def cmip5DrakePassage(nc,keys,**kwargs):
    if 'gridfile' not in list(kwargs.keys()):
        raise AssertionError("drakePassage:\t Needs an `gridFile` kwarg to run calculation.")    
    gridFile_u =     kwargs['gridfile']
    try:    
        drakexsectArea     = drakedetails[(gridFile_u,'Drake_A')]
        tmasku        = drakedetails[(gridFile_u,'tmask')]
    except:
        loadDrakeDetails(gridFile_u)
        drakexsectArea     = drakedetails[(gridFile_u,'Drake_A')]
        tmasku        = drakedetails[(gridFile_u,'tmask')]

    drakexsectArea = np.ma.masked_where((drakexsectArea==0.) + (tmasku==1),drakexsectArea)
    velo = nc.variables[keys[0]]#[:]
    out = []
    if velo.ndim==4:
        for t in np.arange(velo.shape[0]):
            #print 'cmip5DrakePassage',t, velo[t].shape, drakexsectArea.shape, tmasku.shape, gridFile_u
            drk = velo[t] *drakexsectArea 
            
            drk = np.ma.masked_where((drk==0.) + (tmasku==1.) +drk.mask, drk)
            xsection = drk.sum(2)*1.e-6
            totaldrk = drk.sum()*1.e-6                        
            plots=0
            if plots:                        
                from matplotlib import pyplot
                print(drk.shape,xsection.shape, velo.shape, xsection.sum(),totaldrk)                        
                ax1 = pyplot.subplot(411)
                pyplot.pcolormesh(xsection[::-1,:])
                pyplot.title('flux: '+str(xsection.sum())+' Sv, or:'+str(totaldrk))
                pyplot.colorbar()
                
                ax4 = pyplot.subplot(412)
                pyplot.plot(xsection[::-1,:].sum(0))
                pyplot.title('total: '+str(xsection.sum())+' Sv')    
                                    
                ax2 = pyplot.subplot(413)
                pyplot.pcolormesh(drakexsectArea.mean(2)[::-1,:]*1.e-6)
                pyplot.title('Crossectional Area')
                pyplot.colorbar()
                
                ax3 = pyplot.subplot(414)
                v = np.ma.masked_where(drakexsectArea.mask,velo[t])
                pyplot.pcolormesh(v.sum(2)[::-1,:])
                pyplot.title('Velocity')
                pyplot.colorbar()
                                                                
                pyplot.show()

            out.append(totaldrk)
    else:
        assert 0

    out = np.ma.array(out)
    print("cmip5DrakePassage:\tdrake",keys, out, out.shape,tmasku.shape,velo.shape,out.mean())
    return out    # should return 1 d time array.
            

amocdetails     = {}
def loadAMOCdetails(gridfn):
    global amocdetails
    print("loadAMOCdetails:\topening grid file", gridfn)    
    nc = dataset(gridfn,'r')
    if 'AMOC_26N_A' in list(nc.variables.keys()): 
        AMOCkey = 'AMOC_26N_A'
    elif 'AMOC' in list(nc.variables.keys()): 
        AMOCkey = 'AMOC'
    AMOC_26N_A = nc.variables[AMOCkey][:]        
    
    #if 'ADRC_16N_A' in nc.variables.keys(): 
    #    ADRC_key = 'ADRC_16N_A'
    #ADRC_16N_A = nc.variables[ADRC_key][:]            
    
    tmaskv     = nc.variables['tmask'][:]
    amocdetails[(gridfn,'AMOC_26N_A')] = AMOC_26N_A
    #amocdetails[(gridfn,'ADRC_16N_A')] = ADRC_16N_A    
        
    amocdetails[(gridfn,'tmask')]     = tmaskv
    nc.close()
    
    
def cmip5AMOC(nc,keys,**kwargs):
    if 'gridfile' not in list(kwargs.keys()):
        raise AssertionError("cmip5AMOC:\t Needs an `gridFile` kwarg to run calculation.")    
    gridFile_v =     kwargs['gridfile']

    try:    
        xsectArea     = amocdetails[(gridFile_v,'AMOC_26N_A')]
        tmaskv        = amocdetails[(gridFile_v,'tmask')]
    except:
        loadAMOCdetails(gridFile_v)
        xsectArea     = amocdetails[(gridFile_v,'AMOC_26N_A')]
        tmaskv        = amocdetails[(gridFile_v,'tmask')]
        
    xsectArea = np.ma.masked_where((xsectArea==0.) + (tmaskv==1),xsectArea)
    print(gridFile_v,xsectArea.shape,  xsectArea.sum())
    #print nc.__filename__
    velo = nc.variables[keys[0]][:]/1.E06#[:]
    out = []
    
    if velo.ndim==4:
        zlevs = velo.shape[1]
        tlen  = velo.shape[0]
        altshape = velo[0,:,:,0].shape
        for t in np.arange(tlen):
            atlmoc = np.zeros(altshape)
            vel = velo[t] * xsectArea
            vel = np.ma.masked_where((vel==0.)+(tmaskv == 1),vel)
            atlmoc -= vel.sum(2)

            for z in np.arange(zlevs-2,-1,-1):
                atlmoc[z,:] = atlmoc[z+1,:] + atlmoc[z,:]        
            out.append(atlmoc.max())
            #print "cmip5AMOC:", t,atlmoc.max()
    else:
        assert 0                    
    out = np.ma.array(out)
    print("cmip5AMOC: ",keys, out, out.shape,tmaskv.shape,velo.shape,out.mean())
    return out    # should return 1 d time array.


def cmip5ADRC(nc,keys,**kwargs):
    if 'gridfile' not in list(kwargs.keys()):
        raise AssertionError("cmip5ADRC:\t Needs an `gridFile` kwarg to run calculation.")    
    gridFile_v =     kwargs['gridfile']

    try:    
        xsectArea     = amocdetails[(gridFile_v,'AMOC_26N_A')]
        tmaskv        = amocdetails[(gridFile_v,'tmask')]
    except:
        loadAMOCdetails(gridFile_v)
        xsectArea     = amocdetails[(gridFile_v,'AMOC_26N_A')]
        tmaskv        = amocdetails[(gridFile_v,'tmask')]
        
    xsectArea = np.ma.masked_where((xsectArea==0.) + (tmaskv==1),xsectArea)

    velo = nc.variables[keys[0]][:]/1.E06#[:]
    print(xsectArea.shape,  xsectArea.sum())
    out = []
    
    if velo.ndim==4:
        zlevs = velo.shape[1]
        tlen  = velo.shape[0]
        altshape = velo[0,:,:,0].shape
        for t in np.arange(tlen):
            atlmoc = np.zeros(altshape)
            vel = velo[t] * xsectArea
            vel = np.ma.masked_where((vel==0.)+(tmaskv == 1),vel)
            atlmoc -= vel.sum(2)

            for z in np.arange(zlevs-2,-1,-1):
                atlmoc[z,:] = atlmoc[z+1,:] + atlmoc[z,:]        
            out.append(atlmoc.min())
            #print "cmip5AMOC:", t,atlmoc.max()
    else:
        assert 0                    
    out = np.ma.array(out)
    print("cmip5ADRC",keys, out, out.shape,tmaskv.shape,velo.shape,out.mean())
    return out    # should return 1 d time array.
   
 
def cmip5ADRC16(nc,keys,**kwargs):
    if 'gridfile' not in list(kwargs.keys()):
        raise AssertionError("cmip5ADRC:\t Needs an `gridFile` kwarg to run calculation.")    
    gridFile_v =     kwargs['gridfile']

    try:    
        xsectArea     = amocdetails[(gridFile_v,'ADRC_16N_A')]
        tmaskv        = amocdetails[(gridFile_v,'tmask')]
    except:
        loadADRCdetails(gridFile_v)
        xsectArea     = amocdetails[(gridFile_v,'ADRC_16N_A')]
        tmaskv        = amocdetails[(gridFile_v,'tmask')]
        
    xsectArea = np.ma.masked_where((xsectArea==0.) + (tmaskv==1),xsectArea)

    velo = nc.variables[keys[0]][:]/1.E06#[:]
    print(xsectArea.shape,  xsectArea.sum())
    out = []
    
    if velo.ndim==4:
        zlevs = velo.shape[1]
        tlen  = velo.shape[0]
        altshape = velo[0,:,:,0].shape
        for t in np.arange(tlen):
            atlmoc = np.zeros(altshape)
            vel = velo[t] * xsectArea
            vel = np.ma.masked_where((vel==0.)+(tmaskv == 1),vel)
            atlmoc -= vel.sum(2)

            for z in np.arange(zlevs-2,-1,-1):
                atlmoc[z,:] = atlmoc[z+1,:] + atlmoc[z,:]        
            out.append(atlmoc.min())
            #print "cmip5AMOC:", t,atlmoc.max()
    else:
        assert 0                    
    out = np.ma.array(out)
    print("cmip5ADRC16",keys, out, out.shape,tmaskv.shape,velo.shape,out.mean())
    return out    # should return 1 d time array.
        
                
def TwentySixNorth(nc,keys,**kwargs):
    """
    This function loads the AMOC/ADRC array that is used for eORCA
    
    nc: a netcdf openned as a dataset.
    keys: a list of keys to use in this function.
    
    """
    
    if 'areafile' not in list(kwargs.keys()):
        raise AssertionError("AMOC26N:\t Needs an `areafile` kwarg to run calculation.")    
    try:
        maskname = kwargs['maskname']
    except:
        maskname = 'tmask'
    if not loadedArea:
        loadDataMask(kwargs['areafile'], maskname)
    
    if 'altmaskfile' not in list(kwargs.keys()):
        raise AssertionError("AMOC26N:\t Needs an `altmaskfile` kwarg to run calculation.")    
    try:
        altmaskfile = kwargs['altmaskfile']
    except:
        altmaskfile = 'data/basinlandmask_eORCA1.nc'

    if not loadedAltMask: 
        loadAtlanticMask(altmaskfile, maskname='tmaskatl')
    
    zv = np.ma.array(nc.variables[keys[0]][..., latslice26Nnm, :]) # m/s
    atlmoc = np.array(np.zeros_like(zv[0, :, :, 0]))
    e2vshape = e3v_AMOC26N.shape
    xsectArea = (e1v_AMOC26N * e3v_AMOC26N)
    print('TwentySixNorth:', e1v_AMOC26N.shape, e3v_AMOC26N.shape, xsectArea.shape, xsectArea.sum())
    TotalXsection = 0
    for la in range(e2vshape[1]):           #ji, y
        for lo in range(e2vshape[2]):         #jj , x,
            if int(alttmask_AMOC26N[la, lo]) == 0:
                continue
            for z in range(e2vshape[0]):        # jk
                if int(tmask_AMOC26N[z, la, lo]) == 0:
                    continue
                if np.ma.is_masked(zv[0, z, la, lo]):
                    continue
                atlmoc[z, la] = atlmoc[z, la] - e1v_AMOC26N[la, lo]*e3v_AMOC26N[z, la, lo]*zv[0, z, la, lo]/1.E06
                TotalXsection += e1v_AMOC26N[la, lo]*e3v_AMOC26N[z, la, lo]

    print("TotalXsection:", TotalXsection)
    # Cumulative sum from the bottom up.
    for z in range(73, 1, -1):
        atlmoc[z, :] = atlmoc[z+1, :] + atlmoc[z, :]
    return atlmoc
       
 
def AMOC26N(nc, keys, **kwargs):
    atlmoc = TwentySixNorth(nc, keys, **kwargs)
    return atlmoc.max()


def ADRC26N(nc,keys,**kwargs):
    atlmoc = TwentySixNorth(nc, keys, **kwargs)
    return atlmoc.min()
        
                
                        
