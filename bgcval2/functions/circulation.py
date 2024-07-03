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
# You should have received a copy of the revised BSD license aeORCA1_drake_LONg with bgc-val.
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
import os
import numpy as np
import math
#import itertools

from bgcval2.bgcvaltools.dataset import dataset
from bgcval2.bgcvaltools.bv2tools import maenumerate
from bgcval2.functions.get_kwarg_file import get_kwarg_file

# coordinates of Drake Passage in eORCA1
# drake passage is:
# -68.5 (W), -67.2 to -52.6 South


eORCA1_drake_LON=219
eORCA1_drake_LAT0=79
eORCA1_drake_LAT1=109

eORCA025_drake_LON=875
eORCA025_drake_LAT0=317
eORCA025_drake_LAT1=436 

eORCA025_AEU_LON=1057
eORCA025_AEU_LAT0=664
eORCA025_AEU_LAT1=704


# AMOC:
# Coordinates of AMOC calc: 
# 26.5 N
# florida = -81. # W
# Western Sahara = -11. #W


#eORCA1_latslice26N = slice(227,228)
eORCA1_latslice26Nnm = slice(228,229)
eORCA1_latslice32S = slice(137,138)

eORCA025_latslice26Nnm = slice(794,795)
# eORCA025_lonslice26Nnm = slice(825,1105) # Florida
eORCA025_lonslice26Nnm = slice(841, 1105) # Bahamas (removes Florida Straight)

umask_drake = 0
e2u_drake = 0
e2u_aeu = 0
umask_aeu = 0
e3v_AMOC26N = 0
e1v_AMOC26N = 0
tmask_AMOC26N = 0
alttmask_AMOC26N = 0
alttmask = 0
loadedArea = False
loadedAltMask = False
loadedAltMask_full = False
loaded_AEU = False




def maenumerate(marr):
    mask = ~marr.mask.ravel()
    for i, m in zip(np.ndenumerate(marr), mask):
        if m: yield i

def myhaversine(lon1, lat1, lon2, lat2):
    """
            Calculate the great circle distance between two points
            on the earth (specified in decimal degrees)
        """
    # convert decimal degrees to radians
    [lon1, lat1, lon2, lat2] = [math.radians(l) for l in [lon1, lat1, lon2, lat2]]

    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = math.sin(dlat / 2.)**2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2.)**2
    c = 2. * math.asin(math.sqrt(a))
    dist = 6367000. * c

    return dist



def loadDataMask(gridfn, maskname, grid):
    """
    Load files and some key fields:
    """
    global umask_drake
    global e2u_drake
    global e3v_AMOC26N
    global e1v_AMOC26N
    global tmask_AMOC26N
    global loadedArea
    global loadedAltMask

    if grid == 'eORCA1':
        LON = eORCA1_drake_LON
        LAT0 = eORCA1_drake_LAT0
        LAT1 = eORCA1_drake_LAT1
        latslice26Nnm = eORCA1_latslice26Nnm
    elif grid == 'eORCA025':
        LON = eORCA025_drake_LON
        LAT0 = eORCA025_drake_LAT0
        LAT1 = eORCA025_drake_LAT1
        latslice26Nnm = eORCA025_latslice26Nnm
        lonslice26N = eORCA025_lonslice26Nnm
    else:
        assert 0
    print('circulation, loadDataMask, loading:', gridfn)
    nc = dataset(gridfn, 'r')       
    e2u_drake = nc.variables['e2u'][..., LAT0:LAT1, LON]
    umask_drake = nc.variables['umask'][..., LAT0:LAT1, LON]
    print('circulation loadDataMask:', gridfn, nc.variables['e2u'].shape, nc.variables['umask'].shape, e2u_drake.shape, umask_drake.shape)

    if 'e3v_0' in nc.variables.keys():
        # print('amoc:', nc.variables['e3v_0'].shape, nc.variables['e3v_0'].min(), nc.variables['e3v_0'].max()) 
        e3v_AMOC26N = nc.variables['e3v_0'][..., latslice26Nnm, lonslice26N].squeeze()   # z level height 3D 
        e1v_AMOC26N = nc.variables['e1v'][..., latslice26Nnm, lonslice26N]     #
        tmask_AMOC26N = nc.variables['tmask'][..., latslice26Nnm, lonslice26N]
    else:
        e3v_AMOC26N = nc.variables['e3v'][..., latslice26Nnm, :]   # z level height 3D
        e1v_AMOC26N = nc.variables['e1v'][..., latslice26Nnm, :]     #
        tmask_AMOC26N = nc.variables['tmask'][..., latslice26Nnm, :]

    print('e3v_AMOC26N:', e3v_AMOC26N, latslice26Nnm, e3v_AMOC26N.shape)
    nc.close()
    loadedArea = True


def load_AEU_masks(gridfn, grid):
    """
    Loads the AEU cell lenght in the y (j, Meridional, North-South) direction. 
    """
    global e2u_aeu
    global loaded_AEU
    if grid == 'eORCA025':
        LON = eORCA025_AEU_LON
        LAT0 = eORCA025_AEU_LAT0
        LAT1 = eORCA025_AEU_LAT1
    else: assert 0
 
    nc = dataset(gridfn, 'r')
    e2u_aeu = nc.variables['e2u'][..., LAT0:LAT1, LON]
    nc.close()
    loaded_AEU=True
    


def loadAtlanticMask(altmaskfile, maskname='tmaskatl', grid = 'eORCA1'):
    """
    Load the atlantic ocean mask.
    """
    global alttmask_AMOC26N
    global loadedAltMask
    if grid == 'eORCA1':
        latslice26Nnm = eORCA1_latslice26Nnm
    else:
        raise ValueError("Grid not recognised in this calculation: %s", grid)
    nc = dataset(altmaskfile, 'r')        
    alttmask_AMOC26N = nc.variables[maskname][latslice26Nnm, :]
    nc.close()
    loadedAltMask = True


def loadAtlanticMask_full(altmaskfile, maskname='tmaskatl', grid = 'eORCA1'):
    """
    Load the atlantic ocean mask.
    """
    global alttmask
    global loadedAltMask_full
    if grid == 'eORCA1':
        latslice26Nnm = eORCA1_latslice26Nnm
    else:
        raise ValueError("Grid not recognised in this calculation: %s", grid)
    nc = dataset(altmaskfile, 'r')        
    alttmask = nc.variables[maskname][:]
    nc.close()
    loadedAltMask_full = True


def find_keys_in_nc(nc, keys):
    """
    Find some keys in a netcdf.
    """
    intersection_keys = list(set(keys) & set(nc.variables.keys()))
    if not intersection_keys:
        raise KeyError(f'Unable to find keys {keys} in {nc.filename}')
    if len(intersection_keys) > 1:
        print('WARNING: found more than one key {intersection_keys} in {nc.filename}')
    return intersection_keys


def drakePassage(nc, keys, **kwargs):
    """
    This function calculates drake passage for eORCA1. 
    
    nc: a netcdf openned as a dataset.
    keys: a list of keys to use in this function.
    
    """
    areafile = get_kwarg_file(kwargs, 'areafile')
    maskname = kwargs.get('maskname', 'tmask')
    grid = kwargs.get('grid', 'eORCA1')

    if not loadedArea:
        loadDataMask(areafile, maskname, grid)

    if grid == 'eORCA1':
        LON = eORCA1_drake_LON
        LAT0 = eORCA1_drake_LAT0
        LAT1 = eORCA1_drake_LAT1
    elif grid == 'eORCA025':
        LON = eORCA025_drake_LON
        LAT0 = eORCA025_drake_LAT0
        LAT1 = eORCA025_drake_LAT1
        latslice26Nnm = eORCA025_latslice26Nnm
    else:
        assert 0
    print('drakePassage:', grid, 'LON', LON, 'LAT0',LAT0, 'LAT1', LAT1)

    all_e3u_keys = ['thkcello', 'e3u']
    e3u_keys = find_keys_in_nc(nc, all_e3u_keys)
    e3u = nc.variables[e3u_keys[0]][0, :, LAT0:LAT1, LON]


    all_velo_keys = ['uo', 'u3d']
    velo_keys = find_keys_in_nc(nc, all_velo_keys)
    velo = nc.variables[velo_keys[0]][0, :, LAT0:LAT1, LON]
    print('drakePassage:', grid, velo.shape, e3u.shape, e2u_drake.shape, umask_drake.shape)
    drake = np.sum(velo * e3u * e2u_drake * umask_drake) * 1.e-6
    return drake


def TwentySixNorth(nc,keys,**kwargs):
    """
    This function loads the AMOC/ADRC array that is used for eORCA
    
    nc: a netcdf openned as a dataset.
    keys: a list of keys to use in this function.
    
    """
    areafile = get_kwarg_file(kwargs, 'areafile')
    maskname = kwargs.get('maskname', 'tmask')
    grid = kwargs.get('grid', 'eORCA1')

    if not loadedArea:
        loadDataMask(areafile, maskname, grid)

    if grid == 'eORCA1':
        latslice26Nnm = eORCA1_latslice26Nnm

        altmaskfile = get_kwarg_file(kwargs, 'altmaskfile', default = 'bgcval2/data/basinlandmask_eORCA1.nc')
        if not loadedAltMask:
             loadAtlanticMask(altmaskfile, maskname='tmaskatl', grid=grid)
    elif grid == 'eORCA025':
        latslice26Nnm = eORCA025_latslice26Nnm

    else:
        # grid not recognised.
        raise ValueError('TwentySixNorth: grid not recognised: %s', grid)

    if not loadedAltMask:
        # Atlantic Mask not loaded
        raise ValueError('TwentySixNorth: Mask not loaded: ed: %s', grid)

        assert 0 


    zv = np.ma.array(nc.variables[keys[0]][..., latslice26Nnm, :]) # m/s
    atlmoc = np.array(np.zeros_like(zv[0, :, :, 0]))
    print('TwentySixNorth:', e3v_AMOC26N.shape, atlmoc.shape, zv.shape)

    for (z, la, lo), _ in np.ndenumerate(e3v_AMOC26N):
        if not alttmask_AMOC26N[la, lo]:
            continue
        if not tmask_AMOC26N[z, la, lo]:
            continue
        if np.ma.is_masked(zv[0, z, la, lo]):
            continue
        atlmoc[z, la] = atlmoc[z, la] - e1v_AMOC26N[la, lo] * e3v_AMOC26N[z, la, lo] * zv[0, z, la, lo] / 1.E06

    for z in range(e3v_AMOC26N.shape[0] -2, 1, -1): # add from the bottom up
        atlmoc[z, :] = atlmoc[z+1, :] + atlmoc[z, :]
    return atlmoc


def twentysixnorth025(nc,keys,**kwargs):
    """
    This function loads the AMOC array that is used for eORCA025

    nc: a netcdf openned as a dataset.
    keys: a list of keys to use in this function.

    """
    areafile = get_kwarg_file(kwargs, 'areafile')
    maskname = kwargs.get('maskname', 'tmask')
    grid = kwargs.get('grid', 'eORCA025')

    if not loadedArea:
        loadDataMask(areafile, maskname, grid)

    if grid != 'eORCA025':
        assert 0
    
    latslice26N = eORCA025_latslice26Nnm
    lonslice26N = eORCA025_lonslice26Nnm
    vo =  np.ma.array(nc.variables[keys[0]][..., latslice26N, lonslice26N]) # #vo in m/s
    thkcello = np.ma.array(nc.variables['thkcello'][..., latslice26N, lonslice26N]) # #thickness
    depths = np.ma.cumsum(thkcello, axis=1)

    depths = np.ma.masked_where(thkcello.mask + np.abs(depths)<500., depths) # masked above 500m depth.

    e1v = e1v_AMOC26N[:,None, :, :]
    flux = vo * depths * e1v_AMOC26N[:,None, :, :]/1.E06
    
    moc=np.ma.zeros_like(flux)
    np.ma.cumsum(flux[:,::-1], axis=1, out=moc ) # sum floor to surface
    return moc.max()

 
def AMOC26N(nc, keys, **kwargs):
    if kwargs.get('grid',None) == 'eORCA025':
        return twentysixnorth025(nc, keys, **kwargs)
    else:
        atlmoc = TwentySixNorth(nc, keys, **kwargs)
    return atlmoc.max()


def fov_sa(nc, keys, **kwargs):
    # Fov/Mov defined in Jackson 2023 as:
    # We also use diagnostics of the overturning component of the Atlantic freshwater transport (Fov).
    # This is calculated as "equation", where 
    # vbar is the zonal mean of the meridional velocity, 
    # sbar is the zonal mean of the salinity, 
    # and S0 is a reference salinity, 35PSU (Rahmstorf, 1996; Hawkins et al., 2011; Weaver et al., 2012).
    # This is calculated with monthly mean velocity and salinity fields, 
    # which ignores the impacts of the higher-frequency covariances of v and S;
    # however, previous studies have found these to be small (Mecking et al., 2016; Jackson and Wood, 2018a).
    # We use a reference salinity of 35 PSU, except in the case of CESM2, 
    # for which a reference salinity of 34.7 PSU is used, 
    # although the implied difference in transports from these different reference salinities
    # is again very small (Mecking et al., 2017)."
    # https://gmd.copernicus.org/articles/16/1975/2023/

    # Check Grid
    grid = kwargs.get('grid', 'eORCA1')

    # Reference salinity, S0, default is 35. 
    sal_ref = kwargs.get('sal_ref', 35.)

    # Load lats,. lons, cell thickness
    lats = nc.variables['nav_lat'][:]
    lons = nc.variables['nav_lon'][:]
    thkcello = nc.variables['thkcello'][:].squeeze()

    # mask lats and lons to south Atlantic (30S-34S)
    lats = np.ma.masked_outside(lats, -30., -34.)
    lons = np.ma.masked_outside(lons, -65., 20.)
    mask_2d =  lats.mask + lons.mask

    # Check to make sure Longitude boundaries are even here.  
    lons_bounds_0 = nc.variables['bounds_lon'][:].min(2)
    lons_bounds_1 = nc.variables['bounds_lon'][:].max(2)
    lons_bounds_0 = np.ma.masked_where(mask_2d, lons_bounds_0)
    lons_bounds_1 = np.ma.masked_where(mask_2d, lons_bounds_1)
    lons_diff = lons_bounds_1 - lons_bounds_0
    if lons_diff.min() != lons_diff.max():
        print('Can not assume that longitude grid is even')
        assert 0

    # Load Atlantic Mask
#    if not loadedAltMask_full:
#        altmaskfile = get_kwarg_file(kwargs, 'altmaskfile', default = 'bgcval2/data/basinlandmask_eORCA1.nc')
#        loadAtlanticMask_full(altmaskfile, maskname='tmaskatl', grid=grid)
    
    # Load and mask vo and vso (vo * salinity)
    vso =  np.ma.array(nc.variables['vso'][:]).squeeze() # #vso in PSU m/s
    vo =  np.ma.array(nc.variables['vo'][:]).squeeze() # #vso in PSU m/s
    
    # Calculate salinity and subtract reference salininty.
    sal0 = vso/vo - sal_ref

    # Calculate zonal cell length. 
    # Lon grid is evenly spaced, so only need one cell length per latitude.
    unique_lats = {la:True for la in np.ma.masked_where(mask_2d, lats).compressed()}
    zonal_distances = {la: myhaversine(0, la, 1., la) for la in unique_lats.keys()}
 
    # create 3d mask
    mask_3d = [mask_2d for _ in range(75)]
    mask_3d = np.stack(mask_3d, axis=0)

    # check sizes
    if vso.shape != mask_3d.shape:
        print('FOV: Shapes don\'t match')
        assert 0
    
    # Apply masks to 3d data
    vo = np.ma.masked_where(vso.mask + mask_3d + (vo == 0.), vo) # shape alignment?
    sal0 = np.ma.masked_where(vso.mask + mask_3d + (sal0 == 0.), sal0) # shape alignment?

#    from matplotlib import pyplot
#    for name, dat in zip(
#        ['vo', 'sal0', 'xarea', 'lats', 'lons', 'mask2d', 'mask3d', 'alttmask',],
#        [vo, sal0, xarea, lats, lons, mask_2d, mask_3d, alttmask]):
#        print('plotting', name)
#        if name in ['lats', 'lons', 'mask2d', 'alttmask']:
#            plot = dat
#        else:
#           plot = dat.mean(0)
#        pyplot.pcolormesh(plot)
#        pyplot.title(name)
#        pyplot.colorbar()
#        pyplot.savefig('images/'+name+'.png')
#        pyplot.close()

    # calculate cross sectional area by multiplying zonal cell length by cell depth thickness
    xarea = np.ma.masked_where(sal0.mask, thkcello)

    for (z, y, x), thk in maenumerate(xarea):
        la = lats[y, x]
        xarea[z, y, x] = thk * zonal_distances[la]
    #xarea_sum = xarea.sum(axis=(0,2))
    xarea_sum2 = xarea.sum(axis=2)

    # Vbar = SUM( vo(x)*thickcello(x)*dx)  / SUM(  thickcello*dx)

    vobar = (vo* xarea).sum(2) / xarea_sum2

    print('vobar', vobar.shape, vobar.min(), vobar.max())

    sobar = (sal0* xarea).sum(2) / xarea_sum2
    print('sobar', sobar.shape, sobar.min(), sobar.max())


    vsbar = (vobar*sobar * xarea_sum2).sum(1)
    print('vsbar', vsbar.shape, vsbar.min(), vsbar.max())

    total = vsbar.mean()

    # Take the zonal sum of the meridional velocity, the normalised salinity and the cross sectional area 
    #total =  vo * sal0 * xarea 

    # Calculate the cross sectional total, then divide by the total cross section area
    #total = total.sum(axis=(0, 2))/xarea_sum

    #print('total', {f:True for f in total.compressed()}.keys())
    #pyplot.pcolormesh(vo[0] * sal0[0] * xarea[0])
    #pyplot.title('total')
    #pyplot.colorbar()
    #pyplot.savefig('images/total.png')
    #pyplot.close()

    # Take the mean in the meridional area
    #total = total.mean()

    # Apply factors from paper.
    output = (-1./sal_ref) * total
    print(output)
    assert 0
    return output


def AEU(nc, keys, **kwargs):
    """
    Calculate the Atlantic Equatorial Undercurrent
    23 W, 5 south to 5 north
    """
    # Get cross section
    grid = kwargs.get('grid', 'eORCA025')
    gridfn = get_kwarg_file(kwargs, 'areafile')

    if not loaded_AEU:
        load_AEU_masks(gridfn, grid)

    if grid == 'eORCA025':
        AEU_LON=eORCA025_AEU_LON
        lat_slice = slice(eORCA025_AEU_LAT0,eORCA025_AEU_LAT1)
        AEU_LAT0=eORCA025_AEU_LAT0
        AEU_LAT1=eORCA025_AEU_LAT1
    else: assert 0
    
    max_depth = kwargs.get('max_depth', 500.)

    thkcello = nc.variables['thkcello'][:, :, lat_slice, AEU_LON]
    depth = np.ma.abs(np.cumsum(thkcello[:], axis=1))

    # calculate cross section area
    cross_section = e2u_aeu * thkcello

    # Multiply current speed by area for total flux
    uo = nc.variables[keys[0]][:,:,lat_slice, AEU_LON]

    # mask deep water and Southbound current (below zero)
    uo = np.ma.masked_where(uo.mask + (depth> max_depth) + (uo<0.), uo)

    flux = uo * cross_section
    flux = flux.sum() / 1.E06
    return flux

