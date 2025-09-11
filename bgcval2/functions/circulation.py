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


eORCA1_davis_LON=300
eORCA1_davis_LAT0=229
eORCA1_davis_LAT1=250


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
eORCA1_latslice40N = slice(245,246)
eORCA1_latslice55N = slice(272,273)

eORCA1_lonslice_GS = slice(228,229)


eORCA1_latslice32S = slice(137,138)

eORCA025_latslice26Nnm = slice(794,795)
# eORCA025_lonslice26Nnm = slice(825,1105) # Florida
eORCA025_lonslice26Nnm = slice(841, 1105) # Bahamas (removes Florida Straight)

umask_drake = 0
e2u_drake = 0
e2u_aeu = 0
umask_aeu = 0

e3v_AMOC26N = 0
e3v_AMOC40N = 0
e3v_AMOC55N = 0

e1v_AMOC26N = 0
e1v_AMOC40N = 0
e1v_AMOC55N = 0
e1v_davis = 0
tmask_AMOC26N = 0
tmask_AMOC40N = 0
tmask_AMOC55N = 0

alttmask_AMOC26N = 0
alttmask_AMOC40N = 0
alttmask_AMOC55N = 0

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
    global e3v_AMOC40N
    global e3v_AMOC55N
    global e1v_davis    
    global e1v_AMOC26N    
    global e1v_AMOC40N    
    global e1v_AMOC55N
    global tmask_AMOC26N
    global tmask_AMOC40N
    global tmask_AMOC55N
    
    global loadedArea
#    global loadedAltMask

    if grid == 'eORCA1':
        LON = eORCA1_drake_LON
        LAT0 = eORCA1_drake_LAT0
        LAT1 = eORCA1_drake_LAT1
        latslice26Nnm = eORCA1_latslice26Nnm
        latslice40N = eORCA1_latslice40N
        latslice55N = eORCA1_latslice55N
        # no lon slice for eORCA1

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

        e3v_AMOC40N = nc.variables['e3v_0'][..., latslice40N, lonslice40N].squeeze()   # z level height 3D 
        e1v_AMOC40N = nc.variables['e1v'][..., latslice40N, lonslice40N]     #
        tmask_AMOC40N = nc.variables['tmask'][..., latslice40N, lonslice40N]

        e3v_AMOC55N = nc.variables['e3v_0'][..., latslice55Nnm, lonslice55N].squeeze()   # z level height 3D 
        e1v_AMOC55N = nc.variables['e1v'][..., latslice55Nnm, lonslice55N]     #
        tmask_AMOC55N = nc.variables['tmask'][..., latslice55Nnm, lonslice55N]

        e1v_davis = nc.variables['e1v'][eORCA1_davis_LON, eORCA1_davis_LAT0:eORCA1_davis_LAT1]  

    else:
        e3v_AMOC26N = nc.variables['e3v'][..., latslice26Nnm, :]   # z level height 3D
        e1v_AMOC26N = nc.variables['e1v'][..., latslice26Nnm, :]     #
        tmask_AMOC26N = nc.variables['tmask'][..., latslice26Nnm, :]

        e3v_AMOC40N = nc.variables['e3v'][..., latslice40N, :]   # z level height 3D
        e1v_AMOC40N = nc.variables['e1v'][..., latslice40N, :]     #
        tmask_AMOC40N = nc.variables['tmask'][..., latslice40N, :]        

        e3v_AMOC55N = nc.variables['e3v'][..., latslice55N, :]   # z level height 3D
        e1v_AMOC55N = nc.variables['e1v'][..., latslice55N, :]     #
        tmask_AMOC55N = nc.variables['tmask'][..., latslice55N, :]

        e1v_davis = nc.variables['e1v'][eORCA1_davis_LON, eORCA1_davis_LAT0:eORCA1_davis_LAT1]

    #print('e3v_AMOC26N: loaded')#e3v_AMOC26N, latslice26Nnm, e3v_AMOC26N.shape)
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
    global alttmask_AMOC40N
    global alttmask_AMOC55N

    global loadedAltMask
    if grid == 'eORCA1':
        latslice26Nnm = eORCA1_latslice26Nnm
        latslice40N = eORCA1_latslice40N
        latslice55N = eORCA1_latslice55N

    else:
        raise ValueError("Grid not recognised in this calculation: %s", grid)
    nc = dataset(altmaskfile, 'r')        
    alttmask_AMOC26N = nc.variables[maskname][latslice26Nnm, :]
    alttmask_AMOC40N = nc.variables[maskname][latslice40N, :]
    alttmask_AMOC55N = nc.variables[maskname][latslice55N, :]
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


def davisstraightsalt(nc, keys, **kwargs):
    """
    This function calculates the salt flux through the davis straight in eORCA1. 
    
    nc: a netcdf openned as a dataset.
    keys: a list of keys to use in this function.
    
    """
    areafile = get_kwarg_file(kwargs, 'areafile')
    maskname = kwargs.get('maskname', 'tmask')
    grid = kwargs.get('grid', 'eORCA1')

    if not loadedArea:
        loadDataMask(areafile, maskname, grid)

    if grid == 'eORCA1':    
        LON = eORCA1_davis_LON
        LAT0 = eORCA1_davis_LAT0
        LAT1 = eORCA1_davis_LAT1
    else:
        assert 0

    print('Davis straight:', grid, 'LON', LON, 'LAT0',LAT0, 'LAT1', LAT1)

    vso = nc.variables['vso'][0, :, LAT0:LAT1, LON]
    thkcello = nc.variables['thkcello'][0, :, LAT0:LAT1, LON]
    e1v_4d = np.broadcast_to(e1v_davis[np.newaxis, :, :], vso.shape[:])

    vso = np.ma.masked_where(vso==0., vso)

    if vso.shape == thkcello.shape == e1v_4d.shape :
        pass
    else:
        print('Shapes do not match', vso.shape, thkcello.shape, e1v_4d.shape)
        assert 0

    davis = np.ma.sum(vso * e1v_4d * thkcello) * 1.e-6  # PSU m s-1 * m * m or PSU Sv

    return davis


def TwentySixNorth(nc, keys, lat='26N', return_max_depth=False, **kwargs):
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

    altmaskfile = get_kwarg_file(kwargs, 'altmaskfile', default = 'bgcval2/data/basinlandmask_eORCA1.nc')
    if not loadedAltMask:
         loadAtlanticMask(altmaskfile, maskname='tmaskatl', grid=grid)

    if grid == 'eORCA1':
        if lat == '26N':
            latslice = eORCA1_latslice26Nnm
            e1v_AMOC = e1v_AMOC26N
            alttmask_AMOC = alttmask_AMOC26N[:]
            tmask_AMOC = tmask_AMOC26N
            e3v_AMOC = e3v_AMOC26N

        elif lat == '40N':
            latslice = eORCA1_latslice40N
            e1v_AMOC = e1v_AMOC40N
            alttmask_AMOC = alttmask_AMOC40N[:]
            tmask_AMOC = tmask_AMOC40N
            e3v_AMOC = e3v_AMOC40N

        elif lat == '55N':
            latslice = eORCA1_latslice55N
            e1v_AMOC = e1v_AMOC55N
            alttmask_AMOC = alttmask_AMOC55N[:]
            tmask_AMOC = tmask_AMOC55N
            e3v_AMOC = e3v_AMOC55N

        else:
            raise ValueError('Region not recognised', lat)

    elif grid == 'eORCA025':
        assert 0 
        latslice = eORCA025_latslice26Nnm

    else:
        # grid not recognised.
        raise ValueError('TwentySixNorth: grid not recognised: %s', grid)

    if not loadedAltMask:
        # Atlantic Mask not loaded
        raise ValueError('TwentySixNorth: Mask not loaded: ed: %s', grid)

    zv = np.ma.array(nc.variables[keys[0]][..., latslice, :]) # m/s
    zv = np.ma.masked_where(zv.mask + (zv == 0.), zv)

    atlmoc = np.array(np.zeros_like(zv[0, :, :, 0]))

    if 'thkcello' in nc.variables.keys():
        thkcello = nc.variables['thkcello'][0, :, latslice, :]
        thkcello = np.ma.masked_where(thkcello.mask + zv[0].mask, thkcello)
    else:
        thkcello = e3v_AMOC[:]

    depths = np.ma.abs(np.cumsum(thkcello, axis=0))

    for (z, la, lo), _ in np.ndenumerate(thkcello):
        if not alttmask_AMOC[la, lo]:
            continue
        if not tmask_AMOC[z, la, lo]:
            continue
        if np.ma.is_masked(zv[0, z, la, lo]):
            continue
        atlmoc[z, la] = atlmoc[z, la] - e1v_AMOC[la, lo] * thkcello[z, la, lo] * zv[0, z, la, lo] / 1.E06

    for z in range(thkcello.shape[0] -2, 1, -1): # add from the bottom up
        atlmoc[z, :] = atlmoc[z+1, :] + atlmoc[z, :]
    if return_max_depth:
        max_depth_index = np.argmax(atlmoc, keepdims=True)
        max_depth2d = depths[max_depth_index[0], :] # extract the depth layer
        max_depth2d = np.ma.masked_where( max_depth2d.mask + (alttmask_AMOC==0.), max_depth2d) # mask the non-atlantic

#       print('max depth ', atlmoc.max(), atlmoc.shape, max_depth_index, atlmoc[max_depth_index], depths.shape)
#       print('max_depth2d', max_depth2d.min(), max_depth2d.mean(), max_depth2d.max())
        return max_depth2d.mean()
    else:
        # AMOC calculation
        return atlmoc


def AMOC_depth_2(nc,keys,**kwargs):
    """
    Calculates the depth of the amoc maxima.
    
    nc: a netcdf openned as a dataset.
    keys: a list of keys to use in this function.
    """
    return TwentySixNorth(nc, keys, **kwargs, return_max_depth=True)

def AMOC_depth(nc,keys,**kwargs):
    """
    Calculates the depth of the Gulf Stream AMOC border.
    
    nc: a netcdf openned as a dataset.
    keys: a list of keys to use in this function.

    """

    # this is wrong agin
    # change strategy:
    # use V.
    # call it GF_depth
    # extract V window at 26N
    # Find longitude with highest Northern current
    # Extract V and depth in that spot.

    depthw = nc.variables['depthw'][:]
    deptht = nc.variables['deptht'][:].squeeze()
    depthw_bounds = nc.variables['depthw_bounds'][:]

    grid = kwargs.get('grid', 'eORCA1')
    if grid == 'eORCA1':
        eORCA1_lat_26N = 228

    # Load Stream function at 26N
    sf26 = nc.variables[keys[0]][0, :, eORCA1_lat_26N, :].squeeze()
    lats = nc.variables['nav_lat'][ eORCA1_lat_26N, :]

    print('sf36:', sf26, sf26.shape, 'lats', lats)
    non_cumsum = [ ]
    for i, (z, sf) in enumerate(zip(deptht, sf26)):
        if i == 0:
            diff = 0
            non_cumsum.append(sf)
        else:
            diff = sf - non_cumsum[i-1]
            non_cumsum.append(diff)
        print(i, z, sf, diff)
    deptht = np.ma.masked_outside(deptht, 250., 3000)
    #sf26 = np.ma.masked_where((sf26==0.) + sf26.mask + deptht.mask, sf26)
    
    # indexes where streamfunction goes positive to negative:
    changes = np.where(sf26[:-1]*sf26[1:]<0)[0]
    changes = changes[sf26[changes]<0] 

    print(changes)

    assert 0
    return TwentySixNorth(nc, keys, **kwargs, return_max_depth=True)


def gulfstream_depth(nc, keys, **kwargs):
    """
    This function calculates the depth of the bottom of 
    Gulf Stream at 26N.
    This is the sum of Northbound current between

    nc: a netcdf openned as a dataset.
    keys: a list of keys to use in this function.

    """
    areafile = get_kwarg_file(kwargs, 'areafile')
    maskname = kwargs.get('maskname', 'tmask')
    grid = kwargs.get('grid', 'eORCA1')
    maxdepth = np.abs(kwargs.get('maxdepth', 2000. ))

    if not loadedArea:
        loadDataMask(areafile, maskname, grid)

    if grid == 'eORCA1':
        latslice26Nnm = eORCA1_latslice26Nnm
        #data=[-80.5011659 , -79.50119298, -78.50121829, -77.50124181,
        #           -76.50126349, -75.50128329, -74.50130118, -73.50131712,
        #           -72.50133107, -71.50134301, -70.50135293, -69.50136079,
        #           -68.50136658],
        lonslice_70W = slice(211, 217)

        altmaskfile = get_kwarg_file(kwargs, 'altmaskfile', default = 'bgcval2/data/basinlandmask_eORCA1.nc')
    elif grid == 'eORCA025':
        latslice26Nnm = eORCA025_latslice26Nnm
    else:
        # grid not recognised.
        raise ValueError('gulfstream: grid not recognised: %s', grid)

    lats = nc.variables['nav_lat'][latslice26Nnm, lonslice_70W]
    lons = nc.variables['nav_lon'][latslice26Nnm, lonslice_70W]

    print(lats, lons)

    vo = nc.variables[keys[0]][0, :, latslice26Nnm, lonslice_70W].squeeze() # m/s
    thickness = nc.variables['thkcello'][0,:,latslice26Nnm, lonslice_70W].squeeze()
    depth = np.abs(np.cumsum(thickness, axis=0))# depth array

    vo = np.ma.masked_where(vo.mask + (vo == 0.), vo)


    vo_max_index = np.argmax(vo,  keepdims=True)
    index = np.unravel_index(vo.argmax(), vo.shape)
    print('GS depth:', vo.shape, depth.shape, vo_max_index, vo.max(), vo_max_index.shape, depth.max(), index)
    index = np.unravel_index(vo.argmax(), vo.shape)

    gs_max_depth = depth[index]

    print('gs_max_depth:', gs_max_depth, vo[index])
    z1, z2 = 1., 1.
    v1, v2 = 1., 1.
    
    col_depth = depth[:, index[1]]
    col_vo = vo[:, index[1]] 

    found = False
    for i, z in enumerate(col_depth):
        if found:
            continue
        if np.abs(z) < np.abs(gs_max_depth): 
            continue

        v1 = col_vo[i-1]
        v2 = col_vo[i]

        z1 = col_depth[i-1]
        z2 = col_depth[i]
        print('column:', i, (z1, v1), (z2, v2))
        if v2 * v1 <= 0.:
            found = True
            print('Found inflection point:', i, (z1, v1), (z2, v2))

    if not found: 
        print('Failure!, is there no gulf stream here?')
        assert 0 

    m = (z2 - z1)/(v2 - v1)
    c = z2 - m*v2
    zdepth = c #-1.*c/m
    print('Zero V depth', zdepth, ('v = ',m, '* z +', c))
    return zdepth
 

def gulfstream(nc, keys, **kwargs):
    """
    This function calculates the Gulf Stream at 26N.
    This is the sum of Northbound current between 

    nc: a netcdf openned as a dataset.
    keys: a list of keys to use in this function.

    """
    areafile = get_kwarg_file(kwargs, 'areafile')
    maskname = kwargs.get('maskname', 'tmask')
    grid = kwargs.get('grid', 'eORCA1')
    maxdepth = np.abs(kwargs.get('maxdepth', 2000. ))
 
    if not loadedArea:
        loadDataMask(areafile, maskname, grid)

    if grid == 'eORCA1':
        latslice26Nnm = eORCA1_latslice26Nnm
        #data=[-80.5011659 , -79.50119298, -78.50121829, -77.50124181,
        #           -76.50126349, -75.50128329, -74.50130118, -73.50131712,
        #           -72.50133107, -71.50134301, -70.50135293, -69.50136079,
        #           -68.50136658],
        lonslice_70W = slice(207, 220) 

        altmaskfile = get_kwarg_file(kwargs, 'altmaskfile', default = 'bgcval2/data/basinlandmask_eORCA1.nc')
        if not loadedAltMask:
             loadAtlanticMask(altmaskfile, maskname='tmaskatl', grid=grid)
    elif grid == 'eORCA025':
        latslice26Nnm = eORCA025_latslice26Nnm
    else:
        # grid not recognised.
        raise ValueError('gulfstream: grid not recognised: %s', grid)

    if not loadedAltMask:
        # Atlantic Mask not loaded
        raise ValueError('gulfstream: Mask not loaded: ed: %s', grid)
        assert 0

    lats = nc.variables['nav_lat'][latslice26Nnm, lonslice_70W]
    lons = nc.variables['nav_lon'][latslice26Nnm, lonslice_70W]
    vo = np.ma.array(nc.variables[keys[0]][0, :, latslice26Nnm, lonslice_70W]) # m/s
    vo = np.ma.masked_where(vo.mask + (vo <= 0.), vo) 

    thickness = nc.variables['thkcello'][0,:,latslice26Nnm, lonslice_70W] 
    depth = np.abs(np.cumsum(thickness, axis=0))# depth array
    #print(vo.shape, thickness.shape, e1v_AMOC26N.shape)
    gs = 0.
    for (z, la, lo), v in np.ndenumerate(vo):
        if depth[z, la,lo] > maxdepth:
            continue
        if v <= 0:
            continue
        #print((z, la, lo),'depth:', depth[z, la,lo], (lats[la, lo],'N', lons[la, lo], 'E'),  'v:', v, 'thickness:', thickness[z, la, lo], 'width:', e1v_AMOC26N[la, lo])
        gs += v * thickness[z, la, lo] * e1v_AMOC26N[la, lo] / 1.E06

    print('Gulf Stream:', gs) # expecting a value of 32Sv ish.
    # https://www.sciencedirect.com/science/article/pii/S0079661114001694
    return gs




def twentysixnorth025(nc, keys, **kwargs):
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
        atlmoc = TwentySixNorth(nc, keys, lat='26N', **kwargs)
    return atlmoc.max()

def AMOC40N(nc, keys, **kwargs):
    if kwargs.get('grid',None) == 'eORCA025':
        return twentysixnorth025(nc, keys, **kwargs)
    else:
        atlmoc = TwentySixNorth(nc, keys,lat='40N', **kwargs)
    return atlmoc.max()

def AMOC55N(nc, keys, **kwargs):
    if kwargs.get('grid',None) == 'eORCA025':
        return twentysixnorth025(nc, keys, **kwargs)
    else:
        atlmoc = TwentySixNorth(nc, keys, lat='55N', **kwargs)
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
    #vso =  np.ma.array(nc.variables['vso'][:]).squeeze() # #vso in PSU m/s
    vo =  np.ma.array(nc.variables['vo'][:]).squeeze() # #vso in PSU m/s
   
    #print(nc.variables['vso'], '\n', nc.variables['vo']) 
    # Calculate salinity and subtract reference salininty.
    fn2 = nc.filename.replace('grid_V', 'grid_T').replace('grid-V', 'grid-T')
    nc2 = dataset(fn2, 'r')
    sal0 = nc2.variables['so'][:].squeeze()
    nc2.close()

    # Calculate zonal cell length. 
    # Lon grid is evenly spaced, so only need one cell length per latitude.
    unique_lats = {la:True for la in np.ma.masked_where(mask_2d, lats).compressed()}
    zonal_distances = {la: myhaversine(0, la, 1., la) for la in unique_lats.keys()}
 
    # create 3d mask
    mask_3d = [mask_2d for _ in range(75)]
    mask_3d = np.stack(mask_3d, axis=0)

    # check sizes
    if sal0.shape != mask_3d.shape:
        print('FOV: Shapes don\'t match')
        assert 0
    
    # Apply masks to 3d data
    vo = np.ma.masked_where(vo.mask + mask_3d + (vo == 0.), vo) # shape alignment?
    sal0 = np.ma.masked_where(sal0.mask + mask_3d + (sal0 == 0.), sal0) # shape alignment?

    #print('a sal0', sal0.shape, sal0.min(), sal0.max())

    sal0 = sal0 - sal_ref
    #print('b sal0', sal0.shape, sal0.min(), sal0.max())

    #print('vo:', vo.shape, vo.min(), vo.max())
    #print('c sal0', sal0.shape, sal0.min(), sal0.max())
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
    xarea_sum2 = xarea.sum(axis=2)
    #print('xarea_sum2:', xarea_sum2.shape, xarea_sum2.min(), xarea_sum2.max())


    # Calculate vobar and sobar. 
    # Vbar = SUM( vo(x)*thickcello(x)*dx)  / SUM(  thickcello*dx)
    vobar = (vo* xarea).sum(2) / xarea_sum2
    #obar = vo.mean(2)
    #print('vobar', vobar.shape, vobar.min(), vobar.max())

    #sobar = sal0.mean(2)
    sobar = (sal0* xarea).sum(2) / xarea_sum2
    #print('sobar', sobar.shape, sobar.min(), sobar.max())

    vsbar = (vobar*sobar * xarea_sum2).sum(0)
    #vsbar = (vobar*sobar * xarea.mean(2)).sum(0)

    #print('vsbar', vsbar.shape, vsbar.min(), vsbar.max(), vsbar.compressed())

    total = vsbar.mean()

    # Take the zonal sum of the meridional velocity, the normalised salinity and the cross sectional area 
    #total =  vo * sal0 * xarea    # m/s * PSU * m2

    # Calculate the cross sectional total, then divide by the total cross section area
    #total = total.sum(axis=(0, 2))/xarea_sum  # PSU m3/s /m2

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
    output = (-1./sal_ref) * total   # 1/PSU * PSU m/s
    output = output * 1e-6 # Convert m3 s-1 to Sv. 
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

