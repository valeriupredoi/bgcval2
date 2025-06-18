#!/usr/bin/ipython
#
# Copyright 2024, Plymouth Marine Laboratory
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
.. module:: generic_map_legend 
   :platform: Unix
   :synopsis: Tool to make a plot showing a regions.
.. moduleauthor:: Lee de Mora <ledm@pml.ac.uk>
.. active:: No
"""

# implement correct import of params if module in use
# from ..Paths.paths import orcaGridfn, WOAFolder_annual
import matplotlib
matplotlib.use('Agg')
import os
from netCDF4 import Dataset
import numpy as np
from bgcval2.bgcvaltools import bv2tools as bvt
from bgcval2.bgcvaltools.pftnames import getLongName
from bgcval2.bgcvaltools.makeMask import makeMask
from matplotlib import pyplot
import cartopy
import cartopy.crs as ccrs
from cartopy import img_transform, feature as cfeature
from bgcval2._runtime_config import get_run_configuration
from bgcval2.Paths.paths import paths_setter




# Functions#
# Make a single plot for each region.
# one pane for global map centered on the middle of the region.
# One pane zoomed in on center of region
# One pane global map.


def plot_globe(ax):
    pyplot.sca(ax)     
       
    # if quick:
    ax.add_feature(cfeature.OCEAN, zorder=0)
    ax.add_feature(cfeature.LAND, zorder=0, edgecolor='black')
    # else:
        # nc = Dataset(bathy_fn, 'r')
        # lats = nc.variables['lat'][::binning]
        # lons = nc.variables['lon'][::binning]

        # data = nc.variables['elevation'][::binning, ::binning]
        # nc.close()
 
        # data = np.ma.masked_where(data>0., data)

        # pyplot.pcolormesh(
        #             lons, 
        #             lats,
        #             data,
        #             #transform=proj,
        #             transform=ccrs.PlateCarree(),
        #             cmap=cmap,
        #             vmin=vmin, vmax=vmax, 
        #             )
        # ax.coastlines()
        # ax.add_feature(cfeature.LAND, edgecolor='black', facecolor=land_color, linewidth=0.5, zorder=9)

    ax.set_global()
    ax.gridlines()
    return ax


def plot_platcarre(ax):
    pyplot.sca(ax)     
       
    # if quick:
    ax.add_feature(cfeature.OCEAN, zorder=0)
    # else:
        # nc = Dataset(bathy_fn, 'r')
        # lats = nc.variables['lat'][::binning]
        # lons = nc.variables['lon'][::binning]

        # data = nc.variables['elevation'][::binning, ::binning]
        # nc.close()
 
        # data = np.ma.masked_where(data>0., data)

        # pyplot.pcolormesh(
        #             lons, 
        #             lats,
        #             data,
        #             #transform=proj,
        #             transform=ccrs.PlateCarree(),
        #             cmap=cmap,
        #             vmin=vmin, vmax=vmax, 
        #             )
        # ax.coastlines()
        # ax.add_feature(cfeature.LAND, edgecolor='black', facecolor=land_color, linewidth=0.5, zorder=9)

    ax.set_global()
    ax.gridlines()
    return ax


def add_region(fig, ax, lons, lats, data):
    #im = ax.scatter(lons, lats, c=data)
    pyplot.sca(ax)
    #im = ax.contourf(lons, lats, data, zorder=1000)
    im = ax.pcolormesh(lons, lats, data, zorder=1, transform=ccrs.PlateCarree(),)

    #pyplot.colorbar()

    ax.add_feature(cfeature.LAND, zorder=10, edgecolor='black')
    return fig, ax, im


def make_figure(region):
    """
    Make a figure for this region.
    """
    fig_fn = bvt.folder('images/regions')+region+'.png'

    paths_dict, config_user = get_run_configuration("defaults")
    # filter paths dict into an object that's usable below
    paths = paths_setter(paths_dict)   
    orcaGridfn = paths.orcaGridfn

    nc = Dataset(orcaGridfn, 'r')
    dat = nc.variables['mbathy'][:].squeeze()
    lats = nc.variables['nav_lat'][:].squeeze()
    lons = nc.variables['nav_lon'][:].squeeze()
    lons = bvt.makeLonSafeArr(lons)
    nc.close()

    old_mask = np.ma.masked_where(dat.mask + dat ==0, dat).mask

    xd = np.ma.masked_where(old_mask, dat).flatten()
    xt = np.ones_like(xd)
    xz = xt
    xy = np.ma.masked_where(old_mask, lats).flatten()
    xx = np.ma.masked_where(old_mask, lons).flatten()
    old_mask_flat = old_mask.flatten()

    region_mask = makeMask('bathy', region, xt, xz, xy, xx, xd, debug=True)
    #print('done makeMask')
    #assert 0
    new_dat = np.ma.masked_where(region_mask + old_mask_flat, xd)
    new_lon = np.ma.masked_where(region_mask+ old_mask_flat, xx)
    new_lat = np.ma.masked_where(region_mask+ old_mask_flat, xy)

    new_dat = new_dat.reshape(dat.shape)
    #new_lat = lats # new_lat.reshape(lats.shape)
    #new_lon = lons # new_lon.reshape(lons.shape)

    fig = pyplot.figure()
    fig.set_size_inches(12, 8)
    widths = [1, 1, 1]
    heights = [1, 1.75]       
    spec2 = matplotlib.gridspec.GridSpec(
        ncols=len(widths), 
        nrows=len(heights), 
        figure=fig, 
        width_ratios=widths,
        height_ratios=heights,
        hspace=0.30,
        wspace=0.30,)

    print('\''+region+'\': {\'lats\'', new_lat.mean(), 'lon:', new_lon.mean()) #data:', new_dat.min(), new_dat.max())

    ortho_pro=ccrs.Orthographic(new_lon.mean(), new_lat.mean(),)    
    ax_globe = fig.add_subplot(spec2[0, 0], projection=ortho_pro)
    ax_globe = plot_globe(ax_globe)
    fig, ax_globe, im = add_region(fig, ax_globe, lons, lats, new_dat)

    ortho_pro=ccrs.Orthographic(new_lon.mean()+120., new_lat.mean(),)    
    ax_globe1 = fig.add_subplot(spec2[0, 1], projection=ortho_pro)
    ax_globe1 = plot_globe(ax_globe1)
    fig, ax_globe1, im1 = add_region(fig, ax_globe1, lons, lats, new_dat)


    ortho_pro=ccrs.Orthographic(new_lon.mean()-120., new_lat.mean(),)    
    ax_globe2 = fig.add_subplot(spec2[0, 2], projection=ortho_pro)
    ax_globe2 = plot_globe(ax_globe2)
    fig, ax_globe2, im2 = add_region(fig, ax_globe2, lons, lats, new_dat)
    
    pc_proj=cartopy.crs.PlateCarree(central_longitude=new_lon.mean())
    ax_pc = fig.add_subplot(spec2[1, :], projection=pc_proj)    
    ax_pc = plot_platcarre(ax_pc)
    fig, ax_pc, im3 = add_region(fig, ax_pc, lons, lats, new_dat)
    #cbar = pyplot.colorbar(ax=ax_pc, cax=im3)

    pyplot.suptitle(region+': '+getLongName(region))
    print('saving:', fig_fn)
    pyplot.savefig(fig_fn,dpi=300.)
    pyplot.savefig(fig_fn.replace('.png', '_trans.png'), transparent=True)
    pyplot.close()


def main():
#    paths_dict, config_user = get_run_configuration("defaults")

    regions = [
#                'LIseas',
#                'LIGINseas',
#                'GLINseas',
#               'Ascension',
#               'ITCZ',
#               'TristandaCunha',
#               'Pitcairn',
#               'Cornwall',
#                'SubtropicNorthAtlantic',
#                'SPNA',
#                'STNA',        
                'SouthernOcean',
                'subpolar',
                'NorthEastAtlantic',
#               'ArcticOcean',
#               'Equator10',
#               'NorthPacificOcean',
                'SouthAtlanticOcean',
#               'SouthPacificOcean',
#                'NorthAtlanticOcean',
#               'SouthAtlanticOcean',
#                'GINseas',
#                'LabradorSea',
#                'IrmingerSea',
                'EquatorialAtlanticOcean',
#              'Global',
#              'ignoreInlandSeas',
    ]
    for region in regions[:]:
        make_figure(region) 

if __name__ == "__main__":
    main()
