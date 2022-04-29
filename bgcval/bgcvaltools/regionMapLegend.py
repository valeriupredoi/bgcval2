#!/usr/bin/ipython
#
# Copyright 2014, Plymouth Marine Laboratory
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
.. module:: regionMapLegend 
   :platform: Unix
   :synopsis: Tool to make a plot showing various regions.
.. moduleauthor:: Lee de Mora <ledm@pml.ac.uk>

"""

from paths import orcaGridfn, WOAFolder_annual
from netCDF4 import Dataset
import numpy as np
import UKESMpython as ukp
from matplotlib import pyplot
import cartopy.crs as ccrs
from cartopy import img_transform, feature as cfeature
from bgcvaltools.pftnames import getLongName
from bgcvaltools.makeMask import makeMask


def robinPlotCustom(lons,
                    lats,
                    data,
                    filename,
                    title,
                    regionList,
                    zrange=[-100, 100],
                    drawCbar=True,
                    cbarlabel='',
                    doLog=False,
                    dpi=100,
                    cmapname='default',
                    crude=0):
    ####
    # Based on robinplotSingle

    fig = pyplot.figure()
    fig.set_size_inches(10, 5)

    lons = np.array(lons)
    lats = np.array(lats)
    data = np.ma.array(data)

    rbmi = min([
        data.min(),
    ])
    rbma = max([
        data.max(),
    ])

    if rbmi * rbma > 0. and rbma / rbmi > 100.: doLog = True

    print(lons.shape, lats.shape, data.shape)
    lon0 = 0.  #lons.mean()
    if crude:
        ax = pyplot.subplot(
            111)  #,projection=ccrs.PlateCarree(central_longitude=lon0, ))
        im = pyplot.scatter(
            lats,
            lons,
            c=data,
            lw=0.,
            s=3,
            cmap='viridis',
            vmin=rbmi,
            vmax=rbma,
        )
        pyplot.axvline(80.)
        pyplot.axvline(83.5)

        pyplot.axvline(-100.)
        #pyplot.axvline(-105.)
        pyplot.axvline(-96.)

        #pyplot.colorbar(im)
        #title, zrange=[rbmi,rbma],lon0=lon0,drawCbar=False,cbarlabel=cbarlabel,doLog=doLog,cmap = cmapname)
    else:
        ax = pyplot.subplot(111,
                            projection=ccrs.PlateCarree(
                                central_longitude=lon0, ))
        fig, ax, im = ukp.makemapplot(fig,
                                      ax,
                                      lons,
                                      lats,
                                      data,
                                      title,
                                      zrange=[rbmi, rbma],
                                      lon0=lon0,
                                      drawCbar=False,
                                      cbarlabel=cbarlabel,
                                      doLog=doLog,
                                      cmap=cmapname)

    cmap = im.get_cmap()
    #for i in [0,10,100,1000,10000,1000000,100000]:
    #	print i, cmap(i), data.min(),data.max()
    leg = True
    if leg:
        # Shrink current axis's height by 10% on the bottom
        box = ax.get_position()
        ax.set_position(
            [box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])

        for i, r in enumerate(regionList):
            c = cmap((i) / (len(regionList) - 1.))
            print(i, r, c)
            pyplot.plot([], [], lw=8, color=c, label=getLongName(r))

        # Put a legend below current axis
        leg = ax.legend(loc='upper center',
                        bbox_to_anchor=(0.5, -0.05),
                        ncol=3,
                        prop={'size': 9})

        #leg = pyplot.legend(loc='lower center',ncol=3, )
        leg.draw_frame(False)
        leg.get_frame().set_alpha(0.)

    print("robinPlotSingle.py:\tSaving:", filename)
    pyplot.savefig(filename, dpi=dpi)
    pyplot.close()


def robinPlotTransects(lons,
                       lats,
                       data,
                       filename,
                       title,
                       legends=[],
                       zrange=[-100, 100],
                       drawCbar=True,
                       cbarlabel='',
                       doLog=False,
                       dpi=100,
                       cmapname='jet',
                       proj='Arctic'):
    ####
    # Based on robinplotSingle

    lons = np.array(lons)
    lats = np.array(lats)
    data = np.ma.array(data)

    rbmi = min([
        data.min(),
    ])
    rbma = max([
        data.max(),
    ])

    if rbmi * rbma > 0. and rbma / rbmi > 100.: doLog = True

    print(lons.shape, lats.shape, data.shape)
    lon0 = 0.  #lons.mean()
    if proj in ['Arctic', 'Antartic']:
        fig = pyplot.figure()
        print(lons.shape, lats.shape, data.shape, [rbmi, rbma])
        lat2d, lon2d = np.meshgrid(lats, lons)

        if proj in [
                'Arctic',
        ]:
            ax = pyplot.axes(projection=ccrs.Orthographic(
                -10, 70))  #central_longitude=lon0, ))
        if proj in [
                'Antartic',
        ]:
            ax = pyplot.axes(projection=ccrs.Orthographic(
                -10, -70))  #central_longitude=lon0, ))

        m = data.mask
        lon2d = np.ma.masked_where(m, lon2d).compressed()
        lat2d = np.ma.masked_where(m, lat2d).compressed()
        data = np.ma.masked_where(m, data).compressed()
        im = ax.scatter(
            lat2d,
            lon2d,
            c=data,
            lw=0.,
            s=8.,
            marker='s',
            cmap=cmapname,
            transform=ccrs.PlateCarree())  #vmin=rbmi,vmax=rbma,s=3,

        #im = ax.pcolormesh(lat2d,lon2d,data,cmap=cmapname,transform=ccrs.PlateCarree())#vmin=rbmi,vmax=rbma,s=3,

        ax.add_feature(cfeature.LAND, facecolor='0.85')
    if proj in ['Both']:
        fig = pyplot.figure()
        fig.set_size_inches(10, 5)
        print(lons.shape, lats.shape, data.shape, [rbmi, rbma])
        lat2d, lon2d = np.meshgrid(lats, lons)

        m = data.mask
        lon2d = np.ma.masked_where(m, lon2d).compressed()
        lat2d = np.ma.masked_where(m, lat2d).compressed()
        data = np.ma.masked_where(m, data).compressed()

        ax1 = pyplot.subplot(121, projection=ccrs.Orthographic(
            -10, 70))  #central_longitude=lon0, ))
        im1 = ax1.scatter(
            lat2d,
            lon2d,
            c=data,
            lw=0.,
            s=10.,
            marker='s',
            cmap=cmapname,
            transform=ccrs.PlateCarree())  #vmin=rbmi,vmax=rbma,s=3,
        ax1.add_feature(cfeature.LAND, facecolor='0.85')

        ax = pyplot.subplot(122, projection=ccrs.Orthographic(
            -10, -70))  #central_longitude=lon0, ))
        im = ax.scatter(
            lat2d,
            lon2d,
            c=data,
            lw=0.,
            s=10.,
            marker='s',
            cmap=cmapname,
            transform=ccrs.PlateCarree())  #vmin=rbmi,vmax=rbma,s=3,

        #im = ax.pcolormesh(lat2d,lon2d,data,cmap=cmapname,transform=ccrs.PlateCarree())#vmin=rbmi,vmax=rbma,s=3,

        ax.add_feature(cfeature.LAND, facecolor='0.85')

    if proj == 'robin':
        fig = pyplot.figure()
        fig.set_size_inches(10, 5)
        ax = pyplot.subplot(111,
                            projection=ccrs.PlateCarree(
                                central_longitude=lon0, ))
        fig, ax, im = ukp.makemapplot(fig,
                                      ax,
                                      lons,
                                      lats,
                                      data,
                                      title,
                                      zrange=[rbmi, rbma],
                                      lon0=lon0,
                                      drawCbar=False,
                                      cbarlabel=cbarlabel,
                                      doLog=doLog,
                                      cmap=cmapname)

    leg = True
    if leg:
        cmap = im.get_cmap()

        # Shrink current axis's height by 10% on the bottom
        box = ax.get_position()
        ax.set_position(
            [box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
        if proj == 'Both':
            box1 = ax1.get_position()
            ax1.set_position([
                box1.x0, box1.y0 + box1.height * 0.1, box1.width,
                box1.height * 0.9
            ])

        for i, r in enumerate(legends):
            c = cmap((i) / (len(legends) - 1.))
            print('making transect legend:', i, r, c)
            pyplot.plot([], [], color=c, lw=8, label=getLongName(r))

        # Put a legend below current axis

        if proj in [
                'robin',
        ]:
            leg = ax.legend(loc='upper center',
                            bbox_to_anchor=(0.5, -0.05),
                            ncol=4,
                            prop={'size': 10})
        elif proj in ['Both']:
            leg = ax.legend(loc='upper center',
                            bbox_to_anchor=(0.0, -0.05),
                            ncol=4,
                            prop={'size': 10})
        else:
            leg = ax.legend(loc='upper center',
                            bbox_to_anchor=(0.5, -0.05),
                            ncol=4,
                            prop={'size': 9})
        #leg = pyplot.legend(loc='lower center',ncol=3, )
        leg.draw_frame(False)
        leg.get_frame().set_alpha(0.)

    print("robinPlotTransects.py:\tSaving:", filename)
    pyplot.savefig(filename, dpi=dpi)
    pyplot.close()


def makeRegionMap(regionList):

    plotAll = 0  #True	# make plots for all regions
    imageFold = ukp.folder('images/maps')
    #####
    # Load data.
    nc = Dataset(orcaGridfn, 'r')
    bathy = nc.variables['mbathy'][:]
    xy = np.ma.masked_where(bathy == 0,
                            nc.variables['nav_lat'][:]).compressed()
    xx = np.ma.masked_where(bathy == 0,
                            nc.variables['nav_lon'][:]).compressed()
    nc.close()

    cbathy = np.ma.masked_where(bathy == 0, bathy).compressed()
    xt = np.ones_like(cbathy)
    xz = np.ones_like(cbathy)

    ####
    # Calculate masks, based on lat/lon.
    masks = {}
    for r in regionList:
        masks[r] = ~makeMask('', r, xt, xz, xy, xx, cbathy, debug=True)

    #####
    # Turn mask into one field.
    data = np.zeros_like(cbathy)
    for i, r in enumerate(regionList):
        data += (i + 1) * masks[r]
        if plotAll:
            fn = imageFold + 'Region_Legend_' + r + '.png'
            ukp.robinPlotSingle(
                xy,
                xx,
                masks[r],
                fn,
                r,
                drawCbar=True,
                cbarlabel='',
                doLog=False,
                dpi=100,
            )
    data = np.ma.masked_where(data == 0, data)

    #####
    # Send it to the plotting tool.
    colourmaps = [
        'default',
    ]  #'rainbow','jet','gist_earth','terrain','ocean','hsv','gist_rainbow','nipy_spectral',]
    for c in colourmaps:
        fn = imageFold + 'Region_Legend.png'
        robinPlotCustom(xy,
                        xx,
                        data,
                        fn,
                        '',
                        regionList,
                        drawCbar=False,
                        cbarlabel='',
                        doLog=False,
                        dpi=200,
                        cmapname=c)


def makeRegionMapNA(regionList):

    plotAll = 0  #True	# make plots for all regions
    imageFold = ukp.folder('images/maps')
    #####
    # Load data.
    nc = Dataset(orcaGridfn, 'r')
    bathy = nc.variables['mbathy'][:]
    xy = np.ma.masked_where(bathy == 0,
                            nc.variables['nav_lat'][:]).compressed()
    xx = np.ma.masked_where(bathy == 0,
                            nc.variables['nav_lon'][:]).compressed()
    nc.close()

    cbathy = np.ma.masked_where(bathy == 0, bathy).compressed()
    xt = np.ones_like(cbathy)
    xz = np.ones_like(cbathy)

    ####
    # Calculate masks, based on lat/lon.
    masks = {}
    for r in regionList:
        masks[r] = ~makeMask('', r, xt, xz, xy, xx, cbathy, debug=True)

    #####
    # Turn mask into one field.
    data = np.zeros_like(cbathy)
    for i, r in enumerate(regionList):
        data += (i + 1) * masks[r]
        if plotAll:
            fn = imageFold + 'Region_Legend_NA_' + r + '_robin.png'
            ukp.robinPlotSingle(
                xy,
                xx,
                masks[r],
                fn,
                r,
                drawCbar=True,
                cbarlabel='',
                doLog=False,
                dpi=100,
            )

    data = np.ma.masked_where(data == 0, data)

    #####
    # Send it to the plotting tool.
    colourmaps = [
        'default',
    ]  #'rainbow','jet','gist_earth','terrain','ocean','hsv','gist_rainbow','nipy_spectral',]
    for c in colourmaps:
        fn = imageFold + 'Region_Legend_NorthAtlantic.png'
        robinPlotCustom(xy,
                        xx,
                        data,
                        fn,
                        '',
                        regionList,
                        drawCbar=False,
                        cbarlabel='',
                        doLog=False,
                        dpi=200,
                        cmapname=c)


def makeRegionMapPierce():

    PierceRegions = [
        'Enderby',
        'Wilkes',
        'Ross',
        'Amundsen',
        'Weddel',
    ]
    plotAll = 1  #True	# make plots for all regions
    imageFold = ukp.folder('images/maps')
    #####
    # Load data.
    nc = Dataset(orcaGridfn, 'r')
    bathy = nc.variables['mbathy'][:]
    xy = np.ma.masked_where(bathy == 0,
                            nc.variables['nav_lat'][:]).compressed()
    xx = np.ma.masked_where(bathy == 0,
                            nc.variables['nav_lon'][:]).compressed()
    nc.close()

    cbathy = np.ma.masked_where(bathy == 0, bathy).compressed()
    xt = np.ones_like(cbathy)
    xz = np.ones_like(cbathy)

    ####
    # Calculate masks, based on lat/lon.
    masks = {}
    for i, r in enumerate(PierceRegions):
        masks[r] = ~makeMask('', r, xt, xz, xy, xx, cbathy, debug=True)

    #####
    # Turn mask into one field.
    data = np.zeros_like(cbathy)
    for i, r in enumerate(PierceRegions):
        print(i, r, ':', i)
        data += (i + 1) * masks[r]
        data = np.clip(data, 0, i + 1)
        if plotAll:
            fn = imageFold + 'Region_Legend_Pierce_' + r + '_robin.png'
            ukp.robinPlotSingle(
                xy,
                xx,
                masks[r],
                fn,
                r,
                drawCbar=True,
                cbarlabel='',
                doLog=False,
                dpi=100,
            )

            fn = imageFold + 'Region_Legend_Pierce_' + r + '_polar.png'
            td = np.ma.masked_where(masks[r] == 0., cbathy).compressed()
            tx = np.ma.masked_where(masks[r] == 0., xx).compressed()
            ty = np.ma.masked_where(masks[r] == 0., xy).compressed()
            #robinPlotTransects(ty, tx, td,fn,r, dpi=100,cmapname='jet',proj='Antartic')

    data = np.ma.masked_where(data == 0, data)

    #####
    # Send it to the plotting tool.
    colourmaps = [
        'default',
    ]  #'rainbow','jet','gist_earth','terrain','ocean','hsv','gist_rainbow','nipy_spectral',]
    for c in colourmaps:
        fn = imageFold + 'Pierce_SouthernOceanRegions.png'
        robinPlotCustom(xy,
                        xx,
                        data,
                        fn,
                        '',
                        PierceRegions,
                        drawCbar=True,
                        cbarlabel='',
                        doLog=False,
                        dpi=200,
                        cmapname=c)


def makeRegionMapSouthAtlantic():

    PierceRegions = [
        'Wilkes',
        'Weddel',
        'AtlanticSOcean',
        'SouthernOcean',
    ]
    plotAll = 0  # 1#True	# make plots for all regions
    imageFold = ukp.folder('images/maps')

    #####
    # Load data.
    nc = Dataset(orcaGridfn, 'r')
    bathy = nc.variables['mbathy'][:]
    xy = np.ma.masked_where(bathy == 0,
                            nc.variables['nav_lat'][:]).compressed()
    xx = np.ma.masked_where(bathy == 0,
                            nc.variables['nav_lon'][:]).compressed()
    nc.close()

    cbathy = np.ma.masked_where(bathy == 0, bathy).compressed()
    xt = np.ones_like(cbathy)
    xz = np.ones_like(cbathy)

    ####
    # Calculate masks, based on lat/lon.
    masks = {}
    for i, r in enumerate(PierceRegions):
        # masks[r] = ~makeMask('',r, xt,xz,xy,xx,cbathy,debug=True)
        masks[r] = makeMask('', r, xt, xz, xy, xx, cbathy, debug=True)

    #####
    # Turn mask into one field.
    data = np.zeros_like(cbathy)
    for i, r in enumerate(PierceRegions):
        i += 1
        #data += (i+1)* masks[r]
        for itr, d in enumerate(masks[r]):
            if data[itr] == 0 and not d:
                data[itr] += i
        #data = np.clip(data, 0, i+1)

        if plotAll:
            fn = imageFold + 'Region_Legend_SouthAtlantic_' + r + '_robin.png'
            ukp.robinPlotSingle(
                xy,
                xx,
                masks[r],
                fn,
                r,
                drawCbar=True,
                cbarlabel='',
                doLog=False,
                dpi=100,
            )

            fn = imageFold + 'Region_Legend_SouthAtlantic_' + r + '_polar.png'
            td = np.ma.masked_where(masks[r] == 0., cbathy).compressed()
            tx = np.ma.masked_where(masks[r] == 0., xx).compressed()
            ty = np.ma.masked_where(masks[r] == 0., xy).compressed()
            #robinPlotTransects(ty, tx, td,fn,r, dpi=100,cmapname='jet',proj='Antartic')

    data = np.ma.masked_where(data == 0, data)

    #####
    # Send it to the plotting tool.
    colourmaps = [
        'default', 'jet', 'viridis_r'
    ]  #'gist_earth','terrain','ocean','hsv','gist_rainbow','nipy_spectral','rainbow',]
    for c in colourmaps:
        fn = imageFold + 'SouthAtlantic_' + c + '.png'
        robinPlotCustom(xy,
                        xx,
                        data,
                        fn,
                        '',
                        PierceRegions,
                        drawCbar=True,
                        cbarlabel='',
                        doLog=False,
                        dpi=200,
                        cmapname=c)


def makeRegionMapYevgeny():

    regionList = [
        'IrmingerSea', 'YevgenyLabradorSea', 'YevgenyGreenlandIcelandicSeas'
    ]

    plotAll = 0  #True	# make plots for all regions
    imageFold = ukp.folder('images/maps')
    #####
    # Load data.
    nc = Dataset(orcaGridfn, 'r')
    bathy = nc.variables['mbathy'][:]
    xy = np.ma.masked_where(bathy == 0,
                            nc.variables['nav_lat'][:]).compressed()
    xx = np.ma.masked_where(bathy == 0,
                            nc.variables['nav_lon'][:]).compressed()
    nc.close()

    cbathy = np.ma.masked_where(bathy == 0, bathy).compressed()
    xt = np.ones_like(cbathy)
    xz = np.ones_like(cbathy)

    ####
    # Calculate masks, based on lat/lon.
    masks = {}
    masks['IrmingerSea'] = ~makeMask(
        '', 'NordicSea', xt, xz, xy, xx, cbathy, debug=True)

    #: (i=223,j=251) - (i=247,j=251) - (i=247,j=226) - (i=223,j=226)
    lab = np.zeros_like(bathy)
    lab[226 + 40:251 + 40, 223:247] = 1
    masks['YevgenyLabradorSea'] = np.ma.masked_where((bathy == 0),
                                                     lab).compressed()

    norw = np.zeros_like(bathy)
    norw[244 + 40:272 + 40, 249:282] = 1
    masks['YevgenyGreenlandIcelandicSeas'] = np.ma.masked_where(
        (bathy == 0), norw).compressed()

    #lab = np.zeros_like(bathy)
    #lab[226:251,223:247] = 1
    #masks['YevgenyLabradorSea'] = np.ma.masked_where((bathy==0) ,lab).compressed()

    #norw = np.zeros_like(bathy)
    #norw[244:272,249:282] = 1
    #masks['YevgenyGreenlandIcelandicSeas'] = np.ma.masked_where((bathy==0),norw).compressed()

    #####

    #####
    # Turn mask into one field.
    data = np.zeros_like(cbathy)
    for i, r in enumerate(regionList):
        print(i, r, masks[r].shape, data.shape)
        data += (i + 1) * masks[r]
        if plotAll:
            fn = imageFold + 'Region_Legend_NA_' + r + '.png'
            ukp.robinPlotSingle(
                xy,
                xx,
                masks[r],
                fn,
                r,
                drawCbar=True,
                cbarlabel='',
                doLog=False,
                dpi=100,
            )
    data = np.ma.masked_where(data == 0, data)

    #####
    # Send it to the plotting tool.
    colourmaps = [
        'default',
    ]  #'rainbow','jet','gist_earth','terrain','ocean','hsv','gist_rainbow','nipy_spectral',]
    for c in colourmaps:
        fn = imageFold + 'Region_Legend__YevgenyNorthAtlantic.png'
        robinPlotCustom(xy,
                        xx,
                        data,
                        fn,
                        '',
                        regionList,
                        drawCbar=False,
                        cbarlabel='',
                        doLog=False,
                        dpi=200,
                        cmapname=c)


def makeTransectsMap(proj='robin'):
    """
	Makes a plot of the transect lines.
	"""

    plotAll = 0  #True	# make plots for all regions
    imageFold = ukp.folder('images/maps/')  #Transects')

    #####
    # Load data.
    nc = Dataset(WOAFolder_annual + '/woa13_all_o00_01.nc', 'r')  #Oxy
    oxy = nc.variables['o_mn'][0, 0]
    lat = nc.variables['lat'][:]
    lon = nc.variables['lon'][:]
    nc.close()

    maps = np.zeros(oxy.shape)
    transects = [
        'ArcTransect',
        'AntTransect',
        'Transect',
        'CanRusTransect',
        'Equator',
        'PTransect',
        'SOTransect',
    ]
    fn = imageFold + 'Transects_legend_' + proj + '.png'

    for i, transect in enumerate(transects):
        i += 1
        single_map = np.zeros(oxy.shape)

        if transect in ['Transect', 'PTransect']:
            if transect == 'Transect': x = -28.
            if transect == 'PTransect': x = 200.
            k = ukp.getclosestlon(x, lon, debug=True)
            maps[:, k] = i
            single_map[:, k] = i

        if transect in ['SOTransect', 'Equator']:
            if transect == 'SOTransect': y = -60.
            if transect == 'Equator': y = 0.
            k = ukp.getclosestlat(y, lat, debug=True)
            maps[k, :] = i
            single_map[k, :] = i
        if transect in [
                'ArcTransect',
                'AntTransect',
                'CanRusTransect',
        ]:
            numpoints = 300
            if transect in [
                    'ArcTransect',
            ]:
                longi = 0.
                minlat = 50.
                maxlat = 90.
                transectcoords = [(minlat + a * (maxlat - minlat) / numpoints,
                                   longi)
                                  for a in np.arange(numpoints)]  # lat,lon

                longi = -165.
                minlat = 60.
                maxlat = 90.
                transectcoords.extend([
                    (minlat + a * (maxlat - minlat) / numpoints, longi)
                    for a in np.arange(numpoints)
                ])  # lat,lon

            if transect == 'AntTransect':
                longi = 0.
                minlat = -90.
                maxlat = -40.
                transectcoords = [(minlat + a * (maxlat - minlat) / numpoints,
                                   longi)
                                  for a in np.arange(numpoints)]  # lat,lon

            if transect == 'CanRusTransect':
                longi = 83.5
                minlat = 65.
                maxlat = 90.
                transectcoords = [(minlat + a * (maxlat - minlat) / numpoints,
                                   longi)
                                  for a in np.arange(numpoints)]  # lat,lon

                longi = -96.
                minlat = 60.
                maxlat = 90.
                transectcoords.extend([
                    (minlat + a * (maxlat - minlat) / numpoints, longi)
                    for a in np.arange(numpoints)
                ])  # lat,lon

            lon2d, lat2d = np.meshgrid(lon, lat)
            for (ilat, ilon) in sorted(transectcoords):
                la, lo = ukp.getOrcaIndexCC(
                    ilat,
                    ilon,
                    lat2d,
                    lon2d,
                    debug=True,
                )
                maps[la, lo] = i
                single_map[la, lo] = i
        #####
        # Make plot
        singles = False
        if singles:
            fn = ukp.folder(imageFold + 'Transects') + transect + '.png'
            single_map = np.clip(single_map, 0, i)
            single_map = np.ma.masked_where(single_map == 0, single_map)
            ukp.robinPlotSingle(
                lat,
                lon,
                single_map,
                fn,
                transect,
                drawCbar=False,
                cbarlabel='',
                doLog=False,
                dpi=100,
            )
    #maps = np.clip(maps,0,i)

    maps = np.ma.masked_where(maps == 0, maps)  #.compessed()
    robinPlotTransects(lat,
                       lon,
                       maps,
                       fn,
                       '',
                       legends=transects,
                       drawCbar=False,
                       cbarlabel='',
                       doLog=False,
                       dpi=200,
                       proj=proj)


def main():
    makeRegionMapSouthAtlantic()
    return

    makeRegionMapPierce()

    makeRegionMapYevgeny()

    makeRegionMapNA(['NordicSea', 'LabradorSea', 'NorwegianSea'])
    regionList = [  #'Global', 'ignoreInlandSeas',
        'SouthernOcean',
        'Remainder',
        'Equator10',
        'NorthernSubpolarAtlantic',
        'NorthernSubpolarPacific',
        'ArcticOcean',
    ]
    for proj in [
            'Both',
            'robin',
    ]:  #'Arctic', 'Antartic', ]:
        makeTransectsMap(proj=proj)
    makeRegionMap(regionList)


if __name__ == "__main__":
    main()
