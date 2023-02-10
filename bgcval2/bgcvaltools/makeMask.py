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
.. module:: makeMask
   :platform: Unix
   :synopsis: A function that produces a mask after a region.
.. moduleauthor:: Lee de Mora <ledm@pml.ac.uk>

"""
import numpy as np
from calendar import month_name
from shelve import open as shOpen
import os
from bgcval2.Paths.paths import paths
from netCDF4 import Dataset
import bgcval2.bgcvaltools.bv2tools as bvt


def makeMask(name, newSlice, xt, xz, xy, xx, xd, debug=False):
    """
	:param name: The name of the data. (useful for debugging)
	:param newSlice: The name of the regional cut (or slice)
	:param xt: A one-dimensional array of the dataset times.
	:param xz: A one-dimensional array of the dataset depths.
	:param xy: A one-dimensional array of the dataset latitudes.
	:param xx: A one-dimensional array of the dataset longitudes.
	:param xd: A one-dimensional array of the data.
	
	This function produces a mask to hides all points that are not in the requested region.
	
	Note that xt,xz,xy,xx,xd should all be the same shape and size. 
	
	This functional can call itself, if two regional masks are needed.
	
	Please add your own regions, at the bottom of the list, if needed.
	"""
    if debug:
        print("makeMask:\tmakeMask:\tinitialise:\t", name, '\t', newSlice)

    nmask = np.zeros(len(xd))  # nothing masked

    #####
    # Equavalent to keeping the same mask.
    if newSlice in ['All', 'Global', 'regionless', 'layerless']:

        for a in [xt, xz, xy, xx, xd]:
            try:
                nmask += a.mask
            except:
                pass
        return nmask

    #####
    # Simple masks
    if newSlice == 'nonZero': return np.ma.masked_where(xd == 0., nmask).mask
    if newSlice == 'aboveZero': return np.ma.masked_where(xd <= 0., nmask).mask
    if newSlice == 'belowZero': return np.ma.masked_where(xd >= 0., nmask).mask

    #####
    # Simple Regional masks
    if newSlice == 'NorthHemisphere':
        return np.ma.masked_where(xy < 0., nmask).mask
    if newSlice == 'SouthHemisphere':
        return np.ma.masked_where(xy > 0., nmask).mask

    if newSlice == 'Tropics':
        return np.ma.masked_where(abs(xy) > 23., nmask).mask
    if newSlice == 'Equatorial':
        return np.ma.masked_where(abs(xy) > 7., nmask).mask
    if newSlice == 'Temperate':
        return np.ma.masked_where((abs(xy) < 23.) + (abs(xy) > 60.),
                                  nmask).mask
    if newSlice == 'NorthTropics':
        return np.ma.masked_where((xy > 23.) + (xy < 7.), nmask).mask
    if newSlice == 'SouthTropics':
        return np.ma.masked_where((xy < -23.) + (xy > -7.), nmask).mask
    if newSlice == 'NorthTemperate':
        return np.ma.masked_where((xy < 23.) + (xy > 60.), nmask).mask
    if newSlice == 'SouthTemperate':
        return np.ma.masked_where((xy > -23.) + (xy < -60.), nmask).mask

    if newSlice == 'AtlanticTransect':
        return np.ma.masked_where((xx > -26.) + (xx < -30.), nmask).mask
    if newSlice == 'PacificTransect':
        return np.ma.masked_where((xx > -139.) + (xx < -143.), nmask).mask
    if newSlice == '10N':
        return np.ma.masked_where((xy > 12.) + (xy < 8.), nmask).mask
    if newSlice == '10S':
        return np.ma.masked_where((xy > -8.) + (xy < -12.), nmask).mask
    if newSlice == 'SouthernTransect':
        return np.ma.masked_where((xy > -55.) + (xy < -59.), nmask).mask

    if newSlice == 'Arctic':
        return np.ma.masked_where(abs(xy) < 60., nmask).mask
    if newSlice == 'Antarctic':
        return np.ma.masked_where(xy > -60., nmask).mask
    if newSlice == 'NorthArctic':
        return np.ma.masked_where(xy < 60., nmask).mask
    if newSlice == 'SalArtifact':
        return np.ma.masked_where((xd > 15.) + (xd < 10.), nmask).mask
    if newSlice == 'NitArtifact':
        return np.ma.masked_where((xd > 6.) + (xd < 2.), nmask).mask

    #####
    # Complex Regional masks
    if newSlice == 'BlackSea':
        mx = np.ma.masked_outside(xx, 25.9, 41.7).mask
        my = np.ma.masked_outside(xy, 39.8, 48.1).mask
        return np.ma.masked_where(mx + my, nmask).mask

    if newSlice == 'ignoreBlackSea':
        mx = np.ma.masked_inside(xx, 25.9, 41.7).mask
        my = np.ma.masked_inside(xy, 39.8, 48.1).mask
        return np.ma.masked_where(mx * my, nmask).mask

    if newSlice == 'BalticSea':
        mx = np.ma.masked_outside(xx, 12.5, 30.7).mask
        my = np.ma.masked_outside(xy, 53.0, 66.4).mask
        return np.ma.masked_where(mx + my, nmask).mask

    if newSlice == 'ignoreBalticSea':
        mx = np.ma.masked_inside(xx, 12.5, 30.7).mask
        my = np.ma.masked_inside(xy, 53.0, 66.4).mask
        return np.ma.masked_where(mx * my, nmask).mask

    if newSlice == 'RedSea':
        mx = np.ma.masked_outside(xx, 30.0, 43.0).mask
        my = np.ma.masked_outside(xy, 12.4, 30.4).mask
        return np.ma.masked_where(mx + my, nmask).mask

    if newSlice == 'ignoreRedSea':
        mx = np.ma.masked_inside(xx, 30.0, 43.0).mask
        my = np.ma.masked_inside(xy, 12.4, 30.4).mask
        return np.ma.masked_where(mx * my, nmask).mask

    if newSlice == 'PersianGulf':
        mx = np.ma.masked_outside(xx, 47.5, 56.8).mask
        my = np.ma.masked_outside(xy, 22.3, 32.1).mask
        return np.ma.masked_where(mx + my, nmask).mask

    if newSlice == 'ignorePersianGulf':
        mx = np.ma.masked_inside(xx, 47.5, 56.8).mask
        my = np.ma.masked_inside(xy, 22.3, 32.1).mask
        return np.ma.masked_where(mx * my, nmask).mask

    if newSlice == 'ignoreCaspian':
        mx = np.ma.masked_inside(xx, 45.0, 55.0).mask * np.ma.masked_inside(
            xy, 35., 48.).mask  # caspian
        return np.ma.masked_where(mx, nmask).mask

    if newSlice == 'ignoreMediteranean':
        mx = np.ma.masked_inside(xx, -5.8, 42.5).mask  #E
        my = np.ma.masked_inside(xy, 30., 43.).mask  #N
        mx2 = np.ma.masked_inside(xx, 0., 20.).mask  #E
        my2 = np.ma.masked_inside(xy, 32., 47.).mask  #N
        m = mx * my + mx2 * my2
        return np.ma.masked_where(m, nmask).mask

    if newSlice in ['ignoreInlandSeas', 'IndianOcean']:
        mx = np.ma.masked_inside(xx, 47.5, 56.8).mask * np.ma.masked_inside(
            xy, 22.3, 32.1).mask
        mx += np.ma.masked_inside(xx, 30.0, 43.0).mask * np.ma.masked_inside(
            xy, 12.4, 30.4).mask
        mx += np.ma.masked_inside(xx, 12.5, 30.7).mask * np.ma.masked_inside(
            xy, 53.0, 66.4).mask
        mx += np.ma.masked_inside(xx, 25.9, 41.7).mask * np.ma.masked_inside(
            xy, 39.8, 48.1).mask
        mx += np.ma.masked_inside(xx, -5.8, 42.5).mask * np.ma.masked_inside(
            xy, 30., 43.).mask
        mx += np.ma.masked_inside(xx, 0.0, 20.0).mask * np.ma.masked_inside(
            xy, 32., 47.).mask
        mx += np.ma.masked_inside(xx, 45.0, 55.0).mask * np.ma.masked_inside(
            xy, 35., 52.).mask  # caspian
        if newSlice == 'ignoreInlandSeas':
            return np.ma.masked_where(mx, nmask).mask
        mx += np.ma.masked_outside(xx, 25., 100.).mask
        my = np.ma.masked_outside(xy, -50., 30.).mask
        if newSlice == 'IndianOcean':
            return np.ma.masked_where(mx + my, nmask).mask

    if newSlice == 'AMM':
        return np.ma.masked_outside(bvt.makeLonSafeArr(xx), -20., 13.).mask + np.ma.masked_outside(xy, 40., 65.).mask


    if newSlice == 'SouthernOcean':
        return np.ma.masked_where(xy > -40., nmask).mask
    if newSlice == 'AntarcticOcean':
        return np.ma.masked_where(xy > -50., nmask).mask
    if newSlice == 'ignoreArtics':
        return np.ma.masked_outside(xy, -70., 70.).mask
    if newSlice == 'ignoreMidArtics':
        return np.ma.masked_outside(xy, -65., 65.).mask
    if newSlice == 'ignoreMoreArtics':
        return np.ma.masked_outside(xy, -60., 60.).mask
    if newSlice == 'ignoreExtraArtics':
        return np.ma.masked_outside(xy, -50., 50.).mask
    if newSlice == 'NorthAtlanticOcean':
        return np.ma.masked_outside(bvt.makeLonSafeArr(xx), -80.,
                                    0.).mask + np.ma.masked_outside(
                                        xy, 10., 60.).mask
    if newSlice == 'SouthAtlanticOcean':
        return np.ma.masked_outside(bvt.makeLonSafeArr(xx), -65.,
                                    20.).mask + np.ma.masked_outside(
                                        xy, -50., -10.).mask
    if newSlice == 'EquatorialAtlanticOcean':
        return np.ma.masked_outside(bvt.makeLonSafeArr(xx), -65.,
                                    20.).mask + np.ma.masked_outside(
                                        xy, -15., 15.).mask
    if newSlice == 'ITCZ': #Inter‐Tropical Convergence Zone (johns 2020 Sargassum) in the region 0-15N, 15-55W
        return np.ma.masked_outside(bvt.makeLonSafeArr(xx), -55.,
                                    15.).mask + np.ma.masked_outside(
                                        xy, 0., 15.).mask

    if newSlice == 'ArcticOcean':
        mx = np.ma.masked_where(xy < 60., nmask).mask
        mx += np.ma.masked_inside(xx, -45., 15.).mask * np.ma.masked_inside(
            xy, 50., 80.).mask

        return np.ma.masked_where(mx, nmask).mask

    if newSlice == 'NorthernSubpolarAtlantic':
        mx = np.ma.masked_outside(xx, -74., -3.).mask + np.ma.masked_outside(
            xy, 40., 60.).mask
        mx *= np.ma.masked_outside(xx, -45., 15.).mask + np.ma.masked_outside(
            xy, 60., 80.).mask
        return mx

    if newSlice == 'AtlanticSOcean':
        mx = np.ma.masked_outside(xx, -40., 20.).mask + np.ma.masked_outside(
            xy, -50., -75.).mask
        return mx

    if newSlice == 'NordicSea':
        mx = np.ma.masked_outside(xx, -44., -5.).mask
        mx += np.ma.masked_outside(xy, 53., 65.).mask
        return mx

    if newSlice == 'LabradorSea':
        mx = np.ma.masked_outside(xx, -69., -45.).mask
        mx += np.ma.masked_outside(xy, 53., 67.).mask
        return mx

    if newSlice == 'NorwegianSea':
        mx = np.ma.masked_outside(xx, -15., 10.).mask
        mx += np.ma.masked_outside(xy, 67., 76.).mask
        return mx

    if newSlice == 'Cornwall':
        mx = np.ma.masked_outside(xx, -8., -4.).mask
        mx += np.ma.masked_outside(xy, 49., 52.).mask
        return mx

    if newSlice == 'WeddelSea':
        mx = np.ma.masked_outside(xx, -60., -20.).mask
        mx += np.ma.masked_outside(xy, -80., -64.).mask
        return mx


    # Regions from Pierce 1995 - https://doi.org/10.1175/1520-0485(1995)025<2046:CROHAF>2.0.CO;2
    if newSlice == 'Enderby':
        mx = np.ma.masked_outside(xx, 0., 97.5).mask
        mx += np.ma.masked_outside(xy, -80., -60.).mask
        return mx
    if newSlice == 'Wilkes':
        mx = np.ma.masked_outside(xx, 97.5, 172.5).mask
        mx += np.ma.masked_outside(xy, -80., -60.).mask
        return mx
    if newSlice == 'Ross':
        mx = np.ma.masked_inside(xx, -137.5, 172.5).mask
        mx += np.ma.masked_outside(xy, -80., -60.).mask
        return mx
    if newSlice == 'Amundsen':
        mx = np.ma.masked_outside(xx, -137.5, -72.5).mask
        mx += np.ma.masked_outside(xy, -80., -60.).mask
        return mx
    if newSlice == 'Weddel':
        mx = np.ma.masked_outside(xx, -72.5, 0.).mask
        mx += np.ma.masked_outside(xy, -80., -60.).mask
        return mx

    if newSlice == 'YevgenyNordicSea':
        mx = np.ma.masked_outside(xx, -44., -5.).mask
        mx += np.ma.masked_outside(xy, 53., 65.).mask
        return mx

    if newSlice == 'YevgenyLabradorSea':
        mx = np.ma.masked_outside(xx, -69., -45.).mask
        mx += np.ma.masked_outside(xy, 53., 67.).mask
        return mx

    if newSlice == 'YevgenyNorwegianSea':  # same
        mx = np.ma.masked_outside(xx, -15., 10.).mask
        mx += np.ma.masked_outside(xy, 67., 76.).mask
        return mx

    if newSlice == 'NorthernSubpolarPacific':
        mx = np.ma.masked_inside(xx, -100., 120.).mask
        mx += np.ma.masked_inside(xx, 260., 365.).mask
        mx += np.ma.masked_outside(xy, 40., 60.).mask
        return np.ma.masked_where(mx, nmask).mask

    if newSlice == 'Remainder':
        mx = makeMask(name, 'ignoreInlandSeas', xt, xz, xy, xx, xd)
        mx += np.ma.masked_inside(xy, -10., 10.).mask
        mx += np.ma.masked_outside(abs(xy), -40., 40.).mask
        return np.ma.masked_where(mx, nmask).mask

    if newSlice == 'Equator10':
        mx = makeMask(name, 'ignoreInlandSeas', xt, xz, xy, xx, xd)
        mx += np.ma.masked_outside(xy, -10., 10.).mask
        return mx

    if newSlice == 'NorthPacificOcean':
        mx = np.ma.masked_inside(xx, -100., 120.).mask
        mx += np.ma.masked_inside(xx, 260., 365.).mask
        mx += np.ma.masked_outside(xy, 10., 60.).mask
        return mx

    if newSlice == 'EquatorialPacificOcean':
        mx = np.ma.masked_inside(xx, -83., 120.).mask
        mx += np.ma.masked_inside(xx, 260., 365.).mask
        mx += np.ma.masked_outside(xy, -15., 15.).mask
        return mx

    if newSlice == 'SouthPacificOcean':
        mx = np.ma.masked_inside(xx, -70., 140.).mask
        mx += np.ma.masked_inside(xx, 290., 365.).mask
        my = np.ma.masked_outside(xy, -10., -50.).mask
        return np.ma.masked_where(mx + my, nmask).mask

    if newSlice == 'HighLatWinter':
        NHwinter = np.ma.masked_where(
            ~((xt == months['January']) + (xt == months['February']) +
              (xt == months['March'])), nmask).mask
        SHwinter = np.ma.masked_where(
            ~((xt == months['July']) + (xt == months['August']) +
              (xt == months['September'])), nmask).mask
        print("HighLatWinter masking:\tNHwinter:", NHwinter.sum(), 'SH:',
              SHwinter.sum(), 'of', nmask.sum())

        mnhw = np.ma.masked_where((xy < 45.) + NHwinter, nmask).mask
        mshw = np.ma.masked_where((xy > -45.) + SHwinter, nmask).mask
        return np.ma.masked_where(mnhw * mshw, nmask).mask

    if newSlice == 'CCI_JJA':
        mx = np.ma.masked_where(xy < -53., nmask).mask
        mx += np.ma.masked_where(xy > 80., nmask).mask
        return np.ma.masked_where(mx, nmask).mask

    if newSlice == 'CCI_SON':
        mx = np.ma.masked_where(np.abs(xy) > 80., nmask).mask
        return np.ma.masked_where(mx, nmask).mask

    if newSlice == 'CCI_DJF':
        mx = np.ma.masked_where(xy > 53., nmask).mask
        mx += np.ma.masked_where(xy > 80., nmask).mask
        return np.ma.masked_where(mx, nmask).mask

    if newSlice in ['Depth_700m', 'Depth_2000m', 'Depth_700-2000m']:
        print(newSlice, xz.min(), xz.mean(), xz.max(), len(xz))
        print('Depth_700m',
              len(np.ma.masked_where(abs(xz) > 700., nmask).mask))
        print('Depth_2000m',
              len(np.ma.masked_where(abs(xz) > 2000., nmask).mask))
        print(
            'Depth_700-2000m',
            len(
                np.ma.masked_where((abs(xz) < 700.) * (abs(xz) > 2000.),
                                   nmask).mask))
        assert 0

    #####
    # Depths masks
    if newSlice == 'Depth_0-10m':
        return np.ma.masked_where(abs(xz) > 10., nmask).mask
    if newSlice == 'Depth_10-20m':
        return np.ma.masked_where((abs(xz) < 10.) + (abs(xz) > 20.),
                                  nmask).mask
    if newSlice == 'Depth_20-50m':
        return np.ma.masked_where((abs(xz) > 20.) + (abs(xz) > 50.),
                                  nmask).mask
    if newSlice == 'Depth_50-100m':
        return np.ma.masked_where((abs(xz) < 50.) + (abs(xz) > 100.),
                                  nmask).mask
    if newSlice == 'Depth_100-500m':
        return np.ma.masked_where((abs(xz) < 100.) + (abs(xz) > 500.),
                                  nmask).mask
    if newSlice == 'Depth_500m':
        return np.ma.masked_where(abs(xz) < 500., nmask).mask
    if newSlice == 'Depth_700m':
        return np.ma.masked_where(abs(xz) > 700., nmask).mask
    if newSlice == 'Depth_0-50m':
        return np.ma.masked_where(abs(xz) > 50., nmask).mask
    if newSlice == 'Depth_50-100m':
        return np.ma.masked_where((abs(xz) < 50.) + (abs(xz) > 100.),
                                  nmask).mask
    if newSlice == 'Depth_100-200m':
        return np.ma.masked_where((abs(xz) < 100.) + (abs(xz) > 200.),
                                  nmask).mask
    if newSlice == 'Depth_200-500m':
        return np.ma.masked_where((abs(xz) < 200.) + (abs(xz) > 500.),
                                  nmask).mask
    if newSlice == 'Depth_500-1000m':
        return np.ma.masked_where((abs(xz) < 500.) + (abs(xz) > 1000.),
                                  nmask).mask
    if newSlice == 'Depth_700-2000m':
        return np.ma.masked_where((abs(xz) < 700.) * (abs(xz) > 2000.),
                                  nmask).mask
    if newSlice == 'Depth_1000-2000m':
        return np.ma.masked_where((abs(xz) < 1000.) + (abs(xz) > 2000.),
                                  nmask).mask
    if newSlice == 'Depth_1000m':
        return np.ma.masked_where(abs(xz) < 1000., nmask).mask
    if newSlice == 'Depth_2000m':
        return np.ma.masked_where(abs(xz) > 2000., nmask).mask
    if newSlice == 'Shallow': return np.ma.masked_where(xz > 200., nmask).mask
    if newSlice == 'Depth': return np.ma.masked_where(xz < 200., nmask).mask

    #####
    # Masks that require bathymetry:
    if newSlice in [
            'maskBelowBathy',
            'OnShelf',
            'OffShelf',
            'AMM_Shelf',
            'AMM_OffShelf',
    ]:
        if newSlice in ['AMM_Shelf', 'AMM_OffShelf']:
            bathync = Dataset(paths.orcaGridfn,'r')
            print(paths.orcaGridfn)
            latcc = bathync.variables["nav_lat"][:]
            loncc = bathync.variables["nav_lon"][:]
            deptht = np.abs(bathync.variables["nav_lev"][:])
            bathy = bathync.variables["hbatt"][:].squeeze()
            bathync.close()
            shelfDepth = 200.
        else:
            bathync = Dataset(paths.orca1bathy, 'r')
            bathy = np.abs(bathync.variables["bathymetry"][:])
            latcc, loncc = bathync.variables["lat"][:], bathync.variables["lon"][:]
            bathync.close()
            shelfDepth = 250.

        shelveFn = bvt.folder(os.path.join(paths.shelvedir, "MatchingMasks/"))+ newSlice+"_diag_maskMask.shelve"

        try:
            s = shOpen(shelveFn)
            lldict = s['lldict']
            s.close()
        except:
            lldict = {}

        print("Bathy mask: before mask:", newSlice, nmask.sum(), 'of',
              len(nmask))
        i = 0
        for i, z in enumerate(xz):
            try:
                la, lo = lldict[(xy[i], xx[i])]
            except:
                la, lo = bvt.getOrcaIndexCC(xy[i],
                                        xx[i],
                                        latcc,
                                        loncc,
                                        debug=False)
                lldict[(xy[i], xx[i])] = la, lo
            if la == lo == -1:
                print("Corner case:", la, lo, bathy[la, lo])
                nmask[i] = 1
            if newSlice == "maskBelowBathy":
                if (bathy[la, lo] - 10.) > abs(z): 
                    nmask[i] = 1
            elif newSlice in ["OnShelf", 'AMM_Shelf']:
                if bathy[la, lo] >= shelfDepth: 
                    nmask[i] = 1
            elif newSlice in ["OffShelf", 'AMM_OffShelf']:
                if bathy[la, lo] < shelfDepth: 
                    nmask[i] = 1
        if i > 0:
            try:
                s = shOpen(shelveFn)
                s['lldict'] = lldict
                s.close()
            except:
                print(
                    "makeMask:\tWARNING:\tUnable to save lldict at this time")
        print("Bathy mask:", newSlice, nmask.sum(), 'of', len(nmask))
        return nmask

    #####
    # Time masks:
    months = {month_name[i + 1]: i for i in range(0, 12)}
    if newSlice in list(months.keys()):
        print("masking a month:", newSlice, xt[0], xt[-1])
        return np.ma.masked_where(xt != months[newSlice], nmask).mask
    if newSlice == 'JFM':
        return np.ma.masked_where(
            ~(xt == months['January']) + (xt == months['February']) +
            (xt == months['March']), nmask).mask
    if newSlice == 'AMJ':
        return np.ma.masked_where(
            ~(xt == months['April']) + (xt == months['May']) +
            (xt == months['June']), nmask).mask
    if newSlice == 'JAS':
        return np.ma.masked_where(
            ~(xt == months['July']) + (xt == months['August']) +
            (xt == months['September']), nmask).mask
    if newSlice == 'OND':
        return np.ma.masked_where(
            ~(xt == months['October']) + (xt == months['November']) +
            (xt == months['December']), nmask).mask

    #####
    # Old rarely used masks
    if newSlice in ['1-99pc', '5-95pc', '0-99pc'] or newSlice in [
            '0-1pc',
            '1-5pc',
            '5-25pc',
            '25-40pc',
            '40-60pc',
            '60-75pc',
            '75-95pc',
            '95-99pc',
            '99-100pc',
    ]:
        if newSlice in [
                '0-1pc',
                '1-5pc',
                '5-25pc',
                '25-40pc',
                '40-60pc',
                '60-75pc',
                '75-95pc',
                '95-99pc',
                '99-100pc',
        ]:
            tmp = newSlice.replace('pc', '').split('-')
            pcmin, pcmax = float(tmp[0]), float(tmp[1])
            print(newSlice, pcmin, pcmax)
            if pcmin == 0: ymin = yd.min()
            else: ymin = scoreatpercentile(yd, pcmin)
            if pcmax == 100: ymax = yd.max()
            else: ymax = scoreatpercentile(yd, pcmax)

        if newSlice in [
                '1-99pc',
        ]:
            ymin = scoreatpercentile(xd, 1)
            ymax = scoreatpercentile(xd, 99)
        if newSlice in [
                '5-95pc',
        ]:
            ymin = scoreatpercentile(xd, 5)
            ymax = scoreatpercentile(xd, 95)

        if newSlice in [
                '0-99pc',
        ]:
            ymin = xd.min()
            ymax = scoreatpercentile(xd, 99)
        print("makeMask:\t", newSlice, ymin, ymax)
        return np.ma.masked_outside(xd, ymin, ymax).mask

    if newSlice == "0.1": return np.ma.masked_where(xd == 0.1, xd).mask
    if newSlice == "0.2": return np.ma.masked_where(xd == 0.2, xd).mask
    if newSlice == "0.01": return np.ma.masked_where(xd == 0.01, xd).mask
    if newSlice == 'Zoom': return np.ma.masked_where(xd > 10., nmask).mask
    if newSlice == 'Zoom5': return np.ma.masked_where(xd > 5., nmask).mask
    if newSlice == 'Zoom2': return np.ma.masked_where(xd > 2., nmask).mask
    if newSlice == 'TypicalIron':
        return np.ma.masked_where((xd <= 0.) * (xd <= 4.), nmask).mask

    if newSlice in [
            'OffAxis',
            'Overestimate',
            'Underestimate',
            'Matched',
            'Overestimate_2sig',
            'Underestimate_2sig',
            'Overestimate_3sig',
            'Underestimate_3sig',
    ]:
        print("makeMask:\tSlice", newSlice,
              "requires both datasets, and you should never see this")
        assert False

    ####################
    # Please Add your own masks here:

    # End of custom maps.
    ###################

    print("Mask region not accepted:", newSlice)
    assert False
