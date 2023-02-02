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

# This is a long list of possible values for new Slice used in bgcvaltools.makeMask.
# these are all included in the UKESMPython.slicesDict dictonairy. #
# This list is deprecated now, as it's largely unused.
#
import numpy as np
from calendar import month_name


def getSlicesDict():
    """
	Produces a dictionary of slices. lices is another name for the regional cuts.
	"""
    slicesDict = {}
    standardCuts = [
        '5-95pc',
        'ignoreInlandSeas',
        'OffShelf',
        'ignoreExtraArtics',
        'aboveZero',
    ]
    months = [month_name[i + 1] for i in range(0, 12)
              ]
    depthRanges = [
        'Depth_0-50m', 'Depth_50-100m', 'Depth_100-200m', 'Depth_200-500m',
        'Depth_500-1000m', 'Depth_1000-2000m', 'Depth_2000m'
    ]
    percentiles = [
        '0-1pc',
        '1-5pc',
        '5-25pc',
        '25-40pc',
        '40-60pc',
        '60-75pc',
        '75-95pc',
        '95-99pc',
        '99-100pc',
    ]
    latregions = [
        'NorthTemperate', 'SouthTemperate', 'NorthTropics', 'Equatorial',
        'SouthTropics', 'Antarctic', 'NorthArctic', 'Arctic', 'Tropics',
        'Temperate'
    ]
    Hemispheres = [
        'NorthHemisphere',
        'SouthHemisphere',
    ]
    Seas = [
        'ignoreMediteranean',
        'BlackSea',
        'ignoreBlackSea',
        'RedSea',
        'BalticSea',
        'PersianGulf',
        'ignoreInlandSeas',
        'ignoreRedSea',
        'ignoreBalticSea',
        'ignorePersianGulf',
    ]
    Oceans = [
        'SouthPacificOcean',
        'ArcticOcean',
        'AntarcticOcean',
        'NorthAtlanticOcean',
        'SouthAtlanticOcean',
        'NorthPacificOcean',
        'IndianOcean',
        'EquatorialPacificOcean',
        'EquatorialAtlanticOcean',
    ]  #'ignoreExtraArtics','ignoreMidArtics','ignoreArtics','ignoreMoreArtics',]
    QualityCuts = [
        'Overestimate',
        'Underestimate',
        'Overestimate_2sig',
        'Underestimate_2sig',
        'Overestimate_3sig',
        'Underestimate_3sig',
        'Matched',
        'OffAxis',
        '1-99pc',
        '5-95pc',
        '0-99pc',
    ]
    Seasons = ['JFM', 'AMJ', 'JAS', 'OND']
    Transects = [
        'AtlanticTransect', 'PacificTransect', 'SouthernTransect', '10N', '10S'
    ]
    Misc = [
        'HighLatWinter',
    ]
    BGCVal = [
        'Global',
        'ignoreInlandSeas',
        'Equator10',
        'ArcticOcean',
        'NorthernSubpolarAtlantic',
        'NorthernSubpolarPacific',
        'SouthernOcean',
        'Remainder',
        'AMM',
    ]
    AMM = ['AMM','AMM_Shelf', 'AMM_OffShelf']

    OceanMonths = {o: [(o, m) for m in months] for o in Oceans}
    OceanSeasons = {o: [(o, m) for m in Seasons] for o in Oceans}
    HemispheresMonths = {o: [(o, m) for m in months] for o in Hemispheres}
    HemispheresSeasons = {o: [(o, m) for m in Seasons] for o in Hemispheres}

    newSlices = [
        'All',
        'Standard',
    ]
    newSlices.extend(months)
    newSlices.extend(depthRanges)
    newSlices.extend(percentiles)
    newSlices.extend(latregions)
    newSlices.extend(QualityCuts)
    newSlices.extend(Seas)
    newSlices.extend(Oceans)
    newSlices.extend(Hemispheres)
    newSlices.extend(Seasons)
    newSlices.extend(OceanSeasons)
    newSlices.extend(OceanMonths)
    newSlices.extend(HemispheresMonths)
    newSlices.extend(HemispheresSeasons)
    newSlices.extend(Transects)
    newSlices.extend(Misc)
    newSlices.extend(BGCVal)
    newSlices.extend(AMM)

    for om, keys in list(OceanMonths.items()):
        newSlices.extend(keys)
    for om, keys in list(OceanSeasons.items()):
        newSlices.extend(keys)
    for om, keys in list(HemispheresMonths.items()):
        newSlices.extend(keys)
    for om, keys in list(HemispheresSeasons.items()):
        newSlices.extend(keys)

    slicesDict['Default'] = [
        'All',
        'Standard',
    ]
    slicesDict['StandardCuts'] = standardCuts
    slicesDict['AllSlices'] = newSlices
    slicesDict['Months'] = months
    slicesDict['Hemispheres'] = Hemispheres
    slicesDict['Oceans'] = Oceans
    slicesDict['Seasons'] = Seasons
    slicesDict['depthRanges'] = depthRanges
    slicesDict['Transects'] = Transects
    slicesDict['Misc'] = Misc
    slicesDict['BGCVal'] = BGCVal
    slicesDict['AMM'] = AMM
    for om, keys in list(OceanMonths.items()):
        slicesDict[om + 'Months'] = keys
    for om, keys in list(OceanSeasons.items()):
        slicesDict[om + 'Seasons'] = keys
    for om, keys in list(HemispheresMonths.items()):
        slicesDict[om + 'Months'] = keys
    for om, keys in list(HemispheresSeasons.items()):
        slicesDict[om + 'Seasons'] = keys
    return slicesDict

slicesDict = getSlicesDict()


def populateSlicesList(  #plotallcuts = False,
    plotDefaults=True,
    plotMonths=0,  #True
    plotdepthRanges=0,  #True	
    plotpercentiles=0,  #True	
    plotLatRegions=0,  #True
    plotQualityCuts=0,  #True
    plotSeas=0,  #True		 
    plotOceans=0,  #True	
    plotHemispheres=0,
    plotSeasons=0,  # True
    plotOceanSeasons=0,  # True		 
    plotOceanMonths=0,
    plotHemispheresMonths=0,
    plotHemispheresSeasons=0,
    plotTransects=0,
    plotMisc=0,
):
    """
	Takes a series of boolean flags, then creates a list of slices for point to point analyses.
	"""

    if plotDefaults: newSlices = [
            'All',
            'Standard',
    ]  # Defaults
    else: newSlices = []
    if plotMonths: newSlices.extend(slicesDict['Months'])
    if plotdepthRanges: newSlices.extend(slicesDict['depthRanges'])
    if plotpercentiles: newSlices.extend(slicesDict['percentiles'])
    if plotLatRegions: newSlices.extend(slicesDict['latregions'])
    if plotQualityCuts: newSlices.extend(slicesDict['QualityCuts'])
    if plotSeas: newSlices.extend(slicesDict['Seas'])
    if plotOceans: newSlices.extend(slicesDict['Oceans'])
    if plotHemispheres: newSlices.extend(slicesDict['Hemispheres'])
    if plotSeasons: newSlices.extend(slicesDict['Seasons'])
    if plotTransects: newSlices.extend(slicesDict['Transects'])
    if plotMisc: newSlices.extend(slicesDict['Misc'])

    if plotOceanMonths:
        for o in slicesDict['Oceans']:
            newSlices.extend(slicesDict[o + 'Months'])
    if plotOceanSeasons:
        for o in slicesDict['Oceans']:
            newSlices.extend(slicesDict[o + 'Seasons'])
    if plotHemispheresMonths:
        for o in slicesDict['Hemispheres']:
            newSlices.extend(slicesDict[o + 'Months'])
    if plotHemispheresSeasons:
        for o in slicesDict['Hemispheres']:
            newSlices.extend(slicesDict[o + 'Seasons'])

    return newSlices
