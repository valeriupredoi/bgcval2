#
# Copyright 2015, Plymouth Marine Laboratory
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
.. module:: patternAnalyses
   :platform: Unix
   :synopsis: A toolkit containing a selection of ways to make patterns plots, according to the circumstances.
.. moduleauthor:: Lee de Mora <ledm@pml.ac.uk>

"""

from . import makePatternStatsPlots, makeTargets
from bgcval2.bgcvaltools.UKESMpython import folder, reducesShelves, listShelvesContents
from ..bgcvaltools.pftnames import months, Ocean_names, SouthHemispheresMonths, NorthHemispheresMonths
from .slicesDict import populateSlicesList, slicesDict


def InterAnnualPatterns(shelvesAV, jobID, years, Grid):

    #####
    # InterAnnual analysis
    allshelves = listShelvesContents(shelvesAV)

    for model in allshelves.models:
        for name in allshelves.names:
            for depthLevel in allshelves.depthLevels:
                print('InterAnnual analysis:', model, name, depthLevel)

                #####
                # Each year is a different colour.
                yearslyShelves = {}
                for y in years:
                    yearslyShelves[y] = reducesShelves(
                        shelvesAV,
                        models=[model],
                        names=[
                            name,
                        ],
                        years=[
                            y,
                        ],
                        depthLevels=[
                            depthLevel,
                        ],
                        sliceslist=months,
                    )

                filenamebase = folder(
                    'images/' + model + '-' + jobID +
                    '/Patterns/Interannual/' + name +
                    depthLevel) + model + jobID + name + '_' + years[
                        0] + '-' + years[-1] + depthLevel
                makePatternStatsPlots(
                    yearslyShelves,  # {legend, shelves}
                    model + ' ' + name,  # title
                    months,  # xkeysLabels
                    filenamebase,  # filename base
                    grid=Grid,
                )

                #####
                # One long line covering multiple years - monthly. (looks best with up to 4 years.)
                longMonths, longShelves = [], []
                for y in years:
                    shelves = yearslyShelves[y]
                    longMonths.extend([(mn, y) for mn in months])
                    if len(shelves) == 12: longShelves.extend(shelves)
                    else:
                        print("Shelves are not 12 long:", len(shelves))
                        continue
                        assert False
                if len(longMonths) != len(longShelves):
                    print("len(longMonths) != len(longShelves):",
                          len(longMonths), ' != ', len(longShelves))
                    continue
                    assert False

                filenamebase = folder(
                    'images/' + model + '-' + jobID +
                    '/Patterns/Interannual/' + name +
                    depthLevel) + model + jobID + name + '_' + years[
                        0] + '-' + years[-1] + depthLevel + '_longtimeseries'
                print("makePatternStatsPlots:", {name: longShelves},
                      model + ' ' + name, longMonths, filenamebase, Grid)
                makePatternStatsPlots(
                    {name: longShelves},  # {legend, shelves}
                    model + ' ' + name + ' ' + depthLevel,  # title
                    longMonths,  # xkeysLabels
                    filenamebase,  # filename base
                    grid=Grid,
                )

                #####
                # One long line showing multiple years - annual.
                for sl in ['All', 'Standard']:
                    yearslyShelves = reducesShelves(
                        shelvesAV,
                        models=[model],
                        names=[
                            name,
                        ],
                        years=years,
                        depthLevels=[
                            depthLevel,
                        ],
                        sliceslist=[
                            sl,
                        ],
                    )

                    if len(yearslyShelves) != len(years):
                        print("len(yearslyShelves) != len(years):",
                              len(yearslyShelves), ' != ', len(years))
                        assert False

                    filenamebase = folder(
                        'images/' + model + '-' + jobID +
                        '/Patterns/Interannual/' + name +
                        depthLevel) + model + jobID + name + '_' + years[
                            0] + '-' + years[-1] + depthLevel + '_annual' + sl
                    print("makePatternStatsPlots - Annual :", sl,
                          {name: sorted(yearslyShelves)}, years,
                          model + ' ' + name, filenamebase, Grid)
                    makePatternStatsPlots(
                        {name: sorted(yearslyShelves)},  # {legend, shelves}
                        model + ' ' + name + ' ' + depthLevel,  # title
                        sorted(years),  # xkeysLabels
                        filenamebase,  # filename base
                        grid=Grid,
                    )


def BGCvsPhysics(shelvesAV, jobID, Grid, physicsModel='NEMO'):
    allshelves = listShelvesContents(shelvesAV)
    if len(allshelves.models) < 2:
        print("BGCvsPhysics:\tNot enough models")

    print(allshelves)
    if physicsModel not in allshelves.models:
        print("BGCvsPhysics:\tNo NEMO model")

    model = allshelves.models
    if len(model) == 2:
        model.remove('NEMO')
    if len(model) == 1 and type(model) == type([
            'a',
    ]):
        model = model[0]

    bgcnames = allshelves.names
    for b in ['MLD', 'mld', 'temperature', 'salinity']:
        try:
            bgcnames.remove(b)
        except:
            pass

    print("BGCvsPhysics:\tStarting:", model, jobID, bgcnames)
    for year in allshelves.years:
        for depthLevel in allshelves.depthLevels:
            for name in bgcnames:
                for slicekey, slices in zip([
                        'months',
                        'oceans',
                        'SHmonths',
                        'NHmonths',
                ], [
                        months, Ocean_names, SouthHemispheresMonths,
                        NorthHemispheresMonths
                ]):
                    monthlyShelves = {}
                    monthlyShelves['MLD'] = reducesShelves(
                        shelvesAV,
                        models=[
                            physicsModel,
                        ],
                        names=[
                            'mld',
                        ],
                        years=[
                            year,
                        ],
                        sliceslist=slices,
                    )
                    monthlyShelves[name] = reducesShelves(
                        shelvesAV,
                        models=allshelves.models,
                        names=[
                            name,
                        ],
                        years=[
                            year,
                        ],
                        sliceslist=slices,
                        depthLevels=[
                            depthLevel,
                        ])
                    for nm in ['temperature', 'salinity']:
                        monthlyShelves[nm] = reducesShelves(
                            shelvesAV,
                            models=[
                                physicsModel,
                            ],
                            names=[
                                nm,
                            ],
                            years=[
                                year,
                            ],
                            depthLevels=[
                                depthLevel,
                            ],
                            sliceslist=slices,
                        )

                    filenamebase = folder(
                        'images/' + model + '-' + jobID + '/Patterns/' + year +
                        '/' + name + '_vsPhysics'
                    ) + model + jobID + name + '_' + year + depthLevel + '_' + slicekey
                    print("BGCvsPhysics:\t", year, model + ' ' + name,
                          filenamebase, Grid)
                    print("monthlyShelves:", monthlyShelves)

                    makePatternStatsPlots(
                        monthlyShelves,  # {legend, shelves}
                        model + ' ' + name + ' ' + depthLevel + ' ' +
                        year,  # title
                        slices,  # xkeysLabels
                        filenamebase,  # filename base
                        grid=Grid,
                    )


def onePatternAtATime(allshelves, ):
    ### produces a simple pattern plot with one line. and a taylor/target lot
    groups = {
        'Oceans': [],
        'Months': [],
        'Seasons': [],
        'NorthHemisphereMonths': [],
        'SouthHemisphereMonths': [],
        'depthRanges': []
    }
    allshelves = listShelvesContents(shelvesAV)
    for year in allshelves.years:
        for name in allshelves.names:
            for depthLevel in allshelves.depthLevels:
                #####
                # Produce a set of pattern and a target plots for each of the groups here.

                for g in groups:
                    groups[g] = reducesShelves(shelvesAV,
                                               models=[
                                                   model,
                                               ],
                                               depthLevels=[
                                                   depthLevel,
                                               ],
                                               names=[
                                                   name,
                                               ],
                                               sliceslist=slicesDict[g])
                    print(g, groups[g])

                    if len(groups[g]) == 0: continue

                    #####
                    # makeTargets:
                    # Make a target diagram of the shelves of this group.
                    filename = folder(
                        imageFolder + '/Targets/' + year + '/' + name +
                        depthLevel + '/' + g
                    ) + model + '-' + jobID + '_' + year + '_' + name + depthLevel + '_' + g + '.png'
                    makeTargets(
                        groups[g],
                        filename,
                        legendKeys=[
                            'newSlice',
                        ],
                    )
                    #####
                    # makePattern plots:
                    # Make a pattern  diagram of all matches for this particular dataset.
                    xkeys = ''
                    for o in ['Oceans', 'Months', 'depthRanges']:
                        if g.find(o) >= 0: xkeys = o
                    if xkeys == '':
                        print("Could no find x axis keys!", g, 'in',
                              ['Oceans', 'Months'])

                    filenamebase = folder(
                        imageFolder + '/Patterns/' + year + '/' + name +
                        depthLevel + '/' + g
                    ) + 'Months-' + model + '-' + jobID + '_' + year + '_' + name + depthLevel
                    makePatternStatsPlots(
                        {
                            name: groups[g],
                        },  # {legend, shelves}
                        name + ' ' + g,  #xkeysname
                        slicesDict[xkeys],  #xkeysLabels=
                        filenamebase,  # filename base	
                        grid=grid,
                    )
                    #####
                    # After finding all the shelves, we can plot them on the same axis.
                    filenamebase = folder(
                        imageFolder + '/Patterns/' + year + '/' + name +
                        depthLevel + '/ANSH'
                    ) + 'ANSH-Months-' + model + '-' + jobID + '_' + year + '_' + name + depthLevel

                    makePatternStatsPlots(
                        {
                            'North Hemisphere':
                            groups['NorthHemisphereMonths'],
                            'South Hemisphere':
                            groups['SouthHemisphereMonths'],
                            'Global': groups['Months'],
                        },  # {legend, shelves}
                        name + ' Months',  #xkeysname
                        slicesDict['Months'],  #xkeysLabels=
                        filenamebase,  # filename base	
                        grid=grid,
                    )


def modelIntercomparisonAnnual(shelvesAV, imageFolder):

    allshelves = listShelvesContents(shelvesAV)
    models = allshelves.models
    names = allshelves.names
    dls = allshelves.depthLevels
    years = allshelves.years
    plots = ['Oceans', 'depthRanges', 'Transects']

    for year in years:
        for name in names:
            for dl in dls:
                for plot in plots:
                    data = {}
                    for model in models:
                        print("modelIntercomparisonAnnual:", year, name, dl,
                              plot, model, slicesDict[plot])
                        data[model] = reducesShelves(
                            shelvesAV,
                            models=[
                                model,
                            ],
                            depthLevels=[
                                dl,
                            ],
                            names=[
                                name,
                            ],
                            sliceslist=slicesDict[plot])

                    filenamebase = folder(
                        imageFolder + '/Patterns/' + year + name + dl + '/' +
                        plot
                    ) + 'modelIntercomparisonAnnual-' + year + '_' + name + dl + plot
                    makePatternStatsPlots(
                        data,  # {legend, shelves}
                        plot,  # xkeysname
                        slicesDict[plot],  # xkeysLabels
                        filenamebase,  # filename base
                        grid='Flat1deg',
                    )
