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
.. module:: testsuite_p2p
   :platform: Unix
   :synopsis: The tool that does the legwork for the point to point analysis.
.. moduleauthor:: Lee de Mora <ledm@pml.ac.uk>

"""
#Standard Python modules:
from sys import argv, exit
from os.path import exists
from calendar import month_name

#Specific local code:
from .. import UKESMpython as ukp
from . import matchDataAndModel, makeTargets, makePatternStatsPlots
from .p2pPlots import makePlots
from .slicesDict import populateSlicesList, slicesDict
from ..bgcvaltools.pftnames import MaredatTypes, WOATypes, Ocean_names, OceanMonth_names, months, Seasons, Hemispheres, HemispheresMonths, OceanSeason_names

###	Potential problems?
###		Reliance on ORCA1 grid


def testsuite_p2p(
    model='ERSEM',  #'MEDUSA','ERSEM','NEMO'],
    year='1997',
    jobID='xhonc',
    av={},
    plottingSlices=[],
    workingDir='',
    imageFolder='',
    noPlots=False,
    gridFile='',
    annual='',
    noTargets=False,
):
    """
	This analysis package performs the point to point analysis for a single model, for one year, for one job ID.
	
	Arguments:
	    model: Model name
	    
	    year: year (4 digit string), doesn't need to be exact. For climatologies, we typically use 2525 or something absurd.
	    
	    jobID: 5 letter jobID as used on monsoon. 
	    
	    av:
		The AutoVivification (av) is crucial. It controls the analysis. It locates the files. It identifies the field names in the netcdf.
	
		The av has a few very specific requiments in terms of structure and key words.
		Here is an example, of ERSEM chlorophyll:
			av['chl']['Data']['File'] 	= Observation_Filename_path_.netcdf	
			av['chl']['Data']['Vars'] 	= ['Chlorophylla',]		
			av['chl']['ERSEM']['File'] 	= Model_Filename_path_.netcdf
			av['chl']['ERSEM']['Vars'] 	= ['chl',]
			av['chl']['ERSEM']['grid']	= 'ORCA1'				
			av['chl']['layers'] 	= ['',]
		where: 
			'File' is the file path
			'Vars' is the variable as it is call in the netcdf.
			'grid' is the model grid name. These grids are linked to a grid mesh file for calculating cell volume, surface area and land masks.
			'layers' a list of depth levels. This is needed because some WOA files are huges and desktop computers may run the p2p analysis of that data.
				layers options are ['', 'Surface','100m','200m','500m',]
				layers = ['',] indicates look at no depth slicing. (Not recommended for big (>500,000 points) datasets! Don't say I didn't warn you.)

	    plottingSlices:
		plottingSlices is a list of regional, temporal, or statistical slices to be given to the analysis for plotting.
		ie:	plottingSlices = ['All', 				# plots everything,
					  'NorthAtlantic', 			# plots North Atlantic
					  'February',				# plots February
					  ('NorthAtlantic', 'February'),	# plots North Atlantic in February
					  ]
			plottingSlices can be made automatically with UKESMpthon.populateSlicesList()
			plottingSlices can also be added to the av:
				av[name]['plottingSlices'] = a list of slices
	    workingDir: 
	    	workingDir is a location for the working files that are produced during the analysis. 
	    	if no working directory is provided, the default is: ~/WorkingFiles/model-jobID-yyear
	    		
	    imageFolder: 
	    	imageFolder is a location for all the images  that are produced during the analysis. 
	    	if no working directory is provided, the default is: ~/images/model-jobID  	

	    noPlots:
	    	noPlots is a boolean value to turn off the production of images.
	    	This can streamline the analysis routine, if plots are not needed.
	    	
	Returns:
		shelvesAV:
		another AutoVivification with the following structure:
		shelvesAV[model][name][depthLevel][newSlice][xkey] = shelvePath
	
	testsuite_p2p is not the place for intercomparisons of models, years, or jobID. 
	This can be done after calling testsuite_p2p.py, and by using the 

	"""

    print("#############################")
    print("testsuite_p2p:  ")
#    print("models:        ", model)
    print("year:          ", year)
    print("jobID:         ", jobID)
    print("av keys:	      ", sorted(av.keys()))
    print("#############################")

    if len(list(av.keys())) == 0:
        print(
            "No autovivification nested dictionary given. - See testsuite_p2p documentation or a working example."
        )
        exit(0)

    # Location of processing files
    if len(workingDir) == 0:
        workingDir = ukp.folder(''.join(["WorkingFiles/", jobID, '-', year]))
        print("No working directory provided, creating default:", workingDir)

    # Location of image Output files
    if noPlots is False:
        if len(imageFolder) == 0:
            imageFolder = ukp.folder(''.joib(['images/', jobID]))
            print("No image directory provided, creating default:",
                  imageFolder)

    #####
    # Start analysis here:
    shelvesAV = []  #AutoVivification()

    for name in sorted(av.keys()):
#       #####
#       # Start with some tests of the av.
#
#       #####
#       # Testing av for presence of model keyword
#       print("testsuite_p2p: \t", model, jobID, year,
#             name)  #, av[name][model]
#       try:
#           if not isinstance(av[name][model], dict):
#               print("testsuite_p2p: \tWARNING:", model, ' not in av',
#                     list(av[name].keys()))
#               continue
#           if len(list(av[name][model].keys())) == 0:
#               print("testsuite_p2p: \tWARNING:", model, ' not in av',
#                     list(av[name].keys()))
#               continue
#       except KeyError:
#           print("testsuite_p2p: \tWARNING:\tNo ", name, 'in ', model)
#           continue

    #####
    # Testing av for presence of data keyword
#        try:
#            if not isinstance(av[name]['Data'], dict):
#                print("testsuite_p2p: \tWARNING:", 'Data', ' not in av',
#                      list(av[name].keys()))
#                continue
#            if len(list(av[name]['Data'].keys())) == 0:
#                print("testsuite_p2p: \tWARNING:", 'Data', ' not in av',
#                      list(av[name].keys()))
#                continue
#        except KeyError:
#            print("testsuite_p2p: \tWARNING:\tNo ", 'Data', 'in ', jobID)
#            continue

    #####
#    # Testing av for presence of data/obs files.
#        try:
#            if not exists(av[name]['Data']['File']):
#                print("testsuite_p2p.py:\tWARNING:\tFile does not exist",
#                      av[name]['Data']['File'])
#                continue
#        except:
#            print(
#                "testsuite_p2p.py:\tWARNING:\tDict entry does not exist\tav[",
#                name, "][", jobID, '][File]')
#            continue
#        try:
#            if not exists(av[name][jobID]['File']):
#                print("testsuite_p2p.py:\tWARNING:\tFile does not exist",
#                      av[name][jobID]['File'])
#                continue
#        except:
#            print(
#                "testsuite_p2p.py:\tWARNING:\tDict entry does not exist:\tav[",
#                name, "][", jobID, '][File] :', av[name][jobID]['File'])
#            continue

        #####
    # Testing av for presence of grid
#        grid = av[name][model]['grid']
#        if grid in ['', [], {}, None]:
#            print("testsuite_p2p.py:\tERROR:\tgrid not found:\tav[", name,
#                  "][", model, '][grid]: ', grid)
#            assert False
#
#        #####
#    # Testing av for presence of layers
#        if len(av[name]['layers']) == 0:
#            av[name]['layers'] = [
#                '',
#            ]
#            print(
#                "testsuite_p2p: \tWARNING: no 'layers' provided in av, using defaults: ['',]"
#            )

    #####
    # Made it though the initial tests. Time to start the analysis.
        print(
            "\n\n\ntestsuite_p2p.py:\tINFO:\tMade it though initial tests. Running:",
            jobID, year, name, av[name]['layers'])
        for depthLevel in av[name]['layers']:
            depthLevel = str(depthLevel)

            #####
            # matchDataAndModel:
            # Match observations and model.
            # Does not produce and plots.
            b = matchDataAndModel(av[name]['Data']['File'],
                                  av[name][model]['File'],
                                  dataType=name,
                                  modelcoords=av[name][model]['coords'],
                                  modeldetails=av[name][model]['details'],
                                  datacoords=av[name]['Data']['coords'],
                                  datadetails=av[name]['Data']['details'],
                                  datasource=av[name]['Data']['source'],
                                  model=av[name][model]['source'],
                                  jobID=jobID,
                                  year=year,
                                  workingDir=ukp.folder(workingDir + name),
                                  depthLevel=depthLevel,
                                  grid=grid,
                                  gridFile=gridFile)

            #####
            # makePlots:
            # Make some plots of the point to point datasets.
            # MakePlot runs a series of analysis, comparing every pair in DataVars and ModelVars
            #	 under a range of different masks. For instance, only data from Antarctic Ocean, or only data from January.
            # The makePlot produces a shelve file in workingDir containing all results of the analysis.
            if len(plottingSlices) == 0:
                if len(av[name]['plottingSlices']) == 0:
                    nplottingSlices = populateSlicesList()
                    print("No plotting slices provided, using defaults",
                          nplottingSlices)
                else:

                    nplottingSlices = av[name]['plottingSlices']
                    print("Plotting slices provided, using ", nplottingSlices)
            else:
                nplottingSlices = plottingSlices

            imageDir = ukp.folder(imageFolder + 'P2Pplots/' + year + '/' +
                                  name + depthLevel)
            m = makePlots(b.MatchedDataFile,
                          b.MatchedModelFile,
                          name,
                          newSlices=nplottingSlices,
                          model=av[name][model]['source'],
                          datasource=av[name]['Data']['source'],
                          jobID=jobID,
                          depthLevel=depthLevel,
                          year=year,
                          modelcoords=av[name][model]['coords'],
                          modeldetails=av[name][model]['details'],
                          datacoords=av[name]['Data']['coords'],
                          datadetails=av[name]['Data']['details'],
                          shelveDir=ukp.folder(workingDir + name + depthLevel),
                          imageDir=imageDir,
                          compareCoords=True,
                          noPlots=noPlots)

            #shelvesAV[model][name][depthLevel] = m.shelvesAV
            shelvesAV.extend(m.shelvesAV)

            #####
            # no plots doesn't produce any plots, but does produce the list of shelves which can be used in Taylor/Target/Pattern diagrams.
            if noPlots: continue

            if noTargets: continue
            #csvFile = ukp.folder(workingDir+'/CSV')+'summary_file.csv'
            #print "attempting csvFromShelves:",m.shelves, csvFile
            #c = csvFromShelves.csvFromShelves(m.shelves, csvFile ,['check',])

            #####
            # makeTargets:
            # Make a target diagram of all matches for this particular dataset.
            #filename = ukp.folder(imageFolder+'/Targets/'+year+'/AllSlices')+model+'-'+jobID+'_'+year+'_'+name+depthLevel+'.png'
            #t = makeTargets(	m.shelves,
            #			filename,
            #			#name=name,
            #			legendKeys = ['newSlice','ykey',],
            #			debug=True)
            #			#imageDir='',

            #####
            # Produce a set of pattern and a target plots for each of the groups here.
            if annual:
                groups = {
                    'Oceans': [],
                    'depthRanges': [],
                    'BGCVal': [],
                }
            else:
                groups = {
                    'Oceans': [],
                    'Months': [],
                    'Seasons': [],
                    'NorthHemisphereMonths': [],
                    'SouthHemisphereMonths': [],
                    'depthRanges': [],
                    'BGCVal': [],
                }
            for g in groups:
                groups[g] = ukp.reducesShelves(shelvesAV,
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
                filename = ukp.folder(
                    imageFolder + '/Targets/' + year + '/' + name +
                    depthLevel + '/' + g
                ) + model + '_' + jobID + '_' + year + '_' + name + depthLevel + '_' + g + '.png'
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
                for o in ['Oceans', 'Months', 'depthRanges', 'BGCVal']:
                    if g.find(o) >= 0: xkeys = o
                if xkeys == '':
                    print("Could no find x axis keys!", g, 'in',
                          ['Oceans', 'Months', 'BGCVal'])

                filenamebase = ukp.folder(
                    imageFolder + '/Patterns/' + year + '/' + name +
                    depthLevel + '/' + g
                ) + 'Months_' + model + '_' + jobID + '_' + year + '_' + name + depthLevel
                makePatternStatsPlots(
                    {
                        name: groups[g],
                    },  # {legend, shelves}
                    name + ' ' + g,  #xkeysname
                    slicesDict[xkeys],  #xkeysLabels=
                    filenamebase,  # filename base	
                    grid=grid,
                    gridFile=gridFile)
            if not annual:
                #####
                # After finding all the shelves, we can plot them on the same axis.
                filenamebase = ukp.folder(
                    imageFolder + '/Patterns/' + year + '/' + name +
                    depthLevel + '/ANSH'
                ) + 'ANSH-Months_' + model + '_' + jobID + '_' + year + '_' + name + depthLevel

                makePatternStatsPlots(
                    {
                        'North Hemisphere': groups['NorthHemisphereMonths'],
                        'South Hemisphere': groups['SouthHemisphereMonths'],
                        'Global': groups['Months'],
                    },  # {legend, shelves}
                    name + ' Months',  #xkeysname
                    slicesDict['Months'],  #xkeysLabels=
                    filenamebase,  # filename base	
                    grid=grid,
                    gridFile=gridFile)

        if noPlots: continue
        if noTargets: continue
        #####
        # And now by depth levels:
        if annual: groups = [
                'Oceans',
                'depthRanges',
                'BGCVal',
        ]
        else:
            groups = [
                'Oceans',
                'Months',
                'Seasons',
                'depthRanges',
                'BGCVal',
            ]  #'NorthHemisphereMonths':[],'SouthHemisphereMonths':[]}
        for g in groups:
            if len(av[name]['layers']) <= 1: continue
            outShelves = {}
            for dl in av[name]['layers']:
                outShelves[dl] = ukp.reducesShelves(shelvesAV,
                                                    models=[
                                                        model,
                                                    ],
                                                    depthLevels=[
                                                        dl,
                                                    ],
                                                    names=[
                                                        name,
                                                    ],
                                                    sliceslist=slicesDict[g])
            filenamebase = ukp.folder(
                imageFolder + '/Patterns/' + year + '/' + name + 'AllDepths/'
            ) + 'AllDepths_' + g + '_' + model + '_' + jobID + '_' + year + '_' + name
            makePatternStatsPlots(outShelves,
                                  name + ' ' + g,
                                  slicesDict[g],
                                  filenamebase,
                                  grid=grid,
                                  gridFile=gridFile)

    return shelvesAV


if __name__ == "__main__":
    print('The end.')
