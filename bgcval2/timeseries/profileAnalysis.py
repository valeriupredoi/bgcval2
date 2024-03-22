#!/usr/bin/ipython

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
"""
.. module:: profileAnalysis
   :platform: Unix
   :synopsis: A tool for running a depth-profile time series analysis.
.. moduleauthor:: Lee de Mora <ledm@pml.ac.uk>
"""

import numpy as np
from shelve import open as shOpen
from netCDF4 import num2date
import os
import shutil

#Specific local code:
from bgcval2.bgcvaltools import bv2tools as bvt
from . import timeseriesTools as tst
from . import timeseriesPlots as tsp
from ..bgcvaltools.makeEORCAmasks import makeMaskNC
from ..bgcvaltools.pftnames import getLongName
from ..bgcvaltools.dataset import dataset


class profileAnalysis:

    def __init__(
        self,
        modelFiles,
        dataFile,
        dataType='',
        modelcoords='',
        modeldetails='',
        datacoords='',
        datadetails='',
        datasource='',
        model='',
        jobID='',
        layers='',
        regions='',
        metrics='',
        workingDir='',
        imageDir='',
        grid='',
        gridFile='',
        clean=True,
        debug=True,
    ):

        #####
        #	This is the class that does most of the legwork.
        #	First we save all the initialisation settings as class attributes.

        if debug: print("profileAnalysis:\t init.")
        self.modelFiles = modelFiles
        self.dataFile = dataFile
        self.dataType = dataType
        self.modelcoords = modelcoords
        self.modeldetails = modeldetails
        self.datacoords = datacoords
        self.datadetails = datadetails
        self.datasource = datasource
        self.model = model
        self.jobID = jobID
        self.layers = layers
        self.regions = regions
        self.metrics = metrics
        self.grid = grid
        self.gridFile = gridFile
        self.workingDir = workingDir
        self.imageDir = imageDir
        self.debug = debug
        self.clean = clean

        self.gridmaskshelve = bvt.folder(self.workingDir) + '_'.join([
            self.grid,
        ]) + '_masks.shelve'
        self.shelvefn = bvt.folder(self.workingDir) + '_'.join([
            'profile',
            self.jobID,
            self.dataType,
        ]) + '.shelve'
        self.shelvefn_insitu = bvt.folder(self.workingDir) + '_'.join([
            'profile',
            self.jobID,
            self.dataType,
        ]) + '_insitu.shelve'

        self._masksLoaded_ = False
        self.doHov = False
        #####
        # Load Data file
        self.__madeDataArea__ = False
        self.loadData()
        #assert 0

        #####
        # Load Model File
        self.loadModel()

        #####
        # Make the plots:
        self.makePlots()

        if self.debug:
            print("profileAnalysis:\tsafely finished ", self.dataType, (
                self.modeldetails['name']))

    def setmlayers(self):
        """	From the first model netcdf,
  		determine the number of depth layers.
  	"""
        if self.layers not in [
                'All',
                'Every2',
                'Every5',
                'Every10',
        ]:
            self.mlayers = self.layers
            return
        mlayers = self.calclayers(self.modelFiles[0], self.modelcoords['z'])
        if self.layers == 'All': self.mlayers = mlayers
        if self.layers == 'Every2': self.mlayers = mlayers[::2]
        if self.layers == 'Every5': self.mlayers = mlayers[::5]
        if self.layers == 'Every10': self.mlayers = mlayers[::10]
        print(self.mlayers)

    def setdlayers(self):
        """	From the data netcdf,
  		determine the number of depth layers.
  	"""
        if self.layers not in [
                'All',
                'Every2',
                'Every5',
                'Every10',
        ]:
            self.dlayers = self.layers
            return
        dlayers = self.calclayers(self.dataFile, self.datacoords['z'])
        if self.layers == 'All': self.dlayers = dlayers
        if self.layers == 'Every2': self.dlayers = dlayers[::2]
        if self.layers == 'Every5': self.dlayers = dlayers[::5]
        if self.layers == 'Every10': self.dlayers = dlayers[::10]

    def calclayers(self, fn, depthkey):

        nc = dataset(fn, 'r')
        depth = nc.variables[depthkey]
        if depth.ndim == 1:
            layers = np.arange(len(depth))
        if depth.ndim == 3:
            layers = np.arange(len(depth[:, 0, 0]))
        nc.close()
        return list(layers)

    def loadModel(self):
        if self.debug: print("profileAnalysis:\tloadModel.")
        ####
        # load and calculate the model info
        self.setmlayers()

        try:
            if self.clean:
                print("profileAnalysis:\tloadModel:\tUser requested clean run. Wiping old data.")
                assert 0
            with shOpen(self.shelvefn) as sh:
                readFiles = sh['readFiles']
                modeldataD = sh['modeldata']
            # sh = shOpen(self.shelvefn)
            # readFiles = sh['readFiles']
            # modeldataD = sh['modeldata']
            # sh.close()
            print("OprofileAnalysis:\tloadModel:\tpened shelve:", self.shelvefn, '\tread', len(
                readFiles))
        except:
            readFiles = []
            modeldataD = {}
            #		self.setmlayers()
            for r in self.regions:
                for l in self.mlayers:
                    for m in self.metrics:
                        modeldataD[(r, l, m)] = {}

            print("profileAnalysis:\tloadModel:\tCould not open shelve:", self.shelvefn, '\tread', len(
                readFiles))

        ###############
        # Check whethere there has been a change in what was requested:
        for r in self.regions:
            for l in self.mlayers:
                for m in self.metrics:
                    if self.debug:
                        print("profileAnalysis:\tloadModel:\tChecking: ", [
                            r,
                            l,
                            m,
                        ], '\t', end=' ')
                    try:
                        if self.debug:
                            print('has ', len(list(modeldataD[(r, l,
                                                          m)].keys())), 'keys')
                    except:
                        readFiles = []
                        modeldataD[(r, l, m)] = {}
                        if self.debug: print('has no keys')
                    try:
                        if len(list(modeldataD[(r, l, m)].keys())) == 0:
                            print("profileAnalysis:\tloadModel:\tmodeldataD[", (
                                r, l, m), "] has no keys")
                            readFiles = []
                            assert 0

                    except:
                        pass

        #####
        # Summarise checks
        if self.debug:
            print("profileAnalysis:\tloadModel:\tloadModel:post checks:")
            #print "modeldataD:",modeldataD
            print("profileAnalysis:\tloadModel:\tshelveFn:", self.shelvefn)
            print("profileAnalysis:\tloadModel:\treadFiles:", end=' ')
            try:
                print(readFiles[-1])
            except:
                print('...')

        ###############
        # Load files, and calculate fields.
        openedFiles = 0
        for fn in self.modelFiles:
            if fn in readFiles: continue

            if not self._masksLoaded_:
                self.loadMasks()

            print("profileAnalysis:\tloadModel:\tloading new file:", self.dataType, fn, end=' ')
            nc = dataset(fn, 'r')
            ts = tst.getTimes(nc, self.modelcoords)
            meantime = np.mean(ts)
            print("\ttime:", meantime)

            #DL = tst.DataLoader(fn,nc,self.modelcoords,self.modeldetails, regions = self.regions, layers = self.layers,)
            nc = dataset(fn, 'r')
            dataAll = bvt.extractData(nc, self.modeldetails).squeeze()

            for r in self.regions:
                for m in self.metrics:
                    if m == 'mean':
                        data = bvt.mameanaxis(np.ma.masked_where(
                            (self.modelMasks[r] != 1) + dataAll.mask, dataAll),
                                              axis=(1, 2))

                        if self.debug:
                            print("profileAnalysis:\tloadModel.", r, m, self.dataType, '\tyear:', int(
                                meantime), 'mean:', data.mean())
                        #if self.debug:print "profileAnalysis:\tloadModel.",self.dataType, data.shape, data.min(),data.max(), dataAll.shape ,self.modelMasks[r].shape, dataAll.min(),dataAll.max()

                        alllayers = []
                        for l, d in enumerate(data):
                            #print "Saving model data profile",r,m,l,d
                            modeldataD[(r, l, m)][meantime] = d
                            alllayers.append(l)

                        #####
                        # Add a masked value in layers where there is no data.

                        for l in self.mlayers:
                            if l in alllayers: continue
                            modeldataD[(r, l, m)][meantime] = np.ma.masked

                    else:
                        print('ERROR:', m, "not implemented in profile")
                        assert 0

            readFiles.append(fn)
            openedFiles += 1

            nc.close()
            if openedFiles:
                print("Saving shelve:", self.dataType, self.shelvefn, '\tread', len(
                    readFiles))
                with shOpen(self.shelvefn) as sh:
                    sh['readFiles'] = readFiles
                    sh['modeldata'] = modeldataD
                # sh = shOpen(self.shelvefn)
                # sh['readFiles'] = readFiles
                # sh['modeldata'] = modeldataD
                # sh.close()                
                openedFiles = 0
        if openedFiles:
            print("Saving shelve:", self.dataType, self.shelvefn, '\tread', len(
                readFiles))
            with shOpen(self.shelvefn) as sh:
                sh['readFiles'] = readFiles
                sh['modeldata'] = modeldataD

            # sh = shOpen(self.shelvefn)
            # sh['readFiles'] = readFiles
            # sh['modeldata'] = modeldataD
            # sh.close()
        self.modeldataD = modeldataD
        if self.debug:
            print("profileAnalysis:\tloadModel.\t Model loaded:", list(self.modeldataD.keys(
            ))[:3], '...', len(list(self.modeldataD.keys())))

    def loadMasks(self):
        #####
        # Here we load the masks file.
        self.maskfn = 'data/' + self.grid + '_masks.nc'

        if not os.path.exists(self.maskfn):
            print("Making mask file", self.maskfn)

            makeMaskNC(self.maskfn, self.regions, self.grid)

        self.modelMasks = {}

        ncmasks = dataset(self.maskfn, 'r')

        for r in self.regions:
            if r in list(ncmasks.variables.keys()):
                print("Loading mask", r)
                self.modelMasks[r] = ncmasks.variables[r][:]

            else:
                newmask = 'data/' + self.grid + '_masks_' + r + '.nc'

                if not os.path.exists(newmask):
                    makeMaskNC(newmask, [
                        r,
                    ], self.grid)
                nc = dataset(newmask, 'r')
                self.modelMasks[r] = nc.variables[r][:]
                nc.close()

        print("Loaded masks", list(self.modelMasks.keys()))

        ncmasks.close()
        self._masksLoaded_ = True

    def loadData(self):

        if self.debug: print("profileAnalysis:\t loadData.", self.dataFile)

        if not self.dataFile:
            if self.debug:
                print("profileAnalysis:\t No data File provided:", self.dataFile)
            self.dataD = {}
            return

        if not os.path.exists(self.dataFile):
            if self.debug:
                print("profileAnalysis:\tWARNING:\t No such data File:", self.dataFile)
            self.dataD = {}
            return

        ###############
        # load and calculate the real data info
        try:
            if self.clean:
                print("profileAnalysis:\t loadData\tUser requested clean run. Wiping old data.")
                assert 0
            with shOpen(self.shelvefn_insitu) as sh:
                dataD = sh['dataD']
            # sh = shOpen(self.shelvefn_insitu)
            # dataD = sh['dataD']
            # sh.close()            
            print("profileAnalysis:\t loadData\tOpened shelve:", self.shelvefn_insitu)
            self.dataD = dataD
        except:
            dataD = {}
            print("profileAnalysis:\t loadData\tCould not open shelve:", self.shelvefn_insitu)

        ###############
        # Test to find out if we need to load the netcdf, or if we can just return the dict as a self.object.
        needtoLoad = False
        self.setdlayers()
        for r in self.regions:
            #if needtoLoad:continue
            #    for l in self.dlayers:
            for l in sorted(self.dlayers)[:]:
                #if needtoLoad:continue
                try:
                    dat = self.dataD[(r, l)]
                    #test = (len(),self.dataD[(r,l)].shape)
                    print("profileAnalysis:\t loadData\t", (r, l))  #,dat
                except:
                    needtoLoad = True
                    print("profileAnalysis:\t loadData\tUnable to load", (r, l))

        if needtoLoad: pass
        else:
            self.dataD = dataD
            return

        ###############
        # Loading data for each region.
        print("profileAnalysis:\t loadData,\tloading ", self.dataFile)
        nc = dataset(self.dataFile, 'r')
        data = tst.loadData(nc, self.datadetails)

        if not self.__madeDataArea__: self.AddDataArea()

        ###############
        # Loading data for each region.
        dl = tst.DataLoader(
            self.dataFile,
            '',
            self.datacoords,
            self.datadetails,
            regions=self.regions,
            layers=self.dlayers[:],
        )

        #	#for r in self.regions:
        #	 #   for l in self.dlayers:
        #	    	dataD[(r,l)] = dl.load[(r,l,)]
        #	    	dataD[(r,l,'lat')] = dl.load[(r,l,'lat')]
        #	    	dataD[(r,l,'lon')] = dl.load[(r,l,'lon')]
        #		if len(dataD[(r,l)])==0  or np.ma.is_masked(dataD[(r,l)]):
        #			dataD[(r,l)]  = np.ma.masked
        #			dataD[(r,l,'lat')]  = np.ma.masked
        #			dataD[(r,l,'lon')]  = np.ma.masked

        AreaNeeded = len(
            bvt.intersection([
                'mean',
                'median',
                'sum',
            ], self.metrics))

        maskedValue = np.ma.masked  # -999.# np.ma.array([-999.,],mask=[True,])
        #maskedValue = np.ma.array([-999.,],mask=[True,])
        #maskedValue  = -999 #np.ma.array([-999.,],mask=[True,])

        for l in sorted(self.dlayers)[:]:
            for r in self.regions:
                dataD[(r, l)] = dl.load[(
                    r,
                    l,
                )]
                try:
                    meandatad = dataD[(r, l)].mean()
                    datadmask = (~np.ma.array(dataD[(r, l)]).mask).sum()
                except:
                    meandatad = False
                    datadmask = False

                if np.isnan(meandatad) or np.isinf(meandatad) or dataD[(
                        r, l)].mask.all() or np.ma.is_masked(meandatad):
                    meandatad = False
                    datadmask = False

    #print "profileAnalysis:\t load in situ data,\tloaded ",(r,l),  'mean:',meandatad

                if False in [meandatad, datadmask]:
                    dataD[(r, l)] = maskedValue
                    dataD[(r, l, 'lat')] = maskedValue
                    dataD[(r, l, 'lon')] = maskedValue
                    if AreaNeeded: dataD[(r, l, 'area')] = maskedValue
                    print("profileAnalysis:\t loadData\tproblem with ", (
                        r, l), 'data:\t', dataD[(r, l)])

                else:

                    dataD[(r, l, 'lat')] = dl.load[(r, l, 'lat')]
                    dataD[(r, l, 'lon')] = dl.load[(r, l, 'lon')]
                    if AreaNeeded:
                        dataD[(r, l, 'area')] = self.loadDataAreas(
                            dataD[(r, l, 'lat')], dataD[(r, l, 'lon')])
                    else:
                        dataD[(r, l, 'area')] = np.ones_like(dataD[(r, l,
                                                                    'lon')])
                    print("profileAnalysis:\t loadData,\tloading ", (
                        r, l), 'mean:\t', meandatad)

            print("profileAnalysis:\t loadData.\tSaving shelve: (layer", l, ")", self.shelvefn_insitu)
            with shOpen(self.shelvefn_insitu) as sh:
                sh['dataD'] = dataD
            # sh = shOpen(self.shelvefn_insitu)
            # sh['dataD'] = dataD
            # sh.close()
        ###############
        # Savng shelve
        print("profileAnalysis:\t loadData.\tSaving shelve:", self.shelvefn_insitu)
        try:
            with shOpen(self.shelvefn_insitu) as sh:
                sh['dataD'] = dataD
            # sh = shOpen(self.shelvefn_insitu)
            # sh['dataD'] = dataD
            # sh.close()
            print("profileAnalysis:\t loadData.\tSaved shelve:", self.shelvefn_insitu)

        except:
            print("profileAnalysis:\t WARNING.\tSaving shelve failed, trying again.:", self.shelvefn_insitu)
            print("Data is", list(dataD.keys()))

            for key in sorted(dataD.keys()):

                print(key, ':\t', dataD[key])
                with shOpen(bvt.folder('./tmpshelves') + 'tmshelve.shelve') as sh:
                    sh['dataD'] = dataD[key]
                # sh = shOpen(bvt.folder('./tmpshelves') + 'tmshelve.shelve')
                # sh['dataD'] = dataD[key]
                # sh.close()                
                print("saved fine:\t./tmpshelves/tmshelve.shelve")

            shutil.move(self.shelvefn_insitu, self.shelvefn_insitu + '.broken')

            with  shOpen(self.shelvefn_insitu) as sh:
                sh['dataD'] = dataD

#		except:
#			print "profileAnalysis:\t WARNING.\tUnable to Save in situ shelve.\tYou'll have to input it each time.",self.shelvefn_insitu

        self.dataD = dataD

    def AddDataArea(self, ):
        """
  	Adding Area dictionany
  	"""
        if not self.dataFile:
            self.dataAreaDict = {}
            return
        area = tst.makeArea(self.dataFile, self.datacoords)
        nc = dataset(self.dataFile, 'r')
        lats = nc.variables[self.datacoords['lat']][:]
        lons = nc.variables[self.datacoords['lon']][:]
        nc.close()
        print("timeseriesAnalysis:\tAddDataArea:\t", area.shape, lats.shape, lons.shape)
        self.dataAreaDict = {}
        if lats.ndim == 2:
            for (i, j), a in np.ndenumerate(area):
                #if np.ma.is_masked(a):continue
                self.dataAreaDict[(lats[i, j], lons[i, j])] = a
        if lats.ndim == 1:
            for (i, j), a in np.ndenumerate(area):
                #if np.ma.is_masked(a):continue
                self.dataAreaDict[(lats[i], lons[j])] = a
        self.__madeDataArea__ = True

    def loadDataAreas(self, lats, lons):
        """
  	Adding Area for each region.
  	"""
        areas = []
        for la, lo in zip(lats, lons):
            try:
                areas.append(self.dataAreaDict[(la, lo)])
            except:
                areas.append(0.)
        return np.ma.array(areas)

    def mapplotsRegionsLayers(self, ):
        """	Makes a map plot of model vs data for each string-named layer (not numbered layers). 
  	"""
        newlayers = [
            l for l in self.layers
            if type(l) not in [type(0), type(0.)]
        ]
        mDL = tst.DataLoader(
            self.modelFiles[-1],
            '',
            self.modelcoords,
            self.modeldetails,
            regions=self.regions,
            layers=newlayers,
        )
        for r in self.regions:
            for l in self.layers:
                if type(l) in [type(0), type(0.)]: continue
                mapfilename = bvt.folder(self.imageDir + '/' +
                                         self.dataType) + '_'.join([
                                             'map',
                                             self.jobID,
                                             self.dataType,
                                             str(l),
                                             r,
                                         ]) + '.png'
                modeldata = mDL.load[(r, l)]
                modellat = mDL.load[(r, l, 'lat')]
                modellon = mDL.load[(r, l, 'lon')]

                if not len(modeldata): continue

                print("mapplotsRegionsLayers:\t", r, l, "model contains", len(
                    modeldata), 'model data')
                print("mapplotsRegionsLayers:\t", r, l, "model lat:", modellat.min(
                ), modellat.mean(), modellat.max())
                print("mapplotsRegionsLayers:\t", r, l, "model lon:", modellon.min(
                ), modellon.mean(), modellon.max())

                if self.dataFile:
                    datadata = self.dataD[(r, l)]
                    datalat = self.dataD[(r, l, 'lat')]
                    datalon = self.dataD[(r, l, 'lon')]

                else:
                    datadata = np.ma.array([
                        -1000,
                    ], mask=[
                        True,
                    ])
                    datalat = np.ma.array([
                        -1000,
                    ], mask=[
                        True,
                    ])
                    datalon = np.ma.array([
                        -1000,
                    ], mask=[
                        True,
                    ])

                print("mapplotsRegionsLayers:\t", r, l, "contains", len(
                    datadata), 'in situ data')
                print("mapplotsRegionsLayers:\t", r, l, "data lat:", len(
                    datalat), datalat.min(), datalat.mean(), datalat.max())
                print("mapplotsRegionsLayers:\t", r, l, "data lon:", len(
                    datalon), datalon.min(), datalon.mean(), datalon.max())

                titles = [
                    ' '.join([
                        getLongName(t) for t in [
                            self.model, '(' + self.jobID + ')',
                            str(l), self.modeldetails['name']
                        ]
                    ]), ' '.join([
                        getLongName(t) for t in
                        [self.datasource,
                         str(l), self.datadetails['name']]
                    ])
                ]

                tsp.mapPlotPair(
                    modellon,
                    modellat,
                    modeldata,
                    datalon,
                    datalat,
                    datadata,
                    mapfilename,
                    titles=titles,
                    lon0=0.,
                    drawCbar=True,
                    cbarlabel='',
                    dpi=100,
                )

    def makePlots(self):
        if self.debug: print("profileAnalysis:\t makePlots.")

        #####
        # create a dictionary of model and data depths and layers.
        mnc = dataset(self.modelFiles[-1], 'r')
        modelZcoords = {
            i: z
            for i, z in enumerate(mnc.variables[self.modelcoords['z']][:])
        }
        mnc.close()

        if self.dataFile:
            dnc = dataset(self.dataFile, 'r')
            if self.debug:
                print("profileAnalysis:\t makePlots\tOpening", self.dataFile)
            dataZcoords = {
                i: z
                for i, z in enumerate(dnc.variables[self.datacoords['z']][:])
            }
            if self.debug:
                print("profileAnalysis:\t makePlots\tOpened", self.dataFile)

            if self.debug:
                print("profileAnalysis:\t makePlots\tloaded", dataZcoords)
            dnc.close()
        else:
            dataZcoords = {}

        #####
        # Hovmoeller plots
        for r in self.regions:
            for m in self.metrics:
                if m not in [
                        'mean',
                        'median',
                        'min',
                        'max',
                ]: continue
                if self.debug: print("profileAnalysis:\t makePlots\t", r, m)
                #####
                # Load data layers:
                data = {}
                if self.dataFile:
                    for l in self.dlayers:
                        #print "Hovmoeller plots:",r,m,l

                        if type(l) == type('str'):
                            continue  # no strings, only numbered layers.
                        if l > max(dataZcoords.keys()): continue
                        #####
                        # Test for presence/absence of in situ data.
                        try:
                            dataslice = self.dataD[(r, l)]
                            dataslice = dataslice.compressed()
                        except:
                            dataslice = np.ma.array([
                                -1000,
                            ], mask=[
                                True,
                            ])

                        if m == 'mean':
                            try:
                                data[l] = np.ma.mean(dataslice)
                            except:
                                data[l] = np.ma.array([
                                    -1000,
                                ], mask=[
                                    True,
                                ])
                        elif m == 'median':
                            try:
                                data[l] = np.ma.median(dataslice)
                            except:
                                data[l] = np.ma.array([
                                    -1000,
                                ], mask=[
                                    True,
                                ])
                        elif m == 'min':
                            try:
                                data[l] = np.ma.min(dataslice)
                            except:
                                data[l] = np.ma.array([
                                    -1000,
                                ], mask=[
                                    True,
                                ])
                        elif m == 'max':
                            try:
                                data[l] = np.ma.max(dataslice)
                            except:
                                data[l] = np.ma.array([
                                    -1000,
                                ], mask=[
                                    True,
                                ])

        #		if self.debug: print "profileAnalysis:\tmakePlots:\tHovmoeller plots:",r,m,l,'\tdata'#,data[l]

        #####
        # Load model layers:
                modeldata = {}
                for l in self.mlayers:
                    if type(l) == type('str'):
                        continue  # no strings, only numbered layers.
                    if l > max(modelZcoords.keys()): continue
                    modeldata[l] = self.modeldataD[(r, l, m)]
                if self.debug:
                    print("profileAnalysis:\tmakePlots:\tHovmoeller plots:", r, m, '\tloaded model data')

                #####
                # check that multiple layers were requested.
                #if len(data.keys())<1: continue
                if len(list(modeldata.keys())) < 1: continue

                title = ' '.join(
                    [getLongName(t) for t in [r, m, self.dataType]])
                profilefn = bvt.folder(self.imageDir + '/' +
                                       self.dataType) + '_'.join([
                                           'profile',
                                           self.jobID,
                                           self.dataType,
                                           r,
                                           m,
                                       ]) + '.png'
                axislabel = getLongName(
                    self.modeldetails['name']) + ', ' + getLongName(
                        self.modeldetails['units'])
                if bvt.shouldIMakeFile([self.shelvefn, self.shelvefn_insitu],
                                       profilefn,
                                       debug=False):
                    tsp.profilePlot(
                        modeldata,
                        data,
                        profilefn,
                        modelZcoords=modelZcoords,
                        dataZcoords=dataZcoords,
                        xaxislabel=axislabel,
                        title=title,
                    )

                if self.doHov:
                    hovfilename = bvt.folder(self.imageDir + '/' +
                                             self.dataType) + '_'.join([
                                                 'profilehov',
                                                 self.jobID,
                                                 self.dataType,
                                                 r,
                                                 m,
                                             ]) + '.png'
                    if bvt.shouldIMakeFile(
                        [self.shelvefn, self.shelvefn_insitu],
                            hovfilename,
                            debug=False):
                        tsp.hovmoellerPlot(modeldata,
                                           data,
                                           hovfilename,
                                           modelZcoords=modelZcoords,
                                           dataZcoords=dataZcoords,
                                           title=title,
                                           zaxislabel=axislabel,
                                           diff=False)

                    hovfilename_diff = bvt.folder(self.imageDir + '/' +
                                                  self.dataType) + '_'.join([
                                                      'profileDiff',
                                                      self.jobID,
                                                      self.dataType,
                                                      r,
                                                      m,
                                                  ]) + '.png'
                    if bvt.shouldIMakeFile(
                        [self.shelvefn, self.shelvefn_insitu],
                            hovfilename_diff,
                            debug=False):
                        tsp.hovmoellerPlot(modeldata,
                                           data,
                                           hovfilename_diff,
                                           modelZcoords=modelZcoords,
                                           dataZcoords=dataZcoords,
                                           title=title,
                                           zaxislabel=axislabel,
                                           diff=True)
