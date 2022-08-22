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
.. module:: timeseriesAnalysis
   :platform: Unix
   :synopsis: A tool for running a time series analysis.
.. moduleauthor:: Lee de Mora <ledm@pml.ac.uk>
"""

import numpy as np
from shelve import open as shOpen
from netCDF4 import num2date
import os
import shutil
import glob
import errno


#Specific local code:
from bgcval2 import UKESMpython as ukp
from bgcval2.bgcvaltools.pftnames import getLongName
from bgcval2.bgcvaltools.dataset import dataset
from bgcval2.timeseries import timeseriesTools as tst
from bgcval2.timeseries import timeseriesPlots as tsp

#getTimes, loadData


class timeseriesAnalysis:

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
            clean=False,
            debug=True,
            noNewFiles=False,  # stops loading new files
    ):

        #####
        #	This is the class that does most of the legwork.
        #	First we save all the initialisation settings as class attributes.

        if debug: print("timeseriesAnalysis:\t init.")
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
        self.noNewFiles = noNewFiles


        # workingdir set else to match       shelvedir = ukp.folder(paths.shelvedir + "/timeseries/" + jobID)

        self.shelvefn = self.workingDir + '_'.join([
            self.jobID,
            self.dataType,
        ]) + '.shelve'
        self.shelvefn_insitu = self.workingDir + '_'.join([
            self.jobID,
            self.dataType,
        ]) + '_insitu.shelve'

        #####
        # Load Data file
        self.__madeDataArea__ = False
        if self.noNewFiles: pass
        else: self.loadData()

        #####
        # Load Model File
        self.loadModelWeightsDict()
        if 'wcvweighted' in self.metrics: 
            self.loadModelwcvDict()
        self.loadModel()
        #assert 0

        #####
        # return Model data without making new images
        if self.noNewFiles: 
            return

        #####
        # Make the plots:
        self.makePlots()

        if self.debug:
            print("timeseriesAnalysis:\tsafely finished ", self.dataType,
                  (self.modeldetails['name']))

    def loadModel(self):
        if self.debug: print("timeseriesAnalysis:\tloadModel.")
        ####
        # load and calculate the model info
        #if os.path.exists(self.shelvefn):
        #    sh = shOpen(self.shelvefn)
        #    print('seems fine:', self.shelvefn)
        #else:
        #    print('Does not exist', self.shelvefn)
        #    sh = shOpen(self.shelvefn)
        #    print (sh.keys())
        #    readFiles       = sh['readFiles']
        #    modeldataD      = sh['modeldata']
        #    print(readFiles)
        readFiles = []
        modeldataD = {}
        for r in self.regions:
            for l in self.layers:
                for m in self.metrics:
                    modeldataD[(r, l, m)] = {}

        if self.clean:
            print(
                "timeseriesAnalysis:\tloadModel\tUser requested clean run. Wiping old data."
            )
        elif glob.glob(''.join([self.shelvefn, '*'])):
            print("timeseriesAnalysis:\tloadModel\tOpening shelve:",
                   self.shelvefn)

            # explicitly open in read only:
            with shOpen(self.shelvefn, flag='r', writeback=False, protocol=5,) as sh:
                print(sh.keys())
#                if 'readFiles' in sh.keys() and 'modeldata' in sh.keys():
                if sh.get('readFiles', False) and sh.get('modeldata', False):
                     readFiles = sh['readFiles']
                     modeldataD = sh['modeldata']

            if not len(readFiles):
                print("timeseriesAnalysis:\tloadModel\tOpened shelve - but it's empty:",
                      self.shelvefn, '\tread', len(readFiles))

        else:          
            print("timeseriesAnalysis:\tloadModel\tCould not open model shelve:",
                  self.shelvefn, '\tread', len(readFiles))
        #exit(0)
        #assert 0 
        ###############
        # Check whether there has been a change in what was requested:
        for r in self.regions:
            for l in self.layers:
                for m in self.metrics:
                    if self.debug:
                        print("timeseriesAnalysis:\tloadModel\tChecking: ", [
                            r,
                            l,
                            m,
                        ],
                              '\t',
                              end=' ')
                    try:
                        if self.debug:
                            print('has ',
                                  len(list(modeldataD[(r, l, m)].keys())),
                                  'keys')
                    except:
                        readFiles = []
                        modeldataD[(r, l, m)] = {}
                        if self.debug: print('has no keys')
                    try:
                        if len(list(modeldataD[(r, l, m)].keys())) == 0:
                            readFiles = []
                    except:
                        pass
        #####
        # Check if the Input file has changed since the shelve file last changed.
        reDoFiles = []
        #for fn in sorted(readFiles):
        #        #if self.debug:print "timeseriesAnalysis:\tloadModel\tChecking: ",fn

        #	if ukp.shouldIMakeFile(fn, self.shelvefn+'*',debug=False):
        #	#		print "timeseriesAnalysis:\tloadModel\t:this file should be re-analysed: because shelve is older thna input file.", fn
        #		readFiles.remove(fn)
        #	        reDoFiles.append(fn)

        #####
        # Check if the Input file has changed since the shelve file last changed.
        #for fn in sorted(readFiles):
        #	if ukp.shouldIMakeFile(fn, self.shelvefn+'*',debug=False):
        #		print "timeseriesAnalysis:\tloadModel\t:this file should be re-analysed: because input file is newer)", fn
        #		readFiles.remove(fn)

        #####
        # Summarise checks
        if self.debug:
            print("timeseriesAnalysis:\tloadModel:\tpost checks...")
            #print "modeldataD:",modeldataD
            print("timeseriesAnalysis:\tloadModel\tshelveFn:", self.shelvefn)
            print("timeseriesAnalysis:\tloadModel\treadFiles: contains ",
                  len(readFiles),
                  end=' ')
            try:
                print("files.\tUp to ", sorted(readFiles)[-1])
            except:
                print("files.")

        #####
        # No New Files checks - to save time and avoid double work.
        if self.noNewFiles:
            self.modeldataD = modeldataD
            if self.debug:
                print(
                    "timeseriesAnalysis:\tloadModel.\tno New Files requested. Loaded: ",
                    len(list(modeldataD.keys())), 'Model data')
            return

        percentiles = {}
        for m in self.metrics:
            if m.find('pc') > -1:
                pc = float(m.replace('pc', ''))
                percentiles[pc] = True
            if m == 'median': percentiles[50.] = True
        percentiles = sorted(percentiles.keys())

        ###############
        # Load files, and calculate fields.
        openedFiles = 0
        save_every = 10
        for fn in self.modelFiles:
            if fn in readFiles: continue
            print("timeseriesAnalysis:\tloadModel:\tloading new file:",
                  fn,
                  end=' ')
            nc = dataset(fn, 'r')
            ts = tst.getTimes(nc, self.modelcoords)
            meantime = np.mean(ts)
            print("\ttime:", meantime)

            DL = tst.DataLoader(
                fn,
                nc,
                self.modelcoords,
                self.modeldetails,
                regions=self.regions,
                layers=self.layers,
            )

            for l in self.layers:
                for r in self.regions:

                    #####
                    # Check wherether you can skip loading this metric,region,layer
                    skip = True
                    for m in self.metrics:
                        if skip == False: continue
                        try:
                            a = modeldataD[(r, l, m)][meantime]
                            print(
                                "timeseriesAnalysis:\tloadModel\tAlready created ",
                                int(meantime), ':\t', (r, l, m), '\t=', a)
                        except:
                            skip = False
                            print(
                                "timeseriesAnalysis:\tloadModel\tNeed to create:\t",
                                int(meantime), '\t', (r, l, m))
                    if fn in reDoFiles:
                        print(
                            "timeseriesAnalysis:\tloadModel\tNeed to re-load ",
                            int(meantime), ':\t', (r, l, m))
                        skip = False
                    if skip: continue

                    #####
                    # can't skip it, need to load it.
                    layerdata = DL.load[(r, l)]
                    #####
                    # get Weights:
                    volumeWeightedLayers = ['All', 'Transect']

                    if len(
                            ukp.intersection([
                                'mean',
                                'median',
                                'sum',
                            ], self.metrics)):
                        lats = DL.load[(r, l, 'lat')]
                        lons = DL.load[(r, l, 'lon')]
                        if l in volumeWeightedLayers:
                            depths = DL.load[(r, l, 'z')]
                            weights = np.array([
                                self.weightsDict[(la, lo, z)]
                                for la, lo in zip(lats, lons, depths)
                            ])
                        else:
                            #weights = np.array([self.weightsDict[(la,lo)] for la,lo in zip(lats,lons)])
                            weights = []
                            for la, lo, da in zip(lats, lons, layerdata):
                                try:
                                    weights.append(self.weightsDict[(la, lo)])
                                except:
                                    #print "timeseriesAnalysis:\tloadModel\tunable to load area",la,lo,da, 'adding 0. weight'
                                    weights.append(0.)

                    else:
                        weights = np.ones_like(layerdata)

                    if 'wcvweighted' in self.metrics:
                        lats = DL.load[(r, l, 'lat')]
                        lons = DL.load[(r, l, 'lon')]
                        wcvweights = []
                        for la, lo, da in zip(lats, lons, layerdata):
                            try:
                                wcvweights.append(self.wcvDict[(la, lo)])
                            except:
                                continue  #wcvweights.append(0.)
                        wcvweights = np.array(wcvweights)

                        print("Loaded Water Column Volume Weights",
                              wcvweights.min(), wcvweights.mean(),
                              wcvweights.max())

                    if type(layerdata) == type(
                            np.ma.array([
                                1,
                                -999,
                            ], mask=[
                                False,
                                True,
                            ])):
                        weights = np.ma.array(weights)
                        #print weights.mean(),weights.min(),weights.max()
                        weights = np.ma.masked_where(
                            (weights == 0.) + weights.mask + layerdata.mask,
                            weights)  #.compressed()
                        layerdata = np.ma.masked_where(
                            (weights == 0.) + weights.mask + layerdata.mask,
                            layerdata)  #.compressed()
                        weights = weights.compressed()
                        layerdata = layerdata.compressed()
                        if len(layerdata) != len(weights):
                            print("1.b len(	layerdata)!= len(weights)",
                                  len(layerdata), '!=', len(weights))
                            assert 0

                    if len(layerdata) == 0:
                        for m in self.metrics:
                            modeldataD[(
                                r, l, m
                            )][meantime] = np.ma.masked  #np.ma.array([-999,],mask=[True,])
                        continue


#		  	for m in self.metrics:
#		  		try:
#		  			a = modeldataD[(r,l,m)][meantime]
#		  			continue
#		  		except:pass
#				if m == 'mean':   	modeldataD[(r,l,m)][meantime] = np.ma.average(layerdata,weights=weights)
#				if m == 'median':   	modeldataD[(r,l,m)][meantime] = np.ma.median(layerdata)
#				if m == 'sum':   	modeldataD[(r,l,m)][meantime] = np.ma.sum(layerdata)
#				if m == 'metricless':  	modeldataD[(r,l,m)][meantime] = np.ma.sum(layerdata)	# same as sum
#				if m == 'min':   	modeldataD[(r,l,m)][meantime] = np.ma.min(layerdata)
#				if m == 'max':   	modeldataD[(r,l,m)][meantime] = np.ma.max(layerdata)
#				if m.find('pc')>-1:
#					pc = int(m.replace('pc',''))
#					modeldataD[(r,l,m)][meantime] = np.percentile(layerdata,pc)
#
#		  		print "timeseriesAnalysis:\tloadModel\tLoaded metric:", int(meantime),'\t',[(r,l,m)], '\t',modeldataD[(r,l,m)][meantime]

                    if 'mean' in self.metrics:
                        modeldataD[(r, l, 'mean')][meantime] = np.ma.average(
                            layerdata, weights=weights)
                    if 'sum' in self.metrics:
                        modeldataD[(r, l,
                                    'sum')][meantime] = np.ma.sum(layerdata)
                    if 'min' in self.metrics:
                        modeldataD[(r, l,
                                    'min')][meantime] = np.ma.min(layerdata)
                    if 'max' in self.metrics:
                        modeldataD[(r, l,
                                    'max')][meantime] = np.ma.max(layerdata)
                    if 'metricless' in self.metrics:
                        modeldataD[(
                            r, l,
                            'metricless')][meantime] = np.ma.sum(layerdata)

                    if 'wcvweighted' in self.metrics:
                        #print 'wcvweights', wcvweights.shape, wcvweights.min(), wcvweights.mean(),wcvweights.max()
                        #print 'layerdata:',layerdata.shape, layerdata.min(), layerdata.mean(),layerdata.max()
                        modeldataD[(r, l,
                                    'wcvweighted')][meantime] = np.ma.average(
                                        layerdata, weights=wcvweights)

                    if len(percentiles) == 0: continue
                    out_pc = ukp.weighted_percentiles(layerdata,
                                                      percentiles,
                                                      weights=weights)

                    for pc, dat in zip(percentiles, out_pc):
                        modeldataD[(r, l,
                                    ukp.mnStr(pc) + 'pc')][meantime] = dat
                        if pc == 50.:
                            modeldataD[(r, l, 'median')][meantime] = dat

                    print("timeseriesAnalysis:\tloadModel\tLoaded metric:\t",
                          int(meantime), '\t', [(r, l, m)], '\t',
                          modeldataD[(r, l, m)][meantime])

            readFiles.append(fn)
            openedFiles += 1

            nc.close()
            if openedFiles and openedFiles%save_every==0:
                print("timeseriesAnalysis:\tloadModel\tSaving shelve:",
                      self.shelvefn, '\tread', len(readFiles))
                sh = shOpen(self.shelvefn, protocol=5)
                sh['readFiles'] = readFiles
                sh['modeldata'] = modeldataD
                sh.close()
                openedFiles = 0
        if openedFiles:
            print("timeseriesAnalysis:\tloadModel\tSaving shelve - last time:",
                  self.shelvefn, '\tread', len(readFiles))
            sh = shOpen(self.shelvefn, protocol=5)
            sh['readFiles'] = readFiles
            sh['modeldata'] = modeldataD
            sh.close()

        self.modeldataD = modeldataD
        if self.debug:
            print("timeseriesAnalysis:\tloadModel.\t Model loaded:",
                  list(self.modeldataD.keys())[:3], '...',
                  len(list(self.modeldataD.keys())))

    def loadModelWeightsDict(self, ):
        """
  	Adding Area dictionany for Model.
  	"""
        if not os.path.exists(self.gridFile):
            raise FileNotFoundError(
                errno.ENOENT, os.strerror(errno.ENOENT), self.gridFile)
        nc = dataset(self.gridFile, 'r')
        tmask = nc.variables['tmask'][:]
        try:
            area = nc.variables['area'][:]
        except:
            area = nc.variables['e2t'][:] * nc.variables['e1t'][:]
        if tmask.ndim == 3: area = np.ma.masked_where(tmask[0] == 0, area)
        if tmask.ndim == 2: area = np.ma.masked_where(tmask == 0, area)

        #if 'layerless' in self.layers:
        #try:
        #	pvol  = nc.variables['pvol' ][:]
        #except:
        #	pvol = nc.variables['e3t'][:] *area
        #pvol = np.ma.masked_where(tmask==0,pvol)

        self.weightsDict = {}
        if self.modelcoords['lat'] == self.modelcoords['lon'] == False:
            ####
            # A scalar field with no lat or lon coordinates.
            self.weightsDict[(0, 0)] = 1.
            self.weightsDict[(False, False)] = 1.
            return

        lats = nc.variables[self.modelcoords['lat']][:]
        lons = nc.variables[self.modelcoords['lon']][:]
        nc.close()

        if lats.ndim == 2:
            for (i, j), a in np.ndenumerate(area):
                #if np.ma.is_masked(a):continue
                self.weightsDict[(lats[i, j], lons[i, j])] = a

        if lats.ndim == 1:
            for (i, j), a in np.ndenumerate(area):
                #if np.ma.is_masked(a):continue
                self.weightsDict[(lats[i], lons[j])] = a

        if self.debug:
            print("timeseriesAnalysis:\t loadModelWeightsDict.",
                  list(self.weightsDict.keys())[0])

    def loadModelwcvDict(self, ):
        """
  	Adding Water Column Volume (WCV) dictionainy for Model.
  	"""

        nc = dataset(self.gridFile, 'r')
        tmask = nc.variables['tmask'][:]
        try:
            area = nc.variables['area'][:]
        except:
            area = nc.variables['e2t'][:] * nc.variables['e1t'][:]

        #####
        #        mbathy  = (tmask*nc('e3t')).sum(0)

        wcv = (tmask * nc('e3t')).sum(0) * area

        lats = nc.variables['nav_lat'][:]
        lons = nc.variables['nav_lon'][:]
        nc.close()

        self.wcvDict = {}
        if lats.ndim == 2:
            for (i, j), a in np.ndenumerate(wcv):
                if np.ma.is_masked(a): continue
                self.wcvDict[(lats[i, j], lons[i, j])] = a

        if lats.ndim == 1:
            for (i, j), a in np.ndenumerate(wcv):
                if np.ma.is_masked(a): continue
                self.wcvDict[(lats[i], lons[j])] = a

        if self.debug:
            print("timeseriesAnalysis:\t loadModelwcvDict.",
                  list(self.wcvDict.keys())[0])

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
        print("timeseriesAnalysis:\tAddDataArea:\t", area.shape, lats.shape,
              lons.shape)
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

    def loadData(self):

        if self.debug: print("timeseriesAnalysis:\t loadData.", self.dataFile)

        if not self.dataFile:
            if self.debug:
                print("timeseriesAnalysis:\t No data (obs) File provided:",
                      self.dataFile)
            self.dataD = {}
            return

        if not os.path.exists(self.dataFile):
            if self.debug:
                print("timeseriesAnalysis:\tWARNING:\t No such data File:",
                      self.dataFile)
            self.dataD = {}
            return

        ###############
        # load and calculate the real data info
        try:
            if self.clean:
                print(
                    "timeseriesAnalysis:\t loadData\tUser requested clean run. Wiping old data."
                )
                assert 0
            sh = shOpen(self.shelvefn_insitu)
            dataD = sh['dataD']
            sh.close()
            print("timeseriesAnalysis:\t loadData\tOpened shelve:",
                  self.shelvefn_insitu)
            self.dataD = dataD
        except:
            dataD = {}
            print("timeseriesAnalysis:\t loadData\tCould not open in situ shelve:",
                  self.shelvefn_insitu)

        ###############
        # Test to find out if we need to load the netcdf, or if we can just return the dict as a self.object.
        needtoLoad = False
        for r in self.regions:
            if needtoLoad: continue
            for l in self.layers:
                if needtoLoad: continue
                try:
                    print(
                        "timeseriesAnalysis:\t loadData\tChecking if need to Load:",
                        needtoLoad, (r, l), 'len:',
                        len(self.dataD[(r, l)]), self.dataD[(r, l)].shape,
                        np.ma.mean(self.dataD[(r, l)]))
                except:
                    needtoLoad = True
                    print("timeseriesAnalysis:\t loadData\tUnable to Load:",
                          needtoLoad, (r, l))

        if needtoLoad: pass
        else:
            print(
                "timeseriesAnalysis:\t loadData\tDon't need to Load from scratch",
                list(dataD.keys()))
            self.dataD = dataD
            return

        ###############
        # Loading data for each region.
        print("timeseriesAnalysis:\t loadData,\tloading ", self.dataFile)
        #nc = dataset(self.dataFile,'r')
        #data = tst.loadData(nc, self.datadetails)

        ###############
        # Loading data for each region.
        dl = tst.DataLoader(
            self.dataFile,
            '',
            self.datacoords,
            self.datadetails,
            regions=self.regions,
            layers=self.layers,
        )
        if not self.__madeDataArea__: self.AddDataArea()
        for r in self.regions:
            for l in self.layers:
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

                print("timeseriesAnalysis:\t load in situ data,\tloaded ",
                      (r, l), 'mean:', meandatad)
                dataD[(r, l, 'lat')] = dl.load[(r, l, 'lat')]
                dataD[(r, l, 'lon')] = dl.load[(r, l, 'lon')]
                if len(
                        ukp.intersection([
                            'mean',
                            'median',
                            'sum',
                        ], self.metrics)):
                    dataD[(r, l, 'area')] = self.loadDataAreas(
                        dataD[(r, l, 'lat')], dataD[(r, l, 'lon')])
                else:
                    dataD[(r, l, 'area')] = np.ones_like(dataD[(r, l, 'lon')])

                if not meandatad and not datadmask:  #np.ma.is_masked(dataD[(r,l)]):
                    dataD[(r, l)] = np.ma.array([
                        -999,
                    ], mask=[
                        True,
                    ])
                    dataD[(r, l, 'lat')] = np.ma.array([
                        -999,
                    ], mask=[
                        True,
                    ])
                    dataD[(r, l, 'lon')] = np.ma.array([
                        -999,
                    ], mask=[
                        True,
                    ])
                    dataD[(r, l, 'area')] = np.ma.array([
                        -999,
                    ],
                                                        mask=[
                                                            True,
                                                        ])
        #if meandatad and dataD[(r,l)]  == np.ma.array([-999,],mask=[True,]):
        #	print "Massive failiure here:",meandatad, dataD[(r,l)] ,dl.load[(r,l,)]
        #	assert 0
                print("timeseriesAnalysis:\t loadData,\tloading ", (r, l),
                      'mean:', meandatad)

        ###############
        # Savng shelve
        print("timeseriesAnalysis:\t loadData.\tSaving shelve:",
              self.shelvefn_insitu)
        try:
            sh = shOpen(self.shelvefn_insitu)
            sh['dataD'] = dataD
            sh.close()
        except:
            print(
                "timeseriesAnalysis:\t WARNING.\tSaving shelve failed, trying again.:",
                self.shelvefn_insitu)
            shutil.move(self.shelvefn_insitu, self.shelvefn_insitu + '.broken')
            sh = shOpen(self.shelvefn_insitu)
            sh['dataD'] = dataD
            sh.close()

        self.dataD = dataD

    def mapplotsRegionsLayers(self, ):
        """	
        Makes a map plot of model vs data for each string-named layer (not numbered layers). 
  	"""
        newlayers = [
            l for l in self.layers
            if type(l) not in [type(0), type(0.)]
        ]
        fn = self.modelFiles[-1]
        print("mapplotsRegionsLayers:\tLoading file:", fn)
        mDL = tst.DataLoader(
            fn,
            '',
            self.modelcoords,
            self.modeldetails,
            regions=self.regions,
            layers=newlayers,
        )
        nc = dataset(fn, 'r')
        ts = tst.getTimes(nc, self.modelcoords)

        if len(ts) > 1:
            maxtime = ts.max()
            maxtime_index = mDL.timedict_ti[maxtime]
            timestr = str(int(maxtime))
        else:
            timestr = str(int(ts[0]))

    #meantime = np.mean(ts)
        cbarlabel = getLongName(
            self.modeldetails['name']) + ', ' + getLongName(
                self.modeldetails['units'])
        for r in self.regions:
            for l in self.layers:
                if type(l) in [type(0), type(0.)]: continue
                mapfilename = ukp.folder(self.imageDir + '/' +
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
                modelt = mDL.load[(r, l, 't')]

                if len(ts) > 1:
                    timemask = np.ma.masked_where(modelt != maxtime_index,
                                                  modelt).mask
                    modeldata = np.ma.masked_where(timemask,
                                                   modeldata).compressed()
                    modellat = np.ma.masked_where(timemask,
                                                  modellat).compressed()
                    modellon = np.ma.masked_where(timemask,
                                                  modellon).compressed()

                if not len(modeldata): continue
                if modellat.mean() == 0. and modellon.mean() == 0.: continue

                print("mapplotsRegionsLayers:\t", r, l, "model contains",
                      len(modeldata), 'model data')
                print("mapplotsRegionsLayers:\t", r, l, "model lat:",
                      len(modellat), modellat.min(), modellat.mean(),
                      modellat.max())
                print("mapplotsRegionsLayers:\t", r, l, "model lon:",
                      len(modellon), modellon.min(), modellon.mean(),
                      modellon.max())

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

                print("mapplotsRegionsLayers:\t", r, l, "contains",
                      len(datadata), 'in situ data')
                print("mapplotsRegionsLayers:\t", r, l, "data lat:",
                      len(datalat), datalat.min(), datalat.mean(),
                      datalat.max())
                print("mapplotsRegionsLayers:\t", r, l, "data lon:",
                      len(datalon), datalon.min(), datalon.mean(),
                      datalon.max())

                for title_string in [
                        self.model, '(' + self.jobID + ')',
                        str(l), self.modeldetails['name'], timestr,
                        self.datasource,
                        str(l), self.datadetails['name']
                ]:
                    print("mapplotsRegionsLayers:", title_string, end=' ')
                    print(getLongName(title_string))

                titles = [
                    ' '.join([
                        getLongName(t) for t in [
                            self.model, '(' + self.jobID + ')',
                            str(l), self.modeldetails['name'], timestr
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
                    cbarlabel=cbarlabel,
                    dpi=100,
                )

    def makePlots(self):
        if self.debug: print("timeseriesAnalysis:\t makePlots.")

        #####
        # Trafficlight and percentiles plots:
        for r in self.regions:
            for l in self.layers:
                #####
                # Don't make pictures for each integer or float layer, only the ones that are strings.
                if type(l) in [type(0), type(0.)]: continue
                if self.debug:
                    print("timeseriesAnalysis:\t makePlots.\t", r, l)

                #####
                # Test for presence/absence of in situ data.
                try:
                    dataslice = self.dataD[(r, l)]
                    dataweights = self.dataD[(r, l, 'area')]
                    print(
                        "timeseriesAnalysis:\t makePlots, \tLoaded In situ data:",
                        (r, l))
                except:
                    dataslice = []
                    dataweights = []
                    print(
                        "timeseriesAnalysis:\t makePlots, \tNo In situ data:",
                        (r, l))
                try:
                    dataslice = np.ma.masked_where(
                        dataslice.mask + dataweights.mask, dataslice)
                    dataweights = np.ma.masked_where(
                        dataslice.mask + dataweights.mask, dataweights)
                    dataslice = dataslice.compressed()
                    dataweights = dataweights.compressed()
                except:
                    print(
                        "timeseriesAnalysis:\t makePlots, \tCan not compress In situ data:",
                        (r, l))

                if len(dataslice) != len(dataweights):
                    print(
                        "timeseriesAnalysis:\t makePlots, \tlen(dataslice) != len(dataweights)"
                        + str(len(dataslice)) + ' != ' + str(len(dataweights)))
                    print(type(dataslice), type(dataweights))
                    assert 0

                #####
                # Percentiles plots.
                if '20pc' in self.metrics:  #continue
                    modeldataDict = {}
                    timesDict = {}
                    for m in self.metrics:
                        timesDict[m] = sorted(self.modeldataD[(r, l,
                                                               m)].keys())
                        #modeldataDict[m] = [self.modeldataD[(r,l,m)][t] for t in timesDict[m]]
                        modeldataDict[m] = []
                        for t in sorted(timesDict[m]):
                            v = self.modeldataD[(r, l, m)][t]
                            if np.ma.is_masked(v): modeldataDict[m].append(0.)
                            else: modeldataDict[m].append(v)

                    title = ' '.join([
                        getLongName(t)
                        for t in [r, str(l), self.datasource, self.dataType]
                    ])
                    for greyband in [
                            '10-90pc',
                    ]:  #'MinMax',
                        filename = ukp.folder(self.imageDir + '/' +
                                              self.dataType) + '_'.join([
                                                  'percentiles', self.jobID,
                                                  self.dataType, r,
                                                  str(l), greyband
                                              ]) + '.png'
                        if self.debug:
                            print(
                                "timeseriesAnalysis:\t makePlots.\tInvestigating:",
                                filename)

                        if not ukp.shouldIMakeFile(
                            [self.shelvefn, self.shelvefn_insitu],
                                filename,
                                debug=False):
                            continue
                        tsp.percentilesPlot(timesDict,
                                            modeldataDict,
                                            dataslice,
                                            dataweights=dataweights,
                                            title=title,
                                            filename=filename,
                                            units=self.modeldetails['units'],
                                            greyband=greyband)
                else:
                    print(
                        "timeseriesAnalysis:\t makePlots, \tNo percentile plots:",
                        (r, l))

            #####
            # simpletimeseries plots.
                for m in self.metrics:
                    if m not in [
                            'mean',
                            'metricless',
                            'sum',
                            'wcvweighted',
                    ]:
                        continue
                    filename = ukp.folder(self.imageDir + '/' +
                                          self.dataType) + '_'.join([
                                              m,
                                              self.jobID,
                                              self.dataType,
                                              r,
                                              str(l),
                                              m,
                                          ]) + '.png'
                    if self.debug:
                        print(
                            "timeseriesAnalysis:\t makePlots.\tInvestigating simpletimeseries: ",
                            filename)
                    if not ukp.shouldIMakeFile(
                        [self.shelvefn + '*', self.shelvefn_insitu + '*'],
                            filename,
                            debug=False):
                        continue

                    modeldataDict = self.modeldataD[(r, l, m)]
                    times = sorted(modeldataDict.keys())
                    modeldata = [modeldataDict[t] for t in times]
                    title = ' '.join([
                        getLongName(t)
                        for t in [r, str(l), m, self.dataType]
                    ])

                    if dataslice == []:
                        datamean = -999.
                    elif m in [
                            'mean',
                            'metricless',
                    ]:
                        if len(dataweights) != 0 and dataweights.sum() != 0.:
                            datamean = np.average(dataslice,
                                                  weights=dataweights)
                        else:
                            datamean = np.mean(dataslice)
                    elif m in [
                            'sum',
                    ]:
                        if len(dataweights) != 0 and dataweights.sum() != 0.:
                            datamean = np.sum(dataslice, weights=dataweights)
                        else:
                            datamean = np.sum(dataslice)

                    tsp.simpletimeseries(times,
                                         modeldata,
                                         datamean,
                                         title=title,
                                         filename=filename,
                                         units=self.modeldetails['units'],
                                         greyband=False)

        #####
        # map plots for specific regions:
        runmapplots = False
        for r in self.regions:
            for l in self.layers:
                if not runmapplots: continue
                mapfilename = ukp.folder(self.imageDir + '/' +
                                         self.dataType) + '_'.join([
                                             'map',
                                             self.jobID,
                                             self.dataType,
                                             str(l),
                                             r,
                                         ]) + '.png'
                if not ukp.shouldIMakeFile(self.modelFiles[-1],
                                       mapfilename,
                                       debug=False):
                    runmapplots = False
        if runmapplots:
            self.mapplotsRegionsLayers()
