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
.. module:: p2pPlots
   :platform: Unix
   :synopsis: A tool to make plots for point to point analysis.
.. moduleauthor:: Lee de Mora <ledm@pml.ac.uk>

"""

from netCDF4 import Dataset, num2date
from datetime import datetime
from sys import argv
from os.path import exists, split, getmtime, basename
from glob import glob
from shelve import open as shOpen
from matplotlib.colors import LogNorm
from matplotlib import pyplot, ticker
from calendar import month_name
from itertools import product
from scipy.stats import linregress
from scipy.stats.mstats import scoreatpercentile
import numpy as np

#local imports
from ..bgcvaltools.StatsDiagram import StatsDiagram
from ..bgcvaltools.robust import StatsDiagram as robustStatsDiagram
from .. import bgcvaltools.unbiasedSymmetricMetrics as usm
from .. import UKESMpython as ukp
from ..bgcvaltools.pftnames import getLongName, AutoVivification, fancyUnits, CMIP5models  # getmt
from ..bgcvaltools.makeMask import makeMask
from .slicesDict import populateSlicesList, slicesDict

#from bgcvaltools.pftnames import MaredatTypes,IFREMERTypes,WOATypes,GEOTRACESTypes

#import seaborn as sb
"""	This code makes matched plots, hexbins, scatter plots, and so on.

"""

#BioLogScales 	= ['bac','mesozoo','diatoms','picophyto','microzoo','PP','Seawifs', 'iron']

noXYLogs = [
    'pCO2',
    #'nitrateSurface', 	'nitrateAll',	'nitrateTransect',
    #'phosphateSurface',	'phosphateAll',	'phosphateTransect',
    'silicateSurface',
    'silicateAll',
    'silicateTransect',
    'silicate100m',
    'silicate200m',
    'silicate500m',
    'tempSurface',
    'tempAll',
    'tempTransect',
    'temp100m',
    'temp200m',
    'temp1000m',
    'temp500m',
    'salSurface',
    'salAll',
    'salTransect',
    'sal100m',
    'sal200m',
    'sal1000m',
    'sal500m',
]

transectSlices = [
    'All',
    'Global',
]


class makePlots:

    def __init__(self,
                 matchedDataFile,
                 matchedModelFile,
                 name,
                 datasource='',
                 model='',
                 jobID='',
                 year='',
                 depthLevel='',
                 modelcoords='',
                 modeldetails='',
                 datacoords='',
                 datadetails='',
                 shelveDir='',
                 imageDir='',
                 newSlices=['All', 'Standard'],
                 compareCoords=True,
                 noPlots=False,
                 dpi=100):  #xfilename,yfilename,saveShelve=True,

        self.xfn = matchedModelFile
        self.yfn = matchedDataFile
        self.name = name
        self.newSlices = newSlices
        self.depthLevel = depthLevel

        self.xtype = model
        self.ytype = datasource

        self.model = model
        self.jobID = jobID
        self.year = year
        self.shelveDir = shelveDir
        self.compareCoords = compareCoords
        self.months = {month_name[i + 1]: i for i in range(0, 12)}
        self.noPlots = noPlots

        # details about coordinates and what to load.
        self.modelcoords = modelcoords
        self.modeldetails = modeldetails
        self.datacoords = datacoords
        self.datadetails = datadetails
        self.dpi = dpi
        #self.mt = getmt()
        #Models = [m.upper() for m in ['Diat-HadOCC', 'ERSEM','HadOCC', 'MEDUSA','PlankTOM6','PlankTOM10','NEMO','IMARNET','CMIP5',]] # skip these to find in situ data types.
        #Models.extend(['IMARNET_' +m.upper() for m in ['Diat-HadOCC', 'ERSEM','HadOCC', 'MEDUSA','PlankTOM6','PlankTOM10','NEMO',]])
        #Models.extend(['CMIP5_' +m.upper() for m in CMIP5models])
        #ytypes = []
        #for dk in self.mt.keys():
        #	if dk.upper() in Models:
        #		#print "
        #		continue
        #	if self.name in self.mt[dk].keys():ytypes.append(dk)
        #if len(ytypes)==1:
        #	self.ytype = ytypes[0]
        #else:
        #    if len(ytypes)>1:	print "ERROR:\t The same name,(",self.name,") appears in multiple datasets:",ytypes
        #    if len(ytypes)<1:	print "ERROR:\t The name,(",self.name,") not appears in any datasets"
        #  	    print "THis job will probably fa

        #if self.name in MaredatTypes:  	self.ytype = 'Maredat'
        #if self.name in WOATypes:  	self.ytype = 'WOA'
        #if self.name in IFREMERTypes:  	self.ytype = 'IFREMER'
        #if self.name in GEOTRACESTypes: self.ytype = 'GEOTRACES'

        if self.shelveDir == '':
            self.shelveDir = ukp.folder([
                'shelves', self.xtype, self.year, self.ytype, 'Slices',
                self.name + self.depthLevel
            ])
        else:
            self.shelveDir = ukp.folder(self.shelveDir)

        if imageDir == '':

            self.imageDir = ukp.folder([
                'images', self.xtype, 'P2P_plots', self.year,
                self.name + self.depthLevel
            ])
            print("Using default image folder:", self.imageDir)
        else:
            self.imageDir = ukp.folder(imageDir)

        self.run()

    def run(self, ):

        self.xnc = Dataset(self.xfn, 'r')
        self.ync = Dataset(self.yfn, 'r')

        if self.compareCoords: self.CompareCoords()
        #self.defineSlices(self.plotallcuts)

        self.plotWithSlices()

        self.xnc.close()
        self.ync.close()

    def plotWithSlices(self):  #,newSlice,):
        print("plotWithSlices:\txtype:", self.xtype, "\tytype:", self.ytype,
              "\tname:", self.name, self.depthLevel)

        #####
        # Test if any of the plots exist.

        xkeys = []
        ykeys = []

        #nx = self.mt[self.xtype][self.name]

        #if type(nx) == type(['a',]):	xkeys = self.mt[self.xtype][self.name]
        #else:				xkeys.append(self.mt[self.xtype][self.name]['name'])
        #ny = self.mt[self.ytype][self.name]
        #if type(ny) == type(['a',]):	ykeys = self.mt[self.ytype][self.name]
        #else:				ykeys.append(self.mt[self.ytype][self.name]['name'])
        xkeys = [
            self.modeldetails['name'],
        ]
        ykeys = [
            self.datadetails['name'],
        ]

        print("plotWithSlices:\txkeys:", xkeys, '\tykeys:', ykeys)
        if [{}] in [xkeys, ykeys]:
            print(
                "plotWithSlices:\tERROR\t This data type  is not defined in pftnames.py getmt()"
            )
            print(
                "plotWithSlices:\tx:\t['" + self.xtype + "'\t]['" + self.name +
                "'] = ", xkeys)
            print(
                "plotWithSlices:\ty:\t['" + self.ytype + "'\t]['" + self.name +
                "'] = ", ykeys)
            assert False

        #####
        # This section of code is a bit of a time saver.
        # It checks to see if the image and the output shelve exist.
        # If they both exist and and are older than the input netcdfs, the rest of this function is skipped.
        # If one is missing, or the input files are newer than the old image, the function runs as normal.
        # Caveat: if some image can not be made, ie the data makes the mask cover 100% of the data, then the code will run as normal (no skipping).
        self.shelvesAV = []  #AutoVivification()

        plotsToMake = 0
        for newSlice in self.newSlices:
            for xk, yk in product(xkeys, ykeys):
                print('plotWithSlices:\t', newSlice, '\tlisting plotpairs:\tX',
                      xk, ': [', self.xtype, '][', self.name, ']')
                print('plotWithSlices:\t', newSlice, '\tlisting plotpairs:\tY',
                      yk, ': [', self.ytype, '][', self.name, ']')
                print(xk, yk, self.xtype, self.ytype, self.name)

                if type(newSlice) in [type([
                        'a',
                        'b',
                ]), type((
                        'a',
                        'b',
                ))]:
                    ns = ''.join(newSlice)
                else:
                    ns = newSlice

                try:
                    fn = ns + '_' + xk + 'vs' + yk
                except:
                    print(
                        "ERROR:\tcan\'t add ", newSlice, ns, xk, yk,
                        'together as strings. the problem is probably in your mt dictionary in pftnames.'
                    )
                    assert False

        #####
        # Don't make Plots for transects with two spacial cuts.
                if self.depthLevel.lower().find(
                        'transect') > -1 and newSlice not in transectSlices:
                    continue

                #####
                # Does the image exist?
                filename = self.getFileName(newSlice, xk, yk)
                if ukp.shouldIMakeFile([self.xfn, self.yfn],
                                       filename,
                                       debug=False):
                    plotsToMake += 1

        #####
        #Does the shelve file exist?
                shelveName = self.shelveDir + self.name + '_' + ns + '_' + xk + 'vs' + yk + '.shelve'
                if ukp.shouldIMakeFile([self.xfn, self.yfn],
                                       shelveName,
                                       debug=False):
                    plotsToMake += 1

        #####
        # Make a list of shelve meta data, to aid post processing.
                she = ukp.shelveMetadata(model=self.model,
                                         name=self.name,
                                         year=self.year,
                                         depthLevel=self.depthLevel,
                                         newSlice=newSlice,
                                         xkey=xk,
                                         ykey=yk,
                                         shelve=shelveName)
                self.shelvesAV.append(she)  #shelveName
                #self.shelvesAV[newSlice][xk][yk] = shelveName
                try:
                    self.shelves.append(shelveName)
                except:
                    self.shelves = [
                        shelveName,
                    ]

        if plotsToMake == 0:
            print('plotWithSlices:\tAll plots and shelve files already made',
                  self.name, newSlice, xkeys, ykeys)
            return

        #####
        # Load Coordinates
        #time and depth
        self.xt = np.ma.array(self.xnc.variables[self.modelcoords['t']][:])
        self.yt = np.ma.array(self.ync.variables[self.datacoords['t']][:])
        self.xz = np.ma.array(self.xnc.variables[self.modelcoords['z']][:])
        self.yz = np.ma.array(self.ync.variables[self.datacoords['z']][:])

        #lat and lon
        self.xy = np.ma.array(self.xnc.variables[self.modelcoords['lat']][:])
        self.yy = np.ma.array(self.ync.variables[self.datacoords['lat']][:])
        self.xx = np.ma.array(self.xnc.variables[self.modelcoords['lon']][:])
        self.yx = np.ma.array(self.ync.variables[self.datacoords['lon']][:])

        for newSlice in self.newSlices:

            #####
            # Don't make Plots for transects with two spacial cuts.
            if self.depthLevel.lower().find(
                    'transect') > -1 and newSlice not in transectSlices:
                continue

            for xkey, ykey in product(xkeys, ykeys):
                print("plotWithSlices:\t", newSlice, xkey, ykey)
                self.plotsFromKeys(newSlice, xkey, ykey)

    def plotsFromKeys(self, newSlice, xkey, ykey):

        #####
        # check that the plot and shelve should be made
        if type(newSlice) in [type([
                'a',
                'b',
        ]), type((
                'a',
                'b',
        ))]:
            ns = ''.join(newSlice)
        else:
            ns = newSlice
        self.shelveName = self.shelveDir + self.name + '_' + ns + '_' + xkey + 'vs' + ykey + '.shelve'
        filename = self.getFileName(newSlice, xkey, ykey)

        print("plotWithSlices:\tINFO:\tinvestigating:", (newSlice), filename)
        if not ukp.shouldIMakeFile([self.xfn,self.yfn],self.shelveName,debug=False) \
         and not ukp.shouldIMakeFile([self.xfn,self.yfn],filename,debug=False):
            return

        #####
        # Extract remaining data (already know lat,lon,time,depth)
        xd = ukp.extractData(self.xnc, self.modeldetails, key=xkey)
        yd = ukp.extractData(self.ync, self.datadetails, key=ykey)

        #####
        # Build mask
        fullmask = xd.mask.astype(int) + yd.mask.astype(
            int) + np.ma.masked_invalid(xd).mask.astype(
                int) + np.ma.masked_invalid(yd).mask.astype(int)

        if type(newSlice) in [type([
                'a',
        ]), type(('a', ))]:  # newSlice is actaully a list of multiple slices.
            for n in newSlice:
                fullmask += makeMask(self.name, n, self.xt, self.xz, self.xy,
                                     self.xx, xd).astype(int)
                fullmask += makeMask(self.name, n, self.yt, self.yz, self.yy,
                                     self.yx, yd).astype(int)

        elif newSlice == 'Standard':  # Standard is a shorthand for my favourite cuts.
            for stanSlice in slicesDict['StandardCuts']:
                if self.name in ['tempSurface', 'tempTransect', 'tempAll'
                                 ] and stanSlice in [
                                     'aboveZero',
                                 ]:
                    continue

                fullmask += makeMask(self.name, stanSlice, self.xt, self.xz,
                                     self.xy, self.xx, xd).astype(int)
                fullmask += makeMask(self.name, stanSlice, self.yt, self.yz,
                                     self.yy, self.yx, yd).astype(int)

        else:  # newSlice is a simple slice.
            fullmask += makeMask(self.name, newSlice, self.xt, self.xz,
                                 self.xy, self.xx, xd).astype(int)
            fullmask += makeMask(self.name, newSlice, self.yt, self.yz,
                                 self.yy, self.yx, yd).astype(int)
            print('plotWithSlices:\t', fullmask.sum())

        if self.name in ['mld', 'mld_DT02', 'mld_DR003', 'mld_DReqDTm02']:
            mldMask = self.ync.variables['mask'][:]
            fullmask += np.ma.masked_where(mldMask == 0., mldMask).mask

        N = len(self.xt)

        maskcoverpc = 100. * np.clip(fullmask, 0, 1).sum() / float(N)
        if maskcoverpc == 100.:
            print("plotWithSlices:\tNew Mask,", newSlice,
                  ", covers entire dataset.", maskcoverpc, '%', N)
            try:
                self.shelves[newSlice][xk][yk] = ''
            except:
                pass
            return
        print("plotWithSlices:\tNew Mask,", newSlice, ", covers ", maskcoverpc,
              '% of ', N, 'data')

        #####
        # Apply mask to all data.
        nmxx = np.ma.masked_where(fullmask, self.xx).compressed()
        nmxy = np.ma.masked_where(fullmask, self.xy).compressed()
        nmxz = np.ma.masked_where(fullmask, self.xz).compressed()
        nmxt = np.ma.masked_where(fullmask, self.xt).compressed()
        nmyx = np.ma.masked_where(fullmask, self.yx).compressed()
        nmyy = np.ma.masked_where(fullmask, self.yy).compressed()
        nmyz = np.ma.masked_where(fullmask, self.yz).compressed()
        nmyt = np.ma.masked_where(fullmask, self.yt).compressed()
        datax = np.ma.masked_where(fullmask, xd).compressed()
        datay = np.ma.masked_where(fullmask, yd).compressed()

        print("plotWithSlices:\tlenghts", [len(datax), len(datay)], 'x:\t',
              [len(nmxx), len(nmxy)], 'y:\t', [len(nmxz), len(nmyx)], 'z:\t',
              [len(nmyy), len(nmyz)])
        if 0 in [
                len(datax),
                len(datay),
                len(nmxx),
                len(nmxy),
                len(nmxz),
                len(nmyx),
                len(nmyy),
                len(nmyz)
        ]:
            print('plotWithSlices:\tWARNING:\tslice:', newSlice,
                  'There is a zero in one of the fields.')
            #try:	self.shelvesAV[newSlice][xk][yk] = ''
            #except:	pass
            return

        dmin = min([datax.min(), datay.min()])
        dmax = max([datax.max(), datay.max()])
        if dmin == dmax:
            print("plotWithSlices:\tWARNING:\tminimum == maximum,\t (", dmin,
                  ' == ', dmax, ')')
            #try:	self.shelvesAV[newSlice][xk][yk] = ''
            #except:	pass
            return

        #####
        # Prepare units, axis labels and titles.
        if 'units' in list(self.modeldetails.keys()):
            xunits = fancyUnits(self.modeldetails['units'])
        else:
            try:
                xunits = fancyUnits(
                    self.xnc.variables[self.modeldetails['vars'][0]].units,
                    debug=True)
            except:
                print(
                    "plotWithSlices:\tWARNING:\tno units provided for model ",
                    self.modeldetails['name'], 'in details dictionairy')
                xunits = ''

        if 'units' in list(self.datadetails.keys()):
            yunits = fancyUnits(self.datadetails['units'])
        else:
            try:
                yunits = fancyUnits(
                    self.ync.variables[self.datadetails['vars'][0]].units,
                    debug=True)
            except:
                print("plotWithSlices:\tWARNING:\tno units provided for data ",
                      self.datadetails['name'], 'in details dictionairy')
                yunits = ''

        labelx = getLongName(self.xtype) + ' ' + getLongName(
            self.name) + ', ' + xunits
        labely = getLongName(self.ytype) + ' ' + getLongName(
            self.name) + ', ' + yunits

        try:
            title = getLongName(newSlice) + ' ' + getLongName(
                self.name + self.depthLevel)  #+getLongName(self.name)
        except:
            title = newSlice + ' ' + xkey + ' vs ' + ykey

        scatterfn = filename.replace('.png', '_scatter.png')
        robfnxy = filename.replace('.png', '_xyrobin.png')
        robfnquad = filename.replace('.png', '_robinquad.png')
        robfncartopy = filename.replace('.png', '_robinquad-cartopy.png')
        transectquadfn = filename.replace('.png', '_transect.png')
        histfnxy = filename.replace('.png', '_hist.png')
        histsfnxy = filename.replace('.png', '_hists.png')

        #####
        # Can turn off plots to run analysis faster.
        if self.noPlots:
            print("plotWithSlices:\tSkipping plots ...")
        else:
            #####
            # Robinson projection plots - Basemap
            mptbasemap = True  # Don't need both.
            if mptbasemap:
                if ukp.shouldIMakeFile([self.xfn, self.yfn],
                                       robfnquad,
                                       debug=False):
                    ti1 = getLongName(self.xtype)
                    ti2 = getLongName(self.ytype)
                    cbarlabel = xunits
                    if self.name in noXYLogs or dmin * dmax <= 0.:
                        doLog = False
                    else:
                        doLog = True
                    print("plotWithSlices:\tROBIN QUAD:", [ti1, ti2], False,
                          dmin, dmax)
                    ukp.robinPlotQuad(
                        nmxx,
                        nmxy,
                        datax,
                        datay,
                        robfnquad,
                        titles=[ti1, ti2],
                        title=' '.join([
                            getLongName(newSlice),
                            getLongName(self.name),
                            getLongName(self.depthLevel), self.year
                        ]),
                        cbarlabel=cbarlabel,
                        doLog=doLog,
                        vmin=dmin,
                        vmax=dmax,
                    )

            # Robinson projection plots - Cartopy
            #makeCartopy = True	# Don't need both.
            if newSlice == 'Global' and self.depthLevel in [
                    'Surface',
                    '100m',
                    '200m',
                    '500m',
                    '1000m',
            ]:
                # ####
                # Global, as we have interpollation turned on here.
                if ukp.shouldIMakeFile([self.xfn, self.yfn],
                                       robfncartopy,
                                       debug=False):
                    ti1 = getLongName(self.xtype)
                    ti2 = getLongName(self.ytype)
                    cbarlabel = xunits
                    if self.name in noXYLogs or dmin * dmax <= 0.:
                        doLog = False
                    else:
                        doLog = True
                    print("plotWithSlices:\tROBIN QUAD:", [ti1, ti2], False,
                          dmin, dmax)
                    try:
                        ukp.robinPlotQuad(nmxx,
                                          nmxy,
                                          datax,
                                          datay,
                                          robfncartopy,
                                          titles=[ti1, ti2],
                                          title=' '.join([
                                              getLongName(self.name),
                                              getLongName(self.depthLevel),
                                              self.year
                                          ]),
                                          cbarlabel=cbarlabel,
                                          doLog=doLog,
                                          vmin=dmin,
                                          vmax=dmax,
                                          maptype='Cartopy',
                                          scatter=False)
                    except:
                        print("Cartopy is broken again, can't make: ",
                              robfncartopy)

            if self.depthLevel not in [
                    'Surface',
                    '100m',
                    '200m',
                    '500m',
                    '1000m',
            ]:  # No point in making these.
                if ukp.shouldIMakeFile([self.xfn, self.yfn],
                                       transectquadfn,
                                       debug=False):
                    ti1 = getLongName(self.xtype)
                    ti2 = getLongName(self.ytype)
                    cbarlabel = xunits
                    if self.name in noXYLogs or dmin * dmax <= 0.:
                        doLog = False
                    else:
                        doLog = True
                    print("plotWithSlices:\ttransect quad:", [ti1, ti2], False,
                          dmin, dmax)
                    if self.depthLevel in [
                            'ArcTransect',
                            'AntTransect',
                            'CanRusTransect',
                    ]:
                        ukp.ArcticTransectPlotQuad(
                            nmxx,
                            nmxy,
                            nmxz,
                            datax,
                            datay,
                            transectquadfn,
                            titles=[ti1, ti2],
                            title=' '.join([
                                getLongName(self.name),
                                getLongName(self.depthLevel), self.year
                            ]),
                            cbarlabel=cbarlabel,
                            doLog=doLog,
                            vmin=dmin,
                            vmax=dmax,
                            scatter=False,
                            logy=True,
                            transectName=self.depthLevel,
                        )
                    else:
                        ukp.HovPlotQuad(
                            nmxx,
                            nmxy,
                            nmxz,
                            datax,
                            datay,
                            transectquadfn,
                            titles=[ti1, ti2],
                            title=' '.join([
                                getLongName(self.name),
                                getLongName(self.depthLevel), self.year
                            ]),
                            cbarlabel=cbarlabel,
                            doLog=doLog,
                            vmin=dmin,
                            vmax=dmax,
                            scatter=False,
                            logy=True,
                        )

            #####
            # Simultaneous histograms plot	- single
            if ukp.shouldIMakeFile([self.xfn, self.yfn], histfnxy,
                                   debug=False):
                xaxislabel = getLongName(self.name) + ', ' + xunits
                labelx = self.xtype
                labely = self.ytype
                histtitle = title
                histxaxis = xaxislabel
                if self.ytype in ['LANA', 'LANA_p']:
                    labelx = getLongName(self.name)
                    labely = getLongName(self.ytype)
                    histtitle = getLongName(
                        newSlice) + ' DMS: ' + labelx + ' vs ' + labely
                    histxaxis = 'DMS, ' + xunits

                if self.name in noXYLogs or dmin * dmax <= 0.:
                    ukp.histPlot(datax,
                                 datay,
                                 histfnxy,
                                 Title=histtitle,
                                 labelx=labelx,
                                 labely=labely,
                                 dpi=self.dpi,
                                 xaxislabel=histxaxis)
                else:
                    ukp.histPlot(
                        datax,
                        datay,
                        histfnxy,
                        Title=histtitle,
                        labelx=labelx,
                        labely=labely,
                        dpi=self.dpi,
                        xaxislabel=histxaxis,
                        logx=True,
                    )

            # Simultaneous histograms plot	- triple
            if ukp.shouldIMakeFile([self.xfn, self.yfn],
                                   histsfnxy,
                                   debug=False):
                xaxislabel = getLongName(self.name) + ', ' + xunits
                if self.name in noXYLogs or dmin * dmax <= 0.:
                    ukp.histsPlot(datax,
                                  datay,
                                  histsfnxy,
                                  Title=title,
                                  labelx=self.xtype,
                                  labely=self.ytype,
                                  xaxislabel=xaxislabel)
                else:
                    ukp.histsPlot(
                        datax,
                        datay,
                        histsfnxy,
                        Title=title,
                        labelx=self.xtype,
                        labely=self.ytype,
                        xaxislabel=xaxislabel,
                        logx=True,
                    )

            #####
            # Scatter  (hexbin) plot
            if ukp.shouldIMakeFile([self.xfn, self.yfn],
                                   scatterfn,
                                   debug=False):
                gs = 50
                scattitle = title
                slabelx = labelx
                slabely = labely
                if self.ytype in ['LANA', 'LANA_p']:
                    slabelx = getLongName(self.name) + ' DMS, ' + xunits
                    slabely = getLongName(self.ytype) + ' DMS, ' + xunits
                    scattitle = getLongName(newSlice) + ' DMS: ' + getLongName(
                        self.name) + ' vs ' + getLongName(self.ytype)

                    pass

                if self.name in noXYLogs or dmin * dmax <= 0.:
                    ukp.scatterPlot(datax,
                                    datay,
                                    scatterfn,
                                    Title=scattitle,
                                    labelx=slabelx,
                                    labely=slabely,
                                    dpi=self.dpi,
                                    bestfitLine=True,
                                    gridsize=gs)
                else:
                    ukp.scatterPlot(
                        datax,
                        datay,
                        scatterfn,
                        Title=scattitle,
                        labelx=slabelx,
                        labely=slabely,
                        dpi=self.dpi,
                        bestfitLine=True,
                        gridsize=gs,
                        logx=True,
                        logy=True,
                    )

        #####
        # Save fit in a shelve file.
        s = shOpen(self.shelveName)
        print("plotWithSlices:\tSaving ", self.shelveName)
        b1, b0, rValue, pValue, stdErr = linregress(datax, datay)
        print("plotWithSlices:\tlinear regression: \n\tb1:", b1, "\n\tb0:", b0,
              "\n\trValue:", rValue, "\n\tpValue:", pValue, "\n\tstdErr:",
              stdErr)
        s['b1'] = b1
        s['b0'] = b0
        s['rValue'] = rValue
        s['pValue'] = pValue
        s['stdErr'] = stdErr
        s['N'] = len(datax)

        mtaylor = StatsDiagram(datax, datay)
        s['Taylor.E0'] = mtaylor.E0
        s['Taylor.E'] = mtaylor.E
        s['Taylor.R'] = mtaylor.R
        s['Taylor.p'] = mtaylor.p
        s['Taylor.gamma'] = mtaylor.gamma

        s['MNAFE'] = usm.MNAFE(datax, datay)
        s['MNFB'] = usm.MNFB(datax, datay)
        s['NMAEF'] = usm.NMAEF(datax, datay)
        s['NMBF'] = usm.NMBF(datax, datay)

        mrobust = robustStatsDiagram(datax, datay, 0.01)
        print(
            "makePlots.py:\tWARNING: robustStatsDiagram CALCULATED WITH DEFAULT PRECISION (0.01)"
        )
        s['robust.E0'] = mrobust.E0
        s['robust.E'] = mrobust.E
        s['robust.R'] = mrobust.R
        s['robust.p'] = mrobust.p
        s['robust.gamma'] = mrobust.gamma

        s['datax'] = datax
        s['datay'] = datay

        s['x_lon'] = nmxx
        s['x_lat'] = nmxy
        s['x_depth'] = nmxz
        s['x_time'] = nmxt
        s['y_lon'] = nmyx
        s['y_lat'] = nmyy
        s['y_depth'] = nmyz
        s['y_time'] = nmyt

        s['title'] = title
        s['labelx'] = labelx
        s['labely'] = labely
        s['name'] = self.name
        s['depthLevel'] = self.depthLevel
        s['region'] = self.depthLevel
        s['year'] = self.year
        s['xtype'] = self.xtype
        s['ytype'] = self.ytype
        s['xfn'] = self.xfn
        s['yfn'] = self.yfn
        s['slice'] = newSlice
        s['newSlice'] = ns
        s['xkey'] = xkey
        s['ykey'] = ykey
        s.close()

    def CompareCoords(self, ):
        """	This routine plots the coordinates of the data against the coordinates of the model.
		This should produce a straight line plot, ensuring that the matching has been performed correctly.
	"""
        #import seaborn as sb
        #sb.set(style="ticks")
        xcoords = [
            self.modelcoords[k] for k in [
                't',
                'lat',
                'lon',
                'z',
                'lon',
            ]
        ]
        ycoords = [
            self.datacoords[k] for k in [
                't',
                'lat',
                'lon',
                'z',
                'lat',
            ]
        ]

        for xkey, ykey in zip(xcoords, ycoords):
            if xkey not in list(self.xnc.variables.keys()): continue
            if ykey not in list(self.ync.variables.keys()): continue
            filename = self.imageDir + 'CompareCoords' + self.name + xkey + 'vs' + ykey + '.png'
            if not ukp.shouldIMakeFile(
                [self.xfn, self.yfn], filename, debug=False):
                continue
            print("CompareCoords:\tx:", xkey, "\ty:", ykey)
            if xkey not in list(self.xnc.variables.keys()):
                print(xkey, "not in xnc")
                assert False
            if ykey not in list(self.ync.variables.keys()):
                print(ykey, "not in ync")
                assert False

            mask = np.ma.array(self.xnc.variables[xkey][:]).mask + np.ma.array(
                self.ync.variables[ykey][:]).mask
            dx = np.ma.masked_where(
                mask, np.ma.array(self.xnc.variables[xkey][:])).compressed()
            dy = np.ma.masked_where(
                mask, np.ma.array(self.ync.variables[ykey][:])).compressed()

            print("CompareCoords:\t", xkey, ':', len(dx), "\t:", ykey, ':',
                  len(dy), dx.min(), dx.max(), dy.min(), dy.max())

            fig = pyplot.figure()
            fig.set_size_inches(8, 12)
            ax = pyplot.subplot(411)

            rects1 = pyplot.hist((dx, dy),
                                 label=[xkey, ykey],
                                 histtype='bar',
                                 bins=72 / 2)  #,alpha=0.2)
            pyplot.legend()
            ax.set_yscale('log')

            ax.set_title(xkey + ' and ' + ykey)
            ax = pyplot.subplot(412)
            rects3 = pyplot.hist(dx - dy, bins=72, label=[xkey + ' - ' + ykey])
            pyplot.legend()
            ax.set_yscale('log')

            ax = pyplot.subplot(212)
            pyplot.hexbin(dx,
                          dy,
                          bins='log',
                          gridsize=72,
                          cmap=pyplot.get_cmap('gist_yarg'),
                          mincnt=0)  #extent=plotrange,
            cb = pyplot.colorbar()
            mmax = max(dx.max(), dy.max())
            mmin = min(dx.min(), dy.min())

            try:
                fx = np.arange(mmin, mmax, (mmax - mmin) / 20.)
            except:
                print("CompareCoords: unable to make plot for \t", xkey,
                      'np.arange(', mmin, mmax, (mmax - mmin) / 20.,
                      ') fails.')
                continue
            pyplot.plot(fx, fx, 'k--')
            ax.set_aspect("equal")

            pyplot.xlabel(self.xtype + ' ' + xkey)
            pyplot.ylabel(self.ytype + ' ' + ykey)

            print("\tSaving: " + filename)
            pyplot.savefig(
                filename,
                dpi=self.dpi,
            )
            pyplot.close()

    def getFileName(self, newSlice, xkey, ykey):
        #####
        # This needs some work.
        file_prefix = self.imageDir  #ukp.folder(['images',self.xtype,'P2P_plots',self.year,self.name+self.region,])

        file_suffix = '_' + self.xtype + '_' + self.year + '.png'

        for dictkey, dictlist in list(slicesDict.items()):
            if dictkey == 'AllSlices': continue
            if newSlice not in dictlist: continue
            if type(newSlice) in [type([
                    'a',
                    'b',
            ]), type((
                    'a',
                    'b',
            ))]:
                newSlice = list(newSlice)
                for i, n in enumerate(newSlice):
                    if n in slicesDict['Months']:
                        newSlice[i] = ukp.mnStr(self.months[n] + 1) + n
                newSlice = ''.join(newSlice)
            if newSlice in slicesDict['Months']:
                newSlice = ukp.mnStr(self.months[newSlice] + 1) + newSlice
            if dictkey == 'Default': dictkey = ''
        filename = ukp.folder(
            [file_prefix, dictkey]
        ) + self.name + self.depthLevel + '_' + newSlice + '_' + xkey + 'vs' + ykey + file_suffix

        return filename


if __name__ == "__main__":
    print("makePlots isn't written to be run as a __main__")
    print("Look at testsuite_p2p.py for examples on how to run this.")
    print('The end.')
