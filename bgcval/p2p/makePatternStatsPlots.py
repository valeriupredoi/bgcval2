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
.. module:: makePatternStatsPlots
   :platform: Unix
   :synopsis: A tool to make pattern statistics plots.
.. moduleauthor:: Lee de Mora <ledm@pml.ac.uk>

"""

from matplotlib import pyplot
from glob import glob

from bgcvaltools.StatsDiagram import rmsds
import bgcvaltools.unbiasedSymmetricMetrics as usm
from netCDF4 import Dataset

import UKESMpython as ukp
from itertools import product
import os
import numpy as np

try:
    from pyproj import Proj
except:
    print("Unable to import proj")
from shapely.geometry import shape
from os.path import basename, exists
from sys import argv
from shelve import open as shOpen
from bgcvaltools.pftnames import getLongName

purple = [125. / 256., 38. / 256., 205. / 256.]


class makePatternStatsPlots:

    def __init__(self,
                 shelveDict,
                 plotTitle,
                 xkeys,
                 filenamebase,
                 grid='ORCA1',
                 gridFile=''):
        """	makePatternStatsPlots:
  		
		Similarly to makeTargets.py, this routine produces a set of plots based on the shelves made by makePlots.py.
  		
  		It takes:
  		  shelveDict: a dictionairy showing the line label and the shelves.
  		  plotTitle: a name for the set of points on the x-axis. ie months, oceans
  		  xkeys: are the tick lables on the x axis. Can be strings or tuple. These need to appear in the shelve filename.
  		  filenamebase: the base name of the image to be saved.
		  grid: which grid to be used (ORCA1, ORCA025, Flat1deg.) 
		  	the grid needs to be in UKESMpython.
  		The labels of the different lines are taken from he keys of the shelve dictionairy.
  		
  	"""
        #if key_self.keys =='dmspmetrics': 	self.keys = ['dms_p_and','dms_p_ara','dms_p_hal','dms_p_sim','In situ',]
        self.keys = sorted(shelveDict.keys())
        print("makePatternStatsPlots:")
        print("		Title:", plotTitle)
        print("		x axis keys:", xkeys)
        print("		grid:", grid)
        print("		legend keys:", self.keys)
        self.shelveDict = shelveDict

        self.plotTitle = plotTitle
        self.filenamebase = filenamebase
        self.xkeys = xkeys
        self.grid = grid
        self.gridFile = gridFile
        self.setDictionaries()

        if not len(self.AllShelves):
            print("makePatternStatsPlots: Found no shelves.", plotTitle, xkeys,
                  shelveDict)
            return

        self.calculateVolume()
        self.loadMetrics()
        if not self.i:
            print("makePatternStatsPlots: Found no shelves.", plotTitle, xkeys,
                  shelveDict)
            return
        self.plotMetrics()

    def calculateVolume(self):
        # makes a netcdf with the volumes of each pixel in the flat map.
        #
        if self.gridFile == '': self.gridFile = ukp.getGridFile(self.grid)
        volumeDict = {}
        if self.grid == 'Flat1deg':
            nc = Dataset(self.gridFile, 'r')
            lat = nc.variables['lat'][:]
            lon = nc.variables['lon'][:]
            print('calculating volume: ')
            for i, la in enumerate(lat):  #(u'lat', u'lon')
                co = {
                    "type":
                    "Polygon",
                    "coordinates": [[
                        (1., la + 0.5),  #('lon', 'lat')
                        (1., la - 0.5),
                        (0., la - 0.5),
                        (0., la + 0.5)
                    ]]
                }
                clon, clat = list(zip(*co['coordinates'][0]))

                try:
                    pa = Proj("+proj=aea +lat_1=" + str(la - 0.5) +
                              " +lat_2=" + str(la + 0.5) + "+lat_0=" +
                              str(la) + " +lon_0=0.5")
                    x, y = pa(clon, clat)
                    cop = {"type": "Polygon", "coordinates": [list(zip(x, y))]}
                    volume = shape(
                        cop).area * 1.  # Assume thickness of 1m for top layer.
                    for lo in lon:
                        volumeDict[(round(la, 3), round(lo, 3))] = volume
                except:
                    print("Unable to use proj", end=' ')
                    for lo in lon:
                        volumeDict[(round(la, 3), round(lo, 3))] = 0.
            nc.close()
            self.volumeDict = volumeDict
            return
        if self.grid in [
                'ORCA1',
                'ORCA025',
        ]:
            nc = Dataset(self.gridFile, 'r')
            pvol = nc.variables['pvol'][0, :, :]
            lat = nc.variables['nav_lat'][:]
            lon = nc.variables['nav_lon'][:]
            m = np.ma.array(pvol).mask + np.ma.array(lat).mask + np.ma.array(
                lon).mask
            pvol = np.ma.masked_where(m, pvol).compressed()
            lat = np.ma.masked_where(m, lat).compressed()
            lon = np.ma.masked_where(m, lon).compressed()

            for la, lo, v in zip(lat, lon, pvol):
                volumeDict[(round(la, 3), round(lo, 3))] = v
            nc.close()
            self.volumeDict = volumeDict
            return

        if self.grid in [
                'eORCA1',
        ]:
            #####
            # No pvol in this mesh file
            nc = Dataset(self.gridFile, 'r')
            pvol = nc.variables['e1t'][:] * nc.variables['e2t'][:]
            lat = nc.variables['nav_lat'][:]
            lon = nc.variables['nav_lon'][:]
            m = np.ma.array(pvol).mask + np.ma.array(lat).mask + np.ma.array(
                lon).mask
            pvol = np.ma.masked_where(m, pvol).compressed()
            lat = np.ma.masked_where(m, lat).compressed()
            lon = np.ma.masked_where(m, lon).compressed()

            for la, lo, v in zip(lat, lon, pvol):
                volumeDict[(round(la, 3), round(lo, 3))] = v
            nc.close()
            self.volumeDict = volumeDict
            return

        print("calculateVolume: Error. Can't do the", self.grid,
              "Grid.   This is not the case for the ORCA1 and ORCA025 grids")
        assert False

    def calcTotalConcentration(self, sh, key='datax'):
        ####
        # Calculates the total quatity in the dataset.
        totalConc = 0.
        #print 'calcTotalConcentration:',sh['xtype'],key,self.grid,self.gridFile
        #print self.volumeDict
        for la, lo, data in zip(self.loadFromshelve(sh, 'x_lat'),
                                self.loadFromshelve(sh, 'x_lon'), sh[key]):
            #print la,lo,data, [(round(la,3),round(lo,3))]
            #print self.volumeDict[(round(la,3),round(lo,3))]
            totalConc += self.volumeDict[(round(la, 3), round(lo, 3))] * data
        return totalConc

    def calcMeanConcentration(self, sh, key='datax'):
        ####
        # Calculates the weighted by cell volume mean quatity in the dataset.
        totalConc = 0.
        totalVolume = 0.

        for la, lo, data in zip(self.loadFromshelve(sh, 'x_lat'),
                                self.loadFromshelve(sh, 'x_lon'), sh[key]):
            vol = self.volumeDict[(round(la, 3), round(lo, 3))]
            totalConc += vol * data
            totalVolume += vol
        if totalVolume == 0.: return 0.
        return totalConc / vol

    def setDictionaries(self, ):
        #####
        # setDictionaries: set up some dictionaries that are used throughout makePatternStatslots
        #

        # metricDicts: List of metrics to load for each shelve.
        self.metricsDict = {
            'rgam': 'Scale Ratio (Data/Reference)',
            'rE': 'Norm. Diff. Scale',
            'rE0': 'Norm. Bias',
            'rR': 'Spearman Corr.',
            'rP': 'p-value (robust)',
            #'rsig':'sign',
            'tgam': 'STD Ratio (Data/Reference)',
            'tE': 'Norm. Unbiased RMSD',
            'tE0': 'Norm. Bias',
            'tR': 'Pearson  Corr.',
            'tP': 'p-value (Taylor)',
            #'tsig':'sign',
            'N': 'Number',
            'b0': 'Intersect',
            'b1': 'Slope',
            'TotalModel': 'Total Model',
            'TotalInSitu': 'Total In Situ',
            'Model:In situ': 'Model / In situ ratio',
            'MeanModel': 'Weighted Model Mean',
            'MeanInSitu': 'Weighted InSitu Mean',
            'MNFB': 'MNFB',  # 'Mean norm. factor bias',
            'MNAFE': 'MNAFE',  # 'Mean norm. abs. factor error',
            'NMBF': 'NMBF',  # 'Norm. mean bias factor',		# robust
            'NMAEF': 'NMAEF',  # 'Norm. mean abs. error factor',	# robust
            'MedianModel': 'Median',
            'Model/obs. median': 'Model Median: data median'
        }
        # plot Styles: dictionary shoing Plot title and a list of things to plot in a subplot.
        self.plotStyles = {
            'Total N': [
                'TotalModel',
                'N',
            ],
            #'Total':   		['TotalModel','Model:In situ'],
            #'Weighted Mean':   	['MeanModel','Model:In situ'],
            'Median': ['MedianModel', 'Model/obs. median'],
            'Robust': [
                'rE',
                'rE0',
                'rR',
            ],
            'Taylor': [
                'tE',
                'tE0',
                'tR',
            ],
            #'Robust v Taylor Correlation':['rR','tR',],
            #'Robust v Taylor Correlation and N':['rR','tR','N',],
            #'Robust v Taylor Gamma':['rgam','tgam',],
            #'Robust v Taylor E':['rE','tE',],
            #'Robust v Taylor E0':['rE0','tE0',],
            'Linear Regression': [
                'b1',
                'b0',
                'tR',
            ],
            #'Yu MN metrics':['MNFB', 'MNAFE'],
            'Yu NM metrics': ['NMBF', 'NMAEF'],
            #'Yu metrics':['MNFB', 'MNAFE', 'NMBF', 'NMAEF'],
        }

        self.keyslongname = {m: getLongName(m) for m in self.keys}
        self.keyscolours = {
            m: i
            for m, i in zip(self.keys, [
                'r',
                'b',
                purple,
                'g',
                'orange',
                'CornflowerBlue',
                'DarkOrchid',
                'DarkTurquoise',
                'FireBrick',
                'LightSeaGreen',
                'Orchid',
            ])
        }

        self.AllShelves = []
        for key, shelves in list(self.shelveDict.items()):
            self.AllShelves.extend(shelves)

    def loadFromshelve(self, sh, sh_key):
        #####
        # Simple code to return mask if can't load the shelve
        try:
            return sh[sh_key]
        except:
            return np.ma.masked

    def loadMetrics(self, ):
        #####
        # Here, we load the metrics from disk into memory.
        # It produces a nested dictionary which is read out to make the plots.
        #
        print('loadMetrics')
        metrics = {me: {} for me in list(self.metricsDict.keys())}

        i = 0
        for xkey in self.xkeys:
            #	print 'first loop (xkeys):', xkey
            for me in list(self.metricsDict.keys()):
                metrics[me][xkey] = {}

            for key, shelves in list(self.shelveDict.items()):
                #print 'second loop (shelveDict):', 'xkey:',xkey,'key:',key ,'# shelves:',len(shelves)
                if not len(shelves):
                    print('No third loop as empty shelves:', xkey, key)
                    for met in list(self.metricsDict.keys()):
                        metrics[met][xkey][key] = np.ma.masked
                    continue

                for shelve in shelves:
                    #print 'third loop (shelves):', xkey, key, shelve

                    if type(xkey) == type('string'):
                        if shelve.find(xkey) < 0: continue
                    elif type(xkey) in [type(('', )),
                                        type([
                                            '',
                                        ])]:
                        continues = [shelve.find(xk) for xk in xkey]
                        if np.min(continues) == -1: continue

                    print("Found:", xkey, key, shelve)
                    if not os.path.exists(shelve):
                        for met in list(self.metricsDict.keys()):
                            metrics[met][xkey][key] = np.ma.masked
                        continue
                    sh = shOpen(shelve)

                    sig = self.loadFromshelve(sh,
                                              'robust.gamma') > 1 and 1 or -1

                    metrics['rgam'][xkey][key] = self.loadFromshelve(
                        sh, 'robust.gamma')
                    metrics['rE0'][xkey][key] = self.loadFromshelve(
                        sh, 'robust.E0')
                    metrics['rE'][xkey][key] = self.loadFromshelve(
                        sh, 'robust.E') * sig
                    metrics['rR'][xkey][key] = self.loadFromshelve(
                        sh, 'robust.R')
                    metrics['rP'][xkey][key] = self.loadFromshelve(
                        sh, 'robust.p')

                    sig = self.loadFromshelve(sh,
                                              'Taylor.gamma') > 1 and 1 or -1
                    metrics['tgam'][xkey][key] = self.loadFromshelve(
                        sh, 'Taylor.gamma')
                    metrics['tE0'][xkey][key] = self.loadFromshelve(
                        sh, 'Taylor.E0')
                    metrics['tE'][xkey][key] = self.loadFromshelve(
                        sh, 'Taylor.E') * sig
                    metrics['tR'][xkey][key] = self.loadFromshelve(
                        sh, 'Taylor.R')
                    metrics['tP'][xkey][key] = self.loadFromshelve(
                        sh, 'Taylor.p')

                    for k in ['N', 'b1', 'b0']:
                        metrics[k][xkey][key] = self.loadFromshelve(sh, k)

                    metrics['TotalModel'][xkey][
                        key] = self.calcTotalConcentration(sh,
                                                           key='datax') / 1.e12
                    metrics['TotalInSitu'][xkey][
                        key] = self.calcTotalConcentration(sh,
                                                           key='datay') / 1.e12
                    metrics['Model:In situ'][xkey][
                        key] = metrics['TotalModel'][xkey][key] / metrics[
                            'TotalInSitu'][xkey][key]

                    metrics['MeanModel'][xkey][
                        key] = self.calcMeanConcentration(sh, key='datax')
                    metrics['MeanInSitu'][xkey][
                        key] = self.calcMeanConcentration(sh, key='datax')

                    model = self.loadFromshelve(sh, 'datax')
                    obs = self.loadFromshelve(sh, 'datay')
                    if 'MNAFE' in list(sh.keys()):
                        for k in ['MNAFE', 'MNFB', 'NMAEF', 'NMBF']:
                            metrics[k][xkey][key] = self.loadFromshelve(sh, k)
                    else:
                        metrics['MNAFE'][xkey][key] = usm.MNAFE(model, obs)
                        metrics['MNFB'][xkey][key] = usm.MNFB(model, obs)
                        metrics['NMAEF'][xkey][key] = usm.NMAEF(model, obs)
                        metrics['NMBF'][xkey][key] = usm.NMBF(model, obs)

                    metrics['MedianModel'][xkey][key] = np.median(model)
                    metrics['Model/obs. median'][xkey][key] = np.median(model /
                                                                        obs)
                    sh.close()
                    i += 1
        self.i = i
        self.metrics = metrics

    def plotMetrics(self, ):
        ####
        # plotMetrics makes the plots.
        #
        return

        metrics = self.metrics
        title = getLongName(self.plotTitle)

        for plotType, yAxisKeys in list(self.plotStyles.items()):
            plotStyle = 'Lines'
            fn = self.filenamebase + self.plotTitle.replace(
                ' ', '') + '_' + plotType + '.png'
            fn = fn.replace(' ', '')
            if not ukp.shouldIMakeFile(self.AllShelves, fn, debug=False):
                continue
            #if os.path.exists(fn):continue

            fig = pyplot.figure()
            for d, metric in enumerate(yAxisKeys):
                ax = fig.add_subplot(len(yAxisKeys), 1, d + 1)
                xticks = []
                xticklabels = []

                linesDict = {
                    'x': [],
                    'xticks': [],
                    'xticklabels': [],
                    'TotalInSitu': [],
                    'MeanInSitu': [],
                }
                for key in self.keys:
                    linesDict[key] = []
                for o, xkey in enumerate(self.xkeys):

                    linesDict['x'].append(o)
                    linesDict['xticks'].append(o)
                    if type(xkey) in [type(('', )), type([
                            '',
                    ])]:
                        xkeystring = ' '.join(xkey)
                    else:
                        xkeystring = str(xkey)
                    xkeystring = xkeystring.replace('Ocean', '').replace(
                        'North',
                        'North ').replace('South',
                                          'South ')  #.replace(' ','\n')
                    linesDict['xticklabels'].append(xkeystring)

                    if not len(list(metrics[metric][xkey].keys())): continue
                    for m, val in sorted(metrics[metric][xkey].items()):
                        #if plotStyle=='Lines':
                        if val in [
                                None,
                        ]:
                            print('Warning:', 'metrics[', metric, '][', xkey,
                                  '] is empty:', [m, ':', val])
                            val = np.ma.masked
                        #print 'plotting:\t[d,metric]:', [d,metric],'\t[o,xkey]:',[o,xkey], '\t[m,val]:',[m,val]
                        linesDict[m].append(val)
                    linesDict['TotalInSitu'].append(
                        metrics['TotalInSitu'][xkey][m])
                    linesDict['MeanInSitu'].append(
                        metrics['MeanInSitu'][xkey][m])

                #print metric,linesDict.keys(),['x'], linesDict['TotalInSitu']
                if metric in [
                        'TotalModel',
                ]:
                    pyplot.plot(
                        linesDict['x'],
                        linesDict['TotalInSitu'],
                        c='k',
                        lw=2,
                    )  # assume in situ is the same for all
                if metric in [
                        'MeanModel',
                ]:
                    pyplot.plot(
                        linesDict['x'],
                        linesDict['MeanInSitu'],
                        c='k',
                        lw=2,
                    )  # assume in situ is the same for all
                if metric in ['Model:In situ', 'Model/obs. median']:
                    pyplot.plot(
                        linesDict['x'],
                        np.ones_like(linesDict['x']),
                        c='k',
                        lw=2,
                    )  # assume in situ is the same for all

                for key in self.keys:
                    #print linesDict['x'],linesDict[key],self.keyscolours[key]
                    print(key, linesDict['x'], linesDict[key])
                    pyplot.plot(linesDict['x'],
                                linesDict[key],
                                c=self.keyscolours[key],
                                lw=2)  # s = 40,lw=0.)#marker=modelmarkers[m]

                pyplot.ylabel(self.metricsDict[metric])

                if metric == yAxisKeys[-1]:
                    pyplot.xticks(
                        linesDict['xticks'],
                        linesDict['xticklabels'],
                        rotation='vertical',
                    )
                else:
                    pyplot.xticks(linesDict['xticks'], [
                        '',
                    ])

                if d + 1 == 1: ax.set_title(plotType + ' ' + title)

                if metric in [
                        'b0', 'rE0', 'tE0', 'tE', 'rO', 'rE', 'NMBF', 'NMAEF'
                ]:
                    pyplot.axhline(y=0., c='k', ls='--')
                if metric in [
                        'b1',
                        'rgam',
                        'tgam',
                        'rR',
                        'tR',
                ]:
                    pyplot.axhline(y=1., c='k', ls='--')

                labelCount = 0
                for key in self.keys:
                    pyplot.plot([], [],
                                label=self.keyslongname[key],
                                c=self.keyscolours[key],
                                lw=2)  #s = 40,lw=0.marker=modelmarkers[m],)
                    labelCount += 1
                if metric in [
                        'TotalModel',
                        'MeanModel',
                        'Model:In situ',
                        'N',
                        'Model/obs. median',
                ]:
                    pyplot.plot([], [], label='In Situ', c='k', lw=2)
                    labelCount += 1
            pyplot.subplots_adjust(bottom=0.20)
            if labelCount <= 4: ncol = 4
            if labelCount > 4: ncol = 3
            legend = pyplot.legend(
                loc='lower center',
                ncol=ncol,
                borderaxespad=0.,
                numpoints=1,
                scatterpoints=1,
                prop={'size': 8},
            )
            legend.draw_frame(False)
            legend.get_frame().set_alpha(0.)

            print('saving', fn)
            pyplot.savefig(
                fn,
                dpi=300,
            )


#def main():
#	kys=[]
#	kys.extend([o + 'Months' for o in hemispheres])
#kys.extend([o + 'Months' for o in oceans])
#kys.extend([o + 'Seasons' for o in oceans])
#kys.extend([o + 'Seasons' for o in oceans])
#	kys.extend(['months', 'oceans',])#'OceansSeasons','OceanMonths','hemispheres','seasons',
#	for xkeys in kys: #[,]:
#	    for self.keys in ['dmsmetrics', 'dmspmetrics',]: #
#		run(self.keys, xkeys,)

if __name__ == "__main__":
    #main()

    print('The end.')
