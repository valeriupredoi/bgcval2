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
.. module:: comstrucFit
   :platform: Unix
   :synopsis: This code produces a communistruture fit to data, using least squares regression.
.. moduleauthor:: Lee de Mora <ledm@pml.ac.uk>

"""

from scipy.optimize import curve_fit
from numpy import ma, exp, arange, logspace, isnan, isinf


class comstrucFit:
    """	This code produces a fit to data, using least squares regression.
  	Typically run in the following way:
  	fy = comstrucFit(totalchl,planktonFunctionTypePercent, pft='pico').plotForX(chlAxis)
  """

    def Cp(self, tchl, Cmp, Sp):
        return ma.array(Cmp * (1. - exp(-Sp * tchl)) / tchl)  #pico

    def Cpn(self, tchl, Cmpn, Spn):
        return ma.array(Cmpn * (1. - exp(-Spn * tchl)) / tchl)  #piconano

    #def pnmp(self,tchl, Cmp,  Sp, Cmpn, Spn):	return self.pn(tchl, Cmpn, Spn) - self.p(tchl, Cmp,  Sp)#nano
    #def cmpn(self,tchl, Cmpn, Spn):		return tchl - self.pn(	tchl, Cmpn, Spn)  	  	#micro

    def microHirata2011(self, tchl, a0, a1, a2):
        return ma.aray(1. / (a0 + exp(a1 * log10(tchl) + a2)))

    def picoHirata2011(self, tchl, a0, a1, a2, a3, a4):
        return ma.array(1. / (-1. * (a0 + exp(a1 * log10(c) + a2))) +
                        a3 * log10(c) + a4)

    #  def microHirata2011 = lambda c: 1./(0.912 + exp(-2.733 * log10(c) + 0.4)))
    # def picoHirata2011  = lambda c: (1./(-1.*(0.153 + exp(1.031 * log10(c) -1.558))) - 1.860*log10(c)+2.995))
    # def nanoHirata2011  = lambda c: 1. - microHirata2011(c) - picoHirata2011(c))
    #def piconanoHirata2011 = lambda c: 1. - microHirata2011(c))

    def __init__(
        self,
        tchl,
        pftfrac,
        pft='pico',
        chlRange=[0, 0],
        pcRange=[-4., 104.],
        fitTypes='Continuum',
    ):

        if chlRange[0] == chlRange[1]:
            chlRange = [tchl.min(), tchl.max()]
            print("Community Structure Fit:\tSetting range:", chlRange)
        self.chlRange = chlRange
        self.pcRange = pcRange
        self.fitTypes = fitTypes
        self.pft = pft

        tchl = ma.array(tchl)
        pftfrac = ma.array(pftfrac)
        tchl = ma.masked_outside(tchl, chlRange[0], chlRange[1])
        pftfrac = ma.masked_outside(pftfrac, pcRange[0], pcRange[1])
        tchl = ma.masked_where(tchl.mask + pftfrac.mask, tchl)
        pftfrac = ma.masked_where(tchl.mask + pftfrac.mask, pftfrac)

        self.tchl = tchl.compressed()
        self.pftfrac = pftfrac.compressed()
        if fitTypes == 'runningMean':
            self.tchl, self.pftfrac = self.doRunningMean()

        if len(self.tchl) != len(self.pftfrac):
            print("Community Structure Fit:\tERROR:\tTotal Chl and ", self.pft,
                  " % do not match", len(self.tchl), '!=', len(self.pftfrac))
            assert False

        if self.fitTypes == 'Continuum':
            if self.pft == 'pico': self.func = self.Cp
            if self.pft == 'piconano': self.func = self.Cpn
        if self.fitTypes == 'Hirata':
            print("WARNING: This is not tested.")
            if self.pft == 'pico': self.func = self.picoHirata2011
            #if self.pft=='piconano':self.func = self.piconanoHirata2011
            #if self.pft=='nano':	self.func = self.nanoHirata2011
            if self.pft == 'micro': self.func = self.microHirata2011

        #if self.pft not in ['pico', 'piconano',]:#'micro', 'nano',]:
        #	print "Community Structure Fit:\t not performing fit to ",self.pft

        if 'NoFit' not in [
                fitTypes,
        ]: self.fit()

    def doRunningMean(self, nbins=1001):
        """ Bin the two 1D arrays, to simplyfy the fit.
  	"""

        print("Running Mean:", self.chlRange[0], ma.log10(self.chlRange[0]),
              self.chlRange[1], ma.log10(self.chlRange[1]))
        x_axis = logspace(ma.log10(self.chlRange[0]),
                          ma.log10(self.chlRange[1]), nbins)

        xbins = list(zip(x_axis[:-1], x_axis[1:]))

        data = {}
        for x in xbins:
            y = ma.masked_where((self.tchl <= x[0]) + (self.tchl > x[1]),
                                self.pftfrac).compressed()
            x = ma.mean(x)
            if y is ma.masked or not len(y):
                data[x] = ma.masked
                continue
            #data[x] = ma.masked
            print(x, y.mean(), len(y))
            if len(y) < len(self.tchl) / 1000.:
                data[x] = ma.masked
            else:
                data[x] = ma.mean(y)

        tchl = ma.array([i for i in sorted(data.keys())])
        pftfrac = ma.array([data[i] for i in sorted(data.keys())])

        return tchl, pftfrac
        #print "tchl:", self.tchl[0], self.tchl.mean(), self.tchl.min(), self.tchl.max()
        #print "pftfrac:", self.pftfrac[0], self.pftfrac.mean(), self.pftfrac.min(), self.pftfrac.max()

    def fit(self):
        print("Community Structure Fit:\t performing fit", self.pft,
              self.fitTypes)

        try:
            popt, pcov = curve_fit(
                self.func,
                self.tchl,
                self.pftfrac,
            )
        except:
            print("Unable to solve and find a fit")
            popt = ma.masked_all(2)
            pcov = ma.masked_all((2, 2))

#if self.pft=='pico':	popt, pcov = curve_fit(self.p,    self.tchl, self.pftfrac,)

#if self.pft=='piconano':popt, pcov = curve_fit(self.pn,  self.tchl, self.pftfrac,)
#if self.pft=='nano':	popt, pcov = curve_fit(self.pnmp, self.tchl, self.pftfrac,)
#if self.pft=='micro':	popt, pcov = curve_fit(self.cmpn, self.tchl, self.pftfrac,)

        print("Community Structure Fit:\t fit results:", popt, pcov)
        if ma.is_masked(pcov) or isnan(pcov).any() or isinf(pcov).any():
            print(
                "Community Structure Fit:\tWARNING:\tfit results contains a masked Value:",
                popt, pcov)

        self.popt, self.pcov = popt, pcov

        return popt, pcov

    def plotForX(self, x):
        if self.pft not in [
                'pico',
                'piconano',
                'micro',
                'nano',
        ]:
            print("Community Structure Fit:\t not performing fit to ",
                  self.pft)
            z = ma.zeros(len(x))
            return ma.masked_where(z == 0, z)
        if ma.is_masked(self.pcov) or isnan(self.pcov).any() or isinf(
                self.pcov).any():
            print(
                "Community Structure Fit:\tWARNING:\tfit results contains a masked Value:",
                self.popt, self.pcov)
            z = ma.zeros(len(x))
            return ma.masked_where(z == 0, z)

        print("Community Structure Fit:\t plot for x", self.pft, self.popt)
        y = []
        x = ma.array(x)
        if x.min() < self.chlRange[0] or x.max() > self.chlRange[1]:
            print("Community Structure Fit:\tWarning:\tplot for x ", self.pft,
                  " x ouside Chl Range", [x.min(), x, max()], '<',
                  self.chlRange)

    #if:
    #	print "Community Structure Fit:\tWarning:\tplot for x ", self.pft," x > Chl Range", x.max(), '>',self.chlRange

        if self.pft in ['pico', 'piconano', 'micro']:
            y = self.func(x, self.popt[0], self.popt[1])
        if self.pft in [
                'nano',
        ]:
            y = self.func(x, self.popt[0], self.popt[1], self.popt[2],
                          self.popt[3])

        #if self.pft=='pico':
        #if self.pft=='piconano':y = 	self.pn(x, self.popt[0],self.popt[1])
        #if self.pft=='nano':	y = 	self.p(x, self.popt[0],self.popt[1])
        #if self.pft=='micro':	y = 	self.p(x, self.popt[0],self.popt[1])

        if len(y) != len(x):
            print("Community Structure Fit:\tERROR:\tplot for x ", self.pft,
                  " % do not match", len(y), '!=', len(x))
            assert False

        if not len(y):
            print(
                "Community Structure Fit:\tERROR:\tCould not produce y in plotForX"
            )

        else:
            return y

        #if self.pft=='piconano':	popt, pcov = curve_fit(self.ppn,   self.tchl, self.pftfrac,)
        #if self.pft=='nano':		popt, pcov = curve_fit(self.ppnmp, self.tchl, self.pftfrac,)
        assert False

if __name__ == "__main__":
    #test()
    #def test():
    from ncdfView import ncdfView
    fn = "/data/euryale7/scratch/ledm/iMarNet/xhonp/MEANS/xhonpo_20011231m01P.nc"
    print("loading", fn)
    nc = ncdfView(fn, Quiet=True)
    chl = nc('Chl3')[0, 0]
    tchl = nc('Chl1')[0, 0] + nc('Chl2')[0, 0] + nc('Chl3')[0,
                                                            0] + nc('Chl4')[0,
                                                                            0]
    pftfrac = chl / tchl
    pft = 'pico'

    popt, pcov = comstrucFit(tchl, pftfrac, pft=pft)

    print('The end.')
