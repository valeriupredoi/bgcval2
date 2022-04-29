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
# modm@pml.ac.uk
#
"""
.. module:: robust 
   :platform: Unix
   :synopsis: Tool to make a robust statistics plots.
.. moduleauthor:: Momme Butenschon <momm@pml.ac.uk>

"""
from numpy import sqrt, arange, sin, cos, pi, abs, arccos, array, median, max
from scipy.stats import spearmanr
from matplotlib.pyplot import plot, axis, scatter, xlabel, ylabel, clabel, colorbar, text, subplot, tick_params, xticks, yticks, axhline, axvline
from .RobustStatistics import MAE, Qn, unbiasedMAE, IQR, Sn, MID
from matplotlib.ticker import MultipleLocator, NullFormatter

scaleFun = IQR
scaleDiff = MAE
scaleLable = "IQR"
diffLable = "MAE"


class StatsDiagram:

    def __init__(self, data, refdata, precision, *opts, **keys):
        """
        In the robust case a precision of the reference data is required.
        This represents the measurment precision and ensures that the scale measure can not be lesser than the precision. 
        In this case a warning is issued and the scale measure of the reference replaced by the precision.
        """
        self._stats(data, refdata, precision, *opts, **keys)

    def __call__(self, data, refdata, precision, *opts, **keys):
        self._stats(data, refdata, precision, *opts, **keys)

    def _stats(self, data, refdata, precision, *opts, **keys):
        dat = array(data).ravel()
        ref = array(refdata).ravel()
        self.scale = scaleFun(ref)
        if self.scale < precision:
            print("Reference scale lesser than measurment precision!")
            self.scale = precision
            print("\treplaced by precision", precision)
        bias = median(dat - ref)
        self.E0 = (median(dat - ref)) / self.scale
        self.gamma = scaleFun(dat) / self.scale
        self.R, self.p = spearmanr(dat, ref)
        self.E = scaleDiff(data - bias, ref) / self.scale
        print(self)

    def __str__(self):
        return "\tNormalised Bias: "+str(self.E0)+\
     "\n\tNormalised Difference Scale: "+str(self.E)+\
     "\n\tCorrelation Coefficient: "+str(self.R)+\
     "\n\t\t with p-value: "+str(self.p)+\
     "\n\tScale Ratio (Data/Reference): "+str(self.gamma)+\
     "\n\tReference Scale: "+str(self.scale)+'\n'


class Stats:

    def __init__(self, gam, E0, E, rho):
        self.E0 = E0
        self.R = rho
        self.gamma = gam
        self.E = E
        print(self)

    def __call__(self, gam, E0, E, rho):
        self.E0 = E0
        self.R = rho
        self.gamma = gam
        self.E = E
        print(self)

    def __str__(self):
        return "\tNormalised Bias: "+str(self.E0)+\
     "\n\tNormalised Difference Scale: "+str(self.E)+\
     "\n\tCorrelation Coefficient: "+str(self.R)+\
     "\n\tScale Ratio (Data/Reference): "+str(self.gamma)


class Target:

    def drawTargetGrid(self):
        #majL=MultipleLocator(1)
        #minL=MultipleLocator(.5)
        #majF=NullFormatter()
        ax = subplot(111)
        #plot(0.,0.,'ko',markersize=12)
        plot([-1., 1.], [0., 0.], c='k', ls='-')
        plot([0., 0.], [-1., 1.], c='k', ls='-')
        #axhline(y=0.,xmin=-1.,xmax=1.,c='k',ls='-',)
        #axvline(x=0.,ymin=-1.,ymax=1.,c='k',ls='-',)
        #axhline(y=0.,xmin=-2.,xmax=2.,c='k',ls='-',alpha=0.001)
        #axvline(x=0.,ymin=-2.,ymax=2.,c='k',ls='-',alpha=0.001)
        #plot(0.,0.,'ko',markersize=12)
        #ax.spines['left'].set_position('center')
        #ax.spines['bottom'].set_position('center')
        #ax.spines['right'].set_position('center')
        #ax.spines['top'].set_position('center')
        #tick_params(which='major',length=5,width=1.5)
        #tick_params(which='minor',length=3,width=1.2)
        #ax.xaxis.set_ticklabels([])
        #ax.yaxis.set_ticklabels([])
        #ax.xaxis.set_major_formatter(majF)
        #ax.yaxis.set_major_formatter(majF)
        #ax.xaxis.set_major_locator(majL)
        #ax.xaxis.set_minor_locator(minL)
        #ax.yaxis.set_major_locator(majL)
        #ax.yaxis.set_minor_locator(minL)
        #xticks((-1,1),("-1","1"))
        #yticks((-1,1),("-1","1"))
        #radius at minimum RMSD' given R:
        #  Mr=sqrt(1+R^2-2R^2)
        #for R=0.7:
        #plot(.71414284285428498*sin(a),.71414284285428498*cos(a),'k--')
        #radius at observ. uncertainty:
        #plot(.5*sin(a),.5*cos(a),'k:')
        #plot((0,),(0,),'k+')
        text(-.05,
             .5,
             '${Bias}/{IQR}_{ref}$',
             fontsize=16,
             transform=ax.transAxes,
             rotation=90,
             verticalalignment='center',
             horizontalalignment='center')
        text(.5,
             -.05,
             "${sign}({IQR}-{IQR}_{ref})*{MAE'}/{IQR}_{ref}$",
             fontsize=16,
             transform=ax.transAxes,
             verticalalignment='center',
             horizontalalignment='center')


class TargetDiagram(Target, Stats):

    def __init__(self,
                 gam,
                 E0,
                 E,
                 rho,
                 marker='o',
                 s=40,
                 antiCorrelation=False,
                 *opts,
                 **keys):
        Stats.__init__(self, gam, E0, E, rho)
        self.drawTargetGrid()
        if antiCorrelation:
            self._cmin = -1.
        else:
            self._cmin = 0.
        self._cmax = 1.
        self._lpos = []
        self.add(self.gamma,
                 self.E0,
                 self.E,
                 self.R,
                 marker=marker,
                 s=s,
                 *opts,
                 **keys)
        self.cbar = colorbar()
        self.cbar.set_label('Correlation Coefficient')

    def __call__(self, gam, E0, E, rho, marker='o', s=40, *opts, **keys):
        Stats.__call__(self, gam, E0, E, rho)
        self.add(self.gamma,
                 self.E0,
                 self.E,
                 self.R,
                 marker=marker,
                 s=s,
                 *opts,
                 **keys)

    def add(self, gam, E0, E, R, marker='o', s=40, *opts, **keys):
        sig = gam > 1 and 1 or -1
        scatter(sig * E,
                E0,
                c=R,
                vmin=self._cmin,
                vmax=self._cmax,
                marker=marker,
                s=s,
                *opts,
                **keys)
        self._lpos.append((sig * E, E0))
        rmax = abs(array(axis('scaled'))).max()
        rmax = max([rmax, 1.5])
        #plot((0,0),(-rmax,rmax),'k-')
        #plot((rmax,-rmax),(0,0),'k-')
        axis(xmin=-rmax, xmax=rmax, ymax=rmax, ymin=-rmax)

    def labels(self, lstr, *opts, **keys):
        rmax = abs(array(axis())).max()
        rmax = max([rmax, 1.5])
        for n, p in enumerate(self._lpos):
            text(p[0] + .025 * rmax, p[1] + .025 * rmax, lstr[n], *opts,
                 **keys)


class TargetStatistics(StatsDiagram, TargetDiagram):

    def __init__(self,
                 data,
                 refdata,
                 precision,
                 marker='o',
                 s=40,
                 antiCorrelation=False,
                 *opts,
                 **keys):
        StatsDiagram.__init__(self, data, refdata, precision, *opts, **keys)
        self.drawTargetGrid()
        if antiCorrelation:
            self._cmin = -1.
        else:
            self._cmin = 0.
        self._cmax = 1.
        self._lpos = []
        self.add(self.gamma,
                 self.E0,
                 self.E,
                 self.R,
                 marker=marker,
                 s=s,
                 *opts,
                 **keys)
        self.cbar = colorbar()
        self.cbar.set_label('Correlation Coefficient')

    def __call__(self,
                 data,
                 refdata,
                 precision,
                 marker='o',
                 s=40,
                 *opts,
                 **keys):
        StatsDiagram.__call__(self, data, refdata, precision, *opts, **keys)
        self.add(self.gamma,
                 self.E0,
                 self.E,
                 self.R,
                 marker=marker,
                 s=s,
                 *opts,
                 **keys)

    def add(self, gam, E0, E, R, marker='o', s=40, *opts, **keys):
        sig = gam > 1 and 1 or -1
        scatter(sig * E,
                E0,
                c=R,
                vmin=self._cmin,
                vmax=self._cmax,
                marker=marker,
                s=s,
                *opts,
                **keys)
        self._lpos.append((sig * E, E0))
        rmax = abs(array(axis('scaled'))).max()
        rmax = max([rmax, 1.5])
        axis(xmin=-rmax, xmax=rmax, ymax=rmax, ymin=-rmax)

    def labels(self, lstr, *opts, **keys):
        rmax = abs(array(axis())).max()
        for n, p in enumerate(self._lpos):
            text(p[0] + .025 * rmax, p[1] + .025 * rmax, lstr[n], *opts,
                 **keys)


#Examples:

#from numpy.random import randn
#a=randn(10)
#b=randn(10)
#ref=randn(10)
#subplot(221)
#TD=TargetDiagram(a,ref)
#TD(b,ref)
#subplot(222)
#TD=TaylorDiagram(a,ref)
#TD(b,ref)

#R1,p=pearsonr(a,ref)
#E1=a.mean()-ref.mean()
#scale1=a.scale()
#refscale1=ref.scale()
#R2,p=pearsonr(b,ref)
#E2=b.mean()-ref.mean()
#scale2=b.scale()
#refscale2=ref.scale()

#subplot(223)
#TD=TargetPlot(R1,E1,scale1,refscale1)
#TD(R2,E2,scale2,refscale2)
#subplot(224)
#TD=TaylorPlot(R1,E1,scale1,refscale1)
#TD(R2,E2,scale2,refscale2)
