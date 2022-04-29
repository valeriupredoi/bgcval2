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
.. module:: StatsDiagrams 
   :platform: Unix
   :synopsis: Tool to make a Statistical diagrams.
.. moduleauthor:: Momme Butenschon <momm@pml.ac.uk>

"""

from numpy import sqrt, arange, sin, cos, pi, abs, arccos, array, max as npmax, abs as npabs
from scipy.stats import pearsonr
from matplotlib.pyplot import plot, axis, scatter, xlabel, ylabel, clabel, colorbar, text, subplot

rmsds = lambda gamma, R: sqrt(1. + gamma**2 - 2. * gamma * R)


class StatsDiagram:
    """Base class for statistical summary diagrams based on two arrays of
    the same size that will be compared on a point to point base.
    The first array is considered the data to be evaluated, the second is 
    the reference data. It computes all the basic metrics using the 
    ``self._stats`` function:
    
          :self.std: the standard deviation of the reference data
          :self.E0: the mean bias of the two data sets
          :self.gama: the ratio of the std of data over reference data
          :self.R, self.p: the Pearson correlation with p-value
          :self.E: the root-mean square difference of the two dataset
    """

    def __init__(self, data, refdata, *opts, **keys):
        self._stats(data, refdata, *opts, **keys)

    def __call__(self, data, refdata, *opts, **keys):
        """Recomputes metrics for new data."""
        self._stats(data, refdata, *opts, **keys)

    def _stats(self, data, refdata, *opts, **keys):
        dat = array(data).ravel()
        ref = array(refdata).ravel()
        self.std = ref.std()
        self.E0 = (dat - ref).mean() / self.std
        self.gamma = dat.std() / self.std
        self.R, self.p = pearsonr(dat, ref)
        self.E = rmsds(self.gamma, self.R)
        print(self)

    def __str__(self):
        return "\tNormalised Bias: "+str(self.E0)+\
     "\n\tNormalised Unbiased RMSD: "+str(self.E)+\
     "\n\tNormalised RMSD: "+str(sqrt(self.E0**2+self.E**2))+\
     "\n\tCorrelation Coefficient: "+str(self.R)+\
     "\n\t\t with p-value: "+str(self.p)+\
     "\n\tSTD Ratio (Data/Reference): "+str(self.gamma)+\
     "\n\tReference STD: "+str(self.std)+'\n'


class Stats:
    """Base class for statistical summary diagrams, using precalculated metrics rather than the full datasets as input."""

    def __init__(self, gam, E0, rho):
        """Loading the necessary metrics."""
        self.E0 = E0
        self.R = rho
        self.gamma = gam
        self.E = rmsds(gam, rho)
        print(self)

    def __call__(self, gam, E0, rho):
        """Reloading the necessary metrics."""
        self.E0 = E0
        self.R = rho
        self.gamma = gam
        self.E = rmsds(gam, rho)
        print(self)

    def __str__(self):
        return "\tNormalised Bias: "+str(self.E0)+\
     "\n\tNormalised Unbiased RMSD: "+str(self.E)+\
     "\n\tNormalised RMSD: "+str(sqrt(self.E0**2+self.E**2))+\
     "\n\tCorrelation Coefficient: "+str(self.R)+\
     "\n\tSTD Ratio (Data/Reference): "+str(self.gamma)


class Target:
    """Base class providing a function to draw the grid for Target Diagrams.    """

    def drawTargetGrid(self):
        a = arange(0, 2.01 * pi, .02 * pi)
        plot(sin(a), cos(a), 'k')
        #radius at minimum RMSD' given R:
        #  Mr=sqrt(1+R^2-2R^2)
        #for R=0.7:
        plot(.71414284285428498 * sin(a), .71414284285428498 * cos(a), 'k--')
        #radius at observ. uncertainty:
        #plot(.5*sin(a),.5*cos(a),'k:')
        plot((0, ), (0, ), 'k+')
        ylabel('${Bias}/\sigma_{ref}$', fontsize=16)
        xlabel('${sign}(\sigma-\sigma_{ref})*{RMSD}\'/\sigma_{ref}$',
               fontsize=16)


class Taylor:
    """Base class providing a function to draw the grid for Taylor Diagrams.    """

    def drawTaylorGrid(self, R, dr, antiCorrelation):
        if antiCorrelation:
            a0 = -1.
        else:
            a0 = 0.
#Draw circles:
        a = arange(a0, 1.01, .05) * .5 * pi
        self.ax = plot(sin(a), cos(a), 'k-')
        plot(.5 * sin(a), .5 * cos(a), 'k:')
        n = R / .5
        for m in arange(3, n + 1):
            plot(m * .5 * sin(a), m * .5 * cos(a), 'k:')

#Draw rays for correlations at .99,.75,.5,.25 steps:
        rays = list(arange(a0, 1.05, .25))
        rays.append(.99)
        if a0 == -1.: rays.append(-.99)
        for rho in rays:
            plot((R * rho, 0), (R * sin(arccos(rho)), 0), 'k:')
            d = rho >= 0. and 1.02 or 1.25
            text(d * R * rho,
                 1.01 * R * sin(arccos(rho)),
                 str(rho),
                 fontsize=8)
        text(1.01 * R * sin(pi * .25),
             1.02 * R * cos(pi * .25),
             r'$\rho$',
             rotation=-45,
             fontsize=16)
        #text(0.,1.02*R*cos(0.),'0')
        #text(1.03*R*sin(.5*pi),0.,'1')
        goOn = True
        r = dr
        a = arange(-1., 1.01, .05) * .5 * pi
        plot((1, ), (0, ), 'ko')
        while goOn:
            xx = []
            yy = []
            for p in a:
                x = r * sin(p) + 1.
                y = r * cos(p)
                if antiCorrelation:
                    if x**2 + y**2 < R**2:
                        xx.append(x)
                        yy.append(y)
                else:
                    if x**2 + y**2 < R**2:
                        if x > 0.:
                            xx.append(x)
                            yy.append(y)
            if len(xx) > 0:
                plot(xx, yy, 'k--')
            else:
                goOn = False
            r += dr
        xlabel('$\sigma/\sigma_{ref}$', fontsize=16)
        ylabel('$\sigma/\sigma_{ref}$', fontsize=16)


class TaylorDiagram(Taylor, Stats):
    """Class for drawing a Taylor diagram using the pre-calculated metrics."""

    def __init__(self,
                 gam,
                 E0,
                 rho,
                 R=2.5,
                 dr=.5,
                 antiCorrelation=True,
                 marker='o',
                 s=40,
                 *opts,
                 **keys):
        """Initialises the class given the pre-calculated metrics and draws 
        the diagram grid with the first point.
        Markers in the diagram are colour-coded using the mean bias.

             :antiCorrelation: show negative correlations
             :marker: shape used to show points, should be hollow shape to
                      be filled with colour code.
        """
        print("StatsDiagram:\tTaylorDiagram:\t", gam, E0, rho, R)
        Stats.__init__(self, gam, E0, rho)
        R = max(int(2. * self.gamma + 1.) / 2., 1.5)
        self.drawTaylorGrid(R, dr, antiCorrelation)
        self._cmax = max(1., abs(self.E0))
        self._cmin = -self._cmax
        self._lpos = []
        if antiCorrelation:
            self._axis = {
                'xmin': -1.3 * R,
                'xmax': 1.3 * R,
                'ymin': -.1 * R,
                'ymax': 1.1 * R
            }
        else:
            self._axis = {
                'xmin': -.1 * R,
                'xmax': 1.3 * R,
                'ymin': -.1 * R,
                'ymax': 1.1 * R
            }
        self.add(self.gamma,
                 self.E0,
                 self.R,
                 marker=marker,
                 s=s,
                 *opts,
                 **keys)
        self.cbar = colorbar()
        self.cbar.set_label('${Bias}/\sigma_{ref}$')

    def __call__(self, gam, E0, rho, marker='o', s=40, *opts, **keys):
        """Adds points to the diagram adjusting the colour codes.

             :antiCorrelation: show negative correlations
             :marker: shape used to show points, should be hollow shape to
                      be filled with colour code.
        """
        self._cmax = max(abs(E0), self._cmax)
        self._cmin = -self._cmin
        Stats.__call__(self, gam, E0, rho)
        self.add(gam, E0, rho, marker=marker, s=s, *opts, **keys)

    def add(self, gam, E0, R, marker='o', s=40, *opts, **keys):
        """Function to add additional points to the diagram, using invoked 
        by means of the ``__call__`` function."""
        self._cmax = max(1., abs(E0))
        E = rmsds(gam, R)
        scatter(gam * R,
                gam * sin(arccos(R)),
                c=E0,
                vmin=self._cmin,
                vmax=self._cmax,
                marker=marker,
                s=s,
                *opts,
                **keys)
        self._lpos.append((gam * R, gam * sin(arccos(R))))
        #      rmax=abs(array(axis('scaled'))).max()
        # axis(xmin=-rmax,xmax=rmax,ymax=rmax,ymin=-0.1)
        axis(**self._axis)

    def labels(self, lstr, *opts, **keys):
        """Adds labels ``lstr``` to the points in the diagram"""
        yrange = axis()[2:]
        rmax = abs(yrange[1] - yrange[0])
        for n, p in enumerate(self._lpos):
            text(p[0] + .025 * rmax, p[1] + .025 * rmax, lstr[n], *opts,
                 **keys)


absmax = lambda x: npmax(npabs(x))


class TaylorDiagramMulti(Taylor, Stats):

    def __init__(self,
                 gam_s,
                 E0_s,
                 rho_s,
                 R=2.5,
                 dr=.5,
                 antiCorrelation=True,
                 markers='o',
                 s=40,
                 *opts,
                 **keys):
        print("StatsDiagram:\tTaylorDiagram:\t", gam_s, E0_s, rho_s, R)
        Stats.__init__(self, absmax(gam_s), absmax(E0_s), absmax(rho_s))
        R = max(int(2. * self.gamma + 1.) / 2., 1.5)
        self.drawTaylorGrid(R, dr, antiCorrelation)
        self._cmax = max(1., int(absmax(E0_s)) + 1)
        self._cmin = -self._cmax
        self._lpos = []
        if antiCorrelation:
            self._axis = {
                'xmin': -1.3 * R,
                'xmax': 1.3 * R,
                'ymin': -.1 * R,
                'ymax': 1.1 * R
            }
        else:
            self._axis = {
                'xmin': -.1 * R,
                'xmax': 1.3 * R,
                'ymin': -.1 * R,
                'ymax': 1.1 * R
            }
        for i, gam in enumerate(gam_s):
            self.add(gam_s[i],
                     E0_s[i],
                     rho_s[i],
                     marker=markers[i],
                     s=s,
                     *opts,
                     **keys)
        self.cbar = colorbar()
        self.cbar.set_label('${Bias}/\sigma_{ref}$')

    def __call__(self, gam, E0, rho, marker='o', s=40, *opts, **keys):
        #self._cmax=max(abs(E0),self._cmax)
        self._cmin = -self._cmin
        Stats.__call__(self, gam, E0, rho)
        self.add(gam, E0, rho, marker=marker, s=s, *opts, **keys)

    def add(self, gam, E0, R, marker='o', s=40, *opts, **keys):
        #self._cmax=max(1.,abs(E0))
        E = rmsds(gam, R)
        scatter(gam * R,
                gam * sin(arccos(R)),
                c=E0,
                vmin=self._cmin,
                vmax=self._cmax,
                marker=marker,
                s=s,
                *opts,
                **keys)
        self._lpos.append((gam * R, gam * sin(arccos(R))))
        #      rmax=abs(array(axis('scaled'))).max()
        # axis(xmin=-rmax,xmax=rmax,ymax=rmax,ymin=-0.1)
        axis(**self._axis)

    def labels(self, lstr, *opts, **keys):
        yrange = axis()[2:]
        rmax = abs(yrange[1] - yrange[0])
        for n, p in enumerate(self._lpos):
            text(p[0] + .025 * rmax, p[1] + .025 * rmax, lstr[n], *opts,
                 **keys)


class TargetDiagram(Target, Stats):

    def __init__(self,
                 gam,
                 E0,
                 rho,
                 marker='o',
                 s=40,
                 antiCorrelation=False,
                 *opts,
                 **keys):
        """Initialises the class given the pre-calculated metrics and draws 
        the diagram grid with the first point.
        Markers in the diagram are colour-coded using the mean bias.

             :antiCorrelation: show negative correlations
             :marker: shape used to show points, should be hollow shape to
                      be filled with colour code.
        """
        Stats.__init__(self, gam, E0, rho)
        self.drawTargetGrid()
        if antiCorrelation:
            self._cmin = -1.
        else:
            self._cmin = 0.
        self._cmax = 1.
        self._lpos = []
        self.add(self.gamma,
                 self.E0,
                 self.R,
                 marker=marker,
                 s=s,
                 *opts,
                 **keys)
        self.cbar = colorbar()
        self.cbar.set_label('Correlation Coefficient')

    def __call__(self, gam, E0, rho, marker='o', s=40, *opts, **keys):
        """Adds points to the diagram adjusting the colour codes.

             :antiCorrelation: show negative correlations
             :marker: shape used to show points, should be hollow shape to
                      be filled with colour code.
        """
        Stats.__call__(self, gam, E0, rho)
        self.add(self.gamma,
                 self.E0,
                 self.R,
                 marker=marker,
                 s=s,
                 *opts,
                 **keys)

    def add(self, gam, E0, R, marker='o', s=40, *opts, **keys):
        """Function to add additional points to the diagram, using invoked 
        by means of the ``__call__`` function."""
        sig = gam > 1 and 1 or -1
        E = rmsds(gam, R)
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
        #print 'rmax: 1',rmax, array(axis('scaled'))
        #plot((0,0),(-rmax,rmax),'k-')
        #print 'rmax: 2',rmax , array(axis('scaled'))
        #plot((rmax,-rmax),(0,0),'k-')
        #print 'rmax: 3',rmax  , array(axis('scaled'))
        #print "add: axis:", rmax , array(axis('scaled')), abs(array(axis('scaled'))).max(),'cmax:=',[self._cmin,self._cmax]
        axis(xmin=-rmax, xmax=rmax, ymax=rmax, ymin=-rmax)
        #print 'rmax: 4',rmax, array(axis('scaled'))
    def labels(self, lstr, *opts, **keys):
        """Adds labels ``lstr``` to the points in the diagram"""
        rmax = abs(array(axis())).max()
        for n, p in enumerate(self._lpos):
            text(p[0] + .025 * rmax, p[1] + .025 * rmax, lstr[n], *opts,
                 **keys)


class TargetStatistics(StatsDiagram, TargetDiagram):

    def __init__(self,
                 data,
                 refdata,
                 marker='o',
                 s=40,
                 antiCorrelation=False,
                 *opts,
                 **keys):
        """Initialises the class computing all necessary metrics and draws 
        the diagram grid with the first point.
        Markers in the diagram are colour-coded using the mean bias.

             :antiCorrelation: show negative correlations
             :marker: shape used to show points, should be hollow shape to
                      be filled with colour code.
        """
        StatsDiagram.__init__(self, data, refdata, *opts, **keys)
        self.drawTargetGrid()
        if antiCorrelation:
            self._cmin = -1.
        else:
            self._cmin = 0.
        self._cmax = 1.
        self._lpos = []
        self.add(self.gamma,
                 self.E0,
                 self.R,
                 marker=marker,
                 s=s,
                 *opts,
                 **keys)
        self.cbar = colorbar()
        self.cbar.set_label('Correlation Coefficient')

    def __call__(self, data, refdata, marker='o', s=40, *opts, **keys):
        """Adds points to the diagram adjusting the colour codes.

             :antiCorrelation: show negative correlations
             :marker: shape used to show points, should be hollow shape to
                      be filled with colour code.
        """
        StatsDiagram.__call__(self, data, refdata, *opts, **keys)
        self.add(self.gamma,
                 self.E0,
                 self.R,
                 marker=marker,
                 s=s,
                 *opts,
                 **keys)

    def add(self, gam, E0, R, marker='o', s=40, *opts, **keys):
        sig = gam > 1 and 1 or -1
        E = sqrt(1. + gam**2 - 2. * gam * R)
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
        plot((0, 0), (-rmax, rmax), 'k-')
        plot((rmax, -rmax), (0, 0), 'k-')
        axis(xmin=-rmax, xmax=rmax, ymax=rmax, ymin=-rmax)

    def labels(self, lstr, *opts, **keys):
        rmax = abs(array(axis())).max()
        for n, p in enumerate(self._lpos):
            text(p[0] + .025 * rmax, p[1] + .025 * rmax, lstr[n], *opts,
                 **keys)


class TaylorStatistics(StatsDiagram, TaylorDiagram):

    def __init__(self,
                 data,
                 refdata,
                 R=2.5,
                 dr=.5,
                 antiCorrelation=True,
                 marker='o',
                 s=40,
                 *opts,
                 **keys):
        """Initialises the class computing all necessary metrics and draws 
        the diagram grid with the first point.
        Markers in the diagram are colour-coded using the mean bias.

             :antiCorrelation: show negative correlations
             :marker: shape used to show points, should be hollow shape to
                      be filled with colour code.
        """
        StatsDiagram.__init__(self, data, refdata, *opts, **keys)
        R = max(int(2. * self.gamma + 1.) / 2., 1.5)
        self.drawTaylorGrid(R, dr, antiCorrelation)
        self._cmax = max(1., abs(self.E0))
        self._cmin = -self._cmax
        self._lpos = []
        if antiCorrelation:
            self._axis = {
                'xmin': -1.3 * R,
                'xmax': 1.3 * R,
                'ymin': -.1 * R,
                'ymax': 1.1 * R
            }
        else:
            self._axis = {
                'xmin': -.1 * R,
                'xmax': 1.3 * R,
                'ymin': -.1 * R,
                'ymax': 1.1 * R
            }
        self.add(self.gamma,
                 self.E0,
                 self.R,
                 marker=marker,
                 s=s,
                 *opts,
                 **keys)
        self.cbar = colorbar()
        self.cbar.set_label('${Bias}/\sigma_{ref}$')

    def __call__(self, data, refdata, marker='o', s=40, *opts, **keys):
        """Adds points to the diagram adjusting the colour codes.

             :antiCorrelation: show negative correlations
             :marker: shape used to show points, should be hollow shape to
                      be filled with colour code.
        """
        StatsDiagram.__call__(self, data, refdata, *opts, **keys)
        self._cmax = max(abs(self.E0), self._cmax)
        self._cmin = -self._cmax
        self.add(self.gamma,
                 self.E0,
                 self.R,
                 marker=marker,
                 s=s,
                 *opts,
                 **keys)


#Examples:

#from StatsDiagram import *
#from numpy.random import randn
#from matplotlib.pyplot import show,subplot

#a=randn(10)
#b=randn(10)
#ref=randn(10)
#subplot(221)
#TD=TargetStatistics(a,ref)
#TD(b,ref)
#subplot(222)
#TD=TaylorStatistics(a,ref)
#TD(b,ref)

#std1=a.std()
#std2=b.std()
#refstd=ref.std()
#R1,p=pearsonr(a,ref)
#E1=(a.mean()-ref.mean())/refstd
#G1=std1/refstd
#R2,p=pearsonr(b,ref)
#E2=(b.mean()-ref.mean())/refstd
#G2=std2/refstd

#subplot(223)
#TayD=TargetDiagram(G1,E1,R1,)
#TayD(G2,E2,R2,)
#subplot(224)
#TarD=TaylorDiagram(G1,E1,R1,)
#TarD(G2,E2,R2,)

#show()
