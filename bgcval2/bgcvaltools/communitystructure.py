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
.. module:: UKESMpython
   :platform: Unix
   :synopsis: Toolkit for calculating the communicty structure.
.. moduleauthor:: Lee de Mora <ledm@pml.ac.uk>

"""

from numpy import zeros, logspace, log10, exp, min, max, arange, clip
from .comstrucFit import comstrucFit

# Taka's CSparameterisations:
pcclip = lambda a: clip(a * 100., 0., 100.)

microHirata2011 = lambda c: clip(
    100. / (0.9117 + exp(-2.733 * log10(c) + 0.4003)), 0., 100.)
picoHirata2011 = lambda c: clip(
    100. * (1. / (-1. * (0.1529 + exp(1.0306 * log10(c) - 1.5576))) - 1.8587 *
            log10(c) + 2.9954), 0., 100.)
nanoHirata2011 = lambda c: clip(100. - microHirata2011(c) - picoHirata2011(c),
                                0., 100.)
piconanoHirata2011 = lambda c: clip(100. - microHirata2011(c), 0., 100.)
CSparameterisations = [
    'hirata2011', "devred2011", "brewin2012", "brewin2011a", "brewin2010",
    "brewin2014"
]


def continuumModel(
    chl,
    pft,
    fit="brewin2010",
):
    if fit == "brewin2010":
        Cmpn = 1.057
        Spn = 0.851
        Cmp = 0.107
        Sp = 6.801
        CmpmSpm = 0.9
        CmpSp = 0.728

    if fit == "brewin2011a":
        Cmpn = 0.775
        Spn = 1.152
        Cmp = 0.146
        Sp = 5.118
        CmpmSpm = 0.8933
        CmpSp = 0.747

    if fit == "brewin2012":
        Cmpn = 0.937
        Spn = 1.033
        Cmp = 0.170
        Sp = 4.804
        CmpmSpm = 0.968
        CmpSp = 0.817

    if fit == "devred2011":
        Cmpn = 0.546
        Spn = 1.830
        Cmp = 0.148
        Sp = 6.765
        CmpmSpm = 1.
        CmpSp = 1.

    if fit == "brewin2014":
        Cmpn = 2.611
        Spn = 0.364
        Cmp = 0.730
        Sp = 1.047
        CmpmSpm = 0.951
        CmpSp = 0.763

    p = Cmp * (1. - exp(-Sp * chl))
    if pft == 'pico': return p
    pn = Cmpn * (1. - exp(-Spn * chl))
    if pft == 'piconano': return pn
    if pft == 'nano': return pn - p
    if pft == 'micro': return chl - pn
    return False


def chlPercent(
    chl,
    pft,
    fit='Hirata 2011',
):
    """ Takes Chl in [mg Chl / m^3], a pft, and a fit, then returns the %.
	"""
    p = pft.lower()
    fit = fit.lower().replace(' ', '')
    pfts = ['pico', 'nano', 'micro', 'piconano']

    if p not in pfts:
        print("chlPercent:\tERROR:\tWrong functional name", p, 'not in ', pfts)
        return False

    if fit not in CSparameterisations:
        print("chlPercent:\tERROR:\tWrong fit name", fit, 'not in ',
              CSparameterisations)
        return False

    if fit == 'hirata2011':
        if p == 'micro': return microHirata2011(chl)
        if p == 'pico': return picoHirata2011(chl)
        if p == 'nano': return nanoHirata2011(chl)
        if p == 'piconano': return nanoHirata2011(chl) + picoHirata2011(chl)
        return False
    else:
        return clip(100. * continuumModel(chl, pft, fit=fit) / chl, 0., 100.)


def test():
    """ Test the term chlPercent
	"""
    from itertoolsmodule import product
    pfts = ['pico', 'nano', 'micro', 'piconano']
    #CSparameterisations = ['Hirata2011',"Devred2011","Brewin2012","Brewin2011a","Brewin2010", "Brewin 2014",]
    chls = [0.01, 0.1, 0.5, 1., 1.5, 2., 10.]
    for f in CSparameterisations:
        for p, c in product(pfts, chls):
            print(f, p, c, ':\t', chlPercent(c, p, f))
