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
.. module:: paths
   :platform: Unix
   :synopsis: A list of paths to data files.
.. moduleauthor:: Lee de Mora <ledm@pml.ac.uk>
"""

import importlib
from socket import gethostname
#import UKESMpython as ukp
from getpass import getuser
import os
import sys


def folder(name):
    """ This snippet takes a string, makes the folder and the string.
            It also accepts lists of strings.
        """
    if type(name) == type(['a', 'b', 'c']):
        name = '/'.join(name, )
    if name[-1] != '/':
        name = name + '/'
    if os.path.exists(name) is False:
        os.makedirs(name)
        print('makedirs ', name)
    return name


# JASMIN
if gethostname().find('ceda.ac.uk') > -1 or gethostname().find(
        'jasmin') > -1 or gethostname().find('jc.rl.ac.uk') > -1:
    print("analysis-timeseries.py:\tBeing run at CEDA on ", gethostname())
    paths = importlib.import_module('bgcval2.Paths.paths_jasmin')
# MONSOON
elif gethostname().find('monsoon') > -1:
    print("Being run at the Met Office on ", gethostname())
    paths = importlib.import_module('bgcval2.Paths.paths_monsoon')
# PML
elif gethostname().find('pmpc') > -1:
    print("Paths.py:\tBeing run at PML on ", gethostname())
    paths = importlib.import_module('bgcval2.Paths.paths_pml')
# local
elif gethostname().find('valeriu-PORTEGE-Z30-C') > -1:
    print("Paths.py:\tBeing run at V laptop on ", gethostname())
    paths = importlib.import_module('bgcval2.Paths.paths_local')
# github actions
elif gethostname().find('-az') > -1:
    print("Paths.py:\tBeing run at GA machine ", gethostname())
    paths = importlib.import_module('bgcval2.Paths.paths_local')  # bogus
else:
    print("Got host name: ", gethostname())
    raise ValueError("Unidentified host. Run at either JASMIN, MONSOON or PML.")

machinelocation = paths.machinelocation
root_dir = paths.root_dir
#####
# Post processed Data location
shelvedir = paths.shelvedir

#####
# Post processed p2p Data location
p2p_ppDir = paths.p2p_ppDir

######
# Output location for plots.
imagedir = paths.imagedir

#####
# Location of model files.
ModelFolder_pref = paths.ModelFolder_pref

#####
# eORCA1 grid
orcaGridfn = paths.orcaGridfn

#####
# Location of data files.
ObsFolder = paths.ObsFolder
Dustdir = ObsFolder + "/MahowaldDust/"
WOAFolder_annual = ObsFolder + "WOA/annual/"
WOAFolder = ObsFolder + "WOA/"
DMSDir = ObsFolder + "/DMS_Lana2011nc/"
MAREDATFolder = ObsFolder + "/MAREDAT/MAREDAT/"
GEOTRACESFolder = ObsFolder + "/GEOTRACES/GEOTRACES_PostProccessed/"
GODASFolder = ObsFolder + "/GODAS/"
TakahashiFolder = ObsFolder + "/Takahashi2009_pCO2/"
MLDFolder = ObsFolder + "/IFREMER-MLD/"
iMarNetFolder = ObsFolder + "/LestersReportData/"
GlodapDir = ObsFolder + "/GLODAP/"
GLODAPv2Dir = ObsFolder + "/GLODAPv2/GLODAPv2_Mapped_Climatologies/"
OSUDir = ObsFolder + "OSU/"
CCIDir = ObsFolder + "CCI/"
icFold = ObsFolder + "/InitialConditions/"
