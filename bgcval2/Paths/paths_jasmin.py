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

from socket import gethostname
#import UKESMpython as ukp
from getpass import getuser
import os


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


machinelocation = 'JASMIN'
root_dir = "/gws/nopw/j04/ukesm/"
#####
# Post processed Data location
shelvedir = folder(root_dir + "BGC_data/" + getuser() + "/shelves/")

#####
# Post processed p2p Data location
p2p_ppDir = folder(root_dir + "BGC_data/ukesm_postProcessed/")

######
# Output location for plots.
imagedir = folder(root_dir + 'BGC_data/' + getuser() + '/images/')

#####
# Location of model files.
ModelFolder_pref = folder(root_dir + "BGC_data/")

#####
# eORCA1 grid
orcaGridfn = '/gws/nopw/j04/esmeval/bgc-val/mesh_mask_eORCA1_wrk.nc'

#####
# Location of data files.
ObsFolder = "/gws/nopw/j04/esmeval/example_data/bgc/"
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
