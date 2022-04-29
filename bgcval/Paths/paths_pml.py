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
import UKESMpython as ukp
from getpass import getuser

#####
# PML
machinelocation = ''

if gethostname().find('pmpc') > -1:
    print("Paths.py:\tBeing run at PML on ", gethostname())
    machinelocation = 'PML'

    #####
    # Post processed shelve Data location
    shelvedir = ukp.folder('shelves/')

    #####
    # Post processed p2p Data location
    p2p_ppDir = "/data/euryale7/scratch/ledm/ukesm_postProcessed/"

    ######
    # Output location for plots.
    imagedir = ukp.folder('images/')

    #####
    # Location of model files.
    ModelFolder_pref = "/data/euryale7/scratch/ledm/UKESM/MEDUSA/"

    #####
    # eORCA1 grid
    orcaGridfn = '/data/euryale7/scratch/ledm/UKESM/MEDUSA/mesh_mask_eORCA1_wrk.nc'

    #####
    # Location of data files.
    ObsFolder = "/data/euryale7/backup/ledm/Observations/"
    CCIDir = ObsFolder + "CCI/"
    Dustdir = ObsFolder + "/MahowaldDust/"
    DMSDir = ObsFolder + "/DMS_Lana2011nc/"
    GEOTRACESFolder = ObsFolder + "/GEOTRACES/GEOTRACES_PostProccessed/"
    GlodapDir = ObsFolder + "/GLODAP/"
    GLODAPv2Dir = ObsFolder + "/GLODAPv2/GLODAPv2_Mapped_Climatologies/"
    GODASFolder = ObsFolder + "/GODAS/clim/"
    iMarNetFolder = ObsFolder + "/LestersReportData/"
    MAREDATFolder = ObsFolder + "/MAREDAT/MAREDAT/"
    MLDFolder = ObsFolder + "/IFREMER-MLD/"
    OSUDir = ObsFolder + "OSU/"
    TakahashiFolder = ObsFolder + "/Takahashi2009_pCO2/"
    WOAFolder_annual = ObsFolder + "WOA/annual/"
    WOAFolder = ObsFolder + "WOA/"

    icFold = "/data/euryale7/backup/ledm/UKESM/InitialConditions/"
