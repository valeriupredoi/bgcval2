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

from socket import gethostname
import UKESMpython as ukp
from getpass import getuser

#####
# LOCAL_MACHINE_NAME
machinelocation = ''
"""
	This is a set of paths that is used by BGC_Val to locate various datasets.
	The idea is that it is set once per machine, and then used in multiple analyses.
	It also sets the paths for storing the post processed data. 
	
	The specific paths are:
	:shelvedir: A location for post processed data, saved in pythons shelve format.
	:esmvalFolder: A path for the model data.
	:ObsFolder: The base folder for many observational datasets. These will vary depending on your path location.
	
"""

#####
# LOCAL_MACHINE_NAME
if gethostname().find('LOCAL_MACHINE_NAME') > -1:
    print("analysis-timeseries.py:\tBeing run at LOCAL_MACHINE_NAME on ", gethostname(
    ))
    machinelocation = 'LOCAL_MACHINE_NAME'

    #####
    # Post processed Data location
    shelvedir = ukp.folder("/Path/To/Model/Data/BGC_val_data/" + getuser() +
                           "/shelves/")

    #####
    # Location of model files.
    esmvalFolder = "/Path/To/Model/Data/BGC_data/"
    ModelFolder_pref = ukp.folder(esmvalFolder)

    #####
    # eORCA1 grid
    orcaGridfn = '/Path/To/Model/Data/mesh_mask_eORCA1_wrk.nc'

    #####
    # Location of data files.
    ObsFolder = "/Path/To/Observation/Data/"
    WOAFolder_annual = ObsFolder + "WOA/annual/"
    WOAFolder = ObsFolder + "WOA/"
    MAREDATFolder = ObsFolder + "/MAREDAT/MAREDAT/"
    GEOTRACESFolder = ObsFolder + "/GEOTRACES/GEOTRACES_PostProccessed/"
    TakahashiFolder = ObsFolder + "/Takahashi2009_pCO2/"
    MLDFolder = ObsFolder + "/IFREMER-MLD/"
    iMarNetFolder = ObsFolder + "/LestersReportData/"
    GlodapDir = ObsFolder + "/GLODAP/"
    GLODAPv2Dir = ObsFolder + "/GLODAPv2/GLODAPv2_Mapped_Climatologies/"
    OSUDir = ObsFolder + "OSU/"
    CCIDir = ObsFolder + "CCI/"
    icFold = ObsFolder + "/InitialConditions/"
