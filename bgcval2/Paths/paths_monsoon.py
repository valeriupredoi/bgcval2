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

#####
# MONSOON
machinelocation = ''
print("Being run at the Met Office on ", gethostname())
machinelocation = 'MONSOON'

ObsFolder = "/projects/ukesm/" + getuser() + "/BGC-data/"
ModelFolder = "/projects/ukesm/" + getuser() + "/UKESM"
#####
# Location of model files.
MEDUSAFolder_pref = ukp.folder(ModelFolder)
NEMOFolder_pref = ukp.folder(ModelFolder)

#####
# Location of data files.
if annual: WOAFolder = ukp.folder(ObsFolder + "WOA/annual")
else: WOAFolder = ukp.folder(ObsFolder + "WOA/")

MAREDATFolder = ObsFolder + "/MAREDAT/MAREDAT/"
GEOTRACESFolder = ObsFolder + "/GEOTRACES/GEOTRACES_PostProccessed/"
TakahashiFolder = ObsFolder + "/Takahashi2009_pCO2/"
MLDFolder = ObsFolder + "/IFREMER-MLD/"
iMarNetFolder = ObsFolder + "/LestersReportData/"
GlodapDir = ObsFolder + "/GLODAP/"
GLODAPv2Dir = ObsFolder + "/GLODAPv2/GLODAPv2_Mapped_Climatologies/"
OSUDir = ObsFolder + "OSU/"
CCIDir = ObsFolder + "CCI/"
if jobID in [
        "xkrus",
]:
    # Old school ORCA1 grid
    orcaGridfn = ModelFolder + '/mesh_mask_ORCA1_75.nc'
else:
    # New eORCA1 grid
    orcaGridfn = ModelFolder + '/mesh_mask_eORCA1_wrk.nc'
