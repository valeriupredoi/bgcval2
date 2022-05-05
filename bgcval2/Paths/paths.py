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

class paths():
    """Empty class to hold paths."""


def paths_setter(paths_dict):
    """Grab all paths and return them."""
    # [general] paths
    paths.machinelocation = paths_dict["general"]["machinelocation"]
    paths.root_dir = paths_dict["general"]["root_dir"]
    paths.shelvedir = paths_dict["general"]["shelvedir"]
    paths.p2p_ppDir = paths_dict["general"]["p2p_ppDir"]
    paths.imagedir = paths_dict["general"]["imagedir"]
    paths.ModelFolder_pref = paths_dict["general"]["ModelFolder_pref"]
    paths.orcaGridfn = paths_dict["general"]["orcaGridfn"]
    paths.ObsFolder = paths_dict["general"]["ObsFolder"]

    # [tata-files] paths
    paths.Dustdir = paths_dict["data-files"]["Dustdir"]
    paths.WOAFolder_annual = paths_dict["data-files"]["WOAFolder_annual"]
    paths.WOAFolder = paths_dict["data-files"]["WOAFolder"]
    paths.DMSDir = paths_dict["data-files"]["DMSDir"]
    paths.MAREDATFolder = paths_dict["data-files"]["MAREDATFolder"]
    paths.GEOTRACESFolder = paths_dict["data-files"]["GEOTRACESFolder"]
    paths.GODASFolder = paths_dict["data-files"]["GODASFolder"]
    paths.TakahashiFolder = paths_dict["data-files"]["TakahashiFolder"]
    paths.MLDFolder = paths_dict["data-files"]["MLDFolder"]
    paths.iMarNetFolder = paths_dict["data-files"]["iMarNetFolder"]
    paths.GlodapDir = paths_dict["data-files"]["GlodapDir"]
    paths.GLODAPv2Dir = paths_dict["data-files"]["GLODAPv2Dir"]
    paths.OSUDir = paths_dict["data-files"]["OSUDir"]
    paths.CCIDir = paths_dict["data-files"]["CCIDir"]
    paths.icFold = paths_dict["data-files"]["icFold"]

    return paths
