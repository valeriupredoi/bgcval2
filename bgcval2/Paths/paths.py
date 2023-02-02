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
import os

class paths():
    """Empty class to hold paths."""


def paths_setter(paths_dict):
    """Grab all paths and return them."""

    paths_dir = os.path.dirname(os.path.realpath(__file__))
    paths.bgcval2_repo = os.path.dirname(os.path.dirname(paths_dir))

    # [general] paths; mandatory everywhere
    paths.machinelocation = paths_dict["general"]["machinelocation"]
    paths.root_dir = paths_dict["general"]["root_dir"]
    paths.shelvedir = paths_dict["general"]["shelvedir"]
    paths.p2p_ppDir = paths_dict["general"]["p2p_ppDir"]
    paths.shared_mass_scripts = paths_dict["general"]["shared_mass_scripts"]
    paths.imagedir = paths_dict["general"]["imagedir"]
    paths.ModelFolder_pref = paths_dict["general"]["ModelFolder_pref"]
    paths.orcaGridfn = paths_dict["general"].get("orcaGridfn", None)
    paths.amm7Gridfn = paths_dict["general"].get("amm7Gridfn", None)
    paths.orca1bathy = paths_dict["general"].get("orca1bathy", None)
    paths.ObsFolder = paths_dict["general"]["ObsFolder"]

    # [data-files] paths; optional depending on site
    if "Dustdir" in paths_dict["data-files"]:
        paths.Dustdir = paths_dict["data-files"]["Dustdir"]
    if "WOAFolder_annual" in paths_dict["data-files"]:
        paths.WOAFolder_annual = paths_dict["data-files"]["WOAFolder_annual"]
    if "WOAFolder" in paths_dict["data-files"]:
        paths.WOAFolder = paths_dict["data-files"]["WOAFolder"]
    if "DMSDir" in paths_dict["data-files"]:
        paths.DMSDir = paths_dict["data-files"]["DMSDir"]
    if "MAREDATFolder" in paths_dict["data-files"]:
        paths.MAREDATFolder = paths_dict["data-files"]["MAREDATFolder"]
    if "GEOTRACESFolder" in paths_dict["data-files"]:
        paths.GEOTRACESFolder = paths_dict["data-files"]["GEOTRACESFolder"]
    if "GODASFolder" in paths_dict["data-files"]:
        paths.GODASFolder = paths_dict["data-files"]["GODASFolder"]
    if "TakahashiFolder" in paths_dict["data-files"]:
        paths.TakahashiFolder = paths_dict["data-files"]["TakahashiFolder"]
    if "MLDFolder" in paths_dict["data-files"]:
        paths.MLDFolder = paths_dict["data-files"]["MLDFolder"]
    if "iMarNetFolder" in paths_dict["data-files"]:
        paths.iMarNetFolder = paths_dict["data-files"]["iMarNetFolder"]
    if "GlodapDir" in paths_dict["data-files"]:
        paths.GlodapDir = paths_dict["data-files"]["GlodapDir"]
    if "GLODAPv2Dir" in paths_dict["data-files"]:
        paths.GLODAPv2Dir = paths_dict["data-files"]["GLODAPv2Dir"]
    if "OSUDir" in paths_dict["data-files"]:
        paths.OSUDir = paths_dict["data-files"]["OSUDir"]
    if "CCIDir" in paths_dict["data-files"]:
        paths.CCIDir = paths_dict["data-files"]["CCIDir"]
    if "icFold" in paths_dict["data-files"]:
        paths.icFold = paths_dict["data-files"]["icFold"]

    # OPTIONALS: site-specific
    # special paths for MONSOON
    if "MEDUSAFolder_pref" in paths_dict["general"]:
        paths.MEDUSAFolder_pref = paths_dict["general"]["MEDUSAFolder_pref"]
    if "NEMOFolder_pref" in paths_dict["general"]:
        paths.NEMOFolder_pref = paths_dict["general"]["NEMOFolder_pref"]
    

    return paths
