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


def _establish_hostname():
    """Return the hostname where the run is done."""
    if gethostname().find('ceda.ac.uk') > -1 or gethostname().find(
            'jasmin') > -1 or gethostname().find('jc.rl.ac.uk') > -1:
        hostname = "jasmin"
    elif gethostname().find('monsoon') > -1:
        hostname = "monsoon"
    elif gethostname().find('pmpc') > -1:
        hostname = "pml"
    elif gethostname().find('-az') > -1:
        hostname = "github-actions"  # for testing on GA machine
    else:
        print("Got host name: ", gethostname())
        raise ValueError("Unidentified host.
                          Run at either JASMIN, MONSOON or PML.")

    return hostname


class paths(paths_dict):
    """Grab all paths and return them."""
    hostname = _establish_hostname()

    # [general] paths
    paths.machinelocation = paths_dict["standard-paths"][hostname]["machinelocation"]
    paths.root_dir = paths_dict["standard-paths"][hostname]["root_dir"]
    paths.shelvedir = paths_dict["standard-paths"][hostname]["shelvedir"]
    paths.p2p_ppDir = paths_dict["standard-paths"][hostname]["p2p_ppDir"]
    paths.imagedir = paths_dict["standard-paths"][hostname]["imagedir"]
    paths.ModelFolder_pref = paths_dict["standard-paths"][hostname]["ModelFolder_pref"]
    paths.orcaGridfn = paths_dict["standard-paths"][hostname]["orcaGridfn"]

    # [data-files] paths
    paths.ObsFolder = paths_dict["standard-paths"][hostname]["ObsFolder"]
    paths.Dustdir = paths_dict["data-files"][hostname]["Dustdir"]
    paths.WOAFolder_annual = paths_dict["data-files"][hostname]["WOAFolder_annual"]
    paths.WOAFolder = paths_dict["data-files"][hostname]["WOAFolder"]
    paths.DMSDir = paths_dict["data-files"][hostname]["DMSDir"]
    paths.MAREDATFolder = paths_dict["data-files"][hostname]["MAREDATFolder"]
    paths.GEOTRACESFolder = paths_dict["data-files"][hostname]["GEOTRACESFolder"]
    paths.GODASFolder = paths_dict["data-files"][hostname]["GODASFolder"]
    paths.TakahashiFolder = paths_dict["data-files"][hostname]["TakahashiFolder"]
    paths.MLDFolder = paths_dict["data-files"][hostname]["MLDFolder"]
    paths.iMarNetFolder = paths_dict["data-files"][hostname]["iMarNetFolder"]
    paths.GlodapDir = paths_dict["data-files"][hostname]["GlodapDir"]
    paths.GLODAPv2Dir = paths_dict["data-files"][hostname]["GLODAPv2Dir"]
    paths.OSUDir = paths_dict["data-files"][hostname]["OSUDir"]
    paths.CCIDir = paths_dict["data-files"][hostname]["CCIDir"]
    paths.icFold = paths_dict["data-files"][hostname]["icFold"]
