import os
import sys

from socket import gethostname
from getpass import getuser
from pathlib import Path

import yaml


def _read_yaml(config_file):
    """Read yaml file."""
    with open(config_file, 'r') as file:
        cfg = yaml.safe_load(file)

    return cfg


def get_run_configuration(config_file):
    """
    Get runtime configuration settings

    Read config user and default config file
    and store settings in dictionaries.
    """
    # read defaults first
    basedir = os.path.dirname(__file__)
    default_config_file = os.path.join(basedir,
                                       'default-bgcval2-config.yml')
    defaults = _read_yaml(default_config_file)

    # return defaults if no user config specified
    if config_file == "defaults":
        paths = _get_paths(defaults)
        return paths, defaults
    # use user's config file to replace default values
    else:
        # Read user config file
        config_file = _normalize_path(config_file)
        if not os.path.exists(config_file):
            print(f"Specified config file doesnt exist {config_file}")
        print(f"Will use parameters from {config_file} for runtime.")
        user_cfg = _read_yaml(config_file)

        # and replace all user-specific values
        # treat paths separately though
        defaults = dict(defaults)
        paths = _get_paths(defaults, user_cfg)
        # look for anything that's not [paths, ]
        for elem in user_cfg:
            if elem not in ["standard-paths", ]:
                defaults[elem] = user_cfg[elem]


    return paths, defaults


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
        raise ValueError(f"Unidentified host."
                         f"Run at either JASMIN, MONSOON or PML.")

    return hostname


def _expand_paths(paths_dict, hostname):
    """Expand paths to correct full abspaths depending run host."""
    if hostname == "jasmin":
        jasmin_paths = dict(paths_dict)
        root_dir = _normalize_path(jasmin_paths["general"]["root_dir"])
        user = getuser()
        jasmin_paths["general"]["shelvedir"] = os.path.join(root_dir,
                                                            "BGC_data",
                                                            user, "shelves")
        jasmin_paths["general"]["p2p_ppDir"] = os.path.join(
            root_dir,
            "BGC_data/ukesm_postProcessed"
        )
        jasmin_paths["general"]["imagedir"] = os.path.join(
            root_dir,
            'BGC_data', user, 'images'
        )
        jasmin_paths["general"]["ModelFolder_pref"] = os.path.join(root_dir,
                                                                   "BGC_data")
        obs_folder = jasmin_paths["general"]["ObsFolder"]
        for obsdir in jasmin_paths["data-files"]:
            jasmin_paths["data-files"][obsdir] = os.path.join(obs_folder, obsdir)

        return jasmin_paths


def _get_paths(default_config, user_config=None):
    """Assemble the paths object containing all needed runtime paths."""
    hostname = _establish_hostname()
    if hostname not in default_config["standard-paths"]:
        raise ValueError(f"No section standard-paths forund for {hostname} "
                         f"found in {default_config}.")
    default_paths = default_config["standard-paths"][hostname]

    # dict populated with fully working default paths
    paths = _expand_paths(default_paths, hostname)
    paths = dict(paths)

    # replace with user specifics, if any
    if user_config is not None:
        # return immediately if user has no paths
        if "standard-paths" not in user_config:
            return paths
        # if they have, check and replace what they have
        else:
            if "general" in user_config["standard-paths"]:
                for pth_name in user_config["standard-paths"]["general"]:
                    if pth_name in paths["general"]:
                        norm_path = _normalize_path(user_config[pth_name])
                        paths["general"][pth_name] = norm_path
            if "data-files" in user_config["standard-paths"]:
                for pth_name in user_config["standard-paths"]["data-files"]:
                    if pth_name in paths["data-files"]:
                        norm_path = _normalize_path(user_config[pth_name])
                        paths["data-files"][pth_name] = norm_path

    return paths


def _normalize_path(path):
    """Normalize paths.

    Expand ~ character and environment variables and convert path to absolute.

    Parameters
    ----------
    path: str
        Original path

    Returns
    -------
    str:
        Normalized path
    """
    if path is None:
        return None
    return os.path.abspath(os.path.expanduser(os.path.expandvars(path)))
