import os
from getpass import getuser
from socket import gethostname

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
        host = gethostname()
        print("Got host name: ", host)
        raise ValueError(f"Unidentified hostname {host}"
                         f"Run at either JASMIN, MONSOON or PML.")

    return hostname


def _set_jasmin_paths(paths_dict):
    """Fix paths for when running on JASMIN."""
    jasmin_paths = dict(paths_dict)

    # normalize root dir in case user didnt specify abspath
    root_dir = _normalize_path(jasmin_paths["general"]["root_dir"])

    user = getuser()
    shelves_dir = jasmin_paths["general"]["shelvedir"]
    data_dir = jasmin_paths["general"]["ModelFolder_pref"]
    jasmin_paths["general"]["shelvedir"] = os.path.join(
        root_dir,
        data_dir,
        user,
        shelves_dir
    )
    p2p_dir = jasmin_paths["general"]["p2p_ppDir"]
    jasmin_paths["general"]["p2p_ppDir"] = os.path.join(
        root_dir,
        data_dir,
        p2p_dir
    )
    images_dir = jasmin_paths["general"]["imagedir"]
    jasmin_paths["general"]["imagedir"] = os.path.join(
        root_dir,
        data_dir, user,
        images_dir
    )
    jasmin_paths["general"]["ModelFolder_pref"] = os.path.join(
        root_dir,
        data_dir
    )
        
    # normalize obs forlder in case user didnt specify abspath
    obs_folder = _normalize_path(jasmin_paths["general"]["ObsFolder"])

    for obsdir in jasmin_paths["data-files"]:
        jasmin_paths["data-files"][obsdir] = os.path.join(
            obs_folder,
            jasmin_paths["data-files"][obsdir]
        )

    return jasmin_paths


def _set_pml_paths(paths_dict):
    """Fix paths for running on PML."""
    # import UKESM module
    from . import UKESMpython as ukp

    pml_paths = dict(paths_dict)

    # normalize root dir in case user didnt specify abspath
    root_dir = _normalize_path(pml_paths["general"]["root_dir"])

    user = getuser()
    shelves_dir = ukp.folder(pml_paths["general"]["shelvedir"])
    data_dir = pml_paths["general"]["ModelFolder_pref"]
    pml_paths["general"]["shelvedir"] = os.path.join(
        root_dir,
        data_dir,
        user,
        shelves_dir
    )
    p2p_dir = pml_paths["general"]["p2p_ppDir"]
    pml_paths["general"]["p2p_ppDir"] = os.path.join(
        root_dir,
        data_dir,
        p2p_dir
    )
    images_dir = ukp.folder(pml_paths["general"]["imagedir"])
    pml_paths["general"]["imagedir"] = os.path.join(
        root_dir,
        data_dir, user,
        images_dir
    )
    pml_paths["general"]["ModelFolder_pref"] = os.path.join(
        root_dir,
        data_dir
    )

    # normalize obs forlder in case user didnt specify abspath
    obs_folder = _normalize_path(pml_paths["general"]["ObsFolder"])

    for obsdir in pml_paths["data-files"]:
        pml_paths["data-files"][obsdir] = os.path.join(
            obs_folder,
            pml_paths["data-files"][obsdir]
        )

    return pml_paths


def _set_monsoon_paths(paths_dict):
    """Fix runtime paths when running on MONSOON."""

    # import UKESM module
    from . import UKESMpython as ukp

    mons_paths = dict(paths_dict)

    # normalize root dir in case user didnt specify abspath
    root_dir = _normalize_path(mons_paths["general"]["root_dir"])

    user = getuser()
    mons_paths["general"]["ObsFolder"] = os.path.join(
        root_dir,
        user,
        mons_paths["general"]["ModelFolder_pref"]
    )
    mons_paths["general"]["ModelFolder"] = os.path.join(
        root_dir,
        user,
        "UKESM"
    )
    mons_paths["general"]["MEDUSAFolder_pref"] = \
        ukp.folder(mons_paths["general"]["ModelFolder"])
    mons_paths["general"]["NEMOFolder_pref"] = \
        ukp.folder(mons_paths["general"]["ModelFolder"])
    shelves_dir = mons_paths["general"]["shelvedir"]
    data_dir = mons_paths["general"]["ModelFolder_pref"]
    mons_paths["general"]["shelvedir"] = os.path.join(
        root_dir,
        data_dir,
        user,
        shelves_dir
    )
    p2p_dir = mons_paths["general"]["p2p_ppDir"]
    mons_paths["general"]["p2p_ppDir"] = os.path.join(
        root_dir,
        data_dir,
        p2p_dir
    )
    images_dir = mons_paths["general"]["imagedir"]
    mons_paths["general"]["imagedir"] = os.path.join(
        root_dir,
        data_dir, user,
        images_dir
    )
    mons_paths["general"]["ModelFolder_pref"] = os.path.join(
        root_dir,
        data_dir
    )

    # special case for Orca grid file
    # if jobID in ["xkrus"]: mesh_mask_ORCA1_75.nc
    # else mesh_mask_eORCA1_wrk.nc
    if "orcaGridfn_file" in mons_paths["general"]:
        if not os.path.exists(mons_paths["general"]["orcaGridfn_file"]):
            mons_paths["general"]["orcaGridfn"] = os.path.join(
                mons_paths["general"]["ModelFolder"],
                mons_paths["general"]["orcaGridfn_file"]
            )
        else:
            mons_paths["general"]["orcaGridfn"] = \
                mons_paths["general"]["orcaGridfn_file"]
        
    # normalize obs forlder in case user didnt specify abspath
    obs_folder = _normalize_path(mons_paths["general"]["ObsFolder"])

    for obsdir in mons_paths["data-files"]:
        mons_paths["data-files"][obsdir] = os.path.join(
            obs_folder,
            mons_paths["data-files"][obsdir]
        )

    return mons_paths


def _check_paths(paths_dict):
    """Check if each of the paths in paths dict really exists."""
    # don't check for paths to thsese dirs, they will get created later
    
    not_check_dirs = [
        "machinelocation",
        "shelvedir",
        "p2p_ppDir",
        "imagedir",
    ]
    for key, pth in paths_dict["general"].items():
        if key not in not_check_dirs:
            if not os.path.exists(pth):
                raise ValueError(f"Path in [general] {pth} does not exist for "
                                 f"specified path parameter {key}."
                )
        else:
            print(f"Configuration: WARNING [standard-paths][general] "
                  f"Not checking {key}: {pth}")
        print(f"Configuration: [standard-paths][general] {key}: {pth}")
    for key, pth in paths_dict["data-files"].items():
        if not os.path.exists(pth):
            raise ValueError(f"Path in [data-files] {pth} does not exist for "
                             f"specified path parameter {key}."
            )
        print(f"Configuration: [standard-paths][data-files] {key}: {pth}")


def _expand_paths(paths_dict, hostname):
    """Expand paths to correct full abspaths depending run host."""
    if hostname == "jasmin":
        runtime_paths = _set_jasmin_paths(paths_dict)
        _check_paths(runtime_paths)
    elif hostname == "monsoon":
        runtime_paths = _set_monsoon_paths(paths_dict)
        _check_paths(runtime_paths)

    return runtime_paths


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
            if hostname not in user_config["standard-paths"]:
                raise ValueError(f"Running on {hostname} but user config"
                                 f"file does not have {hostname} section.")
            _expand_paths(
                user_config["standard-paths"][hostname],
                hostname
            )
            if "general" in user_config["standard-paths"][hostname]:
                user_generals = user_config["standard-paths"][hostname]["general"]
                for pth_name in user_generals:
                    if pth_name in paths["general"]:
                        paths["general"][pth_name] = user_generals[pth_name]
            if "data-files" in user_config["standard-paths"][hostname]:
                user_obses = user_config["standard-paths"][hostname]["data-files"]
                for pth_name in user_obses:
                    if pth_name in paths["data-files"]:
                        paths["data-files"][pth_name] = user_obses[pth_name]

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
