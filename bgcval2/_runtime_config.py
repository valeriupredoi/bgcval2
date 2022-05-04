import os
import sys
from pathlib import Path

import yaml


def _read_yaml(config_file):
    """Read yaml file."""
    with open(config_file, 'r') as file:
        cfg = yaml.safe_load(file)

    return cfg


def read_config_user_file(config_file):
    """Read config user file and store settings in a dictionary."""
    # read defaults first
    basedir = os.path.dirname(__file__)
    default_config_file = os.path.join(basedir,
                                       'default-bgcval2-config.yml')
    defaults = _read_yaml(default_config_file)

    # return defaults if no user config specified
    if config_file == "defaults":
        return defaults
    # use user's config file to replace default values
    else:
        # Read user config file
        config_file = _normalize_path(config_file)
        if not os.path.exists(config_file):
            print(f"Specified config file doesnt exist {config_file}")
        print(f"Will use parameters from {config_file} for runtime.")
        user_cfg = _read_yaml(config_file)

        # and replace all user-specific values
        for elem in user_cfg:
            defaults[elem] = user_cfg[elem]


    return defaults


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
