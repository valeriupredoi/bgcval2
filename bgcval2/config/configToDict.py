#!/usr/bin/ipython -i
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
from configparser import ConfigParser
import os


def checkConfig(Config, debug=False):
    """
	If it's a string, it opens it as a ConfigParser.
	"""
    if type(Config) == type('string'):
        if debug: print("Reading", Config)
        Config1 = ConfigParser()
        Config1.read(Config)
        Config = Config1
    return Config


def configToDict(fn):
    Config = checkConfig(fn)
    sections = Config.sections()

    if len(sections) == 1:
        section = sections[0]
        options = Config.options(section)
        d = {o: Config.get(section, o) for o in options}
        return d

    d = {}
    for section in sections:
        options = Config.options(section)
        d[section] = {o: Config.get(section, o) for o in options}
    return d


if __name__ == "__main__":
    from sys import argv
    try:
        fn = argv[1]
    except:
        print("Please provide a config file")
        exit()

    d = configToDict(fn)

    #
