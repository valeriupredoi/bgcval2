# Copyright 2014, Plymouth Marine Laboratory
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
.. module:: summaryTargets
   :platform: Unix
   :synopsis: Produces some straightforward summary metrics, using a list of shelves.
.. moduleauthor:: Lee de Mora <ledm@pml.ac.uk>

"""

from .. import UKESMpython as ukp
from . import makeTargets


def summaryTargets(shelvesAV, imageFold, year):
    """	Produces some straightforward summary metrics, using a list of shelves.
	"""

    #####
    # makeTargets:
    # Make a target diagram of the shelves of this group.
    BGCVALregions = [
        'Global',
        'ignoreInlandSeas',
        'Equator10',
        'ArcticOcean',
        'NorthernSubpolarAtlantic',
        'NorthernSubpolarPacific',
        'SouthernOcean',
        'Remainder',
    ]

    for r in BGCVALregions:
        shelves = ukp.reducesShelves(shelvesAV,
                                     models=[
                                         'NEMO',
                                         'MEDUSA',
                                     ],
                                     sliceslist=[
                                         r,
                                     ])

        filename = imageFold + 'Summary_' + year + '_' + r + '.png'
        print("Summary Target", shelves, '\nfilename:', filename)
        makeTargets(
            shelves,
            filename,
            legendKeys=[
                'name',
            ],
        )
    #####
