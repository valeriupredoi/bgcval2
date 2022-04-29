#
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
# momm@pml.ac.uk
#
"""
.. module:: C2Chl
   :platform: Unix
   :synopsis: A class for relationships of carbon to chlorophyll.
.. moduleauthor:: Momme Butenschon <momm@pml.ac.uk>

"""


class C2Chl:
    """class for relationships of carbon to chlorophyll following
    Sathyendranath et al. 2009.
    The class includes the loglog functions of the from log10(C)=m+p*log10(Chl) and their transformation C=m*Chl^p as well as the carbon to chlorophyll ratio C2Chl=m*Chl^(p-1)."""

    def __init__(self, m, p):
        self.loglog = lambda x: self._loglog(x, m, p)
        self.C = lambda x: self._C(x, m, p)
        self.C2Chl = lambda x: self._C2Chl(x, m, p)

    _loglog = lambda self, lchl, m, p: m + p * lchl

    _C = lambda self, chl, m, p: 10**m * chl**p

    _C2Chl = lambda self, chl, m, p: 10**m * chl**(p - 1)

    def __call__(self, lchl):
        return self.loglog(lchl)


#Sathyendranath HPLC (C_P):
SH = C2Chl(1.9, .65)
#Sathyendranath Turner (C_P):
ST = C2Chl(1.81, .63)
#Buck Turner:
B = C2Chl(1.92, .69)

#Sathyendranath HPLC (C_T):
SH_T = C2Chl(2.25, 0.48)
#Sathyendranath Turner (C_T):
ST_T = C2Chl(2.19, 0.45)
