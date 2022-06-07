# -*- coding: utf-8 -*-
#
#  This file is part of PyDeVeny.
#
#   This Source Code Form is subject to the terms of the Mozilla Public
#   License, v. 2.0. If a copy of the MPL was not distributed with this
#   file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
#  Created on 01-Feb-2021
#
#  @author: tbowers

"""Init File
"""

# Boilerplate variables
__author__ = 'Timothy P. Ellsworth Bowers'
__copyright__ = 'Copyright 2022'
__credits__ = ['Lowell Observatory']
__license__ = 'MPL-2.0'
__version__ = '0.1.0'
__email__ = 'tbowers@lowell.edu'
__status__ = 'Development Status :: 4 - Beta'


# Import the user-facing functions to make them available under LDTObserverTools
from .obstools.celestial_time import *
from .obstools.deveny_grangle import main as deveny_grangle
from .obstools.dfocus import dfocus
from .obstools.etc_calc import (
    exptime_given_snr_mag,
    exptime_given_peak_mag,
    snr_given_exptime_mag,
    mag_given_snr_exptime,
    peak_counts
)
