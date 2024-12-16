# pylint: disable=missing-function-docstring
# -*- coding: utf-8 -*-
#
#  This file is part of LDTObserverTools.
#
#   This Source Code Form is subject to the terms of the Mozilla Public
#   License, v. 2.0. If a copy of the MPL was not distributed with this
#   file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
#  Created on 04-Nov-2024
#
#  @author: tbowers

"""DeVeny Grating Angle Calculator TEST Module
"""

import numpy as np

from obstools import deveny_grangle

# NOTE: No unit tests for deveny_grangle_cli() and deveny_grangle_gui()
#       functions.  Need to find another way to test those.


def test_compute_grangle():
    grangle, amag = deveny_grangle.compute_grangle(400, 8000)
    # Ensure returns are floats and the correct values
    assert isinstance(grangle, float)
    assert isinstance(amag, float)
    assert np.isclose(grangle, 27.8919, atol=1.0e-4)
    assert np.isclose(amag, 0.8257, atol=1.0e-4)


def test_grangle_eqn():
    zero = deveny_grangle.grangle_eqn(np.deg2rad(22.54), 300, 5195)
    # Ensure zero is a float and of the correct value
    assert isinstance(zero, float)
    assert np.isclose(zero, 0, atol=1)
    # Check another grating
    zero = deveny_grangle.grangle_eqn(np.deg2rad(21.00), 150, 7220)
    assert np.isclose(zero, 0, atol=1)


def test_lambda_at_angle():
    cenwave = deveny_grangle.lambda_at_angle(22.54, 300)
    # Ensure cenwave is a float and of the correct value
    assert isinstance(cenwave, float)
    assert np.isclose(cenwave, 5195, atol=1)
    # Check another grating
    cenwave = deveny_grangle.lambda_at_angle(25.60, 500)
    assert np.isclose(cenwave, 5000, atol=1)
    # Test the radians option
    cenwave = deveny_grangle.lambda_at_angle(np.deg2rad(40.96), 831, radians=True)
    assert np.isclose(cenwave, 8499, atol=1)


def test_deveny_amag():
    amag = deveny_grangle.deveny_amag(22.54)
    # Ensure amag is a float of the correct value
    assert isinstance(amag, float)
    assert np.isclose(amag, 0.9122, atol=1.0e-4)
    # Check another value
    amag = deveny_grangle.deveny_amag(27.04)
    assert np.isclose(amag, 0.8391, atol=1.0e-4)
