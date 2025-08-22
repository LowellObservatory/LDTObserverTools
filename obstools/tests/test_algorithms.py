# -*- coding: utf-8 -*-
#
#  This file is part of LDTObserverTools.
#
#   This Source Code Form is subject to the terms of the Mozilla Public
#   License, v. 2.0. If a copy of the MPL was not distributed with this
#   file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
#  Created on 06-Nov-2024
#
#  @author: tbowers
# pylint: disable=missing-function-docstring

"""Algorithms TEST Module
"""

import numpy as np
import pytest

from obstools import algorithms


def test_compute_visibility():
    pass

def test_lst_midnight():
    # Test the inputs / outputs
    utdates = ["2020-01-24", "2021-03-02", "2022-04-19", "2023-06-21", "2024-11-29"]
    lsts = algorithms.lst_midnight(utdates)
    assert isinstance(lsts, np.ndarray)
    assert isinstance(lsts[0], str)
    assert len(lsts) == len(utdates)
    assert lsts[0] == "07:46:36"  # LST at midnight on UT 2020-01-24
    assert lsts[2] == "13:23:46"  # LST at midnight on UT 2022-04-19

    # Make sure the routine handles errors properly...
    # Bad Month
    with pytest.raises(ValueError):
        _ = algorithms.lst_midnight(["2022-13-17"])
    # Bad Day
    with pytest.raises(ValueError):
        _ = algorithms.lst_midnight(["2022-03-83"])
    # No "-" separator
    with pytest.raises(ValueError):
        _ = algorithms.lst_midnight(["20210409"])
    # Bad inputs (not list)
    with pytest.raises(ValueError):
        _ = algorithms.lst_midnight("2024-11-07")
