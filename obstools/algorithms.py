# -*- coding: utf-8 -*-
#
#  This file is part of LDTObserverTools.
#
#   This Source Code Form is subject to the terms of the Mozilla Public
#   License, v. 2.0. If a copy of the MPL was not distributed with this
#   file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
#  Created on 12-Nov-2021
#
#  @author: tbowers

"""Algorithms Module

LDTObserverTools contains python ports of the various DeVeny IDL routines

Lowell Discovery Telescope (Lowell Observatory: Flagstaff, AZ)
https://lowell.edu

This file contains various algorithms needed for observing calculations.
"""

# Built-In Libraries

# 3rd-Party Libraries
import astropy.coordinates
import astropy.time
import numpy as np

# Local Libraries


def lst_midnight(utdates: list[str]) -> np.ndarray:
    """Compute the LST at midnight for LDT on a list of UT dates

    The "LST at Midnight" is a helpful guide for determining what objects may
    be observable on a given night.  This routine returns this value for the
    LDT given the input UT dates.

    Parameters
    ----------
    utdates : :obj:`list`
        List of (:obj:`str`) UT dates for which to compute the LST at midnight.
        Each date must be of the form "YYYY-MM-DD"

    Returns
    -------
    :obj:`~numpy.ndarray`
        Array of the output LST in HH:MM:SS format
    """
    # Check input
    if not isinstance(utdates, list):
        raise ValueError(f"Incorrect input type {type(utdates)}")
    midnights = [f"{date}T07:00:00" for date in utdates]
    times = astropy.time.Time(
        midnights,
        format="isot",
        scale="utc",
        location=astropy.coordinates.EarthLocation.of_site("DCT"),
    )
    return times.sidereal_time("apparent").to_string(precision=0, sep=":", pad=True)
