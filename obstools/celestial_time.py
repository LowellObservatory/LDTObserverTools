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

"""LDTObserverTools contains python ports of the various DeVeny IDL routines

Lowell Discovery Telescope (Lowell Observatory: Flagstaff, AZ)
http://www.lowell.edu

This file contains various celestial time utilities for LDT.
"""

# Built-In Libraries

# 3rd-Party Libraries
import astropy.coordinates
import astropy.time

# Local Libraries


def lst_midnight(utdates):
    """lst_midnight Compute the LST at midnight for LDT on a list of UT dates

    [extended_summary]

    Parameters
    ----------
    utdates : `list` of `str`
        List of UT dates for which to compute the LST at midnight.  Each
        date must be of the form "YYYY-MM-DD"

    Returns
    -------
    `list` of `str`
        List of the output LST in HH:MM:SS format
    """
    midnights = [f"{date}T07:00:00" for date in utdates]
    times = astropy.time.Time(
        midnights,
        format="isot",
        scale="utc",
        location=astropy.coordinates.EarthLocation.of_site("DCT"),
    )
    return times.sidereal_time("apparent").to_string(precision=0)
