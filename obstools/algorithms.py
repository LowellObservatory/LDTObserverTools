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
import datetime

# 3rd-Party Libraries
import astropy.coordinates
import astropy.time
import numpy as np

# Local Libraries
from obstools import utils


def compute_visibility(
    obj_id: str,
    ut_times: np.typing.ArrayLike,
    ra: np.typing.ArrayLike,
    dec: np.typing.ArrayLike,
    loc: astropy.coordinates.EarthLocation = None,
) -> utils.Visibility:
    """Compute the visibility of an object

    _extended_summary_

    Parameters
    ----------
    obj_id : :obj:`str`
        The target identification to be placed in the
        :class:`~obstools.utils.Visibility` object
    ut_times : :obj:`~numpy.typing.ArrayLike`
        Array of :obj:`~datetime.datetime` objects for which the visibility is
        to be computed
    ra : :obj:`~numpy.typing.ArrayLike`
        The Right Ascension(s) in DECIMAL DEGREES for which to compute
        visibility.  If an array of R.A.'s is passed in, it must be of equal
        length as the ``ut_times`` array, and is used primarily for
        non-sidereal objects.
    dec : :obj:`~numpy.typing.ArrayLike`
        The Declinations in DECIMAL DEGREES for which to compute visibility.
        If an array of dec's is passed in, it must be of equal length as the
        ``ut_times`` array, and is used primarily for non-sidereal objects.
    loc : :obj:`~astropy.coordinates.EarthLocation`, optional
        The location for which to compute visibilities.  If ``None``, then
        visibilities are computed for LDT.

    Returns
    -------
    :class:`~obstools.utils.Visibility`
        The computed visibility data
    """
    # Error checking and input parsing
    if loc is None:
        loc = astropy.coordinates.EarthLocation.of_site("LDT")
    if not isinstance(loc, astropy.coordinates.EarthLocation):
        raise utils.ObstoolsError(f"Incorrect type of location provided: {type(loc)}")

    # Convert to array things
    try:
        ut_times = np.asarray(ut_times, dtype=datetime.datetime)
        ra = np.asarray(ra, dtype=float)
        dec = np.asarray(dec, dtype=float)
    except ValueError:
        print("Error here, returning")
        return (np.empty(0), np.empty(0), np.empty(0), np.empty(0))

    # Check lengths:
    if ra.size != 1 or dec.size != 1:
        # TODO: But does this ensure what I want?
        if ut_times.size not in (ra.size, dec.size):
            raise utils.ObstoolsError(
                f"Incorrect array lengths for RA/Dec: {ra.size}/{dec.size} "
                f"-- should be {ut_times.size}"
            )
        # TODO: Are there other cases?

    # coord = astropy.coordinates.SkyCoord(ra=ra, dec=dec, frame=)

    az = np.zeros(ut_times.size, dtype=float)
    el = np.zeros(ut_times.size, dtype=float)
    se = np.zeros(ut_times.size, dtype=float)
    le = np.zeros(ut_times.size, dtype=float)

    return utils.Visibility(
        objname=obj_id,
        ut_time=ut_times,
        azimuth=az,
        elevation=el,
        solar_elon=se,
        lunar_elon=le,
    )


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
