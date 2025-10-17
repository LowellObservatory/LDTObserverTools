# -*- coding: utf-8 -*-
#
#  This file is part of LDTObserverTools.
#
#   This Source Code Form is subject to the terms of the Mozilla Public
#   License, v. 2.0. If a copy of the MPL was not distributed with this
#   file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
#  Created on 23-Nov-2022
#
#  @author: tbowers

"""NEO Confirmation Page Object Ephemeris Generator Module

LDTObserverTools contains python ports of various LDT Observer Tools

Lowell Discovery Telescope (Lowell Observatory: Flagstaff, AZ)
https://lowell.edu

This file contains a tool for querying the JPL Scout database for short-shelf-
life NEOs that have not yet been assigned a Horizons identifier and turning the
returned ephemeris into a file that can be ingested into the LDT TCS for
observations.

The NEOCP (NEO Confirmation Page at the Minor Planet Center):
    https://www.minorplanetcenter.net/iau/NEO/toconfirm_tabular.html

The JPL Scout API Definition:
    https://ssd-api.jpl.nasa.gov/doc/scout.html

The result from a query of JPL Scout is a dictionary with the following
structure::

    count: <str>         # Number of ephemeris time points
    object: <dict>       # Various information about the object
    signature: <dict>    # Signature of JPL Scout, including the version #
    eph: <list>          # The list of ephemeris points
    data-fields: <list>  # List of the data field names for each time point / orbit
    orbit-count: <int>   # The number of Monte Carlo orbits used to compute medians

The ``eph`` member is a list of ``count`` ephemeris time points.  Each point
is a dictionary with the following structure::

    sigma-pos: <str>     # The 1-sigma plane-of-sky uncertainty (arcmin)
    limits: <dict>       # Minimum / Maximum results from the Monte Carlo orbits
    time: <str>          # Time for this ephemeris position
    data: <list>         # List of the ``orbit-count`` individual Monte Carlo vals
    sun-flag: <str>      # Flag of where the sun is (null, a, n, c, *)
    median: <dict>       # Medial results from the Monte Carlo orbits
    sigma-limits: <dict> # The 1-sigma minimum/maximum results

It is likely that the pieces of this that we really need are the median values
for each timestamp to convert into something the LDT TCS can ingest.  The
format of this dictionary is as follows::

    ra: <str>        # J2000 RA (degrees)
    dec: <str>       # J2000 Dec (degrees)
    dra: <str>       # Change in RA (arcsec/min)  (accounts for cos(Dec) factor)
    ddec: <str>      # Change in Dec (arcsec/min)
    rate: <str>      # Plane of sky motion (arcsec/min)
    pa: <str>        # Position angle of motion (computed from dra/ddec)
    vmag: <str>      # Visual magnitude
    elong: <str>     # Solar elongation (degrees)
    moon: <str>      # Lunar separation (degrees)
    el: <str>        # Elevation above the local horizon

Assuming::

    result = requests.get(f"https://ssd-api.jpl.nasa.gov/scout.api", params=query).json()

then the pieces we need to be concerned about are::

    for point in result['eph']:
        time = datetime.datetime.fromisoformat(point['time'])
        coord = astropy.coordinates.SkyCoord(
            point['median']['ra']*u.deg,
            point['median']['dec']*u.deg,
            frame='fk5'
        )

        <write appropriate versions to the output file>

The output format for LDT TCS is::

    yyyy mm dd hh mm ss αh αm αs.sss ±δd δm δs.ss

.. warning::

    This module is not yet functional!

"""

# Built-In Libraries
import argparse
import datetime

# 3rd-Party Libraries
import astropy.coordinates
import astropy.units as u
import requests

# Local Libraries
from obstools import utils


def neocp_ephem(neocp_id):
    """NEO Confirmation Page Ephemeris Generator

    _extended_summary_

    Parameters
    ----------
    neocp_id : :obj:`~typing.Any`
        The NEOCP ID for the object
    """
    now = datetime.datetime.fromisoformat("2022-11-23 22:00:00")
    now_p1d = now + datetime.timedelta(days=1)

    # Build the API query
    query = {
        "tdes": neocp_id,
        "obs-code": "G37",
        "eph-start": now.isoformat("T", timespec="seconds"),
        "eph-stop": now_p1d.isoformat("T", timespec="seconds"),
        "eph-step": "1h",
        # "orbits": "true"
    }

    r = requests.get("https://ssd-api.jpl.nasa.gov/scout.api", params=query, timeout=10)

    print(r.url)

    result = r.json()

    # Playing with the output data
    print(f"Count = {result['count']}")

    for k, v in result.items():
        print(f"{k}: {type(v)}")

    print(f"eph: {len(result['eph'])}")

    print("Type and length of each element in eph:")
    for thing in result["eph"]:
        print(f"{thing['time']}  {len(thing['data'])}")
        print(f"{type(thing)}  {list(thing.keys())}")

        for k, v in thing.items():
            print(f"  {k}: {type(v)}")

        for k, v in thing["median"].items():
            print(f"    {k}: {type(v)}")

    with open(f"{neocp_id}_LDT.eph", "w", encoding="utf-8") as f_obj:
        for point in result["eph"]:
            time = datetime.datetime.fromisoformat(point["time"])
            coord = astropy.coordinates.SkyCoord(
                float(point["median"]["ra"]) * u.deg,
                float(point["median"]["dec"]) * u.deg,
                frame="fk5",
            )

            timestamp = time.strftime("%Y %m %d %H %M %S")
            c_str = (
                coord.to_string("hmsdms")
                .replace("h", " ")
                .replace("m", " ")
                .replace("s", " ")
                .replace("d", " ")
                .replace("  ", " ")
                .replace("  ", " ")
            )
            c = c_str.split(" ")
            c_str = (
                f"{int(c[0]):02d} {int(c[1]):02d} "
                f"{float(c[2]):04.1f} {int(c[3]):+02d} "
                f"{int(c[4]):02d} {float(c[5]):02.0f}"
            )

            f_obj.write(f"{timestamp} {c_str}\n")

        f_obj.write("FK5 J2000.0 2000.0\n")


# Command Line Script Infrastructure (borrowed from PypeIt) ==================#
class NeocpEphem(utils.ScriptBase):
    """Script class for ``neocp_ephem`` tool

    Script structure borrowed from :class:`pypeit.scripts.scriptbase.ScriptBase`.
    """

    @classmethod
    def get_parser(
        cls,
        description: str = None,
        width: int = None,
        formatter: argparse.HelpFormatter = argparse.ArgumentDefaultsHelpFormatter,
    ):
        """Construct the command-line argument parser.

        Parameters
        ----------
        description : :obj:`str`, optional
            A short description of the purpose of the script.
        width : :obj:`int`, optional
            Restrict the width of the formatted help output to be no longer
            than this number of characters, if possible given the help
            formatter.  If None, the width is the same as the terminal
            width.
        formatter : :obj:`~argparse.HelpFormatter`
            Class used to format the help output.

        Returns
        -------
        :obj:`~argparse.ArgumentParser`
            Command-line interpreter.
        """

        parser = super().get_parser(
            description="Generate LDT Ephemeris Files for NEOCP objects", width=width
        )
        parser.add_argument(
            "obj_id",
            action="store",
            type=str,
            help="The NEOCP temporary designation ID (e.g., 'P10vY9r')",
        )
        return parser

    @staticmethod
    def main(args):
        """Main Driver

        Simple function that calls the main driver function.
        """
        # Giddy up!
        neocp_ephem(args.obj_id)
