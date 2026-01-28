# -*- coding: utf-8 -*-
#
#  This file is part of LDTObserverTools.
#
#   This Source Code Form is subject to the terms of the Mozilla Public
#   License, v. 2.0. If a copy of the MPL was not distributed with this
#   file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
#  Created on 15-Jan-2026
#
#  @author: tbowers
# pylint: disable=c-extension-no-member

"""Plot stellar visibilities

LDTObserverTools contains python ports of various LDT Observer Tools

Lowell Discovery Telescope (Lowell Observatory: Flagstaff, AZ)
https://lowell.edu
"""

# Built-In Libraries
import argparse
import datetime
import pathlib

# 3rd-Party Libraries
import astropy.coordinates
import astropy.table
import astropy.time
import astropy.units as u
import numpy as np

# Local Libraries
from obstools import plotting_functions
from obstools import utils


def make_visibility(
    objname: str, ra: str, dec: str, utstart: str = None, utend: str = None
):
    """Make the visibility

    _extended_summary_

    Parameters
    ----------
    objname : :obj:`str`
        Object name (to be used for the plot)
    ra : :obj:`str`
        Sexagesimal formatted RA with ':'
    dec : str
        Sexagesimal formatted declination with ':'
    """

    if utstart is None:
        utstart = datetime.datetime.now(datetime.UTC).isoformat()
    if utend is None:
        utend = (
            datetime.datetime.fromisoformat(utstart) + datetime.timedelta(days=0.5)
        ).isoformat()
        npts = 100
    else:
        npts = 1000

    # Set up AstroPy objects
    coord = astropy.coordinates.SkyCoord(ra=ra, dec=dec, frame="icrs", unit="hour,deg")
    loc = astropy.coordinates.EarthLocation.of_site("LDT")

    # 1. Initialize start and end as Astropy Time objects
    # (No need for np.datetime64 conversion)
    t_start = astropy.time.Time(datetime.datetime.fromisoformat(utstart))
    t_end = astropy.time.Time(datetime.datetime.fromisoformat(utend))

    # 2. Generate a linear range from 0 to 1
    # This replaces np.linspace(t_start, t_end)
    steps = np.linspace(0, 1, num=npts)

    # 3. Use Astropy arithmetic to calculate the array
    # (t_end - t_start) creates a TimeDelta object
    time_array = t_start + steps * (t_end - t_start)

    # Verify first and last
    print(f"Start: {time_array[0].iso}")
    print(f"End:   {time_array[-1].iso}")

    # Convert to AltAz
    altaz = coord.transform_to(
        astropy.coordinates.AltAz(obstime=time_array, location=loc)
    )

    sun_coord = astropy.coordinates.get_sun(time_array)
    moon_coord = astropy.coordinates.get_body("moon", time_array)

    print(moon_coord[-1])
    print(coord)
    print(coord.transform_to('gcrs'))

    solar_elon = sun_coord.separation(coord)
    lunar_elon = moon_coord.separation(coord)

    data = []
    for t, a, e, s, l in zip(time_array, altaz.az, altaz.alt, solar_elon, lunar_elon):
        data.append(
            {
                "datetime": t.to_datetime(),
                "azimuth": a,
                "elevation": e,
                "solar_elon": s.to("degree").value,
                "lunar_elon": l.to("degree").value,
            }
        )

    data = astropy.table.Table(data)

    plot_fn = plot_visibility(objname, data)


def plot_visibility(objname: str, data: astropy.table.Table) -> pathlib.Path:
    """Creates a visibility plot and returns the pathname thereto

    Visibility plots allow for planning of when to optimally observe
    objects.  This is even more crucial with the moving targets that
    require ephemerides.

    Returns
    -------
    :obj:`~pathlib.Path`
        The path to the created visibility plot
    """
    # If ``self.data`` is empty, do not try to plot visibility
    if not data:
        return None

    # TODO: Check that the azimuth -> lunar_elon columns have been filled
    vis = utils.Visibility(
        objname=f"{objname}",
        ut_time=data["datetime"].data,
        azimuth=data["azimuth"].data,
        elevation=data["elevation"].data,
        solar_elon=data["solar_elon"].data,
        lunar_elon=data["lunar_elon"].data,
    )
    return plotting_functions.plot_visibility(vis)


# Command Line Script Infrastructure (borrowed from PypeIt) ==================#
class SiderealVisibility(utils.ScriptBase):
    """Script class for ``sidereal_visibility`` tool

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
            description="Sidereal Visibility Plotter", width=width
        )
        parser.add_argument("name", type=str, help="Object name")
        parser.add_argument(
            "ra", type=str, help="Object Right Ascension (Sexagesimal with ':')"
        )
        parser.add_argument(
            "dec", type=str, help="Object Declination (Sexagesimal with ':')"
        )
        parser.add_argument(
            "--utstart",
            type=str,
            default=None,
            help="UT Start of Visibility in ISO format",
        )
        parser.add_argument(
            "--utend",
            type=str,
            default=None,
            help="UT End of Visibility in ISO format",
        )
        return parser

    @staticmethod
    def main(args):
        """Main Driver

        Simple function that calls the primary function.
        """
        # Giddy Up!
        make_visibility(
            args.name, args.ra, args.dec, utstart=args.utstart, utend=args.utend
        )
