# -*- coding: utf-8 -*-
#
#  This file is part of LDTObserverTools.
#
#   This Source Code Form is subject to the terms of the Mozilla Public
#   License, v. 2.0. If a copy of the MPL was not distributed with this
#   file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
#  Created on 08-Aug-2025
#
#  @author: tbowers
# pylint: disable=c-extension-no-member

"""Multiple-Source LDT Ephemeris Generator Module

LDTObserverTools contains python ports of various LDT Observer Tools

Lowell Discovery Telescope (Lowell Observatory: Flagstaff, AZ)
https://lowell.edu

This file contains a GUI tool for querying various online databases to generate
ephemerides for non-sidereal objects that can be ingested into the LDT TCS for
observations.

Possible databse sources include:
    * JPL Horizons (included in TCS, but also here for completeness)
    * JPL Scout / NEOCP
    * NORAD (artifical satellites)
    * Minor Planet Center
    * AstOrb (Lowell)
    * IMCCE (France)

The goal is to use the database-specific API to query ephemeris information and
produce a TCS-compliant file that can be FTP'd to the TCS computer for
ingestion.
"""

# Built-In Libraries
import argparse
import dataclasses
import datetime
import enum
import io
import pathlib
import re

# 3rd-Party Libraries
import astropy.coordinates
import astropy.table
import astropy.time
import astropy.units as u
from bs4 import BeautifulSoup
import numpy as np
import paramiko
from PyQt6 import QtCore, QtGui, QtWidgets
import requests

# Local Libraries
from obstools import plotting_functions
from obstools import utils
from obstools.UI.EphemMainWindow import Ui_EphemMainWindow

LDT_OBSCODE = "G37"

# Define the API as dataclasses, computation routines, and the external script
__all__ = [
    "EphemObj",
    "astorb_ephem",
    "celestrak_ephem",
    "horizons_ephem",
    "imcce_ephem",
    "mpc_ephem",
    "neocp_ephem",
    "projpluto_ephem",
    "EphemerisGenerator",
]


class EphemError(Exception):
    """Base Epehemeris Generator Module Exception"""


class EphemQueryError(EphemError):
    """Ephemeris Query Error"""


class RefFrame(enum.Enum):
    """Reference Frame enumeration class

    This class contains the ``enum`` values used for identifying epoch.
    We do not use the value of the enumerations at the moment, so simply
    auto-set each value.
    """

    FK5 = enum.auto()  # J2000.0 FK5 system
    ICRS = enum.auto()  # International Celestial Reference System
    APPT = enum.auto()  # Local apparent coordinates


@dataclasses.dataclass
class EphemObj:
    """Ephemeris Object data class

    The output format for LDT TCS is::

        yyyy mm dd hh mm ss αh αm αs.sss ±δd δm δs.ss
    """

    source: str
    objname: str
    obj_id: str
    utstart: datetime.datetime
    utend: datetime.datetime
    stepsize: datetime.timedelta
    ref_frame: RefFrame
    data: astropy.table.Table = dataclasses.field(default_factory=astropy.table.Table)

    def __post_init__(self):
        # Structure of the table with predefined columns
        if not self.data:
            # Input Columns from the ephemeris-query methods
            self.data.add_column(
                astropy.table.Column(dtype=datetime.datetime, name="datetime")
            )
            self.data.add_column(
                astropy.table.Column(dtype=astropy.coordinates.SkyCoord, name="coords")
            )
            # Extra Columns for Visibility Plots
            self.data.add_column(astropy.table.Column(dtype=float, name="azimuth"))
            self.data.add_column(astropy.table.Column(dtype=float, name="elevation"))
            self.data.add_column(astropy.table.Column(dtype=float, name="solar_elon"))
            self.data.add_column(astropy.table.Column(dtype=float, name="lunar_elon"))

    def __str__(self) -> str:
        """Print a string representation of the object

        Includes database source, object ID, UT start time, ephemeris step
        size, and total number of points in the table.

        Returns
        -------
        :obj:`str`
            A brief description of the object
        """
        return (
            f"{self.source}  ID: {self.obj_id}  UT: {self.utstart}  "
            f"ΔT: {self.stepsize.total_seconds()/60.:.1f} min  N: {len(self.data)}"
        )

    @property
    def tcs_format(self) -> str:
        """Produce the TCS format, ready to write to disk

        Reminder to self about LDT ephem format:

        yyyy mm dd hh mm ss Ah Am As.ssssss +/-Dd Dm Ds.sssss
        FK5 J2000.0 2000.0

        Data Format
        The 12 data columns are: Year, Month, Day, Hour, Minute, Seconds, RA
        (hours min sec), and Dec(deg min sec), with the following format:

        yyyy mm dd hh mm ss Ah Am As.ssssss +/-Dd Dm Ds.sssss

        Note: All dates / times are UT.

        The final line of the ephemeris file may include the reference frame
        for the ephemeris:
        For J2000 coordinates:
        FK5 J2000.0 2000.0
        For ICRF:
        ICRS ICRF 2000.0
        For apparent coordinates:  "APPT J20xx.xx 20xx.xx" where the 'x' must
        be replaced by the epoch of your observing night.  (e.g. the night of
        31 October 2020 was J2020.83)

        Note: If this line is not included, you need to inform your TO of the
        reference frame.

        Returns
        -------
         :obj:`str`
             The write-ready string ready for sending to the TCS
        """
        # Create a string buffer to capture the table output
        output_buffer = io.StringIO()

        # Write the table to the string buffer in a specified ASCII format
        tcs_ephem = self.parse_datetime_coords()
        tcs_ephem.remove_columns(["azimuth", "elevation", "solar_elon", "lunar_elon"])
        tcs_ephem.write(output_buffer, format="ascii.no_header")

        # Get the string content from the buffer
        ascii_string = output_buffer.getvalue()
        # Append the epoch information -- NO NEWLINE (TCS chokes on it)
        match self.ref_frame:
            case RefFrame.FK5:
                ascii_string += "FK5 J2000.0 2000.0"
            case RefFrame.ICRS:
                ascii_string += "ICRS ICRF 2000.0"
            case RefFrame.APPT:
                t = astropy.time.Time(self.utstart)
                t.format = "jyear_str"
                ascii_string += f"APPT {t.value} {t.value.strip('J')}"

        # Return the string object
        return ascii_string

    @property
    def visibility(self) -> pathlib.Path:
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
        if not self.data:
            return None

        # TODO: Check that the azimuth -> lunar_elon columns have been filled
        vis = utils.Visibility(
            objname=f"{self.source}.{self.objname}",
            ut_time=self.data["datetime"].data,
            azimuth=self.data["azimuth"].data,
            elevation=self.data["elevation"].data,
            solar_elon=self.data["solar_elon"].data,
            lunar_elon=self.data["lunar_elon"].data,
        )
        return plotting_functions.plot_visibility(vis)

    def parse_datetime_coords(self) -> astropy.table.Table:
        """Parse the datetime and SkyCoord into columns

        This class takes complex objects (:obj:`~datetime.datetime` and
        :obj:`~astropy.coordinates.SkyCoord`) as input, but must return
        simple text for the TCS ephemeris.  This method parses out those
        objects into the needed columns so that the individual ephemeris
        source routines don't have to.

        Returns
        -------
        :obj:`~astropy.table.Table`
            The parsed table ready for printing as a TCS ephemeris
        """
        # Make a copy of the data table
        table = self.data.copy()

        # Add table columns
        table["year"] = [dt.year for dt in table["datetime"]]
        table["month"] = [dt.month for dt in table["datetime"]]
        table["day"] = [dt.day for dt in table["datetime"]]
        table["hour"] = [dt.hour for dt in table["datetime"]]
        table["minute"] = [dt.minute for dt in table["datetime"]]
        table["second"] = [dt.second for dt in table["datetime"]]
        table["ra_h"] = [
            int(c.ra.to_string(unit=u.hour, sep=":").split(":")[0])
            for c in table["coords"]
        ]
        table["ra_m"] = [
            int(c.ra.to_string(unit=u.hour, sep=":").split(":")[1])
            for c in table["coords"]
        ]
        table["ra_s"] = [
            float(c.ra.to_string(unit=u.hour, sep=":").split(":")[2])
            for c in table["coords"]
        ]
        table["dec_d"] = [
            int(c.dec.to_string(unit=u.deg, sep=":").split(":")[0])
            for c in table["coords"]
        ]
        table["dec_m"] = [
            int(c.dec.to_string(unit=u.deg, sep=":").split(":")[1])
            for c in table["coords"]
        ]
        table["dec_s"] = [
            float(c.dec.to_string(unit=u.deg, sep=":").split(":")[2])
            for c in table["coords"]
        ]

        # Remove the object columns and return
        table.remove_columns(["datetime", "coords"])
        return table


# Ephemeris Query Functions for Online Databases =============================#
def astorb_ephem(
    astorb_id: str,
    utstart: datetime.datetime,
    utend: datetime.datetime,
    stepsize: datetime.timedelta,
) -> EphemObj:
    """astorb_ephem _summary_

    Lowell Observatory AstOrb Database

    Parameters
    ----------
    astorb_id : :obj:`str`
        The Astorb identification of the object for which to generate ephemeris
    utstart : :obj:`~datetime.datetime`
        The starting UT time for the ephemeris
    utend : :obj:`~datetime.datetime`
        The ending UT time for the ephemeris
    stepsize : :obj:`~datetime.timedelta`
        The stepsize of ephemeris points

    Returns
    -------
    :class:`EphemObj`
        The Ephemeris Object class containing the requested data
    """
    raise NotImplementedError


def celestrak_ephem(
    norad_id: str,
    utstart: datetime.datetime,
    utend: datetime.datetime,
    stepsize: datetime.timedelta,
) -> EphemObj:
    """Get artificial satellite ephemeris from Project Pluto

    .. code-block::

        FORM METHOD=GET ACTION="https://www.projectpluto.com/cgi-bin/sat_id/sat_cgi"

        https://www.projectpluto.com/cgi-bin/sat_id/sat_cgi?
        obj_name=43435&
        time=2025-07-31+00%3A00%3A00&
        round_step=on&
        num_steps=24&
        step_size=1h&
        obs_code=G37&
        show_separate_motions=on

        obj_name is the NORAD id #
        time is UT
        step_size takes units h, m, s etc.
        obs_code is MPC code (G37 for LDT)

    Parameters
    ----------
    norad_id : :obj:`str`
        The NORAD identification of the object for which to generate ephemeris
    utstart : :obj:`~datetime.datetime`
        The starting UT time for the ephemeris
    utend : :obj:`~datetime.datetime`
        The ending UT time for the ephemeris
    stepsize : :obj:`~datetime.timedelta`
        The stepsize of ephemeris points

    Returns
    -------
    :class:`EphemObj`
        The Ephemeris Object class containing the requested data
    """
    # Get the object name and COSPAR-ID from the input ``norad_id``
    common_name, cospar_id = get_name_from_norad(norad_id)

    # Get the TLE via a query to CelesTrak
    query = {"CATNR": norad_id}
    try:
        r = requests.get(
            "https://celestrak.org/NORAD/elements/gp.php",
            params=query,
            timeout=10,
        )
    except requests.exceptions.ReadTimeout as err:
        raise EphemQueryError(
            "Timeout while waiting for https://celestrak.org query"
        ) from err

    # Check if HTTP request was okay
    if not r.ok:
        raise EphemQueryError(
            f"Got code {r.status_code} from https://celestrak.org query"
        )

    # Parse out the return text
    lines = [l.strip() for l in r.text.split("\n")]

    # If no General Perturbations data found for this object, return empty
    if "No GP data found" in lines[0]:
        print(
            f"No TLE found at CelesTrak for {common_name.strip()}({cospar_id.strip()})\n"
        )
        return EphemObj(
            objname=f"{common_name.strip()}({cospar_id.strip()})",
            obj_id=norad_id,
            source="CelesTrak",
            utstart=utstart,
            utend=utend,
            stepsize=stepsize,
            ref_frame=RefFrame.FK5,
        )

    tle_str = "\n".join(lines)
    print(f"This is the TLE:\n{tle_str}")

    query = {
        "format": "text",
        "MAKE_EPHEM": "'YES'",
        "COMMAND": "'TLE'",
        "EPHEM_TYPE": "'OBSERVER'",
        "CENTER": f"'{LDT_OBSCODE}@399'",
        "TLE": "\n".join(lines),
        "START_TIME": utstart.isoformat(timespec="seconds"),
        "STOP_TIME": utend.isoformat(timespec="seconds"),
        "STEP_SIZE": f"'{stepsize.total_seconds()/60.:.0f} MINUTES'",
        "QUANTITIES": "'1,4,9,23,25'",
        "REF_SYSTEM": "'ICRF'",
        "CAL_FORMAT": "'CAL'",
        "CAL_TYPE": "'M'",
        "TIME_DIGITS": "'SECONDS'",
        "ANG_FORMAT": "'HMS'",
        "APPARENT": "'REFRACTED'",
        "RANGE_UNITS": "'KM'",
        "SUPPRESS_RANGE_RATE": "'NO'",
        "SKIP_DAYLT": "'NO'",
        "SOLAR_ELONG": "'0,180'",
        "EXTRA_PREC": "'NO'",
        "R_T_S_ONLY": "'NO'",
        "CSV_FORMAT": "'YES'",
        "OBJ_DATA": "'YES'",
    }

    try:
        r = requests.get(
            "https://ssd.jpl.nasa.gov/api/horizons.api",
            params=query,
            timeout=10,
        )
    except requests.exceptions.ReadTimeout as err:
        raise EphemQueryError(
            "Timeout while waiting for https://ssd.jpl.nasa.gov query"
        ) from err

    # Check if HTTP request was okay
    if not r.ok:
        raise EphemQueryError(
            f"Got code {r.status_code} from https://ssd.jpl.nasa.gov query"
        )

    # Parse out the return text
    lines = r.text.split("\n")

    api_version = lines[0]
    api_source = lines[1]
    about = lines[4]

    # Loop through the lines and get the ephemeris section
    hephem = []
    soe = False
    for line in lines:
        # Find the Start of Ephemeris
        if not soe:
            if "$$SOE" in line:
                soe = True
            continue
        # Stop when we get to the End of Ephemeris
        if "$$EOE" in line:
            break
        # Oherwise, place the line in the hephem list
        hephem.append(line)

    # Break apart each line into columns
    table = []
    for line in hephem:
        cols = line.split(",")

        table.append(
            {
                "UT_TIME": cols[0],
                "RA": cols[3],
                "DEC": cols[4],
                "azimuth": float(cols[5]),
                "elevation": float(cols[6]),
                "MAGNITUDE": cols[7],
                "solar_elon": float(cols[9]),
                "lunar_elon": float(cols[11]),
            }
        )
    table = astropy.table.Table(table)

    # Stuff the information into python data structures for parsing
    table["coords"] = astropy.coordinates.SkyCoord(
        table["RA"], table["DEC"], unit=(u.hour, u.deg), frame="fk5"
    )
    try:
        table["datetime"] = [
            datetime.datetime.strptime(t.strip(), "%Y-%b-%d %H:%M:%S")
            for t in table["UT_TIME"]
        ]
    except ValueError as err:
        # Sometimes JPL just wants to send minutes...
        table["datetime"] = [
            datetime.datetime.strptime(t.strip(), "%Y-%b-%d %H:%M")
            for t in table["UT_TIME"]
        ]

    # JPL Horizons data are J2000
    return EphemObj(
        objname=f"{common_name.strip()}({cospar_id.strip()})",
        obj_id=norad_id,
        source="CelesTrak",
        utstart=utstart,
        utend=utend,
        stepsize=stepsize,
        ref_frame=RefFrame.FK5,
        # TODO: This overwrites my nicely defined table above, but need to
        #       figure out a way to insert `table` columns into `EphemObj.data`
        #       columns.
        data=table[
            "datetime",
            "coords",
            "azimuth",
            "elevation",
            "solar_elon",
            "lunar_elon",
        ],
    )


def horizons_ephem(
    horizons_id: str,
    utstart: datetime.datetime,
    utend: datetime.datetime,
    stepsize: datetime.timedelta,
) -> EphemObj:
    """horizons_ephem _summary_

    Use AstroQuery

    JPL Horizons Queries (astroquery.jplhorizons/astroquery.solarsystem.jpl.horizons)

    Parameters
    ----------
    horizons_id : :obj:`str`
        The JPL/Horizons identification of the object for which to generate
        ephemeris
    utstart : :obj:`~datetime.datetime`
        The starting UT time for the ephemeris
    utend : :obj:`~datetime.datetime`
        The ending UT time for the ephemeris
    stepsize : :obj:`~datetime.timedelta`
        The stepsize of ephemeris points

    Returns
    -------
    :class:`EphemObj`
        The Ephemeris Object class containing the requested data
    """
    common_name, cospar_id = get_name_from_norad(horizons_id)

    # Build the Horizons query
    query = {
        "format": "text",
        "MAKE_EPHEM": "'YES'",
        "COMMAND": f"'{cospar_id}'",
        "EPHEM_TYPE": "'OBSERVER'",
        "CENTER": f"'{LDT_OBSCODE}@399'",
        "START_TIME": utstart.isoformat(timespec="seconds"),
        "STOP_TIME": utend.isoformat(timespec="seconds"),
        "STEP_SIZE": f"'{stepsize.total_seconds()/60.:.0f} MINUTES'",
        "QUANTITIES": "'1,4,9,23,25'",
        "REF_SYSTEM": "'ICRF'",
        "CAL_FORMAT": "'CAL'",
        "CAL_TYPE": "'M'",
        "TIME_DIGITS": "'SECONDS'",
        "ANG_FORMAT": "'HMS'",
        "APPARENT": "'REFRACTED'",
        "RANGE_UNITS": "'KM'",
        "SUPPRESS_RANGE_RATE": "'NO'",
        "SKIP_DAYLT": "'NO'",
        "SOLAR_ELONG": "'0,180'",
        "EXTRA_PREC": "'NO'",
        "R_T_S_ONLY": "'NO'",
        "CSV_FORMAT": "'YES'",
        "OBJ_DATA": "'YES'",
    }
    try:
        r = requests.get(
            "https://ssd.jpl.nasa.gov/api/horizons.api",
            params=query,
            timeout=10,
        )
    except requests.exceptions.ReadTimeout as err:
        raise EphemQueryError(
            "Timeout while waiting for https://ssd.jpl.nasa.gov query"
        ) from err

    # Check if HTTP request was okay
    if not r.ok:
        raise EphemQueryError(
            f"Got code {r.status_code} from https://ssd.jpl.nasa.gov query"
        )

    # Parse out the return text
    lines = r.text.split("\n")

    # api_version = lines[0]
    # api_source = lines[1]
    # about = lines[4]

    # Loop through the lines and get the ephemeris section
    hephem = []
    soe = False
    objname = ""
    for line in lines:
        # Find the Start of Ephemeris
        if not soe:
            if "$$SOE" in line:
                soe = True
            continue
        # Stop when we get to the End of Ephemeris
        if "$$EOE" in line:
            break
        # Oherwise, place the line in the hephem list
        hephem.append(line)

    if not hephem:
        print(f"Nothing returned from JPL/Horizons for {common_name}({cospar_id})\n")
        return EphemObj(
            objname=f"{common_name.strip()}({cospar_id.strip()})",
            obj_id=horizons_id,
            source="JPLHorizons",
            utstart=utstart,
            utend=utend,
            stepsize=stepsize,
            ref_frame=RefFrame.FK5,
        )
    print(f"Building ephemeris for {common_name}({cospar_id})...")

    # Break apart each line into columns
    table = []
    for line in hephem:
        cols = line.split(",")

        table.append(
            {
                "UT_TIME": cols[0],
                "RA": cols[3],
                "DEC": cols[4],
                "azimuth": float(cols[5]),
                "elevation": float(cols[6]),
                "MAGNITUDE": cols[7],
                "solar_elon": float(cols[9]),
                "lunar_elon": float(cols[11]),
            }
        )
    table = astropy.table.Table(table)

    # Stuff the information into python data structures for parsing
    table["coords"] = astropy.coordinates.SkyCoord(
        table["RA"], table["DEC"], unit=(u.hour, u.deg), frame="fk5"
    )
    try:
        table["datetime"] = [
            datetime.datetime.strptime(t.strip(), "%Y-%b-%d %H:%M:%S")
            for t in table["UT_TIME"]
        ]
    except ValueError as err:
        # Sometimes JPL just wants to send minutes...
        table["datetime"] = [
            datetime.datetime.strptime(t.strip(), "%Y-%b-%d %H:%M")
            for t in table["UT_TIME"]
        ]

    # JPL Horizons data are J2000
    return EphemObj(
        objname=f"{common_name.strip()}({cospar_id.strip()})",
        obj_id=horizons_id,
        source="JPLHorizons",
        utstart=utstart,
        utend=utend,
        stepsize=stepsize,
        ref_frame=RefFrame.FK5,
        # TODO: This overwrites my nicely defined table above, but need to
        #       figure out a way to insert `table` columns into `EphemObj.data`
        #       columns.
        data=table[
            "datetime",
            "coords",
            "azimuth",
            "elevation",
            "solar_elon",
            "lunar_elon",
        ],
    )


def imcce_ephem(
    imcce_id: str,
    utstart: datetime.datetime,
    utend: datetime.datetime,
    stepsize: datetime.timedelta,
) -> EphemObj:
    """imcce_ephem _summary_

    Use AstroQuery

    Institut de Mécanique Céleste et de Calcul des Éphémérides (IMCCE) Solar
    System Services (astroquery.imcce/astroquery.solarsystem.imcce)

    Parameters
    ----------
    imcce_id : :obj:`str`
        The IMCCE identification of the object for which to generate ephemeris
    utstart : :obj:`~datetime.datetime`
        The starting UT time for the ephemeris
    utend : :obj:`~datetime.datetime`
        The ending UT time for the ephemeris
    stepsize : :obj:`~datetime.timedelta`
        The stepsize of ephemeris points

    Returns
    -------
    :class:`EphemObj`
        The Ephemeris Object class containing the requested data
    """
    raise NotImplementedError


def mpc_ephem(
    mpc_id: str,
    utstart: datetime.datetime,
    utend: datetime.datetime,
    stepsize: datetime.timedelta,
) -> EphemObj:
    """mpc_ephem _summary_

    Use AstroQuery

    Minor Planet Center Queries (astroquery.mpc/astroquery.solarsystem.MPC)

    Parameters
    ----------
    mpc_id : :obj:`str`
        The MPC identification of the object for which to generate ephemeris
    utstart : :obj:`~datetime.datetime`
        The starting UT time for the ephemeris
    utend : :obj:`~datetime.datetime`
        The ending UT time for the ephemeris
    stepsize : :obj:`~datetime.timedelta`
        The stepsize of ephemeris points

    Returns
    -------
    :class:`EphemObj`
        The Ephemeris Object class containing the requested data
    """
    raise NotImplementedError


def neocp_ephem(
    neocp_id: str,
    utstart: datetime.datetime,
    utend: datetime.datetime,
    stepsize: datetime.timedelta,
) -> EphemObj:
    """NEO Confirmation Page Ephemeris Generator

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


    Parameters
    ----------
    neocp_id : :obj:`str`
        The MPC identification of the object for which to generate ephemeris
    utstart : :obj:`~datetime.datetime`
        The starting UT time for the ephemeris
    utend : :obj:`~datetime.datetime`
        The ending UT time for the ephemeris
    stepsize : :obj:`~datetime.timedelta`
        The stepsize of ephemeris points

    Returns
    -------
    :class:`EphemObj`
        The Ephemeris Object class containing the requested data
    """
    # Note that the API limits the number of ephemeris output records to 500.
    # Specify appropriate values for eph-start, eph-stop, and eph-step to
    # ensure this limit is not exceeded.
    if (n_steps := (utend - utstart) / stepsize) > 500:
        raise EphemQueryError(
            f"Requesting too many steps ({n_steps}); "
            "max = 500 -- adjust start/stop/step!"
        )

    # Build the API query
    query = {
        "tdes": neocp_id,
        "obs-code": LDT_OBSCODE,
        "eph-start": utstart.isoformat(timespec="seconds"),
        "eph-stop": utend.isoformat(timespec="seconds"),
        "eph-step": f"{np.round(stepsize.total_seconds()/60,0):.0f}m",
        "fov-diam": 13.5,  # LMI FOV
        # "orbits": "true"
    }

    try:
        # This service can take a while...
        r = requests.get(
            "https://ssd-api.jpl.nasa.gov/scout.api", params=query, timeout=60
        )
    except requests.exceptions.ReadTimeout as err:
        raise EphemQueryError(
            "Timeout while waiting for https://ssd-api.jpl.nasa.gov query"
        ) from err

    # Check if HTTP request was okay
    if not r.ok:
        raise EphemQueryError(
            f"Got code {r.status_code} from https://ssd-api.jpl.nasa.gov query"
        )

    # The results from this source are in JSON format
    result = r.json()

    # Check that the API version has not changed since development
    #   Version: 1.3 (2021 September)
    if result["signature"]["version"] != "1.3":
        raise EphemQueryError(
            "Server API version different from this client!  "
            f"Expecting 1.3, got {result['signature']['version']}"
        )

    # Parse the returned object into an AstroPy Table
    table = []
    for row in result["eph"]:
        row_dict = {"time": row["time"]}
        row_dict.update(row["median"])
        # Ensure things are floats
        for col in ["elong", "moon", "el"]:
            row_dict[col] = float(row_dict[col])
        table.append(row_dict)

    table = astropy.table.Table(table)
    # Just grab the columns I want, in the order I want them (bad JSON...)
    table = table["time", "ra", "dec", "elong", "moon", "el"]
    table.rename_columns(
        table.colnames,
        ["UT_TIME", "RA", "DEC", "solar_elon", "lunar_elon", "elevation"],
    )
    # JPL Scout doesn't provide azimuith directly.
    table["azimuth"] = np.zeros_like(table["elevation"].data, dtype=float)

    # Stuff the information into python data structures for parsing
    table["coords"] = astropy.coordinates.SkyCoord(
        table["RA"], table["DEC"], unit=u.deg, frame="fk5"
    )
    table["datetime"] = [datetime.datetime.fromisoformat(t) for t in table["UT_TIME"]]

    # JPL Scout data are J2000
    return EphemObj(
        objname=neocp_id,
        obj_id=neocp_id,
        source="JPLScout",
        utstart=utstart,
        utend=utend,
        stepsize=stepsize,
        ref_frame=RefFrame.FK5,
        # TODO: This overwrites my nicely defined table above, but need to
        #       figure out a way to insert `table` columns into `EphemObj.data`
        #       columns.
        data=table[
            "datetime",
            "coords",
            "azimuth",
            "elevation",
            "solar_elon",
            "lunar_elon",
        ],
    )


def projpluto_ephem(
    norad_id: str,
    utstart: datetime.datetime,
    utend: datetime.datetime,
    stepsize: datetime.timedelta,
) -> EphemObj:
    """Get artificial satellite ephemeris from Project Pluto

    .. code-block::

        FORM METHOD=GET ACTION="https://www.projectpluto.com/cgi-bin/sat_id/sat_cgi"

        https://www.projectpluto.com/cgi-bin/sat_id/sat_cgi?
        obj_name=43435&
        time=2025-07-31+00%3A00%3A00&
        round_step=on&
        num_steps=24&
        step_size=1h&
        obs_code=G37&
        show_separate_motions=on

        obj_name is the NORAD id #
        time is UT
        step_size takes units h, m, s etc.
        obs_code is MPC code (G37 for LDT)

    Parameters
    ----------
    norad_id : :obj:`str`
        The NORAD identification of the object for which to generate ephemeris
    utstart : :obj:`~datetime.datetime`
        The starting UT time for the ephemeris
    utend : :obj:`~datetime.datetime`
        The ending UT time for the ephemeris
    stepsize : :obj:`~datetime.timedelta`
        The stepsize of ephemeris points

    Returns
    -------
    :class:`EphemObj`
        The Ephemeris Object class containing the requested data
    """
    # Get the object name and COSPAR-ID from the input ``norad_id``
    common_name, cospar_id = get_name_from_norad(norad_id)

    # Build the ProjectPluto query
    query = {
        "obj_name": norad_id,
        "time": utstart.isoformat(timespec="seconds"),
        "round_step": "on",
        "num_steps": f"{(utend - utstart) / stepsize:.0f}",
        "step_size": f"{stepsize.total_seconds():.0f}s",
        "obs_code": LDT_OBSCODE,
        "show_separate_motions": "on",
    }
    try:
        r = requests.get(
            "https://www.projectpluto.com/cgi-bin/sat_id/sat_cgi",
            params=query,
            timeout=10,
        )
    except requests.exceptions.ReadTimeout as err:
        raise EphemQueryError(
            "Timeout while waiting for https://www.projectpluto.com query"
        ) from err

    # Check if HTTP request was okay
    if not r.ok:
        raise EphemQueryError(
            f"Got code {r.status_code} from https://www.projectpluto.com query"
        )

    print(r.url)

    # Parse out the HTML data, which is inside the <pre></pre> tags
    soup = BeautifulSoup(r.text, "html.parser")
    body = soup.find("pre").get_text()

    # Remove the preamble and ps, leaving just the table; parse with AstroPy
    # Column names in lowercase match directly those in EphemObj.data, others
    # require parsing/processing or are not used
    colnames = [
        "year",
        "month",
        "day",
        "TIME",
        "ra_h",
        "ra_m",
        "ra_s",
        "dec_d",
        "dec_m",
        "dec_s",
        "azimuth",
        "elevation",
        "solar_elon",
        "lunar_elon",
        "D_KM",
        "AS_SEC",
        "PA",
        "DRA",
        "DDEC",
        "MAG",
    ]
    dtypes = [
        int,
        int,
        int,
        str,
        int,
        int,
        float,
        int,
        int,
        float,
        float,
        float,
        float,
        float,
        str,
        str,
        str,
        str,
        str,
        str,
    ]

    # Stuff the thing into an AstroPy table
    try:
        table_array = np.array([r.split() for r in body.split("\n")[8:-5]])
    except ValueError as err:
        raise EphemQueryError(
            "Check that you are requesting valid NORAD / ProjectPluto Objects"
        ) from err

    try:
        _, ncol = table_array.shape
    except ValueError as err:
        print(
            f"Table array from ProjectPluto appears empty for object {common_name}({cospar_id})\n"
        )
        return EphemObj(
            objname=f"{common_name.strip()}({cospar_id.strip()})",
            obj_id=norad_id,
            source="PP",
            utstart=utstart,
            utend=utend,
            stepsize=stepsize,
            ref_frame=RefFrame.FK5,
        )
    table = astropy.table.Table(table_array, names=colnames[:ncol], dtype=dtypes[:ncol])

    # Merge info into object columns
    table["datetime"] = [
        datetime.datetime.fromisoformat(
            f"{row['year']:04d}-{row['month']:02d}-{row['day']:02d}T{row['TIME']}"
        )
        for row in table
    ]
    table["coords"] = [
        astropy.coordinates.SkyCoord(
            ra=f"{row['ra_h']}:{row['ra_m']}:{row['ra_s']}",
            dec=f"{row['dec_d']}:{row['dec_m']}:{row['dec_s']}",
            unit=(u.hour, u.deg),
            frame="fk5",
        )
        for row in table
    ]

    # ProjectPluto data are J2000
    return EphemObj(
        objname=f"{common_name.strip()}({cospar_id.strip()})",
        obj_id=norad_id,
        source="PP",
        utstart=utstart,
        utend=utend,
        stepsize=stepsize,
        ref_frame=RefFrame.FK5,
        # TODO: This overwrites my nicely defined table above, but need to
        #       figure out a way to insert `table` columns into `EphemObj.data`
        #       columns.
        data=table[
            "datetime",
            "coords",
            "azimuth",
            "elevation",
            "solar_elon",
            "lunar_elon",
        ],
    )


# Utility Functions ==========================================================#
def get_name_from_norad(norad_id: str) -> tuple[str, str]:
    """Pull the name and COSPAR ID based on the NORAD ID

    _extended_summary_

    Parameters
    ----------
    norad_id : :obj:`str`
        The input NORAD ID

    Returns
    -------
    objname : :obj:`str`
        The satellite's common name
    cospar_id : :obj:`str`
        The COSPAR ID for the satellite
    """
    table = astropy.table.Table.read(utils.DATA / "satcat.csv", format="ascii.csv")
    print(f"Looking up NORAD ID: {norad_id}")
    idx = table["NORAD_CAT_ID"] == int(norad_id)
    objname = table["OBJECT_NAME"][idx].value[0]
    cospar_id = table["OBJECT_ID"][idx].value[0]
    print(f"Found {objname} = {cospar_id}")
    return (
        remove_parentheses_and_contents(objname).replace(" R/B", "").replace(" ", "_"),
        cospar_id,
    )


def remove_parentheses_and_contents(text: str) -> str:
    """Removes parentheses and their contents from a string.

    Parameters
    ----------
    text : :obj:`str`
        The input string

    Returns
    -------
    :obj:`str`
        The string with parentheses and their contents removed
    """
    # The regex r"\s*\(.*?\)\s*" matches:
    # \s*: zero or more whitespace characters (optional, for cleaning up spaces around parentheses)
    # \( : the literal opening parenthesis (escaped because '(' is a special regex character)
    # .*?: any character (except newline) zero or more times, non-greedily
    # \) : the literal closing parenthesis (escaped)
    # \s*: zero or more whitespace characters (optional)
    return re.sub(r"\s*\(.*?\)\s*", "", text)


# GUI Classes ================================================================#
class EphemWindow(utils.ObstoolsGUI, Ui_EphemMainWindow):
    """Ephemeris Generator Main Window Class

    The UI is defined in EphemMainWindow.ui and translated (via pyuic6) into python
    in EphemMainWindow.py.  This class inherits the UI and defines the various
    actions needed to generate ephemerides from the GUI inputs.
    """

    def __init__(self):
        super().__init__()
        self.setupUi(self)
        self.show()

        # # Connect the table view window
        # self.tableWindow = TableWindow(self.table_colnames)

        # Connect buttons to actions
        self.exitButton.pressed.connect(self.exit_button_clicked)
        self.addObjectsButton.pressed.connect(self.add_objects_button_clicked)
        self.removeObjectsButton.pressed.connect(self.remove_objects_button_clicked)
        self.generateSelectedButton.pressed.connect(
            self.generate_selected_button_clicked
        )
        self.generateAllButton.pressed.connect(self.generate_all_button_clicked)
        self.saveGeneratedButton.pressed.connect(self.save_generated_button_clicked)
        self.visibilitySelectedButton.pressed.connect(
            self.selected_visibility_button_clicked
        )
        self.visibilityAllButton.pressed.connect(self.all_visibility_button_clicked)
        self.inputUTStart.dateTimeChanged.connect(self.check_legal_time)
        self.inputUTEnd.dateTimeChanged.connect(self.check_legal_time)

        # self.radioExpSnMag.toggled.connect(
        #     lambda checked: self.set_stacked_page(0, checked)
        # )
        # self.radioExpPkMag.clicked.connect(
        #     lambda checked: self.set_stacked_page(1, checked)
        # )

        self.set_fonts_and_logo()

        # Set the UT start time to the top of the next hour
        now = datetime.datetime.now(datetime.UTC)
        self.inputUTStart.setDateTime(
            QtCore.QDateTime(
                QtCore.QDate(now.year, now.month, now.day),
                QtCore.QTime(now.hour + 1, 0, 0),
                spec=QtCore.Qt.TimeSpec.UTC,
            )
        )
        # Set the UT end time to 24 hours further on
        now = datetime.datetime.now(datetime.UTC)
        self.inputUTEnd.setDateTime(
            QtCore.QDateTime(
                QtCore.QDate(now.year, now.month, now.day + 1),
                QtCore.QTime(now.hour + 1, 0, 0),
                spec=QtCore.Qt.TimeSpec.UTC,
            )
        )

        # Set the database source buttons as a `QButtonGroup`
        self.sourceGroup = QtWidgets.QButtonGroup(self)
        for i, source in enumerate(
            [
                self.sourceHorizons,
                self.sourceScout,
                self.sourceMPC,
                self.sourceProjPluto,
                self.sourceCelestrak,
                self.sourceAstorb,
                self.sourceIMCCE,
            ],
            1,
        ):
            self.sourceGroup.addButton(source, i)

        # Disable currently unimplemented data sources
        self.sourceHorizons.setDisabled(False)
        self.sourceScout.setDisabled(False)
        self.sourceMPC.setDisabled(True)
        self.sourceProjPluto.setDisabled(False)
        self.sourceCelestrak.setDisabled(False)
        self.sourceAstorb.setDisabled(True)
        self.sourceIMCCE.setDisabled(True)
        # Set the default source
        self.sourceCelestrak.setChecked(True)

        # Set default values
        self.pulled_ephems = []
        self.built_visibilities = []
        # self.last_aux_data = AuxData()
        # self.etc_table = astropy.table.Table()

    # Button Click Callbacks =============================#
    def add_objects_button_clicked(self):
        """The user clicked the "Add Objects" button

        Open a `getMultiLineText` input dialog, then parse the lines input into
        the QListWidget object for display in the main window.
        """
        input_string, ok = QtWidgets.QInputDialog.getMultiLineText(
            self, "Enter Objects", "Enter object IDs, one per line", ""
        )

        # Once the dialog closes, "uncheck" the button
        self.addObjectsButton.setChecked(False)
        if not ok:
            # User cancelled the input
            return

        # Otherwise, fill in the `objectList` object with any entries
        object_list = [
            item.strip() for item in input_string.split("\n") if item.strip()
        ]

        if object_list:
            self.objectList.addItems(object_list)
            n_obj = len(object_list)
            self.labelStatus.setText(
                f"Added {n_obj} object{'s' if n_obj > 1 else ''} to the list"
            )

    def remove_objects_button_clicked(self):
        """The user clicked the "Remove Objects" button

        Check that items are selected, open a "Confirm" dialog, then remove the
        selected lines from the QListWidget.
        """
        if not self.objectList.count() or not self.objectList.selectedIndexes():
            # Set the `labelStatus` to an error message and return
            self.labelStatus.setText("ERROR: No object(s) selected for removal")
            return

        # Ask for confirmation
        axed_lines = sorted(l.row() for l in self.objectList.selectedIndexes())
        button = QtWidgets.QMessageBox.question(
            self,
            "Remove Lines?",
            f"Are you sure you want to remove the {len(axed_lines)} "
            "selected lines from the Objects list?",
        )
        # Once the dialog closes, "uncheck" the button
        self.removeObjectsButton.setChecked(False)
        if button == QtWidgets.QMessageBox.StandardButton.No:
            return

        # Remove the selected items from the ``objectList`` from the bottom
        for line in axed_lines[::-1]:
            item = self.objectList.takeItem(line)
            # Delete the item to free memory
            del item
        self.labelStatus.setText(
            f"Removed {len(axed_lines)} object{'s' if len(axed_lines) > 1 else ''}"
            " from the list"
        )

    def generate_selected_button_clicked(self):
        """The user clicked the "Generate Selected Ephemeris" button

        Poll the Databse Source radio button and the Ephemeris Options, and the
        currently selected item in the Objects list, then pass everything along
        to the appropriate query function.
        """
        # First, check that there are items in the `objectList`, and that at
        #   least one is selected.
        if not self.objectList.count() or not self.objectList.selectedIndexes():
            # Set the `labelStatus` to an error message and return
            self.labelStatus.setText(
                "ERROR: No object(s) selected for ephemeris generation!"
            )
            return

        # Get IDs and Giddy Up!
        selected_ids = [item.text() for item in self.objectList.selectedItems()]
        self.query_ephems(selected_ids)

    def generate_all_button_clicked(self):
        """The user clicked the "Generate All Ephemerides" button

        Poll the Databse Source radio button and the Ephemeris Options, and the
        entire Objects list, then pass everything along to the appropriate
        query function.
        """
        # First, check that there are items in the `objectList`
        if not self.objectList.count():
            # Set the `labelStatus` to an error message and return
            self.labelStatus.setText(
                "ERROR: No object(s) available for ephemeris generation!"
            )
            return

        # Get IDs and Giddy Up!
        selected_ids = [
            self.objectList.item(i).text() for i in range(self.objectList.count())
        ]
        self.query_ephems(selected_ids)

    def selected_visibility_button_clicked(self):
        """The user clicked the "View Selected Visibility" button

        Generate the visibility plot for the selected object(s) and show
        it/them in a separate window with the option to save to disk.

        .. todo::

            Right now, this just saves the selected ones to disk.  Need to
            spawn a seperate window and show the matplotlib object in it...
            that necessitates changes to the underlying plotting routine.
        """
        # First, check that there are items in the `objectList`, and that at
        #   least one is selected.
        if not self.objectList.count() or not self.objectList.selectedIndexes():
            # Set the `labelStatus` to an error message and return
            self.labelStatus.setText(
                "ERROR: No object(s) selected for ephemeris visibility!"
            )
            return

        # Get IDs and Giddy Up!
        selected_ids = [item.text() for item in self.objectList.selectedItems()]
        self.build_visibilities(selected_ids)

    def all_visibility_button_clicked(self):
        """The user clicked the "Save All Visibilities" button

        Generate the visibility plots for all objects and save them to disk.
        """
        # First, check that there are items in the `objectList`
        if not self.objectList.count():
            # Set the `labelStatus` to an error message and return
            self.labelStatus.setText(
                "ERROR: No object(s) available for ephemeris visibility!"
            )
            return

        # Get IDs and Giddy Up!
        selected_ids = [
            self.objectList.item(i).text() for i in range(self.objectList.count())
        ]
        self.build_visibilities(selected_ids)

    def save_generated_button_clicked(self):
        """The user clicked the "Save Generated to TCS" button

        Save the ephemerides in the ``self.pulled_ephems`` attribute to disk /
        SFTP them to TCS in some way.
        """
        # Write the status
        sources = sorted({e.source for e in self.pulled_ephems})
        obj_ids = sorted({e.obj_id for e in self.pulled_ephems})
        self.labelStatus.setText(f"Writing TCS files for {sources} objects {obj_ids}")

        # Loop over the list
        for ephem in self.pulled_ephems:
            fn = (
                utils.EPHEMS
                / f"{ephem.source}.{ephem.objname}.{ephem.utstart.strftime('%Y%m%d')}.eph"
            )
            with open(fn, "w", encoding="utf-8") as f_obj:
                f_obj.write(ephem.tcs_format)

    # GUI Functionality methods ==========================#
    def check_legal_time(self):
        """Check if the start/end times are legal

        Show the end time in RED if it is NOT LATER THAN the start time.
        """
        # Get the current palette
        palette = self.inputUTEnd.palette()
        if (
            self.inputUTEnd.dateTime().toPyDateTime()
            <= self.inputUTStart.dateTime().toPyDateTime()
        ):
            # Set the text color to red
            palette.setColor(QtGui.QPalette.ColorRole.Text, QtGui.QColor("red"))
        else:
            # Set the text color to default
            palette.setColor(QtGui.QPalette.ColorRole.Text, QtGui.QColor(""))
        # Apply the modified palette
        self.inputUTEnd.setPalette(palette)

    def build_visibilities(self, selected_ids: list[str]):
        """Build the visibility plots

        _extended_summary_

        Parameters
        ----------
        selected_ids : :obj:`list`
            List of selected Object IDs for which to pull ephemerides
        """
        # Build up the list of previously pulled ephemerides
        selected_source = self.sourceGroup.checkedButton().text()
        pulled_eph = [
            e.obj_id for e in self.pulled_ephems if e.source == selected_source
        ]

        # Loop through the ID's and make sure ephemerides have been generated
        get_new_ephems = []
        for sel_id in selected_ids:
            # If not already pulled, add it to the list to pull
            if sel_id not in pulled_eph:
                get_new_ephems.append(sel_id)
        if get_new_ephems:
            self.query_ephems(get_new_ephems, clear_pulled=False)

        # Loop through all pulled ephems and construct the visibility objects
        for ephem in self.pulled_ephems:
            self.built_visibilities.append(ephem.visibility)

    def query_ephems(self, selected_ids: list[str], clear_pulled: bool = True):
        """Query the database source

        Pass everything along to the appropriate query function.

        Parameters
        ----------
        selected_ids : :obj:`list`
            List of selected Object IDs for which to pull ephemerides
        clear_pulled : :obj:`bool`, optional
            Clear the list of pulled ephemerides before requesting more?
        """
        # Check that the times are legal
        if (
            self.inputUTEnd.dateTime().toPyDateTime()
            <= self.inputUTStart.dateTime().toPyDateTime()
        ):
            self.labelStatus.setText("ERROR: End time is NOT LATER THAN start time!")
            return

        # Get the selected source and set the status message
        selected_source = self.sourceGroup.checkedButton().text()
        self.labelStatus.setText(
            f"Running {selected_source} ephemeris for object(s) {selected_ids}"
        )

        # Assign the ephem method
        match selected_source:
            case "JPL Horizons":
                ephem_method = horizons_ephem
            case "JPL Scout":
                ephem_method = neocp_ephem
            case "MPC":
                ephem_method = mpc_ephem
            case "NORAD / PP":
                ephem_method = projpluto_ephem
            case "NORAD / CT+H":
                ephem_method = celestrak_ephem
            case "Astorb":
                ephem_method = astorb_ephem
            case "IMCCE":
                ephem_method = imcce_ephem
            case _:
                raise ValueError(
                    f"Unknown database source indicated: {selected_source}"
                )

        # Loop through selected IDs
        if clear_pulled:
            self.pulled_ephems = []
        for sel_id in selected_ids:
            # Pull the ephemeris from the source
            try:
                ephem = ephem_method(
                    sel_id,
                    utstart=self.inputUTStart.dateTime().toPyDateTime(),
                    utend=self.inputUTEnd.dateTime().toPyDateTime(),
                    stepsize=datetime.timedelta(minutes=self.inputStepSize.value()),
                )
            except EphemQueryError as err:
                self.labelStatus.setText(f"ERROR: Query error ({err})")
                return

            # Do something with the ephem object
            self.pulled_ephems.append(ephem)
            # print(ephem.tcs_format)


# Command Line Script Infrastructure (borrowed from PypeIt) ==================#
class EphemerisGenerator(utils.ScriptBase):
    """Script class for ``ephemeris_generator`` tool

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
            description="Generate LDT Ephemeris Files for non-sidereal objects",
            width=width,
        )
        parser.add_argument(
            "--neocp_objid",
            action="store",
            type=str,
            help="The NEOCP temporary designation ID (e.g., 'P10vY9r')",
        )
        return parser

    @staticmethod
    def main(args):
        """Main Driver

        Set up the top-level PyQt6 objects and start the event loop
        """
        # Create the QApplication object and main Qt window
        app = QtWidgets.QApplication([])
        _ = EphemWindow()

        # Giddy up!
        app.exec()
