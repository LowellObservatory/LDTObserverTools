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

The output format for LDT TCS is::

    yyyy mm dd hh mm ss αh αm αs.sss ±δd δm δs.ss

.. warning::

    This module is not yet functional!
    
"""

# Built-In Libraries
import argparse
import dataclasses
import datetime
import sys

# 3rd-Party Libraries
import astropy.coordinates
import astropy.table
import astropy.units as u
import numpy as np
from PyQt6 import QtCore, QtGui, QtWidgets
import requests

# Local Libraries
from obstools import utils
from obstools.UI.EphemMainWindow import Ui_EphemMainWindow

LDT_OBSCODE = "G37"


@dataclasses.dataclass
class EphemObj:
    """Ephemeris Object

    _extended_summary_
    """

    epoch: str
    data: astropy.table.Table


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
    astorb_id : str
        _description_
    utstart : datetime.datetime
        _description_
    utend : datetime.datetime
        _description_
    stepsize : datetime.timedelta
        _description_

    Returns
    -------
    EphemObj
        _description_
    """


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


def neocp_ephem(neocp_id):
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


def norad_ephem(
    norad_id: str,
    utstart: datetime.datetime,
    utend: datetime.datetime,
    stepsize: datetime.timedelta,
) -> EphemObj:
    """norad_ephem _summary_

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


# GUI Classes ================================================================#
class EphemWindow(QtWidgets.QMainWindow, Ui_EphemMainWindow):
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
        self.generateSelectedButton.pressed.connect(
            self.generate_selected_button_clicked
        )
        self.generateAllButton.pressed.connect(self.generate_all_button_clicked)

        # self.computeButton.pressed.connect(self.compute_button_clicked)
        # self.radioExpSnMag.toggled.connect(
        #     lambda checked: self.set_stacked_page(0, checked)
        # )
        # self.radioExpPkMag.clicked.connect(
        #     lambda checked: self.set_stacked_page(1, checked)
        # )
        # self.radioSnExpMag.clicked.connect(
        #     lambda checked: self.set_stacked_page(2, checked)
        # )
        # self.radioMagSnExp.clicked.connect(
        #     lambda checked: self.set_stacked_page(3, checked)
        # )
        # self.buttonAdd2Table.clicked.connect(self.add_data_button_clicked)
        # self.buttonShowTable.clicked.connect(self.show_table_button_clicked)

        # Correctly point to the Lowell Logo
        self.LowellLogo.setPixmap(
            QtGui.QPixmap(str(utils.UI / "lowelllogo_horizontal_web.png"))
        )
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

        # Fix the font sizes
        if sys.platform.startswith("linux"):
            # Reset the font size to system
            font = QtGui.QFont()
            self.centralwidget.setFont(font)
            # Make the title label bigger
            font.setPointSize(int(np.round(font.pointSize() * 13 / 7, 0)))
            self.labelTitle.setFont(font)

        # Disable currently unimplemented data sources
        self.sourceHorizons.setDisabled(True)
        self.sourceScout.setDisabled(True)
        self.sourceMPC.setDisabled(True)
        self.sourceAstorb.setDisabled(True)
        self.sourceIMCCE.setDisabled(True)

        # Set default values
        # self.last_input_data = ETCData()
        # self.last_aux_data = AuxData()
        # self.etc_table = astropy.table.Table()

    def exit_button_clicked(self):
        """The user clicked the "Exit" button

        Display a confirmation dialog, and quit if "Yes"
        """
        button = QtWidgets.QMessageBox.question(
            self,
            "",
            "Are you sure you want to quit?",
            buttons=QtWidgets.QMessageBox.StandardButton.Ok
            | QtWidgets.QMessageBox.StandardButton.Cancel,
        )

        if button == QtWidgets.QMessageBox.StandardButton.Ok:
            QtWidgets.QApplication.quit()

    def generate_selected_button_clicked(self):
        """The user clicked the "Generate Selected Ephemeris" button

        Poll the Databse Source radio button and the Ephemeris Options, and the
        currently selected item in the Objects list, then pass everything along
        to the appropriate query function.
        """

    def generate_all_button_clicked(self):
        """The user clicked the "Generate All Ephemerides" button

        Poll the Databse Source radio button and the Ephemeris Options, and the
        entire Objects list, then pass everything along to the appropriate
        query function.
        """


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
