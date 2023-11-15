# -*- coding: utf-8 -*-
#
#  This file is part of LDTObserverTools.
#
#   This Source Code Form is subject to the terms of the Mozilla Public
#   License, v. 2.0. If a copy of the MPL was not distributed with this
#   file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
#  Created on 25-Jan-2021
#
#  @author: tbowers

"""DeVeny Grating Angle Calculator Module

LDTObserverTools contains python ports of various LDT Observer Tools

Lowell Discovery Telescope (Lowell Observatory: Flagstaff, AZ)
https://lowell.edu

This file contains the ``deveny_grangle`` routine for computing the needed grating
tilt angle to be set in order to center the desired wavelength on the CCD.
Both a CLI (direct copy of the IDL version) and GUI version are included here.
"""

# Built-In Libraries
import os
import sys

# 3rd-Party Libraries
import numpy as np
import scipy.optimize
import PySimpleGUI as sg

# Local Libraries
from obstools import utils


# CONSTANTS
PIXSCALE = 2.94  # Base pixels per arcsec: 1 / (0.34 arcsec / pixel)
CAMCOL = np.deg2rad(55.00)  # DeVeny Optical Angle -- Camera-to-Collimator
COLL = np.deg2rad(10.00)  # DeVeny Optical Angle -- Collimator-to-Grating
TGOFFSET = 0.0  # Mechanical Offset in the grating angle


def deveny_grangle_cli():
    """Compute the desired grating angle given grating and central wavelength

    Command line version of ``deveny_grangle``, direct port of the IDL version.
    Takes no arguments, returns nothing, and prints output to screen.
    """
    # Get input from user
    print(" Enter grating resolution (g/mm):")
    gpmm = float(input())
    print(" Enter central wavelength (A):")
    wavelen = float(input())

    # Compute the grating angle and anamorphic demagnification
    grangle, amag = compute_grangle(gpmm, wavelen)

    print(f"\n Grating: {gpmm:.0f} g/mm")
    print(f" Central Wavelength: {wavelen} A")
    print(f" DeVeny grating tilt = {grangle+TGOFFSET:.2f} deg")
    print(
        " Slit demagnification (pixels/arcsec, 0.34 arcsec/pixel): "
        + f"{PIXSCALE*amag:.2f}\n"
    )


def deveny_grangle_gui(max_gui: bool = False):
    r"""Main Driver for the DeVeny Grangle GUI

    Compute the desired grating angle given grating and central wavelength

    GUI version of ``deveny_grangle`` built using ``PySimpleGUI``.  The GUI
    includes a drop-down menu for available gratings, and checks for a requested
    wavelength between 3000A and 11,000A.  This interface uses the same
    subroutines as the CLI version and produces the same results.

    This version optionally allows for the calcuation of the central wavelength
    given a grating angle

    Parameters
    ----------
    max_gui : :obj:`bool`, optional
        Display the MAX GUI (forward and backward calculations)  (Default: False)
    """
    # Define the gratings for the drop-down menu
    gratings = [
        "DV1 - 150 g/mm, 5000 Å",
        "DV2 - 300 g/mm, 4000 Å",
        "DV3 - 300 g/mm, 6750 Å",
        "DV4 - 400 g/mm, 8500 Å",
        "DV5 - 500 g/mm, 5500 Å",
        "DV6 - 600 g/mm, 4900 Å",
        "DV7 - 600 g/mm, 6750 Å",
        "DV8 - 831 g/mm, 8000 Å",
        "DV9 - 1200 g/mm, 5000 Å",
        "DV10 - 2160 g/mm, 5000 Å",
    ]

    # Define the color scheme for the GUI
    sg.theme(utils.SG_THEME)

    # Define the window layout
    row1 = [
        sg.Text("Select Grating:"),
        sg.Drop(
            values=(gratings),
            auto_size_text=True,
            default_value=gratings[1],
            key="Grat",
        ),
    ]
    row2 = [
        sg.Text("Enter Central Wavelength:"),
        sg.Input(key="-WAVEIN-", size=(6, 1)),
        sg.Text("Å"),
    ]
    if max_gui:
        row8 = [
            sg.Text("Enter Grating Tilt:"),
            sg.Input(key="-TILTIN-", size=(6, 1)),
            sg.Text("º"),
        ]
        row3 = [
            sg.Button("Compute Tilt"),
            sg.Button("Compute Wavelength"),
            sg.Button("Done"),
        ]
    else:
        row3 = [sg.Button("Compute"), sg.Button("Done")]
    row4 = [sg.Text("             Grating: "), sg.Text(size=(20, 1), key="-GRATOUT-")]
    row5 = [sg.Text("  Central Wavelength: "), sg.Text(size=(20, 1), key="-WAVEOUT-")]
    row6 = [sg.Text("DeVeny Grating Tilt = "), sg.Text(size=(20, 1), key="-TILTOUT-")]
    row7 = [
        sg.Text('Slit demagnification (0.34"/pixel): '),
        sg.Text(size=(15, 1), key="-DEMAGOUT-"),
    ]

    # Define the rows based on which GUI we're making
    rows = (
        [row1, row2, row8, row3, row4, row5, row6, row7]
        if max_gui
        else [row1, row2, row3, row4, row5, row6, row7]
    )

    # Make the pysimplegui "Error performing wm_overrideredirect" go away
    old_stdout = sys.stdout
    with open(os.devnull, "w", encoding="utf8") as f_null:
        sys.stdout = f_null

        # Create the Window
        window = sg.Window(
            "DeVeny Grating Angle Calculator",
            rows,
            location=(10, 10),
            finalize=True,
            element_justification="center",
            font="Helvetica 18",
        )

    # Return the STDOUT to the command line
    sys.stdout = old_stdout

    # Wait for events
    while True:
        event, values = window.read()
        if event in [sg.WIN_CLOSED, "Done"]:
            break

        if event in ["Compute", "Compute Tilt"]:
            # Check for non-numeric entries for Central Wavelength
            if not values["-WAVEIN-"].isnumeric():
                window["-GRATOUT-"].update("")
                window["-WAVEOUT-"].update("Please Enter a Number")
                window["-TILTOUT-"].update("")
                window["-DEMAGOUT-"].update("")
                if max_gui:
                    window["-TILTIN-"].update("")
                # Wait for next event
                continue

            # Convert wavelen to float, and check for valid range
            wavelen = float(values["-WAVEIN-"])
            if wavelen < 3000 or wavelen > 11000:
                window["-GRATOUT-"].update("")
                window["-WAVEOUT-"].update("Wavelength out of range")
                window["-TILTOUT-"].update("")
                window["-DEMAGOUT-"].update("")
                if max_gui:
                    window["-TILTIN-"].update("")
                # Wait for next event
                continue

            # Compute the grating angle and anamorphic demagnification
            gpmm = float(values["Grat"].split(" - ")[1].split(" g/mm")[0])
            grangle, amag = compute_grangle(gpmm, wavelen)

            # Update the window with the calculated values
            window["-GRATOUT-"].update(values["Grat"])
            window["-WAVEOUT-"].update(f"{values['-WAVEIN-']} Å")
            window["-TILTOUT-"].update(f"{grangle+TGOFFSET:.2f}º")
            window["-DEMAGOUT-"].update(f"{PIXSCALE*amag:.2f} pixels/arcsec")
            if max_gui:
                window["-TILTIN-"].update(f"{grangle+TGOFFSET:.2f}")

        elif event == "Compute Wavelength":
            # Check for non-numeric entries for Grating Tilt
            if not utils.check_float(values["-TILTIN-"]):
                window["-GRATOUT-"].update("")
                window["-WAVEOUT-"].update("")
                window["-TILTOUT-"].update("Please Enter a Number")
                window["-DEMAGOUT-"].update("")
                window["-WAVEIN-"].update("")
                # Wait for next event
                continue

            # Convert wavelen to float, and check for valid range
            tilt = float(values["-TILTIN-"])
            if tilt < 0 or tilt > 48:
                window["-GRATOUT-"].update("")
                window["-WAVEOUT-"].update("")
                window["-TILTOUT-"].update("Tilt angle out of range")
                window["-DEMAGOUT-"].update("")
                window["-WAVEIN-"].update("")
                # Wait for next event
                continue

            # Compute the grating angle and anamorphic demagnification
            gpmm = float(values["Grat"].split(" - ")[1].split(" g/mm")[0])
            wavelen = lambda_at_angle(tilt, gpmm)
            amag = deveny_amag(tilt)

            # Update the window with the calculated values
            window["-GRATOUT-"].update(values["Grat"])
            window["-WAVEOUT-"].update(f"{wavelen:.0f} Å")
            window["-TILTOUT-"].update(f"{tilt+TGOFFSET:.2f}º")
            window["-DEMAGOUT-"].update(f"{PIXSCALE*amag:.2f} pixels/arcsec")
            window["-WAVEIN-"].update(f"{wavelen:.0f}")
        else:
            print("Something funny happened... should never print.")

    # All done, close window
    window.close()


# Computation Utility Functions ==============================================#
def compute_grangle(lpmm: float, wavelen: float):
    """Compute the needed grating angle

    Given the grating's line density and the desired central wavelength,
    compute the required grating angle.  Uses :func:`scipy.optimize.newton` as
    the root-solver.

    Parameters
    ----------
    lpmm : :obj:`float`
        The line density of the grating in g/mm
    wavelen : :obj:`float`
        The central wavelength in angstroms for which to compute the tilt

    Returns
    -------
    grangle : :obj:`float`
        The desired grating angle
    amag : :obj:`float`
        The anamorphic demagnification of the spectrograph at this grangle
    """
    # Initial guess: 20º
    theta = np.deg2rad(20.0)

    # Call the newton method from scipy.optimize to solve the grating equation
    grangle = np.rad2deg(
        scipy.optimize.newton(grangle_eqn, x0=theta, args=(lpmm, wavelen))
    )
    amag = deveny_amag(grangle)
    return grangle, amag


def grangle_eqn(theta: float, lpmm: float, wavelen: float) -> float:
    """The grating equation used to find the angle

    This is the equation for which :func:`scipy.optimize.newton` is finding
    the root.

    Parameters
    ----------
    theta : :obj:`float`
        The grating angle being tested for (in radians)
    lpmm : :obj:`float`
        The line density of the grating in g/mm
    wavelen : :obj:`float`
        The central wavelength in angstroms for which to compute the tilt

    Returns
    -------
    :obj:`float`
        The portion of the grating equation to be set to zero.
    """
    return lambda_at_angle(theta, lpmm, radians=True) - wavelen


def lambda_at_angle(theta: float, lpmm: float, radians: bool = False) -> float:
    """Compute the central wavelength given theta

    Use the grating equation to compute the central wavelength given theta

    Parameters
    ----------
    theta : :obj:`float`
        Grating Angle
    lpmm : :obj:`float`
        The line density of the grating in g/mm
    radians : :obj:`bool`, optional
        The input angle is in radians (Default: False)

    Returns
    -------
    :obj:`float`
        The computed central wavelength
    """
    # Condition the inputs
    if not isinstance(theta, float):
        theta = float(theta)
    if not isinstance(lpmm, float):
        lpmm = float(lpmm)
    if not radians:
        theta = np.deg2rad(theta)

    return (np.sin(COLL + theta) + np.sin(COLL + theta - CAMCOL)) * 1.0e7 / lpmm


def deveny_amag(grangle: float) -> float:
    r"""Compute the anamorphic demagnification of the slit

    The rays hitting the grating in the plane of α and β diffract to the camera
    in such a way that the beam width changes as a function of α and β, whereas
    rays incident on the grating in the perpendicular plane have the same beam
    width upon incidence and reflection. Because there is a difference between
    the beam widths for the two planes, there will be different magnification
    levels (Schweizer, 1979). Whenever perpendicular planes have different
    magnifications, this is called "anamorphic" (de)magnification. Schweizer
    (1979), however, thinks the term "anamorphic magnification" is somewhat
    inaccurate, and prefers "grating magnification". Historically, the DeVeny
    manuals and associated code (`e.g.`, ``deveny_grangle``) use "anamorphic",
    so we continue that there. The resulting magnification in the direction of
    dispersion due to the grating, `r`, arising from differentiation of the
    grating equation, is:

    .. math::

        r \equiv \frac{dβ}{dα} = \frac{cos α}{cosβ}

    (Schweizer, 1979). Practical spectrograph design aligns the slit
    perpendicular to the dispersion direction, and so the change in
    magnification is in the direction of the slit width, hence our quoted
    "anamorphic demagnification of slit width".

    Parameters
    ----------
    grangle : :obj:`float`
        The desired grating angle (in degrees)

    Returns
    -------
    :obj:`float`
        The anamorphic demagnification factor
    """
    alpha = np.deg2rad(grangle) + COLL
    beta = CAMCOL - alpha

    return np.cos(alpha) / np.cos(beta)


# Command Line Script Infrastructure (borrowed from PypeIt) ==================#
class DevenyGrangle(utils.ScriptBase):
    """Script class for ``deveny_grangle`` tool

    Script structure borrowed from :class:`pypeit.scripts.scriptbase.ScriptBase`.
    """

    @classmethod
    def get_parser(cls, width=None):
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
            description="DeVeny Grating Angle Calculator", width=width
        )
        parser.add_argument(
            "--cli",
            action="store_true",
            help="Use the command-line version of this tool",
        )
        parser.add_argument(
            "--max",
            action="store_true",
            help="Use the MAX version of the GUI (compute wavelength from angle)",
        )
        return parser

    @staticmethod
    def main(args):
        """Main Driver

        Simple function that calls the appropriate driver function.
        """
        # Giddy up!
        if args.cli:
            deveny_grangle_cli()
        else:
            deveny_grangle_gui(max_gui=args.max)
