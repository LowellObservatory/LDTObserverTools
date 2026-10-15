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
import argparse

# 3rd-Party Libraries
import numpy as np
import scipy.optimize

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
            description="DeVeny Grating Angle Calculator", width=width
        )
        parser.add_argument(
            "--cli",
            action="store_true",
            help="Use the command-line version of this tool",
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
            from obstools.deveny_gui import run_gui

            run_gui(tab=0)
