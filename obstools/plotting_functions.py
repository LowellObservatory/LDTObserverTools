# -*- coding: utf-8 -*-
#
#  This file is part of LDTObserverTools.
#
#   This Source Code Form is subject to the terms of the Mozilla Public
#   License, v. 2.0. If a copy of the MPL was not distributed with this
#   file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
#  Created on 22-Aug-2025
#
#  @author: tbowers

"""Plotting Function Module

LDTObserverTools contains python ports of the various DeVeny IDL routines

Lowell Discovery Telescope (Lowell Observatory: Flagstaff, AZ)
https://lowell.edu

This file contains various plotting functions needed by other routines in this
package.
"""

# Built-In Libraries
import pathlib

# 3rd-Party Libraries
import matplotlib.dates as mdates
import matplotlib.pyplot as plt

# Local Libraries
from obstools import utils

# CONSTANTS


# Plotting Routines ==========================================================#
def plot_visibility(
    vis: utils.Visibility, tsz: int | float = 8, out_fn: pathlib.Path = None
) -> pathlib.Path:
    """Create the visibility plot for an object

    _extended_summary_

    Parameters
    ----------
    vis : :class:`~obstools.utils.Visibility`
        Input object for which to make the plot
    tsz : :obj:`int` or :obj:`float`, optional
        Base fontsize to be used for the plot
    out_fn : :obj:`~pathlib.Path`, optional
        Output filename in which to place the plot.  If ``None``, then the plot
        will be placed in the ``visibilities/`` directory tagged with the
        ``obj_id`` and UT Date at the start of the plot.

    Returns
    -------
    :obj:`~pathlib.Path`
        Path to the file where the plot was saved
    """
    # Set the filename of the output
    if not out_fn:
        utdate = vis.ut_time[0]
        out_fn = utils.VISIBILITIES / f"{vis.objname}.{utdate.strftime('%Y%m%d')}.pdf"

    # Set up the plotting environment
    _, axis = plt.subplots()

    # Begin plotting!
    axis.plot(vis.ut_time, vis.elevation, "k-", label="Object Elevation")
    axis.set_xlabel("UT Time", fontsize=tsz)
    axis.set_ylabel("Elevation (deg)", fontsize=tsz)
    set_std_tickparams(axis, tsz)
    axis.tick_params(right=False)
    axis.set_ylim(0, 90)
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter("%b-%d %H:%M"))

    # Create secondary axis
    axis2 = axis.twinx()  # Create a secondary y-axis sharing the x-axis
    axis2.plot(vis.ut_time, vis.solar_elon, "r--", label="Solar Elongation", alpha=0.5)
    axis2.plot(vis.ut_time, vis.lunar_elon, "g:", label="Lunar Elongation", alpha=0.5)
    axis2.set_ylabel("Elongation (deg)", color="C0", fontsize=tsz)
    axis2.tick_params(
        which="both", direction="in", labelsize=tsz, axis="y", labelcolor="C0"
    )
    axis2.set_ylim(0, 180)

    axis.legend(fontsize=tsz, loc="upper left")
    axis2.legend(fontsize=tsz, loc="upper right")
    axis.set_title(f"Visibility for {vis.objname}", fontsize=tsz + 2)

    # Clean up
    plt.tight_layout()
    plt.savefig(out_fn)
    plt.close()

    # Return the output filename
    return out_fn


# Plotting Utility Functions =================================================#
def set_std_tickparams(axis: plt.axis, tsz: int | float):
    """Set standard tick parameters for a plot

    These are my own "standards", based on plots I used to make in IDL.

    Parameters
    ----------
    axis : :obj:`~matplotlib.pyplot.axis`
        PyPlot axis for whom the tick parameters must be set
    tsz : :obj:`int` or :obj:`float`
        TypeSiZe
    """
    axis.tick_params(
        axis="both",
        which="both",
        direction="in",
        top=True,
        right=True,
        labelsize=tsz,
    )
