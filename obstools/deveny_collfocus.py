# -*- coding: utf-8 -*-
#
#  This file is part of LDTObserverTools.
#
#   This Source Code Form is subject to the terms of the Mozilla Public
#   License, v. 2.0. If a copy of the MPL was not distributed with this
#   file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
#  Created on 31-May-2023
#
#  @author: bshafransky, tbowers

"""DeVeny Collimator Focus Range Estimator GUI Module

LDTObserverTools contains python ports of various LDT Observer Tools

Lowell Discovery Telescope (Lowell Observatory: Flagstaff, AZ)
https://lowell.edu

This file contains the ``deveny_collfocus`` routine for computing the estimated
collimator focus value and range for the LOUI Focus Sequence tab.  The GUI can
gather current information (mount temperature and grating angle) from the
ActiveMQ broker, if available, but the user may also manually enter these
values when the tool is used away from the telescope.

This calculator was originally written by B. Shafransky (LDT TO), and has been
incorporated into the LDTObserverTools package.
"""

# Built-In Libraries
import argparse

# 3rd-Party Libraries
import numpy as np

# Local Libraries
from obstools import utils

# CONSTANTS
MIN_COLLFOC = 7.5  # Smallest collimator focus value (in mm)
STEPSIZE = 0.5  # Default step size (in mm)


# Computation Utility Functions ==============================================#
def calculate_collimator_focus(
    mount_temp: float, grating_tilt: float, order_blocker: bool
) -> float:
    """Calculate Collimator Focus Value

    Using the formula presented in the DeVeny manual (built from collimator
    focus values from 2017-2020), compute the approximate focus value given
    the current conditions.

    The smallest value returned is 7.75mm -- a physical limit due to the limit
    switch on the collimator focus stage.

    Parameters
    ----------
    mount_temp : :obj:`float`
        The current mount temperature in degrees Celsius
    grating_tilt : :obj:`float`
        The current grating tilt angle in degrees
    order_blocker : :obj:`bool`
        Whether an order-blocking filter (`e.g.`, OG570) is in place

    Returns
    -------
    :obj:`float`
        The estimated collimator focus value
    """

    return max(
        11.0
        - 0.08 * mount_temp
        - 0.14 * (grating_tilt - 25.0)
        + 0.7 * int(order_blocker),
        7.75,
    )


def calculate_focus_sequence(estimated_focus: float) -> tuple[float, float, int]:
    """Calculate Estimated Focus Sequence

    The DeVeny LOUI Focus Sequence tab requests three pieces of information in
    order to execute a focus run: starting position, step size, and number of
    steps.  This function computes these three values based on the estimated
    focus value (from temperature and tilt) and some prior information about
    the behavior of the collimator.

    Parameters
    ----------
    estimated_focus : :obj:`float`
        The estimated focus value from the equation

    Returns
    -------
    startpos : :obj:`float`
        The initial position of the Focus Sequence
    stepsize : :obj:`float`
        The stepsize to be used
    n_steps : :obj:`int`
        The number of steps for the Focus Sequence
    """
    # Compute the range from `estimated_focus` to the minumum (no negatives!)
    dist_above_min = max(estimated_focus - MIN_COLLFOC, 0)
    base_steps = (dist_above_min / STEPSIZE) * 2
    # Make n_steps ODD, with a minumum of 5 steps, and a maximum of 9
    n_steps = min(max(np.ceil(base_steps) // 2 * 2 + 1, 5), 9)
    # Determine initial position (est focus minus 1/2 the range, with limit)
    startpos = max(
        np.round(estimated_focus * 2, 0) / 2 - (n_steps - 1) * STEPSIZE / 2, 7.5
    )
    # Return
    return startpos, STEPSIZE, int(n_steps)


def extract_broker_values(status_dict: dict) -> tuple[str, str]:
    """Get Current Values from the Broker

    When running at LDT, we can use the ActiveMQ broker to get the current
    mount temperature and grating tilt.  This is a fancy little trick that
    will make the software developer very happy.

    Parameters
    ----------
    status_dict : :obj:`dict`
        The status dictionary from the broker listener

    Returns
    -------
    :obj:`str`
        The string value (or ``"~ Not Found ~"``) corresponding to this
        dictionary
    """
    # Grating Tilt Angle
    if "GratingTilt" in status_dict:
        return f"{status_dict['GratingTilt']:.2f}"

    # Mount Temperature
    if "MountTemperature" in status_dict:
        return f"{status_dict['MountTemperature']:.1f}"

    # Else, nothing found
    return "~ Not Found ~"


# Command Line Script Infrastructure (borrowed from PypeIt) ==================#
class DevenyCollfocus(utils.ScriptBase):
    """Script class for ``deveny_collfocus`` tool

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
            description="DeVeny Collimator Focus Sequence Estimator", width=width
        )
        parser.add_argument(
            "-d", "--debug", action="store_true", default=False, help=argparse.SUPPRESS
        )
        return parser

    @staticmethod
    def main(args):
        """Main Driver

        Simple function that calls the main driver function.
        """
        # Giddy up!
        from obstools.deveny_gui import run_gui

        run_gui(tab=1)
