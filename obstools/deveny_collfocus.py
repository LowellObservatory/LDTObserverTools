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

"""LDTObserverTools contains python ports of various LDT Observer Tools

Lowell Discovery Telescope (Lowell Observatory: Flagstaff, AZ)
http://www.lowell.edu

This file contains the ``deveny_collfocus`` routine for computing the estimated
collimator focus value and range for the LOUI Focus Sequence tab.  The GUI will
gather current information (mount temperature and grating angle) from the
ActiveMQ broker, if available, but the user can also manually enter these
values when the tool is used away from the telescope.

This calculator was originally written by B. Shafransky (TO), and has been
incorporated into the LDT Observer Tools package.
"""

# Built-In Libraries
import argparse
import os
import sys

# 3rd-Party Libraries
import numpy as np
import PySimpleGUI as sg

try:
    import stomp

    USE_STOMP = True
except ImportError:
    USE_STOMP = False

# Local Libraries
from obstools import utils

# CONSTANTS
MIN_COLLFOC = 7.5  # Smallest collimator focus value (in mm)
STEPSIZE = 0.5  # Default step size (in mm)


def deveny_collfocus():
    """Main Driver for the DeVeny Collimator Focus Sequence Estimator GUI

    Compute the estimated focus and LOUI Focus Sequence range given the mount
    temperature, grating tilt angle, and presence of an order-blocking filter.
    """
    # Define the color scheme for the GUI
    sg.theme(utils.SG_THEME)

    # Define the window layout
    row1 = [
        sg.Text("Enter Mount Temperature:", font="courier 18"),
        sg.Input(key="-TEMPIN-", size=(6, 1)),
        sg.Text("ºC"),
    ]
    row2 = [
        sg.Text("     Enter Grating Tilt:", font="courier 18"),
        sg.Input(key="-TILTIN-", size=(6, 1)),
        sg.Text("º"),
    ]
    row3 = [
        sg.Text("     Select Rear Filter:", font="courier 18"),
        sg.Drop(
            values=(["CLEAR", "GG420", "GG495", "OG570"]),
            auto_size_text=True,
            default_value="CLEAR",
            key="-FILTER-",
        ),
    ]
    row4 = [sg.Button("Compute"), sg.Button("Done")]
    row5 = [
        sg.Text("     Estimated Focus Value:"),
        sg.Text(size=(15, 1), key="-FOCUSOUT-"),
    ]
    row6 = [sg.Text("Suggested Focus Sequence: ", font="Helvetica 18 underline")]
    row7 = [
        sg.Text("  Initial Position (mm):", font="courier 18"),
        sg.Text(
            size=(3, 1),
            key="-STARTPOS-",
            text_color=sg.theme_input_text_color(),
            background_color=sg.theme_input_background_color(),
        ),
        sg.Text("Step Size (mm):", font="courier 18"),
        sg.Text(
            size=(3, 1),
            key="-STEPSIZE-",
            text_color=sg.theme_input_text_color(),
            background_color=sg.theme_input_background_color(),
        ),
    ]
    row8 = [
        sg.Text("        Number of Steps:", font="courier 18"),
        sg.Text(
            size=(3, 1),
            key="-NSTEPS-",
            text_color=sg.theme_input_text_color(),
            background_color=sg.theme_input_background_color(),
        ),
    ]

    # Define the rows based on which GUI we're making
    rows = [row1, row2, row3, row4, row5, row6, row7, row8]

    # Make the pysimplegui "Error performing wm_overrideredirect" go away
    old_stdout = sys.stdout
    with open(os.devnull, "w", encoding="utf8") as f_null:
        sys.stdout = f_null

        # Create the Window
        window = sg.Window(
            "DeVeny Collimator Focus Sequence Estimator",
            rows,
            location=(10, 10),
            finalize=True,
            element_justification="left",
            font="Helvetica 18",
        )

    # Return the STDOUT to the command line
    sys.stdout = old_stdout

    # Wait for events
    while True:
        event, values = window.read()
        if event in [sg.WIN_CLOSED, "Done"]:
            break

        if event == "Compute":
            # Check for non-numeric values in entries
            try:
                tempin = float(values["-TEMPIN-"])
                tiltin = float(values["-TILTIN-"])
            except ValueError:
                window["-FOCUSOUT-"].update("Please Enter a Number")
                window["-STARTPOS-"].update("")
                window["-STEPSIZE-"].update("")
                window["-NSTEPS-"].update("")
                window["-TEMPIN-"].update("")
                window["-TILTIN-"].update("")
                # Wait for next event
                continue

            # Compute the estimated focus and focus sequence range
            est_focus = calculate_collimator_focus(
                tempin, tiltin, values["-FILTER-"] != "CLEAR"
            )
            focseq_range = calculate_focus_sequence(est_focus)

            # Update the window with the calculated values
            window["-FOCUSOUT-"].update(f"{np.round(est_focus,1)} mm")
            window["-STARTPOS-"].update(f"{np.round(focseq_range[0],1)}")
            window["-STEPSIZE-"].update(f"{np.round(focseq_range[1],1)}")
            window["-NSTEPS-"].update(f"{focseq_range[2]}")

        elif event == "Get from Broker":
            # Retrieve the current values from the broker
            pass

    # All done, close window
    window.close()


def calculate_collimator_focus(
    mount_temp: float, grating_tilt: float, order_blocker: bool
) -> float:
    """Calculate Collimator Focus Value

    Using the formula presented in the DeVeny manual (built from collimator
    focus values from 2017-2020), compute the approximate focus value given
    the current conditions.

    The smallest value returned is 7.75mm -- a physical limit due to the limit
    switch on the collimator focus stage.

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

    _extended_summary_

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


# Command Line Entry Point ===================================================#
def entry_point():
    """Command-Line Entry Point"""
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        prog="deveny_collfocus",
        description="DeVeny Collimator Focus Sequence Estimator",
    )
    args = parser.parse_args()

    # Giddy Up!
    sys.exit(deveny_collfocus())
