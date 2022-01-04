# -*- coding: utf-8 -*-
#
#  This file is part of PyDeVeny.
#
#   This Source Code Form is subject to the terms of the Mozilla Public
#   License, v. 2.0. If a copy of the MPL was not distributed with this
#   file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
#  Created on 25-Jan-2021
#
#  @author: tbowers

"""LDTObserverTools contains python ports of various LDT Observer Tools

Lowell Discovery Telescope (Lowell Observatory: Flagstaff, AZ)
http://www.lowell.edu

This file contains the deveny_grangle routine for computing the needed grating
tilt angle to be set in order to center the desired wavelength on the CCD.
Both a CLI (direct copy of the IDL version) and GUI version are included here.
"""

# Built-In Libraries
import os
import sys

# 3rd-Party Libraries
import numpy as np
from scipy import optimize
import PySimpleGUI as sg

# Local Libraries

# CONSTANTS
PIXSCALE = 2.94             # Base pixels per arcsec: 1 / (0.34 arcsec / pixel)
CAMCOL = np.deg2rad(55.00)  # DeVeny Optical Angle -- Camera-to-Collimator
COLL = np.deg2rad(10.00)    # DeVeny Optical Angle -- Collimator-to-Grating
TGOFFSET = 0.0              # Mechanical Offset in the grating angle


def deveny_grangle_cli():
    """Compute the desired grating angle given grating and central wavelength

    Command line version of deveny_grangle, direct copy of the IDL version.
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
    print(" Slit demagnification (pixels/arcsec, 0.34 arcsec/pixel): " + \
          f"{PIXSCALE*amag:.2f}\n")


def deveny_grangle_gui(max_gui=False):
    """Compute the desired grating angle given grating and central wavelength

    GUI version of deveny_grangle.  Uses PySimpleGUI.  Includes a drop-down
    menu for available gratings, and checks for a wavelength between 3000 and
    11,000 angstroms.  Uses the same subroutines as the CLI version and
    produces the same results.

    This version optionally allows for the calcuation of the central wavelength
    given a grating angle
    """
    # Define the gratings for the drop-down menu
    gratings = ["DV1 - 150 g/mm, 5000 Å",
                "DV2 - 300 g/mm, 4000 Å",
                "DV3 - 300 g/mm, 6750 Å",
                "DV4 - 400 g/mm, 8500 Å",
                "DV5 - 500 g/mm, 5500 Å",
                "DV6 - 600 g/mm, 4900 Å",
                "DV7 - 600 g/mm, 6750 Å",
                "DV8 - 831 g/mm, 8000 Å",
                "DV9 - 1200 g/mm, 5000 Å",
                "DV10 - 2160 g/mm, 5000 Å"]

    # Define the color scheme for the GUI
    sg.theme('light grey 1')

    # Define the window layout
    row1 = [sg.Text("Select Grating:"),
            sg.Drop(values=(gratings), auto_size_text=True,
                    default_value=gratings[1],key="Grat")]
    row2 = [sg.Text("Enter Central Wavelength:"),
            sg.Input(key="-WAVEIN-", size=(6,1)), sg.Text("Å")]
    if max_gui:
        row8 = [sg.Text("Enter Grating Tilt:"),
                sg.Input(key="-TILTIN-", size=(6,1)), sg.Text("º")]
        row3 = [sg.Button("Compute Tilt"), sg.Button("Compute Wavelength"),
                sg.Button("Done")]
    else:
        row3 = [sg.Button("Compute"), sg.Button("Done")]
    row4 = [sg.Text("             Grating: "),
            sg.Text(size=(20,1), key="-GRATOUT-")]
    row5 = [sg.Text("  Central Wavelength: "),
            sg.Text(size=(20,1), key="-WAVEOUT-")]
    row6 = [sg.Text("DeVeny Grating Tilt = "),
            sg.Text(size=(20,1), key="-TILTOUT-")]
    row7 = [sg.Text("Slit demagnification (0.34\"/pixel): "),
            sg.Text(size=(15,1), key="-DEMAGOUT-")]

    # Define the rows based on which GUI we're making
    rows = [row1, row2, row8, row3, row4, row5, row6, row7] if max_gui \
        else [row1, row2, row3, row4, row5, row6, row7]

    # Make the pysimplegui "Error performing wm_overrideredirect" go away
    old_stdout = sys.stdout
    with open(os.devnull, 'w', encoding='utf8') as f_null:
        sys.stdout = f_null

    # Create the Window
    window = sg.Window(
        "DeVeny Grating Angle Calculator",
        rows,
        location=(10, 10),
        finalize=True,
        element_justification="center",
        font="Helvetica 18")

    # Return the STDOUT to the command line
    sys.stdout = old_stdout

    # Wait for events
    while True:
        event, values = window.read()
        if event in [sg.WIN_CLOSED, 'Done']:
            break

        if event in ['Compute', 'Compute Tilt']:
            # Check for non-numeric entries for Central Wavelength
            if not values['-WAVEIN-'].isnumeric():
                window['-GRATOUT-'].update("")
                window['-WAVEOUT-'].update("Please Enter a Number")
                window['-TILTOUT-'].update("")
                window['-DEMAGOUT-'].update("")
                if max_gui:
                    window['-TILTIN-'].update("")
                # Wait for next event
                continue

            # Convert wavelen to float, and check for valid range
            wavelen = float(values['-WAVEIN-'])
            if wavelen < 3000 or wavelen > 11000:
                window['-GRATOUT-'].update("")
                window['-WAVEOUT-'].update("Wavelength out of range")
                window['-TILTOUT-'].update("")
                window['-DEMAGOUT-'].update("")
                if max_gui:
                    window['-TILTIN-'].update("")
                # Wait for next event
                continue

            # Compute the grating angle and anamorphic demagnification
            gpmm = float(values['Grat'].split(' - ')[1].split(' g/mm')[0])
            grangle, amag = compute_grangle(gpmm, wavelen)

            # Update the window with the calculated values
            window['-GRATOUT-'].update(values['Grat'])
            window['-WAVEOUT-'].update(f"{values['-WAVEIN-']} Å")
            window['-TILTOUT-'].update(f"{grangle+TGOFFSET:.2f}º")
            window['-DEMAGOUT-'].update(f"{PIXSCALE*amag:.2f} pixels/arcsec")
            if max_gui:
                window['-TILTIN-'].update(f"{grangle+TGOFFSET:.2f}")

        elif event == 'Compute Wavelength':
            # Check for non-numeric entries for Grating Tilt
            if not check_float(values['-TILTIN-']):
                window['-GRATOUT-'].update("")
                window['-WAVEOUT-'].update("")
                window['-TILTOUT-'].update("Please Enter a Number")
                window['-DEMAGOUT-'].update("")
                window['-WAVEIN-'].update("")
                # Wait for next event
                continue

            # Convert wavelen to float, and check for valid range
            tilt = float(values['-TILTIN-'])
            if tilt < 0 or tilt > 48:
                window['-GRATOUT-'].update("")
                window['-WAVEOUT-'].update("")
                window['-TILTOUT-'].update("Tilt angle out of range")
                window['-DEMAGOUT-'].update("")
                window['-WAVEIN-'].update("")
                # Wait for next event
                continue

            # Compute the grating angle and anamorphic demagnification
            gpmm = float(values['Grat'].split(' - ')[1].split(' g/mm')[0])
            wavelen = lambda_at_angle(tilt, gpmm)
            amag = deveny_amag(tilt)

            # Update the window with the calculated values
            window['-GRATOUT-'].update(values['Grat'])
            window['-WAVEOUT-'].update(f"{wavelen:.0f} Å")
            window['-TILTOUT-'].update(f"{tilt+TGOFFSET:.2f}º")
            window['-DEMAGOUT-'].update(f"{PIXSCALE*amag:.2f} pixels/arcsec")
            window['-WAVEIN-'].update(f"{wavelen:.0f}")
        else:
            print("Something funny happened... should never print.")

    # All done, close window
    window.close()


def compute_grangle(lpmm, wavelen):
    """compute_grangle Compute the needed grating angle

    [extended_summary]

    Parameters
    ----------
    lpmm : `float`
        The line density of the grating in g/mm
    wavelen : `float`
        The central wavelength in angstroms for which to compute the tilt

    Returns
    -------
    grangle : `float`
        The desired grating angle
    amag : `float`
        The anamorphic demagnification of the spectrograph at this grangle
    """
    # Initial guess: 20º
    theta = np.deg2rad(20.)

    # Call the newton method from scipy.optimize to solve the grating equation
    grangle = np.rad2deg( optimize.newton(grangle_eqn, theta,
                                          args=(lpmm, wavelen)) )
    amag = deveny_amag(grangle)
    return grangle, amag


def grangle_eqn(theta, lpmm, wavelen):
    """grangle_eqn The grating equation used to find the angle

    The scipy.optimize.newton() function looks for where this equation
    equals zero.

    Parameters
    ----------
    theta : `float`
        The grating angle being tested for
    lpmm : `float`
        The line density of the grating in g/mm
    wavelen : `float`
        The central wavelength in angstroms for which to compute the tilt

    Returns
    -------
    `float`
        The portion of the grating equation to be set to zero.
    """
    return (np.sin((COLL + theta)) + np.sin(COLL + theta - CAMCOL)) * \
         1.e7 / lpmm - wavelen


def lambda_at_angle(theta, lpmm, radians=False):
    """lambda_at_angle Compute the central wavelength given theta

    Use the grating equation to compute the central wavelength given theta

    Parameters
    ----------
    theta : `float`
        Grating Angle
    lpmm : `float`
        The line density of the grating in g/mm
    radians : `bool`, optional
        The input angle is in radians [Default: False]

    Returns
    -------
    `float``
        The computed central wavelength
    """
    # Condition the inputs
    if not isinstance(theta, float):
        theta = float(theta)
    if not isinstance(lpmm, float):
        lpmm = float(lpmm)
    if not radians:
        theta = np.deg2rad(theta)

    return ( np.sin(COLL + theta) + np.sin(COLL + theta - CAMCOL) ) * \
               1.e7 / lpmm


def deveny_amag(grangle):
    """deveny_amag Compute the anamorphic demagnification of the slit

    Computes the anamorphic demagnification of the slit given grangle

    Parameters
    ----------
    grangle : `float`
        The desired grating angle (in degrees)

    Returns
    -------
    `flaot`
        The anamorphic demagnification factor
    """
    alpha = np.deg2rad(grangle) + COLL
    mbeta = CAMCOL - np.deg2rad(np.rad2deg(alpha))

    return np.cos(alpha) / np.cos(mbeta)


def check_float(potential_float):
    """Simple funtion to check whether something is a float"""
    try:
        float(potential_float)
        return True
    except ValueError:
        return False


#=========================================================#
def main(cli=False, max_gui=False):
    """main Main driver for calling the appropriate functions

    Parameters
    ----------
    cli : `bool`, optional
        Run the command-line version of the tool [Default: False]
    max_gui : `bool`, optional
        Run the max-GUI version of the tool [Default: False]
    """
    # If CLI, do this regardless of the MAX option
    if cli:
        deveny_grangle_cli()
    else:
        deveny_grangle_gui(max_gui=max_gui)


if __name__ == "__main__":

    # Set up the environment to import the program
    import argparse

    # Parse command line arguments
    parser = argparse.ArgumentParser(prog='deveny_grangle',
                        description='DeVeny Grating Angle Calculator')
    parser.add_argument('--cli', action='store_true',
                        help='Use the command-line version of this tool')
    parser.add_argument('--max', action='store_true',
                        help='Use the MAX version of the GUI (compute wavelength from angle)')
    args = parser.parse_args()

    # Giddy Up!
    main(cli=args.cli, max_gui=args.max)
