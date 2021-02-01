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

"""PyDeVeny contains python ports of the various DeVeny IDL routines

Lowell Discovery Telescope (Lowell Observatory: Flagstaff, AZ)
http://www.lowell.edu

This file contains the deveny_grangle routine for computing the needed grating
tilt angle to be set in order to center the desired wavelength on the CCD.
Both a CLI (direct copy of the IDL version) and GUI version are included here.
"""

# Built-In Libraries
import sys

# Numpy & SciPy & PySimpleGui
import numpy as np
from scipy import optimize
import PySimpleGUI as sg

# CONSTANTS
PIXSCALE = 2.94     # Base pixels per arcsec: 1 / (0.34 arcsec / pixel)


def deveny_grangle_cli():
    """Compute the desired grating angle given grating and central wavelength
    
    Command line version of deveny_grangle, direct copy of the IDL version.
    Takes no arguments, returns nothing, and prints output to screen.
    """

    # Global variables -- available to grating_eqn function
    global gpmm, wavelen

    # Mechanical Offset
    tgoffset = 0.0

    # Get input from user
    print(" Enter grating resolution (g/mm):")    
    gpmm = float(input())
    print(" Enter central wavelength (A):")
    wavelen = float(input())

    # Compute the grating angle and anamorphic demagnification
    grangle, amag = compute_grangle(wavelen, gpmm)

    print(f"\n Grating: {gpmm:.0f} g/mm")
    print(f" Central Wavelength: {wavelen} A")
    print(f" DeVeny grating tilt = {grangle+tgoffset:.2f} deg")
    print(f" Slit demagnification (pixels/arcsec, 0.34 arcsec/pixel): " + \
          f"{PIXSCALE*amag:.2f}\n")


def deveny_grangle_gui():
    """Compute the desired grating angle given grating and central wavelength
    
    GUI version of deveny_grangle.  Uses PySimpleGUI.  Includes a drop-down
    menu for available gratings, and checks for a wavelength between 3000 and
    11,000 angstroms.  Uses the same subroutines as the CLI version and 
    produces the same results.
    """

    # Global variables -- available to grating_eqn function
    global gpmm, wavelen

    # Mechanical Offset
    tgoffset = 0.0

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
    row3 = [sg.Button("Compute"), sg.Button("Done")]
    row4 = [sg.Text("             Grating: "), 
            sg.Text(size=(20,1), key="-GRATOUT-")]
    row5 = [sg.Text("  Central Wavelength: "), 
            sg.Text(size=(20,1), key="-WAVEOUT-")]
    row6 = [sg.Text("DeVeny Grating Tilt = "), 
            sg.Text(size=(20,1), key="-TILTOUT-")]
    row7 = [sg.Text("Slit demagnification (0.34\"/pixel): "),
            sg.Text(size=(15,1), key="-DEMAGOUT-")]

    # Create the Window
    window = sg.Window(
        "DeVeny Grating Angle Calculator",
        [row1, row2, row3, row4, row5, row6, row7],
        location=(0, 0),
        finalize=True,
        element_justification="center",
        font="Helvetica 18")

    # Wait for events
    while True:
        event, values = window.read()
        if event == sg.WIN_CLOSED or event == 'Done':
            break
        if event == 'Compute':
            # Check for non-numeric entries for Central Wavelength
            if not values['-WAVEIN-'].isnumeric():
                window['-GRATOUT-'].update("")
                window['-WAVEOUT-'].update(f"Please Enter a Number")
                window['-TILTOUT-'].update("")
                window['-DEMAGOUT-'].update("")
                # Wait for next event
                continue

            # Convert wavelen to float, and check for valid range
            wavelen = float(values['-WAVEIN-'])
            if wavelen < 3000 or wavelen > 11000:
                window['-GRATOUT-'].update("")
                window['-WAVEOUT-'].update(f"Wavelength out of range")
                window['-TILTOUT-'].update("")
                window['-DEMAGOUT-'].update("")
                # Wait for next event
                continue

            # Compute the grating angle and anamorphic demagnification
            gpmm = float(values['Grat'].split(' - ')[1].split(' g/mm')[0])
            grangle, amag = compute_grangle(wavelen, gpmm)

            # Update the window with the calculated values
            window['-GRATOUT-'].update(values['Grat'])
            window['-WAVEOUT-'].update(f"{values['-WAVEIN-']} Å")
            window['-TILTOUT-'].update(f"{grangle+tgoffset:.2f} deg")
            window['-DEMAGOUT-'].update(f"{PIXSCALE*amag:.2f} pixels/arcsec")

    # All done, close window        
    window.close()


def deveny_grangle_maxgui():
    """Compute the desired grating angle given grating and central wavelength
    
    GUI version of deveny_grangle.  Uses PySimpleGUI.  Includes a drop-down
    menu for available gratings, and checks for a wavelength between 3000 and
    11,000 angstroms.  Uses the same subroutines as the CLI version and 
    produces the same results.

    This version also allows for the calcuation of the central wavelength given
    a grating angle
    """

    # Global variables -- available to grating_eqn function
    global gpmm, wavelen

    # Mechanical Offset
    tgoffset = 0.0

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
    row8 = [sg.Text("Enter Grating Tilt:"), 
            sg.Input(key="-TILTIN-", size=(6,1)), sg.Text("º")]
    row3 = [sg.Button("Compute Tilt"), sg.Button("Compute Wavelength"), 
            sg.Button("Done")]
    row4 = [sg.Text("             Grating: "), 
            sg.Text(size=(20,1), key="-GRATOUT-")]
    row5 = [sg.Text("  Central Wavelength: "), 
            sg.Text(size=(20,1), key="-WAVEOUT-")]
    row6 = [sg.Text("DeVeny Grating Tilt = "), 
            sg.Text(size=(20,1), key="-TILTOUT-")]
    row7 = [sg.Text("Slit demagnification (0.34\"/pixel): "),
            sg.Text(size=(15,1), key="-DEMAGOUT-")]

    # Create the Window
    window = sg.Window(
        "DeVeny Grating Angle Calculator",
        [row1, row2, row8, row3, row4, row5, row6, row7],
        location=(0, 0),
        finalize=True,
        element_justification="center",
        font="Helvetica 18")

    # Wait for events
    while True:
        event, values = window.read()
        if event == sg.WIN_CLOSED or event == 'Done':
            break

        elif event == 'Compute Tilt':
            # Check for non-numeric entries for Central Wavelength
            if not values['-WAVEIN-'].isnumeric():
                window['-GRATOUT-'].update("")
                window['-WAVEOUT-'].update(f"Please Enter a Number")
                window['-TILTOUT-'].update("")
                window['-DEMAGOUT-'].update("")
                window['-TILTIN-'].update("")
                # Wait for next event
                continue

            # Convert wavelen to float, and check for valid range
            wavelen = float(values['-WAVEIN-'])
            if wavelen < 3000 or wavelen > 11000:
                window['-GRATOUT-'].update("")
                window['-WAVEOUT-'].update(f"Wavelength out of range")
                window['-TILTOUT-'].update("")
                window['-DEMAGOUT-'].update("")
                window['-TILTIN-'].update("")
                # Wait for next event
                continue

            # Compute the grating angle and anamorphic demagnification
            gpmm = float(values['Grat'].split(' - ')[1].split(' g/mm')[0])
            grangle, amag = compute_grangle(wavelen, gpmm)

            # Update the window with the calculated values
            window['-GRATOUT-'].update(values['Grat'])
            window['-WAVEOUT-'].update(f"{values['-WAVEIN-']} Å")
            window['-TILTOUT-'].update(f"{grangle+tgoffset:.2f}º")
            window['-DEMAGOUT-'].update(f"{PIXSCALE*amag:.2f} pixels/arcsec")
            window['-TILTIN-'].update(f"{grangle+tgoffset:.2f}")
       
        elif event == 'Compute Wavelength':
            # Check for non-numeric entries for Grating Tilt
            if not check_float(values['-TILTIN-']):
                window['-GRATOUT-'].update("")
                window['-WAVEOUT-'].update("")
                window['-TILTOUT-'].update(f"Please Enter a Number")
                window['-DEMAGOUT-'].update("")
                window['-WAVEIN-'].update("")
                # Wait for next event
                continue

            # Convert wavelen to float, and check for valid range
            tilt = float(values['-TILTIN-'])
            if tilt < 0 or tilt > 48:
                window['-GRATOUT-'].update("")
                window['-WAVEOUT-'].update("")
                window['-TILTOUT-'].update(f"Tilt angle out of range")
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
            window['-TILTOUT-'].update(f"{tilt+tgoffset:.2f}º")
            window['-DEMAGOUT-'].update(f"{PIXSCALE*amag:.2f} pixels/arcsec")
            window['-WAVEIN-'].update(f"{wavelen:.0f}")
        else:
            print("Something funny happened... should never print.")

    # All done, close window        
    window.close()


def compute_grangle(wavelen, gpmm):
    """Compute the needed grating angle

    This function does the heavy lifting for computing the grating angle
    :param wavelen: The central wavelength in angstroms
    :param gpmm: The groove density of the grating in g/mm 
    :return: The computed grating angle
    """
    # Initial guess: 20º
    theta = np.deg2rad(20.) 

    # Call the newton method from scipy.optimize to solve the grating equation
    grangle = np.rad2deg(optimize.newton(grangle_eqn, theta))
    amag = deveny_amag(grangle)
    return (grangle, amag)


def grangle_eqn(theta):
    """The grating equation used to find the angle
    
    The scipy.optimize.newton() function looks for where this equation
    equals zero.
    """
    # DeVeny optical angles
    camcol = np.deg2rad(55.00)
    coll = np.deg2rad(10.00)

    gx = (np.sin((coll + theta)) + np.sin(coll + theta - camcol)) * \
         1.e7 / gpmm - wavelen
    return gx


def lambda_at_angle(theta, gpmm, radians=False):
    """Use the grating equation to compute the central wavelength given theta
    
    :param theta: The specified grating angle
    :param gpmm: The groove density of the grating in g/mm 
    :param radians: The input angle is in radians [Default: False]
    :return: The computed central wavelength
    """
   # DeVeny optical angles
    camcol = np.deg2rad(55.00)
    coll = np.deg2rad(10.00)

    if not radians:
        theta = np.deg2rad(theta)

    wavelen = (np.sin((coll + theta)) + np.sin(coll + theta - camcol)) * \
               1.e7 / gpmm
    return wavelen


def deveny_amag(grangle):
    """Computes the anamorphic demagnification of the slit given grangle"""
    # DeVeny optical angles
    collang = 10.0
    camcollang = 55.0
    alpha = np.deg2rad(grangle + collang)
    mbeta = np.deg2rad(camcollang - np.rad2deg(alpha))

    return(np.cos(alpha) / np.cos(mbeta))


def check_float(potential_float):
    """Simple funtion to check whether something is a float"""
    try:
        float(potential_float)
        return True
    except ValueError:
        return False


if __name__ == "__main__":
    """File is run from the command line

    Direct run takes either 0 (default:GUI) or 1 (choose: CLI or GUI) command
    line arguments.  Anything else gets an error message and end.
    """
    # Check for command line arguments
    if len(sys.argv) == 1:
        # If no command line argument, assume GUI
        deveny_grangle_gui()
    elif len(sys.argv) == 2:
        # If exactly one command line argument, check what it is.
        if sys.argv[1].lower() == 'cli':
            deveny_grangle_cli()
        elif sys.argv[1].lower() == 'gui':
            deveny_grangle_gui()
        elif sys.argv[1].lower() == 'max':
            deveny_grangle_maxgui()
        else:
            print("This routine only accepts 'CLI', 'GUI', " + \
                "or 'MAX' as arguments.")
    else:
        # Can't deal with more than one command line argument
        print("This routine accepts zero or one command line arguments.")
