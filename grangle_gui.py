# -*- coding: utf-8 -*-
#
#  This file is part of PyDeVeny.
#
#   This Source Code Form is subject to the terms of the Mozilla Public
#   License, v. 2.0. If a copy of the MPL was not distributed with this
#   file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
#  Created on 26-Jan-2021
#
#  @author: tbowers

"""PyDeVeny contains python ports of the various DeVeny IDL routines

Lowell Discovery Telescope (Lowell Observatory: Flagstaff, AZ)
http://www.lowell.edu
"""

# Built-In Libraries

# Numpy & SciPy & PySimpleGUI
import numpy as np
from scipy import optimize
import PySimpleGUI as sg

# Mechanical Offset
tgoffset = 0.0

def get_gpmm(grat_str):
    return float(grat_str.split(' - ')[1].split(' g/mm')[0])

def grangle_eqn(theta):
    """The grating equation used to find the angle
    """
    # DeVeny optical angles
    camcol = np.deg2rad(55.00)
    coll = np.deg2rad(10.00)

    gx = (np.sin((coll + theta)) + np.sin(coll + theta - camcol)) * 1.e7 / gpmm - wavelen
    return gx


def deveny_amag(grangle):
    """Computes the anamorphic demagnification of the slit given grangle
    """
    # DeVeny optical angles
    collang = 10.0
    camcollang = 55.0
    alpha = np.deg2rad(grangle + collang)
    mbeta = np.deg2rad(camcollang - np.rad2deg(alpha))

    return(np.cos(alpha) / np.cos(mbeta))


#sg.theme('dark purple 6')
sg.theme('light grey 1')

gratings = ["DV1 - 150 g/mm, 5000 Å",
            "DV2 - 300 g/mm, 4000 Å",
            "DV3 - 300 g/mm, 6750 Å",
            "DV4 - 400 g/mm, 8500 Å",
            "DV5 - 500 g/mm, 5500 Å",
            "DV6 - 600 g/mm, 4900 Å",
            "DV7 - 600 g/mm, 6750 Å",
            "DV8 - 831 g/mm, 8000 Å",
            "DV9 - 1200 g/mm, 5000 Å",
            "DV10 - 2160 g/mm, 5000 Å"
            ]


# Define the window layout
layout = [
    [sg.Text("Select Grating:"), 
    sg.Drop(values=(gratings), auto_size_text=True, 
            default_value=gratings[1],key="Grat")],
    [sg.Text("Enter Central Wavelength:"), sg.Input(key="-WAVEIN-", size=(6,1)), sg.Text("Å")],
    [sg.Button("Compute"), sg.Button("Done")],
    [sg.Text("             Grating: "), sg.Text(size=(20,1), key="-GRATOUT-")],
    [sg.Text("  Central Wavelength: "), sg.Text(size=(20,1), key="-WAVEOUT-")],
    [sg.Text("DeVeny Grating Tilt = "), sg.Text(size=(20,1), key="-TILTOUT-")],
    [sg.Text("Slit demagnification (0.34\"/pixel): "),
    sg.Text(size=(15,1), key="-DEMAGOUT-")]
]

# Create the form and show it without the plot
window = sg.Window(
    "Compute DeVeny Grating Angle",
    layout,
    location=(0, 0),
    finalize=True,
    element_justification="center",
    font="Helvetica 18",
)

while True:
    event, values = window.read()
    if event == sg.WIN_CLOSED or event == 'Done':
        break
    if event == 'Compute':
        if not values['-WAVEIN-'].isnumeric():
            #print("Central Wavelength is not a number.  Try again.")
            window['-GRATOUT-'].update("")
            window['-WAVEOUT-'].update(f"Please Enter a Number")
            window['-TILTOUT-'].update("")
            window['-DEMAGOUT-'].update("")
            continue

        wavelen = float(values['-WAVEIN-'])
        gpmm = get_gpmm(values['Grat'])

        # Call the newton method from scipy.optimize to solve the grating equation
        theta = np.deg2rad(20.)
        grangle = optimize.newton(grangle_eqn, theta)
        grangle = np.rad2deg(grangle)
        amag = deveny_amag(grangle)

        window['-GRATOUT-'].update(values['Grat'])
        window['-WAVEOUT-'].update(f"{values['-WAVEIN-']} Å")
        window['-TILTOUT-'].update(f"{grangle+tgoffset:.2f} deg")
        window['-DEMAGOUT-'].update(f"{2.94*amag:.2f} pixels/arcsec")
        #print(values["Grat"])
        #print("We're computing!")



window.close()