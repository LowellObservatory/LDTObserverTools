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
"""

# Built-In Libraries

# Numpy & SciPy
import numpy as np
from scipy import optimize


def deveny_grangle():
    """Compute the desired grating angle given grating and central wavelength
    
    Takes no arguments, and prints output to screen.
    """
    # Global variables
    global gpmm, wavelen

    # Mechanical Offset
    tgoffset = 0.0

    # Get input from user
    print(" Enter grating resolution (g/mm):")    
    gpmm = float(input())
    print(" Enter central wavelength (A):")
    wavelen = float(input())

    # Call the newton method from scipy.optimize to solve the grating equation
    theta = np.deg2rad(20.)
    grangle = optimize.newton(grangle_eqn, theta)
    grangle = np.rad2deg(grangle)
    amag = deveny_amag(grangle)

    print(f"\n Grating: {gpmm:.0f} g/mm")
    print(f" Central Wavelength: {wavelen} A")
    print(f" DeVeny grating tilt = {grangle+tgoffset:.2f} deg")
    print(f" Slit demagnification (pixels/arcsec, 0.34 arcsec/pixel): {2.94*amag:.2f}\n")


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


if __name__ == "__main__":
    deveny_grangle()
