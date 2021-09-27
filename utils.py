# -*- coding: utf-8 -*-
#
#  This file is part of LDTObserverTools.
#
#   This Source Code Form is subject to the terms of the Mozilla Public
#   License, v. 2.0. If a copy of the MPL was not distributed with this
#   file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
#  Created on 01-Feb-2021
#
#  @author: tbowers

"""LDTObserverTools contains python ports of the various DeVeny IDL routines

Lowell Discovery Telescope (Lowell Observatory: Flagstaff, AZ)
http://www.lowell.edu

This file contains the dfocus routine for computing the required collimator
focus for the DeVeny Spectrograph based on a focus sequence completed by the
DeVeny LOUI.
"""

# Built-In Libraries
import warnings

# 3rd-Party Libraries
import numpy as np
from scipy import optimize

# Local Libraries

# CONSTANTS


def gaussfit(x, y, nterms=3, estimates=None):
    """gaussfit Function similar to IDL's GAUSSFIT

    Big caveat: as implemented, can only estimate the initial parameters for
    POSITIVE gaussians (emission), and cannot correctly estimate parameters
    for negative (absorption) gaussians.  The function will still happily fit
    a negative gaussian if given the proper estimates.

    Utilizes scipy.optimize.curvefit and the helper function gaussfit_func
    below.

    Parameters
    ----------
    x : [type]
        [description]
    y : [type]
        [description]
    nterms : int, optional
        [description], by default 3
    estimates : [type], optional
        [description], by default None

    Returns
    -------
    [type]
        [description]

    Raises
    ------
    ValueError
        [description]
    ValueError
        [description]
    """
    global nterms_gf
    nterms_gf = nterms

    if nterms < 3 or nterms > 6:
        raise ValueError(f"{nterms} is an invalid number of terms.")

    # This block is for estimating the parameters if none given.
    if estimates is None:
        # Subtract a linear term if nterm == 5 or 6 or constant for nterm == 4
        if nterms > 3:
            p = np.polyfit(x, y, 0 if nterms == 4 else 1)
            y_modified = y - np.polyval(p, x)
        # Do nothing if nterm == 3
        else:
            y_modified = y

        # Find the estimates of a0, a1, a2:
        dx = np.diff(x)[0]
        a0 = np.max(y_modified)
        a1 = x[0] + first_moment_1d(y_modified)*dx
        # Use points where value > (1/e)*max to estimate width
        s_idx = np.where(y_modified > a0/np.e)
        a2 = np.abs(x[s_idx][-1] - x[s_idx][0])/2.
        # Need to estimate the width at 1/e from the peak...

        # Construct the estimates list
        estimates = [a0, a1, a2]
        if nterms > 3:
            estimates = estimates + list(np.flip(p))
        if nterms == 6:
            estimates.append(0.)

    # Else, make sure the number of estimate values equals nterms
    else:
        if len(estimates) != nterms:
            raise ValueError("Estimate array must contain NTERMS elements.")

    # Pad out the estimates list to make the fitting function happy
    if nterms == 3:
        estimates = estimates + [0., 0., 0.]
    elif nterms == 4:
        estimates = estimates + [0., 0.]
    elif nterms == 5:
        estimates = estimates + [0.]
        
    # print(estimates)

    return  optimize.curve_fit(gaussfit_func, x, y, p0=estimates)



def gaussfit_func(x, a0, a1, a2, a3, a4, a5):
    """gaussfit_func Gaussian Function

    [extended_summary]

    Parameters
    ----------
    x : `array`
        X values over which to compute the gaussian
    a0 : `float`
        Gaussian amplitude
    a1 : `float`
        Gaussian mean (mu)
    a2 : `float`
        Gaussian width (sigma)
    a3 : `float`
        Baseline atop which the Gaussian sits
    a4 : `float`
        Slope of the baseline atop which the Gaussian sits
    a5 : `float`
        Quadratic term of the baseline atop which the Gaussian sits

    Returns
    -------
    `array`
        The Y values of the Gaussian corresponding to X
    """
    global nterms_gf

    # Silence RuntimeWarning for overflow, this function only
    warnings.simplefilter('ignore', RuntimeWarning)
    z = (x - a1) / a2

    if nterms_gf == 3:
        return a0 * np.exp(-z**2 / 2.)
    if nterms_gf == 4:
        return a0 * np.exp(-z**2 / 2.) + a3
    if nterms_gf == 5:
        return a0 * np.exp(-z**2 / 2.) + a3 + a4*x
    if nterms_gf == 6:
        return a0 * np.exp(-z**2 / 2.) + a3 + a4*x + a5*x**2


def first_moment_1d(line):
    """first_moment_1d Returns the 1st moment of line

    [extended_summary]

    Parameters
    ----------
    line : `array`
        1-dimensional array to find the 1st moment of

    Returns
    -------
    `float`
        The first moment of the input array relative to element #
    """
    # Only use positive values -- set negative values to zero
    line[np.where(line < 0)] = 0

    # Make the counting array
    yy = np.arange(len(line))

    # Return the first moment
    return np.sum(yy * line) / np.sum(line)

