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

This file contains various utility functions needed by other routines in this
package.
"""

# Built-In Libraries
import pathlib
import warnings

# 3rd-Party Libraries
from pkg_resources import resource_filename
import numpy as np
import scipy.optimize


# Local Libraries

# CONSTANTS
SG_THEME = "light grey 1"
CONFIG = pathlib.Path(resource_filename("obstools", "config"))


def check_float(potential_float):
    """Simple funtion to check whether something is a float

    Parameters
    ----------
    potential_float : :obj:`~typing.Any`
        Value to check for float

    Returns
    -------
    :obj:`bool`
        Whether it am or it ain't a :obj:`float`.
    """
    try:
        float(potential_float)
        return True
    except ValueError:
        return False


def gaussfit(x, y, nterms=3, estimates=None, bounds=None, debug=False):
    """Function similar to IDL's GAUSSFIT

    Big caveat: as implemented, can only estimate the initial parameters for
    POSITIVE gaussians (emission), and cannot correctly estimate parameters
    for negative (absorption) gaussians.  The function will still happily fit
    a negative gaussian if given the proper estimates.

    Utilizes :func:`scipy.optimize.curve_fit` and the helper function
    :func:`gaussian_function` below.

    Parameters
    ----------
    x : :obj:`~numpy.ndarray`
        Input abscissa values for the fitting
    y : :obj:`~numpy.ndarray`
        Input ordinate values for the fitting
    nterms : :obj:`int`, optional
        Number of terms to use in the Gaussian fit (Default: 3)
    estimates : :obj:`~numpy.ndarray`, optional
        Estimated values for the ``nterms`` parameters
        (see :func:`~scipy.optimize.curve_fit`; Default: None)
    bounds : :obj:`~numpy.ndarray`, optional
        Bounds on the ``nterms`` parameters
        (see :func:`~scipy.optimize.curve_fit`; Default: None)
    debug : :obj:`bool`, optional
        Print debugging statements.  (Defualt: False)

    Returns
    -------
    popt : :obj:`~numpy.ndarray`
        Optimal values for the parameters (see :func:`~scipy.optimize.curve_fit`)

    pcov : :obj:`~numpy.ndarray`
        The estimated covariance of ``popt`` (see :func:`~scipy.optimize.curve_fit`)
    """

    if nterms < 3 or nterms > 6:
        raise ValueError(f"{nterms} is an invalid number of terms.")

    if bounds is not None and len(bounds[0]) != nterms:
        raise ValueError("Bounds array must contain NTERMS elements.")

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
        a1 = x[0] + first_moment_1d(y_modified) * dx
        # Use points where value > (1/e)*max to estimate width
        s_idx = np.where(y_modified > a0 / np.e)
        a2 = np.abs(x[s_idx][-1] - x[s_idx][0]) / 2.0
        if debug:
            print(f"Estimated values: a0={a0:.1f}  a1={a1:.2f}  a2={a2:.2f}")

        # Construct the estimates list
        estimates = [a0, a1, a2]

        # Check estimate against bounds
        if bounds is not None:
            for i, est in enumerate(estimates):
                est = bounds[0][i] if est < bounds[0][i] else est
                est = bounds[1][i] if est > bounds[1][i] else est
                estimates[i] = est

        if nterms > 3:
            estimates = estimates + list(np.flip(p))
        if nterms == 6:
            estimates.append(0.0)

    # Else, make sure the number of estimate values equals nterms
    else:
        if len(estimates) != nterms:
            raise ValueError("Estimate array must contain NTERMS elements.")

    if bounds is None:
        bounds = (-np.inf, np.inf)
    if debug:
        print(bounds)

    popt, pcov = scipy.optimize.curve_fit(
        gaussian_function, x, y, p0=estimates, bounds=bounds, ftol=1e-6
    )

    if debug:
        print(f"Estimated/Fit Width: {a2} / {popt[2]}")

    return popt, pcov


def gaussian_function(
    x: np.ndarray,
    a0: float,
    a1: float,
    a2: float,
    a3: float = 0.0,
    a4: float = 0.0,
    a5: float = 0.0,
):
    """Gaussian Function

    Construct a basic Gaussian using at least 3, but up to 6 parameters.

    Parameters
    ----------
    x : :obj:`~numpy.ndarray`
        X values over which to compute the gaussian
    a0 : :obj:`float`
        Gaussian amplitude
    a1 : :obj:`float`
        Gaussian mean (mu)
    a2 : :obj:`float`
        Gaussian width (sigma)
    a3 : :obj:`float`, optional
        Baseline atop which the Gaussian sits.  (Default: 0)
    a4 : :obj:`float`, optional
        Slope of the baseline atop which the Gaussian sits.  (Default: 0)
    a5 : :obj:`float`, optional
        Quadratic term of the baseline atop which the Gaussian sits.
        (Default: 0)

    Returns
    -------
    :obj:`~numpy.ndarray`
        The Y values of the Gaussian corresponding to X
    """
    # Silence RuntimeWarning for overflow, this function only
    warnings.simplefilter("ignore", RuntimeWarning)
    z = (x - a1) / a2

    return a0 * np.exp(-(z**2) / 2.0) + a3 + a4 * x + a5 * x**2


def first_moment_1d(line):
    """Returns the 1st moment of line

    Parameters
    ----------
    line : :obj:`~numpy.ndarray`
        1-dimensional array to find the 1st moment of

    Returns
    -------
    :obj:`float`
        The first moment of the input array relative to element #
    """
    # Only use positive values -- set negative values to zero
    line[np.where(line < 0)] = 0

    # Make the counting array
    yy = np.arange(len(line))

    # Return the first moment
    return np.sum(yy * line) / np.sum(line)


def good_poly(x, y, order, thresh, return_full=False):
    """Robust fitting of a polynomial to data

    This is a python port of an IDL routine written years ago by M. Buie.

    This is a multi-pass fitting routine that fits a fixed order polynomial
    to the input data.  After each pass, the scatter of the fit relative
    to the fitted line is computed.  Each point is examined to see if it
    falls beyond THRESH sigma from the line.  If is does, it is removed
    from the data and the fit is tried again.  This will make up to two
    attempts to remove bad data.

    Written in IDL 1991-1998, Marc W. Buie, Lowell Observatory

    Parameters
    ----------
    x : :obj:`~numpy.ndarray`
        Input dataset, independant values.
    y : :obj:`~numpy.ndarray`
        Input dataset, dependant values.
    order : :obj:`int`
        Order of the polynomial fit (linear = 1).
    thresh : :obj:`float`
        Sigma threshold for removing outliers.
    return_full : :obj:`bool`, optional
        If True, also return:
            yfit : Fitted values for y that match the input vector.
            newx : X values from input that were considered good.
            newy : Y values from input that were considered good.

    Returns
    -------
    :obj:`~numpy.ndarray`
        Array of fit parameters, as in :func:`numpy.polyfit`.
    Also, optionally, the ``return_full`` bits
    """
    # Make copies to not mess up the inputs
    xx = x
    yy = y

    # Filter out NaNs
    if False in (good_idx := np.logical_and(~np.isnan(x), ~np.isnan(y))):
        xx = xx[good_idx]
        yy = yy[good_idx]

    if (array_length := len(xx)) == 0:
        return warn_and_return_zeros(return_full, x, xx, yy, order)

    # Check for fewer data points than the requested polynomial order
    if array_length < order:
        coeff = np.zeros(order + 1)
        if array_length != 1:
            sigma = np.std(yy)
            coeff[0] = np.mean(yy)
        else:
            coeff[0] = yy[0]
            sigma = yy[0]
            sigma = 1.0 if sigma == 0.0 else sigma
        print("Not enough data to support even a non-robust polynomial fit.")
        if return_full:
            yfit = [coeff[0]] * len(x)
            return coeff, yfit, xx, yy
        return coeff

    # Initial fit with all the data.
    coeff = np.polyfit(xx, yy, order)
    yfit = np.polyval(coeff, xx)
    flat = (yy - yfit) + np.sum(yfit) / array_length
    mean, sigma = np.mean(flat), np.std(flat)

    # Remove all points beyond threshold sigma
    good = np.where(np.abs(flat - mean) < thresh * sigma)
    nbad = array_length - len(good)
    xx, yy = xx[good], yy[good]
    if (array_length := len(xx)) == 0:
        return warn_and_return_zeros(return_full, x, xx, yy, order)

    # Do a second pass if there were any bad points removed
    if nbad != 0:
        coeff = np.polyfit(xx, yy, order)
        yfit = np.polyval(coeff, xx)
        flat = (yy - yfit) + np.sum(yfit) / array_length
        mean, sigma = np.mean(flat), np.std(flat)

        # Remove all points beyond threshold sigma
        good = np.where(np.abs(flat - mean) < thresh * sigma)
        nbad = array_length - len(good)
        xx, yy = xx[good], yy[good]
        if (array_length := len(xx)) == 0:
            return warn_and_return_zeros(return_full, x, xx, yy, order)

    # Do a third pass if there were any more bad points removed
    if nbad != 0:
        coeff = np.polyfit(xx, yy, order)
        yfit = np.polyval(coeff, xx)
        flat = (yy - yfit) + np.sum(yfit) / array_length
        mean, sigma = np.mean(flat), np.std(flat)

    # Check that the fit coefficients are finite:
    if not np.all(np.isfinite(coeff)):
        return warn_and_return_zeros(return_full, x, xx, yy, order)

    if return_full:
        return coeff, yfit, xx, yy
    return coeff


def warn_and_return_zeros(return_full: bool, x, xx, yy, order, raise_warn=False):
    """Set warning and return zeroes from :func:`good_poly`

    This function is a DRY.  Since this block is used several times in
    :func:`good_poly`, separate out as a function.

    Parameters
    ----------
    return_full : :obj:`bool`
        Return a bunch of stuff
    x :  :obj:`~numpy.ndarray`
        [description]
    xx :  :obj:`~numpy.ndarray`
        [description]
    yy :  :obj:`~numpy.ndarray`
        [description]
    order : :obj:`int`
        The order of the polynomial fit
    raise_warn : :obj:`bool`, optional
        Actually raise the warning this function is meant to  [Default: False]

    Returns
    -------
    :obj:`~numpy.ndarray`
        An array of zeros of the proper length
    """
    if raise_warn:
        warnings.warn("No good values to fit, return zeros.", UserWarning)
    if return_full:
        yfit = [0] * len(x)
        return [0] * (order + 1), yfit, xx, yy
    return [0] * (order + 1)
