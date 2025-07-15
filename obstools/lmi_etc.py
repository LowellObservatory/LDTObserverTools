# -*- coding: utf-8 -*-
#
#  This file is part of LDTObserverTools.
#
#   This Source Code Form is subject to the terms of the Mozilla Public
#   License, v. 2.0. If a copy of the MPL was not distributed with this
#   file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
#  Created on 29-Dec-2021
#  GUI Added on 11-Jul-2025
#
#  @author: tbowers
# pylint: disable=no-name-in-module

"""LMI Exposure Time Calculator Module

LDTObserverTools contains python ports of various LDT Observer Tools

Lowell Discovery Telescope (Lowell Observatory: Flagstaff, AZ)
https://lowell.edu

This file contains the LMI exposure time calculator routine (ported from PHP)
from Phil Massey's webpage.
http://www2.lowell.edu/users/massey/LMI/etc_calc.php

These routines are used for computing the required exposure times for LMI based
on various requirements.

The values for ``Star20`` are the count rates in `e-`/`sec`/`image` at airmass=0
for a 20th magnitude star measured with a radius = 1.4 x FWHM in pixels.

.. warning::

    The LMI-specific pixel scale, gain, and read noise are hard-coded into
    this module.

    
.. note::

    A GUI wrapper for these functions is forthcoming.
"""

# Built-In Libraries
import argparse
import dataclasses
import sys

# 3rd-Party Libraries
import astropy.table
import numpy as np
from PyQt6.QtWidgets import QApplication, QMainWindow

# Local Libraries
from obstools import utils

# Constants
SCALE = 0.12  # "/pix
READ_NOISE = 6.0  # e-
GAIN = 2.89  # e-/ADU
BIAS = 1050  # ADU (approx) for 2x2 binning


# Define the API as the "User-Interface Computation Routines"
# __all__ = [
#     "exptime_given_snr_mag",
#     "exptime_given_peak_mag",
#     "snr_given_exptime_mag",
#     "mag_given_snr_exptime",
#     "peak_counts",
# ]


@dataclasses.dataclass
class ETCData:
    """ETC Data Class

    snr : :obj:`float`
        Desired signal-to-noise ratio
    mag : :obj:`float`
        Magnitude in the band of the star desired
    exptime : :obj:`float`
        User-defined exposure time (seconds)
    peak : :obj:`float`
        Desired peak count level on the CCD (e-)
    airmass : :obj:`float`
        Airmass at which the observation will take place
    band : :obj:`str`
        The LMI filter for which to perform the calculation
    phase : :obj:`float`
        Moon phase (0-14)
    seeing : :obj:`float`
        Size of the seeing disk (arcsec)
    binning : :obj:`int`, optional
        Binning of the CCD  (Default: 2)
    """

    snr: float = None
    mag: float = None
    exptime: float = None
    peak: float = None
    airmass: float = None
    band: str = None
    phase: float = None
    seeing: float = None
    binning: int = 2


@dataclasses.dataclass
class BandData:
    """Band-Specific Data Class

    filter : :obj:`str`
        Name of the filter
    star20 : :obj:`float`
        Expected counts/sec from a 20th magnitude star
    extinction : :obj:`float`
        The in-band extinction in magnitudes per airmass
    sky0 : :obj:`float`
        Constant sky brightness term (magnitudes)
    sky1 : :obj:`float`
        Linear in moon phase sky brightness term (magnitudes)
    sky2 : :obj:`float`
        Quadratic in moon phase sky brightness term (magnitudes)
    """

    filter: str = None
    star20: float = None
    extinction: float = None
    sky0: float = None
    sky1: float = None
    sky2: float = None


# User-Interface Computation Routines ========================================#
def exptime_given_snr_mag(input_data: ETCData) -> float:
    """Compute the exposure time given SNR and magnitude

    Given a desired signal-to-noise ratio and stellar magnitude, compute the
    exposure time required for a particular LMI Filter, moon phase, seeing and
    CCD binning.

    Parameters
    ----------
    input_data : :class:`ETCData`
        The input data needed for this calculation, placed in an
        :class:`ETCData` class.

    Returns
    -------
    :obj:`float`
        The desired exposure time in seconds
    """
    # Check the inputs
    check_etc_inputs(input_data)

    # Load the necessary counts information
    star_counts = star_counts_per_sec(input_data)
    sky_counts = sky_counts_per_sec_in_aperture(input_data)
    read_counts = read_noise_in_aperture(input_data)

    # Do the computation
    k_a = star_counts**2
    k_b = -input_data.snr**2 * (star_counts + sky_counts)
    k_c = -input_data.snr**2 * read_counts**2
    return (-k_b + np.sqrt(k_b * k_b - 4.0 * k_a * k_c)) / (2.0 * k_a)


def exptime_given_peak_mag(input_data: ETCData) -> float:
    """Compute the exposure time given peak and mag

    Given a desired peak count level on the CCD and stellar magnitude, compute
    the exposure time required for a particular LMI Filter, moon phase, seeing
    and CCD binning.

    Parameters
    ----------
    input_data : :class:`ETCData`
        The input data needed for this calculation, placed in an
        :class:`ETCData` class.

    Returns
    -------
    :obj:`float`
        The desired exposure time in seconds
    """
    # Check the inputs
    check_etc_inputs(input_data)

    # Load the necessary counts information
    star_counts = star_counts_per_sec(input_data)
    sky_counts = sky_counts_per_sec_in_aperture(input_data)

    # Do the computation
    fwhm = input_data.seeing / (SCALE * input_data.binning)
    sky_count_per_pixel_per_sec = sky_counts / pixels_in_aperture(input_data)
    return (input_data.peak - BIAS * GAIN) / (
        star_counts / (1.13 * fwhm**2) + sky_count_per_pixel_per_sec
    )


def snr_given_exptime_mag(input_data: ETCData) -> float:
    """Compute the SNR given exposure time and magnitude

    Given a desired exposure time and stellar magnitude, compute the resulting
    signal-to-noise ratio for a particular LMI Filter, moon phase, seeing and
    CCD binning.

    Parameters
    ----------
    input_data : :class:`ETCData`
        The input data needed for this calculation, placed in an
        :class:`ETCData` class.

    Returns
    -------
    :obj:`float`
        The desired signal-to-noise ratio
    """
    # Check the inputs
    check_etc_inputs(input_data=input_data)

    # Load the necessary counts information
    star_counts = star_counts_per_sec(input_data)
    sky_counts = sky_counts_per_sec_in_aperture(input_data)
    read_counts = read_noise_in_aperture(input_data)

    # Do the computation
    signal = star_counts * input_data.exptime
    noise = np.sqrt(
        star_counts * input_data.exptime
        + sky_counts * input_data.exptime
        + read_counts**2
    )
    return signal / noise


def mag_given_snr_exptime(input_data: ETCData) -> float:
    """Compute the magnitude given SNR and exposure time

    Given a desired signal-to-noise ratio and exposure time, compute the
    limiting magnitude for a particular LMI Filter, moon phase, seeing and
    CCD binning.

    Parameters
    ----------
    input_data : :class:`ETCData`
        The input data needed for this calculation, placed in an
        :class:`ETCData` class.

    Returns
    -------
    :obj:`float`
        The limiting stellar magnitude
    """
    # Check the inputs
    check_etc_inputs(input_data)

    # Load the necessary counts information
    band_info = get_band_values(input_data.band)
    sky_counts = sky_counts_per_sec_in_aperture(input_data)
    read_counts = read_noise_in_aperture(input_data)

    # Do the computation
    k_a = input_data.exptime**2
    k_b = -input_data.snr**2 * input_data.exptime
    k_c = -input_data.snr**2 * (
        sky_counts * input_data.exptime + read_counts * read_counts
    )
    cts_from_star_per_sec = (-k_b + np.sqrt(k_b**2 - 4.0 * k_a * k_c)) / (2.0 * k_a)
    mag_raw = -2.5 * np.log10(cts_from_star_per_sec / band_info.star20) + 20.0
    return mag_raw - band_info.extinction * input_data.airmass


def peak_counts(input_data: ETCData) -> float:
    """Compute the peak counts on the CCD for an exptime and mag

    Given a desired exposure time and stellar magnitude, compute the resulting
    counts on the CCD for a particular LMI Filter, moon phase, seeing and
    CCD binning.

    Parameters
    ----------
    input_data : :class:`ETCData`
        The input data needed for this calculation, placed in an
        :class:`ETCData` class.

    Returns
    -------
    :obj:`float`
        The desired counts on the CCD (e-)
    """
    # Load the necessary counts information
    star_counts = star_counts_per_sec(input_data)
    sky_counts = sky_counts_per_sec_in_aperture(input_data)

    # Do the computation
    fwhm = input_data.seeing / (SCALE * input_data.binning)
    sky_count_per_pixel_per_sec = sky_counts / pixels_in_aperture(input_data)
    return (
        star_counts / (1.13 * fwhm**2) + sky_count_per_pixel_per_sec
    ) * input_data.exptime + BIAS * GAIN


# Helper Routines (Alphabetical) =============================================#
def check_etc_inputs(input_data: ETCData):
    """Check the ETC inputs for valid values

    Does a cursory check on the ETC inputs for proper range, etc.  These are
    not exhaustive checks (i.e. checking for proper type on all values), but
    a good starting point nonetheless.

    Parameters
    ----------
    input_data : :class:`ETCData`
        The input data needed for this calculation, placed in an
        :class:`ETCData` class.

    Raises
    ------
    :obj:`ValueError`
        If any of the inputs are deemed to be out of range
    """
    if input_data.airmass < 1 or input_data.airmass > 3:
        raise ValueError(f"Invalid airmass specified: {input_data.airmass}")
    if input_data.phase < 0 or input_data.phase > 14:
        raise ValueError(f"Invalid moon phase specified: {input_data.phase}")
    if input_data.seeing < 0.5 or input_data.seeing > 3.0:
        raise ValueError(f"Invalid seeing specified: {input_data.seeing}")
    if (
        not isinstance(input_data.binning, int)
        or input_data.binning < 1
        or input_data.binning > 4
    ):
        raise ValueError(f"Invalid binning specified: {input_data.binning}")
    if input_data.exptime and (input_data.exptime < 0.001 or input_data.exptime > 1200):
        raise ValueError(f"Invalid exposure time specified: {input_data.exptime}")
    if input_data.mag and (input_data.mag < -1 or input_data.mag > 28):
        raise ValueError(f"Invalid stellar magnitude specified: {input_data.mag}")
    if input_data.snr and input_data.snr < 0.1:
        raise ValueError(f"Invalid signal-to-noise specified: {input_data.snr}")


def get_band_values(band: str) -> BandData:
    """Return the band-specific star and sky values

    Pull the correct row from ``etc_filter_info.ecsv`` containing the star count
    and sky parameters (brightness and extinction).

    Parameters
    ----------
    band : :obj:`str`
        The LMI filter to use

    Returns
    -------
    :class:`BandData`
        The class representation of the table row.  If the band is
        improperly specified (i.e. is not in the table), raise an error.
    """
    # Read in the table, and index the filter column
    table = astropy.table.Table.read(utils.DATA / "etc_filter_info.ecsv")
    table.add_index("Filter")

    band_info = BandData()
    try:
        # Extract the row, and return it as a BandData class
        row = table.loc[band]
        for key, val in zip(row.colnames, row):
            setattr(band_info, key, val)
    except KeyError as err:
        # Raise an error
        raise ValueError(f"Improper LMI band provided: {band}") from err
    return band_info


def pixels_in_aperture(input_data: ETCData) -> float:
    """Number of pixels in the measuring aperture

    Counts the number of pixels in the measuring aperture, based on the seeing
    and CCD binning scheme.

    .. note::
        The minimum value of the return value ``N_pix`` is 9, which corresponds
        to a FWHM of 2.54 pixels.  For 2x2 binning, this occurs at a seeing of
        0.61" (3x3 binning = 0.91", 4x4 binning = 1.22").

    Parameters
    ----------
    input_data : :class:`ETCData`
        The input data needed for this calculation, placed in an
        :class:`ETCData` class.

    Returns
    -------
    :obj:`float`
        Equivalent number of pixels within the measuring aperture
    """
    fwhm = input_data.seeing / (SCALE * input_data.binning)
    return np.max([1.4 * fwhm**2, 9.0])


def read_noise_in_aperture(input_data: ETCData) -> float:
    """Calculate read-noise contribution in the aperture

    Compute the read-noise contribution to the measuring aperture by
    multiplying the read noise per pixel by the square root of the number
    of pixels.

    Parameters
    ----------
    input_data : :class:`ETCData`
        The input data needed for this calculation, placed in an
        :class:`ETCData` class.

    Returns
    -------
    :obj:`float`
        The read-noise contribution to the photometry aperture
    """
    return READ_NOISE * np.sqrt(pixels_in_aperture(input_data))


def sky_counts_per_sec_in_aperture(input_data: ETCData) -> float:
    """Determine sky counts per aperture per second

    [extended_summary]

    Parameters
    ----------
    input_data : :class:`ETCData`
        The input data needed for this calculation, placed in an
        :class:`ETCData` class.

    Returns
    -------
    :obj:`float`
        Sky counts per second in the aperture
    """
    band_info = get_band_values(input_data.band)
    sky_brightness_per_arcsec2 = (
        band_info.sky0
        + band_info.sky1 * input_data.phase
        + band_info.sky2 * input_data.phase**2
    )
    sky_count_per_arcsec2_per_sec = band_info.star20 * np.power(
        10, -((sky_brightness_per_arcsec2 - 20) / 2.5)
    )
    rscale = SCALE * input_data.binning
    sky_count_per_pixel_per_sec = sky_count_per_arcsec2_per_sec * rscale * rscale
    return pixels_in_aperture(input_data) * sky_count_per_pixel_per_sec


def star_counts_per_sec(input_data: ETCData) -> float:
    """Compute the counts per second from a star

    Compute the counts per second from a star given a band, magnitude, and
    airmass.

    Parameters
    ----------
    input_data : :class:`ETCData`
        The input data needed for this calculation, placed in an
        :class:`ETCData` class.

    Returns
    -------
    :obj:`float`
        Number of counts per second for the described star
    """
    band_info = get_band_values(input_data.band)
    mag_corrected = input_data.mag + band_info.extinction * input_data.airmass
    return band_info.star20 * np.power(10, -((mag_corrected - 20) / 2.5))


# Command Line Script Infrastructure (borrowed from PypeIt) ==================#
class LmiEtc(utils.ScriptBase):
    """Script class for ``lmi_etc`` tool

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
            description="LMI Exposure Time Calculator", width=width
        )
        return parser

    @staticmethod
    def main(args):
        """Main Driver

        Simple function that calls the main driver function.
        """
        # Giddy up!

        # You need one (and only one) QApplication instance per application.
        # Pass in sys.argv to allow command line arguments for your app.
        # If you know you won't use command line arguments QApplication([]) works too.
        app = QApplication(sys.argv)

        # Create a Qt widget, which will be our window.
        window = QMainWindow()
        window.show()  # IMPORTANT!!!!! Windows are hidden by default.

        # Start the event loop.
        app.exec()
