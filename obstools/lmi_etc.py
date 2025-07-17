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
# pylint: disable=c-extension-no-member

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
from PyQt6 import QtCore
from PyQt6 import QtGui
from PyQt6 import QtWidgets

# Local Libraries
from obstools import utils
from obstools.ETCWindow import Ui_MainWindow

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

    Attributes
    ----------
    snr : :obj:`float`
        Desired signal-to-noise ratio
    mag : :obj:`float`
        Magnitude in the band of the star desired
    exptime : :obj:`float`
        User-defined exposure time (seconds)
    peak : :obj:`float`
        Desired peak count level on the CCD (e-)
    airmass : :obj:`float`
        Airmass at which the observation will take place  (Default: 1.0)
    band : :obj:`str`
        The LMI filter for which to perform the calculation  (Default: "V")
    phase : :obj:`float`
        Moon phase (0-14)  (Default: 0)
    seeing : :obj:`float`
        Size of the seeing disk (arcsec)  (Default: 1.0")
    binning : :obj:`int`, optional
        Binning of the CCD  (Default: 2)
    """

    snr: float = None
    mag: float = None
    exptime: float = None
    peak: float = None
    airmass: float = 1.0
    band: str = "V"
    phase: float = 0
    seeing: float = 1.0
    binning: int = 2


@dataclasses.dataclass
class AuxData:
    """Auxillary Data Class

    Attributes
    ----------
    n_pix : :obj:`float`
        Number of pixels in the measuing aperture
    sky_bright : :obj:`float`
        The in-band sky brightness (e- / pix)
    num_e_star : :obj:`float`
        Total number of electrons from the star in the aperture
    peak_e_star : :obj:`float`
        Peak number of electrons in the center of the image
    noise_star : :obj:`float`
        Noise contribution (in e-) from the star (Poisson)
    noise_sky : :obj:`float`
        Noise contribution (in e-) from the sky (Poisson)
    noise_ccd : :obj:`float`
        Noise contribution (in e-) from the CCD (readout)
    """

    n_pix: float = 0
    sky_bright: float = 0
    num_e_star: float = 0
    peak_e_star: float = 0
    noise_star: float = 0
    noise_sky: float = 0
    noise_ccd: float = 0


@dataclasses.dataclass
class BandData:
    """Band-Specific Data Class

    Attributes
    ----------
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

    # Do the computation, quadratic style
    k_a = star_counts**2
    k_b = -input_data.snr**2 * (star_counts + sky_counts)
    k_c = -input_data.snr**2 * read_counts**2
    return (-k_b + np.sqrt(k_b**2 - 4.0 * k_a * k_c)) / (2.0 * k_a)


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
    noise = np.sqrt(signal + sky_counts * input_data.exptime + read_counts**2)
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

    # Do the computation, quadratic style
    k_a = input_data.exptime**2
    k_b = -input_data.snr**2 * input_data.exptime
    k_c = -input_data.snr**2 * (sky_counts * input_data.exptime + read_counts**2)
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
        for name, value in zip(row.colnames, row):
            setattr(band_info, name.lower(), value)
    except (KeyError, AttributeError) as err:
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
    print(band_info)
    mag_corrected = input_data.mag + band_info.extinction * input_data.airmass
    return band_info.star20 * np.power(10, -((mag_corrected - 20) / 2.5))


# GUI Class ==================================================================#
class ETCWindow(QtWidgets.QMainWindow, Ui_MainWindow):
    """Exposure Time Calculator Main Window Class

    The UI is defined in ETCWindow.ui and translated (via pyuic6) into python
    in ETCWindow.py.  This class inherits the UI and defines the various
    actions needed to compute ETCs from the GUI inputs.
    """

    def __init__(self):
        super().__init__()
        self.setupUi(self)
        self.show()

        # Connect buttons to actions
        self.exitButton.pressed.connect(self.exit_button_clicked)
        self.computeButton.pressed.connect(self.compute_button_clicked)
        self.radioExpSnMag.toggled.connect(
            lambda checked: self.set_stacked_page(0, checked)
        )
        self.radioExpPkMag.clicked.connect(
            lambda checked: self.set_stacked_page(1, checked)
        )
        self.radioSnExpMag.clicked.connect(
            lambda checked: self.set_stacked_page(2, checked)
        )
        self.radioMagSnExp.clicked.connect(
            lambda checked: self.set_stacked_page(3, checked)
        )

        # Set default values
        self.last_input_data = ETCData()
        self.last_aux_data = AuxData()

    def exit_button_clicked(self):
        """The user clicked the "Exit" button

        Display a confirmation dialog, and quit if "Yes"
        """
        button = QtWidgets.QMessageBox.question(
            self, "", "Are you sure you want to quit?"
        )

        if button == QtWidgets.QMessageBox.StandardButton.Yes:
            QtWidgets.QApplication.quit()

    def compute_button_clicked(self):
        """The user clicked the "Compute" button

        Poll the ETC mode, and read in the mode-specific values along with all
        of the current Telescope Setup values to an :class:`ETCData` object.

        Run the mode-specific function to compute the needed pieces, and fill
        in the "Output" text labels in the GUI.
        """
        # Create an ETCData object, and load in the Telescope Setup information
        input_data = ETCData(
            band=self.inputBandpass.currentText(),
            binning=self.inputBinning.value(),
            seeing=self.inputSeeing.value(),
            phase=self.inputLunarphase.currentIndex(),
            airmass=self.inputAirmass.value(),
        )
        # Load in the proper mode-specific information and compute
        match self.stackedWidget.currentIndex():
            case 0:
                input_data.snr = float(self.input_1_Snr.text())
                input_data.mag = float(self.input_1_Magnitude.text())
                input_data.exptime = exptime_given_snr_mag(input_data)

            case 1:
                input_data.peak = float(self.input_2_PkCts.text())
                input_data.mag = float(self.input_2_Magnitude.text())
                input_data.exptime = exptime_given_peak_mag(input_data)

            case 2:
                input_data.exptime = float(self.input_3_Exptime.text())
                input_data.mag = float(self.input_3_Magnitude.text())
                input_data.snr = snr_given_exptime_mag(input_data)

            case 3:
                input_data.snr = float(self.input_4_Snr.text())
                input_data.exptime = float(self.input_4_Exptime.text())
                input_data.mag = mag_given_snr_exptime(input_data)

            case _:
                raise ValueError(
                    f"Incorrect index provided: {self.stackedWidget.currentIndex()}"
                )

        # Compute the auxillary data
        aux_data = self.compute_aux_data(input_data)

        # Display the results in the GUI
        self.outputBandpass.setText(input_data.band)
        self.outputMagnitude.setText(f"{input_data.mag:.1f}")
        self.outputSnratio.setText(f"{input_data.snr:.1f}")
        self.outputExptime.setText(f"{input_data.exptime:.1f}")
        self.outputBinning.setText(f"{input_data.binning:.0f}")

        self.outputNpix.setText(f"{aux_data.n_pix:.1f}")
        self.outputSkybright.setText(f"{aux_data.sky_bright:.1f}")
        self.outputNestar.setText(f"{aux_data.num_e_star:.1f}")
        self.outputPeakestar.setText(f"{aux_data.peak_e_star:.1f}")
        self.outputNoisestar.setText(f"{aux_data.noise_star:.1f}")
        self.outputNoisesky.setText(f"{aux_data.noise_sky:.1f}")
        self.outputNoiseccd.setText(f"{aux_data.noise_ccd:.1f}")

        # Put all results into instance attributes for possible saving
        self.last_input_data = input_data
        self.last_aux_data = aux_data

    def compute_aux_data(self, inuput_data: ETCData) -> AuxData:
        """compute_aux_data _summary_

        _extended_summary_

        Parameters
        ----------
        inuput_data : ETCData
            _description_

        Returns
        -------
        AuxData
            _description_
        """
        return AuxData()

    def set_stacked_page(self, index: int, checked: bool):
        """Show the proper ETC Mode page

        Parameters
        ----------
        index : :obj:`int`
            The stackedWidget index to show
        checked : :obj:`bool`
            Is the event a "checked" or "unchecked"?
        """
        # Only change the page if the radio button was checked (not unchecked)
        if checked:
            self.stackedWidget.setCurrentIndex(index)


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
        app = QtWidgets.QApplication(sys.argv)

        # Create a Qt widget, which will be our window.
        window = ETCWindow()
        window.show()  # IMPORTANT!!!!! Windows are hidden by default.

        # Start the event loop.
        app.exec()
