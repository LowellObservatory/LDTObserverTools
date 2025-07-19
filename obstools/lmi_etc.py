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
"""

# Built-In Libraries
import argparse
import dataclasses

# 3rd-Party Libraries
import astropy.table
import numpy as np
from PyQt6 import QtWidgets

# Local Libraries
from obstools.UI.ETCMainWindow import Ui_MainWindow
from obstools.UI.ETCDataWindow import Ui_ETCDataWindow
from obstools import utils

# Constants
SCALE = 0.12  # "/pix
READ_NOISE = 6.0  # e-
GAIN = 2.89  # e-/ADU
BIAS = 1050  # ADU (approx) for 2x2 binning


# Define the API as dataclasses, computation routines, and the external script
__all__ = [
    "ETCData",
    "AuxData",
    "exptime_given_snr_mag",
    "exptime_given_peak_mag",
    "snr_given_exptime_mag",
    "mag_given_snr_exptime",
    "LmiEtc",
]


# Dataclass objects ==========================================================#
@dataclasses.dataclass
class ETCData:
    """ETC Data Class

    Attributes
    ----------
    exptime : :obj:`float`
        User-defined exposure time (seconds)
    band : :obj:`str`
        The LMI filter for which to perform the calculation
    mag : :obj:`float`
        Magnitude in the band of the star desired
    snr : :obj:`float`
        Desired signal-to-noise ratio
    binning : :obj:`int`, optional
        Binning of the CCD
    seeing : :obj:`float`
        Size of the seeing disk (arcsec)
    airmass : :obj:`float`
        Airmass at which the observation will take place
    phase : :obj:`float`
        Moon phase (0-14)
    peak : :obj:`float`
        Desired peak count level on the CCD (e-)
    """

    exptime: float = None
    band: str = "V"
    mag: float = None
    snr: float = None
    binning: int = 2
    seeing: float = 1.0
    airmass: float = 1.0
    phase: float = 0
    peak: float = None


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


def peak_counts(input_data: ETCData, star_only: bool = False) -> float:
    """Compute the peak counts on the CCD for an exptime and mag

    Given a desired exposure time and stellar magnitude, compute the resulting
    counts on the CCD for a particular LMI Filter, moon phase, seeing and
    CCD binning.

    Parameters
    ----------
    input_data : :class:`ETCData`
        The input data needed for this calculation, placed in an
        :class:`ETCData` class.
    star_only : :obj:`bool`, optional
        Return the peak counts from the star only (not counting sky and detector bias)

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
    if star_only:
        return star_counts / (1.13 * fwhm**2) * input_data.exptime
    return (
        star_counts / (1.13 * fwhm**2) + sky_count_per_pixel_per_sec
    ) * input_data.exptime + BIAS * GAIN


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
    sky_count_per_pixel_per_sec = sky_count_per_arcsec2_per_sec * rscale**2
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


# GUI Classes ================================================================#
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

        # Connect the table view window
        self.tableWindow = TableWindow(self.table_colnames)

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
        self.buttonAdd2Table.clicked.connect(self.add_data_button_clicked)
        self.buttonShowTable.clicked.connect(self.show_table_button_clicked)

        # Set default values
        self.last_input_data = ETCData()
        self.last_aux_data = AuxData()
        self.etc_table = astropy.table.Table()

    def exit_button_clicked(self):
        """The user clicked the "Exit" button

        Display a confirmation dialog, and quit if "Yes"
        """
        button = QtWidgets.QMessageBox.question(
            self,
            "",
            "Are you sure you want to quit?",
            buttons=QtWidgets.QMessageBox.StandardButton.Ok
            | QtWidgets.QMessageBox.StandardButton.Cancel,
        )

        if button == QtWidgets.QMessageBox.StandardButton.Ok:
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
        self.outputExptime.setText(f"{input_data.exptime:.2f}")
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

    def show_table_button_clicked(self):
        """The user clicked the "Show Table" button

        The table of saved ETC values can be shown in a separate window.  By
        pressing the "Show Table" button, that window is displayed and the
        button text is changed to allow the user to hide the window upon next
        click of the button.
        """
        if self.tableWindow.isVisible():
            # Hide the window and make the button say "Show"
            self.tableWindow.hide()
            self.buttonShowTable.setText("Show Table")
        else:
            # Show the window and make the button say "Hide"
            self.tableWindow.show()
            self.buttonShowTable.setText("Hide Table")

    def add_data_button_clicked(self):
        """The user clicked the "Add to Table" button

        Add the current contents of :attr:`last_input_data` and
        :attr:`last_aux_data` to a :obj:`~astropy.table.Table` saved as an
        instance attribute.
        """
        save_dict = dataclasses.asdict(self.last_input_data)
        save_dict.update(dataclasses.asdict(self.last_aux_data))
        new_row = self.set_format(astropy.table.Table([save_dict]))
        self.etc_table = astropy.table.vstack([self.etc_table, new_row])
        # Send the new row to the TableWindow
        self.tableWindow.add_row_to_table(new_row)

    def compute_aux_data(self, input_data: ETCData) -> AuxData:
        """Compute auxillary output data

        These are the ancillary data computed from the ``input_data`` requested
        for the exposure time calculator.

        Parameters
        ----------
        input_data : :class:`ETCData`
            The input data needed for this calculation, placed in an
            :class:`ETCData` class.

        Returns
        -------
        :class:`AuxData`
            The auxillary data object associated with this ``input_data``
        """
        n_pix = pixels_in_aperture(input_data)
        total_star_cts = star_counts_per_sec(input_data) * input_data.exptime
        total_sky_cts = sky_counts_per_sec_in_aperture(input_data) * input_data.exptime
        return AuxData(
            n_pix=n_pix,
            sky_bright=total_sky_cts / n_pix,
            num_e_star=total_star_cts,
            peak_e_star=peak_counts(input_data, star_only=True),
            noise_star=np.sqrt(total_star_cts),
            noise_sky=np.sqrt(total_sky_cts),
            noise_ccd=read_noise_in_aperture(input_data),
        )

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

    @property
    def table_colnames(self) -> list:
        """Column names for the data table

        This method builds the column names from the dataclass attributes for
        both the :class:`ETCData` class containing the inputs to the
        calculation and the :class:`AuxData` class containing the detailed
        outputs.

        Returns
        -------
        :obj:`list`
            The column names for the output data table
        """
        return [f.name for f in dataclasses.fields(ETCData)] + [
            f.name for f in dataclasses.fields(AuxData)
        ]

    @staticmethod
    def set_format(table: astropy.table.Table) -> astropy.table.Table:
        """Set column formats for the data table

        This primarily deals with rounding and ease of display

        Parameters
        ----------
        table : :obj:`~astropy.table.Table`
            Input (unclean) table

        Returns
        -------
        :obj:`~astropy.table.Table`
            Output (cleaned) table
        """
        # TODO: Implement this!!!!!
        return table


class TableWindow(QtWidgets.QMainWindow, Ui_ETCDataWindow):
    """Exposure Time Calculator Saved Values Table Window Class

    The UI is defined in ETCWindow.ui and translated (via pyuic6) into python
    in ETCWindow.py.  This class inherits the UI and defines the various
    actions needed to compute ETCs from the GUI inputs.

    Parameters
    ----------
    table_colnames : :obj:`list`
        List of the column names for the data table
    """

    def __init__(self, table_colnames: list):
        self.table_colnames = table_colnames
        super().__init__()
        self.setupUi(self)

        # Connect buttons to actions
        self.buttonSavetable.clicked.connect(self.save_table_button_clicked)
        self.buttonCleartable.clicked.connect(self.clear_table_button_clicked)
        self.buttonRemoverow.clicked.connect(self.remove_row_button_clicked)

        # Looks prettier with this stuff
        self.update_table_cols()
        self.tableDatalog.resizeColumnsToContents()
        self.tableDatalog.resizeRowsToContents()

        # Actually show the table
        self.tableDatalog.show()

        # Init things
        self.last_row = 0

    def save_table_button_clicked(self):
        """The user clicked the "Save Table" button

        _extended_summary_
        """
        # Open the save file dialog, including format question
        print("Saving Table!!!")
        # TODO: Implement this!!!!!

    def clear_table_button_clicked(self):
        """The user clicked the "Clear Table" button

        _extended_summary_
        """
        # Open an "are you sure" dialog, then clear the table
        button = QtWidgets.QMessageBox.question(
            self,
            "",
            "Are you sure you want to clear the data table?",
            buttons=QtWidgets.QMessageBox.StandardButton.Ok
            | QtWidgets.QMessageBox.StandardButton.Cancel,
        )
        if button == QtWidgets.QMessageBox.StandardButton.Cancel:
            return

        # Clear the table, reset column/row N to zero, update table columns
        self.tableDatalog.clear()
        self.tableDatalog.setColumnCount(0)
        self.tableDatalog.setRowCount(0)
        self.update_table_cols()

        # Looks prettier with this stuff, and show the table
        self.tableDatalog.resizeColumnsToContents()
        self.tableDatalog.resizeRowsToContents()
        self.tableDatalog.show()

    def remove_row_button_clicked(self):
        """The user clicked the "Remove Row" button

        _extended_summary_
        """
        # Remove the indexed row
        bad = self.tableDatalog.currentRow()
        # -1 means we didn't select anything
        if bad != -1:
            # Clear the data we don't need anymore
            self.tableDatalog.removeRow(self.tableDatalog.currentRow())

            # Redraw
            self.tableDatalog.setVerticalHeaderLabels(
                f"{i}" for i in range(1, self.tableDatalog.rowCount() + 2)
            )
            self.tableDatalog.show()

    def update_table_cols(self):
        """Update the number and labels of table columns

        Update based on the current content of self.headers
        """
        # Add the number of columns we'll need for the header keys given
        for _ in self.table_colnames:
            col_position = self.tableDatalog.columnCount()
            self.tableDatalog.insertColumn(col_position)
        self.tableDatalog.setHorizontalHeaderLabels(self.table_colnames)

    def add_row_to_table(self, new_row: astropy.table.Table):
        """Add a row to the displayed data table

        Take the latest calculation from the ETC and add it to the TableWidget
        object holding the data.

        Parameters
        ----------
        new_row : :obj:`~astropy.table.Table`
            The new data row as an AstroPy table object.
        """
        # Disable fun stuff while we update
        self.tableDatalog.setSortingEnabled(False)
        self.tableDatalog.horizontalHeader().setSectionsMovable(False)
        self.tableDatalog.horizontalHeader().setDragEnabled(False)
        self.tableDatalog.horizontalHeader().setDragDropMode(
            QtWidgets.QAbstractItemView.DragDropMode.NoDragDrop
        )

        # Capture the last row position so we know where to start
        self.last_row = self.tableDatalog.rowCount()
        self.tableDatalog.insertRow(self.last_row)

        # Actually set the labels for rows
        self.tableDatalog.setVerticalHeaderLabels(
            [f"{i}" for i in range(1, self.last_row + 2)]
        )

        # Create the data table items and populate things
        #   Note! This is for use with headerDict style of grabbing stuff
        for n, row in enumerate(new_row):
            for m, hkey in enumerate(self.table_colnames):
                newitem = QtWidgets.QTableWidgetItem(str(row[hkey]))
                self.tableDatalog.setItem(n + self.last_row, m, newitem)

        # Resize to minimum required, then display
        # self.tableDatalog.resizeColumnsToContents()
        self.tableDatalog.resizeRowsToContents()

        # Seems to be more trouble than it's worth, so keep this commented
        # self.tableDatalog.setSortingEnabled(True)

        # Reenable fun stuff
        self.tableDatalog.horizontalHeader().setSectionsMovable(True)
        self.tableDatalog.horizontalHeader().setDragEnabled(True)
        self.tableDatalog.horizontalHeader().setDragDropMode(
            QtWidgets.QAbstractItemView.DragDropMode.InternalMove
        )

        # Looks prettier with this stuff
        self.tableDatalog.resizeColumnsToContents()
        self.tableDatalog.resizeRowsToContents()

        self.tableDatalog.show()

        # Should add this as a checkbox option to always scroll to bottom
        #   whenever a new file comes in...
        self.tableDatalog.scrollToBottom()


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

        Set up the top-level PyQt6 objects and start the event loop
        """
        # Create the QApplication object and main Qt window
        app = QtWidgets.QApplication([])
        _ = ETCWindow()

        # Giddy up!
        app.exec()
