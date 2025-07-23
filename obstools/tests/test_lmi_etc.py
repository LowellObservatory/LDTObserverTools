# -*- coding: utf-8 -*-
#
#  This file is part of LDTObserverTools.
#
#   This Source Code Form is subject to the terms of the Mozilla Public
#   License, v. 2.0. If a copy of the MPL was not distributed with this
#   file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
#  Created on 06-Nov-2024
#
#  @author: tbowers
# pylint: disable=missing-function-docstring

"""LMI Exposure Time Calculator TEST Module
"""

import numpy as np
import pytest

from obstools import lmi_etc


def test_ETCData():
    # Test building `input_data` for the 4 GUI mode; calculated value is init'd at zero
    input_data = lmi_etc.ETCData(
        band="V", binning=2, seeing=1.2, airmass=1.65, phase=0, mag=18, snr=100
    )
    assert isinstance(input_data.exptime, float)
    assert input_data.exptime == 0.0
    input_data = lmi_etc.ETCData(
        band="B", binning=1, seeing=1.1, airmass=1.0, phase=7, peak=15000, mag=16.2
    )
    assert isinstance(input_data.exptime, float)
    assert input_data.exptime == 0.0
    input_data = lmi_etc.ETCData(
        band="u'", binning=3, seeing=0.7, airmass=1.12, phase=14, exptime=0.5, mag=10.2
    )
    assert isinstance(input_data.snr, float)
    assert input_data.snr == 0.0
    input_data = lmi_etc.ETCData(
        band="Hα-on", binning=4, seeing=2.5, airmass=2.3, phase=1, snr=10, exptime=600
    )
    assert isinstance(input_data.mag, float)
    assert input_data.mag == 0.0


def test_exptime_given_snr_mag():
    input_data = lmi_etc.ETCData(
        band="V", binning=2, seeing=1.2, airmass=1.65, phase=0, mag=18, snr=100
    )
    exptime = lmi_etc.exptime_given_snr_mag(input_data)
    assert isinstance(exptime, float)
    assert np.isclose(exptime, 3.459, atol=0.001)


def test_exptime_given_peak_mag():
    input_data = lmi_etc.ETCData(
        band="B", binning=1, seeing=1.1, airmass=1.0, phase=7, peak=15000, mag=16.2
    )
    exptime = lmi_etc.exptime_given_peak_mag(input_data)
    assert isinstance(exptime, float)
    assert np.isclose(exptime, 55.396, atol=0.001)


def test_snr_given_exptime_mag():
    input_data = lmi_etc.ETCData(
        band="u'", binning=3, seeing=0.7, airmass=1.12, phase=14, exptime=0.5, mag=10.2
    )
    snr = lmi_etc.snr_given_exptime_mag(input_data)
    assert isinstance(snr, float)
    assert np.isclose(snr, 21.627, atol=0.001)


def test_mag_given_snr_exptime():
    input_data = lmi_etc.ETCData(
        band="Hα-on", binning=4, seeing=2.5, airmass=2.3, phase=1, snr=10, exptime=600
    )
    mag = lmi_etc.mag_given_snr_exptime(input_data)
    assert isinstance(mag, float)
    assert np.isclose(mag, 21.424, atol=0.001)


def test_check_etc_inputs():
    # Input a reasonable set of values
    input_data = lmi_etc.ETCData(
        exptime=10.0,
        band="VR",
        mag=21.4,
        snr=75,
        binning=2,
        seeing=1.2,
        airmass=1.15,
        phase=2,
    )
    lmi_etc.check_etc_inputs(input_data)

    # Try bad values for each input
    input_data.airmass = 0.5
    with pytest.raises(ValueError):
        lmi_etc.check_etc_inputs(input_data)
    input_data.airmass = 1.15
    input_data.phase = 21
    with pytest.raises(ValueError):
        lmi_etc.check_etc_inputs(input_data)
    input_data.phase = 2
    input_data.seeing = 0.1
    with pytest.raises(ValueError):
        lmi_etc.check_etc_inputs(input_data)
    input_data.seeing = 1.2
    input_data.binning = 10
    with pytest.raises(ValueError):
        lmi_etc.check_etc_inputs(input_data)
    input_data.binning = 2.5
    with pytest.raises(ValueError):
        lmi_etc.check_etc_inputs(input_data)
    input_data.binning = 2
    input_data.exptime = 1.0e-5
    with pytest.raises(ValueError):
        lmi_etc.check_etc_inputs(input_data)
    input_data.exptime = 2000
    with pytest.raises(ValueError):
        lmi_etc.check_etc_inputs(input_data)
    input_data.exptime = 10.0
    input_data.mag = -26.74  # Looking at the Sun
    with pytest.raises(ValueError):
        lmi_etc.check_etc_inputs(input_data)
    input_data.mag = 30.0
    with pytest.raises(ValueError):
        lmi_etc.check_etc_inputs(input_data)
    input_data.mag = 21.4
    input_data.snr = 0.05
    with pytest.raises(ValueError):
        lmi_etc.check_etc_inputs(input_data)


def test_get_band_values():
    # Make sure good inputs yield good outputs
    band_info = lmi_etc.get_band_values("V")
    assert isinstance(band_info, lmi_etc.BandData)
    assert band_info.extinction == 0.14
    assert band_info.star20 == 670

    # Try a bad one
    with pytest.raises(ValueError):
        band_info = lmi_etc.get_band_values("F502N")  # No HST/WFC3 filters


def test_peak_counts():
    input_data = lmi_etc.ETCData(
        exptime=10.0, band="VR", mag=21.4, binning=2, seeing=1.2, airmass=1.15, phase=2
    )
    # This includes the sky and detector background
    counts = lmi_etc.peak_counts(input_data)
    assert isinstance(counts, float)
    assert np.isclose(counts, 3282, atol=1)
    # This is the star only
    counts = lmi_etc.peak_counts(input_data, star_only=True)
    assert np.isclose(counts, 85, atol=1)


def test_pixels_in_aperture():
    input_data = lmi_etc.ETCData(
        exptime=10.0, band="VR", mag=21.4, binning=2, seeing=1.2, airmass=1.15, phase=2
    )
    pix = lmi_etc.pixels_in_aperture(input_data)
    assert isinstance(pix, float)
    assert np.isclose(pix, 35.00, atol=0.01)
    # Decrease the seeing
    small_pix = 11.91
    input_data.seeing = 0.7
    pix = lmi_etc.pixels_in_aperture(input_data)
    assert np.isclose(pix, small_pix, atol=0.01)
    # Change the filter (should have same pixels)
    input_data.band = "u'"
    pix = lmi_etc.pixels_in_aperture(input_data)
    assert np.isclose(pix, small_pix, atol=0.01)
    # Halve the binning, 4x the pixels
    input_data.binning /= 2
    pix = lmi_etc.pixels_in_aperture(input_data)
    assert np.isclose(pix, small_pix * 4, atol=0.01)
    # Test minimum pix per aperture
    input_data.binning = 2
    input_data.seeing = 0.5
    pix = lmi_etc.pixels_in_aperture(input_data)
    assert np.isclose(pix, 9.00, atol=0.01)


def test_read_noise_in_aperture():
    input_data = lmi_etc.ETCData(
        exptime=10.0, band="VR", mag=21.4, binning=2, seeing=1.2, airmass=1.15, phase=2
    )
    readnoise = lmi_etc.read_noise_in_aperture(input_data)
    def_rn = 35.50  # Define the readnoise
    assert isinstance(readnoise, float)
    assert np.isclose(readnoise, def_rn, atol=0.01)
    # Readnoise is unaffected by band or airmass
    input_data.band = "yish"
    readnoise = lmi_etc.read_noise_in_aperture(input_data)
    assert np.isclose(readnoise, def_rn, atol=0.01)
    input_data.airmass = 2.34
    readnoise = lmi_etc.read_noise_in_aperture(input_data)
    assert np.isclose(readnoise, def_rn, atol=0.01)


def test_sky_counts_per_sec_in_aperture():
    input_data = lmi_etc.ETCData(
        exptime=10.0, band="VR", mag=21.4, binning=2, seeing=1.2, airmass=1.15, phase=2
    )
    sky_counts = lmi_etc.sky_counts_per_sec_in_aperture(input_data)
    assert isinstance(sky_counts, float)
    def_sc = 567.42  # Define the sky counts
    assert np.isclose(sky_counts, def_sc, atol=0.01)
    # More moon is more counts
    input_data.phase = 12
    sky_counts = lmi_etc.sky_counts_per_sec_in_aperture(input_data)
    assert sky_counts > def_sc
    # Smaller binning is less counts
    input_data.phase = 2
    input_data.binning = 1
    sky_counts = lmi_etc.sky_counts_per_sec_in_aperture(input_data)
    assert sky_counts < def_sc


def test_star_counts_per_sec():
    input_data = lmi_etc.ETCData(
        exptime=10.0, band="VR", mag=21.4, binning=2, seeing=1.2, airmass=1.15, phase=2
    )
    star_counts = lmi_etc.star_counts_per_sec(input_data)
    assert isinstance(star_counts, float)
    def_sc = 240.12  # Define the sky counts
    assert np.isclose(star_counts, def_sc, atol=0.01)
    # Brighter star, more counts
    input_data.mag -= 0.1
    star_counts = lmi_etc.star_counts_per_sec(input_data)
    assert star_counts > def_sc
    # Highser airmass, less counts
    input_data.mag += 0.1
    input_data.airmass += 0.1
    star_counts = lmi_etc.star_counts_per_sec(input_data)
    assert star_counts < def_sc


# Test GUI staticmethods
def test_ETCWindow_clean_etc_data():
    input_data = lmi_etc.ETCData(
        exptime=14.26135, airmass=2.3414, mag=21.43714, peak=7983.168
    )
    assert input_data.exptime == 14.26135
    # Clean the input_data
    input_data = lmi_etc.ETCWindow.clean_etc_data(input_data)
    assert input_data.exptime != 14.26135
    assert input_data.exptime == 14.26
    assert input_data.airmass != 2.3414
    assert input_data.airmass == 2.34
    assert input_data.mag != 21.43714
    assert input_data.mag == 21.44
    assert input_data.peak != 7983.168
    assert input_data.peak == 7983


def test_ETCWindow_compute_aux_data():
    input_data = lmi_etc.ETCData(
        exptime=10.0, band="VR", mag=21.4, binning=2, seeing=1.2, airmass=1.15, phase=2
    )
    aux_data = lmi_etc.ETCWindow.compute_aux_data(input_data)
    assert isinstance(aux_data, lmi_etc.AuxData)
    # Check all the outputs
    assert isinstance(aux_data.n_pix, float)
    assert isinstance(aux_data.sky_bright, float)
    assert isinstance(aux_data.num_e_star, float)
    assert isinstance(aux_data.peak_e_star, float)
    assert isinstance(aux_data.noise_star, float)
    assert isinstance(aux_data.noise_sky, float)
    assert isinstance(aux_data.noise_ccd, float)
    assert aux_data.n_pix == 35.00
    assert aux_data.sky_bright == 162.1
    assert aux_data.num_e_star == 2401.2
    assert aux_data.peak_e_star == 85
    assert aux_data.noise_star == 49.00
    assert aux_data.noise_sky == 75.33
    assert aux_data.noise_ccd == 35.50
