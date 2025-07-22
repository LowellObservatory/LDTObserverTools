# -*- coding: utf-8 -*-
#
#  This file is part of LDTObserverTools.
#
#   This Source Code Form is subject to the terms of the Mozilla Public
#   License, v. 2.0. If a copy of the MPL was not distributed with this
#   file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
#  Created on 18-Nov-2024
#
#  @author: tbowers
# pylint: disable=missing-function-docstring

"""DeVeny Collimator Focus Calculator TEST Module

NOTE: Use unit testing as a guide for rewriting the dfocus module in a more
      pythonic, less IDL style.  In particular, think about breaking down the
      functions into functionality rather than narrative steps.
"""

import pathlib

import astropy.nddata
import ccdproc
import numpy as np
import pytest

from obstools import dfocus
from obstools.tests.tstutlis import TEST_FILES
from obstools.utils import ObstoolsError

DV_TEST_FILES = TEST_FILES / "deveny"


def test_dfocus():
    # Not rreally sure how to test this one
    pass


def test_parse_focus_log():

    # Test parsing "last log", output types, and contents of output
    retval = dfocus.parse_focus_log(DV_TEST_FILES / "focus", "last")
    assert isinstance(retval, tuple)
    focus_cl, focus_id = retval
    assert isinstance(focus_cl, ccdproc.ImageFileCollection)
    assert isinstance(focus_id, str)
    assert focus_id == "20250203.081131"
    assert len(focus_cl.files) == 9

    # Test parsing of specifically named log
    focus_cl, focus_id = dfocus.parse_focus_log(
        DV_TEST_FILES / "focus", "deveny_focus.20250203.081131"
    )
    assert isinstance(focus_cl, ccdproc.ImageFileCollection)

    # Test when there is no log file using 'last'
    with pytest.raises(ObstoolsError) as err:
        retval = dfocus.parse_focus_log(pathlib.Path(".").resolve(), "last")
    assert str(err.value) == "No successful focus run completed in this directory"

    # Test when there is no log file using a specific log file
    with pytest.raises(ObstoolsError) as err:
        retval = dfocus.parse_focus_log(
            pathlib.Path(".").resolve(), "deveny_focus.20250203.081131"
        )
    assert str(err.value) == "Specified focus run not in this directory"


def test_parse_focus_headers():

    # Get the testing ImageFileCollection  -- This works if the above test does not fail
    focus_cl, _ = dfocus.parse_focus_log(
        DV_TEST_FILES / "focus", "deveny_focus.20250203.081131"
    )

    # Test parsing the focus headers
    focus_pars = dfocus.parse_focus_headers(focus_cl)
    assert isinstance(focus_pars, dfocus.FocusParams)

    # Check outputs versus expected
    assert focus_pars.mid_file == DV_TEST_FILES / "20250203.0025.fits"
    assert np.isclose(focus_pars.nominal, 2.76, atol=0.01)
    assert focus_pars.start == 9.0
    assert focus_pars.end == 13.0
    assert focus_pars.delta == 0.5
    assert focus_pars.mnttemp == 7.45
    assert focus_pars.binning == "1x1"

    # Use a focus log with all the same input DeVeny frame
    focus_cl, _ = dfocus.parse_focus_log(
        DV_TEST_FILES / "focus", "deveny_focus.20250202.allsame"
    )
    with pytest.raises(ObstoolsError) as err:
        focus_pars = dfocus.parse_focus_headers(focus_cl)
    assert str(err.value) == "No change in focus over this set of images"


@pytest.mark.filterwarnings("ignore::astropy.wcs.FITSFixedWarning")
def test_centered_trace():

    # Load in a CCD image
    spec2d = astropy.nddata.CCDData.read(DV_TEST_FILES / "20250203.0025.fits").data
    # Test building the trace from the data array
    trace = dfocus.centered_trace(spec2d)
    assert isinstance(trace, np.ndarray)

    # Check that the trace is one element "tall" across the array
    _, nx = spec2d.shape
    assert trace.shape == (1, nx)


@pytest.mark.filterwarnings("ignore::astropy.wcs.FITSFixedWarning")
def test_extract_spectrum():

    # Load in a CCD image, and build the trace -- tested above
    spec2d = astropy.nddata.CCDData.read(DV_TEST_FILES / "20250203.0025.fits").data
    trace = dfocus.centered_trace(spec2d)

    # Test extracting the 1D spectrum
    spec1d = dfocus.extract_spectrum(spec2d, trace, window=11, thresh=100.0)
    assert isinstance(spec1d, np.ndarray)

    # Check that the extracted spectrum is, indeed, one-dimensional
    _, nx = spec2d.shape
    assert spec1d.shape == (nx,)

    # Check raised errors
    traces = np.concatenate([trace, trace, trace, trace])
    with pytest.raises(ObstoolsError) as err:
        spec1d = dfocus.extract_spectrum(spec2d, traces, window=11)
    assert str(err.value) == "Cannot deal with multiple traces"


@pytest.mark.filterwarnings("ignore::astropy.wcs.FITSFixedWarning")
def test_find_lines():

    # Get the 1d spectrum from the trimmed image -- tested above
    ccd = astropy.nddata.CCDData.read(DV_TEST_FILES / "20250203.0025.fits")
    spec2d = ccdproc.trim_image(ccd, ccd.header["TRIMSEC"]).data
    trace = dfocus.centered_trace(spec2d)
    spec1d = dfocus.extract_spectrum(spec2d, trace, window=11, thresh=100.0)

    # Test finding the lines and check the return values
    retval = dfocus.find_lines(spec1d, 100.0)
    assert isinstance(retval, tuple)
    centers, fwhm = retval
    assert isinstance(centers, np.ndarray)
    assert isinstance(fwhm, np.ndarray)
    assert centers.size == fwhm.size

    # Check that the "correct" number of lines were found for this spectrum
    assert len(centers) == 42

    # Check that using a larger threshold produces fewer found lines
    centers, _ = dfocus.find_lines(spec1d, 1000.0)
    assert len(centers) == 29


@pytest.mark.filterwarnings("ignore::astropy.wcs.FITSFixedWarning")
def test_get_lines_from_ccd():

    # Test the wrapper function on the test frame; check return value
    ccd = astropy.nddata.CCDData.read(DV_TEST_FILES / "20250203.0025.fits")
    lines = dfocus.get_lines_from_ccd(ccd, 100.0)
    assert isinstance(lines, dfocus.LineInfo)
    # Test types of all ``LineInfo`` components
    assert isinstance(lines.spec_1d, np.ndarray)
    assert isinstance(lines.trace, np.ndarray)
    assert isinstance(lines.centers, np.ndarray)
    assert isinstance(lines.fwhm, np.ndarray)


@pytest.mark.filterwarnings("ignore::astropy.wcs.FITSFixedWarning")
def test_fit_focus_curves():

    # Get the focus parameters and middle image
    focus_cl, _ = dfocus.parse_focus_log(
        DV_TEST_FILES / "focus", "deveny_focus.20250203.081131"
    )
    focus_pars = dfocus.parse_focus_headers(focus_cl)
    mid_ccd = astropy.nddata.CCDData.read(focus_pars.mid_file)
    mid_lines = dfocus.get_lines_from_ccd(mid_ccd, 100.0)

    # Run the loop over all files to build the linewidth array
    line_width_array = []
    for ccd in focus_cl.ccds():
        this_lines = dfocus.get_lines_from_ccd(
            ccd, 100.0, trace=mid_lines.trace, verbose=False
        )
        line_dx = -4.0 * (ccd.header["COLLFOC"] - mid_ccd.header["COLLFOC"])
        line_widths = []
        for cen in mid_lines.centers:
            idx = np.abs((cen + line_dx) - this_lines.centers) < 3.0
            width = this_lines.fwhm[idx][0] if np.sum(idx) else np.nan
            line_widths.append(width if width > 2.0 else np.nan)
        line_width_array.append(np.array(line_widths))
    line_width_array = np.array(line_width_array)

    # Test the fitting of focus curves
    focus_curves = dfocus.fit_focus_curves(line_width_array, focus_pars)
    assert isinstance(focus_curves, dfocus.FocusCurves)
    assert isinstance(focus_curves.min_focus_values, np.ndarray)
    assert isinstance(focus_curves.optimal_focus_values, np.ndarray)
    assert isinstance(focus_curves.min_linewidths, np.ndarray)
    assert isinstance(focus_curves.fit_pars, np.ndarray)

    # We should be finding 42 lines fit with quadratic functions
    n_lines = 42
    n_fitpars = 3
    assert focus_curves.min_focus_values.shape == (n_lines,)
    assert focus_curves.optimal_focus_values.shape == (n_lines,)
    assert focus_curves.min_linewidths.shape == (n_lines,)
    assert focus_curves.fit_pars.shape == (n_lines, n_fitpars)
    # Check the values of the outputs
    assert np.isclose(np.nanmedian(focus_curves.min_focus_values), 9.477, atol=0.001)
    assert np.isclose(
        np.nanmedian(focus_curves.optimal_focus_values), 10.395, atol=0.001
    )
    assert np.isclose(np.nanmedian(focus_curves.min_linewidths), 2.38, atol=0.01)


def test_plot_lines():
    # Only contains plotting commands -- test at a later point, if need be
    pass


def test_plot_optimal_focus():
    # Only contains plotting commands -- test at a later point, if need be
    pass


def test_plot_focus_curves():
    # Only contains plotting commands -- test at a later point, if need be
    pass


def test_find_lines_in_spectrum():
    # Functionality is already tested -- maybe consider adding an integration test
    pass
