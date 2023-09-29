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

"""DeVeny Collimator Focus Calculator Module

LDTObserverTools contains python ports of various LDT Observer Tools

Lowell Discovery Telescope (Lowell Observatory: Flagstaff, AZ)
https://lowell.edu

This file contains the dfocus routine for computing the required collimator
focus for the DeVeny Spectrograph based on a focus sequence completed by the
DeVeny LOUI.

.. include:: ../include/links.rst
"""

# Built-In Libraries
import argparse
import os
import pathlib
import shutil
import sys
import warnings

# 3rd-Party Libraries
import astropy.io.fits
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np
from pypeit.scripts import scriptbase
import scipy.signal
from tqdm import tqdm

# Local Libraries
from obstools import deveny_grangle
from obstools import utils


# User-facing Function =======================================================#
def dfocus(
    path: pathlib.Path,
    flog: str = "last",
    thresh: float = 100.0,
    debug: bool = False,
    launch_preview: bool = True,
    leave_focus_files: bool = False,
):
    """Find the optimal DeVeny collimator focus value

    This is the user-facing ``dfocus`` function that calls all of the various
    following subroutines.  This function operates identically to the original
    IDL routine ``dfocus.pro``, but with the additional options of debugging and
    whether to launch Preview.app to show the PDF plots generated.

    Parameters
    ----------
    path : :obj:`~pathlib.Path`
        The path in which to find the ``deveny_focus.*`` files.
    flog : :obj:`str`, optional
        Focus log to process.  If unspecified, process the last sequence in
        the directory.  (Default: 'last')
        flog musy be of form: ``deveny_focus.YYYYMMDD.HHMMSS``
    thresh : :obj:`float`, optional
        Line intensity threshold above background for detection (Default: 100.0)
    debug : :obj:`bool`, optional
        Print debug statements  (Default: False)
    launch_preview : :obj:`bool`, optional
        Display the plots by launching Preview  (Default: True)
    leave_focus_files : :obj:`bool`, optional
        Do NOT move the focus frames to the ``focus/`` directory  (Default: False)
    """
    # Make a pretty title for the output of the routine
    n_cols = (os.get_terminal_size()).columns
    print("=" * n_cols)
    print("  DeVeny Collimator Focus Calculator")

    # Initialize a dictionary to hold lots of variables
    focus = initialize_focus_values(path, flog)
    if focus["delta"] == 0:
        print("\n** No successful focus run completed in this directory. **\n")
        sys.exit(1)

    # Process the middle image to get line centers, arrays, trace
    centers, trace, mid_collfoc, mspectra = process_middle_image(
        focus, thresh, debug=debug
    )

    # Run through files, showing a progress bar
    print("\n Processing arc images...")
    prog_bar = tqdm(total=focus["n"], unit="frame", unit_scale=False, colour="yellow")

    line_width_array = []
    for foc_file in focus["files"]:
        # Trim and extract the spectrum
        spectrum, collfoc = trim_deveny_image(foc_file)
        spectra = extract_spectrum(spectrum, trace, win=11)

        # Find FWHM of lines:
        _, these_centers, fwhm = find_lines(spectra, thresh=thresh, verbose=False)

        # Empirical shifts in line location due to off-axis paraboloid
        #  collimator mirror
        line_dx = -4.0 * (collfoc - mid_collfoc)

        # Keep only the lines from `these_centers` that match the
        #  reference image
        line_widths = []
        for cen in centers:
            # Find line within 3 pix of the (adjusted) reference line
            idx = np.where(np.absolute((cen + line_dx) - these_centers) < 3.0)[0]
            # If there's something there wider than 2 piux, use it... else NaN
            width = fwhm[idx][0] if len(idx) else np.nan
            line_widths.append(width if width > 2.0 else np.nan)

        # Append these linewidths to the larger array for processing
        line_width_array.append(np.asarray(line_widths))
        prog_bar.update(1)

    # Close the progress bar, end of loop
    prog_bar.close()
    line_width_array = np.asarray(line_width_array)

    print(
        f"\n  Median value of all linewidths: {np.nanmedian(line_width_array):.2f} pix"
    )

    # Fit the focus curve:
    min_focus_index, optimal_focus_index, min_linewidths, fit_pars = fit_focus_curves(
        line_width_array, fnom=focus["nominal"]
    )

    # Convert the returned indices into actual COLLFOC values, find medians
    print("=" * n_cols)
    min_focus_values = min_focus_index * focus["delta"] + focus["start"]
    optimal_focus_values = optimal_focus_index * focus["delta"] + focus["start"]
    # med_min_focus = np.real(np.nanmedian(min_focus_values))
    med_opt_focus = np.real(np.nanmedian(optimal_focus_values))

    if focus["binning"] != "1x1":
        print(f"*** CCD is operating in binning {focus['binning']} (col x row)")
    print(f"*** Recommended (Median) Optimal Focus Position: {med_opt_focus:.2f} mm")
    print(f"*** Note: Current Mount Temperature is: {focus['mnttemp']:.1f}ºC")

    # =========================================================================#
    # Make the multipage PDF plot
    with PdfPages(pdf_fn := path / f"pyfocus.{focus['id']}.pdf") as pdf:
        #  The plot shown in the IDL0 window: Plot of the found lines
        find_lines(
            mspectra,
            thresh=thresh,
            do_plot=True,
            focus_dict=focus,
            pdf=pdf,
            verbose=False,
        )

        # The plot shown in the IDL2 window: Plot of best-fit fwid vs centers
        plot_optimal_focus(focus, centers, optimal_focus_values, med_opt_focus, pdf=pdf)

        # The plot shown in the IDL1 window: Focus curves for each identified line
        plot_focus_curves(
            centers,
            line_width_array,
            min_focus_values,
            optimal_focus_values,
            min_linewidths,
            fit_pars,
            focus["delta"],
            focus["start"],
            fnom=focus["nominal"],
            pdf=pdf,
        )

    # Print the location of the plots, and move focus frames (if desired)
    print(f"\n  Plots have been saved to: {pdf_fn.name}\n")
    if not leave_focus_files:
        for foc_file in focus["files"]:
            shutil.move(foc_file, path)

    # Try to open with Apple's Preview App... if can't, oh well.
    if launch_preview:
        try:
            os.system(f"/usr/bin/open -a Preview {pdf_fn}")
        except Exception as err:
            print(f"Cannot open Preview.app\n{err}")


# Helper Functions (Chronological) ===========================================#
def initialize_focus_values(path: pathlib.Path, flog: str):
    """Initialize a dictionary of focus values

    Create a dictionary of values (mainly from the header) that can be used by
    subsequent routines.

    Parameters
    ----------
    path : :obj:`~pathlib.Path`
        The path to the current working directory
    flog : :obj:`str`
        Identifier for the focus log to be processed

    Returns
    -------
    :obj:`dict`
        Dictionary of the various needed quantities
    """
    # Parse the log file to obtain file list
    n_files, files, focus_id = parse_focus_log(path, flog)

    # Escape hatch if no files found (e.g., run in the wrong directory)
    if not n_files:
        return {"delta": 0}

    # Pull the spectrograph setup from the first focus file:
    hdr0 = astropy.io.fits.getheader(files[0])
    slitasec = hdr0["SLITASEC"]
    grating = hdr0["GRATING"]
    grangle = hdr0["GRANGLE"]
    lampcal = hdr0["LAMPCAL"]
    mnttemp = hdr0["MNTTEMP"]
    binning = hdr0["CCDSUM"].replace(" ", "x")

    # Compute the nominal line width
    nominal_focus = 2.94 * slitasec * deveny_grangle.deveny_amag(grangle)

    # Pull the collimator focus values from the first and last files
    focus_0 = hdr0["COLLFOC"]
    focus_1 = (astropy.io.fits.getheader(files[-1]))["COLLFOC"]
    # Find the delta between focus values
    try:
        delta_focus = (focus_1 - focus_0) / (n_files - 1)
    except ZeroDivisionError:
        delta_focus = 0

    # Examine the middle image
    mid_file = files[n_files // 2]

    dfl_title = (
        f"{mid_file.name}   Grating: {grating}   GRANGLE: "
        + f"{grangle:.2f}   Lamps: {lampcal}"
    )

    opt_title = (
        f"Grating: {grating}    Slit width: {slitasec:.2f} arcsec"
        + f"    Binning: {binning}    Nominal line width: "
        + f"{nominal_focus:.2f} pixels"
    )

    return {
        "n": n_files,
        "files": files,
        "mid_file": mid_file,
        "id": focus_id,
        "nominal": nominal_focus,
        "start": focus_0,
        "end": focus_1,
        "delta": delta_focus,
        "plot_title": dfl_title,
        "opt_title": opt_title,
        "mnttemp": mnttemp,
        "binning": binning,
    }


def parse_focus_log(path: pathlib.Path, flog: str):
    """Parse the focus log file produced by the DeVeny LOUI

    The DeVeny focus log file consists of filename, collimator focus, and other
    relevant information::

        :  Image File Name  ColFoc    Grating  GrTilt  SltWth     Filter    LampCal  MntTmp
        20230613.0026.fits    7.50   600/4900   27.04    1.20  Clear (C)   Cd,Ar,Hg    9.10
        20230613.0027.fits    8.00   600/4900   27.04    1.20  Clear (C)   Cd,Ar,Hg    9.10
        20230613.0028.fits    8.50   600/4900   27.04    1.20  Clear (C)   Cd,Ar,Hg    9.10
        20230613.0029.fits    9.00   600/4900   27.04    1.20  Clear (C)   Cd,Ar,Hg    9.10
        20230613.0030.fits    9.50   600/4900   27.04    1.20  Clear (C)   Cd,Ar,Hg    9.10
        20230613.0031.fits   10.00   600/4900   27.04    1.20  Clear (C)   Cd,Ar,Hg    9.10
        20230613.0032.fits   10.50   600/4900   27.04    1.20  Clear (C)   Cd,Ar,Hg    9.10

    This function parses out the filenames of the focus images for this run,
    largely discarding the remaining information in the focus log file.

    Parameters
    ----------
    path : :obj:`~pathlib.Path`
        The path to the current working directory
    flog : :obj:`str`
        Identifier for the focus log to be processed

    Returns
    -------
    n_files : :obj:`int`
        Number of files for this focus run
    files: :obj:`list`
        List of :obj:`~pathlib.Path` files associated with this focus run
    focus_id: :obj:`str`
        The focus ID
    """
    # Just to be sure...
    path = path.resolve()

    # Get the correct flog
    if flog.lower() == "last":
        focfiles = sorted(path.glob("deveny_focus*"))
        try:
            flog = focfiles[-1]
        except IndexError:
            # In the event of no files, return empty things
            return 0, [], ""
    else:
        flog = path / flog

    files = []
    with open(flog, "r", encoding="utf8") as file_object:
        # Discard file header
        file_object.readline()
        # Read in the remainder of the file, grabbing just the filenames
        for line in file_object:
            files.append(path.parent / line.strip().split()[0])

    # Return the list of files, and the FocusID
    return len(files), files, flog.name[-15:]


def process_middle_image(focus, thresh, debug=False):
    """Process the middle focus image

    This finds the lines to be measured -- presumably the middle is closest to focus

    Parameters
    ----------
    focus : :obj:`dict`, optional
        Dictionary containing needed variables for plot
    thresh : :obj:`float`
        Line intensity threshold above background for detection
    debug : :obj:`bool`, optional
        Print debug statements  (Default: False)

    Returns
    -------
    centers, :obj:`~numpy.ndarray`
        Centers
    trace, :obj:`~numpy.ndarray`
        Trace
    mid_collfoc, :obj:`float`
        Collimator focus of the middle frame
    mspectra, :obj:`~numpy.ndarray`
        The spectrum from the middle frame (for later plotting)
    """
    print(f"\n Processing center focus image {focus['mid_file']}...")
    spectrum, mid_collfoc = trim_deveny_image(focus["mid_file"])

    # Build the trace for spectrum extraction
    n_y, n_x = spectrum.shape
    trace = np.full(n_x, n_y / 2, dtype=float).reshape((1, n_x))
    mspectra = extract_spectrum(spectrum, trace, win=11)
    if debug:
        print(f"Traces: {trace}")
        print(f"Middle Spectrum: {mspectra}")

    # Find the lines in the extracted spectrum
    n_c, centers, _ = find_lines(mspectra, thresh=thresh)
    if debug:
        print(f"Back in the main program, number of lines: {n_c}")
        print(f"Line Centers: {[f'{cent:.1f}' for cent in centers]}")

    return centers, trace, mid_collfoc, mspectra


def trim_deveny_image(filename):
    """Trim a DeVeny Image

    The IDL code from which this was ported contains a large amount of
    vistigial code from previous versions of the DeVeny camera, including
    instances where the CCD was read out using 2 amplifiers and required
    special treatment in order to balance the two sides of the output image.

    The code below consists of the lines of code that were actually running
    using the keywords passed from current version of ``dfocus.pro``, and the
    pieces of that code that are actually used.

    Specifically, this routine trims off the 50 prescan and 50 postscan pixels,
    as well as several rows off the top and bottom.  (Extracts rows ``[12:512]``)

    Parameters
    ----------
    filename : :obj:`~pathlib.Path`
        Filename of the spectrum to get and trim

    Returns
    -------
    :obj:`~numpy.ndarray`
        The trimmed CCD image
    """

    # Parameters for DeVeny (2015 Deep-Depletion Device):
    nxpix, prepix = 2048, 50

    # Read in the file
    with astropy.io.fits.open(filename) as hdul:
        image = hdul[0].data
        collfoc = hdul[0].header["COLLFOC"]

    # Trim the image (remove top and bottom rows) -- why this particular range?
    # Trim off the 50 prepixels and the 50 postpixels; RETURN
    return image[12:512, prepix : prepix + nxpix], collfoc


def extract_spectrum(spectrum, traces, win):
    """Object spectral extraction routine

    Extract spectra by averaging over the specified window

    Parameters
    ----------
    spectrum : :obj:`~numpy.ndarray`
        The trimmed spectral image
    traces : :obj:`~numpy.ndarray`
        The trace(s) along which to extract spectra
    win : :obj:`int`
        Window over which to average the spectrum

    Returns
    -------
    :obj:`~numpy.ndarray`
        2D or 3D array of spectra of individual orders
    """
    # Spec out the shape, and create an empty array to hold the output spectra
    norders, n_x = traces.shape
    spectra = np.empty((norders, n_x), dtype=float)
    speca = np.empty(n_x, dtype=float)

    # Set extraction window size
    half_window = int(win) // 2

    for order in range(norders):
        # Because of python indexing, we need to "+1" the upper limit in order
        #   to get the full wsize elements for the average
        trace = traces[order, :].astype(int)
        for i in range(n_x):
            speca[i] = np.average(
                spectrum[trace[i] - half_window : trace[i] + half_window + 1, i]
            )
        spectra[order, :] = speca.reshape((1, n_x))

    return spectra


def find_lines(
    image,
    thresh=20.0,
    minsep=11,
    verbose=True,
    do_plot=False,
    focus_dict=None,
    pdf=None,
):
    """Automatically find and centroid lines in a 1-row image

    Uses :func:`scipy.signal.find_peaks` for this task

    Parameters
    ----------
    image : :obj:`~numpy.ndarray`
        Extracted spectrum
    thresh : :obj:`float`, optional
        Threshold above which to indentify lines [Default: 20 DN above bkgd]
    minsep : :obj:`int`, optional
        Minimum line separation for identification [Default: 11 pixels]
    verbose : :obj:`bool`, optional
        Produce verbose output?  [Default: False]
    do_plot : :obj:`bool`, optional
        Create a plot on the provided axes?  [Default: False]
    focus_dict : :obj:`dict`, optional
        Dictionary containing needed variables for plot  [Default: None]

    Returns
    -------
    n_c : :obj:`int`
        Number of lines found and returned
    centers : :obj:`~numpy.ndarray`
        Line centers (pixel #)
    fwhm : :obj:`~numpy.ndarray`
        The computed FWHM for each peak
    """
    # Get size and flatten to 1D
    _, n_x = image.shape
    spec = np.ndarray.flatten(image)

    # Find background from median value of the image:
    bkgd = np.median(spec)
    if verbose:
        print(
            f"  Background level: {bkgd:.1f}"
            + f"   Detection threshold level: {bkgd+thresh:.1f}"
        )

    # Use scipy to find peaks & widths -- no more janky IDL-based junk
    centers, _ = scipy.signal.find_peaks(
        newspec := spec - bkgd, height=thresh, distance=minsep
    )
    fwhm = (scipy.signal.peak_widths(newspec, centers))[0]

    if verbose:
        print(f" Number of lines found: {len(centers)}")

    # Produce a plot for posterity, if directed
    if do_plot:
        # Set up the plot environment
        _, axis = plt.subplots()
        tsz = 8

        # Plot the spectrum, mark the peaks, and label them
        axis.plot(np.arange(len(spec)), newspec)
        axis.set_ylim(0, (yrange := 1.2 * max(newspec)))
        axis.plot(centers, newspec[centers.astype(int)] + 0.02 * yrange, "k*")
        for cen in centers:
            axis.text(
                cen,
                newspec[int(np.round(cen))] + 0.03 * yrange,
                f"{cen:.0f}",
                fontsize=tsz,
            )

        # Make pretty & Save
        axis.set_title(focus_dict["plot_title"], fontsize=tsz * 1.2)
        axis.set_xlabel("CCD Column", fontsize=tsz)
        axis.set_ylabel("I (DN)", fontsize=tsz)
        axis.set_xlim(0, n_x + 2)
        axis.tick_params("both", labelsize=tsz, direction="in", top=True, right=True)
        plt.tight_layout()
        if pdf is None:
            plt.show()
        else:
            pdf.savefig()
        plt.close()

    return len(centers), centers, fwhm


def fit_focus_curves(fwhm, fnom=2.7, norder=2, debug=False):
    """Fit line / star focus curves

    [extended_summary]

    Parameters
    ----------
    fwhm : :obj:`~numpy.ndarray`
        Array of FWHM for all lines as a function of COLLFOC
    fnom : :obj:`float`, optional
        Nominal FHWM of an in-focus line. (Default: 2.7)
    norder : :obj:`int`, optional
        Polynomial order of the focus fit (Default: 2 = Quadratic)
    debug : :obj:`bool`, optional
        Print debug statements  (Default: False)

    Returns
    -------
    min_cf_idx_value : :obj:`~numpy.ndarray`
        Best fit focus values
    optimal_cf_idx_value : :obj:`~numpy.ndarray`
        Best fit linewidths
    min_linewidth : :obj:`~numpy.ndarray`
        Minimum linewidths
    foc_fits : :obj:`~numpy.ndarray`
        Fit parameters (for plotting)
    """
    # Warning Filter -- Polyfit RankWarning, don't wanna hear about it
    warnings.simplefilter("ignore", np.RankWarning)

    # Create the various arrays / lists needed
    n_focus, n_centers = fwhm.shape
    min_linewidth = []
    min_cf_idx_value = []
    optimal_cf_idx_value = []
    foc_fits = []

    # Fitting arrays (these are indices for collimator focus)
    cf_idx_coarse = np.arange(n_focus, dtype=float)
    cf_idx_fine = np.arange(0, n_focus - 1 + 0.1, 0.1, dtype=float)

    # Loop through lines to find the best focus for each one
    for i in range(n_centers):
        # Data are the FWHM for this line at different COLLFOC
        fwhms_of_this_line = fwhm[:, i]

        # Find unphysically large or small FWHM (or NaN) -- set to np.nan
        bad_idx = np.where(
            np.logical_or(fwhms_of_this_line < 1.0, fwhms_of_this_line > 15.0)
        )
        fwhms_of_this_line[bad_idx] = np.nan
        fwhms_of_this_line[np.isnan(fwhms_of_this_line)] = np.nan

        # If more than 3 of the FHWM are bad for this line, skip and go on
        if len(bad_idx) > 3:
            # Add values to the lists for proper indexing
            for flst in [min_linewidth, min_cf_idx_value, optimal_cf_idx_value]:
                flst.append(None)
            continue

        # Do a polynomial fit (norder) to the FWHM vs COLLFOC index
        # fit = np.polyfit(cf_idx_coarse, fwhms_of_this_line, norder)
        fit = utils.good_poly(cf_idx_coarse, fwhms_of_this_line, norder, 2.0)
        foc_fits.append(fit)
        if debug:
            print(f"In fit_focus_curves(): fit = {fit}")

        # If good_poly() returns zeros, deal with it accordingly
        if all(value == 0 for value in fit):
            min_cf_idx_value.append(np.nan)
            min_linewidth.append(np.nan)
            optimal_cf_idx_value.append(np.nan)
            continue

        # Use the fine grid to evaluate the curve miniumum
        focus_curve = np.polyval(fit, cf_idx_fine)  # fitfine
        min_cf_idx_value.append(cf_idx_fine[np.argmin(focus_curve)])  # focus
        min_linewidth.append(np.min(focus_curve))  # minfoc

        # Compute the nominal focus position as the larger of the two points
        #  where the polymonial function crosses fnom
        coeffs = [fit[0], fit[1], fit[2] - fnom]
        if debug:
            print(f"Roots: {np.roots(coeffs)}")
        optimal_cf_idx_value.append(np.max(np.real(np.roots(coeffs))))

    # After looping, return the items as numpy arrays
    return (
        np.asarray(min_cf_idx_value),
        np.asarray(optimal_cf_idx_value),
        np.asarray(min_linewidth),
        np.asarray(foc_fits),
    )


# Plotting Routines ==========================================================#
def plot_optimal_focus(
    focus, centers, optimal_focus_values, med_opt_focus, debug=False, pdf=None
):
    """Make the Optimal Focus Plot (IDL2 Window)

    [extended_summary]

    Parameters
    ----------
    focus : :obj:`dict`
        Dictionary of the various focus-related quantities
    centers : :obj:`~numpy.ndarray`
        Array of the centers of each line
    optimal_focus_values : :obj:`~numpy.ndarray`
        Array of the optimal focus values for each line
    med_opt_focus : :obj:`float`
        Median optimal focus value
    debug : :obj:`bool`, optional
        Print debug statements  (Default: False)
    """
    if debug:
        print("=" * 20)
        print(centers.dtype, optimal_focus_values.dtype, type(med_opt_focus))
    _, axis = plt.subplots()
    tsz = 8
    axis.plot(centers, optimal_focus_values, ".")
    axis.set_xlim(0, 2050)
    axis.set_ylim(focus["start"] - focus["delta"], focus["end"] + focus["delta"])
    axis.set_title(
        "Optimal focus position vs. line position, median =  "
        + f"{med_opt_focus:.2f} mm  "
        + f"(Mount Temp: {focus['mnttemp']:.1f}$^\\circ$C)",
        fontsize=tsz * 1.2,
    )
    axis.hlines(
        med_opt_focus,
        0,
        1,
        transform=axis.get_yaxis_transform(),
        color="magenta",
        ls="--",
    )
    axis.set_xlabel(f"CCD Column\n{focus['opt_title']}", fontsize=tsz)
    axis.set_ylabel("Optimal Focus (mm)", fontsize=tsz)
    axis.grid(which="both", color="#c0c0c0", linestyle="-", linewidth=0.5)

    axis.tick_params("both", labelsize=tsz, direction="in", top=True, right=True)
    plt.tight_layout()
    if pdf is None:
        plt.show()
    else:
        pdf.savefig()
    plt.close()


def plot_focus_curves(
    centers,
    line_width_array,
    min_focus_values,
    optimal_focus_values,
    min_linewidths,
    fit_pars,
    delta_focus,
    focus_0,
    fnom=2.7,
    pdf=None,
):
    """Make the big plot of all the focus curves (IDL1 Window)

    [extended_summary]

    Parameters
    ----------
    centers : :obj:`~numpy.ndarray`
        List of line centers from find_lines()
    line_width_array : :obj:`~numpy.ndarray`
        Array of line widths from each COLLFOC setting for each line
    min_focus_values : :obj:`~numpy.ndarray`
        List of the minimum focus values found from the polynomial fit
    optimal_focus_values : :obj:`~numpy.ndarray`
        List of the optimal focus values found from the polynomial fit
    min_linewidths : :obj:`~numpy.ndarray`
        List of the minumum linewidths found from the fitting
    fit_pars : :obj:`~numpy.ndarray`
        Array of the polynomial fit parameters for each line
    df : :obj:`float`
        Spacing between COLLFOC settings
    focus_0 : :obj:`float`
        Lower end of the COLLFOC range
    fnom : :obj:`float`, optional
        Nominal (optimal) linewidth  (Default: 2.7)
    """
    # Warning Filter -- Matplotlib doesn't like going from masked --> NaN
    warnings.simplefilter("ignore", UserWarning)

    # Set up variables
    n_foc, n_c = line_width_array.shape
    focus_idx = np.arange(n_foc)
    focus_x = focus_idx * delta_focus + focus_0

    # Set the plotting array
    ncols = 6
    nrows = n_c // ncols + 1
    fig, axes = plt.subplots(ncols=ncols, nrows=nrows, figsize=(8.5, 11))
    tsz = 6  # type size

    for i, axis in enumerate(axes.flatten()):
        if i < n_c:
            # Plot the points and the polynomial fit
            axis.plot(focus_x, line_width_array[:, i], "kD", fillstyle="none")
            axis.plot(focus_x, np.polyval(fit_pars[i, :], focus_idx), "g-")

            # Plot vertical lines to indicate minimum and optimal focus
            axis.vlines(min_focus_values[i], 0, min_linewidths[i], color="r", ls="-")
            axis.vlines(optimal_focus_values[i], 0, fnom, color="b", ls="-")

            # Plot parameters to make pretty
            # ax.set_ylim(0, np.nanmax(line_width_array[:,i]))
            axis.set_ylim(0, 7.9)
            axis.set_xlim(np.min(focus_x) - delta_focus, np.max(focus_x) + delta_focus)
            axis.set_xlabel("Collimator Position (mm)", fontsize=tsz)
            axis.set_ylabel("FWHM (pix)", fontsize=tsz)
            axis.set_title(
                f"LC: {centers[i]:.0f}  Fnom: {fnom:.2f} pixels", fontsize=tsz
            )
            axis.tick_params(
                "both", labelsize=tsz, direction="in", top=True, right=True
            )
            axis.grid(which="both", color="#c0c0c0", linestyle="-", linewidth=0.5)
        else:
            # Clear any extra positions if there aren't enough lines
            fig.delaxes(axis)

    plt.tight_layout()
    if pdf is None:
        plt.show()
    else:
        pdf.savefig()
    plt.close()


# Extra Routines =============================================================#
def find_lines_in_spectrum(filename, thresh=100.0):
    """Find the line centers in a spectrum

    This function is not directly utilized in ``dfocus``, but rather is included
    as a wrapper for several functions that can be used by other programs.

    Given the filename of an arc-lamp spectrum, this function returns a list
    of the line centers found in the image.

    Parameters
    ----------
    filename : :obj:`str`
        Filename of the arc frame to find lines in
    thresh : :obj:`float`, optional
        Line intensity threshold above background for detection [Default: 100]

    Returns
    -------
    :obj:`~numpy.ndarray`
        List of line centers found in the image
    """
    # Get the trimmed image
    spectrum, _ = trim_deveny_image(filename)

    # Build the trace for spectrum extraction
    n_y, n_x = spectrum.shape
    traces = np.full(n_x, n_y / 2, dtype=float).reshape((1, n_x))
    spectra = extract_spectrum(spectrum, traces, win=11)

    # Find the lines!
    _, centers, _ = find_lines(spectra, thresh=thresh)

    return centers


# Command Line Script Infrastructure (borrowed from PypeIt) ==================#
class DFocus(scriptbase.ScriptBase):
    """Script class for scrub_deveny_pickup tool

    Script structure borrowed from :class:`pypeit.scripts.sciptbase.ScriptBase`.
    """

    @classmethod
    def name(cls):
        """
        Provide the name of the script.  By default, this is the name of the
        module.
        """
        return f"{cls.__module__.rsplit('.', maxsplit=1)[-1]}"

    @classmethod
    def get_parser(cls, width=None):
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
            description="DeVeny Collimator Focus Calculator", width=width
        )
        parser.add_argument(
            "--flog",
            action="store",
            type=str,
            help="focus log to use (default: last)",
            default="last",
        )
        parser.add_argument(
            "--thresh",
            action="store",
            type=float,
            help="threshold for line detection (default: 100)",
            default=100.0,
        )
        parser.add_argument(
            "--nodisplay",
            action="store_false",
            help="DO NOT launch Preview.app to display plots",
        )
        parser.add_argument(
            "--leave_files",
            action="store_true",
            help="DO NOT move the focus frames to focus/",
        )
        return parser

    @staticmethod
    def main(args):
        """Main Driver

        Simple function that calls the primary function.
        """
        # Giddy Up!
        dfocus(
            pathlib.Path(".").resolve(),
            flog=args.flog,
            thresh=args.thresh,
            launch_preview=args.nodisplay,
            leave_focus_files=args.leave_files,
        )
