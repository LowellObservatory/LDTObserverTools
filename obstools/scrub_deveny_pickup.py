# -*- coding: utf-8 -*-
#
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
#  Created on 05-Dec-2022
#
#  @author: tbowers

"""LDTObserverTools contains python ports of various LDT Observer Tools

Lowell Discovery Telescope (Lowell Observatory: Flagstaff, AZ)
http://www.lowell.edu

This module cleans the RF pickup noise from DeVeny spectral 2D images by
fitting a sinusoid to and then subtracting it from each row of the image.

The frequency of the sinusoid tends to vary from frame to frame but is roughly
CONSTANT within each frame.  Between the non-integer wavelength and the slight
delay between row read-out (caused by the advancement of the parallel register)
the phase shifts from row to row, yielding an ever-changing pattern in the
output images.  In order to identify the frequency in a given image, this
module fist taken the FFT of the flattened (`i.e.`, `time series`) image to
identify the strongest frequency.  (The function that performs the FFT also
modifies a copy of the data in various ways to remove spurious signals in the
FFT.)  From there, it uses a tightly bounded curve fitting (from
:obj:`scipy.optimize`) to identify the sinusoid in each line of the image.

Cosmic rays and other impulsive features in the 2D image can cause spurious
matches for the sinusoid.  So to guard against introducing artifacts into the
cleaned image from the fitting, the row-by-row fit coefficients are each
considered as a function and smoothed with an appropriate median window to
ensure the final subtracted pattern is a smooth function, matching the
underlying pickup noise.

The output of the module is a multi-extension FITS file containing [1] the
cleaned image, [2] the original image, [3] the pattern subtracted from the
original image, and [4] a FITS BinTable containing the sinusoid fit parameters
for each row of the image.
"""

# Built-In Libraries
import argparse
import datetime
import pathlib
import sys
import warnings

# 3rd Party Libraries
import astropy.io.fits
import astropy.nddata
import astropy.table
import astropy.units as u
import astropy.visualization
import astropy.wcs
import matplotlib.pyplot as plt
import numpy as np
from pypeit.spec2dobj import Spec2DObj
import scipy.fft
import scipy.ndimage
import scipy.optimize
import scipy.signal
from tqdm import tqdm

# Internal Imports
from obstools import utils

# Mean pixel readout time (7.25 µs)
PIX_DWELL = 7.25e-6

# Silence a couple known-but-don't-care-about warnings
warnings.simplefilter("ignore", astropy.wcs.FITSFixedWarning)
warnings.simplefilter("ignore", astropy.io.fits.verify.VerifyWarning)


# Narrative Functions ========================================================#
def iterative_pypeit_clean(filename: pathlib.Path):
    """Iteratively use PypeIt to clean pickup noise

    As part of the data-reduction process, PypeIt fits and extracts images as
    well as the sky background.  We can use the reduced ``spec2d`` files to
    get the "source-free" noise background, which _should_ make identifying the
    sinusoidal noise more straightforward.

    Parameters
    ----------
    filename : :obj:`pathlib.Path`
        The name of the file to clean
    """
    # Get the corresponding `spec2d` file for this raw file
    setup_dirs = sorted(filename.parent.glob("ldt_deveny_?"))
    print(setup_dirs)
    spec2d_file = [
        next(sd.joinpath("Science").glob(f"spec2d_{filename.stem}*"))
        for sd in setup_dirs
    ][0]

    # Load in the spec2d file
    spec2d = Spec2DObj.from_file(spec2d_file, "DET01")
    print(spec2d)
    interval = astropy.visualization.ZScaleInterval(n_samples=10000)

    img_gpm = spec2d.select_flag(invert=True)

    resid_thing = (spec2d.sciimg - spec2d.skymodel - spec2d.objmodel) * img_gpm.astype(
        float
    )

    resid_thing = np.rot90(resid_thing, k=-1)
    orig_thing = np.rot90(spec2d.sciimg, k=-1)

    # Get the expected pixel period from the FFT of the flattened array
    orig_pixperiod = flatten_clean_fft(orig_thing.copy(), "OBJECT")
    print(f"This is the FFT-predicted pixel period (orig): {orig_pixperiod:.1f} pix")
    resid_pixperiod = flatten_clean_fft(resid_thing.copy(), "BIAS", fft_plot=True)
    print(f"This is the FFT-predicted pixel period (resid): {resid_pixperiod:.1f} pix")

    # Fit the actual images
    resid_fitc = fit_lines(resid_thing.copy(), resid_pixperiod)
    orig_fitc = fit_lines(orig_thing.copy(), orig_pixperiod)

    for fit_coeffs in [orig_fitc, resid_fitc]:
        # Smooth the fit coefficients as a function of row to remove artifacts due
        #  to cosmic rays or strong sources.
        # NOTE: scipy.signal.medfilt() uses zero-padding for the ends of the signal,
        #       therefore subtract the mean and re-add it to eliminate edge weirdness
        fit_coeffs["smooth_a"] = smooth_coeffs(fit_coeffs["a"], kernel_size=9)
        fit_coeffs["smooth_lambda"] = smooth_coeffs(fit_coeffs["lambda"], kernel_size=9)
        # fit_coeffs["smooth_phi"] = fit_coeffs["phi"]
        fit_coeffs["smooth_phi"] = smooth_coeffs(fit_coeffs["phi"], ftype="savgol")
        # Print a happy statement of fact
        print(
            f"This is the mean fit pixel period: {fit_coeffs['smooth_lambda'].mean():.1f} pix"
        )

    cleaned_orig, pattern_orig = apply_pattern(orig_thing, orig_fitc)
    cleaned_resid, pattern_resid = apply_pattern(resid_thing, resid_fitc)

    # =====================================================#
    # Display the thing
    _, axes = plt.subplots(nrows=6, figsize=(9, 15))
    tsz = 8

    # Panel #1: Draw the original science image
    vmin, vmax = interval.get_limits(spec2d.sciimg)
    axes[0].imshow(np.rot90(spec2d.sciimg, k=-1), vmin=vmin, vmax=vmax, origin="lower")
    axes[0].axis("off")
    axes[0].set_ylabel("sciimg", fontsize=tsz)

    # Panel #2: Draw the sky model
    vmin, vmax = interval.get_limits(spec2d.skymodel)
    axes[1].imshow(
        np.rot90(spec2d.skymodel, k=-1), vmin=vmin, vmax=vmax, origin="lower"
    )
    axes[1].axis("off")
    axes[1].set_ylabel("skymodel", fontsize=tsz)

    # Panel #6: Draw the cleaned science image
    axes[5].imshow(
        np.rot90(spec2d.sciimg, k=-1) - pattern_resid + np.median(pattern_resid),
        vmin=vmin,
        vmax=vmax,
        origin="lower",
    )
    axes[5].axis("off")
    axes[5].set_ylabel("sciimg - pattern", fontsize=tsz)

    # Panel #3: Draw the resid image
    vmin, vmax = interval.get_limits(resid_thing)
    axes[2].imshow(resid_thing, vmin=vmin, vmax=vmax, origin="lower")
    axes[2].axis("off")
    axes[2].set_ylabel("resid", fontsize=tsz)

    # Panel #4: Draw the model
    # vmin, vmax = interval.get_limits(pattern_resid)
    axes[3].imshow(pattern_resid, vmin=vmin, vmax=vmax, origin="lower")
    axes[3].axis("off")
    axes[3].set_ylabel("pattern", fontsize=tsz)

    # Panel #5: Draw the cleaned image
    axes[4].imshow(
        resid_thing - pattern_resid + np.median(pattern_resid),
        vmin=vmin,
        vmax=vmax,
        origin="lower",
    )
    axes[4].axis("off")
    axes[4].set_ylabel("resid - pattern", fontsize=tsz)

    axes[0].set_title(filename.name, fontsize=tsz + 2)

    plt.tight_layout()
    plt.savefig("scrub_thing.pdf")
    plt.close()

    # Make disgnostic plot
    make_sinusoid_fit_plots(resid_fitc, filename, resid_pixperiod)


def clean_pickup(
    filename: pathlib.Path,
    use_hann: bool = False,
    sine_plot: bool = True,
    fft_plot: bool = False,
):
    """Clean the Pickup Noise

    This is the first version of the cleaning code, which operates directly on
    the raw DeVeny image.  It suffers from confusion from bright night-sky
    lines and brighter objects.  These deficiencies led to the development of
    the iterative scrubber.  This code is kept here for historical purposes.

    Parameters
    ----------
    filename : :obj:`pathlib.Path`
        The name of the file to clean
    use_hann : bool
        Use a Hann window when cleaning the image for FFT (Default: False)
    sine_plot : bool
        Create a plot of the sinusoid pattern fits for this image?
        (Default: True)
    fft_plot : bool
        Create a plot of the FFT analysis for this image?
        (Default: False)
    """
    # Read in the image and create copy arrays on which to work
    ccd = astropy.nddata.CCDData.read(filename)
    _, ncol = ccd.data.shape

    # Check if this images has already been post-processed by this routine
    if "postproc" in ccd.header and ccd.header["postproc"] == "deveny_pickup":
        print(
            f"\n   Skipping file {filename.name}; already processed with deveny_pickup"
        )
        return

    # Check image type; only process OBJECT and BIAS
    if (imgtyp := ccd.header["imagetyp"]) in ["OBJECT", "BIAS"]:
        print(f"\n   Processing {imgtyp} file {filename.name}")
    else:
        print(f"\n   Skipping {imgtyp} file {filename.name}")
        return

    # Get the expected pixel period from the FFT of the flattened array
    pixel_period = flatten_clean_fft(
        ccd.data.copy(), imgtyp, use_hann=use_hann, fft_plot=fft_plot
    )
    print(f"This is the FFT-predicted pixel period: {pixel_period:.1f} pix")

    # Compute the fit coefficients
    fit_coeffs = fit_lines(ccd.data, pixel_period)

    # Smooth the fit coefficients as a function of row to remove artifacts due
    #  to cosmic rays or strong sources.
    # NOTE: scipy.signal.medfilt() uses zero-padding for the ends of the signal,
    #       therefore subtract the mean and re-add it to eliminate edge weirdness
    fit_coeffs["smooth_a"] = smooth_coeffs(fit_coeffs["a"])
    fit_coeffs["smooth_lambda"] = smooth_coeffs(fit_coeffs["lambda"])
    fit_coeffs["smooth_phi"] = smooth_coeffs(fit_coeffs["phi"], ftype="savgol")
    # Print a happy statement of fact
    print(
        f"This is the mean fit pixel period: {fit_coeffs['smooth_lambda'].mean():.1f} pix"
    )

    cleaned_array, pattern_array = apply_pattern(ccd.data, fit_coeffs)

    # Make a plot, if desired
    if sine_plot:
        make_sinusoid_fit_plots(fit_coeffs, filename, pixel_period)

    # Package everything into a multiextension FITS file ===========#
    time_str = datetime.datetime.utcnow().isoformat(sep=" ", timespec="seconds")
    history_str = f"Written by package deveny_pickup: {time_str} UTC"
    # Primary: HDU #0
    primary_hdu = astropy.io.fits.PrimaryHDU(header=ccd.header)
    primary_hdu.update_header()
    primary_hdu.header.append(
        (
            "postproc",
            "deveny_pickup",
            "Post-processing algorithm or package applied",
        )
    )
    primary_hdu.header.append(
        (
            "post_ext",
            1,
            "Extension holding cleaned image (zero-indexed)",
        )
    )
    primary_hdu.header["HISTORY"] = history_str

    # For the image HDUs, include a basic header
    img_hdr = astropy.io.fits.Header(
        {
            "BUNIT": "ADU",
            "HISTORY": history_str,
        }
    )

    # Cleaned Image: HDU #1
    clean_hdu = astropy.io.fits.ImageHDU(
        astropy.nddata.CCDData(cleaned_array, unit=u.adu),
        name="CLEANED",
        header=img_hdr,
    )
    primary_hdu.header.append(("EXT0001", "CLEANED"))

    # Original Data: HDU #2
    orig_hdu = astropy.io.fits.ImageHDU(ccd, name="ORIGINAL", header=img_hdr)
    primary_hdu.header.append(("EXT0002", "ORIGINAL"))

    # Removed Pattern: HDU #3
    pattern_hdu = astropy.io.fits.ImageHDU(
        astropy.nddata.CCDData(pattern_array, unit=u.adu),
        name="PATTERN",
        header=img_hdr,
    )
    primary_hdu.header.append(("EXT0003", "PATTERN"))

    # Fit coefficients BinTable: HDU #4
    table_hdu = astropy.io.fits.BinTableHDU(fit_coeffs, name="FIT DATA")
    table_hdu.header.set(
        "fft_per",
        np.around(pixel_period, 1),
        "Pattern periodicity computed from FFT (pixels)",
        before="EXTNAME",
    )
    table_hdu.header.set(
        "mean_per",
        np.around(fit_coeffs["smooth_lambda"].mean(), 1),
        "Mean fit pattern periodicity (pixels)",
        before="EXTNAME",
    )
    table_hdu.header.set(
        "n_lam",
        np.around(ncol / fit_coeffs["smooth_lambda"].mean(), 2),
        "Number of (mean) wavelengths across the image",
        before="EXTNAME",
    )
    table_hdu.header.set(
        "mean_amp",
        np.float64(np.around(np.median(fit_coeffs["smooth_a"]) * 2.0, 2)),
        "Median pattern peak-to-peak amplitude (ADU)",
        before="EXTNAME",
    )
    for col in [c for c in fit_coeffs.colnames if "corr_" in c]:
        table_hdu.header.set(
            col,
            np.float64(np.around(np.median(fit_coeffs[col]), 4)),
            "Median correlation coeff between parameters",
            before="EXTNAME",
        )

    table_hdu.header[
        "COMMENT"
    ] = "Table contains the sinusoid fit coefficients for each row of the image"
    table_hdu.header["HISTORY"] = history_str
    primary_hdu.header.append(("EXT0004", "FIT DATA"))

    # Assemble into HDUList and write to disk
    hdul = astropy.io.fits.HDUList(
        [primary_hdu, clean_hdu, orig_hdu, pattern_hdu, table_hdu]
    )
    fn_parts = filename.name.split(".")
    out_fn = f"{fn_parts[0]}.C{fn_parts[1]}.{'.'.join(fn_parts[2:])}"
    hdul.writeto(out_fn, overwrite=True)


# Task-Oriented Helper Functions =============================================#
def flatten_clean_fft(
    data_array: np.ndarray, imgtyp: str, use_hann: bool = False, fft_plot: bool = False
) -> float:
    """Find the peak frequency of sinusoidal noise

    The function name describes what it does:
        * Flatten the 2D array into a 1D timeseries
        * Clean the ends of the CCD for smoother transition
        * Take the FFT to find the proper frequency to return

    In order to more effectively find the frequencies of persistent signal
    across the detector, start by flattening out the data array into a one-
    dimensional timeseries-like object representing the order in which the
    pixels are read out.

    According to information in Chapter 5 of `Scientific Charge-Coupled
    Devices` by James R. Janesick (2001), the time to advance a row along
    the parallel register is small compared to the individual pixel
    digitization time because the latter requires a capacitive stabilization
    time to minimize readnoise, whereas parallel register moves simply move
    the charge without trying to measure it.  Therefore, there is no spacing
    in between the rows of the flattened arrays to represent the time required
    for the parallel charge transfer.

    Once the array is flattened, we deal with row-to-row edge effects by
    modifying the flattened array to replace the 10 pixels at the end of one
    row and the 10 pixels at the start of the next row with a linearly varying
    series of values.  The first 10 and last 10 pixels of the flattened array
    are set to the next adjacent value.  This removes some of the "ringing"
    power at integer frequencies of an entire row and places it into the
    fundamental row-length frequency.  There are still sharp edges in the
    flattened array caused by night sky lines, but most of the power from these
    lines should be at the fundamental frequency (whole-row).

    The detected pickup noise has a period in the range 100-300 pixels, so
    peaks in the power spectrum at lower frequencies than this (`e.g.`, the
    fundamental row-length frequency) are masked when choosing the return
    value.

    Additionally, the assumption of continuous readout is close but not
    entriely accurate.  The slight pause in the readout caused by the parallel
    register readout causes the continuous readout amplifier pickup signal
    to manifest itself as a Gaussian packet of related frequencies in the FFT.
    By smoothing the power spectrum ( abs(FFT)^2 ) with a Gaussian kernel, narrow
    spikes in the FFT (caused by autocorrelation of night sky lines, cosmic
    rays, edge effects) are diluted compared to the sought signal.  The
    returned period is that of the maximum of the Gaussian-smoothed absolute
    value squared.

    Parameters
    ----------
    data_array : :obj:`~numpy.ndarray`
        _description_
    imgtyp : :obj:`str`
        The frame image type from the FITS header
    use_hann : :obj:`bool`
        Use a Hann window on each row in addition to subtracting the row mean
        (Default: False)
    fft_plot : :obj:`bool`
        Create a debugging plot of the FFT analysis?  (Default: False)

    Returns
    -------
    :obj:`float`
        The pixel period of the sinusoidal oscillation to be removed
    """

    # Compute the shape of the input array
    nrow, ncol = data_array.shape

    # Before flattening, mask out the main object, as fit from vertical profile
    if imgtyp == "OBJECT":
        vert_profile = np.median(data_array, axis=1)
        medval = np.median(vert_profile)
        # Fit a gaussian
        popt, _ = scipy.optimize.curve_fit(
            utils.gaussian_function,
            np.arange(nrow),
            vert_profile,
            p0=[
                np.max(vert_profile) - medval,  # a
                np.argmax(vert_profile),  # mu
                5.0,  # sig
                medval,  # y0
            ],
            bounds=[0, [2**16 - medval, vert_profile.size, 15, np.max(vert_profile)]],
        )
        # Build the mask at ±5σ, if the object is bright enough to matter
        if popt[0] > 20:
            mask = (popt[1] + 5.0 * np.array([-popt[2], popt[2]])).astype(int)
            data_array[mask[0] : mask[1], :] = medval

    # If so desired, apply a Hann window to each row
    if use_hann:
        hann = 0.5 * (1.0 - np.cos(2.0 * np.pi * np.arange(ncol) / ncol))
        data_array = data_array * hann

    # Flatten the array; force Column-Major to keep python from getting ideas
    flat_array = data_array.flatten(order="C").astype(np.float64)

    # Smooth out the beginning of the first and end of the last rows
    flat_array[:10] = flat_array[10]
    flat_array[-10:] = flat_array[-11]

    # For all other row breaks, index by `rownum``
    for rownum in range(1, nrow):
        # Build up a slice
        idx = np.s_[rownum * ncol - 10 : rownum * ncol + 10]
        # Pull the chunk and replace it with something smooth
        chunk = flat_array[idx]
        flat_array[idx] = np.linspace(chunk[0], chunk[-1], num=len(chunk))

    # Onward to the FFT!
    # Subtract the mean to remove power from 0 frequency
    y_fft = scipy.fft.fft(flat_array - np.mean(flat_array))

    # Compute the frequency array; N = N points, T = Spacing (PIX_DWELL)
    x_fft = scipy.fft.fftfreq(n_pts := flat_array.size, PIX_DWELL)

    # Mask out periods longer than ~500 pixels (this also removes negative freq's)
    y_fft[x_fft < pixper_tofrom_hz(500)] = 0.0

    # Compute the power spectrum, smoothed with a 10-Hz Gaussian
    pspec = scipy.ndimage.gaussian_filter1d(np.abs(y_fft) ** 2, 10)

    # Find the peak of the smoothed power spectrum, seeking Gaussin packets
    peak_freq = x_fft[np.argmax(pspec[: n_pts // 2])]

    if fft_plot:
        create_fft_plot(flat_array, y_fft, x_fft, peak_freq)

    # Return the period (in pixels) that corresponds to this frequency
    return pixper_tofrom_hz(peak_freq)


def fit_lines(data_array: np.ndarray, pixel_period: float) -> astropy.table.Table:
    """Fit a sinusoid to each line in the image

    This is like a mini-driver function that fits a sinusoid to each line in
    the image.

    Parameters
    ----------
    data_array : :obj:`~numpy.ndarray`
        The image array to be processed.  This should be HORIZONTAL and in the
        same orientation as images written out by lois directly from the
        instrument.
    pixel_period : :obj:`float`
        The predicted pixel period for this image array, as computed by
        :func:`flatten_clean_fft`.

    Returns
    -------
    :obj:`~astropy.table.Table`
        A table object containing the fit coefficients, one table row per row
        of the input ``data_array``.
    """
    # Array shape
    nrow, _ = data_array.shape
    img_mean = np.nanmean(data_array)

    # TODO: Do some sigma-clipping here to the whole image, maybe?

    # These are the starting guesses and bounds for the sinusoidal fit
    #  [a, lam, phi, y0, lin, quad]
    p0 = [4, pixel_period, 0.5, img_mean]  # + 4 * [0]
    bounds = (
        [0, pixel_period / 2, 0, img_mean - 100],  # + 4 * [-np.inf],
        [10, pixel_period * 2, 1, img_mean + 100],  # + 4 * [np.inf],
    )

    # Loop over the rows in the image to fit the pattern
    fit_coeffs = []
    progress_bar = tqdm(
        total=nrow,
        unit="row",
        unit_scale=False,
        colour="#87CEEB",
    )

    for img_row in range(nrow):
        # Pull this line, minus the last few pixels at each end
        line = data_array[img_row, 5:-5]

        # To mitigate cosmic rays and night sky lines, sigma clip @ 5σ
        sig_clip = np.std(line) * 5.0 + np.median(line)
        line[line > sig_clip] = sig_clip

        # Smooth the line with a median filter at 1/10th the pixel period
        line = smooth_coeffs(line, kernel_size=nearest_odd(pixel_period / 10.0))

        # Perform the curve fit
        try:
            xp = np.arange(line.size)
            popt, pcov = scipy.optimize.curve_fit(
                sinusoid, xp, line, p0=p0, bounds=bounds
            )
        except RuntimeError:
            # Reached max function evaluations; set popt and pcov
            print("Runtime error!")
            popt = [0, pixel_period, 0.5, 2400]
            pcov = np.diag(np.ones(len(popt)))

        diagnostic = False
        if diagnostic:
            # ======
            # Make a diagnostic plot
            _, axis = plt.subplots(figsize=(9, 3))
            axis.plot(xp, line, "k-", linewidth=0.75)
            axis.plot(xp, sinusoid(xp, *popt), "r-")
            plt.show()
            plt.close()

        # Compute standard deviations and correlation matrix
        pstd = np.sqrt(np.diag(pcov))
        # Check for zero stddev, replace with small stddev
        pstd[pstd == 0] = 1.0e-8
        inv_pstd = np.linalg.inv(np.diag(pstd))
        pcor = np.matmul(np.matmul(inv_pstd, pcov), inv_pstd)

        # Make sure the phase shift is a smoothly varying function, even if
        #  it goes outside the range 0 - 1.
        # NOTE: This explicitely assumes the phase shift from one line to the
        #       next is always less than 0.5 (π radians, 180º).
        if fit_coeffs:
            while popt[2] - fit_coeffs[-1]["phi"] > 0.5:
                popt[2] -= 1.0
            while popt[2] - fit_coeffs[-1]["phi"] < -0.5:
                popt[2] += 1.0

        # Add the fit coefficients to the list
        fit_coeffs.append(
            {
                "a": popt[0],
                "a_err": pstd[0],
                "lambda": popt[1],
                "lambda_err": pstd[1],
                "phi": popt[2],
                "phi_err": pstd[2],
                "y0": popt[3],
                "y0_err": pstd[3],
                "corr_a_lambda": pcor[0][1],
                "corr_a_phi": pcor[0][2],
                "corr_a_y0": pcor[0][3],
                "corr_lambda_phi": pcor[1][2],
                "corr_lambda_y0": pcor[1][3],
                "corr_phi_y0": pcor[2][3],
            }
        )
        progress_bar.update(1)

    progress_bar.close()
    # After the image has been fit, convert the list of dict into a Table
    return astropy.table.Table(fit_coeffs)


def smooth_coeffs(
    col: astropy.table.Column, kernel_size: int = 21, ftype: str = "median"
) -> np.ndarray:
    """Smooth out the coefficients from the fit coefficients table

    Attempt to allow for a smooth noise pattern from row to row.

     Parameters
     ----------
     col : :obj:`~astropy.table.Column`
         _description_
     kernel_size : :obj:`int`, optional
         Kernel or window size for the smoothing function (Default: 21)
     ftype : :obj:`str`, optional
         Smoothing function to use.  Options are "median" and "savgol".
         (Default: "median")

     Returns
     -------
     :obj:`~numpy.ndarray`
         The smoothed coefficients
    """
    if ftype == "median":
        return (
            scipy.signal.medfilt(col - (med := np.median(col)), kernel_size=kernel_size)
            + med
        )
    if ftype == "savgol":
        return scipy.signal.savgol_filter(
            col, window_length=kernel_size, polyorder=2, mode="nearest"
        )
    raise ValueError(f"Filter type {ftype} not supported by this function.")


def apply_pattern(
    input_array: np.ndarray, fit_coeffs: astropy.table.Table
) -> tuple[np.ndarray, np.ndarray]:
    """Construct the pattern from the fit coefficients

    _extended_summary_

    Parameters
    ----------
    input_array : :obj:`~numpy.ndarray`
        The input array to be cleaned
    fit_coeffs : `astropy.table.Table`
        _description_

    Returns
    -------
    cleaned_array : :obj:`~numpy.ndarray`
        The cleaned sience array
    pattern_array : :obj:`~numpy.ndarray`
        The pattern removed from ``input_array``
    """
    # Compute the image mean for the pattern array
    img_mean = np.nanmean(input_array)

    # Create the arrays that the processed data will go into
    cleaned_array = input_array.copy().astype(np.float64)
    pattern_array = input_array.copy().astype(np.float64)

    # Loop through the image and use the smoothed sinusoid fit coefficients
    for img_row, table_row in enumerate(fit_coeffs):
        # Apply the adjusted pattern to the entire row
        line = input_array[img_row, :]

        # Compute the pattern
        pattern = sinusoid(
            np.arange(line.size),
            table_row["smooth_a"],
            table_row["smooth_lambda"],
            table_row["smooth_phi"],
            0,
        )

        # Fill in the new arrays with the cleaned data and pattern
        cleaned_array[img_row, :] = line - pattern
        pattern_array[img_row, :] = pattern + img_mean

    # Return when done
    return cleaned_array, pattern_array


# Plotting Functions for QA / Debugging ======================================#
def create_fft_plot(
    flat_array: np.ndarray, y_fft: np.ndarray, x_fft: np.ndarray, peak_freq: float
):
    """Create the FFT analysis plot

    _extended_summary_

    Parameters
    ----------
    flat_array : :obj:`~numpy.ndarray`
        The flattened (1D) pixel array
    y_fft : :obj:`~numpy.ndarray`
        The (complex) FFT values
    x_fft : :obj:`~numpy.ndarray`
        The FFT frequencies for the values in ``y_fft``
    peak_freq : :obj:`float`
        The peak frequency measured from the FFT
    """
    # Compute the power spectrum as |FFT|^2; smoothed with a 10-Hz Gaussian
    pspec = scipy.ndimage.gaussian_filter1d(np.abs(y_fft) ** 2, 10)

    # Set up the plotting environment
    _, axes = plt.subplots(nrows=4, figsize=(6.4, 6.4))
    tsz = 8

    # Make a series of secondary axes
    secaxes = axes.copy()

    # Begin
    axes[0].plot(flat_array, linewidth=0.1, color="k")
    axes[0].set_xlabel("Linear pixel number", fontsize=tsz)

    axes[1].plot(x_fft, np.real(y_fft), linewidth=0.1)
    axes[1].set_ylabel("Re(FFT)", fontsize=tsz)

    axes[2].plot(x_fft, np.imag(y_fft), linewidth=0.1)
    axes[2].set_ylabel("Im(FFT)", fontsize=tsz)

    axes[3].plot(x_fft, pspec, linewidth=0.3)
    ylim = axes[3].get_ylim()
    axes[3].plot(x_fft, np.abs(y_fft) ** 2, zorder=0, linewidth=0.1)
    axes[3].set_ylabel("|FFT|^2", fontsize=tsz)
    axes[3].set_ylim(ylim)

    for i, axis in enumerate(axes[1:], 1):
        secaxes[i] = axis.secondary_xaxis(
            "bottom", functions=(pixper_tofrom_hz, pixper_tofrom_hz)
        )
        secaxes[i].set_xlabel("Sinusoid Period (Pixels)", fontsize=tsz)
        axis.set_xlim(pixper_tofrom_hz(np.array([500, 75])))
        axis.set_xticks([])
        axis.vlines(
            peak_freq,
            0,
            1,
            transform=axis.get_xaxis_transform(),
            color="black",
            linestyle="--",
            zorder=-1,
            linewidth=0.2,
        )

    for axis, secax in zip(axes, secaxes):
        axis.tick_params(
            axis="both",
            which="both",
            direction="in",
            top=True,
            right=True,
            labelsize=tsz,
        )
        secax.tick_params(
            axis="both",
            which="both",
            direction="in",
            top=True,
            right=True,
            labelsize=tsz,
        )
    axes[0].set_title("FFT Analysis", fontsize=tsz + 2)

    # End
    plt.tight_layout()
    plt.savefig("fft_analysis.pdf")
    plt.close()


def make_sinusoid_fit_plots(
    fit_coeffs: astropy.table.Table, filename: pathlib.Path, pixel_period: float
):
    """Create a set of diagnostic plots from the set of sinusoid fits

    _extended_summary_

    Parameters
    ----------
    fit_coeffs : :obj:`~astropy.table.Table`
        The fit coefficients table
    filename : :obj:`~pathlib.Path`
        The filename of the original file, used in the plot title
    pixel_period : :obj:`float`
        The estimated pixel period of the sinusoidal noise from the FFT
    """
    _, axes = plt.subplots(nrows=3, figsize=(6.4, 6.4), gridspec_kw={"hspace": 0})
    tsz = 8

    axes[0].plot(np.arange(len(fit_coeffs)), fit_coeffs["a"])
    axes[0].plot(np.arange(len(fit_coeffs)), fit_coeffs["smooth_a"])
    axes[0].set_ylabel("Sinusoid Amplitude (ADU)")

    axes[1].plot(np.arange(len(fit_coeffs)), fit_coeffs["lambda"])
    axes[1].plot(np.arange(len(fit_coeffs)), fit_coeffs["smooth_lambda"])
    axes[1].hlines(
        pixel_period,
        0,
        1,
        transform=axes[1].get_yaxis_transform(),
        linestyle="--",
        zorder=0,
        color="C2",
    )
    axes[1].set_ylabel("Sinusoid Period (pixels)")

    axes[2].plot(np.arange(len(fit_coeffs)), fit_coeffs["phi"])
    axes[2].plot(np.arange(len(fit_coeffs)), fit_coeffs["smooth_phi"])
    axes[2].set_ylabel("Sinusoid Phase Shift (phase)")
    axes[2].set_xlabel("Image Row Number")

    for axis in axes:
        axis.tick_params(
            axis="both",
            which="both",
            direction="in",
            top=True,
            right=True,
            labelsize=tsz,
        )
    axes[0].set_title(f"Sinusoid Pattern Fits for {filename.name}", fontsize=tsz + 2)

    plt.tight_layout()
    plt.savefig("sinusoid_fits.pdf")
    plt.close()


# Utility Functions (Alphabetical) ===========================================#
def nearest_odd(x: float) -> int:
    """Find the nearest odd integer

    https://www.mathworks.com/matlabcentral/answers/45932-round-to-nearest-odd-integer#accepted_answer_56149

    Parameters
    ----------
    x : :obj:`float`
        Input number

    Returns
    -------
    :obj:`int`
        The nearest odd integer
    """
    return int(2 * np.floor(x / 2) + 1)


def pixper_tofrom_hz(x: np.ndarray) -> np.ndarray:
    """Convert to/from pixel period and Hertz

    _extended_summary_

    Parameters
    ----------
    x : :obj:`~numpy.ndarray`
        Input value(s) to convert

    Returns
    -------
    array-like
        Converted output(s)
    """
    warnings.simplefilter("ignore", RuntimeWarning)
    return 1.0 / (PIX_DWELL * x)


def sinusoid(
    x: np.ndarray,
    a: float,
    lam: float,
    phi: float,
    y0: float,
    lin: float = 0,
    quad: float = 0,
    cube: float = 0,
    quar: float = 0,
) -> np.ndarray:
    """Return a basic sinusoid (for use with :func:`scipy.optimize.curve_fit`)

    _extended_summary_

    Parameters
    ----------
    x : :obj:`~numpy.ndarray`
        The abscissa values for which to return the ordinate
    a : :obj:`float`
        The amplitude of the sinusoid (in units of ordinate)
    lam : :obj:`float`
        The wavelength of the sinusoid (in units of abscissa), equivalent to
        `2π/k` (where `k` is the wavenumber).
    phi : :obj:`float`
        The phase shift of the sinusoid (in units of phase, nominally 0-1)
    y0 : :obj:`float`
        The vertical offset of the sinusoid from zero (in units of ordinate)
    lin : :obj:`float`, optional
        The linear term added to the fit (in units of ordinate/abscissa)
        Default: 0
    quad : :obj:`float`, optional
        The quadratic term added to the fit (in units of ordinate/abscissa**2)
        Default: 0
    cube : :obj:`float`, optional
        The cubic term added to the fit (in units of ordinate/abscissa**3)
        Default: 0
    quar : :obj:`float`, optional
        The quartic term added to the fit (in units of ordinate/abscissa**4)
        Default: 0

    Returns
    -------
    array_like
        The sinusoid ordinate
    """
    return (
        a * np.sin(2.0 * np.pi * x / lam + 2.0 * np.pi * phi)
        + y0
        + lin * x
        + quad * x**2
        + cube * x**3
        + quar * x**4
    )


# Command Line Interface Entry Point =========================================#
def main(files: list | str, use_hann: bool = False, no_plots: bool = False):
    """Main driver

    Main driver function that takes the input list and calls the cleaning
    function for each one.

    Parameters
    ----------
    files : list or str
        File or files to process
    use_hann : bool
        Use a Hann window when cleaning the image for FFT (Default: False)
    no_plots : bool
        Create a plots during the image analysis?  (Default: False)
    """
    # Ensure the input is a list
    if not isinstance(files, list):
        files = [files]

    # Giddy up!
    for file in files:
        iterative_pypeit_clean(pathlib.Path(file).resolve())

        # clean_pickup(
        #     pathlib.Path(file),
        #     use_hann=use_hann,
        #     sine_plot=not no_plots,
        #     fft_plot=not no_plots,
        # )


def entry_point(args=None):
    """Main entry point for the package

    _extended_summary_

    Parameters
    ----------
    args : Any, optional
        Command-line arguments passed in [Defualt: None]
    """

    # Use argparse for the Command-Line Script
    parser = argparse.ArgumentParser(
        prog="scrub_deveny_pickup",
        description="Clean RF pickup noise from DeVeny raw frames",
    )
    parser.add_argument(
        "file",
        nargs="+",
        type=str,
        help="File(s) to clean",
    )
    parser.add_argument(
        "--no_plots",
        default=False,
        action="store_true",
        help="Do not create any plots during the analysis (Default: create plots)",
    )
    res = parser.parse_args(args)

    # Giddy up!
    sys.exit(main(res.file, no_plots=res.no_plots))
