# -*- coding: utf-8 -*-
#
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
#  Created on 05-Dec-2022
#
#  @author: tbowers

"""Clean RF Pickup Noise from DeVeny Spectrographs 2D Images

This module is part of the DeVeny Pickup package, written at
Lowell Observatory.

This module cleans the RF pickup noise from DeVeny spectral 2D images by
fitting a sinusoid to and then subtracting it from each row of the image.

The frequency of the sinusoid tends to vary from frame to frame but is
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
import astropy.wcs
import matplotlib.pyplot as plt
import numpy as np
import scipy.fft
import scipy.ndimage
import scipy.optimize
import scipy.signal
from tqdm import tqdm

# Internal Imports

# Mean pixel readout time (7.25 µs)
PIX_DWELL = 7.25e-6

# Silence a couple known-but-don't-care-about warnings
warnings.simplefilter("ignore", astropy.wcs.FITSFixedWarning)
warnings.simplefilter("ignore", astropy.io.fits.verify.VerifyWarning)


def main(files, use_hann=False, no_plots=False):
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
        clean_pickup(
            pathlib.Path(file),
            use_hann=use_hann,
            sine_plot=not no_plots,
            fft_plot=not no_plots,
        )


def clean_pickup(
    filename: pathlib.Path,
    use_hann=False,
    sine_plot=True,
    fft_plot=False,
):
    """Clean the Pickup Noise

    This is the main functional method of the module.

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
    nrow, ncol = ccd.data.shape

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

    # Compute once, use forever!
    img_mean = np.mean(ccd.data)

    # These are the starting guesses and bounds for the sinusoidal fit
    #  [a, lam, phi, y0, lin, quad]
    p0 = [4, pixel_period, 0.5, 2400]  # + 4 * [0]
    bounds = (
        [0, pixel_period - 10, 0, 2000],  # + 4 * [-np.inf],
        [10, pixel_period + 10, 1, 3000],  # + 4 * [np.inf],
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
        line = ccd.data[img_row, 5:-5]

        # To mitigate cosmic rays and night sky lines, sigma clip @ 5σ
        sig_clip = np.std(line) * 5.0 + np.median(line)
        line[line > sig_clip] = sig_clip

        # Perform the curve fit
        try:
            popt, pcov = scipy.optimize.curve_fit(
                sinusoid, np.arange(line.size), line, p0=p0, bounds=bounds
            )
        except RuntimeError:
            # Reached max function evaluations; set popt and pcov
            popt = p0.copy()
            pcov = np.diag(np.ones(len(popt)))

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
    fit_coeffs = astropy.table.Table(fit_coeffs)

    # Smooth the fit coefficients as a function of row to remove artifacts due
    #  to cosmic rays or strong sources.
    # NOTE: scipy.signal.medfilt() uses zero-padding for the ends of the signal,
    #       therefore subtract the mean and re-add it to eliminate edge weirdness
    fit_coeffs["smooth_a"] = (
        scipy.signal.medfilt(
            fit_coeffs["a"] - (med_a := np.median(fit_coeffs["a"])), kernel_size=21
        )
        + med_a
    )
    fit_coeffs["smooth_lambda"] = (
        scipy.signal.medfilt(
            fit_coeffs["lambda"] - (med_lambda := np.median(fit_coeffs["lambda"])),
            kernel_size=21,
        )
        + med_lambda
    )

    fit_coeffs["smooth_phi"] = scipy.signal.savgol_filter(
        fit_coeffs["phi"], window_length=21, polyorder=2, mode="nearest"
    )
    # Print a happy statement of fact
    print(
        f"This is the mean fit pixel period: {fit_coeffs['smooth_lambda'].mean():.1f} pix"
    )

    # Create the arrays that the processed data will go into
    cleaned_array = ccd.data.copy().astype(np.float64)
    pattern_array = ccd.data.copy().astype(np.float64)

    # Loop through the image and use the smoothed sinusoid fit coefficients
    for img_row, table_row in enumerate(fit_coeffs):
        # Apply the adjusted pattern to the entire row
        line = ccd.data[img_row, :]

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

    # Make a plot, if desired
    if sine_plot:
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
        axes[0].set_title(
            f"Sinusoid Pattern Fits for {filename.name}", fontsize=tsz + 2
        )

        plt.tight_layout()
        plt.savefig("sinusoid_fits.pdf")
        plt.close()

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


def flatten_clean_fft(
    data_array: np.ndarray, imgtyp: str, use_hann=False, fft_plot=False
):
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
            gauss1d,
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


def sinusoid(x, a, lam, phi, y0, lin=0, quad=0, cube=0, quar=0):
    """Return a basic sinusoid (for use with :func:`scipy.optimize.curve_fit`)

    _extended_summary_

    Parameters
    ----------
    x : array_like
        The abscissa values for which to return the ordinate
    a : float
        The amplitude of the sinusoid (in units of ordinate)
    lam : float
        The wavelength of the sinusoid (in units of abscissa), equivalent to
        `2π/k` (where `k` is the wavenumber).
    phi : float
        The phase shift of the sinusoid (in units of phase, nominally 0-1)
    y0 : float
        The vertical offset of the sinusoid from zero (in units of ordinate)
    lin : float, optional
        The linear term added to the fit (in units of ordinate/abscissa)
        Default: 0
    quad : float, optional
        The quadratic term added to the fit (in units of ordinate/abscissa**2)
        Default: 0
    cube : float, optional
        The cubic term added to the fit (in units of ordinate/abscissa**3)
        Default: 0
    quar : float, optional
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


def gauss1d(x, a, mu, sig, y0):
    """Return a basic gaussian (for use with :func:`scipy.optimize.curve_fit`)

    _extended_summary_

    Parameters
    ----------
    x : array_like
        The abscissa values for which to return the ordinate
    a : float
        The amplitude of the gaussian (in units of ordinate)
    mu : float
        The mean of the gaussian (in units of abscissa).
    sig : float
        The width of the gaussian (in units of abscissa).
    y0 : float
        The vertical offset of the gaussian from zero (in units of ordinate)

    Returns
    -------
    array_like
        The gaussian ordinate
    """
    return a * np.exp(-((x - mu) ** 2) / (2.0 * sig**2)) + y0


def create_fft_plot(flat_array, y_fft, x_fft, peak_freq):
    """Create the FFT analysis plot

    _extended_summary_

    Parameters
    ----------
    data_array : _type_
        _description_
    flat_array : _type_
        _description_
    y_fft : _type_
        _description_
    x_fft : _type_
        _description_
    peak_freq : _type_
        _description_
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


def pixper_tofrom_hz(x):
    """Convert to/from pixel period and Hertz

    _extended_summary_

    Parameters
    ----------
    x : array-like
        Input value to convert

    Returns
    -------
    array-like
        Converted output
    """
    warnings.simplefilter("ignore", RuntimeWarning)
    return 1.0 / (PIX_DWELL * x)


# Command Line Interface Entry Point =========================================#
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
