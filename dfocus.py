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
import glob

# 3rd-Party Libraries
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
from tqdm import tqdm

# Local Libraries
from .deveny_grangle import deveny_amag
from .utils import good_poly

# CONSTANTS


def dfocus(flog='last', thresh=100., debug=False):
    """dfocus Find the optimal DeVeny collimator focus value

    [extended_summary]

    Parameters
    ----------
    flog : `str`, optional
        Focus log to process.  If unspecified, process the last sequence in
        the directory.  [Default: 'last']
    thresh : `float`, optional
        Line intensity threshold above background for detection [Default: 100]
    debug : `bool`. optional
        Print debug statements  [Default: False]
    """

    # Parse the log file to obtain file list
    files, focus_id = parse_focus_log(flog)
    nfiles = len(files)

    # Pull the spectrograph setup from the first focus file:
    hdr0 = fits.getheader(f"../{files[0]}" + files[0])
    slitasec = hdr0['SLITASEC']
    grating = hdr0['GRATING']
    grangle = hdr0['GRANGLE']
    lampcal = hdr0['LAMPCAL']

    # Nominal line width
    fnom = 2.94 * slitasec * deveny_amag(grangle)

    # Pull the collimator focus values from the first and last files
    focus_0 = hdr0['COLLFOC']
    focus_1 = (fits.getheader(f"../{files[-1]}"))['COLLFOC']

    # Examine the middle image
    middle_file = f"../{files[int(nfiles/2)]}"

    print(f'\n Processing object image {middle_file}...')
    spectrum, middle_cf = trim_deveny_image(middle_file)

    # Build the trace for spectrum extraction
    n_y, n_x = spectrum.shape
    traces = np.full(n_x, n_y/2, dtype=float).reshape((1,n_x))
    mspectra = extract_spectrum(spectrum, traces, win=11)
    if debug:
        print(f"Traces: {traces}")
        print(f"Middle Spectrum: {mspectra}")

    # Find the lines in the extracted spectrum
    dfl_title = f'{middle_file}   Grating: {grating}   GRANGLE: ' + \
                f'{grangle:.2f}   {lampcal}'

    centers, _ = find_lines(mspectra, thresh=thresh, title=dfl_title)
    n_c = len(centers)
    print(F"Back in the main program, number of lines: {n_c}")
    if debug:
        print(f"Line Centers: {[f'{cent:.1f}' for cent in centers]}")

    # Create an array to hold the FWHM values from all lines from all images
    line_width_array = np.empty((nfiles, n_c), dtype=float)

    # Run through files:
    print("\n Processing arc images...")
    # Show progress bar for processing arcs
    prog_bar = tqdm(total=nfiles, unit='frame', unit_scale=False, colour='yellow')
    for i in range(nfiles):
        spectrum, collfoc = trim_deveny_image(f"../{files[i]}")
        spectra = extract_spectrum(spectrum, traces, win=11)

        # Find FWHM of lines:
        this_cen, fwhm = find_lines(spectra, thresh=thresh, verbose=False)
        line_dx = -4.0 * (collfoc - middle_cf)

        # Keep lines matching this CENTER
        fw = []
        for c in centers:
            d = np.absolute((c+line_dx) - this_cen)
            idx = np.where(d < 3.)[0]
            width = fwhm[idx][0] if len(idx) else np.nan
            fw.append(width if width > 2.0 else np.nan)
        fwhm = np.asarray(fw)
        # fwhm = fit_lines(spectra, centers, boxf=5,
        #                  collfoc=collfoc, middle_cf=middle_cf)
        line_width_array[i,:] = fwhm
        prog_bar.update(1)
    prog_bar.close()


    print(f"Median of the array fw: {np.nanmedian(line_width_array)}")

    # Find the delta between focus values, and the nominal linewidth for focus
    df = (focus_1 - focus_0) / (nfiles - 1)
    fnom = 2.94 * slitasec * deveny_amag(grangle)

    # Fit the focus curve:
    min_focus_index, optimal_focus_index, min_linewidths, fit_pars = \
        fit_focus_curves(line_width_array, fnom=fnom)

    print(optimal_focus_index)

    # Convert the returned indices into actual COLLFOC values, find medians
    print("="*30)
    print(optimal_focus_index.dtype, type(df), type(focus_0))
    min_focus_values = min_focus_index * df + focus_0             # fpos
    optimal_focus_values = optimal_focus_index * df + focus_0     # fwid
    med_min_focus = np.real(np.nanmedian(min_focus_values))       # fwmin
    med_opt_focus = np.real(np.nanmedian(optimal_focus_values))

    print(f"Median Optimal Focus Position: {med_opt_focus:.3f}")

    #=========================================================================#
    # Plot Section

    # The plot shown in the IDL0 window: Plot of the found lines
    _, ax = plt.subplots()
    _ = find_lines(mspectra, thresh=thresh, title=dfl_title, ax=ax)
    plt.tight_layout()
    plt.savefig(f"pyfocus.{focus_id}.eps")

    # The plot shown in the IDL2 window: Plot of best-fit fwid vs centers
    print("="*20)
    print(centers.dtype, optimal_focus_values.dtype, type(med_opt_focus))
    _, ax = plt.subplots()
    tsz = 8
    ax.plot(centers, optimal_focus_values, '.')
    ax.set_xlim(0,2050)
    ax.set_ylim(focus_0, focus_1)
    ax.set_title("Optimal focus position vs. line position, median = " + \
                 f"{med_opt_focus:.3f}")
    ax.hlines(med_opt_focus, 0, 1, transform=ax.get_yaxis_transform(),
              color='magenta', ls='--')

    ax.tick_params('both', labelsize=tsz, direction='in', top=True, right=True)
    plt.tight_layout()
    plt.show()

    # The plot shown in the IDL1 window: Focus curves for each identified line
    plot_focus_curves(centers, line_width_array, min_focus_values,
                      optimal_focus_values, min_linewidths, fit_pars,
                      df, focus_0, fnom=fnom)


def parse_focus_log(flog):
    """parse_focus_log Parse the focus log file produced by the DeVeny LOUI

    [extended_summary]

    Parameters
    ----------
    flog : `str`
        Identifier for the focus log to be processed

    Returns
    -------
    `list`, `str`
        List of files associated with this focus run, and the focus ID
    """
    if flog.lower() == 'last':
        focfiles = sorted(glob.glob('deveny_focus*'))
        flog = focfiles[-1]

    files = []
    with open(flog) as f:
        # Discard file header
        f.readline()
        # Read in the remainder of the file, grabbing just the filenames
        for line in f:
            files.append(line[2:20])

    # Return the list of files, and the FocusID
    return files, flog[13:]


def trim_deveny_image(filename):
    """trim_deveny_image Trim a DeVeny Image

    The IDL code from which this was ported contains a large amount of
    vistigial code from previous versions of the DeVeny camera, including
    instances where the CCD was read out using 2 amplifiers and required
    special treatment in order to balance the two sides of the output image.

    The code below consists of the lines of code that were actually running
    using the keywords passed from current version of dfocus.pro, and the
    pieces of that code that are actually used.

    Specifically, this routine trims off the 50 prescan and 50 postscan pixels,
    as well as several rows off the top and bottom.  (Extracts rows [12:512])

    Parameters
    ----------
    filename : `str`
        Filename of the spectrum to get and trim

    Returns
    -------
    `np.ndarray`
        The trimmed CCD image
    """

    # Parameters for DeVeny (2015 Deep-Depletion Device):
    nxpix, prepix = 2048, 50

    # Read in the file
    with fits.open(filename) as hdul:
        image = hdul[0].data
        collfoc = hdul[0].header['COLLFOC']

    # Trim the image (remove top and bottom rows) -- why this particular range?
    # Trim off the 50 prepixels and the 50 postpixels; RETURN
    return image[12:512,prepix:prepix+nxpix], collfoc


def extract_spectrum(spectrum, traces, win):
    """extract_spectrum Object spectral extraction routine

    Extract spectra by averaging over the specified window

    Parameters
    ----------
    spectrum : `np.ndarray`
        The trimmed spectral image
    traces : `np.ndarray`
        The trace(s) along which to extract spectra
    win : `int`
        Window over which to average the spectrum

    Returns
    -------
    `numpy.ndarray`
        2D or 3D array of spectra of individual orders
    """
    # Spec out the shape, and create an empty array to hold the output spectra
    norders, n_x = traces.shape
    spectra = np.empty((norders, n_x), dtype=float)
    speca = np.empty(n_x, dtype=float)

    # Set extraction window size
    half_window = int(np.floor(win/2))

    for io in range(norders):
        # Because of python indexing, we need to "+1" the upper limit in order
        #   to get the full wsize elements for the average
        trace = traces[io,:].astype(int)
        for i in range(n_x):
            speca[i] = np.average(spectrum[trace[i] - half_window :
                                           trace[i] + half_window + 1, i])
        spectra[io,:] = speca.reshape((1,n_x))

    return spectra


def find_lines(image, thresh=20., findmax=50, minsep=11, fit_window=15,
               verbose=True, ax=None, mark=False, title=''):
    """find_lines Automatically find and centroid lines in a 1-row image

    *** Look into using scipy.signal.find_peaks() for this task! ***

    Parameters
    ----------
    image : `ndarray`
        Extracted spectrum
    thresh : `float`, optional
        Threshold above which to indentify lines [Default: 20 DN above bkgd]
    findmax : `int`, optional
        Maximum number of lines to find [Default: 50]
    minsep : `int`, optional
        Minimum line separation for identification [Default: 11 pixels]
    fit_window : `int`, optional
        Size of the window to fit Gaussian [Default: 15 pixels]
    verbose : `bool`, optional
        Produce verbose output?  [Default: False]
    ax : `pyplot.Axes`, optional
        Create a plot on the provided axes  [Default: None]
    mark : `bool`, optional
        Mark the lines on the output plot?  [Default: False]
    title : `str`
        Title to use for the output plot.  [Default: '']

    Returns
    -------
    centers : `ndarray`
        Line centers (pixel #)
    fwhm : `ndarray`
        The computed FWHM for each peak
    """
    # Get size and flatten to 1D
    _, n_x = image.shape
    spec = np.ndarray.flatten(image)

    # Find background from median value of the image:
    bkgd = np.median(spec)
    if verbose:
        print(f'  Background level: {bkgd:.1f}' + \
            f'   Detection threshold level: {bkgd+thresh:.1f}')

    # Use scipy to find peaks
    centers, _ = signal.find_peaks(spec - bkgd, height=thresh, distance=minsep)
    fwhm = (signal.peak_widths(spec - bkgd, centers))[0]

    if verbose:
        print(f" Number of lines: {len(centers)}")

    # Produce a plot for posterity, if directed
    if ax is not None:
        tsz = 8
        print("At this point the code makes some plots.  Yippee.")
        ax.plot(np.arange(len(spec)), spec)
        ax.set_title(title, fontsize=tsz)
        ax.set_xlabel('CCD Column', fontsize=tsz)
        ax.set_ylabel('I (DN)', fontsize=tsz)
        ax.set_xlim(0, n_x+2)
        ax.set_ylim(0, (yrange := 1.2*max(spec)))
        ax.plot(centers, spec[centers.astype(int)]+0.02*yrange, 'k*')
        for c in centers:
            ax.text(c, spec[int(np.round(c))]+0.03*yrange, f"{c:.3f}",
                    fontsize=tsz)
        ax.tick_params('both', labelsize=tsz, direction='in',
                       top=True, right=True)

    return centers, fwhm


def fit_focus_curves(fwhm, fnom=2.7, norder=2):
    """fit_focus_curves Fit line / star focus curves

    [extended_summary]

    Parameters
    ----------
    fwhm : `np.ndarray`
        Array of FWHM for all lines as a function of COLLFOC
    fnom : `float`, optional
        Nominal FHWM of an in-focus line. [Default: 2.7]
    norder : `int`, optional
        Polynomial order of the focus fit [Default: 2 = Quadratic]

    Returns
    -------
    `1d_array`, `1d_array`, `1d_array`, `2d_array`
        Tuple of best fit focus, best fit linewidth, the minumum linewidth,
        and the actual fit parameters (for plotting)
    """
    # Create the various arrays / lists needed
    n_focus, n_centers = fwhm.shape
    min_linewidth = []
    min_cf_idx_value = []
    optimal_cf_idx_value = []
    foc_fits = []

    #print("Here's the nasty array we're fitting!")
    #print(fwhm)

    # Fitting arrays (these are indices for collimator focus)
    cf_idx_coarse = np.arange(n_focus, dtype=float)
    cf_idx_fine = np.arange(0, n_focus-1 + 0.1, 0.1, dtype=float)

    # Loop through lines to find the best focus for each one
    for i in range(n_centers):

        # Data are the FWHM for this line at different COLLFOC
        fwhms_of_this_line = fwhm[:,i]

        # Find unphysically large or small FWHM (or NaN) -- set to np.nan
        bad_idx = np.where(np.logical_or(fwhms_of_this_line < 1.0,
                           fwhms_of_this_line > 15.0))
        fwhms_of_this_line[bad_idx] = np.nan
        fwhms_of_this_line[np.isnan(fwhms_of_this_line)] = np.nan

        # If more than 3 of the FHWM are bad for this line, skip and go on
        if len(bad_idx) > 3:
            # Add values to the lists for proper indexing
            for flst in [min_linewidth, min_cf_idx_value, optimal_cf_idx_value]:
                flst.append(None)
            continue

        # Do a polynomial fit (norder) to the FWHM vs COLLFOC index
        #fit = np.polyfit(cf_idx_coarse, fwhms_of_this_line, norder)
        fit = good_poly(cf_idx_coarse, fwhms_of_this_line, norder, 2.)
        foc_fits.append(fit)
        print(f"In fit_focus_curves(): fit = {fit}")

        # Use the fine grid to evaluate the curve miniumum
        focus_curve = np.polyval(fit, cf_idx_fine)                   # fitfine
        min_cf_idx_value.append(cf_idx_fine[np.argmin(focus_curve)]) # focus
        min_linewidth.append(np.min(focus_curve))                    # minfoc

        # Compute the nominal focus position as the larger of the two points
        #  where the polymonial function crosses fnom
        a0, a1, a2 = fit[0], fit[1], fit[2]-fnom
        print(f"Roots: {np.roots([a0,a1,a2])}")
        optimal_cf_idx_value.append(np.max(np.real(np.roots([a0,a1,a2]))))             # fwidth
        #print(fits[-1])

    return np.asarray(min_cf_idx_value), np.asarray(optimal_cf_idx_value), \
           np.asarray(min_linewidth), np.asarray(foc_fits)


def plot_focus_curves(centers, line_width_array, min_focus_values,
                      optimal_focus_values, min_linewidths, fit_pars,
                      df, focus_0, fnom=2.7):
    """plot_focus_curves Make the big plot of all the focus curves

    [extended_summary]

    Parameters
    ----------
    centers : `array`
        List of line centers from fine_lines()
    line_width_array : `2d array`
        Array of line widths from each COLLFOC setting for each line
    min_focus_values : `array`
        List of the minimum focus values found from the polynomial fit
    optimal_focus_values : `array`
        List of the optimal focus values found from the polynomial fit
    min_linewidths : `array`
        List of the minumum linewidths found from the fitting
    fit_pars : `2d array`
        Array of the polynomial fit parameters for each line
    df : `float`
        Spacing between COLLFOC settings
    focus_0 : `float`
        Lower end of the COLLFOC range
    fnom : `float`, optional
        Nominal (optimal) linewidth [Default: 2.7]
    """

    # Set up variables
    n_foc, n_c = line_width_array.shape
    focus_idx = np.arange(n_foc)
    fx = focus_idx * df + focus_0

    # Set the plotting array
    ncols = 6
    nrows = np.floor(n_c/ncols).astype(int) + 1
    fig, axs = plt.subplots(ncols=ncols, nrows=nrows, figsize=(8.5,11))
    tsz = 6     # type size

    for i, ax in enumerate(axs.flatten()):
        if i < n_c:
            ax.plot(fx, line_width_array[:,i], 'k^')
            ax.plot(fx, np.polyval(fit_pars[i,:], focus_idx), 'g-')
            ax.vlines(min_focus_values[i], 0, min_linewidths[i], color='r', ls='-')
            ax.vlines(optimal_focus_values[i], 0, fnom, color='b', ls='-')
            ax.set_ylim(0, np.nanmax(line_width_array[:,i]))
            ax.set_xlim(np.min(fx)-df, np.max(fx)+df)
            ax.set_xlabel('Collimator Position (mm)', fontsize=tsz)
            ax.set_ylabel('FWHM (pix)', fontsize=tsz)
            ax.set_title(f"LC: {centers[i]:.2f}  Fnom: {fnom:.2f} pixels", fontsize=tsz)
            ax.tick_params('both', labelsize=tsz, direction='in', top=True, right=True)
            ax.grid(which='both', color='#c0c0c0', linestyle='-', linewidth=0.5)
        else:
            # Clear any extra positions if there aren't enough lines
            fig.delaxes(ax)

    plt.tight_layout()
    plt.show()


def find_lines_in_spectrum(filename, thresh=100.):
    """find_lines_in_spectrum Find the line centers in a spectrum

    This function is not directly utilized in DFOCUS, but rather is included
    as a wrapper for several functions that can be used by other programs.

    Given the filename of an arc-lamp spectrum, this function returns a list
    of the line centers found in the image.

    Parameters
    ----------
    filename : `str`
        Filename of the arc frame to find lines in
    thresh : `float`, optional
        Line intensity threshold above background for detection [Default: 100]

    Returns
    -------
    `array`
        List of line centers found in the image
    """
    # Get the trimmed image
    spectrum, _ = trim_deveny_image(filename)

    # Build the trace for spectrum extraction
    n_y, n_x = spectrum.shape
    traces = np.full(n_x, n_y/2, dtype=float).reshape((1,n_x))
    spectra = extract_spectrum(spectrum, traces, win=11)

    # Find the lines!
    centers, _ = find_lines(spectra, thresh=thresh)

    return centers


#=========================================================#
def main(args):
    """main [summary]

    [extended_summary]

    Parameters
    ----------
    args : [type]
        [description]
    """
    dfocus()


if __name__ == '__main__':
    import sys
    main(sys.argv)
