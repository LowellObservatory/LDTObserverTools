# -*- coding: utf-8 -*-
#
#  This file is part of PyDeVeny.
#
#   This Source Code Form is subject to the terms of the Mozilla Public
#   License, v. 2.0. If a copy of the MPL was not distributed with this
#   file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
#  Created on 01-Feb-2021
#
#  @author: tbowers

"""PyDeVeny contains python ports of the various DeVeny IDL routines

Lowell Discovery Telescope (Lowell Observatory: Flagstaff, AZ)
http://www.lowell.edu

This file contains the dfocus routine for computing the required collimator
focus for the DeVeny Spectrograph based on a focus sequence completed by the
DeVeny LOUI.
"""

# Built-In Libraries
import glob
import warnings

# 3rd-Party Libraries
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize
from scipy import signal
#import PySimpleGUI as sg

# Local Libraries
from .deveny_grangle import deveny_amag

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
    files, _ = parse_focus_log(flog)
    nfiles = len(files)

    # Pull the spectrograph setup from the first focus file:
    hdr0 = fits.getheader(f'../' + files[0])
    slitasec = hdr0['SLITASEC']
    grating = hdr0['GRATING']
    grangle = hdr0['GRANGLE']
    lampcal = hdr0['LAMPCAL']

    # Nominal line width
    fnom = 2.94 * slitasec * deveny_amag(grangle)

    # Pull the collimator focus values from the first and last files
    focus_0 = hdr0['COLLFOC']
    focus_1 = (fits.getheader(f'../' + files[-1]))['COLLFOC']

    # Examine the middle image
    middle_file = '../' + files[int(nfiles/2)]

    print(f'\n Processing object image {middle_file}...')
    spectrum = trim_deveny_image(middle_file)

    # Build the trace for spectrum extraction
    ny, nx = spectrum.shape
    traces = np.full(nx, ny/2, dtype=float).reshape((1,nx))
    mspectra = extract_spectrum(spectrum, traces, win=11)
    if debug:
        print(f"Traces: {traces}")
        print(f"Middle Spectrum: {mspectra}")

    # Find the lines in the extracted spectrum
    dfl_title = f'{middle_file}   Grating: {grating}   GRANGLE: ' + \
                f'{grangle:.2f}   {lampcal}'

    centers, _ = find_lines(mspectra, thresh=thresh, title=dfl_title)
    nc = len(centers)
    if debug:
        print(F"Back in the main program, number of lines: {nc}")
    print(f"Line Centers: {[f'{cent:.1f}' for cent in centers]}")

    # Create an array to hold the FWHM values from all lines from all images
    line_width_array = np.empty((nfiles, nc), dtype=float)

    # Run through files:
    for i in range(nfiles):
        print(f"\n Processing arc image {files[i]} ...")
        spectrum = trim_deveny_image(f"../{files[i]}")
        print(f"  Extracting spectra from image {files[i]}...")
        spectra = extract_spectrum(spectrum, traces, win=11)

        # Find FWHM of lines:
        fwhm = fit_lines(spectra, centers)
        line_width_array[i,:] = fwhm

    print(f"Median of the array fw: {np.nanmedian(line_width_array)}")

    # Find the delta between focus values, and the nominal linewidth for focus
    df = (focus_1 - focus_0) / (nfiles - 1)
    fnom = 2.94 * slitasec * deveny_amag(grangle)

    # Fit the focus curve:
    min_focus_index, optimal_focus_index, _ = \
        fit_focus_curves(line_width_array, fnom=fnom)
    min_focus_value = min_focus_index * df + focus_0
    optimal_focus_value = optimal_focus_index * df + focus_0
    med_min_focus = np.real(np.nanmedian(min_focus_value))
    med_opt_focus = np.real(np.nanmedian(optimal_focus_value))

    print(f"Median focus position: {med_min_focus:.3f}, " + \
          f"Median Full Width: {med_opt_focus:.3f}")

    #=========================================================================#
    # Plot Section

    # The plot shown in the IDL0 window: Plot of the found lines
    fig, ax = plt.subplots()
    _ = find_lines(mspectra, thresh=thresh, title=dfl_title, ax=ax)
    plt.tight_layout()
    plt.show()

    # The plot shoen in the IDL2 window: Plot of best-fit fwid vs centers
    fig, ax = plt.subplots()
    tsz = 8
    ax.plot(centers, optimal_focus_value, '.')
    ax.set_xlim(0,2050)
    ax.set_ylim(focus_0, focus_1)
    ax.set_title("Optimal focus position vs. line position, median = " + \
                 f"{med_opt_focus:.3f}")
    ax.hlines(med_opt_focus, 0, 1, transform=ax.get_yaxis_transform(), 
              color='magenta', ls='--')
    # ax.hlines(9.411, 0, 1, transform=ax.get_yaxis_transform(), 
    #           color='magenta', ls='--')

    ax.tick_params('both', labelsize=tsz, direction='in', top=True, right=True)
    plt.tight_layout()
    plt.show()

    # The plot shoen in the IDL2 window: Plot of best-fit fwid vs centers
    plot_focus_curves()

    """
    dplotfocus, x, fw, fits, flist, fnom=fnom
    plot, centers, fpos, xra=[0,2050], yra=[f0, f1], xsty=1, $
        /ynoz, ticklen=1, psym=5, $
        title='Minimum focus position vs. line position, median = ' + $
        strmid(mfpos, 3, 8), xtitle='CCD column', ytitle='Focus (mm)'
    plot, centers, fwid, xra=[0,2050], yra=[f0, f1], xsty=1, $
        /ynoz, ticklen=1, psym=5, $
        title='Optimal focus position vs. line position, median = ' + $
        strmid(mfwid, 3, 8), xtitle='CCD column', ytitle='Optimal Focus (mm)', $
        subtitle='Grating: ' + grating + $
        '   Slit width:' + strmid(slitasec, 4, 6) + ' arcsec' + $
        '    Nominal line width:' + strmid(fnom,4,7) + ' pixels'
    plots, [0, 2050], [mfwid, mfwid], color=60, thick=3, /data
    setplot,'x'

    window, 0, xsize=750, ysize=450, xpos=50, ypos=725
    centers=dflines(mspectra, fwhm=fwhm, thresh=thresh, $
        title = mfile + '   Grating:  ' + grating + '   GRANGLE:  ' + $
        strtrim(grangle,2) + '   ' + lampcal)
    window, 1, xsize=1500, ysize=650, xpos=50, ypos=50
    dplotfocus, x, fw, fits, flist, fnom=fnom
    window, 2, xsize=750, ysize=450, xpos=805, ypos=725
    plot, centers, fwid, xra=[0,2050], yra=[f0, f1], xsty=1, $
        /ynoz, ticklen=1, psym=5, $
        title='Optimal focus position vs. line position, median = ' + $
        strmid(mfwid, 3, 8), xtitle='CCD column', ytitle='Optimal Focus (mm)', $
        subtitle='Grating: ' + grating + $
        '   Slit width:' + strmid(slitasec, 4, 6) + ' arcsec' + $
        '    Nominal line width:' + strmid(fnom,4,7) + ' pixels'
    plots, [0, 2050], [mfwid, mfwid], color=60, thick=3, /data
    loadct, 0

    return
    end
    """

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
    image = fits.getdata(filename)

    # Trim the image (remove top and bottom rows) -- why this particular range?
    # Trim off the 50 prepixels and the 50 postpixels; RETURN
    return image[12:512,prepix:prepix+nxpix]


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
    norders, nx = traces.shape
    spectra = np.empty((norders, nx), dtype=float)
    speca = np.empty(nx, dtype=float)

    # Set extraction window size
    half_window = int(np.floor(win/2))

    for io in range(norders):
        # Because of python indexing, we need to "+1" the upper limit in order
        #   to get the full wsize elements for the average
        trace = traces[io,:].astype(int)
        for i in range(nx):
            speca[i] = np.average(spectrum[trace[i] - half_window :
                                           trace[i] + half_window + 1, i])
        spectra[io,:] = speca.reshape((1,nx))

    return spectra


def find_lines(image, thresh=20., findmax=50, minsep=11, fit_window=15,
               verbose=False, ax=None, mark=False, title=''):
    """find_lines Automatically find and centroid lines in a 1-row image

    [extended_summary]

    Parameters
    ----------
    image : `np.ndarray`
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
    (`list of float`, `list of float`)
        centers: Line centers (pixel #)
        fwhm: The computed FWHM
    """
    # Silence OptimizeWarning, this function only
    warnings.simplefilter('ignore', optimize.OptimizeWarning)

    # Define the half-window
    fhalfwin = int(np.floor(fit_window/2))

     # Get size and flatten to 1D
    _, nx = image.shape
    spec = np.ndarray.flatten(image)

    # Find background from median value of the image:
    bkgd = np.median(spec)
    print(f'  Background level: {bkgd:.1f}' + \
          f'   Detection threshold level: {bkgd+thresh:.1f}')

    # Create empty lists to fill
    cent, fwhm = ([], [])
    j0 = 0

    # Step through the cut and identify peaks:
    for j in range(nx):

        # If we get too close to the end, skip
        if j > (nx - minsep):
            continue

        # If the spectrum at this pixel is above the THRESH...
        if spec[j] > (bkgd + thresh):

            # Mark this pixel as j1
            j1 = j

            # If this is too close to the last one, skip
            if np.abs(j1 - j0) < minsep:
                continue

            # Loop through 0-FINDMAX...  (find central pixel?)
            for jf in range(findmax):
                itmp0 = spec[jf + j]
                itmp1 = spec[jf + j + 1]
                if itmp1 < itmp0:
                    icntr = jf + j
                    break

            # If central pixel is too close to the edge, skip
            if (icntr < minsep/2) or (icntr > (nx - minsep/2 - 1)):
                continue

            # Set up the gaussian fitting for this line
            xmin, xmax = (icntr - fhalfwin, icntr + fhalfwin + 1)
            xx = np.arange(xmin, xmax, dtype=float)
            temp = spec[xmin : xmax]
            # Filter the SPEC to smooth it a bit for fitting
            temp = signal.medfilt(temp, kernel_size=3)

            # Run the fit, with error checking
            try:
                p0 = [1000, np.mean(xx), 3, bkgd]
                aa, _ = optimize.curve_fit(gaussfit_func, xx, temp, p0=p0)
            except RuntimeError:
                continue  # Just skip this one

            # If the width makes sense, save
            if (fw := aa[2] * 2.355) > 1.0:      # sigma -> FWHM
                cent.append(aa[1])
                fwhm.append(fw)

            # Set j0 to this pixel before looping on
            j0 = jf + j

    # Make list into an array, check again that the centers make sense
    centers = np.asarray(cent)
    c_idx = np.where(np.logical_and(centers > 0, centers <= nx))
    centers = centers[c_idx]

    if verbose:
        print(f" Number of lines: {len(centers)}")

    # Produce a plot for posterity, if directed
    if ax is not None:
        tsz = 8
        print(f"At this point the code makes some plots.  Yippee.")
        ax.plot(np.arange(len(spec)), spec)
        ax.set_title(title, fontsize=tsz)
        ax.set_xlabel('CCD Column', fontsize=tsz)
        ax.set_ylabel('I (DN)', fontsize=tsz)
        ax.set_xlim(0, nx+2)
        ax.set_ylim(0, (yrange := 1.2*max(spec)))
        ax.plot(centers, spec[centers.astype(int)]+0.02*yrange, 'k*')
        for c in centers:
            ax.text(c, spec[int(np.round(c))]+0.03*yrange, f"{c:.1f}", fontsize=tsz)
        ax.tick_params('both', labelsize=tsz, direction='in', top=True, right=True)

    return centers, fwhm


def fit_lines(spectrum, centers, return_amp=False, boxf=None):
    """fit_lines Procedure to fit line profiles in a focus image

    [extended_summary]

    Parameters
    ----------
    spectrum : `np.ndarray`
        Extracted spectrum
    centers : `list of float`
        Line centers based on the 'middle' image
    return_amp : `bool`, optional
        Return the amplitudes of the lines? [Default: False]
    boxf : [type], optional
        [description], by default None

    Returns
    -------
    `np.ndarray`
        List of the FHWM of the lines from this image
    """
    # Process inputs
    if return_amp:
        afit = []
        amax = []
    if boxf is None:
        boxi = 17
        boxf = 11
    else:
        boxi = boxf + 6
    xx = np.arange(2 * boxf + 1).flatten()

    fwhm = []
    # Loop through the lines!
    for center in centers:
        # Check that we're within a valid range
        if (center - boxi) < 0 or (center + boxi) > 2030:
            continue

        # Compute the 1st moment to estimate the window for Gaussian fitting
        box = spectrum[:,int(center) - boxi : int(center) + boxi + 1]
        ccnt = first_moment_1d(box) + center - boxi
        if (ccnt - boxi) < 0 or (ccnt + boxi) > 2030:
            continue

        # This is the box for Gaussian fitting
        box = spectrum[:,int(ccnt) - boxf : int(ccnt) + boxf + 1].flatten()

        # Run the fit, with error checking
        try:
            p0 = [1000, np.mean(xx), 3, 0]
            a, _ = optimize.curve_fit(gaussfit_func, xx, box, p0=p0)
        except RuntimeError:
            # Add "None" to the lists to maintain indexing
            fwhm.append(None)
            if return_amp:
                afit.append(None)
                amax.append(np.max(box))
            continue

        # Add the fit values to the appropriate arrays
        fwhm.append(a[2] * 2.355)     # FWHM = 2.355 * sigma
        if return_amp:
            afit.append(a[0])
            amax.append(np.max(box))

    # Return either the single FWHM or tuple containing FWHM & amplitudes
    if return_amp:
        return np.asarray(fwhm), np.asarray(afit), np.asarray(amax)
    else:
        return np.asarray(fwhm)


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
    `array`, `array`, `array`
        Tuple of best fit focus, best fit linewidth, and the minumum
        linewidth
    """
    # Create the various arrays / lists needed
    n_focus, n_centers = fwhm.shape
    min_linewidth = []
    min_focus_index = []
    optimal_focus_index = []

    print("Here's the nasty array we're fitting!")
    print(fwhm)

    # Fitting arrays
    x_coarse = np.arange(n_focus, dtype=float)
    x_fine = np.arange(0, n_focus-1 + 0.1, 0.1, dtype=float)

    # Loop through lines to find the best focus for each one
    for i in range(n_centers):

        # Data are the FWHM for this line at different COLLFOC
        data = fwhm[:,i]

        # Find unphysically large or small FWHM (or NaN) -- set to 50.
        bad_idx = np.where(np.logical_or(data < 1.0, data > 15.0))
        data[bad_idx] = 50.
        data[np.isnan(data)] = 50.

        # If more than 3 of the FHWM are bad for this line, skip and go on
        if len(bad_idx) > 3:
            # Add values to the lists for proper indexing
            for foc_list in [min_linewidth, min_focus_index, optimal_focus_index]:
                foc_list.append(None)
            continue

        # Do the fit to the FWHM vs COLLFOC
        fit = np.polyfit(x_coarse, data, norder)

        # Use the fine grid to evaluate the curve miniumum
        focus_curve = np.polyval(fit, x_fine)
        min_focus_index.append(x_fine[np.argmin(focus_curve)])
        min_linewidth.append(np.min(focus_curve))

        # Compute the nominal focus position as the larger of the two points
        #  where the polymonial crosses fnom
        fit[2] -= fnom
        focus_roots = np.roots(fit)
        optimal_focus_index.append(np.max(focus_roots))

    return np.asarray(min_focus_index), np.asarray(optimal_focus_index), \
           np.asarray(min_linewidth)


def plot_focus_curves():
    # This will be the Python translation of dplotfocus
    """
    pro dplotfocus, x, fwhm, fits, filelist, fnom=fnom

    common curves, fx, centers, fpos, fwid, fwmin
    
    sz = size(fwhm) & nc = sz(1)
    !p.multi(1)=6
    !p.multi(2)=nc/6+1
    if (nc/4+1) gt 5 then !p.multi(2)=5
    for i=0,nc-1 do begin

        plot, fx, fwhm(i,*), psym=4, yra=[0,max(fwhm(i,*))], $
            ticklen=1,ysty=1, ymargin=[4, 4], $
            xtitle='Collimator Position (mm)', $
            ytitle='FWHM', $
            title='LC: '+strtrim(centers(i),1)+'  Fnom: '+strmid(fnom,6,5)+' pixels'
        oplot, fx, poly(x,fits(*,i))
        plots, [fpos(i), fpos(i)], [0, fwmin(i)], /data, color=90, thick=3
        plots, [fwid(i), fwid(i)], [0, fnom], /data, color=60, thick=3

    endfor

    !p.multi=0

    return
    end
    """
    pass


def gaussfit_func(x, a0, a1, a2, a3):
    """gaussfit_func Gaussian Function

    [extended_summary]

    Parameters
    ----------
    x : `array`
        X values over which to compute the gaussian
    a0 : `float`
        Gaussian amplitude
    a1 : `float`
        Gaussian mean (mu)
    a2 : `float`
        Gaussian width (sigma)
    a3 : `float`
        Baseline atop which the Gaussian sits

    Returns
    -------
    `array`
        The Y values of the Gaussian corresponding to X
    """
    # Silence RuntimeWarning for overflow, this function only
    warnings.simplefilter('ignore', RuntimeWarning)
    z = (x - a1) / a2
    return a0 * np.exp(-z**2 / 2.) + a3


def first_moment_1d(line):
    """first_moment_1d Returns the 1st moment of line

    [extended_summary]

    Parameters
    ----------
    line : `array`
        1-dimensional array to find the 1st moment of

    Returns
    -------
    `float`
        The first moment of the input array relative to element #
    """
    # Only use positive values -- set negative values to zero
    line[np.where(line < 0)] = 0

    # Make the counting array
    yy = np.arange(len(line))

    # Return the first moment
    return np.sum(yy * line) / np.sum(line)


#=========================================================#
def main(args):
    dfocus()


if __name__ == '__main__':
    import sys
    main(sys.argv)
