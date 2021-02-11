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
import sys
import warnings

# Numpy & SciPy & PySimpleGui
import numpy as np
from scipy import optimize
from scipy import signal
import PySimpleGUI as sg

# AstroPy
from astropy.io import fits

# Local Libraries
from .deveny_grangle import deveny_amag

# CONSTANTS


def dfocus(flog='last', thresh=100.):
    """Find the optimal DeVeny collimator focus value

    if log is unspecified, process the last sequence
    """

    #global curves, fx, centers, fpos, fwid, fwmin

    # Parse the log file to obtain file list;
    files, focid = flogparse(flog)
    nfiles = len(files)

    # Pull the spectrograph setup from the first focus file:
    hdr0 = fits.getheader(f'../' + files[0])
    slitasec = hdr0['SLITASEC']
    grating = hdr0['GRATING']
    grangle = hdr0['GRANGLE']
    lampcal = hdr0['LAMPCAL']
    filtrear = hdr0['FILTREAR']
    fnom = 2.94 * slitasec * deveny_amag(grangle)

    # Pull the collimator focus values from the first and last files
    f0 = hdr0['COLLFOC']
    hdr1 = fits.getheader(f'../' + files[-1])
    f1 = hdr1['COLLFOC']

    # Set up arrays for populating the focus values
    df = (f1 - f0) / (nfiles - 1)
    x = np.arange(nfiles, dtype=float)
    fx = np.arange(nfiles, dtype=float) * df + f0

    # Examine the middle image
    mfile = '../' + files[int(nfiles/2)]
    swext=4  # Vestigial extraction pramater
    win = 11 # Something about the window to extract
    print(f'\n Processing object image {mfile}...')
    spectrum = dvtrim(mfile)
    ny, nx = spectrum.shape
    #print(f'The shape of spectrum: {spectrum.shape}')
    traces = np.full(nx, ny/2, dtype=float).reshape((1,nx))
    #print(type(traces.shape))
    #print(f'The shape of traces: {traces.shape}')
    print(traces)
    mspectra = dextract(spectrum, traces, win, swext=swext)
    print(mspectra)

    # Find the lines:
    dfl_title = f'{mfile}   Grating:  {grating}   GRANGLE:  ' + \
                f'{grangle:.0f}   {lampcal}'
    centers, fwhm = dflines(mspectra, thresh=thresh, title = dfl_title)
    nc = len(centers)
    print(F"Back in the main program, number of lines: {nc}")
    print(f"Line Centers: {[f'{cent:.1f}' for cent in centers]}")
    fw = np.empty((nfiles, nc), dtype=float)

    # Run through files:
    for i in range(nfiles):
        print(f"\n Processing arc image {files[i]} ...")
        spectrum = dvtrim(f"../{files[i]}")
        print(f"  Extracting spectra from image {files[i]}...")
        spectra = dextract(spectrum, traces, win, swext=swext)

        # Find FWHM of lines:
        fwhm = dfitlines(spectra, centers)
        fw[i,:] = fwhm

    print(f"Median of the array fw: {np.median(fw)}")

    # Fit the lines:
    (linefits, (focus, fwidth, minfoc)) = dfitcurves(fw, fnom=fnom)
    fpos = focus * df + f0
    fwid = fwidth * df + f0
    fwmin = minfoc
    mfpos = np.median(fpos)
    mfwid = np.median(fwid)

    print(f"mfpos: {mfpos}, mfwid: {mfwid}")


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

def flogparse(flog):
    """Parse the focus log file produced by the DeVeny LOUI

    flog: if 'last', then read in the most recent log file, but
            the user can also pass in a specific log file to use
    """
    if flog.lower() == 'last':
        focfiles = sorted(glob.glob('deveny_focus*'))
        flog = focfiles[-1]
    focid = flog[13:]
 
    files = []
    with open(flog) as f:
        f.readline()   # Discard file header
        for line in f:
            files.append(line[2:20])
    
    return files, focid



def dvtrim(filename):
    """Trim a DeVeny Image

    The IDL code from which this was ported contains a large amount of
    vistigial code from previous versions of the DeVeny camera, including
    instances where the CCD was read out using 2 amplifiers and required
    special treatment in order to balance the two sides of the output image.

    The code below consists of the lines of code that were actually running 
    using the keywords passed from current version of dfocus.pro, and the
    pieces of that code that are actually used.
    """

    # Parameters for DeVeny (2015 Deep-Depletion Device):
    nxpix = 2048
    prepix = 50
    # postpix = 50        # Overscan pixels

    # Read in the file
    image = fits.getdata(filename)

    # Trim the image (remove top and bottom rows) -- why this particular range?
    image = image[12:512,:]

    # Trim off the 50 prepixels and the 50 postpixels; RETURN
    return image[:,prepix:prepix+nxpix]


def dextract(spectrum,traces,nspix,swext=2,npixavg=None):
    """Object spectral extraction routine
   
    Options:
    swext = 1: extract point source spectra by averaging over window
    swext = 2: extract full spatial orders with spatial interpolation
    swext = 3: extract full orders without interpolation 
    swext = 4: extract spectra by averaging over window

    Output:
    spectra: 2 or 3-d array of spectra of individual orders
    """
    
    # Case out the dimensionality of traces... 0 -> return
    if traces.ndim == 0:
        return 0
    elif traces.ndim == 1:
        nx = traces.size
        norders = 1
    else:
        norders, nx = traces.shape

    # Case out the option swext
    if swext == 0:
        return 0

    elif swext == 1:
        if npixavg is None:
            npixavg = 19
        spectra = np.empty((norders,nx), dtype=float)
        for io in range(norders):
            spectra[io,:] = dspecfit(spectrum, traces[io,:], bwidth=nspix, 
                                     extwidth=npixavg)
    
    elif swext == 2:
        spectra = np.empty((nspix, nx, norders), dtype=float)
        fnspix = np.arange(nspix, dtype=float) - nspix/2
        for io in range(norders):
            for ix in range(nx):
                # Interpolation:
                xt = traces[io,ix].astype(int) + fnspix
                ut = traces[io, ix] + fnspix
                vector = spectrum[xt, ix]
                tvector = interpol(vector, xt, ut)
                spectra[:,ix,io] = tvector

    elif swext == 3:
        spectra = np.empty((nspix, nx, norders), dtype=float)
        inspix = np.arange(nspix, dtype=int) - nspix/2
        for io in range(norders):
            for ix in range(ix):
                # Interpolation
                xt = traces[io,ix].astype(int) + inspix
                spectra[:,ix,io] = spectrum[xt, ix]

    elif swext == 4:
        if npixavg is None:
            npixavg = nspix
        spectra = np.empty((norders, nx), dtype=float)

        #print(f'The shape of spectra: {spectra.shape}')
        #print(f'The shape of traces: {traces.shape}')

        for io in range(norders):
            spectra[io,:] = specavg(spectrum, traces[io,:], npixavg)
    
    else:
        print("Silly user, you can't do that.")
        return 0

    return spectra


def specavg(spectrum, trace, wsize):
    """Extract an average spectrum along trace of size wsize
    :param spectrun: input spectrum
    :param trace: the trace along which to extract
    :param wsize: the size of the extraction (usually odd)
    :return:
    """
    # Case out the dimensionality of traces... 0 -> return
    if spectrum.ndim == 0:
        return 0
    elif spectrum.ndim == 1:
        nx = spectrum.size
    else:
        ny, nx = spectrum.shape
    speca = np.empty(nx, dtype=float)

    whalfsize = int(np.floor(wsize/2))

    # Because of python indexing, we need to "+1" the upper limit in order
    #   to get the full wsize elements for the average
    for i in range(nx):
        speca[i] = np.average(spectrum[int(trace[i]) - whalfsize : 
                                      int(trace[i]) + whalfsize + 1, i])
    
    #print(f"The shape of speca: {speca.shape}")
    return speca.reshape((1,nx))


def dspecfit(spec, trace, bwidth=101, extwidth=19):
    """Fit the background to a spectrum and subtract it

    Default values:
    :bwidth: 101
    :extwidth: 19
    """

    """     
    tval=fix(trace)

    nfit=bwidth & win=[bwidth/2-extwidth/2,bwidth/2+extwidth/2] ;  w=0 over win
    sz=size(spec)
    fits=fltarr(nfit,sz(1)) & datas=fits & subs=fits

    for i=0,sz(1)-1 do begin

        ww=fltarr(nfit) & ww(*)=1. & ww(win(0):win(1))=0.
        data=spec(i,tval(i)-nfit/2:tval(i)+nfit/2)
        coef=polyfitw(findgen(nfit),data,ww,1)	; go with linear fit

        fit=poly(findgen(nfit),coef)
        fits(*,i)=fit
        datas(*,i)=data
        subs(*,i)=data-fit

    endfor

    gplott=fltarr(sz(1))

    for i=0,sz(1)-1 do gplott(i)=total(subs(win(0):win(1),i))

    return,gplott

    end
    """
    return 0


def gaussfit_func(x, a0, a1, a2, a3):
 
    # Silence RuntimeWarning for overflow, this function only
    warnings.simplefilter('ignore', RuntimeWarning)
    z = (x - a1) / a2
    y = a0 * np.exp(-z**2 / a2) + a3
    return y

def dflines(image, thresh=20., mark=False, title=''):
    """Automatically find and centroid lines in a 1-row image
 
    :image:
    :thresh: 20 DN above background
    :mark:
    :title:

    Returns:
    :centers: List of line centers (pixel #)
    :fwhm: The computed FWHM
    """

    warnings.simplefilter('ignore', optimize.OptimizeWarning)

    nx, ny = image.shape
    avgj = np.ndarray.flatten(image)

    print(f"Shapes of image: {image.shape}, and avgj: {avgj.shape}")

    peaks = []
    fwhm = []
 
    # Create background from median value of the image:
    bkgd = np.median(image)
    print(f'  Background level: {bkgd:.1f}')

    # Step through the cut and identify peaks:
    fmax = 11     # Minimum separations
    fwin = 15
    fhalfmax = int(np.floor(fmax/2))
    fhalfwin = int(np.floor(fwin/2))
    findmax = 50
    j0 = 0

    if mark:
        centers = [0]
    else:
        for j in range(ny):
            # print(f"avgj[j]: {avgj[j]}, {avgj[j].shape}")
            if j > (ny - fmax):
                continue
            if avgj[j] > (bkgd + thresh):
                j1 = j
                if np.abs(j1 - j0) < fmax:      
                    continue
                for jf in range(findmax):
                    itmp0 = avgj[jf + j]
                    itmp1 = avgj[jf + j + 1]
                    if itmp1 < itmp0:
                        icntr = jf + j
                        break
                if (icntr < fmax/2) or (icntr > (ny - fmax/2 - 1)):
                    continue
                xx = np.arange(fwin, dtype=float) + float(icntr - fhalfwin)
                temp = avgj[(icntr) - fhalfwin : (icntr) + fhalfwin + 1]
                temp = signal.medfilt(temp, kernel_size=3)
                
                # Run the fit, with error checking
                try:
                    p0 = [1000, np.mean(xx), 3, bkgd]
                    aa, cv = optimize.curve_fit(gaussfit_func, xx, temp, p0=p0)
                except RuntimeError:
                    continue  # Just skip this one
                # print(f"Gaussfit parameters: {aa}")
                tempfit = gaussfit_func(xx, *aa)
                center = aa[1]
                fw = aa[2] * 1.177 * 2.
                pmax = aa[0]
                if fw > 1.0:
                    peaks.append(center)
                    fwhm.append(fw)
                j0 = jf + j

        centers = np.asarray(peaks)
        #fwhm = fwhm[1:]
    
    cc = np.where(np.logical_and(centers >=0, centers <=2100))
    centers = centers[cc]

    szc = len(centers)
    print(f" Number of lines: {szc}")

    print(f"At this point the code makes some plots.  Yippee.")
 
    """ 
    plot,avgj,xra=[0,ny+2],xsty=1,yra=[0,max(avgj)+0.2*max(avgj)], $
    title=title, xtitle='CCD column', ytitle='I (DN)'
    for id=0,szc(1)-1 do begin
        plots,centers(id),avgj(centers(id))+20,psym=1
        xyouts,centers(id),avgj(centers(id))+30,strtrim(centers(id),2), $
        /data,orientation=0.
    endfor

    """
 
    return (centers, fwhm)


def dfitlines(spectrum, cnt, cplot=False):
    """Procedure to fit line profiles in a focus image

    """
    
    nc = len(cnt)
    fwhm = np.empty(nc, dtype=float)
    boxi = 17
    boxf = 11
    xx = np.arange(2 * boxf + 1).flatten()

    # Loop through the lines!
    for j in range(nc):
        if (cnt[j] - boxi) < 0 or (cnt[j] + boxi) > 2030:
            continue
        box = spectrum[:,int(cnt[j]) - boxi : int(cnt[j]) + boxi + 1]
        ccnt = com1d(box) + cnt[j] - boxi
        if (ccnt - boxi) < 0 or (ccnt + boxi) > 2030:
            continue
        box = spectrum[:,int(ccnt) - boxf : int(ccnt) + boxf + 1].flatten()
        
        # Run the fit, with error checking
        try:
            p0 = [1000, np.mean(xx), 3, 0]
            a, cv = optimize.curve_fit(gaussfit_func, xx, box, p0=p0)
        except RuntimeError:
            continue  # Just skip this one
        # print(f"Gaussfit parameters: {aa}")
        tmp = gaussfit_func(xx, *a)
        fwhm[j] = a[2] * 2.355     # FWHM = 2.355 * sigma
    
    return fwhm


def com1d(line):
    """Returns the center-of-mass of line
    """

    ltz = np.where(line < 0)
    if len(ltz) != 0:
        line[ltz] = 0
    
    nx, ny = line.shape
    yy = np.arange(ny)

    cy = yy * line
    dy = np.sum(cy) / np.sum(line)

    return dy


def dfitcurves(fwhm, fnom=2.7):
    """Fit line / star focus curves
    """
    
    thresh = 2.0
    norder = 2

    nf, nc = fwhm.shape
    fits = np.empty((nc, norder+1), dtype=float)
    minfoc = np.empty(nc, dtype=float)
    x = np.arange(nf, dtype=float)
    focus = np.empty(nc, dtype=float)
    fwidth = np.empty(nc, dtype=float)
    xfine = np.arange((nf-1)*10.0 + 1.0, dtype=float)/10.0
    wts = np.full(nf, 1.0, dtype=float)

    for i in range(nc):
        data = fwhm[:,i]
        out = np.where(np.logical_or(data < 1.0, data > 15.0))
        if count := len(out) > 3:
            continue
        if count != 0:
            wts[out] = 0.0
            data[out] = 50.0
        fit = np.polyfit(x, data, norder)
        # The numpy function returns coeff's in opposite order to IDL func
        fits[i,:] = np.flip(fit)

        # Curve Miniumum
        fitfine = np.polyval(fit, xfine)
        minfine = np.min(fitfine)
        focus[i] = xfine[np.argmin(fitfine)]
        minfoc[i] = minfine

        # Nominal focus position (the larger of the two roots):
        fit[2] -= fnom
        fpix = np.roots(fit)
        fwidth[i] = np.max(fpix)

    return (fits, (focus, fwidth, minfoc))
