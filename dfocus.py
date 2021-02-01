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
import pathlib
import sys

# Numpy & SciPy & PySimpleGui
import numpy as np
from scipy import optimize
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

    global curves, fx, centers, fpos, fwid, fwmin

    # Parse the log file to obtain file list;
    files, focid = flogparse(flog)
    nfiles = len(files)

    # Pull the spectrograph setup from the first focus file:
    with fits.open(f'../' + files[0]) as hdu:
        hdr0 = hdu[0].header
    
    slitasec = hdr0['SLITASEC']
    grating = hdr0['GRATING']
    grangle = hdr0['GRANGLE']
    lampcal = hdr0['LAMPCAL']
    filtrear = hdr0['FILTREAR']
    fnom = 2.94 * slitasec * deveny_amag(grangle)

    # Pull the collimator focus values from the first and last files
    f0 = hdr0['COLLFOC']
    with fits.open(f'../' + files[-1]) as hdu:
        hdr1 = hdu[0].header
    f1 = hdr1['COLLFOC']

    # Set up arrays for populating the focus values
    df = (f1 - f0) / (nfiles - 1)
    x = np.arrange(nfiles, dtype=float)
    fx = np.arrange(nfiles, dtype=float) * df + f0

    # Examine the middle image
    mfile = '../' + files[int(nfiles/2)]
    swext=4  # Vestigial extraction pramater
    win = 11 # Something about the window to extract
    print(f'\n Processing object image {mfile}...')
    spectrum = dvtrim(mfile)
    # sz = size(spectrum)
    # print(sz)
    # traces = make_array(sz(1),value=sz(2)/2,/float)
    # dextract, spectrum, traces, win, spectra=mspectra, swext=swext





def flogparse(flog):
    """Parse the focus log file produced by the DeVeny LOUI

    flog: if 'last', then read in the most recent log file, but
            the user can also pass in a specific log file to use
    """
    if flog.lower() == 'last':
        path = pathlib.Path('.')
        focfiles = sorted(path.rglob('deveny_focus*'))
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
      """

    # Parameters for DeVeny (2015 Deep-Depletion Device):
    darkcol = 0
    nxpix = 2048
    prepix = 50
    postpix = 50

    with fits.open(filename) as hdu:
        header = hdu[0].header
        image = float(hdu[0].data)
    fimage = image

    return False

    """ image=image(*,12:511)
    sz=size(image)
    ny=sz(2)

    ; set parameters
    hotcol=[2007,2008,1128,1129]-1
    numamps = sxpar(header,'NUMAMP')


    gaina=1. & gainb=1.	; no gain amplification

    jcnt=sz(2)/2			; half length
    

    xx=findgen(sz(2))		; dummy array for bias fitting




    ; setting trim=1 will make a 1-amp size gain-corrected image out of a two amp
    if keyword_set(trim) then begin
    ixb=sz(1)-postpix-nxpix		; left amplifier
    ixa=sz(1)-postpix-nxpix/2		; right amplifier
    for iy=0,sz(2)-1 do begin
    image(ixb:ixb+nxpix/2-1,iy)=(image(ixb:ixb+nxpix/2-1,iy)+bfit(iy))
    image(ixa:ixa+nxpix/2-1,iy)=(image(ixa:ixa+nxpix/2-1,iy)+bfit(iy))
    endfor
    imagetmp=image(prepix/2:prepix/2+nxpix+prepix/2+postpix/2-1,*)
    return,imagetmp
    endif else begin
    imagetmp=image(prepix:prepix+nxpix-1,0:sz(2)-1,*)
    return,imagetmp
    endelse

    end """
