============================
DeVeny Pickup Noise Scrubber
============================

Status: *In development*

Overview
========

The spectral channel CCD has nasty pickup noise.  There is a ground-loop
somewhere on the cube (or in the dewar) that is picking up EMI from the various
equipment on the cube.  The instrument team is working diligently to find the
source of and ameliorate the effects of this pickup.

In the meantime, this tool is designed to identify and fit out the sinusoidal
signal in each line of the CCD image.  It works in an iterative fashion,
requiring the raw data frames to first be processed through the PypeIt
spectroscopic data reduction pipeline before being processed here.

PypeIt fits out the sky background and objects in the 2d spectrum, leaving a
residual image, which contains the sinusoidal noise.  This tool uses this image
to fit the noise line-by-line.

Usage
=====

The output of this tool is a multiextension FITS file that contains the cleaned
raw image, the original raw image, the extracted pattern (with both a zero mean
and with the raw image mean), and a table of the row-by-row sinusoidal fit
coefficients.  This cleaned image may be processed with PypeIt as usual (the
parameters for DeVeny data within PypeIt recognize this post-processed image
format).