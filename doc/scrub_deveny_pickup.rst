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

   - ``scrub_deveny_pickup``: The DeVeny spectrograph suffers from EM ground-loop
     pickup noise in the readout electronics.  While the Instrument Groups is working
     diligently to identify the source of this noise, it is still affecting frames at
     the 5-10 DN level.  This routine will attempt to clean the pattern from the raw
     frames.  It works in an iterative fashion with the spectroscopic data reduciton
     pipeline `PypeIt <https://pypeit.readthedocs.io/en/release/index.html>`_ to
     identify and extract sources and sky signal before attempting to fit and remove
     the remaining pickup noise.  If you wish to use this routine, please follow the
     instructions in :ref:`optional-dependencies` to include PypeIt in your
     LDTObserverTools installation.


Usage
=====

The output of this tool is a multiextension FITS file that contains the cleaned
raw image, the original raw image, the extracted pattern (with both a zero mean
and with the raw image mean), and a table of the row-by-row sinusoidal fit
coefficients.  This cleaned image may be processed with PypeIt as usual (the
parameters for DeVeny data within PypeIt recognize this post-processed image
format).