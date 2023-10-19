============================
DeVeny Pickup Noise Scrubber
============================

Status: *In development*

Overview
========

The spectral channel CCD of the DeVeny Optical Spectrograph has been subject
to ground-loop EMI pickup during most of its tenure on the Lowell Discovery
Telescope (LDT).  This manifests in the science images as sinusoidal pattern
noise added in the readout electronics.  In January 2018, a significant ground
loop was broken, eliminating the formerly prominent "corduroy" noise pattern.
There still exists, however, a low-level (:math:`\pm 5`) DN pickup signal
caused by a still-unidentified ground loop on the instrument cube.  While the
Instrument Group works diligently to identify and remove this remaining ground
loop, the present tool is available for observers whose low-SNR extractions are
affected by this residual signal.

.. _raw_frame:
.. figure:: figures/scrubber_raw_frame.*
   :class: with-shadow
   :alt: Raw DeVeny frame with sinusoidal pickup noise

   -- Raw frame from DeVeny.  This is the frame for which we will
   show the scrubber processing steps in the sections below.  The object is
   the thin horizontal line in the center of the frame, and there are several
   prominent (vertical) night sky lines in addition to the wavy pattern caused
   by AC EMI pickup in the readout electronics of the CCD.

This tool relies on an iterative data reduction approach.  The raw frame is
first processed by a conventional data reduction pipeline (`PypeIt
<https://pypeit.readthedocs.io/en/release/index.html>`_) to extract sky lines
and object spectra, the sinusoidal pickup noise is fit to the residual image
(raw frame minus sky and objects), the resulting pattern is subtracted from the
raw frame, and finally the scrubbed image is reprocessed with PypeIt to extract
the final 1D spectra for analysis.  Because the pickup noise is only a few DN
in amplitude, it is essential to remove bright features (cosmic rays, sky
lines, etc.) before attempting to fit this pattern noise.

Jumping to the punch line, :numref:`spec1d_comps` shows a comparison of the
extracted 1D spectra from both the raw scrubbed images for the raw frame shown
in :numref:`raw_frame` to illustrate the utility of this tool and the extent to
which it can make useful extant DeVeny data that is subject to this ground-loop
EMI pickup.

.. _spec1d_comps:
.. figure:: figures/spec1d_scrub_compare.*
   :class: with-shadow
   :alt: Comparison of raw and scrubbed extracted 1D spectra

   -- Comparison of the PypeIt-extracted 1D spectra from the raw and scrubbed
   frames to illustrate the utility of this tool.



.. _scrub_process:

Data Processing Steps
=====================

1. Before the pickup noise may be fit out, it is necessary to perform a basic
   spectroscopic reduction using the PypeIt pipeline.  The outline of how to
   perform this reduction is given in the `DeVeny User Manual
   <https://confluence.lowell.edu/display/LDTOI/DeVeny+Optical+Spectrograph>`_.
   In particular, for this first pass, it is important to not modify any of the
   object finding or extraction parameters in your PypeIt Reduction File, as
   the scrubber expects certain qualities of the processed images.  The user
   should, however, include the following in the Parameter Block to avoid
   artifacts being introduced into the resulting pattern file from the effects
   of local sky subtraction:[1]_


    .. code-block:: ini

      [reduce]
         [[skysub]]
            local_maskwidth = 200.
            no_poly = True

   .. warning::

      If you are using a version of PypeIt ``< 1.14.1``, then you will
      instead need to add the following to the Parameter Block to ensure the
      traced slit edges do not shrink unreasonably.  Regions outside the marked
      slits will have 0 value in the residual image, and therefore any
      sinusoidal signal there will not be fit out.  (The additional parameters
      included here were added to the DeVeny default set in PypeIt version
      ``1.14.1``.)

      .. code-block:: ini

         [reduce]
            [[skysub]]
               local_maskwidth = 200.
               no_poly = True
            [[findobj]]
               find_trim_edge = 0,0
         [calibrations]
            [[slitedges]]
               minimum_slit_length = 170
               max_nudge = 5
               det_buffer = 0
            [[flatfield]]
               tweak_slits = False



   .. important::

      The scrubber expects the reduced files to be in a particular location
      with respect to the raw files, namely ``<RAW_DIR>/ldt_deveny_?/Science/``.
      It is the author's intent to introduce an option to this tool that allows
      the specification of arbitrary reduced file directories, but that has not
      yet been implemented.

2. Once the PypeIt reduction is complete, you are ready to run the scrubber.
   In the directory containing the raw files, run the scrubber (see
   :ref:`usage`) on the files in the directory.  If an input file is **NOT**
   an object science frame (*e.g.*, a bias or dome flat), it will be skipped
   without being processed.  The scrubber will use both the raw frame and the
   processed 2D spectrum (in the ``ldt_deveny_?/Science`` directory) to fit
   the sinusoidal pickup noise and produce a cleaned image.

3. After the scrubber has completed its work, you will need to re-run the
   ``pypeit_setup`` script to include the scrubbed files in a new PypeIt
   Reduction File.  At this point, you may remove the lines in the ``.pypeit``
   file corresponding to the raw images and add any parameter modifications
   you desire, including those affecting object finding and extraction.


.. [1] We prefer this method over simply turning off local sky subtraction
       because in that case no object model is written to the ``spec2d``
       file by PypeIt, and it is oftentimes important to remove the object
       from the frame before attempting to fit the sinusoidal pickup noise. 


.. _usage:

Usage
=====

Command Line
^^^^^^^^^^^^

The tool usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/scrub_deveny_pickup.rst

The ``--proc_dir`` option may be used to specify the location of the
PypeIt-processed files, if other than as a subdirectory of the raw data
directory.  This works in the same sense but opposite direction as the ``-r
<RAWDIR>`` option to ``pypeit_setup``.

To simply overwrite the raw data file with the output of this tool, the
``--overwrite_raw`` option may be specified.



Terminal Output
^^^^^^^^^^^^^^^

When run, the tool prints to the screen information about the current
processing step and produces progress bars when repeatedly running an algorithm
on multiple rows of the image.  The terminal output for the example frame used
in this documentation is shown below.

.. code-block:: console

   $ scrub_deveny_pickup 20290101.0001.fits

   Processing frame 20290101.0001.fits
      --> FFT-predicted pixel period: 169.1 pix
   Checking the object model for bulk sinusoidal signal...
   100%|██████████████████████████████████████████| 54/54 [00:01<00:00, 31.39row/s]
   * Adding the object model back into the residual image for fitting; extracted
     object contains target sinusoidal signal, amplitude = 5.20.
   Fitting sinusoids to each line in the image...
   100%|████████████████████████████████████████| 508/508 [00:10<00:00, 50.65row/s]
      --> Mean fit pixel period: 172.7 pix
   Refitting lines with poor fit in first pass...
   100%|██████████████████████████████████████████| 74/74 [00:01<00:00, 66.74row/s]
      --> Mean fit pixel period: 173.6 pix
    Writing QA plots to ./ldt_deveny_A/QA/PDFs
   Writing out scrubbed FITS file: 20290101.0001_scrub.fits

Important pieces to note in the terminal output:

* The FFT-predicted and mean sinusoid fit pixel periods are printed for quality
  assurance purposes.

* If the PypeIt-reduced ``spec2d`` file contains a non-zero object model, the
  code will analyze it for the sinusoidal signal at the predicted period.  If
  the object model contains this signal, then the residual image (``sciimg`` - 
  ``skymodel`` - ``objmodel``) will be missing power from the pickup noise
  (*i.e.*, that power was extracted into the object model).  In this case, the
  tool will print a statement indicating the object model will be included in
  the residual image, and the amplitude of the sinusoidal fit (for QA
  purposes).  If the object(s) remains extracted (*i.e.* does not contain bulk
  sinusoidal signal), the tool will print the message: ``Object model appears
  clean of target sinusoidal signal``.

* The first pass at fitting sinusoids to the image should include all rows in
  the (trimmed) image (FITS keyword ``TRIMSEC``).  In this case, the CCD was
  binned ``1x1``, so there were 516 rows in the original image, trimmed down to
  508.

* As an iterative step, the tool assesses the rms residuals of the fit and
  determines if there are outliers that need to be refit using initial guess
  parameters from nearby well-fit rows.  The number of rows indicated in this
  progress bar is less than the total image, and can vary based on the
  particular pattern noise in the image.

* Finally, the location of QA plots and filename of the scrubbed FITS file are
  shown for reference.


Output FITS File
^^^^^^^^^^^^^^^^

The output of this tool is a multiextension FITS with filename identical to the
raw file (including path) except with "_scrub" appended before ".fits".  It
includes the following HDUs:

  * Primary HDU: contains the raw file's complete FITS header
  * Cleaned Image HDU: (``raw`` - ``pattern``)
  * Original Image HDU: (``raw``)
  * Pattern Image HDU: (``pattern``)
  * Pattern about Raw Mean Image HDU: (``pattern`` + ``mean(raw)``)
  * Fit Coefficients BinTable HDU: (``fit_coeffs``)

The AstroPy utility ``fitsinfo`` yields the following:

.. code-block:: console

   $ fitsinfo 20290101.0001_scrub.fits
   Filename: 20290101.0001_scrub.fits
   No.    Name      Ver    Type      Cards   Dimensions   Format
     0  PRIMARY       1 PrimaryHDU     193   ()
     1  CLEANED       1 ImageHDU        10   (2148, 516)   float64
     2  ORIGINAL      1 ImageHDU        12   (2148, 516)   int16 (rescales to uint16)
     3  PATTERN0      1 ImageHDU        10   (2148, 516)   float64
     4  PATTERN1      1 ImageHDU        10   (2148, 516)   float64
     5  FIT DATA      1 BinTableHDU     50   508R x 15C   [D, D, D, D, D, D, D, D, D, D, D, D, D, D, D]

The ``ORIGINAL`` raw data frame is included in the scrubbed output for reference
and posterity.  The fitted sinusoidal pattern is included in two different
formats, both with a zero mean (``PATTERN0``) and the raw image mean
(``PATTERN1``) for use as desired (*e.g.* it is easier to compare ``PATTERN1``
to the raw image in ``ds9``, but the actual signal removed from the raw image
is ``PATTERN0``).  Finally, the sinusoidal fit coefficients are included for
perusal.

Of note, PypeIt will recognize the scrubbed image and use the ``CLEANED`` HDU
for processing without user intervention.


.. _scrub_details:

Algorithmic Details of the Scrubbing
====================================

Here we describe the processing details for the scrubber and illustrate example
QA plots generated by the tool.  All QA plots are placed into the PypeIt
``QA/`` directory for convenience.[2]_

The fitting of the sinusoidal pickup noise is done on the residual
PypeIt-processed image, in which cosmic rays, the sky model, and (frequently)
the object have been removed, leaving behind the underlying noise in the image.
This is done because bright features in the raw science frame such as night sky
lines and objects obfuscate and cover the AC signal we wish to remove.  We take
advantage of the sophisticated modeling algorithms in PypeIt to do this heavy
lifting for us, leaving us to concentrate on the DeVeny-specific issues at
hand.  See :numref:`image_comparisons` below for an illustration of the various
processing steps and frames used.

.. note::


   Include a description of when the object may or may not be included in the
   residual image -- namely if the object model from PypeIt contains
   sinusoidal signal of the period indicated from the FFT.


FFT Analysis
^^^^^^^^^^^^

The EMI pickup noise present in the DeVeny images has an unpredictably variable
apparent frequency (or wavelength in pixel space) from frame to frame.  The
secular variation in the frequency is slow enough, however, that the measured
wavelength of the sinusoidal noise is approximately constant across a single
image given the :math:`\sim 8` second readout time of the spectral channel
CCD.

Estimation of the mean sinusoid period is done with a fast fourier transform
(FFT).  Rather than attempt to pull the (varying) horizontal frequencies out
of a 2D FFT, the image is instead flattened into a timeseries-like array, where
pixels are arranged in the order in which they were read out.  In its
processing, PypeIt trims the raw image to remove overscan sections and unruly
rows at the edges of the CCD.  As such, the last pixel of one row is not
strictly temporally adjacent to the first pixel of the next row in the same way
that pixels within a given row are.  Rather than attempt to fill in a
temporally appropriate gap between rows (which would introduce artifacts in the
resulting FFT), we simply stitch the rows together as-is and rely upon the FFT
to pick out the prominent frequencies in the flattened array.

The QA plot for the FFT analysis is shown in :numref:`fft_analysis`.

.. _fft_analysis:
.. figure:: figures/scrubber_fft_analysis.*
   :class: with-shadow
   :alt: FFT analysis of the flattened frame

   -- FFT QA plot, showing the flattened image array, the real
   (amplitude as a function of frequency) and imaginary (phase as a function
   of frequency) components of the FFT, along with the power spectrum
   (absolute value squared).  The power spectrum is further smoothed with a
   gaussian kernel to help isolate real signal (with variable frequency, as
   discussed above) from artifacts and ringing in the FFT.  The peak at 169.1
   pixels, indicating the frequency with the most power in the flattened array:
   most likely the period of the AC EMI pickup noise.

The power spectrum of the FFT (absolute value squared) is smoothed with a 10-Hz
wide gaussian kernel since there is variation in the frequency of the sinusoid
from row to row.  This has the added benefit of squashing artificial peaks in
the power spectrum caused by ringing or aliasing of the signal across rows.
We use the pixel period of the peak of the power spectrum as the initial guess
for the iterative fitting of a sinusoid (in pixel space) to each row of the
residual image.

Row-by-Row Sinusoid Fits
^^^^^^^^^^^^^^^^^^^^^^^^

:numref:`sinusoid_fits` shows the result of iteratively fitting a sinusoid to
each row of the PypeIt residual image for this example frame.  In addition to a
basic sinusoid (amplitude, period, and phase) we include a quadratic polynomial
to account for secular drift in the background.  The
:func:`~scipy.optimize.curve_fit` function from :mod:`scipy.optimize` is used
to perform a non-linear least-squares fit to each row, with bounds placed on
the sinusoid terms to keep the final fit reasonably close to initial guess
values.

The rms of the fit is computed for each row as an estimate of how well any
particular row's model matches the data.  The mean and standard deviation of
the row-by-row rms values (bottom panel in :numref:`sinusoid_fits`) are used to
identify outliers that likely have a poor fit.  These lines are refit using a
the fit values from nearby "good" rows as the initial guesses and tighter
bounds on the fit values.  As a result of this iterative process, the set of
sinusoidal fits tends to have a narrow range of values from line to line and
produces a pattern image (see :numref:`image_comparisons`) that closely
resembles the unwanted pattern in the raw frame.  Note, however, that the final
pixel periods tend to be slightly larger than that predicted from the FFT
(green dashed line).

.. _sinusoid_fits:
.. figure:: figures/scrubber_sinusoid_fits.*
   :class: with-shadow
   :alt: Sinusoid fits to each row

   -- Sinusoid fit QA plots, showing the row-by-row fit
   parameters for the amplitude, period, and phase in the first three panels,
   and the rms of the fit in the bottom panel.  In the period panel, the
   green dashed line indicates the peak of the FFT power spectrum used as the
   initial guess for the sinusoid fit parameters.


Pickup Noise Pattern Construction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Finally, we use the sinusoid fits to produce a pattern image.  This is the
zero-mean sinusoid (sans quadratic polynomial) that *should* represent only the
EMI pickup noise (as an additive AC-only signal).  Subtract that and make
additional QA plots to show how totally awesome this is!!!


.. _image_comparisons:
.. figure:: figures/scrubber_image_comparisons.*
   :class: with-shadow
   :alt: Image comparison plot showing initial, intermediate, and final images

   -- This QA plot illustrates the image-space effects of the sinusoidal fits.
   Panels are all shown with the IRAF ZScale mapping (maybe include link?).
   The panels are described below.

Panel Description:

1. The science image after the initial PypeIt processing.  Without scrubbing,
   this would be the 2D image from which the object(s) are extracted.

2. The PypeIt sky model, where the local sky modeling has included the entire
   slit thanks to the ``local_maskwidth = 200.`` specified in the PypeIt
   Reduction File.

3. The PypeIt residual image, where the sky model (and sometimes the object
   model) has been subtracted and cosmic rays (and other bad pixels) have been
   masked.  What *should* remain in this image is flat residual noise with rms
   based on the sky model photon statistics.  In the case of DeVeny images
   subject to the EMI pickup noise, that noise should be the dominant feature
   in this frame.

4. The modeled sinusoidal pickup noise.  While the sinusoidal fit to each row
   of the residual image (#3) includes a quadratic polynomial to account for
   secular variation in the background, only the sinusoid is included in this
   pattern image due to the AC nature of the pickup.

5. The scrubbed residual image, where the pattern image (#4) has been
   subtracted off.  This frame should be largely featureless (except if the
   object is still included).

6. This is the scrubbed science frame (#1 minus #4), and it should be clear of
   the sinusoidal noise.  Even in frames where the phase difference of the
   pickup noise from row to row is such that large-scale wavy patterns are not
   visible in the base science image (#1), this scrubbed frame should be less
   noisy since a :math:`\sim 5` DN sinusoid has been removed from each line.

Panels #1, #2, and #6 are shown with common scale limits based on the processed
science image, and panels #3 - #5 are shown with common scale limits based on
the residual image.



.. [2] Within the PypeIt ``QA/`` directory, that pipeline places its QA plots
       into a ``PNG/`` subdirectory (because the plots are PNG format).  To
       keep QA files from this tool separate (yet also in an easy-to-find
       location), QA plots generated here are in the ``PDF\`` subdirectory
       (because, you guessed it, they are in PDF format).


--------

.. _appendix:

Appendix
========

For the adventurous reader, here is some additional information about the
development of the scrubber and various points of interest.

.. _fourier_analysis:

Line-by-Line Fits vs Fourier-Space Filtering
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Given the nearly constant frequency (manifested as a periodic sinusoid in
pixel space) of the AC EMI, one is tempted to attempt filtering the signal in
Fourier space rather than doing a line-by-line fitting.

The FFT computation is done on a flattened 1D version of the image, with the
pixels in time-order of readout.  Now, if the exact time between rows were
included as spacers between the rows of data in the 1D array, a single period
may emerge.  Since this is difficult to guess (and can vary within the readout
electronics), the rows were simply strung together one after another.  The
effect of this is a significant amount of ringing in the FFT caused by
breaks in the sinusoidal pattern from row to row.

We can run the computed EMI pattern image (4th panel in
:numref:`image_comparisons`) through the same FFT analysis routine to yield the
QA plot shown in :numref:`pattern_fft`.

.. _pattern_fft:
.. figure:: figures/pattern_scrubber_fft_analysis.*
   :class: with-shadow
   :alt: FFT analysis of the constructed pattern

   -- FFT QA plot, showing the flattened image array, the real
   (amplitude as a function of frequency) and imaginary (phase as a function
   of frequency) components of the FFT, along with the power spectrum
   (absolute value squared).  The power spectrum is further smoothed with a
   gaussian kernel to help isolate real signal (with variable frequency, as
   discussed above) from artifacts and ringing in the FFT.  The peak at 169.1
   pixels, indicating the frequency with the most power in the flattened array:
   most likely the period of the AC EMI pickup noise.

When compared to :numref:`fft_analysis`, most of the structure apparent in the
analysis of the residual image still appears here.  This is why simply applying
a notch filter in Fourier space around the peak of the FFT does not clean the
EMI noise, but rather amplifies the problem.  It is comforting that the peak
period found in the pattern image is identical to that in the residual image
(169.1 pixels), although both are smaller than the average pixel period found
in the sinusoid fitting (173.6 pixels).

Finally, we can perform the same analysis on the scrubbed residual image (5th
panel in :numref:`image_comparisons`).  The result is shown in
:numref:`cleaned_fft`.

.. _cleaned_fft:
.. figure:: figures/cleaned_scrubber_fft_analysis.*
   :class: with-shadow
   :alt: FFT analysis of the cleaned residual image

   -- FFT QA plot, showing the flattened image array, the real
   (amplitude as a function of frequency) and imaginary (phase as a function
   of frequency) components of the FFT, along with the power spectrum
   (absolute value squared).  The power spectrum is further smoothed with a
   gaussian kernel to help isolate real signal (with variable frequency, as
   discussed above) from artifacts and ringing in the FFT.  The peak at 169.1
   pixels, indicating the frequency with the most power in the flattened array:
   most likely the period of the AC EMI pickup noise.

.. _sky_wobbles:

Introduced Structure in the Sky Model Caused by EMI
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Sometimes, the apparent phase of the sinusoidal EMI lines up such that PypeIt
includes some of that signal in its sky model.  For instance,
:numref:`skysub_1stpass` shows the Image Comparison Plots for a different raw
Deveny frame.  Visual inspection indicates that some of the sinusoidal signal
is attributed to the sky, yielding a not completely uniform EMI pattern in the
residual image.

.. note::

   This is ``20221130.0114.fits``.  Just in case I need to redo the images as
   I improve the underlying routine.

.. _skysub_1stpass:
.. figure:: figures/skysub_1stpass_image_comparisons.*
   :class: with-shadow
   :alt: First-pass scrubbing image comparisons

   -- Image comparison plots (like those in :numref:`image_comparisons`) for a
   different raw DeVeny frame.  Notice the oscillating sky model (2nd panel)
   that is more pronounced that that in :numref:`image_comparisons`.

The cleaned science image (bottom panel in :numref:`skysub_1stpass`) shows the
effect of this oscillating sky pattern because it was removed from the raw
frame before the sinusoidal fitting was done.

To assess how much of a problem this is, we can run the scrubbed file back
through this tool to produce a ``_scrub_scrub.fits`` frame.  The Image
Comparison plots for the 2nd pass through the scrubber are shown in
:numref:`skysub_2ndpass`.

.. _skysub_2ndpass:
.. figure:: figures/skysub_2ndpass_image_comparisons.*
   :class: with-shadow
   :alt: Second-pass scrubbing image comparisons

   -- Image comparison plots for the second pass through this tool of the frame
   shown in :numref:`skysub_1stpass`.  Note that the sky model appears to be
   nearly identical to that from the first pass, and no hint of the oscillating
   sky remains in the residual image (3rd panel).

Since the introduced oscillation is once again pulled out in the sky model, the
residual image only shows the remaining horizontal streaking seen in the bottom
panel of :numref:`skysub_1stpass`, which is also an artifact of the sky model
oscillations.  *Mumble, mumble,* something about flux calibration, but shape of
the spectrum is unaffected.
