==================================
DeVeny Collimator Focus Calculator
==================================

Overview
========

   - ``dfocus``: Compute the needed collimator focus based on a series of arc line
     frames taken at various collimator settings.  Read in the arc lamp frames in
     the current night's focus directory, find the appropriate spectral lines in each
     frame, compute the FHWM (or other measure) of those lines, plot the FHWM as a
     function of collimator position and suggest the optimal focus position.  This
     program is executed identically to the old IDL version.  The python version
     uses :obj:`scipy.signal` processing routines for identifying line peaks and widths,
     resulting in more accurate and consistent estimates of the correct collimator
     focus position.  Rather than separately producing plots to the screen and disk,
     this version writes all plots to a PDF, then launches ``Preview.app`` to display
     the plots on the screen.  Newly added is a readout of the mount temperature so
     the user can determine when/if the collimator needs to be refocused during the
     night.  Online help is available with the ``-h`` option.  [`Completed: 2021-11-19`]


Usage
=====

Yeah, use it!