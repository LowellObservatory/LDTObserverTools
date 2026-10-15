.. include:: include/links.rst

.. |nbsp| unicode:: 0xA0 
   :trim:

.. _deveny_collfocus:

=================================
DeVeny Collimator Focus Estimator
=================================

Status: Completed 2023-07-31

Overview
========

The internal optics of the DeVeny Spectrograph are focused by pistoning the
collimator mirror (a separate task from focusing the telescope onto the 
spectrograph).  The DeVeny LOUI includes a tab for performing a focus sequence
(images of arc lines interleaved with collimator focus moves).  The purpose
of this GUI tool is for estimating the values you should input into the LOUI.

For a discussion of the physics behind the collimator focus procedure
(including the origin of the equation below), see the DeVeny User Manual.  An
empirical analysis of collimator focus values from data taken in 2017-2020
yield the approximate relationship:

.. math::

    \begin{aligned}
    \mathrm{collimator~focus}~\approx&~11.0 - 0.08~T_{\rm m} \mathrm{(^{\circ} C)} - 0.14~(\mathrm{tilt} - 25^{\circ}) + 0.7~\delta_{\rm filter}~~~\mathrm{mm} \\
    \delta_{\rm filter} &=  \begin{cases} \mbox{0,} & \mbox{if no order-blocking filter (REARFILT = Clear (C))}  \\ \mbox{1,} & \mbox{if using an order-blocking filter (e.g., GG420)} \end{cases}
    \end{aligned}

where :math:`T_{\rm m}` is the mount temperature and tilt is the grating tilt
angle.  The RMS of the fit is :math:`\sim 0.6` |nbsp| mm.  This equation is coded into
the GUI tool.

Unlike other focus procedures at LDT, the DeVeny LOUI focus sequence script
does not take the *expected* focus value, but rather one end of the range,
along with the step size and number of steps.  This tool is designed to assist
you in constructing a range that centers on the value computed from the
equation above to most efficiently complete the focus procedure.

After the images have been taken, use the :doc:`dfocus` program to analyze the
images and compute the appropriate focus position.

Usage
=====

Run ``deveny_gui`` and select the **Collimator Focus** tab. The existing
``deveny_collfocus`` command opens the same tab.

Enter the mount temperature, grating tilt, and rear filter, then click
**Compute focus sequence**. The tab shows the estimated focus (with 0.6 mm
uncertainty) and the initial position, step size, and number of steps for
LOUI. When ActiveMQ is configured at the telescope, **Get values from broker**
loads the current temperature and tilt; manual entry is always available.

The lower limit on the collimator focus motorized stage is 7.75 mm.
