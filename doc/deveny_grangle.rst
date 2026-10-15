.. include:: include/links.rst

.. |nbsp| unicode:: 0xA0 
   :trim:

.. _deveny_grangle:

===============================
DeVeny Grating Angle Calculator
===============================

Status: Completed 2021-01-26

Overview
========

.. _grateq:

The Grating Equation
^^^^^^^^^^^^^^^^^^^^

The angles at which the diffracted light reflects off the grating are given
by the grating equation:

.. math::

    \begin{equation}
    m\lambda = d(\sin \theta _{i}-\sin \theta _{m})~,
    \end{equation}

where :math:`d` is the spacing between adjacent grooves on the grating,
:math:`\theta_{i} = \theta_{\rm grangle} + 10^{\circ}` is the incident angle,
and :math:`\theta_{m} = 55^{\circ} - \theta_{i}` is the outgoing angle for light
with wavelength :math:`\lambda` diffracting into order :math:`m` (see the
DeVeny user manual for a description of the physical layout of the
spectrograph).  The DeVeny gratings are typically operated in
1\ :superscript:`st`\  order (:math:`m = 1`), although use of :math:`m=2` would
possible with the proper short-pass order-blocking filters to remove
1\ :superscript:`st`\ -order light.

This tool computes the necessary tilt :math:`\theta_{\rm grangle}` by
numerically solving the equation below for a given grating.  The computation
uses the groove density (in g/mm) in place of :math:`d`, and finds
:math:`\theta_{\rm grangle}` for a specified 1\ :superscript:`st`\ -order
central wavelength (:math:`\lambda_c`) in Angstroms on the spectral CCD.

.. math::

    \begin{equation}
    \lambda_c~({\rm \mathring{A}}) = \frac{\sin(\theta_{\rm grangle} +10^{\circ}) - \sin(45^{\circ} - \theta_{\rm grangle})}{ {\rm groove~density~(g/mm)} }\times 10^7
    \end{equation}


Usage
=====

Run ``deveny_gui`` and select the **Grating Angle** tab. The existing
``deveny_grangle`` command opens the same tab; ``deveny_grangle --cli``
retains the command-line calculator.

Choose a grating, then enter a central wavelength from 3000 to 11000 Å and
click **Compute tilt**. The tab shows the grating tilt and slit demagnification.
To reverse the calculation, enter a tilt from 0 to 48 degrees and click
**Compute wavelength**.
