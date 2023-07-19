.. |License| image:: https://img.shields.io/github/license/LowellObservatory/LDTObserverTools
   :target: https://github.com/LowellObservatory/LDTObserverTools/blob/main/LICENSE

.. |astropy| image:: http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat
    :target: http://www.astropy.org/

.. |forks| image:: https://img.shields.io/github/forks/LowellObservatory/LDTObserverTools?style=social
   :target: https://github.com/LowellObservatory/LDTObserverTools

.. |stars| image:: https://img.shields.io/github/stars/LowellObservatory/LDTObserverTools?style=social
   :target: https://github.com/LowellObservatory/LDTObserverTools

.. |watch| image:: https://img.shields.io/github/watchers/LowellObservatory/LDTObserverTools?style=social
   :target: https://github.com/LowellObservatory/LDTObserverTools

.. |github| image:: https://img.shields.io/badge/GitHub-LDTObserverTools-brightgreen
   :target: https://github.com/LowellObservatory/LDTObserverTools

.. image:: https://raw.githubusercontent.com/LowellObservatory/LDTObserverTools/main/doc/_static/obstools_logo.png
    :target: https://github.com/LowellObservatory/LDTObserverTools
    :width: 500

.. include:: include/links.rst

.. _obstools_main:

LDTObserverTools |forks| |stars| |watch|
========================================

|github| |astropy| |License|


**Version**: |version|

LDTObserverTools is a collection of command-line and GUI tools for observers at
the Lowell Discovery Telescope (LDT) in Happy Jack, AZ.

.. contents:: Table of Contents
    :depth: 1
    :local:

----

.. _programs:

List of Programs
================

.. _extant:

Tools Contained in this Package:
--------------------------------

   - ``deveny_grangle``: Compute the desired grating angle based on selected
     grating and desired central wavelength.  This routine comes with two interfaces.
     The default GUI features a dropdown menu for grating selection and contains error
     checking on the input for central wavelength.  There is a ``MAX_GUI`` option for
     computing central wavelength given the grating angle in addition to the standard
     GUI features.  Also included is a command-line interface, identical to the old
     IDL function.  Online help is available with the ``-h`` option.
     [`Completed: 2021-01-26`]

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

   - ``fix_ldt_header``: Replace / add / update FITS keyword values.  While named for
     the LDT, it can beused with any FITS file (or list of files).  The inspiration for
     this script is the cases when the LDT / LOUI does not input the proper information
     in FITS headers (`e.g.`, ``GRATING = UNKNOWN`` for DeVeny before the drop-down
     menu has been selected).  Online help is available with the ``-h`` option.
     [`Completed: 2022-10-17`]

.. _future:

Future Tools (planned or in development):
-----------------------------------------

   - ``lmi_etc``: Python GUI version of the LMI exposure time calculator
     (http://www2.lowell.edu/users/massey/LMI/etc_calc.php).

   - ``neocp_ephem``: An ephemeris-generating tool for querying the JPL Scout database
     for short-shelf- life NEOs that have not yet been assigned a Horizons identifier.
     This tool will turn the returned data into a file that can be ingested into the
     LDT TCS for observations.

   - ``deveny_collfoc_range``: Use the specified grating angle and mount temperature
     to suggest a range for use with the DeVeny LOUI collimator focus sequence
     function.  This is important because, unlike all other focus routines at LDT,
     this function takes the **starting point**, step, and number of exoposures
     rather than the **expected focus value**, step, and number of exposures.  Having
     a routine to compute the desired values would make this step easier and less
     prone to error (`i.e.`, only searching on one side of the expected focus).

   - ``validate_input_list``: The extant input list validation tool
     (https://confluence.lowell.edu/display/LDTOI/Validate+Input+List) was
     produced in 2015 using an old Java Runtime Environment that is not available on
     most modern operating systems.  As such, the tool is virtually useless.  This
     python program would provide this key functionality in a modern environment.

   - ``observer_target_list``: The extant observer target list tool
     (https://confluence.lowell.edu/display/LDTOI/Observer+Target+List+User+Manual)
     was produced in 2015 using an old Java Runtime Environment that is not available on
     most modern operating systems.  It still runs happily on the LDT observer machines,
     but those will need to be replaced at some point.  This python version will provide
     a future-proof solution.
   
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

==================

.. _installing:

Installing
==========

.. _environment:

Set up a clean python environment
---------------------------------

Because installing a python tool like LDTObserverTools also installs or upgrades
its :ref:`dependencies`, the best course of action is to set up a clean python
environment into which the installation will occur.  This mitigates any possible
dependency conflicts with other packages you use.

The recommended method of setting up a new envrionment is with ``conda``:

.. code-block:: console

    conda create -n obstools python=3.11
    conda activate obstools

See `Managing Environments with Conda
<https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`_
for more details.

.. _pip_install:

Installing LDTObserverTools via ``pip``
---------------------------------------

.. note::
  The commands here assume ``pip`` is associated with Python3.  To check, run
  ``pip --version`` from the command line and check that the associated python
  version is ``>= 3.9``.

  Also, you will need ``git`` installed on your system / in your environment.

To install the latest version of LDTObserverTools and its required dependencies,
execute:

.. code-block:: console

    pip install "obstools @ git+https://github.com/LowellObservatory/LDTObserverTools"

This will download the latest version of the package from GitHub and install it along
with its required dependencies.  (Note: whether or not you need quotation marks depends
on your particular shell -- ``bash``, ``zsh``, etc.)

As the package undergoes continued development, it will be necessary to upgrade
your installation to access the latest features.  The upgrade process should
simply be a matter of executing:

.. code-block:: console

    pip install "obstools @ git+https://github.com/LowellObservatory/LDTObserverTools" --upgrade

.. _optional-dependencies:

Optional Dependencies
^^^^^^^^^^^^^^^^^^^^^

Some of the instrument-specific routines in this package require additional dependencies
that are not otherwise needed by the majority of the routines herein.

   - If you are using the ``deveny_pickup_cleaner`` routine, you will need the
     spectroscopic data reduction pipeline PypeIt for the iterative cleaning of
     the pickup noise.  It can be installed by including it in the optional
     dependencies, `e.g.`:

      .. code-block:: console

        pip install "obstools[pypeit] @ git+https://github.com/LowellObservatory/LDTObserverTools"
