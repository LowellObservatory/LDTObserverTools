.. include:: include/links.rst

.. |nbsp| unicode:: 0xA0 
   :trim:

.. _ephem_generator:

=========================================
NEO Confirmation Page Ephemeris Generator
=========================================

Status: *In development*

Overview
========

   - ``ephem_generator``: An ephemeris-generating tool for querying the JPL Scout database
     for short-shelf- life NEOs that have not yet been assigned a Horizons identifier.
     This tool will turn the returned data into a file that can be ingested into the
     LDT TCS for observations.


Also include MPC and/or AstOrb (Lowell) and/or IMCCE (in France) as sources for this tool.


Usage
=====

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/ephem_generator.rst
