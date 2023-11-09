.. include:: include/links.rst

.. |nbsp| unicode:: 0xA0 
   :trim:

.. _fix_ldt_header:

==============================
Simple FITS Header Fixing Tool
==============================

Status: Completed 2022-10-17

Overview
========

   - ``fix_ldt_header``: Replace / add / update FITS keyword values.  While named for
     the LDT, it can beused with any FITS file (or list of files).  The inspiration for
     this script is the cases when the LDT / LOUI does not input the proper information
     in FITS headers (`e.g.`, ``GRATING = UNKNOWN`` for DeVeny before the drop-down
     menu has been selected).  Online help is available with the ``-h`` option.
     [`Completed: 2022-10-17`]


Usage
=====

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/fix_ldt_header.rst
