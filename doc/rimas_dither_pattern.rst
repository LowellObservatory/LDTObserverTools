.. include:: include/links.rst

.. |nbsp| unicode:: 0xA0
   :trim:

.. _rimas_dither_pattern:

================================
RIMAS Dither Metadata Correction
================================

Status: Completed 2026-08-20

Overview
========

The ``rimas_dither_pattern`` tool creates an ECSV sidecar table that corrects
RIMAS dither metadata.  Use it when a commanded ABBA or ON/OFF sequence is not
represented correctly in the FITS headers, including sequences in which one or
more exposures were commanded but not written to disk.

Each sequence must contain the full ordered list of commanded frame numbers.
The tool uses a frame's position relative to the first frame in its sequence to
assign ``dithnum`` and ``dithpos``.  It writes one correction row for each
selected RIMAS arm and frame number, including frames for which no FITS file was
saved.  By default, corrections are generated for both the YJ and HK arms.

The output table contains the following columns:

``filename``
    Expected RIMAS FITS filename constructed from the observing date, frame
    number, and arm.
``dithpat``
    Canonical dither pattern name, either ``ABBA`` or ``ONOFF``.
``dithnum``
    One-indexed position within the commanded sequence.
``dithpos``
    Dither position corresponding to ``dithnum``.
``dithseq``
    One-indexed sequence identifier, assigned in command-line order.
``source``
    Description of the source used to construct the corrections.

Usage
=====

Repeat ``-s`` or ``--sequence`` once for each commanded dither sequence.  Bash
brace expansion provides a compact way to pass every frame in a continuous
range.  For example, the shell expands each expression below into four frame
arguments before invoking the tool:

.. code-block:: console

    $ rimas_dither_pattern ABBA --date 20260305 \
        -s {86..89} -s {90..93} -s {94..97}
    Wrote 24 corrections to /current/path/rimas_dither_corrections.ecsv

This example creates corrections for three ABBA sequences and both detector
arms.  Include every originally commanded frame in each range, even when its
FITS file was not saved.  Comma-separated lists are also accepted when brace
expansion is unavailable, for example ``-s 86,87,88,89``.

Use ``--arms YJ`` or ``--arms HK`` to restrict output to one arm, ``--source``
to describe where the corrected sequence information came from, and
``--overwrite`` to replace an existing output file.

The complete command-line reference is shown below:

.. include:: help/rimas_dither_pattern.rst
