# -*- coding: utf-8 -*-
#
#  This file is part of LDTObserverTools.
#
#   This Source Code Form is subject to the terms of the Mozilla Public
#   License, v. 2.0. If a copy of the MPL was not distributed with this
#   file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
#  Created on 17-Oct-2022
#
#  @author: tbowers

"""FITS File Utility Module

LDTObserverTools contains python ports of various LDT Observer Tools

Lowell Discovery Telescope (Lowell Observatory: Flagstaff, AZ)
https://lowell.edu

This file contains FITS Header utility routines.
"""

# Built-In Libraries
import argparse
import pathlib

# 3rd-Party Libraries
import ccdproc

# Local Libraries

# CONSTANTS


def fix_ldt_header(files, keyword: str, new_value):
    """Change FITS header keywords

    Sometimes at the telescope, incorrect or incomplete information is placed
    into the FITS header.  This routine is a simple wraper around CCDPROC
    functions for easily making changes to these keywords.

    Parameters
    ----------
    files : :obj:`str` or :obj:`~pathlib.Path` or :obj:`list`
        The file(s) for which to update FITS keywords
    keyword : :obj:`str`
        FITS keyword to update
    new_value : :obj:`~typing.Any`
        New value for the FITS keyword
    """
    if isinstance(files, list):
        files = [pathlib.Path(f).resolve() for f in files]
    else:
        files = [pathlib.Path(files).resolve()]

    # Build the IFC
    icl = ccdproc.ImageFileCollection(filenames=files)

    for hdr in icl.headers(overwrite=True):
        # Attempt to get numerical values as numbers, not strings
        try:
            hdr[keyword] = float(new_value)
        except ValueError:
            hdr[keyword] = new_value


# Command-Line Entry Point ===================================================#
def entry_point():
    """Command-Line Entry Point"""
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        prog="fix_ldt_header", description="Fix a keyword in LDT FITS headers"
    )
    parser.add_argument(
        "file",
        action="store",
        type=str,
        nargs="+",
        help="File(s) on which to operate",
    )
    parser.add_argument(
        "keyword",
        action="store",
        type=str,
        help="FITS keyword to change",
    )
    parser.add_argument(
        "new_value",
        action="store",
        type=str,
        help="New header keyword value to insert",
    )
    args = parser.parse_args()

    # Giddy Up!
    fix_ldt_header(files=args.file, keyword=args.keyword, new_value=args.new_value)
