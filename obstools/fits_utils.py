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

"""LDTObserverTools contains python ports of various LDT Observer Tools

Lowell Discovery Telescope (Lowell Observatory: Flagstaff, AZ)
http://www.lowell.edu

This file contains FITS Header utility routines.
"""

# Built-In Libraries
import argparse

# 3rd-Party Libraries
import ccdproc

# Local Libraries

# CONSTANTS


def fix_ldt_header(files, keyword, new_value):

    icl = ccdproc.ImageFileCollection(files)

    for hdr in icl.headers(overwrite=True):
        hdr[keyword] = new_value


def entry_point():
    """Command-Line Entry Point for fix_ldt_header()

    _extended_summary_
    """
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
