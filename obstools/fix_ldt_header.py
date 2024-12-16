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
import pathlib

# 3rd-Party Libraries
import ccdproc

# Local Libraries
from obstools import utils

# CONSTANTS


def fix_ldt_header(files: str | pathlib.Path | list, keyword: str, new_value):
    """Change FITS header keywords

    Sometimes at the telescope, incorrect or incomplete information is placed
    into the FITS header.  This routine is a simple wrapper around CCDPROC
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


# Command Line Script Infrastructure (borrowed from PypeIt) ==================#
class FixLdtHeader(utils.ScriptBase):
    """Script class for ``fix_ldt_header`` tool

    Script structure borrowed from :class:`pypeit.scripts.scriptbase.ScriptBase`.
    """

    @classmethod
    def get_parser(cls, width=None):
        """Construct the command-line argument parser.

        Parameters
        ----------
        description : :obj:`str`, optional
            A short description of the purpose of the script.
        width : :obj:`int`, optional
            Restrict the width of the formatted help output to be no longer
            than this number of characters, if possible given the help
            formatter.  If None, the width is the same as the terminal
            width.
        formatter : :obj:`~argparse.HelpFormatter`
            Class used to format the help output.

        Returns
        -------
        :obj:`~argparse.ArgumentParser`
            Command-line interpreter.
        """

        parser = super().get_parser(
            description="Fix a keyword in LDT FITS headers", width=width
        )
        parser.add_argument(
            "file",
            action="store",
            type=str,
            nargs="+",
            help="File(s) on which to operate",
        )
        parser.add_argument(
            "keyword", action="store", type=str, help="FITS keyword to change"
        )
        parser.add_argument(
            "new_value",
            action="store",
            type=str,
            help="New header keyword value to insert",
        )
        return parser

    @staticmethod
    def main(args):
        """Main Driver

        Simple function that calls the fixer.
        """
        # Giddy up!
        fix_ldt_header(files=args.file, keyword=args.keyword, new_value=args.new_value)
