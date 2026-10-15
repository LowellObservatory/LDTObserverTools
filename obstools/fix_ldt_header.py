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
import math
import pathlib
import sys

# 3rd-Party Libraries
import astropy.time
import ccdproc

# Local Libraries
from obstools import utils

# CONSTANTS


def fix_ldt_header(
    files: str | pathlib.Path | list,
    keyword: str = None,
    new_value=None,
    rimas_exp: bool = False,
):
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
    rimas_exp : :obj:`bool`, optional
        Correct exposure-time keywords in RIMAS files taken before 2026-03-15.
    """
    if isinstance(files, list):
        files = [pathlib.Path(f).resolve() for f in files]
    else:
        files = [pathlib.Path(files).resolve()]

    # Build the IFC
    icl = ccdproc.ImageFileCollection(filenames=files)

    if rimas_exp:
        fix_old_rimas(icl)
        return

    if keyword is None or new_value is None:
        raise ValueError("keyword and new_value are required for a keyword edit")

    for hdr in icl.headers(overwrite=True):
        # Attempt to get numerical values as numbers, not strings
        try:
            hdr[keyword] = float(new_value)
        except ValueError:
            hdr[keyword] = new_value


def fix_old_rimas(icl: ccdproc.ImageFileCollection):
    """Correct exposure-time keywords in older RIMAS FITS files.

    Validate the entire collection before overwriting any header, so a newer
    or malformed file cannot leave a partially corrected batch.
    """
    cutoff_mjd = 61114.0  # 2026-03-15T00:00:00 UTC
    for hdr in icl.headers():
        date = hdr.get("DATE-OBS", hdr.get("DATE"))
        if date is None:
            raise ValueError("RIMAS header has neither DATE-OBS nor DATE")
        if astropy.time.Time(date).mjd >= cutoff_mjd:
            raise ValueError("RIMAS exposure correction requires data before 2026-03-15")
        try:
            total = float(hdr["EXPTIME"])
            frame = float(hdr["FRTIME"])
        except (KeyError, TypeError, ValueError) as err:
            raise ValueError("RIMAS header needs numeric EXPTIME and FRTIME") from err
        if not (math.isfinite(total) and math.isfinite(frame) and 0 <= frame <= total):
            raise ValueError("RIMAS exposure times must be finite and nonnegative")

    for hdr in icl.headers(overwrite=True):
        total = float(hdr["EXPTIME"])
        frame = float(hdr["FRTIME"])
        hdr["EXPTIME"] = (
            round(total, 2),
            "[s] exposure time of all frames incl. pedestal",
        )
        hdr["EXPTIMEE"] = (
            round(total - frame, 2),
            "[s] effective exposure time of reduced frame",
        )


# Command Line Script Infrastructure (borrowed from PypeIt) ==================#
class FixLdtHeader(utils.ScriptBase):
    """Script class for ``fix_ldt_header`` tool

    Script structure borrowed from :class:`pypeit.scripts.scriptbase.ScriptBase`.
    """

    @classmethod
    def parse_args(cls, options=None):
        """Select the RIMAS or keyword-edit parser from the supplied arguments."""
        argv = list(sys.argv[1:] if options is None else options)
        if "--rimas-exp" not in argv and "--rimas_exp" not in argv:
            return super().parse_args(argv)

        parser = super().get_parser(
            description="Correct exposure-time keywords in older RIMAS FITS files"
        )
        parser.add_argument("--rimas-exp", "--rimas_exp", action="store_true")
        parser.add_argument("file", nargs="+", help="RIMAS FITS files to correct")
        return parser.parse_args(argv)

    @classmethod
    def get_parser(
        cls,
        description: str = None,
        width: int = None,
        formatter: argparse.HelpFormatter = argparse.ArgumentDefaultsHelpFormatter,
    ):
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
            description=description or "Fix a keyword in LDT FITS headers",
            width=width,
            formatter=formatter,
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
        parser.epilog += "\nFor older RIMAS files: fix_ldt_header --rimas-exp FILE [FILE ...]"
        return parser

    @staticmethod
    def main(args):
        """Main Driver

        Simple function that calls the fixer.
        """
        # Giddy up!
        fix_ldt_header(
            files=args.file,
            keyword=getattr(args, "keyword", None),
            new_value=getattr(args, "new_value", None),
            rimas_exp=getattr(args, "rimas_exp", False),
        )
