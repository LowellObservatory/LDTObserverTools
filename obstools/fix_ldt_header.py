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
import typing

# 3rd-Party Libraries
import astropy.time
import ccdproc
import numpy as np

# Local Libraries
from obstools import utils

# CONSTANTS


def fix_ldt_header(
    files: str | pathlib.Path | list,
    keyword: str,
    new_value: typing.Any,
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
        Fix RIMAS exposure time keywords for data prior to mid-March 2026
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

    for hdr in icl.headers(overwrite=True):
        # Attempt to get numerical values as numbers, not strings
        try:
            hdr[keyword] = float(new_value)
        except ValueError:
            hdr[keyword] = new_value


def fix_old_rimas(icl: ccdproc.ImageFileCollection):
    """Fix older RIMAS exposure time keywords

    _extended_summary_

    Parameters
    ----------
    icl : :obj:`~ccdproc.ImageFileCollection`
        The Image File Collection object containing the files to be changed
    """
    for hdr in icl.headers(overwrite=True):

        # Check that the file is old enough
        mjd = astropy.time.Time(hdr["DATE"]).mjd
        cutoff = 61114.0  # 2026-03-15T00:00:00
        if mjd > cutoff:
            print("ERROR: Newer RIMAS file formatl cannot update exposure times!")

        # Compute the new values and update the keywords and comments
        tot_exptime = hdr["EXPTIME"]
        frtime = hdr["FRTIME"]
        hdr["EXPTIME"] = (
            np.round(tot_exptime, 2),
            "[s] exposure time of all frames incl. pedestal",
        )
        hdr["EXPTIMEE"] = (
            np.round(tot_exptime - frtime, 2),
            "[s] effective exposure time of reduced frame",
        )


# Command Line Script Infrastructure (borrowed from PypeIt) ==================#
class FixLdtHeader(utils.ScriptBase):
    """Script class for ``fix_ldt_header`` tool

    Script structure borrowed from :class:`pypeit.scripts.scriptbase.ScriptBase`.
    """

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

        # Step 1: Check for the bypass flag only
        temp_parser = argparse.ArgumentParser(add_help=False)
        temp_parser.add_argument("--rimas_exp", action="store_true")
        args, _ = temp_parser.parse_known_args()

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

        # Step 2: Define full parser based on flag presence
        if not args.rimas_exp:
            parser.add_argument(
                "keyword", action="store", type=str, help="FITS keyword to change"
            )
            parser.add_argument(
                "new_value",
                action="store",
                type=str,
                help="New header keyword value to insert",
            )

        parser.add_argument(
            "--rimas_exp",
            action="store_true",
            help="Fix older RIMAS exposure time keywords",
        )

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
