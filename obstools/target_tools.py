# -*- coding: utf-8 -*-
#
#  This file is part of LDTObserverTools.
#
#   This Source Code Form is subject to the terms of the Mozilla Public
#   License, v. 2.0. If a copy of the MPL was not distributed with this
#   file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
#  Created on 17-Oct-2025
#
#  @author: tbowers
# pylint: disable=c-extension-no-member

"""Observing Target List Creation and Manipulation Tools

LDTObserverTools contains python ports of various LDT Observer Tools

Lowell Discovery Telescope (Lowell Observatory: Flagstaff, AZ)
https://lowell.edu

This file will contain a GUI tool for creating and manipulating observing
target list files.  At the moment, it's first goal is to convert extant
``.tls`` files into RIMAS ``Observation List`` CSV files.

Eventually, this tool should be able to validate user-created ``.tls`` and
LMI slew-dither pattern files, create ``.tls`` files from external sources
(*e.g.*, SIMBAD table output), and any direct user input via the GUI.
"""

# Built-In Libraries
import argparse
import dataclasses
import io
import pathlib
import typing

# 3rd-Party Libraries
import astropy.coordinates
import astropy.table
import astropy.units as u

# Local Libraries
from obstools import utils


@dataclasses.dataclass
class TargetList:
    """Generic Target List data class"""

    source: str
    orig_format: str
    data: astropy.table.Table = dataclasses.field(default_factory=astropy.table.Table)

    def __post_init__(self):
        # Structure of the table with predefined columns
        if not self.data:
            #
            self.data.add_column(astropy.table.Column(dtype=str, name="objname"))
            self.data.add_column(astropy.table.Column(dtype=str, name="objtype"))
            self.data.add_column(astropy.table.Column(dtype=str, name="observer"))
            self.data.add_column(
                astropy.table.Column(dtype=astropy.coordinates.SkyCoord, name="coords")
            )
            self.data.add_column(astropy.table.Column(dtype=str, name="epoch"))
            # target_star.spherical_offsets_by(1.3*u.arcmin, -0.7*u.arcmin)
            self.data.add_column(
                astropy.table.Column(dtype=u.Quantity, name="ra_offset")
            )
            self.data.add_column(
                astropy.table.Column(dtype=u.Quantity, name="dec_offset")
            )
            self.data.add_column(
                astropy.table.Column(dtype=u.Quantity, name="rot_offset")
            )
            self.data.add_column(astropy.table.Column(dtype=u.Quantity, name="mu_ra"))
            self.data.add_column(astropy.table.Column(dtype=u.Quantity, name="mu_dec"))
            self.data.add_column(astropy.table.Column(dtype=u.Quantity, name="dra"))
            self.data.add_column(astropy.table.Column(dtype=u.Quantity, name="ddec"))
            self.data.add_column(astropy.table.Column(dtype=float, name="magnitude"))
            self.data.add_column(astropy.table.Column(dtype=str, name="filter1"))
            self.data.add_column(astropy.table.Column(dtype=str, name="filter2"))
            self.data.add_column(astropy.table.Column(dtype=str, name="filter3"))
            self.data.add_column(astropy.table.Column(dtype=str, name="filter4"))
            self.data.add_column(astropy.table.Column(dtype=int, name="num_exp"))
            self.data.add_column(astropy.table.Column(dtype=float, name="exptime"))
            self.data.add_column(astropy.table.Column(dtype=str, name="dither_type"))
            self.data.add_column(
                astropy.table.Column(dtype=u.Quantity, name="dither_x")
            )
            self.data.add_column(
                astropy.table.Column(dtype=u.Quantity, name="dither_y")
            )
            self.data.add_column(
                astropy.table.Column(dtype=u.Quantity, name="dither_total")
            )
            self.data.add_column(astropy.table.Column(dtype=str, name="rotator_frame"))
            self.data.add_column(astropy.table.Column(dtype=float, name="rotator_pa"))
            self.data.add_column(astropy.table.Column(dtype=str, name="comment1"))
            self.data.add_column(astropy.table.Column(dtype=str, name="comment2"))
            self.data.add_column(astropy.table.Column(dtype=str, name="command_option"))

    def __str__(self) -> str:
        """Print a string representation of the object

        Includes database source, object ID, UT start time, ephemeris step
        size, and total number of points in the table.

        Returns
        -------
        :obj:`str`
            A brief description of the object
        """
        return (
            f"Generic Target List from {self.source}, original format: {self.orig_format} "
            f"N items: {len(self.data)}"
        )

    @property
    def tls_format(self) -> str:
        """Produce the TLS format, ready to write to disk

        #title=true ra=true dec=true epoch=false muRA=true muDec=true magnitude=true dRA=true dDec=true rotatorPA=true rotatorFrame=true comment=true
        #
        "Vega - Alpha Lyrae" 18:36:56.336 +38:47:01.28 200.94 286.23 0.03 1.0 2.0 3.0 Target "line 1"
        "HR 173 - HD 3795" 00:40:32.809 -23:48:17.56 635.18 -363.56 6.14 1.0 2.0 3.0 Fixed "line 2"

        Returns
        -------
         :obj:`str`
             The write-ready string ready for saving to disk
        """
        # Create a string buffer to capture the table output
        output_buffer = io.StringIO()

        # Write the table to the string buffer in a specified ASCII format
        tcs_ephem = self.data.copy()
        tcs_ephem.remove_columns(["azimuth", "elevation", "solar_elon", "lunar_elon"])
        tcs_ephem.write(output_buffer, format="ascii.no_header")

        # Get the string content from the buffer
        ascii_string = output_buffer.getvalue()

        # Return the string object
        return ascii_string

    @property
    def slewdither_format(self) -> str:
        """Produce the LMI slew-dither format, ready to write to disk

        #title=true ra=true dec=true exposureTime=true numExposures=true filter=true muRA=false muDec=false epoch=false dRA=false dDec=false rotatorPA=false rotatorFrame=false xi=true eta=true comment=true commandOption=true
        #
        "object 1" 00:15:00.00 +45:30:00.0  10.0  1 R  0  0 "Object 1, slew, no dither"         Slew
        "object 1" 00:15:00.00 +46:30:00.0  10.0  1 R 30 30 "Object 1, dither 30,30"            Dither
        "object 1" 00:15:00.00 +46:30:00.0  10.0  1 R 60 60 "Object 1, dither 60,60"            Dither
        "object 1" 00:15:00.00 +46:30:00.0   0.0  1 R  0  0 "Object 1, dither 0,0 - clear"      Dither
        "object 2" 00:37:00.00 +47:30:00.0   7.5  2 V  0  0 "Object 2, slew no dither"          Slew
        "object 2" 00:37:00.00 +47:30:00.0   7.5  2 V 30 30 "Object 2, dither 30,30"            Dither
        "object 2" 00:37:00.00 +47:30:00.0   7.5  2 V 60 60 "Object 2, dither 60,60"            Dither
        "object 2" 00:37:00.00 +47:30:00.0   0.0  2 V  0  0 "Object 2, dither 0,0 - clear"      Dither

        Returns
        -------
         :obj:`str`
             The write-ready string ready for saving to disk
        """
        # Create a string buffer to capture the table output
        output_buffer = io.StringIO()

        # Write the table to the string buffer in a specified ASCII format
        tcs_ephem = self.data.copy()
        tcs_ephem.remove_columns(["azimuth", "elevation", "solar_elon", "lunar_elon"])
        tcs_ephem.write(output_buffer, format="ascii.no_header")

        # Get the string content from the buffer
        ascii_string = output_buffer.getvalue()

        # Return the string object
        return ascii_string

    @property
    def rimas_format(self) -> str:
        """Produce the RIMAS CSV format, ready to write to disk

        Priority,BlockID,Observer,ObjectName,ObjectType,RA,DEC,RAoffset,DECoffset,ROToffset,Filter1,Filter2,Filter3,Filter4,DitherType,DitherX,DitherY,DitherTotal,Images,IntegrationTime,Comment1,Comment2,xi,eta
        0,P05315,"joe","test","Science",19:47:29.2,+64:29:12.4,0.00,0.00,0.00,J,H,open,open,Random,20.0,0.0,10,1,15.00,"","",0.000000,0.000000

        Returns
        -------
         :obj:`str`
             The write-ready string ready for saving to disk
        """
        # Create a string buffer to capture the table output
        output_buffer = io.StringIO()

        # Write the table to the string buffer in a specified ASCII format
        rimas_csv = self.data.copy()

        # Do the work of organizing the output
        rimas_csv.rename_columns(["objname", "comment1"], ["ObjectName", "Comment1"])

        rimas_csv["RA"] = rimas_csv["coords"].ra.to_string(
            unit=u.hour, sep=":", precision=2
        )
        rimas_csv["DEC"] = rimas_csv["coords"].dec.to_string(
            unit=u.deg, sep=":", precision=2, alwayssign=True
        )

        # Surround the following columns with quotes
        quote_cols = ["Observer", "ObjectName", "ObjectType", "Comment1", "Comment2"]
        for qc in quote_cols:
            if qc in rimas_csv.colnames:
                rimas_csv[qc] = [f'"{v}"' for v in rimas_csv[qc]]

        rimas_csv.pprint()

        # Get the string content from the buffer
        rimas_csv.write(output_buffer, format="ascii.csv", quotechar="'")
        ascii_string = output_buffer.getvalue()

        # Return the string object
        return ascii_string

    @classmethod
    def read_from_tls(cls, filename: pathlib.Path) -> typing.Self:
        """Return a :class:`TargetList` object read from a ``.tls`` file

        Read in an LDT Target List (``.tls``) file and return an instance
        of this class with the proper data columns filled in.

        Parameters
        ----------
        filename : :obj:`~pathlib.Path`
            ``.tls`` file to read

        Returns
        -------
        :obj:`~typing.Self`
            Returns an instance of the class
        """
        # Do the reading and dropping of things into the proper columns
        with open(filename, "r", encoding="utf-8") as f_obj:
            lines = f_obj.readlines()

        # First line is the format
        formats = lines[0][1:].split()
        if (n_lab := len(formats)) != 12:
            raise ValueError(
                "TLS file does not match TLS format!  "
                f"Should have 12 labels, has {n_lab}"
            )
        # Get the included columns
        cols = []
        for fmt in formats:
            if fmt.split("=")[1].strip().lower() == "true":
                cols.append(fmt.split("=")[0].strip().lower())

        # Load the remainder of the file into a table-like object (2nd line is blank comment)
        table_array = []
        for line in lines[2:]:
            parts = utils.split_preserving_quotes(line)
            table_array.append({c: p.replace('"', "") for c, p in zip(cols, parts)})

        table = astropy.table.Table(table_array)

        # Massage the table to conform to class specification
        table["coords"] = astropy.coordinates.SkyCoord(
            ra=table["ra"], dec=table["dec"], frame="fk5", unit=[u.hour, u.deg]
        )
        table.remove_columns(["ra", "dec"])
        table.rename_columns(
            ["title", "comment", "rotatorpa", "rotatorframe"],
            ["objname", "comment1", "rotator_pa", "rotator_frame"],
        )

        return cls(source="user file", orig_format=".tls file", data=table)


def convert_tls_to_rimas(filename: str | pathlib.Path):
    """Convert an LDT ``.tls`` file to RIMAS Observation List ``.csv``

    Read in the ``.tls`` file and write out the corresponding ``.csv`` file.

    Parameters
    ----------
    filename : :obj:`str` or :obj:`~pathlib.Path`
        Name of the ``.tls`` file to convert
    """
    if not isinstance(filename, pathlib.Path):
        filename = pathlib.Path(filename)

    # Get the RIMAS ``.csv`` format
    rimas_csv = TargetList.read_from_tls(filename).rimas_format

    print(rimas_csv)

    ofile = filename.with_suffix(".rimas.csv")
    with open(ofile, "w", encoding="utf-8") as f_obj:
        f_obj.writelines(rimas_csv)


# Command Line Script Infrastructure (borrowed from PypeIt) ==================#
class TargetTools(utils.ScriptBase):
    """Script class for ``target_tools`` tool

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

        parser = super().get_parser(
            description="Obserrving target list manipulator",
            width=width,
        )
        parser.add_argument(
            "file",
            action="store",
            type=str,
            help="File(s) on which to operate",
        )
        return parser

    @staticmethod
    def main(args):
        """Main Driver"""
        convert_tls_to_rimas(args.file)
