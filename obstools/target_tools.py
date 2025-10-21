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
import re
import sys
import typing

# 3rd-Party Libraries
import astropy.coordinates
import astropy.io.votable
import astropy.table
import astropy.time
import astropy.units as u
import numpy as np

# Local Libraries
from obstools import utils


@dataclasses.dataclass
class TargetList:
    """Generic Target List data class"""

    source: str
    orig_format: str
    tls_collist: list = None
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
        # Write the table to the string buffer in a specified ASCII format
        tls_table = self.data.copy()
        tls_cols = [
            "objname",
            "ra",
            "dec",
            "mu_ra",
            "mu_dec",
            "magnitude",
            "dra",
            "ddec",
            "rotator_pa",
            "rotator_frame",
            "comment1",
        ]

        # Convert 'coords' to string RA/Dec
        present_coords = tls_table["coords"]
        tls_table["ra"] = present_coords.ra.to_string(
            unit=u.hour, sep=":", precision=2, pad=True
        )
        tls_table["dec"] = present_coords.dec.to_string(
            unit=u.deg, sep=":", precision=2, alwayssign=True, pad=True
        )

        # Make blank columns for missing data
        n_rows = len(tls_table)
        for col in tls_cols:
            if col not in tls_table.colnames:
                # Float zeros
                if col in ["mu_ra", "mu_dec", "dra", "ddec", "rotator_pa"]:
                    tls_table[col] = np.zeros(n_rows, dtype=float)
                # Strings
                if col == "rotator_frame":
                    tls_table[col] = np.full(n_rows, "TARGET")
                if col == "comment1":
                    tls_table[col] = np.full(n_rows, "")
                # Magnitude
                if col == "magnitude":
                    tls_table[col] = np.full(n_rows, 100, dtype=float)

        tls_table = tls_table[tls_cols]

        # Surround the following columns with quotes and remove multiple spaces
        quote_cols = ["objname", "comment1"]
        for qc in quote_cols:
            if qc in tls_table.colnames:
                tls_table[qc] = [f'"{re.sub(r'\s+', ' ', v)}"' for v in tls_table[qc]]

        # Set output formats for columns
        for col in ["mu_ra", "mu_dec", "dra", "ddec"]:
            tls_table[col].info.format = ".1f"
        tls_table["magnitude"].info.format = ".2f"
        tls_table["rotator_pa"].info.format = ".0f"
        tls_table.sort("ra")
        longest_name = np.amax([len(name) for name in tls_table["objname"]])
        tls_table["objname"].info.format = f"<{longest_name}"
        # Replace blank comments
        tls_table["comment1"][tls_table["comment1"] == '""'] = '"---"'

        # Create a string buffer to capture the table output
        output_buffer = io.StringIO()
        # Get the string content from the buffer
        tls_table.write(
            output_buffer,
            format="ascii.fixed_width_no_header",
            quotechar="'",
            bookend=False,
            delimiter=None,
        )
        ascii_string = output_buffer.getvalue()

        # Add the TLS header and return the string object
        hdr = (
            "#title=true ra=true dec=true epoch=false muRA=true muDec=true "
            "magnitude=true dRA=true dDec=true rotatorPA=true "
            "rotatorFrame=true comment=true\n#\n"
        )
        return hdr + ascii_string

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
        # Create a copy of the table to re-organize
        rimas_csv = self.data.copy()

        # Do the work of organizing the output
        rimas_csv.rename_columns(["objname", "comment1"], ["ObjectName", "Comment1"])

        # now = astropy.time.Time.now()
        # d_lon = rimas_csv['mu_ra']

        # present_coords = rimas_csv["coords"].spherical_offsets_by(d_lon, d_lat)

        present_coords = rimas_csv["coords"]

        rimas_csv["RA"] = present_coords.ra.to_string(unit=u.hour, sep=":", precision=2)
        rimas_csv["DEC"] = present_coords.dec.to_string(
            unit=u.deg, sep=":", precision=2, alwayssign=True
        )

        # Add blank-ish columns
        n_rows = len(rimas_csv)
        fzeros = ["RAoffset", "DECoffset", "ROToffset", "xi", "eta", "DitherY"]
        izeros = ["Priority"]
        for col in fzeros:
            rimas_csv[col] = np.zeros(n_rows, dtype=float)
        for col in izeros:
            rimas_csv[col] = np.zeros(n_rows, dtype=int)
        # Set the filters
        rimas_csv["Filter1"] = np.full(n_rows, "J")
        rimas_csv["Filter2"] = np.full(n_rows, "H")
        rimas_csv["Filter3"] = np.full(n_rows, "open")
        rimas_csv["Filter4"] = np.full(n_rows, "open")
        # Set the dithers
        rimas_csv["DitherType"] = np.full(n_rows, "None")
        rimas_csv["DitherX"] = np.full(n_rows, 20.0, dtype=float)
        rimas_csv["DitherTotal"] = np.ones(n_rows, dtype=int)
        # Remaining columns
        rimas_csv["BlockID"] = [f"P{v:05d}" for v in np.arange(7000, 7000 + n_rows)]
        rimas_csv["Observer"] = np.full(n_rows, "P. Muirhead / F. Fatmasiefa (BU)")
        rimas_csv["ObjectType"] = np.full(n_rows, "Science")
        rimas_csv["Images"] = np.full(n_rows, 10, dtype=int)
        rimas_csv["IntegrationTime"] = np.full(n_rows, 15.0, dtype=float)
        rimas_csv["Comment2"] = np.full(n_rows, "")

        # Surround the following columns with quotes
        quote_cols = ["Observer", "ObjectName", "ObjectType", "Comment1", "Comment2"]
        for qc in quote_cols:
            if qc in rimas_csv.colnames:
                rimas_csv[qc] = [f'"{v}"' for v in rimas_csv[qc]]

        rimas_csv = rimas_csv[
            "Priority",
            "BlockID",
            "Observer",
            "ObjectName",
            "ObjectType",
            "RA",
            "DEC",
            "RAoffset",
            "DECoffset",
            "ROToffset",
            "Filter1",
            "Filter2",
            "Filter3",
            "Filter4",
            "DitherType",
            "DitherX",
            "DitherY",
            "DitherTotal",
            "Images",
            "IntegrationTime",
            "Comment1",
            "Comment2",
            "xi",
            "eta",
        ]

        # Create a string buffer to capture the table output
        output_buffer = io.StringIO()
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
        print(formats)
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

        # If proper motions not in table, add them as zeros
        if "muRA" not in table.colnames:
            table["muRA"] = np.zeros(len(table))
            table["muDec"] = np.zeros(len(table))

        # Massage the table to conform to class specification
        table["coords"] = astropy.coordinates.SkyCoord(
            ra=table["ra"],
            dec=table["dec"],
            frame="fk5",
            unit=[u.hour, u.deg],
            pm_ra_cosdec=table["muRA"] * u.mas / u.yr,
            pm_dec=table["muDec"] * u.mas / u.yr,
            distance=10 * u.pc,
            obstime=astropy.time.Time("J2000"),
        )

        table.remove_columns(["ra", "dec"])
        table.rename_columns(
            ["title", "comment", "rotatorpa", "rotatorframe"],
            ["objname", "comment1", "rotator_pa", "rotator_frame"],
        )

        return cls(
            source="user file", orig_format=".tls file", tls_collist=cols, data=table
        )

    @classmethod
    def read_from_rimas(cls, filename: pathlib.Path) -> typing.Self:
        """Return a :class:`TargetList` object read from a RIMAS CSV file

        Read in a RIMAS Observation File (``.csv``) file and return an instance
        of this class with the proper data columns filled in.

        Parameters
        ----------
        filename : :obj:`~pathlib.Path`
            RIMAS ``.csv`` file to read

        Returns
        -------
        :obj:`~typing.Self`
            Returns an instance of the class
        """

    @classmethod
    def read_from_simbad(cls, filename: pathlib.Path) -> typing.Self:
        """Return a :class:`TargetList` object read from a Simbad file

        Read in a Simbad catalog file and return an instance of this class with
        the proper data columns filled in.

        Parameters
        ----------
        filename : :obj:`~pathlib.Path`
            Simbad file to read

        Returns
        -------
        :obj:`~typing.Self`
            Returns an instance of the class
        """

    @classmethod
    def read_from_votable(cls, filename: pathlib.Path) -> typing.Self:
        """Return a :class:`TargetList` object read from a VOTable file

        Read in a Virtual Observatory Table catalog file and return an instance
        of this class with the proper data columns filled in.

        Parameters
        ----------
        filename : :obj:`~pathlib.Path`
            VOTable file to read

        Returns
        -------
        :obj:`~typing.Self`
            Returns an instance of the class
        """
        # Read in the VOTable and extract the columns we want
        table = astropy.io.votable.parse_single_table(filename).to_table()
        table = table[
            "MAIN_ID",
            "RA_d",
            "DEC_d",
            "PMRA",
            "PMDEC",
            "FLUX_V",
            "SP_TYPE",
        ]

        # Check for missing values -- set appropriately
        for col in ["PMRA", "PMDEC"]:
            if isinstance(table[col], astropy.table.MaskedColumn):
                table[col] = table[col].filled(0.0)
        if isinstance(table["FLUX_V"], astropy.table.MaskedColumn):
            table["FLUX_V"] = table["FLUX_V"].filled(100.0)
        if isinstance(table["SP_TYPE"], astropy.table.MaskedColumn):
            table["SP_TYPE"] = table["SP_TYPE"].filled("---")

        # Massage the table to conform to class specification
        table["coords"] = astropy.coordinates.SkyCoord(
            ra=table["RA_d"],
            dec=table["DEC_d"],
            frame="fk5",
            unit=[u.deg, u.deg],
            pm_ra_cosdec=table["PMRA"],
            pm_dec=table["PMDEC"],
            obstime=astropy.time.Time("J2000"),
        )
        table.rename_columns(
            ["MAIN_ID", "PMRA", "PMDEC", "FLUX_V", "SP_TYPE"],
            ["objname", "mu_ra", "mu_dec", "magnitude", "comment1"],
        )
        table.remove_columns(["RA_d", "DEC_d"])

        # Return the class instance
        return cls(source="user file", orig_format="VOTable file", data=table)


# Directly callable functions ================================================#
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

    ofile = filename.with_suffix(".rimas.csv")
    with open(ofile, "w", encoding="utf-8") as f_obj:
        f_obj.writelines(rimas_csv)


def convert_votable_to_tls(filename: str | pathlib.Path):
    """Convert a VOTable file to an LDT ``.tls`` file

    Read in the VOTable file and write out the corresponding ``.tls`` file

    Parameters
    ----------
    filename : :obj:`str` or :obj:`~pathlib.Path`
        Name of the VOTable file to convert
    """
    if not isinstance(filename, pathlib.Path):
        filename = pathlib.Path(filename)

    # Get the LDT ``.tls`` format
    ldt_tls = TargetList.read_from_votable(filename).tls_format

    ofile = filename.with_suffix(".tls")
    with open(ofile, "w", encoding="utf-8") as f_obj:
        f_obj.writelines(ldt_tls)


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
        parser.add_argument(
            "--votable",
            action="store_true",
            help="The input is a VOTable to be converted to .tls",
        )
        return parser

    @staticmethod
    def main(args):
        """Main Driver"""

        # Giddy up!
        if args.votable:
            sys.exit(convert_votable_to_tls(args.file))

        sys.exit(convert_tls_to_rimas(args.file))
