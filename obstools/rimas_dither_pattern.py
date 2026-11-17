# -*- coding: utf-8 -*-
#
#  This file is part of LDTObserverTools.
#
#   This Source Code Form is subject to the terms of the Mozilla Public
#   License, v. 2.0. If a copy of the MPL was not distributed with this
#   file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
#  Created on 19-Aug-2026
#
#  @author: tbowers

"""Create ECSV sidecars that correct RIMAS dither metadata

LDTObserverTools contains python ports of various LDT Observer Tools

Lowell Discovery Telescope (Lowell Observatory: Flagstaff, AZ)
https://lowell.edu

The generated tables identify the dither pattern, position, and sequence for
all commanded RIMAS frames, including exposures that were not written to disk.
"""

# Built-In Libraries
import argparse
from collections.abc import Iterable, Sequence
import dataclasses
import pathlib

# 3rd-Party Libraries
import astropy.table

# Local Libraries
from obstools import utils

# Module-level Constants
PATTERNS = {"ABBA": ("A", "B", "B", "A"), "ONOFF": ("On", "Off")}


@dataclasses.dataclass(frozen=True)
class DitherSequence:
    """Represent the recorded frames in one commanded dither sequence.

    Parameters
    ----------
    frames : tuple of int
        Commanded frame numbers in strictly increasing order, including frames
        that were not written to disk.
    sequence_id : int
        Positive identifier for the commanded sequence.

    Attributes
    ----------
    frames : tuple of int
        Recorded frame numbers in the sequence.
    sequence_id : int
        Identifier for the sequence.
    """

    frames: tuple[int, ...]
    sequence_id: int


def normalize_pattern(sequence_type: str) -> str:
    """Return the canonical name of a supported dither pattern.

    Parameters
    ----------
    sequence_type : str
        Pattern name. Matching ignores case and surrounding whitespace.

    Returns
    -------
    str
        Uppercase canonical pattern name.

    Raises
    ------
    ValueError
        If ``sequence_type`` is not supported.
    """
    pattern = sequence_type.strip().upper()
    if pattern not in PATTERNS:
        raise ValueError(
            f"Unsupported sequence type {sequence_type!r}; choose {', '.join(PATTERNS)}."
        )
    return pattern


def parse_frame_list(value: str | Iterable[str]) -> tuple[int, ...]:
    """Parse frame numbers from one or more strings.

    Parameters
    ----------
    value : str or iterable of str
        Frame numbers supplied as individual strings, comma-separated strings,
        or a mixture of both. Whitespace and empty fields are ignored.

    Returns
    -------
    tuple of int
        Parsed frame numbers in input order.

    Raises
    ------
    ValueError
        If a nonempty field is not an integer or no frames are given.
    """
    values = (value,) if isinstance(value, str) else value
    try:
        frames = tuple(
            int(item.strip())
            for entry in values
            for item in entry.split(",")
            if item.strip()
        )
    except ValueError as exc:
        raise ValueError(
            f"Invalid frame list {value!r}; use space- or comma-separated integers."
        ) from exc
    if not frames:
        raise ValueError("A dither sequence must contain at least one frame.")
    return frames


def validate_sequence(sequence: DitherSequence, sequence_type: str) -> None:
    """Validate one sequence against a dither pattern.

    Parameters
    ----------
    sequence : DitherSequence
        Sequence to validate.
    sequence_type : str
        Supported dither pattern name.

    Raises
    ------
    ValueError
        If the pattern or sequence is invalid.
    """
    pattern, frames = normalize_pattern(sequence_type), sequence.frames
    if not frames:
        raise ValueError("A dither sequence must contain at least one frame.")
    if sequence.sequence_id < 1:
        raise ValueError("Sequence IDs must be positive.")
    if any(frame < 0 for frame in frames):
        raise ValueError("Frame numbers must be non-negative.")
    if any(b <= a for a, b in zip(frames, frames[1:])):
        raise ValueError(
            f"Frames in sequence {sequence.sequence_id} must increase strictly."
        )
    if frames[-1] - frames[0] >= len(PATTERNS[pattern]):
        raise ValueError(
            f"Frames in sequence {sequence.sequence_id} span more than one {pattern} pattern."
        )


def build_correction_table(
    sequence_type: str,
    sequences: Sequence[DitherSequence],
    date: str,
    arms: Iterable[str] = ("YJ", "HK"),
    source: str = "observing log",
) -> astropy.table.Table:
    """Build a RIMAS dither metadata correction table.

    Parameters
    ----------
    sequence_type : str
        Supported dither pattern name.
    sequences : sequence of DitherSequence
        Recorded frames grouped by commanded sequence.
    date : str
        Observing date in ``YYYYMMDD`` format.
    arms : iterable of str, optional
        Detector arms to include; supported values are ``"YJ"`` and ``"HK"``.
    source : str, optional
        Description of the correction metadata source.

    Returns
    -------
    astropy.table.Table
        Correction rows for every recorded frame and selected arm.

    Raises
    ------
    ValueError
        If inputs are invalid, identifiers or frames are duplicated, or no rows
        can be produced.
    """
    pattern = normalize_pattern(sequence_type)
    if len(date) != 8 or not date.isdigit():
        raise ValueError("Date must have YYYYMMDD format.")
    use_arms = tuple(str(arm).strip().upper() for arm in arms)
    if not use_arms or any(arm not in ("YJ", "HK") for arm in use_arms):
        raise ValueError("Arms must contain one or both of YJ and HK.")
    if len(set(use_arms)) != len(use_arms):
        raise ValueError("Arms must not contain duplicates.")

    rows, seen_frames, seen_ids = [], set(), set()
    for sequence in sequences:
        validate_sequence(sequence, pattern)
        if sequence.sequence_id in seen_ids:
            raise ValueError(f"Duplicate sequence ID {sequence.sequence_id}.")
        seen_ids.add(sequence.sequence_id)
        first = sequence.frames[0]
        for frame in sequence.frames:
            if frame in seen_frames:
                raise ValueError(f"Frame {frame} occurs in more than one sequence.")
            seen_frames.add(frame)
            dithnum = frame - first + 1
            for arm in use_arms:
                rows.append(
                    (
                        f"{date}.rimas.{frame:04d}.{arm}.fits",
                        pattern,
                        dithnum,
                        PATTERNS[pattern][dithnum - 1],
                        sequence.sequence_id,
                        source,
                    )
                )
    if not rows:
        raise ValueError("At least one dither sequence is required.")
    table = astropy.table.Table(
        rows=rows,
        names=("filename", "dithpat", "dithnum", "dithpos", "dithseq", "source"),
    )
    table.meta["description"] = "RIMAS dither metadata corrections"
    return table


def write_correction_table(
    table: astropy.table.Table, output: str | pathlib.Path, overwrite: bool = False
) -> pathlib.Path:
    """Write a correction table in ECSV format.

    Parameters
    ----------
    table : astropy.table.Table
        Correction table to write.
    output : str or pathlib.Path
        Destination path. Missing parent directories are created.
    overwrite : bool, optional
        Whether to replace an existing file.

    Returns
    -------
    pathlib.Path
        Absolute path of the written file.
    """
    output_path = pathlib.Path(output).expanduser().resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    table.write(output_path, format="ascii.ecsv", overwrite=overwrite)
    return output_path


class RimasDitherPattern(utils.ScriptBase):
    """Provide the CLI for creating RIMAS dither sidecars.

    Notes
    -----
    One ``--sequence`` option is accepted per commanded sequence. Sequence
    identifiers follow the order of those options.
    """

    @classmethod
    def get_parser(
        cls,
        description: str | None = None,
        width: int | None = None,
        formatter: type[
            argparse.HelpFormatter
        ] = argparse.ArgumentDefaultsHelpFormatter,
    ) -> argparse.ArgumentParser:
        """Construct the command-line argument parser.

        Parameters
        ----------
        description : str or None, optional
            Parser description. Retained for base-class compatibility.
        width : int or None, optional
            Maximum width of formatted help output.
        formatter : type of argparse.HelpFormatter, optional
            Help formatter class.

        Returns
        -------
        argparse.ArgumentParser
            Configured command-line parser.
        """
        parser = super().get_parser(
            description="Create an ECSV sidecar correcting RIMAS dither metadata.",
            width=width,
            formatter=formatter,
        )
        parser.add_argument("sequence_type", choices=tuple(PATTERNS))
        parser.add_argument(
            "-s",
            "--sequence",
            action="append",
            nargs="+",
            required=True,
            metavar="FRAME",
            help=(
                "Full frame list for one commanded sequence, including frames not "
                "written to disk. Repeat for each sequence; Bash brace expansion "
                "such as -s {86..89} is supported."
            ),
        )
        parser.add_argument("--date", required=True, help="Observing date (YYYYMMDD)")
        parser.add_argument(
            "--arms", nargs="+", default=("YJ", "HK"), choices=("YJ", "HK")
        )
        parser.add_argument("-o", "--output", default="rimas_dither_corrections.ecsv")
        parser.add_argument("--source", default="observing log")
        parser.add_argument("--overwrite", action="store_true")
        return parser

    @staticmethod
    def main(args: argparse.Namespace) -> None:
        """Create and write a correction table from parsed arguments.

        Parameters
        ----------
        args : argparse.Namespace
            Parsed command-line arguments produced by :meth:`get_parser`.
        """
        sequences = [
            DitherSequence(parse_frame_list(value), index)
            for index, value in enumerate(args.sequence, start=1)
        ]
        table = build_correction_table(
            args.sequence_type, sequences, args.date, args.arms, args.source
        )
        output = write_correction_table(table, args.output, args.overwrite)
        print(f"Wrote {len(table)} corrections to {output}")
