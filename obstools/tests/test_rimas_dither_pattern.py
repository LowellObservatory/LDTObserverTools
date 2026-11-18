"""Tests for the RIMAS dither sidecar generator."""

import pytest
from astropy.table import Table
from obstools.rimas_dither_pattern import (
    DitherSequence,
    RimasDitherPattern,
    build_correction_table,
    parse_frame_list,
    write_correction_table,
)


def test_abba_with_missing_frames():
    sequences = [
        DitherSequence((86, 87, 89), 1),
        DitherSequence((90, 91, 92, 93), 2),
        DitherSequence((94, 96), 3),
    ]
    table = build_correction_table("ABBA", sequences, "20260305", arms=("HK",))
    assert list(table["dithnum"]) == [1, 2, 4, 1, 2, 3, 4, 1, 3]
    assert list(table["dithpos"]) == ["A", "B", "A", "A", "B", "B", "A", "A", "B"]
    assert list(table["dithseq"]) == [1, 1, 1, 2, 2, 2, 2, 3, 3]


def test_both_arms_default():
    table = build_correction_table("ONOFF", [DitherSequence((7, 8), 1)], "20260305")
    assert len(table) == 4
    assert table[0]["filename"] == "20260305.rimas.0007.YJ.fits"


def test_validation():
    with pytest.raises(ValueError, match="span more than one ABBA"):
        build_correction_table("ABBA", [DitherSequence((86, 90), 1)], "20260305")
    sequences = [DitherSequence((86, 87), 1), DitherSequence((87, 88), 2)]
    with pytest.raises(ValueError, match="more than one sequence"):
        build_correction_table("ABBA", sequences, "20260305")


def test_parse_and_write(tmp_path):
    assert parse_frame_list("86, 87,89") == (86, 87, 89)
    table = build_correction_table(
        "ABBA", [DitherSequence((86, 87, 89), 1)], "20260305", arms=("HK",)
    )
    output = write_correction_table(table, tmp_path / "corrections.ecsv")
    reread = Table.read(output, format="ascii.ecsv")
    assert reread.colnames == table.colnames
    assert list(reread["dithnum"]) == [1, 2, 4]


def test_cli_bash_brace_expansion():
    parser = RimasDitherPattern.get_parser()
    args = parser.parse_args(
        [
            "ABBA",
            "--date",
            "20260305",
            "-s",
            "86",
            "87",
            "88",
            "89",
            "-s",
            "90",
            "91",
            "92",
            "93",
            "-s",
            "94",
            "95",
            "96",
            "97",
        ]
    )

    sequences = [parse_frame_list(value) for value in args.sequence]

    assert sequences == [
        (86, 87, 88, 89),
        (90, 91, 92, 93),
        (94, 95, 96, 97),
    ]
