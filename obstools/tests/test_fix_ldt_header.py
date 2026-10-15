"""Tests for FITS header correction and its command-line modes."""

import argparse

from astropy.io import fits
import pytest

from obstools.fix_ldt_header import FixLdtHeader, fix_ldt_header


def make_rimas_file(path, date, exposure=12.345, frame=0.615):
    header = fits.Header()
    header["DATE-OBS"] = date
    header["EXPTIME"] = exposure
    header["FRTIME"] = frame
    fits.PrimaryHDU(header=header).writeto(path)


def test_rimas_exposure_correction(tmp_path):
    image = tmp_path / "old-rimas.fits"
    make_rimas_file(image, "2026-03-14T23:59:59")

    fix_ldt_header(image, rimas_exp=True)

    header = fits.getheader(image)
    assert header["EXPTIME"] == 12.35
    assert header["EXPTIMEE"] == 11.73


def test_existing_keyword_edit(tmp_path):
    image = tmp_path / "ordinary.fits"
    make_rimas_file(image, "2026-03-14")

    fix_ldt_header(image, "OBJECT", "Target")

    assert fits.getheader(image)["OBJECT"] == "Target"
    assert "EXPTIMEE" not in fits.getheader(image)


def test_rimas_batch_is_validated_before_writing(tmp_path):
    old = tmp_path / "old.fits"
    new = tmp_path / "new.fits"
    make_rimas_file(old, "2026-03-14T23:59:59")
    make_rimas_file(new, "2026-03-15T00:00:00")

    with pytest.raises(ValueError, match="before 2026-03-15"):
        fix_ldt_header([old, new], rimas_exp=True)

    assert fits.getheader(old)["EXPTIME"] == 12.345
    assert "EXPTIMEE" not in fits.getheader(old)


def test_rimas_rejects_missing_frame_time(tmp_path):
    image = tmp_path / "missing.fits"
    make_rimas_file(image, "2026-03-14")
    with fits.open(image, mode="update") as hdus:
        del hdus[0].header["FRTIME"]

    with pytest.raises(ValueError, match="EXPTIME and FRTIME"):
        fix_ldt_header(image, rimas_exp=True)
    assert "EXPTIMEE" not in fits.getheader(image)


def test_existing_keyword_cli_and_rimas_cli():
    standard = FixLdtHeader.parse_args(["image.fits", "OBJECT", "Target"])
    assert standard.file == ["image.fits"]
    assert standard.keyword == "OBJECT"
    assert standard.new_value == "Target"

    rimas = FixLdtHeader.parse_args(["--rimas-exp", "image.fits"])
    assert rimas.file == ["image.fits"]
    assert rimas.rimas_exp
    assert not hasattr(rimas, "keyword")
    assert FixLdtHeader.parse_args(["--rimas_exp", "image.fits"]).rimas_exp

    parser = FixLdtHeader.get_parser(formatter=argparse.RawDescriptionHelpFormatter)
    assert isinstance(
        parser.formatter_class("test"), argparse.RawDescriptionHelpFormatter
    )
