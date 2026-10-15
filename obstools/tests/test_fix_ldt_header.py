# SPDX-License-Identifier: MPL-2.0
"""Regression tests for FITS header repair."""

from obstools.fix_ldt_header import fix_old_rimas


def test_fix_old_rimas_skips_new_format():
    """Leave files newer than the legacy format cutoff untouched."""
    old = {"DATE": "2026-03-01T00:00:00", "EXPTIME": 12.345, "FRTIME": 1.0}
    new = {"DATE": "2026-04-01T00:00:00", "EXPTIME": 12.345, "FRTIME": 1.0}

    class Collection:
        """Minimal header collection exposing the ccdproc iterator API."""

        def headers(self, overwrite):
            """Yield both header formats."""
            assert overwrite
            yield old
            yield new

    fix_old_rimas(Collection())
    assert "EXPTIMEE" in old
    assert new == {"DATE": "2026-04-01T00:00:00", "EXPTIME": 12.345, "FRTIME": 1.0}
