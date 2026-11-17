# SPDX-License-Identifier: MPL-2.0
"""Regression tests for NEOCP ephemeris requests."""

import datetime

from obstools.neocp_ephem import neocp_ephem


def test_neocp_ephem_uses_current_utc_date(monkeypatch, tmp_path):
    """Request current UTC dates instead of a historical observing night."""
    captured = {}

    class Response:
        """Minimal Scout API response for an empty ephemeris."""

        url = "https://example.invalid/scout"

        def raise_for_status(self):
            """Accept the mocked response as successful."""

        def json(self):
            """Provide an empty ephemeris."""
            return {"count": 0, "eph": []}

    def fake_get(url, *, params, timeout):
        """Capture the request without using the network."""
        captured.update(params)
        return Response()

    monkeypatch.setattr("obstools.neocp_ephem.requests.get", fake_get)
    monkeypatch.chdir(tmp_path)
    before = datetime.datetime.now(datetime.timezone.utc).replace(tzinfo=None)
    neocp_ephem("TEST")
    after = datetime.datetime.now(datetime.timezone.utc).replace(tzinfo=None)

    start = datetime.datetime.fromisoformat(captured["eph-start"])
    stop = datetime.datetime.fromisoformat(captured["eph-stop"])
    assert before - datetime.timedelta(seconds=1) <= start <= after
    assert stop - start == datetime.timedelta(days=1)
    assert (tmp_path / "TEST_LDT.eph").read_text() == "FK5 J2000.0 2000.0\n"
