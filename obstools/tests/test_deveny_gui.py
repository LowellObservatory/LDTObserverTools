# SPDX-License-Identifier: MPL-2.0
"""Focused checks for the combined DeVeny interface."""

import os
import sys

if sys.platform.startswith("linux"):
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest
from PyQt6 import QtWidgets

from obstools.deveny_gui import DevenyToolsWindow


@pytest.fixture(scope="module")
def window():
    """Create a headless DeVeny window for widget behavior checks."""
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    _ = app
    result = DevenyToolsWindow()
    yield result
    result.close()


def test_tab_calculations(window):
    """The new tabs expose both grating directions and focus sequence values."""
    assert [window.tabWidget.tabText(i) for i in range(4)] == [
        "Grating Angle",
        "Collimator Focus",
        "Focus Analysis",
        "DeVeny ETC",
    ]
    window.wavelengthInput.setText("5000")
    window.compute_tilt()
    assert window.grangleResult.text() == "22.35°"
    assert window.demagResult.text() == "2.69 pixels/arcsec"
    window.compute_wavelength()
    assert window.grangleResult.text() == "5000 Å"

    window.temperatureInput.setText("10")
    window.collTiltInput.setText("25")
    window.compute_collfocus()
    assert window.focusResult.text() == "10.2 ± 0.6"
    assert (
        window.startResult.text(),
        window.stepResult.text(),
        window.stepsResult.text(),
    ) == (
        "8.0",
        "0.5",
        "9",
    )


def test_input_validation(window):
    """Invalid values produce a visible error without leaving stale results."""
    window.wavelengthInput.setText("nan")
    window.compute_tilt()
    assert "between 3000 and 11000" in window.grangleResult.text()
    assert not window.demagResult.text()
    window.focusDirectory.setText("/nonexistent")
    window.run_focus()
    assert "existing data directory" in window.focusOutput.toPlainText()
