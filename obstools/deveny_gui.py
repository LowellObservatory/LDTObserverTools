# SPDX-License-Identifier: MPL-2.0
"""Tabbed PyQt6 interface for the DeVeny observer tools."""

from __future__ import annotations

import argparse
import pathlib
import sys

import numpy as np
from PyQt6 import QtCore, QtGui, QtWidgets

from obstools import deveny_collfocus, deveny_grangle, utils
from obstools.UI.DevenyToolsWindow import Ui_DevenyToolsWindow


class DevenyToolsWindow(QtWidgets.QMainWindow, Ui_DevenyToolsWindow):
    """Show the DeVeny calculators in one window.

    Parameters
    ----------
    tab : int, optional
        Initially selected tab, in display order. Default is the grating tab.
    """

    def __init__(self, tab: int = 0) -> None:
        super().__init__()
        self.setupUi(self)
        self.LowellLogo.setPixmap(
            QtGui.QPixmap(str(utils.UI / "lowelllogo_horizontal_web.png"))
        )
        self.LowellLogo.setScaledContents(False)
        if sys.platform.startswith("linux"):
            font = QtGui.QFont()
            self.centralwidget.setFont(font)
            font.setPointSize(max(font.pointSize() + 5, 14))
            self.labelTitle.setFont(font)
        self.gratingCombo.setCurrentIndex(1)
        self.tabWidget.setCurrentIndex(tab)
        self.focusDirectory.setText(str(pathlib.Path.cwd()))
        self._focus_process: QtCore.QProcess | None = None

        self.computeTiltButton.clicked.connect(self.compute_tilt)
        self.computeWavelengthButton.clicked.connect(self.compute_wavelength)
        self.computeFocusButton.clicked.connect(self.compute_collfocus)
        self.brokerButton.clicked.connect(self.get_broker_values)
        self.browseDirectoryButton.clicked.connect(self.browse_directory)
        self.runFocusButton.clicked.connect(self.run_focus)

    def _grooves(self) -> float:
        """Return the line density of the selected grating."""
        return float(self.gratingCombo.currentText().split(" - ")[1].split(" g/mm")[0])

    def compute_tilt(self) -> None:
        """Compute grating tilt and slit demagnification from wavelength."""
        try:
            wavelength = float(self.wavelengthInput.text())
            if not 3000 <= wavelength <= 11000 or not np.isfinite(wavelength):
                raise ValueError("Wavelength must be between 3000 and 11000 Å.")
            angle, demagnification = deveny_grangle.compute_grangle(
                self._grooves(), wavelength
            )
        except (ValueError, RuntimeError) as error:
            self.grangleResult.setText(str(error))
            self.demagResult.clear()
            return
        self.tiltInput.setText(f"{angle + deveny_grangle.TGOFFSET:.2f}")
        self.grangleResult.setText(f"{angle + deveny_grangle.TGOFFSET:.2f}°")
        self.demagResult.setText(
            f"{deveny_grangle.PIXSCALE * demagnification:.2f} pixels/arcsec"
        )

    def compute_wavelength(self) -> None:
        """Compute central wavelength from grating tilt."""
        try:
            angle = float(self.tiltInput.text())
            if not 0 <= angle <= 48 or not np.isfinite(angle):
                raise ValueError("Grating tilt must be between 0 and 48°.")
            wavelength = deveny_grangle.lambda_at_angle(angle, self._grooves())
            demagnification = deveny_grangle.deveny_amag(angle)
        except ValueError as error:
            self.grangleResult.setText(str(error))
            self.demagResult.clear()
            return
        self.wavelengthInput.setText(f"{wavelength:.0f}")
        self.grangleResult.setText(f"{wavelength:.0f} Å")
        self.demagResult.setText(
            f"{deveny_grangle.PIXSCALE * demagnification:.2f} pixels/arcsec"
        )

    def compute_collfocus(self) -> None:
        """Estimate focus and the LOUI focus sequence."""
        try:
            temperature = float(self.temperatureInput.text())
            tilt = float(self.collTiltInput.text())
            if not np.isfinite(temperature) or not np.isfinite(tilt):
                raise ValueError("Enter finite numeric temperature and tilt values.")
        except ValueError:
            self.focusResult.setText("Enter numeric temperature and tilt values.")
            for label in (self.startResult, self.stepResult, self.stepsResult):
                label.clear()
            return
        focus = deveny_collfocus.calculate_collimator_focus(
            temperature, tilt, self.filterCombo.currentText() != "CLEAR"
        )
        start, step, count = deveny_collfocus.calculate_focus_sequence(focus)
        self.focusResult.setText(f"{focus:.1f} ± 0.6")
        self.startResult.setText(f"{start:.1f}")
        self.stepResult.setText(f"{step:.1f}")
        self.stepsResult.setText(str(count))

    def get_broker_values(self) -> None:
        """Load current telemetry when the optional ActiveMQ service is configured."""
        config = utils.CONFIG / "activemq_joe.yaml"
        if not config.is_file():
            self.focusResult.setText("Broker configuration is unavailable.")
            return
        try:
            from obstools import broker_listener

            listener = broker_listener.ActiveMQ_Listener(config)
            temperature = deveny_collfocus.extract_broker_values(
                listener.mounttemp_from_broker
            )
            tilt = deveny_collfocus.extract_broker_values(listener.grangle_from_broker)
        except (ImportError, OSError, AttributeError, RuntimeError) as error:
            self.focusResult.setText(f"Broker unavailable: {error}")
            return
        self.temperatureInput.setText(temperature)
        self.collTiltInput.setText(tilt)

    def browse_directory(self) -> None:
        """Choose the directory containing a DeVeny focus run."""
        directory = QtWidgets.QFileDialog.getExistingDirectory(
            self, "Choose focus data directory", self.focusDirectory.text()
        )
        if directory:
            self.focusDirectory.setText(directory)

    def run_focus(self) -> None:
        """Run the existing focus analysis without blocking the Qt event loop."""
        directory = pathlib.Path(self.focusDirectory.text()).expanduser()
        if not directory.is_dir():
            self.focusOutput.setPlainText("Choose an existing data directory.")
            return
        try:
            threshold = float(self.focusThreshold.text())
            if not np.isfinite(threshold) or threshold <= 0:
                raise ValueError
        except ValueError:
            self.focusOutput.setPlainText("Enter a positive numeric line threshold.")
            return
        focus_log = self.focusLog.text().strip() or "last"
        if focus_log != "last" and pathlib.Path(focus_log).name != focus_log:
            self.focusOutput.setPlainText("Choose a focus log in the data directory.")
            return
        self.focusOutput.clear()
        self.runFocusButton.setEnabled(False)
        self._focus_process = QtCore.QProcess(self)
        self._focus_process.setWorkingDirectory(str(directory))
        self._focus_process.setProcessChannelMode(
            QtCore.QProcess.ProcessChannelMode.MergedChannels
        )
        self._focus_process.readyReadStandardOutput.connect(self._read_focus_output)
        self._focus_process.finished.connect(self._focus_finished)
        self._focus_process.errorOccurred.connect(self._focus_error)
        self._focus_process.start(
            sys.executable,
            [
                "-c",
                "from obstools.dfocus import DFocus; DFocus.entry_point()",
                "--flog",
                focus_log,
                "--thresh",
                str(threshold),
                "--nodisplay",
            ],
        )

    def _read_focus_output(self) -> None:
        """Append the focus process output to the GUI."""
        if self._focus_process is not None:
            data = bytes(self._focus_process.readAllStandardOutput()).decode(
                "utf-8", errors="replace"
            )
            self.focusOutput.insertPlainText(data)

    def _focus_finished(
        self, exit_code: int, _status: QtCore.QProcess.ExitStatus
    ) -> None:
        """Report completion of the focus process."""
        self._read_focus_output()
        self.runFocusButton.setEnabled(True)
        self.focusOutput.appendPlainText(
            "Focus analysis complete."
            if exit_code == 0
            else f"Focus analysis failed ({exit_code})."
        )
        self._focus_process = None

    def _focus_error(self, _error: QtCore.QProcess.ProcessError) -> None:
        """Show an operating system error from the focus process."""
        if self._focus_process is not None:
            self.focusOutput.appendPlainText(self._focus_process.errorString())
        self.runFocusButton.setEnabled(True)

    def closeEvent(self, event: QtGui.QCloseEvent) -> None:
        """Stop a running focus process before the window closes."""
        if (
            self._focus_process is not None
            and self._focus_process.state() != QtCore.QProcess.ProcessState.NotRunning
        ):
            self._focus_process.kill()
            self._focus_process.waitForFinished(2000)
        super().closeEvent(event)


def run_gui(tab: int = 0) -> None:
    """Launch the DeVeny Tools GUI.

    Parameters
    ----------
    tab : int, optional
        Initially selected tab, in display order.
    """
    app = QtWidgets.QApplication.instance()
    owns_app = app is None
    if owns_app:
        app = QtWidgets.QApplication(sys.argv)
    window = DevenyToolsWindow(tab=tab)
    window.show()
    if owns_app:
        app.exec()


class DevenyGui(utils.ScriptBase):
    """Console entry point for the combined DeVeny GUI."""

    @classmethod
    def get_parser(
        cls,
        description: str | None = None,
        width: int | None = None,
        formatter: type[
            argparse.HelpFormatter
        ] = argparse.ArgumentDefaultsHelpFormatter,
    ) -> argparse.ArgumentParser:
        """Construct the standard command-line parser.

        Parameters
        ----------
        description : str, optional
            Description override (unused).
        width : int, optional
            Maximum help width.
        formatter : type[argparse.HelpFormatter], optional
            Help formatter class.

        Returns
        -------
        argparse.ArgumentParser
            Parser for the application entry point.
        """
        return super().get_parser(
            description="DeVeny Tools GUI", width=width, formatter=formatter
        )

    @staticmethod
    def main(args: argparse.Namespace) -> None:
        """Launch the tabbed application.

        Parameters
        ----------
        args : argparse.Namespace
            Parsed command-line arguments (unused).
        """
        run_gui()
