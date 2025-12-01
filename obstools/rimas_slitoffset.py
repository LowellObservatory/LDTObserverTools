# -*- coding: utf-8 -*-
#
#  This file is part of LDTObserverTools.
#
#   This Source Code Form is subject to the terms of the Mozilla Public
#   License, v. 2.0. If a copy of the MPL was not distributed with this
#   file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
#  Created on 01-Dec-2025
#
#  @author: tbowers
# pylint: disable=c-extension-no-member

"""RIMAS offset-to-slit calculation GUI

LDTObserverTools contains python ports of various LDT Observer Tools

Lowell Discovery Telescope (Lowell Observatory: Flagstaff, AZ)
https://lowell.edu
"""

# Built-In Libraries
import argparse
import dataclasses
import pathlib
import re
import sys

# 3rd-Party Libraries
import astropy.table
import numpy as np
from PyQt6 import QtGui
from PyQt6 import QtWidgets

# Local Libraries
from obstools.UI.RIMAS2SlitMainWindow import Ui_MainWindow
from obstools import utils

# GUI Classes ================================================================#
class RIMASWindow(utils.ObstoolsGUI, Ui_MainWindow):
    """Exposure Time Calculator Main Window Class

    The UI is defined in ETCWindow.ui and translated (via pyuic6) into python
    in ETCWindow.py.  This class inherits the UI and defines the various
    actions needed to compute ETCs from the GUI inputs.
    """

    def __init__(self):
        super().__init__()
        self.setupUi(self)
        self.show()


        # Connect buttons to actions
        self.exitButton.pressed.connect(self.exit_button_clicked)

        self.set_fonts_and_logo()

# Command Line Script Infrastructure (borrowed from PypeIt) ==================#
class RimasSlitoffset(utils.ScriptBase):
    """Script class for ``rimas_slitoffset`` tool

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
            description="RIMAS Slit Offset Calculator", width=width
        )
        return parser

    @staticmethod
    def main(args):
        """Main Driver

        Set up the top-level PyQt6 objects and start the event loop
        """
        # Create the QApplication object and main Qt window
        app = QtWidgets.QApplication([])
        _ = RIMASWindow()

        # Giddy up!
        app.exec()
