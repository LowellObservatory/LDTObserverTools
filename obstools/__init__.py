# -*- coding: utf-8 -*-
# SPDX-License-Identifier: MPL-2.0
#
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
#  Created on 06-Jun-2022
#
#  @author: tbowers

"""Init File"""


# Imports for signal and log handling
import os
from typing import TextIO
import warnings

# Local Imports
from obstools.version import version

# Imports for building help docs
try:
    from obstools import deveny_collfocus
    from obstools import deveny_grangle
    from obstools import dfocus
    from obstools import lmi_etc
    from obstools import fix_ldt_header
    from obstools import neocp_ephem
    from obstools import scrub_deveny_pickup
except ImportError as err:
    pass


def short_warning(
    message: str | Warning,
    category: type[Warning],
    filename: str,
    lineno: int,
    file: TextIO | None = None,
    line: str | None = None,
) -> str:
    """Format a warning with its source filename and line number.

    Parameters
    ----------
    message : str or Warning
        Warning text or instance.
    category : type[Warning]
        Warning category.
    filename : str
        Path to the source of the warning.
    lineno : int
        Line number in the source file.
    file : TextIO, optional
        Unused destination supplied by ``warnings.formatwarning``.
    line : str, optional
        Unused source line supplied by ``warnings.formatwarning``.

    Returns
    -------
    str
        Short formatted warning.
    """
    return f" {category.__name__}: {message} ({os.path.split(filename)[1]}:{lineno})\n"


warnings.formatwarning = short_warning


# Set version
__version__ = version


# Build the list of script classes
def script_classes() -> dict:
    """Build the list of script classes in the package

    Returns
    -------
    :obj:`dict`
        Dictionary of {name:class} for all script classes
    """
    import numpy as np
    from obstools import utils

    # Recursively collect all subclasses
    scr_class = np.array(list(utils.all_subclasses(utils.ScriptBase)))
    scr_name = np.array([c.name() for c in scr_class])
    # Construct a dictionary with the script name and class
    srt = np.argsort(scr_name)
    return dict(zip(scr_name[srt], scr_class[srt]))
