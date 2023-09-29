# -*- coding: utf-8 -*-
#
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
#  Created on 06-Jun-2022
#
#  @author: tbowers

"""Init File
"""


# Imports for signal and log handling
import os
import warnings

# Local Imports
from .version import version
from obstools import dfocus
from obstools import scrub_deveny_pickup


def short_warning(message, category, filename, lineno, file=None, line=None):
    """
    Return the format for a short warning message.
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
    from pypeit.scripts import scriptbase
    from pypeit.utils import all_subclasses

    # Recursively collect all subclasses
    # Since we use PypeIt's ScriptBase, remove all classes from that package
    scr_class = np.array(
        [
            cls
            for cls in list(all_subclasses(scriptbase.ScriptBase))
            if "pypeit" not in cls.name()
        ]
    )
    scr_name = np.array([c.name() for c in scr_class])
    # Construct a dictionary with the script name and class
    srt = np.argsort(scr_name)
    return dict(zip(scr_name[srt], scr_class[srt]))
