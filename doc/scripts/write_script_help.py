"""
Dynamically build the rst documentation with the script help text.

This module was borrowed and updated from PyepIt
    https://github.com/pypeit/PypeIt/
"""

# Built-In Libraries
import os
import pathlib
import time

# 3rd Party Libraries

# Internal Imports
from obstools import script_classes
from obstools import utils

# The repository root is up two levels from here
# NOTE: This is a hack needed for this script to run on GH pages.
OBSTOOLS_ROOT = pathlib.Path(__file__).parents[2]

# -----------------------------------------------------------------------------


def write_help(script_cls: utils.ScriptBase, opath: pathlib.Path, width: int = 80):
    """Write the ``.rst`` help files for all scripts

    Parameters
    ----------
    script_cls : :class:`~pypeit.scripts.scriptbase.ScriptBase`
        The :class:`~pypeit.scripts.scriptbase.ScriptBase` subclass
    opath : :obj:`~pathlib.Path`
        The path into which to place the resulting ``.rst`` file
    width : :obj:`int`, optional
        The width of the help text before wrapping  (Default: 80)
    """
    exe = script_cls.name
    ofile = os.path.join(opath, f"{exe}.rst")
    lines = [".. code-block:: console", ""]
    lines += [f"    $ {exe} -h"]
    parser = script_cls.get_parser(width=width)
    parser.prog = exe
    lines += ["    " + l for l in parser.format_help().split("\n")]
    # Remove ".dev....." junk from the version
    lines[-2] = lines[-2].split('.dev')[0]
    print(f"Writing: {ofile}")
    with open(ofile, "w", encoding="utf-8") as f_obj:
        f_obj.write("\n".join(lines))


if __name__ == "__main__":
    t = time.perf_counter()

    # Define (and create if nonexistant) the directory into which to write
    path = OBSTOOLS_ROOT / "doc" / "help"
    path.mkdir(parents=True, exist_ok=True)

    # Get the list of script names and script classes
    for scr_class in script_classes().values():
        write_help(scr_class, path)

    print(f"Elapsed time: {time.perf_counter() - t} seconds")
