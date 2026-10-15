"""
Construct an rst table with the dependencies

This module was borrowed and updated from PyepIt
    https://github.com/pypeit/PypeIt/
"""

# Built-In Libraries
import pathlib
import tomllib

# 3rd Party Libraries
import numpy as np
from pypeit.utils import string_table

# The repository root is up two levels from here
# NOTE: This is a hack needed for this script to run on GH pages.
OBSTOOLS_ROOT = pathlib.Path(__file__).parents[2]


def write_dependency_table(setup_file: pathlib.Path, path: pathlib.Path):
    """Write the package dependency table

    Parameters
    ----------
    setup_file : :obj:`~pathlib.Path`
        The name and location of pyproject.toml
    path : :obj:`~pathlib.Path`
        The path within the documentation of where to write this table
    """
    ofile = path / "dependencies_table.rst"

    with setup_file.open("rb") as stream:
        project = tomllib.load(stream)["project"]

    user_requires = np.sort(project["dependencies"])
    dev_requires = np.sort(project["optional-dependencies"]["dev"])
    pypeit_requires = np.sort(project["optional-dependencies"]["pypeit"])
    required_python = project["requires-python"]

    data_table = np.empty((4, 2), dtype=object)
    data_table[0, :] = ["Python Version", f"``{required_python}``"]
    data_table[1, :] = [
        "Required for users",
        ", ".join([f"``{u}``" for u in user_requires]),
    ]
    data_table[2, :] = [
        "Optional ``pypeit`` requirements",
        ", ".join([f"``{u}``" for u in pypeit_requires]),
    ]
    data_table[3, :] = [
        "Required for developers",
        ", ".join([f"``{d}``" for d in dev_requires]),
    ]

    lines = string_table(data_table, delimeter="rst", has_header=False)
    with open(ofile, "w", encoding="utf-8") as f_obj:
        f_obj.write("\n".join(line.rstrip() for line in lines.splitlines()) + "\n")
    print(f"Wrote: {ofile}")


def main():
    """Driver"""
    output_root = OBSTOOLS_ROOT / "doc" / "include"
    if not output_root.is_dir():
        raise NotADirectoryError(f"{output_root} does not exist!")

    setup_file = OBSTOOLS_ROOT / "pyproject.toml"
    if not setup_file.is_file():
        raise FileNotFoundError(f"{setup_file} does not exist!")
    write_dependency_table(setup_file, output_root)


if __name__ == "__main__":
    main()
