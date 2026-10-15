# SPDX-License-Identifier: MPL-2.0
"""
Construct an rst table with the dependencies

The table is generated from the package metadata without optional dependencies.
"""

# Built-In Libraries
import pathlib
import tomllib

# The repository root is up two levels from here
# NOTE: This is a hack needed for this script to run on GH pages.
OBSTOOLS_ROOT = pathlib.Path(__file__).parents[2]


def write_dependency_table(pyproject_file: pathlib.Path, path: pathlib.Path) -> None:
    """Write a reStructuredText table of package dependencies.

    Parameters
    ----------
    pyproject_file : pathlib.Path
        Path to the project's ``pyproject.toml`` file.
    path : pathlib.Path
        Directory in which to write the dependency table.
    """
    ofile = path / "dependencies_table.rst"

    with pyproject_file.open("rb") as stream:
        project = tomllib.load(stream)["project"]

    user_requires = sorted(project["dependencies"])
    dev_requires = sorted(project["optional-dependencies"]["dev"])
    pypeit_requires = sorted(project["optional-dependencies"]["pypeit"])
    required_python = project["requires-python"]

    rows = [
        ("Python Version", f"``{required_python}``"),
        ("Required for users", ", ".join(f"``{item}``" for item in user_requires)),
        (
            "Optional ``pypeit`` requirements",
            ", ".join(f"``{item}``" for item in pypeit_requires),
        ),
        ("Required for developers", ", ".join(f"``{item}``" for item in dev_requires)),
    ]

    label_width = max(len(label) for label, _ in rows)
    value_width = max(len(value) for _, value in rows)
    separator = f"{'=' * label_width}  {'=' * value_width}"
    lines = [separator]
    lines.extend(f"{label:<{label_width}}  {value}" for label, value in rows)
    lines.append(separator)
    ofile.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote: {ofile}")


def main() -> None:
    """Generate the dependency table in the documentation directory."""
    output_root = OBSTOOLS_ROOT / "doc" / "include"
    if not output_root.is_dir():
        raise NotADirectoryError(f"{output_root} does not exist!")

    pyproject_file = OBSTOOLS_ROOT / "pyproject.toml"
    if not pyproject_file.is_file():
        raise FileNotFoundError(f"{pyproject_file} does not exist!")
    write_dependency_table(pyproject_file, output_root)


if __name__ == "__main__":
    main()
