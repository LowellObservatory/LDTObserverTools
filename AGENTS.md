# Repository Directives

## Python Environment

- Use the `obstools` conda environment for Python quality and test commands when it is available.
- Run focused `pytest` checks for changed behavior and compile touched modules before broader checks.
- Use `black` for formatting and `pylint` for lint checks when available.

## Library and Path Choices

- Prefer `numpy`, `astropy`, and `scipy` for scientific calculations already using those libraries.
- Use `pathlib` for paths and file operations; reserve `os` for operating-system interfaces.
- Keep optional broker and PypeIt imports behind their existing optional dependency boundaries.

## Packaging and Documentation

- Maintain package metadata, dependencies, extras, and entry points in `pyproject.toml`.
- Keep documentation helpers and CI in sync with packaging changes.
- Preserve the MPL-2.0 license and add `SPDX-License-Identifier: MPL-2.0` to new Python files.
- Add complete type hints and NumPy-style docstrings when creating or modifying public functions and methods. Document public constructor parameters on the class.

## Working in This Repository

- Search narrowly with `rg` and inspect only relevant file ranges.
- Do not modify, stage, or discard unrelated untracked observing data and generated outputs.
- Keep tests and verification proportionate to the changed code.
- Review the final diff and commit only task-related files when Git metadata is writable.
