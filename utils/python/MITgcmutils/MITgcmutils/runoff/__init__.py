"""Sparse runoff forcing files for pkg/exf: schema, integrity checker, examples.

* :mod:`.schema`: constants of sparse-runoff schema 1.0 (names, units,
  calendars, tolerances), as specified in ``docs/runoff_schema.md``.
* :func:`check_files` (module :mod:`.check`): validates one file or a set of
  yearly files and returns a :class:`~.check.Report` of rule findings. The
  command line is ``python -m MITgcmutils.runoff.check FILE [FILE ...]
  [--grid-dir DIR] [--strict] [--json OUT]``.
* :func:`write_example` (module :mod:`.example`): writes a small valid file.

Requires ``netCDF4`` (with ``cftime``), installed by the ``runoff`` extra:
``pip install MITgcmutils[runoff]``.
"""

from .check import Finding, Report, check_files
from .example import write_example
from . import schema

__all__ = ["check_files", "write_example", "Finding", "Report", "schema"]
