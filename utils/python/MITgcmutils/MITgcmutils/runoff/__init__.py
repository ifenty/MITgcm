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

``check_files``, ``Finding``, ``Report`` and ``write_example`` are imported
lazily, on first access (module ``__getattr__``, PEP 562). Importing
:mod:`.check` eagerly here would put it in ``sys.modules`` before
``python -m MITgcmutils.runoff.check`` runs it as ``__main__``, and runpy would
then print a RuntimeWarning about the module being imported twice.
"""

from importlib import import_module

from . import schema

__all__ = ["check_files", "write_example", "Finding", "Report", "schema"]

#: Lazily imported exports: name -> submodule that defines it.
_LAZY = {"check_files": "check", "Finding": "check", "Report": "check",
         "write_example": "example"}


def __getattr__(name):
    """Import a lazy export from its submodule on first access and cache it."""
    module = _LAZY.get(name)
    if module is None:
        raise AttributeError("module {0!r} has no attribute {1!r}".format(__name__, name))
    value = getattr(import_module("." + module, __name__), name)
    globals()[name] = value
    return value


def __dir__():
    return sorted(set(globals()) | set(_LAZY))
