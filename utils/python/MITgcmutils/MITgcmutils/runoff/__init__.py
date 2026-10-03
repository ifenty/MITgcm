"""Sparse runoff forcing files for pkg/exf: schema, integrity checker, examples,
target-table builder, dense-file converter.

* :mod:`.schema`: constants of sparse-runoff schema 1.0 (names, units,
  calendars, tolerances), as specified in ``docs/runoff_schema.md``.
* :func:`check_files` (module :mod:`.check`): validates one file or a set of
  yearly files and returns a :class:`~.check.Report` of rule findings. The
  command line is ``python -m MITgcmutils.runoff.check FILE [FILE ...]
  [--grid-dir DIR] [--strict] [--json OUT] [--tables-only]``.
* :func:`write_example` (module :mod:`.example`): writes a small valid file.
* :func:`build_targets` and :func:`write_targets` (module :mod:`.targets`):
  build the source, alias and target tables from source locations and MITgcm
  grid output (pointwise or spread emission), and write them to a new file or
  into an existing runoff file. The command line is ``python -m
  MITgcmutils.runoff.targets SOURCES --grid-dir DIR -o OUT [options]``.
* :func:`dense_to_sparse` and :func:`sparse_to_dense` (module :mod:`.convert`):
  convert a dense exf ``runoffFile`` (m/s), with the grid and the exf timing
  settings, to a sparse file, and back. The command line is ``python -m
  MITgcmutils.runoff.convert DENSE -o OUT --grid-dir DIR [options]``.

Requires ``netCDF4`` (with ``cftime``), installed by the ``runoff`` extra:
``pip install MITgcmutils[runoff]``. :mod:`.targets` uses ``scipy`` when it
is installed (k-d tree search) and a slower pure-numpy search otherwise.

``check_files``, ``Finding``, ``Report``, ``write_example``,
``build_targets``, ``write_targets``, ``dense_to_sparse`` and
``sparse_to_dense`` are imported lazily, on first access
(module ``__getattr__``, PEP 562). Importing :mod:`.check`, :mod:`.targets`
or :mod:`.convert` eagerly here would put it in ``sys.modules`` before
``python -m MITgcmutils.runoff.check`` (or ``.targets``, ``.convert``) runs it
as ``__main__``, and
runpy would then print a RuntimeWarning about the module being imported twice.
"""

from importlib import import_module

from . import schema

__all__ = ["check_files", "write_example", "build_targets", "write_targets",
           "dense_to_sparse", "sparse_to_dense", "Finding", "Report", "schema"]

#: Lazily imported exports: name -> submodule that defines it.
_LAZY = {"check_files": "check", "Finding": "check", "Report": "check",
         "write_example": "example",
         "build_targets": "targets", "write_targets": "targets",
         "dense_to_sparse": "convert", "sparse_to_dense": "convert"}


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
