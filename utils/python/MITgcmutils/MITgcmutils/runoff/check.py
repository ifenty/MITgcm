"""Integrity checker for sparse-runoff NetCDF files (schema 1.0).

Implements every rule of section 9 of ``docs/runoff_schema.md`` (project
repository). Each rule produces :class:`Finding` objects with the rule id, a
level (``E`` error, ``W`` warning, ``I`` information) and a message that
names the file, the variable and, where applicable, the ``source_id`` and the
time of the offending record.

Library use::

    from MITgcmutils.runoff import check_files
    report = check_files(["runoff_2000.nc", "runoff_2001.nc"], grid_dir="run/")
    report.errors          # list of E findings

Command line::

    python -m MITgcmutils.runoff.check FILE [FILE ...] [--grid-dir DIR]
                                       [--strict] [--json OUT] [--tables-only]

Exit status: 0 no errors, 1 errors (or warnings with ``--strict``), 2 usage
or I/O problem (unreadable file, missing grid file, unwritable JSON).

Implementation notes that the schema leaves to the checker:

* Tables-only mode (``tables_only=True``, ``--tables-only``) checks a file
  that holds only the source, alias and target tables, such as the output of
  :mod:`MITgcmutils.runoff.targets` before its time series are added. It does
  not require the ``time`` dimension (S03) or the ``time`` and ``runoff_flux``
  variables (S04), and it skips the time and time-series rules M01-M06,
  D01-D09 and P01 entirely, even for time variables that are present. With
  several files it also skips the X01 time order and continuity checks; the
  X01 table comparison still runs. Every other rule runs unchanged. One
  ``S10`` (I) finding per file lists what was skipped
  (:data:`schema.TABLES_ONLY_SKIPPED`), and ``Report.stats["tables_only"]``
  records the same list.

* All variables are read with ``set_auto_maskandscale(False)`` and
  ``set_auto_chartostring(False)``: netCDF4 neither masks nor unpacks, so the
  checker sees the stored values, exactly as the Fortran reader does.
  ``scale_factor``/``add_offset`` on a model-read variable (``time``,
  ``time_bnds``, ``source_id``, ``target_source``, ``target_cell``,
  ``target_fraction``, ``target_level``, ``target_cell_area``, ``runoff_*``) is
  an error (S07), and fill and range checks then apply to the stored values.
  ``target_lon``/``target_lat`` and user variables may be packed. A value counts as
  missing when it is NaN, equals the variable's ``_FillValue`` (or, without
  that attribute, the netCDF default fill of the variable's type, which is
  what unwritten records read as), or equals a ``missing_value``.
* S08 covers only the text attributes the model reads (section 1:
  :data:`schema.MODEL_READ_GLOBAL_TEXT_ATTRS`, ``units``/``calendar`` of
  ``time`` and ``time_bnds``, ``units`` of ``runoff_*``): each must be ASCII
  and stored as NC_CHAR. Descriptive ``mitgcm_grid_name``/
  ``mitgcm_grid_description`` may be NC_STRING and UTF-8. The stored type,
  which netCDF4-python does not expose, is read with ``nc_inq_atttype`` from
  the libnetcdf that netCDF4 loaded, through :mod:`ctypes`, on a second
  read-only open of the file. If that library can't be loaded, the types are
  not checked and S08 is reported as a W finding instead (exit status 1 only
  with ``--strict``); the ASCII check still runs.
* S09: ``_FillValue``/``missing_value`` on a numeric model-read variable must
  be a number of the variable's own type. Text values are never used as fill
  values by the other rules. ``source_id`` is text, so it is not checked.
* P02 reads ``netCDF4.Variable.filters()``; only deflate (``zlib``), shuffle
  and fletcher32 are allowed on model-read variables. Filters that
  ``filters()`` doesn't report (unknown plugins) are not detected. If a time
  series with a disallowed filter can't be read here, that is reported under
  P02 instead of stopping the check.
* A month or year edge is recognized within the time tolerance: a bound up to
  :data:`schema.TIME_EQUAL_TOL_SECONDS` before an edge counts as that edge
  (M05).
* Variable-length strings: an unwritten element reads as ``""``, the NC_STRING
  default fill, so an empty string counts as missing. Char arrays are decoded
  after removing trailing NUL padding; ``source_id`` also has trailing blanks
  removed (Fortran blank padding) before its pattern is checked.
* Time series are read in blocks of whole records of at most
  ``max_block_bytes`` (default 50 MiB), so memory stays bounded; static tables
  and the ``time`` axis are read whole.
* Calendar arithmetic (monthly, yearly, annual repeat, yearly file names, file
  ordering) uses :mod:`cftime`. Two times are equal within
  :data:`schema.TIME_EQUAL_TOL_SECONDS` everywhere. ``calendar`` values are
  matched ignoring letter case; across files (X01) calendars compare by their
  MITgcm mapping, so ``standard`` and ``gregorian`` files can be combined.
* ``runoff_temperature`` fill and NaN mean "surface temperature" and are not
  range-checked (D04); ±Inf is never a missing-value marker: every ±Inf value,
  and a ``_FillValue``/``missing_value`` of ±Inf, is an error (D09).
* X01 continuity: each file's first ``time_bnds`` start must equal the previous
  file's last end, so every file of a multi-file set needs ``time_bnds``.
  With ``fixed`` sampling the spacing across a file boundary must equal
  ``mitgcm_time_period``. ``_YYYY`` files must also share one offset of their
  first ``time`` value from 1 January of the year in their name.
* R03 compares unpacked ``target_lon``/``target_lat``, which may be packed
  because the model doesn't read them.
* U01 requires the units listed in :data:`schema.TABLE_UNITS` on those schema
  variables when present, and forbids units on index variables; user-added
  ``*_lon``/``*_lat`` variables are not checked.
* Rules with a potentially large number of occurrences list the first
  :data:`MAX_DETAILS` per rule and variable, then one summary finding with the
  total count.
"""

import argparse
import glob
import json
import os
import sys
from dataclasses import asdict, dataclass, field

import numpy as np

from . import schema as S

__all__ = ["Finding", "Report", "CheckIOError", "check_files", "main",
           "DEFAULT_BLOCK_BYTES", "MAX_DETAILS"]

#: Upper bound on the bytes of one block of time-series records.
DEFAULT_BLOCK_BYTES = 50 * 2**20
#: Findings listed individually per (rule, variable) before summarizing.
MAX_DETAILS = 20


class CheckIOError(OSError):
    """A usage or I/O problem (exit status 2), not a finding about the file."""


@dataclass
class Finding:
    """One rule violation or observation.

    ``source_id``/``source_index`` identify the source, and ``record``/``time``
    the time record (value in the file's time units), where they apply.
    """
    rule: str
    level: str
    message: str
    file: str = None
    variable: str = None
    source_id: str = None
    source_index: int = None
    record: int = None
    time: float = None

    def to_dict(self):
        return asdict(self)


@dataclass
class Report:
    """Result of :func:`check_files`: every finding, plus reading statistics.

    ``stats["blocks_read"]`` maps ``"<file>::<variable>"`` to the number of
    record blocks read for that time series.
    """
    files: list
    grid_dir: str = None
    findings: list = field(default_factory=list)
    stats: dict = field(default_factory=dict)

    def _by_level(self, level):
        return [f for f in self.findings if f.level == level]

    @property
    def errors(self):
        return self._by_level("E")

    @property
    def warnings(self):
        return self._by_level("W")

    @property
    def infos(self):
        return self._by_level("I")

    def rules(self, level=None):
        """Set of rule ids that fired, optionally only at ``level``."""
        return {f.rule for f in self.findings if level is None or f.level == level}

    def exit_code(self, strict=False):
        """0 when acceptable, 1 on errors (or on warnings when ``strict``)."""
        if self.errors or (strict and self.warnings):
            return 1
        return 0

    def to_dict(self):
        return {
            "schema_version": S.SCHEMA_VERSION,
            "files": [str(p) for p in self.files],
            "grid_dir": None if self.grid_dir is None else str(self.grid_dir),
            "counts": {"E": len(self.errors), "W": len(self.warnings),
                       "I": len(self.infos)},
            "findings": [f.to_dict() for f in self.findings],
            "stats": self.stats,
        }

    def format_text(self):
        lines = ["[{0}] {1} {2}".format(f.level, f.rule, f.message)
                 for f in self.findings]
        lines.append("{0} error(s), {1} warning(s), {2} info".format(
            len(self.errors), len(self.warnings), len(self.infos)))
        return "\n".join(lines)

    def _add(self, rule, message, file=None, variable=None):
        self.findings.append(Finding(rule, S.RULES[rule], message,
                                     None if file is None else str(file), variable))


# ---------------------------------------------------------------------------
# Small helpers


def _attr(obj, name):
    """Attribute value, or None when absent."""
    return obj.getncattr(name) if name in obj.ncattrs() else None


def _str_attr(value):
    """A string attribute as ``str``, or None when it isn't a single string."""
    if isinstance(value, str):
        return value
    if isinstance(value, bytes):
        try:
            return value.decode("utf-8")
        except UnicodeDecodeError:
            return None
    return None


def _scalar(value):
    """``(python value, numpy kind)`` of a one-element numeric attribute."""
    if value is None or isinstance(value, (str, bytes)):
        return None, "S"
    a = np.asarray(value)
    if a.size != 1:
        return None, a.dtype.kind
    return a.reshape(()).item(), a.dtype.kind


def _numeric_values(value):
    """The numbers of a numeric attribute value, or ``[]`` for text (S09)."""
    if value is None or isinstance(value, (str, bytes)):
        return []
    a = np.asarray(value)
    if a.dtype.kind not in "fiu":
        return []
    return np.ravel(a).tolist()


def _is_ascii(value):
    """False when a text attribute value has a non-ASCII character (S08)."""
    if isinstance(value, str):
        return value.isascii()
    if isinstance(value, bytes):
        return all(c < 128 for c in value)
    if isinstance(value, (list, tuple, np.ndarray)) and np.asarray(value).dtype.kind in "OUS":
        return all(_is_ascii(v) for v in np.ravel(np.asarray(value, dtype=object)))
    return True


def _unpacked(var):
    """float64 values of ``var`` unpacked as CF does (``packed * scale_factor +
    add_offset``), with fill and missing values as NaN. Used for variables the
    model doesn't read, which may be packed (R03)."""
    raw = np.asarray(var[:])
    a = raw.astype(np.float64)
    if raw.dtype.kind in "fiu":
        a[_fill_mask(raw, _fill_values(var))] = np.nan
    sf, _ = _scalar(_attr(var, "scale_factor"))
    ao, _ = _scalar(_attr(var, "add_offset"))
    if sf is not None:
        a = a * float(sf)
    if ao is not None:
        a = a + float(ao)
    return a


def _type_class(var):
    """``float``, ``int``, ``char``, ``vstring`` or ``other``."""
    dt = var.dtype
    if dt is str:
        return "vstring"
    if getattr(var, "_isvlen", False) or getattr(var, "_iscompound", False):
        return "other"
    try:
        dt = np.dtype(dt)
    except TypeError:
        return "other"
    if dt.kind == "f":
        return "float"
    if dt.kind in "iu":
        return "int"
    if dt.kind == "S" and dt.itemsize == 1:
        return "char"
    return "other"


def _describe_type(var):
    tc = _type_class(var)
    return "string" if tc == "vstring" else "{0} ({1})".format(tc, var.dtype)


def _matches(var, dims, cls):
    tc, vd = _type_class(var), tuple(var.dimensions)
    if cls == "string":
        return ((tc == "vstring" and vd == dims)
                or (tc == "char" and len(vd) == len(dims) + 1 and vd[:-1] == dims))
    if cls == "char":
        return tc == "char" and len(vd) == len(dims) and vd[:-1] == dims[:-1]
    if cls == "double":
        return tc == "float" and np.dtype(var.dtype) == np.float64 and vd == dims
    return tc == cls and vd == dims


def _fill_values(var):
    """Values that mean "missing" for ``var`` (see module notes)."""
    import netCDF4
    try:
        dt = np.dtype(var.dtype)
    except TypeError:
        return []
    if dt.kind not in "fiu":
        return []
    names = var.ncattrs()
    vals = []
    # Text values are not numbers the reader could match; S09 reports them.
    if "_FillValue" in names:
        vals.extend(_numeric_values(var.getncattr("_FillValue")))
    else:
        key = dt.str[1:]
        if key in netCDF4.default_fillvals:
            vals.append(netCDF4.default_fillvals[key])
    if "missing_value" in names:
        vals.extend(_numeric_values(var.getncattr("missing_value")))
    return vals


def _fill_mask(a, fills):
    """Boolean mask of ``a`` equal to any of ``fills`` (compared in a's type)."""
    if not fills:
        return np.zeros(np.shape(a), dtype=bool)
    with np.errstate(over="ignore", invalid="ignore"):
        f = np.array([np.asarray(v).astype(a.dtype) for v in fills], dtype=a.dtype)
    return np.isin(a, f)


def _read_strings(var):
    """Decode a string variable to ``(list of str, missing mask)``."""
    tc = _type_class(var)
    if tc == "vstring":
        raw = np.asarray(var[:], dtype=object).ravel()
        fill = _str_attr(_attr(var, "_FillValue"))
        out = ["" if v is None else str(v) for v in raw]
        missing = np.array([s == "" or (fill is not None and s == fill) for s in out],
                           dtype=bool)
        return out, missing
    raw = np.asarray(var[:])
    out = []
    for row in raw.reshape(raw.shape[0], -1) if raw.ndim else []:
        out.append(row.tobytes().rstrip(b"\x00").decode("utf-8", errors="replace"))
    missing = np.array([s == "" for s in out], dtype=bool)
    return out, missing


def _value_label(v, is_fill):
    if np.isnan(v):
        return "NaN"
    if np.isinf(v):
        return "+Inf" if v > 0 else "-Inf"
    if is_fill:
        return "fill/missing value {0:.9g}".format(v)
    return "{0:.9g}".format(v)


def _model_read(name):
    """True for variables the model reads (schema section 1)."""
    return name in S.MODEL_READ_VARIABLES or name.startswith(S.MODEL_READ_PREFIXES)


_LIBNETCDF = []   # cache: [ctypes library or None]


def _libnetcdf():
    """The libnetcdf shared library netCDF4 uses, via ctypes, or None.

    Prefers the copy already mapped into this process (Linux ``/proc``), then
    ``ctypes.util.find_library("netcdf")``.
    """
    if _LIBNETCDF:
        return _LIBNETCDF[0]
    import ctypes
    import ctypes.util
    import netCDF4  # noqa: F401 - make sure netCDF4's libnetcdf is loaded
    candidates = []
    try:
        with open("/proc/self/maps") as f:
            candidates += sorted({ln.split()[-1] for ln in f
                                  if "libnetcdf" in ln and ln.split()[-1].startswith("/")})
    except OSError:
        pass
    found = ctypes.util.find_library("netcdf")
    if found:
        candidates.append(found)
    lib = None
    for path in candidates:
        try:
            lib = ctypes.CDLL(path)
            for fn in ("nc_open", "nc_close", "nc_inq_varid", "nc_inq_atttype"):
                getattr(lib, fn)
            break
        except (OSError, AttributeError):
            lib = None
    _LIBNETCDF.append(lib)
    return lib


def _attr_types(path, queries):
    """Stored netCDF types of attributes, ``{(variable or None, name): nc_type}``.

    Opens ``path`` read-only a second time through libnetcdf. Returns None when
    the library can't be loaded (including a loader failure) or the file can't
    be opened this way. Attributes that don't exist are left out.
    """
    import ctypes
    try:
        lib = _libnetcdf()
    except Exception:  # noqa: BLE001 - any loader failure means "not checked"
        lib = None
    if lib is None:
        return None
    ncid = ctypes.c_int()
    if lib.nc_open(os.fsencode(path), 0, ctypes.byref(ncid)) != 0:   # NC_NOWRITE
        return None
    out = {}
    try:
        for var, att in queries:
            varid = ctypes.c_int(-1)                                   # NC_GLOBAL
            if var is not None and lib.nc_inq_varid(
                    ncid, var.encode("utf-8"), ctypes.byref(varid)) != 0:
                continue
            xtype = ctypes.c_int()
            if lib.nc_inq_atttype(ncid, varid, att.encode("utf-8"),
                                  ctypes.byref(xtype)) == 0:
                out[(var, att)] = xtype.value
    finally:
        lib.nc_close(ncid)
    return out


def _check_attr_types(ctx):
    """S08: text attributes the model reads (section 1) are ASCII NC_CHAR.

    Only :data:`schema.MODEL_READ_GLOBAL_TEXT_ATTRS`, ``units``/``calendar`` of
    ``time`` and ``time_bnds`` and ``units`` of ``runoff_*`` variables are
    checked; descriptive ``mitgcm_grid_*`` text may be NC_STRING and UTF-8.
    The ASCII check needs no libnetcdf; the NC_STRING check does.
    """
    ds = ctx.ds
    queries = [(None, a) for a in S.MODEL_READ_GLOBAL_TEXT_ATTRS if a in ds.ncattrs()]
    for name in ("time", "time_bnds"):
        if name in ds.variables:
            queries += [(name, a) for a in S.MODEL_READ_TIME_ATTRS
                        if a in ds.variables[name].ncattrs()]
    for name, var in ds.variables.items():
        if name.startswith("runoff_") and "units" in var.ncattrs():
            queries.append((name, "units"))
    for var, att in queries:
        value = (ds if var is None else ds.variables[var]).getncattr(att)
        if not _is_ascii(value):
            what = "global attribute" if var is None else "attribute"
            ctx.add("S08", "{0} '{1}' = {2!r} is not ASCII; the model reads it as "
                    "ASCII text (NF_GET_ATT_TEXT)".format(what, att, value), var)
    types = _attr_types(ctx.path, queries)
    if types is None:
        ctx.add("S08", "attribute types not checked: libnetcdf could not be loaded "
                "through ctypes to read them, so NC_STRING model-read text "
                "attributes can't be detected", level="W")
        return
    for (var, att), xtype in types.items():
        if xtype == S.NC_STRING:
            what = "global attribute" if var is None else "attribute"
            ctx.add("S08", "{0} '{1}' is stored as NC_STRING; the Fortran reader "
                    "(NF_GET_ATT_TEXT) needs a char (NC_CHAR) attribute".format(what, att),
                    var)


def _days_in_month(year, month, calendar):
    """Length of a month in one of the allowed CF calendars."""
    if calendar == "360_day":
        return 30
    days = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)[month - 1]
    if month == 2 and calendar not in ("noleap", "365_day"):
        if calendar in ("standard", "gregorian") and year < 1583:
            days += year % 4 == 0          # Julian part of the mixed calendar
        else:
            days += year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)
    return days


class _TimeAxis:
    """CF time units and calendar, with cftime conversion helpers."""

    def __init__(self, units, match, calendar):
        import cftime
        self.cftime = cftime
        self.units = units.strip()
        unit, y, mo, d, hh, mi, ss, _z = match.groups()
        self.factor = S.TIME_UNIT_SECONDS[unit]
        self.calendar = calendar
        # Validate the reference date in this calendar; raises ValueError.
        y, mo, d = int(y), int(mo), int(d)
        hh, mi, ss = int(hh or 0), int(mi or 0), int(ss or 0)
        if not (1 <= mo <= 12 and 1 <= d <= _days_in_month(y, mo, calendar)
                and hh < 24 and mi < 60 and ss < 60):
            raise ValueError("{0:04d}-{1:02d}-{2:02d} {3:02d}:{4:02d}:{5:02d} is not a "
                             "valid date".format(y, mo, d, hh, mi, ss))
        cftime.datetime(y, mo, d, hh, mi, ss, calendar=calendar)
        self.cf_units = "{0} since {1:04d}-{2:02d}-{3:02d} {4:02d}:{5:02d}:{6:02d}".format(
            unit, y, mo, d, hh, mi, ss)
        self.tol = S.TIME_EQUAL_TOL_SECONDS / self.factor

    def num2date(self, x):
        return self.cftime.num2date(float(x), self.cf_units, calendar=self.calendar,
                                    only_use_cftime_datetimes=True)

    def date2num(self, d):
        return float(self.cftime.date2num(d, self.cf_units, calendar=self.calendar))

    def datetime(self, *args):
        return self.cftime.datetime(*args, calendar=self.calendar)

    def label(self, x):
        try:
            return str(self.num2date(x))
        except Exception:  # non-finite or out of range: no date to show
            return "?"


class _FileContext:
    """Per-file state shared by the rule groups, and the finding collector.

    ``tables_only`` (set by :func:`check_files`) tells the rule groups that the
    time axis and time series are not required in this file.
    """

    def __init__(self, path, ds, report, tables_only=False):
        self.path = str(path)
        self.ds = ds
        self.report = report
        self.tables_only = tables_only   # time and time-series rules skipped
        self.counts = {}
        self.ok = {}             # variable -> passed its structural check
        self.ids = None          # decoded source ids, when readable
        self.nsrc = None
        self.nx = self.ny = None
        self.ntime = None
        self.time_values = None  # float64 time values
        self.time_valid = False  # finite, no fill, strictly increasing
        self.tax = None          # _TimeAxis when units and calendar are valid
        self.calendar = None
        self.bounds = None       # float64 (n, 2) time_bnds when readable
        self.bounds_ok = None    # per-record: bounds finite and not fill
        self.targets = None      # dict of target arrays when readable
        self.bad_filters = {}    # variable -> disallowed HDF5 filters (P02)

    # -- collection -------------------------------------------------------

    def add(self, rule, text, variable=None, source_index=None, record=None,
            time=None, level=None):
        """Add a finding; ``level`` overrides the rule's level (S08 not checked: W)."""
        key = (rule, variable)
        n = self.counts.get(key, 0) + 1
        self.counts[key] = n
        if n > MAX_DETAILS:
            return
        sid = None
        if source_index is not None:
            source_index = int(source_index)
            if self.ids is not None and 0 <= source_index < len(self.ids):
                sid = self.ids[source_index]
        prefix = self.path + ": " + (variable + ": " if variable else "")
        self.report.findings.append(Finding(
            rule, level or S.RULES[rule], prefix + text, self.path, variable, sid,
            source_index, None if record is None else int(record),
            None if time is None else float(time)))

    def add_many(self, rule, variable, items, describe):
        """Add up to the remaining detail budget of ``items``; count the rest.

        ``describe(item)`` returns a dict of :meth:`add` keyword arguments
        including ``text``.
        """
        items = list(items) if not isinstance(items, np.ndarray) else items
        k = len(items)
        if k == 0:
            return
        key = (rule, variable)
        room = max(0, MAX_DETAILS - self.counts.get(key, 0))
        for item in items[:room]:
            kw = describe(item)
            self.add(rule, kw.pop("text"), variable, **kw)
        if k > room:
            self.counts[key] = self.counts.get(key, 0) + (k - room)

    def flush(self):
        for (rule, variable), n in self.counts.items():
            if n > MAX_DETAILS:
                prefix = self.path + ": " + (variable + ": " if variable else "")
                self.report.findings.append(Finding(
                    rule, S.RULES[rule],
                    prefix + "{0} more {1} findings not listed ({2} in total)".format(
                        n - MAX_DETAILS, rule, n),
                    self.path, variable))

    # -- descriptions -----------------------------------------------------

    def src(self, i):
        i = int(i)
        if self.ids is not None and 0 <= i < len(self.ids):
            return "source '{0}' (index {1})".format(self.ids[i], i)
        return "source index {0}".format(i)

    def tdesc(self, r):
        r = int(r)
        if self.time_values is None or not 0 <= r < len(self.time_values):
            return "record {0}".format(r)
        v = self.time_values[r]
        if self.tax is not None:
            return "time={0:.10g} {1} ({2}), record {3}".format(
                v, self.tax.units, self.tax.label(v), r)
        return "time={0:.10g}, record {1}".format(v, r)

    def tval(self, r):
        if self.time_values is None or not 0 <= int(r) < len(self.time_values):
            return None
        return self.time_values[int(r)]

    def cell(self, c):
        c = int(c)
        if self.nx:
            return "cell {0} (i={1}, j={2})".format(c, c % self.nx, c // self.nx)
        return "cell {0}".format(c)

    def target(self, k):
        t = self.targets
        s = int(t["source"][k])
        src = self.src(s) if self.nsrc is not None and 0 <= s < self.nsrc \
            else "source index {0}".format(s)
        return "target row {0} ({1}, {2})".format(int(k), src, self.cell(t["cell"][k]))


# ---------------------------------------------------------------------------
# Rule groups


def _check_structure(ctx):
    """S01-S07, S09, P02 and G01: file-level structure.

    Format (S01), schema version (S02), dimensions (S03), static-table
    variables (S04), undefined reserved names (S05), recommended global
    attributes (S06), no packing of model-read variables (S07), numeric
    missing-value attributes of their own type on model-read variables (S09),
    allowed HDF5 filters on model-read variables (P02) and the grid size
    attributes (G01). Also records ``ctx.bad_filters`` for the P02 read guard
    in :func:`_check_timeseries`. In tables-only mode the ``time`` dimension
    (S03) and the ``time`` and ``runoff_flux`` variables (S04) are not
    required; the form of those present is still checked (S04).
    """
    ds = ctx.ds
    # Tables-only mode: the time axis and time series are not required.
    optional = (S.DIM_TIME, "time") + S.REQUIRED_TIMESERIES if ctx.tables_only else ()
    # S01
    if ds.data_model not in S.NETCDF_FORMATS:
        ctx.add("S01", "file format is {0}; schema {1} requires {2}".format(
            ds.data_model, S.SCHEMA_VERSION, " or ".join(S.NETCDF_FORMATS)))
    attrs = ds.ncattrs()

    # S02
    key = "mitgcm_runoff_schema_version"
    if key not in attrs:
        ctx.add("S02", "required global attribute {0} is missing".format(key))
    else:
        v = _str_attr(ds.getncattr(key))
        parts = None if v is None else v.strip().split(".")
        if (parts is None or len(parts) != 2
                or not all(p.isdigit() for p in parts)):
            ctx.add("S02", "global attribute {0} = {1!r} is not a string "
                    "'MAJOR.MINOR'".format(key, ds.getncattr(key)))
        elif int(parts[0]) not in S.SUPPORTED_MAJOR_VERSIONS:
            ctx.add("S02", "schema version {0} is not supported (this checker reads "
                    "major version(s) {1})".format(v, S.SUPPORTED_MAJOR_VERSIONS))

    # S03
    for d in S.REQUIRED_DIMS:
        if d in optional:
            continue
        if d not in ds.dimensions:
            ctx.add("S03", "required dimension '{0}' is missing".format(d))
        elif len(ds.dimensions[d]) == 0:
            ctx.add("S03", "dimension '{0}' has size 0".format(d))
    if S.DIM_NV in ds.dimensions and len(ds.dimensions[S.DIM_NV]) != 2:
        ctx.add("S03", "dimension 'nv' has size {0}, not 2".format(
            len(ds.dimensions[S.DIM_NV])))
    if "source_id" in ds.variables:
        v = ds.variables["source_id"]
        if len(v.dimensions) == 2:
            dname = v.dimensions[1]
            size = len(ds.dimensions[dname])
            if size > S.MAX_ID_LEN:
                ctx.add("S03", "id string dimension '{0}' has size {1} > {2}".format(
                    dname, size, S.MAX_ID_LEN), variable="source_id")
    if S.DIM_SOURCE in ds.dimensions:
        ctx.nsrc = len(ds.dimensions[S.DIM_SOURCE])
    if S.DIM_TIME in ds.dimensions:
        ctx.ntime = len(ds.dimensions[S.DIM_TIME])

    # S04: static-table variables (time series: presence here, form in D01)
    for name, (dims, cls, req) in S.TABLE_VARIABLES.items():
        required = (req is True or (isinstance(req, str) and req in ds.dimensions)) \
            and name not in optional
        if name not in ds.variables:
            ctx.ok[name] = False
            if required:
                why = "" if req is True else " (dimension '{0}' exists)".format(req)
                ctx.add("S04", "required variable '{0}' is missing{1}".format(name, why))
            continue
        var = ds.variables[name]
        ok = _matches(var, dims, cls)
        ctx.ok[name] = ok
        if not ok:
            want = "({0})".format(", ".join("<strlen>" if d == "*" else d for d in dims))
            if cls == "string":
                want += " string, or char with a trailing string-length dimension"
            else:
                want += " " + cls
            ctx.add("S04", "has dimensions ({0}) and type {1}; expected {2}".format(
                ", ".join(var.dimensions), _describe_type(var), want), variable=name)
    for name in S.REQUIRED_TIMESERIES:
        if name not in ds.variables and name not in optional:
            ctx.add("S04", "required variable '{0}' is missing".format(name))

    # S05
    for name in ds.variables:
        if (name.startswith("runoff_") and name not in S.TIMESERIES_VARIABLES
                and not name.startswith(S.PTRACER_PREFIX)):
            ctx.add("S05", "variable '{0}' is not defined by the schema (reserved "
                    "prefix 'runoff_'; defined: {1}, {2}<NAME>)".format(
                        name, ", ".join(S.TIMESERIES_VARIABLES), S.PTRACER_PREFIX),
                    variable=name)
    for a in attrs:
        if a.startswith(S.RESERVED_ATTR_PREFIX) and a not in S.MITGCM_ATTRS:
            ctx.add("S05", "global attribute '{0}' is not defined by the schema "
                    "(reserved prefix 'mitgcm_')".format(a))

    # S06
    missing = [a for a in S.RECOMMENDED_GLOBAL_ATTRS if a not in attrs]
    if missing:
        ctx.add("S06", "recommended global attributes missing: " + ", ".join(missing))

    # S07: no packing of model-read variables
    for name, var in ds.variables.items():
        if not _model_read(name):
            continue
        packing = [a for a in S.PACKING_ATTRS if a in var.ncattrs()]
        if packing:
            ctx.add("S07", "carries {0}; the model reads stored values without "
                    "unpacking, so model-read variables must not be packed (values "
                    "are checked as stored)".format(" and ".join(
                        "{0} = {1}".format(a, np.asarray(var.getncattr(a)).tolist())
                        for a in packing)),
                    name)

    # S09: missing-value attributes of numeric model-read variables are numbers
    # of the variable's own type (the reader uses NF_GET_ATT_DOUBLE). source_id
    # is text, so its fill is text by definition and not checked here.
    for name, var in ds.variables.items():
        if not _model_read(name) or _type_class(var) not in ("int", "float"):
            continue
        vdt = np.dtype(var.dtype)
        for att in S.MISSING_VALUE_ATTRS:
            if att not in var.ncattrs():
                continue
            val = var.getncattr(att)
            if isinstance(val, (str, bytes)):
                kind = "text"
            else:
                adt = np.asarray(val).dtype
                if adt.kind in "fiu" and adt == vdt:
                    continue
                kind = "type {0}".format(adt)
            ctx.add("S09", "{0} = {1!r} is {2}; on a model-read variable it must be a "
                    "number of the variable's own type ({3}), because the model reads "
                    "it with NF_GET_ATT_DOUBLE".format(att, val, kind, vdt), name)

    # P02: only deflate, shuffle and fletcher32 on model-read variables
    ctx.bad_filters = {}
    for name, var in ds.variables.items():
        if not _model_read(name):
            continue
        try:
            filters = var.filters() or {}
        except Exception:  # noqa: BLE001 - no filter information (e.g. NETCDF3)
            continue
        bad = sorted(k for k, v in filters.items() if v and k not in S.ALLOWED_FILTER_KEYS)
        if bad:
            ctx.bad_filters[name] = bad
            ctx.add("P02", "uses HDF5 filter(s) {0}; model-read variables may use only "
                    "deflate (zlib), shuffle and fletcher32, because the model's netCDF "
                    "build may lack other filter plugins".format(", ".join(bad)), name)

    # G01
    dims = {}
    for key in ("mitgcm_grid_nx", "mitgcm_grid_ny"):
        if key not in attrs:
            ctx.add("G01", "required global attribute {0} is missing".format(key))
            continue
        raw = ds.getncattr(key)
        val, kind = _scalar(raw)
        if kind not in ("i", "u") or val is None or val <= 0:
            ctx.add("G01", "global attribute {0} = {1!r} is not a positive "
                    "integer".format(key, raw))
        else:
            dims[key] = int(val)
    if len(dims) == 2:
        ctx.nx, ctx.ny = dims["mitgcm_grid_nx"], dims["mitgcm_grid_ny"]


def _check_sources(ctx):
    """I01-I04."""
    ds = ctx.ds
    if ctx.ok.get("source_id"):
        var = ds.variables["source_id"]
        raw = np.asarray(var[:])
        ids, bad_utf = [], []
        for row in raw.reshape(raw.shape[0], -1):
            b = row.tobytes().rstrip(b"\x00").rstrip(b" ")
            try:
                ids.append(b.decode("ascii"))
                bad_utf.append(False)
            except UnicodeDecodeError:
                ids.append(b.decode("utf-8", errors="replace"))
                bad_utf.append(True)
        ctx.ids = ids
        for i, sid in enumerate(ids):
            if sid == "":
                ctx.add("I01", "source index {0}: source_id is empty".format(i),
                        "source_id", source_index=i)
            elif bad_utf[i]:
                ctx.add("I01", "{0}: source_id contains non-ASCII bytes".format(
                    ctx.src(i)), "source_id", source_index=i)
            elif len(sid) > S.MAX_ID_LEN:
                ctx.add("I01", "{0}: source_id has {1} characters (> {2})".format(
                    ctx.src(i), len(sid), S.MAX_ID_LEN), "source_id", source_index=i)
            elif not S.SOURCE_ID_RE.match(sid):
                ctx.add("I01", "{0}: source_id {1!r} must use only ASCII letters, "
                        "digits, '_', '-' and '.', and start with a letter or "
                        "digit".format(ctx.src(i), sid), "source_id", source_index=i)
        # I02 / I03
        seen, lower = {}, {}
        for i, sid in enumerate(ids):
            seen.setdefault(sid, []).append(i)
            lower.setdefault(sid.lower(), set()).add(sid)
        for sid, idx in seen.items():
            if len(idx) > 1 and sid != "":
                ctx.add("I02", "source_id '{0}' is used by source indices {1}".format(
                    sid, idx), "source_id", source_index=idx[0])
        for low, variants in lower.items():
            if len(variants) > 1:
                names = sorted(variants)
                ctx.add("I03", "source_ids {0} differ only in letter case".format(
                    names), "source_id", source_index=seen[names[0]][0])
    # I04
    if ctx.ok.get("source_type"):
        types, _missing = _read_strings(ds.variables["source_type"])
        for i, t in enumerate(types):
            if t not in S.SOURCE_TYPES:
                ctx.add("I04", "{0}: source_type {1!r} is not a recommended value "
                        "({2})".format(ctx.src(i), t, ", ".join(S.SOURCE_TYPES)),
                        "source_type", source_index=i)


def _check_aliases(ctx):
    """A01-A02."""
    ds = ctx.ds
    if S.DIM_ALIAS not in ds.dimensions:
        return
    src = names = None
    if ctx.ok.get("alias_source"):
        src = np.asarray(ds.variables["alias_source"][:]).astype(np.int64)
    if ctx.ok.get("alias_name"):
        names, missing = _read_strings(ds.variables["alias_name"])
        for k in np.nonzero(missing)[0]:
            s = None if src is None else src[k]
            where = "" if s is None or ctx.nsrc is None or not 0 <= s < ctx.nsrc \
                else " of " + ctx.src(s)
            ctx.add("A01", "alias row {0}{1}: alias_name is empty or missing".format(
                k, where), "alias_name",
                source_index=None if not where else s)
    if src is not None and ctx.nsrc is not None:
        bad = np.nonzero((src < 0) | (src >= ctx.nsrc))[0]
        ctx.add_many("A01", "alias_source", bad, lambda k: {
            "text": "alias row {0} ({1!r}): alias_source={2} is outside [0, {3})".format(
                k, None if names is None else names[k], src[k], ctx.nsrc)})
    if src is not None and names is not None:
        seen = {}
        for k, (s, n) in enumerate(zip(src.tolist(), names)):
            if n == "":
                continue
            seen.setdefault((s, n), []).append(k)
        for (s, n), rows in seen.items():
            if len(rows) > 1:
                ok = ctx.nsrc is not None and 0 <= s < ctx.nsrc
                ctx.add("A02", "{0} has alias {1!r} {2} times (alias rows {3})".format(
                    ctx.src(s) if ok else "source index {0}".format(s), n,
                    len(rows), rows), "alias_name", source_index=s if ok else None)


def _check_targets(ctx):
    """T01-T08."""
    ds = ctx.ds
    if not all(ctx.ok.get(v) for v in ("target_source", "target_cell", "target_fraction")):
        return
    ts = np.asarray(ds.variables["target_source"][:]).astype(np.int64)
    tc = np.asarray(ds.variables["target_cell"][:]).astype(np.int64)
    fvar = ds.variables["target_fraction"]
    fraw = np.asarray(fvar[:])
    ffill = _fill_mask(fraw, _fill_values(fvar))
    tf = fraw.astype(np.float64)
    if ctx.ok.get("target_level"):
        tl = np.asarray(ds.variables["target_level"][:]).astype(np.int64)
    else:
        tl = np.ones_like(ts)
    ctx.targets = {"source": ts, "cell": tc, "fraction": tf, "level": tl}
    nsrc = ctx.nsrc or 0

    # T01
    src_ok = (ts >= 0) & (ts < nsrc)
    ctx.add_many("T01", "target_source", np.nonzero(~src_ok)[0], lambda k: {
        "text": "target row {0} ({1}): target_source={2} is outside [0, {3})".format(
            k, ctx.cell(tc[k]), ts[k], nsrc)})
    # T02
    if ctx.nx is not None:
        ncell = ctx.nx * ctx.ny
        bad = np.nonzero((tc < 0) | (tc >= ncell))[0]
        ctx.add_many("T02", "target_cell", bad, lambda k: {
            "text": "{0}: target_cell={1} is outside [0, nx*ny={2})".format(
                ctx.target(k), tc[k], ncell),
            "source_index": ts[k] if src_ok[k] else None})
    # T03
    if ts.size:
        triples = np.stack([ts, tc, tl], axis=1)
        uniq, inverse, counts = np.unique(triples, axis=0, return_inverse=True,
                                          return_counts=True)
        inverse = np.ravel(inverse)
        dup_groups = np.nonzero(counts > 1)[0]

        def _dup(g):
            rows = np.nonzero(inverse == g)[0].tolist()
            s, c, lev = uniq[g]
            return {"text": "{0}, {1}, level {2} appears {3} times (target rows {4})"
                            .format(ctx.src(s) if 0 <= s < nsrc else
                                    "source index {0}".format(s),
                                    ctx.cell(c), lev, len(rows), rows),
                    "source_index": s if 0 <= s < nsrc else None}
        ctx.add_many("T03", "target_cell", dup_groups, _dup)
    # T04
    with np.errstate(invalid="ignore"):
        bad = ffill | ~np.isfinite(tf) | (tf < 0) | (tf > 1)
    ctx.add_many("T04", "target_fraction", np.nonzero(bad)[0], lambda k: {
        "text": "{0}: target_fraction={1} is not a finite value in [0, 1]".format(
            ctx.target(k), _value_label(tf[k], ffill[k])),
        "source_index": ts[k] if src_ok[k] else None})
    # T05: float64 sums per source
    if nsrc:
        sums = np.bincount(ts[src_ok], weights=tf[src_ok], minlength=nsrc)
        counts = np.bincount(ts[src_ok], minlength=nsrc)
        bad = []
        for s in range(nsrc):
            if counts[s] == 0 or not abs(sums[s] - 1.0) <= S.FRACTION_SUM_TOL:
                bad.append(s)

        def _sum(s):
            if counts[s] == 0:
                text = "{0} has no targets (its fractions must sum to 1)".format(
                    ctx.src(s))
            else:
                text = ("{0}: fractions of its {1} target(s) sum to {2:.12g} "
                        "(|sum - 1| = {3:.3g} > {4:g})".format(
                            ctx.src(s), counts[s], sums[s], abs(sums[s] - 1.0),
                            S.FRACTION_SUM_TOL))
            return {"text": text, "source_index": s}
        ctx.add_many("T05", "target_fraction", bad, _sum)
    # T06
    ctx.add_many("T06", "target_fraction", np.nonzero(tf == 0)[0], lambda k: {
        "text": "{0}: target_fraction is exactly 0".format(ctx.target(k)),
        "source_index": ts[k] if src_ok[k] else None})
    # T07
    if ctx.ok.get("target_level"):
        ctx.add_many("T07", "target_level", np.nonzero(tl != 1)[0], lambda k: {
            "text": "{0}: target_level={1}; schema {2} allows only 1".format(
                ctx.target(k), tl[k], S.SCHEMA_VERSION),
            "source_index": ts[k] if src_ok[k] else None})
    # T08
    if ts.size > 1:
        unsorted = (ts[1:] < ts[:-1]) | ((ts[1:] == ts[:-1]) & (tc[1:] < tc[:-1]))
        if unsorted.any():
            k = int(np.nonzero(unsorted)[0][0]) + 1
            ctx.add("T08", "targets are not sorted by source and then cell (first "
                    "out of order: target row {0})".format(k), "target_source")


def _check_time(ctx):
    """M01-M04 on ``time`` and ``time_bnds``, then the sampling and file-name checks.

    Units and calendar (M01, M02 warning), finite strictly increasing times
    (M03) and bounds (M04) are checked here; it then calls
    :func:`_check_sampling` (M05) and :func:`_check_yearly_name` (M06).
    """
    ds = ctx.ds
    if not ctx.ok.get("time"):
        return
    tv = ds.variables["time"]
    raw = np.asarray(tv[:])
    tfill = _fill_mask(raw, _fill_values(tv))
    t = raw.astype(np.float64)
    ctx.time_values = t
    n = t.size

    # M01 / M02
    units = _attr(tv, "units")
    su = _str_attr(units)
    match = None
    if units is None:
        ctx.add("M01", "has no units attribute", "time")
    else:
        match = None if su is None else S.TIME_UNITS_RE.match(su.strip())
        if match is None:
            ctx.add("M01", "units {0!r} is not '<days|hours|minutes|seconds> since "
                    "<YYYY-MM-DD>[ hh:mm[:ss]][Z]'".format(units), "time")
    cal_attr = _attr(tv, "calendar")
    if cal_attr is None:
        ctx.add("M02", "calendar attribute is missing; '{0}' is assumed".format(
            S.DEFAULT_CALENDAR), "time")
        cal = S.DEFAULT_CALENDAR
    else:
        cal = _str_attr(cal_attr)
        cal = None if cal is None else cal.strip().lower()   # CF: ignore case
        if cal not in S.CALENDARS:
            ctx.add("M01", "calendar {0!r} is not one of {1} (letter case "
                    "ignored)".format(cal_attr, ", ".join(S.CALENDARS)), "time")
            cal = None
    ctx.calendar = cal
    if match is not None and cal is not None:
        try:
            ctx.tax = _TimeAxis(su, match, cal)
        except (ValueError, TypeError) as e:
            ctx.add("M01", "units {0!r}: invalid reference date in calendar {1!r} "
                    "({2})".format(units, cal, e), "time")
    if "time_bnds" in ds.variables:
        bv = ds.variables["time_bnds"]
        bu = _attr(bv, "units")
        if bu is not None and (_str_attr(bu) or "").strip() != (su or "").strip():
            ctx.add("M01", "units {0!r} differ from time units {1!r}".format(bu, units),
                    "time_bnds")
        bc = _attr(bv, "calendar")
        tcal = _str_attr(cal_attr) if cal_attr is not None else S.DEFAULT_CALENDAR
        if bc is not None and (_str_attr(bc) or "").strip().lower() != \
                (tcal or "").strip().lower():
            ctx.add("M01", "calendar {0!r} differs from the time calendar {1!r}".format(
                bc, cal_attr), "time_bnds")

    # M03
    bad = tfill | ~np.isfinite(t)
    if n == 0:
        ctx.add("M03", "time has no records", "time")
    ctx.add_many("M03", "time", np.nonzero(bad)[0], lambda r: {
        "text": "record {0}: time is {1}".format(r, _value_label(t[r], tfill[r])),
        "record": r})
    if n > 1:
        with np.errstate(invalid="ignore"):
            nonincr = ~(np.diff(t) > 0) & ~bad[:-1] & ~bad[1:]
        ctx.add_many("M03", "time", np.nonzero(nonincr)[0], lambda r: {
            "text": "time is not strictly increasing: {0} is followed by {1}".format(
                ctx.tdesc(r), ctx.tdesc(r + 1)), "record": r + 1, "time": t[r + 1]})
    else:
        nonincr = np.zeros(0, dtype=bool)
    ctx.time_valid = n > 0 and not bad.any() and not nonincr.any()

    # M04
    if "time_bnds" in ds.variables:
        if ctx.ok.get("time_bnds") and len(ds.variables["time_bnds"].dimensions) == 2 \
                and ds.variables["time_bnds"].shape[1] == 2:
            bv = ds.variables["time_bnds"]
            braw = np.asarray(bv[:])
            bfill = _fill_mask(braw, _fill_values(bv))
            b = braw.astype(np.float64)
            rows_bad = (bfill | ~np.isfinite(b)).any(axis=1)
            ctx.bounds, ctx.bounds_ok = b, ~rows_bad
            ctx.add_many("M04", "time_bnds", np.nonzero(rows_bad)[0], lambda r: {
                "text": "{0}: bounds [{1}, {2}] are missing or non-finite".format(
                    ctx.tdesc(r), _value_label(b[r, 0], bfill[r, 0]),
                    _value_label(b[r, 1], bfill[r, 1])), "record": r, "time": ctx.tval(r)})
            both = ~rows_bad & ~bad
            with np.errstate(invalid="ignore"):
                outside = both & ~((b[:, 0] <= t) & (t <= b[:, 1]))
                nonpos = ~rows_bad & ~(b[:, 0] < b[:, 1])
            ctx.add_many("M04", "time_bnds", np.nonzero(outside)[0], lambda r: {
                "text": "{0}: bounds [{1:.10g}, {2:.10g}] do not contain the time".format(
                    ctx.tdesc(r), b[r, 0], b[r, 1]), "record": r, "time": t[r]})
            ctx.add_many("M04", "time_bnds", np.nonzero(nonpos)[0], lambda r: {
                "text": "{0}: bounds [{1:.10g}, {2:.10g}] do not have positive "
                        "length".format(ctx.tdesc(r), b[r, 0], b[r, 1]),
                "record": r, "time": ctx.tval(r)})
            if n > 1:
                tol = ctx.tax.tol if ctx.tax is not None else 0.0
                gap = (~rows_bad[:-1] & ~rows_bad[1:]
                       & (np.abs(b[:-1, 1] - b[1:, 0]) > tol))
                ctx.add_many("M04", "time_bnds", np.nonzero(gap)[0], lambda r: {
                    "text": "bounds are not contiguous: record {0} ends at {1:.10g} but "
                            "{2} starts at {3:.10g}".format(r, b[r, 1], ctx.tdesc(r + 1),
                                                           b[r + 1, 0]),
                    "record": r + 1, "time": ctx.tval(r + 1)})
    elif n > 1:
        ctx.add("M04", "time_bnds is required when there is more than one record "
                "({0} records)".format(n), "time_bnds")

    _check_sampling(ctx)
    _check_yearly_name(ctx)


def _check_sampling(ctx):
    """M05: mitgcm_time_* attributes, and their consistency with time (section 7)."""
    ds, t, n = ctx.ds, ctx.time_values, ctx.time_values.size
    attrs = ds.ncattrs()
    samp = None
    if "mitgcm_time_sampling" not in attrs:
        ctx.add("M05", "required global attribute mitgcm_time_sampling is missing")
    else:
        raw = ds.getncattr("mitgcm_time_sampling")
        samp = _str_attr(raw)
        if samp not in S.TIME_SAMPLINGS:
            ctx.add("M05", "mitgcm_time_sampling = {0!r} is not one of {1}".format(
                raw, ", ".join(S.TIME_SAMPLINGS)))
            samp = None
    period = None
    if "mitgcm_time_period" in attrs:
        raw = ds.getncattr("mitgcm_time_period")
        val, kind = _scalar(raw)
        if kind not in ("i", "u", "f") or val is None or not np.isfinite(val) or val <= 0:
            ctx.add("M05", "mitgcm_time_period = {0!r} is not a positive number of "
                    "seconds".format(raw))
        else:
            period = float(val)
    elif samp == "fixed":
        ctx.add("M05", "mitgcm_time_period is required when mitgcm_time_sampling "
                "is 'fixed'")
    repeat = S.DEFAULT_TIME_REPEAT
    if "mitgcm_time_repeat" in attrs:
        raw = ds.getncattr("mitgcm_time_repeat")
        repeat = _str_attr(raw)
        if repeat not in S.TIME_REPEATS:
            ctx.add("M05", "mitgcm_time_repeat = {0!r} is not one of {1}".format(
                raw, ", ".join(S.TIME_REPEATS)))
            repeat = None

    tax, b, bok = ctx.tax, ctx.bounds, ctx.bounds_ok
    if samp == "constant" and n != 1:
        ctx.add("M05", "mitgcm_time_sampling is 'constant' but time has {0} records "
                "(exactly 1 required)".format(n))
    if samp == "fixed" and period is not None and tax is not None and ctx.time_valid \
            and n > 1:
        dts = np.diff(t) * tax.factor
        with np.errstate(invalid="ignore"):
            bad = ~(np.abs(dts - period) <= S.TIME_EQUAL_TOL_SECONDS)
        ctx.add_many("M05", "time", np.nonzero(bad)[0], lambda r: {
            "text": "spacing from {0} to the next record is {1:.10g} s, not "
                    "mitgcm_time_period = {2:.10g} s".format(ctx.tdesc(r), dts[r], period),
            "record": r, "time": t[r]})
    if samp in ("monthly", "yearly") and tax is not None and b is not None:
        prev = None
        bad = []
        for r in range(n):
            if not bok[r]:
                prev = None
                continue
            # Which month/year the record belongs to: decode the start moved
            # forward by the tolerance, so a start up to the tolerance before
            # an edge (e.g. 86 us early) counts as that edge, not as the last
            # instant of the previous month.
            d0 = tax.num2date(b[r, 0] + tax.tol)
            if samp == "monthly":
                key = (d0.year, d0.month)
                nxt = (d0.year + d0.month // 12, d0.month % 12 + 1)
                exp0 = tax.datetime(d0.year, d0.month, 1)
                exp1 = tax.datetime(nxt[0], nxt[1], 1)
                what = "calendar month {0:04d}-{1:02d}".format(*key)
            else:
                key = d0.year
                exp0 = tax.datetime(d0.year, 1, 1)
                exp1 = tax.datetime(d0.year + 1, 1, 1)
                what = "calendar year {0:04d}".format(key)
            e0, e1 = tax.date2num(exp0), tax.date2num(exp1)
            if abs(b[r, 0] - e0) > tax.tol or abs(b[r, 1] - e1) > tax.tol:
                bad.append((r, "{0}: bounds {1} to {2} are not exactly the {3} "
                               "({4} to {5})".format(ctx.tdesc(r), tax.label(b[r, 0]),
                                                     tax.label(b[r, 1]), what, exp0, exp1)))
                prev = None
                continue
            if prev is not None:
                want = ((prev[0] + prev[1] // 12, prev[1] % 12 + 1)
                        if samp == "monthly" else prev + 1)
                if key != want:
                    bad.append((r, "{0}: {1} does not follow the previous record's "
                                   "period consecutively".format(ctx.tdesc(r), what)))
            prev = key
        ctx.add_many("M05", "time_bnds", bad, lambda item: {
            "text": item[1], "record": item[0], "time": ctx.tval(item[0])})
        # exf interpolates between period midpoints and ignores the file's
        # times, so time must be the midpoint of its bounds (section 3.1).
        mid = 0.5 * (b[:, 0] + b[:, 1])
        with np.errstate(invalid="ignore"):
            off = bok & np.isfinite(t) & ~(np.abs(t - mid) * tax.factor
                                           <= S.TIME_EQUAL_TOL_SECONDS)
        ctx.add_many("M05", "time", np.nonzero(off)[0], lambda r: {
            "text": "{0}: with mitgcm_time_sampling = '{1}', time must be the "
                    "midpoint of its bounds, {2:.10g} ({3}); it differs by "
                    "{4:.6g} s".format(ctx.tdesc(r), samp, mid[r], tax.label(mid[r]),
                                       (t[r] - mid[r]) * tax.factor),
            "record": r, "time": t[r]})
    if repeat == "annual":
        if b is None:
            if n == 1:
                ctx.add("M05", "mitgcm_time_repeat is 'annual' but there is no "
                        "time_bnds to show that the records cover exactly one year",
                        "time_bnds")
        elif tax is not None and bok.all() and n > 0:
            d0 = tax.num2date(b[0, 0])
            try:
                d1 = tax.datetime(d0.year + 1, d0.month, d0.day, d0.hour, d0.minute,
                                  d0.second, d0.microsecond)
                e1 = tax.date2num(d1)
            except ValueError:
                d1, e1 = None, None
            if e1 is None or abs(b[-1, 1] - e1) > tax.tol:
                ctx.add("M05", "mitgcm_time_repeat is 'annual' but the records cover "
                        "{0} to {1}, not exactly one year (expected end {2})".format(
                            tax.label(b[0, 0]), tax.label(b[-1, 1]), d1), "time_bnds")
            # exf monthly records (period -12) run January to December.
            first = tax.num2date(b[0, 0] + tax.tol)
            if samp == "monthly" and first.month != 1:
                ctx.add("M05", "a monthly annual climatology must start in January "
                        "(exf monthly records run January to December), but the first "
                        "record starts {0}".format(tax.label(b[0, 0])), "time_bnds",
                        record=0, time=ctx.tval(0))


def _check_yearly_name(ctx):
    """M06: a ``_YYYY.nc`` file holds only records of year YYYY."""
    m = S.YEARLY_FILE_RE.search(os.path.basename(ctx.path))
    if m is None or ctx.tax is None:
        return
    year = int(m.group(1))
    tax, t = ctx.tax, ctx.time_values
    try:
        lo = tax.date2num(tax.datetime(year, 1, 1))
        hi = tax.date2num(tax.datetime(year + 1, 1, 1))
    except ValueError:
        return
    if ctx.bounds is not None:
        b, ok = ctx.bounds, ctx.bounds_ok
        with np.errstate(invalid="ignore"):
            out = ok & ((b[:, 0] < lo - tax.tol) | (b[:, 1] > hi + tax.tol))
        desc = lambda r: "bounds {0} to {1}".format(tax.label(b[r, 0]), tax.label(b[r, 1]))
    else:
        with np.errstate(invalid="ignore"):
            out = np.isfinite(t) & ((t < lo) | (t >= hi))
        desc = lambda r: "time {0}".format(tax.label(t[r]))
    ctx.add_many("M06", "time", np.nonzero(out)[0], lambda r: {
        "text": "{0}: {1} lie outside year {2} named by the file".format(
            ctx.tdesc(r), desc(r), year), "record": r, "time": ctx.tval(r)})


def _check_timeseries(ctx, max_block_bytes):
    """D01, D07, D09, P01 and P02 on each time series, then its values.

    Form and units (D01), ptracer name and units (D07), an Inf
    ``_FillValue``/``missing_value`` on the temperature (D09), chunking (P01)
    and, for a variable with a disallowed filter, values that this netCDF
    library can't read (P02). The values are then streamed in bounded blocks
    by :func:`_stream_values` (D02-D09).
    """
    ds = ctx.ds
    names = [n for n in ds.variables
             if n in S.TIMESERIES_VARIABLES or n.startswith(S.PTRACER_PREFIX)]
    for name in names:
        var = ds.variables[name]
        is_ptr = name.startswith(S.PTRACER_PREFIX)
        form_ok = True
        if tuple(var.dimensions) != S.TIMESERIES_DIMS:
            ctx.add("D01", "has dimensions ({0}); expected (time, source)".format(
                ", ".join(var.dimensions)), name)
            form_ok = False
        if _type_class(var) != "float":
            ctx.add("D01", "has type {0}; expected float or double".format(
                _describe_type(var)), name)
            form_ok = False
        units = _attr(var, "units")
        su = _str_attr(units)
        if is_ptr:
            tracer = name[len(S.PTRACER_PREFIX):]
            if not S.PTRACER_NAME_RE.match(tracer):
                ctx.add("D07", "ptracer name {0!r} must be letters, digits and '_' "
                        "matching a PTRACERS_names entry".format(tracer), name)
            if su is None or su.strip() == "":
                ctx.add("D07", "units attribute is missing or empty; it must equal "
                        "the ptracer's own units", name)
        elif units is None:
            ctx.add("D01", "has no units attribute; allowed: {0}".format(
                ", ".join(S.UNITS[name])), name)
        elif su is None or su.strip() not in S.UNITS[name]:
            ctx.add("D01", "units {0!r} are not allowed; allowed: {1}".format(
                units, ", ".join(S.UNITS[name])), name)
        # D09: an infinite missing-value marker on the temperature (text: S09)
        if name == S.TEMPERATURE_VAR:
            for att in S.MISSING_VALUE_ATTRS:
                if att not in var.ncattrs():
                    continue
                vals = np.asarray(_numeric_values(var.getncattr(att)), dtype=np.float64)
                if np.isinf(vals).any():
                    ctx.add("D09", "{0} = {1} is infinite; Inf is never a missing-value "
                            "marker (use a finite fill value or NaN)".format(
                                att, vals.tolist()), name)
        # P01
        try:
            chunks = var.chunking()
        except Exception:  # NETCDF3 files have no chunking information
            chunks = "contiguous"
        if isinstance(chunks, (list, tuple)) and chunks and chunks[0] > 1:
            ctx.add("P01", "chunk shape {0} spans {1} records along time; reading "
                    "one record decompresses several (use 1)".format(
                        tuple(chunks), chunks[0]), name)
        if form_ok and name in ctx.bad_filters:
            # The checker's own netCDF may lack the plugin too: report, don't crash.
            try:
                _stream_values(ctx, name, var, max_block_bytes)
            except (RuntimeError, OSError) as e:
                ctx.add("P02", "values could not be read with this netCDF library "
                        "({0}), so D rules were not checked".format(e), name)
        elif form_ok:
            _stream_values(ctx, name, var, max_block_bytes)


def _stream_values(ctx, name, var, max_block_bytes):
    """D02-D09 on the values of one time series, read in blocks of whole records.

    Flux: missing or non-finite (D02), negative (D03). Temperature: ±Inf
    (D09), outside the warning range (D04). Salinity: missing or negative
    (D05), above the warning threshold (D06). Ptracer: missing (D07), negative
    (D08). Blocks hold at most ``max_block_bytes``.
    """
    nt, ns = var.shape
    if nt == 0 or ns == 0:
        return
    itemsize = np.dtype(var.dtype).itemsize
    rows = max(1, int(max_block_bytes // (ns * itemsize)))
    fills = _fill_values(var)
    lo, hi = S.TEMPERATURE_RANGE
    nblocks = 0
    for i0 in range(0, nt, rows):
        i1 = min(nt, i0 + rows)
        x = np.asarray(var[i0:i1, :])
        nblocks += 1
        isfill = _fill_mask(x, fills)
        with np.errstate(invalid="ignore"):
            nan = np.isnan(x)
            missing = isfill | nan | np.isinf(x)
            checks = []
            if name == S.FLUX_VAR:
                checks.append(("D02", missing, "missing or non-finite flux"))
                checks.append(("D03", ~missing & (x < 0), "negative flux"))
            elif name == S.TEMPERATURE_VAR:
                # Fill/NaN mean "surface temperature" and are skipped; Inf is
                # never a missing-value marker, so every ±Inf is an error (D09),
                # even when it equals an (itself invalid) Inf _FillValue.
                inf = np.isinf(x)
                present = ~(isfill | nan | inf)
                checks.append(("D09", inf, "infinite temperature"))
                checks.append(("D04", present & ((x < lo) | (x > hi)),
                               "temperature outside [{0:g}, {1:g}] degC".format(lo, hi)))
            elif name == S.SALINITY_VAR:
                checks.append(("D05", missing, "missing or non-finite salinity"))
                checks.append(("D05", ~missing & (x < 0), "negative salinity"))
                checks.append(("D06", ~missing & (x > S.SALINITY_MAX),
                               "salinity above {0:g}".format(S.SALINITY_MAX)))
            else:
                checks.append(("D07", missing, "missing or non-finite concentration"))
                checks.append(("D08", ~missing & (x < 0), "negative concentration"))
        for rule, mask, what in checks:
            if not mask.any():
                continue
            rr, cc = np.nonzero(mask)
            ctx.add_many(rule, name, np.stack([rr, cc], axis=1),
                         lambda rc, what=what, x=x, isfill=isfill, i0=i0: {
                "text": "{0}: {1} {2} at {3}".format(
                    ctx.src(rc[1]), what, _value_label(x[rc[0], rc[1]],
                                                       isfill[rc[0], rc[1]]),
                    ctx.tdesc(i0 + rc[0])),
                "source_index": rc[1], "record": i0 + rc[0],
                "time": ctx.tval(i0 + rc[0])})
    ctx.report.stats.setdefault("blocks_read", {})[
        "{0}::{1}".format(ctx.path, name)] = nblocks


def _check_units(ctx):
    """U01: units of the schema table variables in :data:`schema.TABLE_UNITS`
    (required when the variable is present) and of index variables (forbidden).
    Time units are checked by :func:`_check_time` and time-series units by
    :func:`_check_timeseries`; user variables are not checked."""
    for name, var in ctx.ds.variables.items():
        units = _attr(var, "units")
        su = None if units is None else _str_attr(units)
        su = None if su is None else su.strip()
        if name in S.INDEX_VARIABLES:
            if units is not None:
                ctx.add("U01", "index variable has units {0!r}; it must have no units "
                        "attribute".format(units), name)
            continue
        allowed = S.TABLE_UNITS.get(name)
        if allowed is None:
            continue
        if units is None:
            ctx.add("U01", "units attribute is missing; allowed: {0}".format(
                ", ".join(allowed)), name)
        elif su not in allowed:
            ctx.add("U01", "units {0!r} are not allowed; allowed: {1}".format(
                units, ", ".join(allowed)), name)


# ---------------------------------------------------------------------------
# Grid checks (R01-R03)


class _Grid:
    """MITgcm grid output read lazily from ``grid_dir`` with :func:`rdmds`."""

    def __init__(self, grid_dir):
        self.dir = str(grid_dir)
        self._cache = {}

    def get(self, name):
        if name not in self._cache:
            self._cache[name] = _read_grid_field(self.dir, name)
        return self._cache[name]


def _read_grid_field(grid_dir, name):
    """Global 2D ``(ny, nx)`` array of ``name`` (level 1 for hFacC).

    Accepts a global ``name.meta/.data`` pair or tiled
    ``name.XXX.YYY.meta/.data`` files, as :func:`MITgcmutils.mds.rdmds` does.
    """
    from ..mds import rdmds, readmeta
    base = os.path.join(grid_dir, name)
    metas = (glob.glob(base + ".[0-9][0-9][0-9].[0-9][0-9][0-9].meta")
             or glob.glob(base + ".meta"))
    if not metas:
        raise CheckIOError("--grid-dir {0}: no {1}.meta or {1}.XXX.YYY.meta "
                           "found".format(grid_dir, name))
    try:
        ndims = len(readmeta(metas[0])[0])
        # Level 1 of a 3D hFacC; a 2D file has no level dimension.
        a = rdmds(base, lev=[0]) if name == "hFacC" and ndims >= 3 else rdmds(base)
    except Exception as e:  # noqa: BLE001 - report any reader failure as I/O
        raise CheckIOError("cannot read {0}: {1}".format(base, e))
    return np.asarray(a, dtype=np.float64)


def _check_grid(ctx, grid):
    """R01-R03, and G01 when the grid's shape doesn't match (nx, ny).

    ``target_lon``/``target_lat`` are not model-read and may be packed; R03
    compares their unpacked values (``scale_factor``/``add_offset`` applied).
    """
    if ctx.nx is None or ctx.targets is None:
        return
    nx, ny = ctx.nx, ctx.ny
    t = ctx.targets
    cells = t["cell"]
    rows = np.nonzero((cells >= 0) & (cells < nx * ny))[0]
    src_ok = (t["source"] >= 0) & (t["source"] < (ctx.nsrc or 0))

    def field2d(name):
        a = grid.get(name)
        if a.shape != (ny, nx):
            ctx.add("G01", "grid field {0} in {1} has shape (ny, nx) = {2}, but the "
                    "file declares mitgcm_grid_nx={3}, mitgcm_grid_ny={4}".format(
                        name, grid.dir, a.shape, nx, ny))
            return None
        return a.ravel()   # C order of (ny, nx): index i + nx*j

    hf = field2d("hFacC")
    if hf is not None:
        land = rows[hf[cells[rows]] == 0]
        ctx.add_many("R01", "target_cell", land, lambda k: {
            "text": "{0} is on land (hFacC = 0 at level 1 in {1})".format(
                ctx.target(k), grid.dir),
            "source_index": t["source"][k] if src_ok[k] else None})
    if ctx.ok.get("target_cell_area"):
        rac = field2d("RAC")
        if rac is not None:
            area = np.asarray(ctx.ds.variables["target_cell_area"][:]).astype(np.float64)
            with np.errstate(invalid="ignore"):
                bad = rows[~(np.abs(area[rows] - rac[cells[rows]])
                             <= S.AREA_RTOL * np.abs(rac[cells[rows]]))]
            ctx.add_many("R02", "target_cell_area", bad, lambda k: {
                "text": "{0}: target_cell_area={1:.9g} m2 differs from RAC={2:.9g} m2 "
                        "by more than {3:g} relative".format(
                            ctx.target(k), area[k], rac[cells[k]], S.AREA_RTOL),
                "source_index": t["source"][k] if src_ok[k] else None})
    for vname, gname in (("target_lon", "XC"), ("target_lat", "YC")):
        if not ctx.ok.get(vname):
            continue
        g = field2d(gname)
        if g is None:
            continue
        # Not model-read, so it may be packed: compare unpacked degrees.
        v = _unpacked(ctx.ds.variables[vname])
        with np.errstate(invalid="ignore"):
            d = v[rows] - g[cells[rows]]
            if vname == "target_lon":
                d = (d + 180.0) % 360.0 - 180.0
            bad = rows[~(np.abs(d) <= S.LONLAT_TOL)]
        ctx.add_many("R03", vname, bad, lambda k, v=v, g=g, gname=gname: {
            "text": "{0}: {1:.9g} differs from {2}={3:.9g} by more than {4:g} "
                    "degrees".format(ctx.target(k), v[k], gname, g[cells[k]],
                                     S.LONLAT_TOL),
            "source_index": t["source"][k] if src_ok[k] else None})


# ---------------------------------------------------------------------------
# Multi-file consistency (X01)


def _normalized(var):
    """Values of a static variable in a form comparable across files."""
    tc = _type_class(var)
    if tc in ("vstring", "char") and (tc == "vstring" or var.ndim >= 2):
        return ("str", tuple(var.dimensions[:1]), _read_strings(var)[0])
    a = np.asarray(var[:])
    if tc == "float":
        return ("float", tuple(var.dimensions), a.astype(np.float64))
    if tc == "int":
        return ("int", tuple(var.dimensions), a.astype(np.int64))
    return (tc, tuple(var.dimensions), a)


def _snapshot(ctx):
    ds = ctx.ds
    static = {}
    for name, var in ds.variables.items():
        if S.DIM_TIME not in var.dimensions:
            static[name] = _normalized(var)
    grid = {a: ds.getncattr(a) for a in ds.ncattrs() if a.startswith("mitgcm_grid_")}
    extent = None
    if ctx.tax is not None and ctx.time_valid:
        t, tax = ctx.time_values, ctx.tax

        def fields(x):   # calendar-free date fields, rebuilt in another file's calendar
            d = tax.num2date(x)
            return (d.year, d.month, d.day, d.hour, d.minute, d.second, d.microsecond)
        first_t, last_t = fields(t[0]), fields(t[-1])
        start = end = None
        if ctx.bounds is not None and ctx.bounds_ok.all():
            start, end = fields(ctx.bounds[0, 0]), fields(ctx.bounds[-1, 1])
        extent = (first_t, last_t, start, end)
    # Yearly files: offset of the first *time value* from 1 January of the
    # year in the file name, in seconds, in the file's own calendar. exf
    # (useExfYearlyFields) uses one fldStartTime, the first record's time
    # offset, for every year (section 7).
    year = offset = None
    m = S.YEARLY_FILE_RE.search(os.path.basename(ctx.path))
    if m is not None:
        year = int(m.group(1))
        if extent is not None:
            try:
                jan1 = ctx.tax.date2num(ctx.tax.datetime(year, 1, 1))
                offset = (ctx.time_values[0] - jan1) * ctx.tax.factor
            except ValueError:
                offset = None
    # Fixed sampling: the period, for the spacing across file boundaries.
    attrs = ds.ncattrs()
    sampling = (_str_attr(ds.getncattr("mitgcm_time_sampling"))
                if "mitgcm_time_sampling" in attrs else None)
    period = None
    if "mitgcm_time_period" in attrs:
        val, kind = _scalar(ds.getncattr("mitgcm_time_period"))
        if kind in ("i", "u", "f") and val is not None and np.isfinite(val) and val > 0:
            period = float(val)
    return {"path": ctx.path, "vars": set(ds.variables), "static": static,
            "grid": grid, "calendar": ctx.calendar, "extent": extent,
            "year": year, "offset": offset, "sampling": sampling, "period": period}


def _same(a, b):
    ka, da, va = a
    kb, db, vb = b
    if ka != kb or da != db:
        return False, None
    if ka == "str":
        if len(va) != len(vb):
            return False, None
        for i, (x, y) in enumerate(zip(va, vb)):
            if x != y:
                return False, (i, x, y)
        return True, None
    va, vb = np.asarray(va), np.asarray(vb)
    if va.shape != vb.shape:
        return False, None
    if ka == "float":
        eq = (va == vb) | (np.isnan(va) & np.isnan(vb))
    else:
        eq = va == vb
    if np.all(eq):
        return True, None
    i = int(np.nonzero(~np.ravel(eq))[0][0])
    return False, (i, np.ravel(va)[i], np.ravel(vb)[i])


def _compare_snapshots(report, ref, cur):
    """X01 between the first file ``ref`` and a later file ``cur``."""
    p, rp = cur["path"], ref["path"]
    if cur["vars"] != ref["vars"]:
        only_cur = sorted(cur["vars"] - ref["vars"])
        only_ref = sorted(ref["vars"] - cur["vars"])
        report._add("X01", "{0}: variable set differs from {1} (only here: {2}; only "
                    "there: {3})".format(p, rp, only_cur, only_ref), p)
    keys = set(ref["grid"]) | set(cur["grid"])
    for k in sorted(keys):
        a, b = ref["grid"].get(k), cur["grid"].get(k)
        if (a is None) != (b is None) or (a is not None and not np.array_equal(
                np.asarray(a), np.asarray(b))):
            report._add("X01", "{0}: global attribute {1} = {2!r} differs from {3} "
                        "({4!r})".format(p, k, b, rp, a), p)
    for name in sorted(set(ref["static"]) & set(cur["static"])):
        same, where = _same(ref["static"][name], cur["static"][name])
        if not same:
            detail = ("" if where is None else
                      " (first difference at index {0}: {1!r} here, {2!r} there)".format(
                          where[0], where[2], where[1]))
            report._add("X01", "{0}: {1}: table differs from {2}{3}; yearly files "
                        "must have identical source, alias and target tables".format(
                            p, name, rp, detail), p, name)


def _check_order(report, snaps):
    """X01: consecutive files, in the order given, form one continuous series.

    Calendars must map to the same MITgcm calendar (``standard`` =
    ``gregorian`` = ``proleptic_gregorian``; ``noleap`` = ``365_day``). Each
    file's first ``time_bnds`` start must equal the previous file's last end;
    without bounds, continuity can't be shown, which is itself an X01 error.
    With ``fixed`` sampling in both files, the first time of a file minus the
    last time of the previous one must equal ``mitgcm_time_period``. Dates are
    compared in the later file's calendar, with the section-7 time tolerance.
    ``_YYYY`` files are then checked for a common offset of their first time
    value from 1 January.
    """
    import cftime
    for prev, cur in zip(snaps[:-1], snaps[1:]):
        if prev["calendar"] is None or cur["calendar"] is None:
            continue  # invalid calendar: reported per file (M01)
        if S.CALENDARS[prev["calendar"]] != S.CALENDARS[cur["calendar"]]:
            report._add("X01", "{0}: calendar {1!r} (MITgcm {2}) differs from {3} "
                        "({4!r}, MITgcm {5})".format(
                            cur["path"], cur["calendar"], S.CALENDARS[cur["calendar"]],
                            prev["path"], prev["calendar"],
                            S.CALENDARS[prev["calendar"]]),
                        cur["path"], "time")
            continue
        if prev["extent"] is None or cur["extent"] is None:
            continue  # time problems are reported per file (M01/M03)

        def dt(fields):
            return cftime.datetime(*fields, calendar=cur["calendar"])
        p_first, p_last, p_start, p_end = prev["extent"]
        c_first, c_last, c_start, c_end = cur["extent"]
        tol = S.TIME_EQUAL_TOL_SECONDS
        broken = False       # a gap or overlap was already reported for this pair
        if p_end is not None and c_start is not None:
            # Continuity: first bound = previous file's last bound.
            gap = (dt(c_start) - dt(p_end)).total_seconds()
            broken = abs(gap) > tol
            if broken:
                kind = ("a gap of {0:.10g} s ({1:.6g} days)".format(gap, gap / 86400.0)
                        if gap > 0 else
                        "an overlap of {0:.10g} s ({1:.6g} days); the files are not "
                        "in time order or overlap (list files in time order)".format(
                            -gap, -gap / 86400.0))
                report._add("X01", "{0}: time_bnds start at {1}, but {2} ends at {3}: "
                            "{4}. Each file's first bound must equal the previous "
                            "file's last bound".format(
                                cur["path"], dt(c_start), prev["path"], dt(p_end), kind),
                            cur["path"], "time_bnds")
        else:
            lacking = [s["path"] for s, e in ((prev, p_end), (cur, c_start)) if e is None]
            report._add("X01", "{0}: continuity with {1} can't be checked because {2} "
                        "has no valid time_bnds; every file of a multi-file set needs "
                        "time_bnds".format(cur["path"], prev["path"],
                                           " and ".join(lacking)),
                        cur["path"], "time_bnds")
            broken = (dt(p_last) - dt(c_first)).total_seconds() >= -tol
            if broken:
                report._add("X01", "{0}: files are not in time order or overlap: last "
                            "time of {1} is {2} but first time of {3} is {4} (list "
                            "files in time order)".format(
                                cur["path"], prev["path"], dt(p_last), cur["path"],
                                dt(c_first)), cur["path"], "time")
        # Fixed sampling: the spacing across the file boundary is the period too
        # (not repeated when a gap or overlap was already reported for the pair).
        if not broken and prev["sampling"] == cur["sampling"] == "fixed" \
                and prev["period"] is not None and cur["period"] is not None:
            if abs(prev["period"] - cur["period"]) > tol:
                report._add("X01", "{0}: mitgcm_time_period = {1:.10g} s differs from "
                            "{2} ({3:.10g} s)".format(cur["path"], cur["period"],
                                                     prev["path"], prev["period"]),
                            cur["path"], "time")
            else:
                spacing = (dt(c_first) - dt(p_last)).total_seconds()
                if abs(spacing - cur["period"]) > tol:
                    report._add("X01", "{0}: first time {1} is {2:.10g} s after the last "
                                "time of {3} ({4}), not mitgcm_time_period = {5:.10g} s; "
                                "with fixed sampling the spacing across files must equal "
                                "the period".format(cur["path"], dt(c_first), spacing,
                                                    prev["path"], dt(p_last),
                                                    cur["period"]),
                                cur["path"], "time")
    _check_yearly_offsets(report, snaps)


def _check_yearly_offsets(report, snaps):
    """X01: ``_YYYY`` files share one offset of their first ``time`` value from
    1 January of the year in their name (exf's single ``fldStartTime``)."""
    yearly = [s for s in snaps if s["year"] is not None and s["offset"] is not None]
    if len(yearly) < 2:
        return
    ref = yearly[0]
    for cur in yearly[1:]:
        if abs(cur["offset"] - ref["offset"]) > S.TIME_EQUAL_TOL_SECONDS:
            report._add("X01", "{0}: first time is {1:.10g} s ({2:.6g} days) after "
                        "1 January {3}, but the first time of {4} is {5:.10g} s ({6:.6g} "
                        "days) after 1 January {7}; yearly files must share one start "
                        "offset (exf uses one fldStartTime for every year)".format(
                            cur["path"], cur["offset"], cur["offset"] / 86400.0,
                            cur["year"], ref["path"], ref["offset"],
                            ref["offset"] / 86400.0, ref["year"]),
                        cur["path"], "time")


def _check_tables_only(ctx, multi):
    """S10: say which rules tables-only mode skipped in this file.

    The skipped rules are :data:`schema.TABLES_ONLY_SKIPPED`, plus
    :data:`schema.TABLES_ONLY_SKIPPED_MULTI` when several files are checked
    together (``multi``). Returns that list, which :func:`check_files` also
    stores in ``Report.stats["tables_only"]``.
    """
    skipped = list(S.TABLES_ONLY_SKIPPED) + (list(S.TABLES_ONLY_SKIPPED_MULTI)
                                             if multi else [])
    ctx.add("S10", "tables-only check: the time axis and time series are not "
            "required, and these rules were not checked: " + "; ".join(skipped))
    return skipped


# ---------------------------------------------------------------------------
# Entry points


def check_files(paths, grid_dir=None, *, max_block_bytes=DEFAULT_BLOCK_BYTES,
                tables_only=False):
    """Check one runoff file, or several files of one data set, against schema 1.0.

    Parameters
    ----------
    paths : str, os.PathLike or list of them
        Files to check. Several files are also compared with each other (X01)
        and must be listed in time order.
    grid_dir : str or os.PathLike, optional
        Directory with MITgcm grid output (``hFacC``, ``RAC``, ``XC``, ``YC``,
        global or tiled ``.meta/.data``); enables R01-R03.
    max_block_bytes : int
        Upper bound on the bytes of one block of time-series records read at once.
    tables_only : bool
        Check only the source, alias and target tables and the global
        attributes, for a file without a time axis or time series (for
        example the output of :func:`MITgcmutils.runoff.targets.write_targets`).
        The time and time-series rules are skipped and each file gets one
        ``S10`` (I) finding listing them (see the module notes).

    Returns
    -------
    Report

    Raises
    ------
    CheckIOError, OSError
        A file or grid field can't be opened or read (CLI exit status 2).
    """
    import netCDF4
    if isinstance(paths, (str, bytes, os.PathLike)):
        paths = [paths]
    paths = [os.fspath(p) for p in paths]
    report = Report(files=paths, grid_dir=None if grid_dir is None else str(grid_dir))
    grid = None if grid_dir is None else _Grid(grid_dir)
    multi = len(paths) > 1
    ref, snaps = None, []
    for path in paths:
        try:
            ds = netCDF4.Dataset(path, "r")
        except OSError as e:
            raise CheckIOError("cannot open {0} as NetCDF: {1}".format(path, e))
        try:
            ds.set_auto_maskandscale(False)   # stored values: no masking, no unpacking
            ds.set_auto_chartostring(False)
            ctx = _FileContext(path, ds, report, tables_only=tables_only)
            _check_structure(ctx)
            _check_attr_types(ctx)
            _check_sources(ctx)
            _check_aliases(ctx)
            _check_targets(ctx)
            if tables_only:
                report.stats["tables_only"] = _check_tables_only(ctx, multi)
            else:
                _check_time(ctx)
                _check_timeseries(ctx, max_block_bytes)
            _check_units(ctx)
            if grid is not None:
                _check_grid(ctx, grid)
            ctx.flush()
            if multi:
                snap = _snapshot(ctx)
                if ref is None:
                    ref = snap
                    snaps.append(dict(snap, static=None))
                else:
                    _compare_snapshots(report, ref, snap)
                    snap.pop("static")   # only the first file's tables stay in memory
                    snaps.append(snap)
        finally:
            ds.close()
    if multi and not tables_only:
        _check_order(report, snaps)
    return report


def main(argv=None):
    """Command-line entry point; returns the exit status."""
    parser = argparse.ArgumentParser(
        prog="python -m MITgcmutils.runoff.check",
        description="Check sparse-runoff NetCDF files against schema {0}. Exit 0: no "
                    "errors; 1: errors (or warnings with --strict); 2: usage or I/O "
                    "problem.".format(S.SCHEMA_VERSION))
    parser.add_argument("files", nargs="+", metavar="FILE",
                        help="runoff file(s); several yearly files in time order")
    parser.add_argument("--grid-dir", metavar="DIR",
                        help="MITgcm grid output (hFacC, RAC, XC, YC) for R01-R03")
    parser.add_argument("--strict", action="store_true", help="fail on warnings too")
    parser.add_argument("--json", metavar="OUT", help="write the report as JSON")
    parser.add_argument("--tables-only", action="store_true",
                        help="check only the source, alias and target tables (a file "
                             "without time series); time and time-series rules are "
                             "skipped and listed in an S10 finding")
    try:
        args = parser.parse_args(argv)
    except SystemExit as e:  # usage error (2) or --help (0)
        return int(e.code or 0)
    try:
        report = check_files(args.files, grid_dir=args.grid_dir,
                             tables_only=args.tables_only)
    except OSError as e:
        print("error: {0}".format(e), file=sys.stderr)
        return 2
    print(report.format_text())
    if args.json:
        try:
            with open(args.json, "w") as f:
                json.dump(report.to_dict(), f, indent=2, default=_json_default)
        except OSError as e:
            print("error: cannot write {0}: {1}".format(args.json, e), file=sys.stderr)
            return 2
    return report.exit_code(strict=args.strict)


def _json_default(o):
    if isinstance(o, np.generic):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


if __name__ == "__main__":
    sys.exit(main())
