"""Convert dense MITgcm runoff files to sparse-runoff NetCDF (schema 1.0), and back.

:func:`dense_to_sparse` reads a dense ``pkg/exf`` ``runoffFile`` (m/s, one
global 2D record per time, big-endian float32 or float64 as ``exf_iprec``
says), the grid (cell area ``rA`` and the surface wet mask, ``hFacC`` at
level 1) and the exf timing settings, and writes a schema-1.0 file
(``docs/runoff_schema.md``). :func:`sparse_to_dense` is the inverse. The
command line is ``python -m MITgcmutils.runoff.convert`` (:func:`main`).

Values
------
* The global cell index is ``cell = i + nx*j`` in the dense record, which has
  shape ``(ny, nx)`` with ``i`` fastest. Nothing else about the grid kind is
  used, so lat-lon, cubed-sphere (exch2) and LLC layouts are handled alike.
* ``runoff_flux`` of a source is ``sum_c dense(c, t) * rA(c)`` in m3/s,
  computed and stored in float64 (``flux_dtype``). :func:`sparse_to_dense`
  computes ``sum_s flux * fraction / rA``. For one-cell sources this round
  trip gives back float32 dense values exactly and float64 dense values to
  one unit in the last place: ``(dense * rA) / rA`` is not always ``dense`` in
  float64, and where it is not, no float64 flux gives ``dense`` back.
* Sources. By default every wet cell that is nonzero in at least one record
  is one source with fraction 1. With a grouping file (``sources``) the cells
  of a group form one source whose fractions are the cells' shares of the
  group's flux. This needs shares that are constant in time: a group whose
  shares change by more than ``share_tol`` between records is **refused**
  with an error that names it, unless ``split_varying`` is set, in which case
  that group is written as one source per cell (``<group>_<cell>``). The
  same holds for a group with no defined fractions in [0, 1]: one whose flux
  sums to zero over all records, or with a cell of negative share (both need
  negative runoff).
* ``runoff_temperature`` (optional dense ``runoftempFile``, degC) is the
  fraction-weighted mean over the source's cells, which equals the
  flux-weighted mean because the shares are constant; for a one-cell source
  it is the cell value. It is recorded for every record, also where the flux
  is zero. A source has one temperature, so cell temperatures that differ
  within a group are replaced by this mean and are not recovered by
  :func:`sparse_to_dense`.
* A nonzero value on a land cell (``hFacC`` = 0 at level 1) and a NaN or Inf
  value are errors that name the cell. Negative values give a warning.
* The wet mask is level 1 of ``hFacC``, the surface level in z coordinates.
  In pressure coordinates the surface is level ``Nr``, and level 1 is dry
  wherever the ocean is shallower than the deepest level, so the converter
  refuses runoff on such cells when it reads the mask from ``grid_dir``. From
  Python, pass the surface mask as the ``hfac`` array.

Timing: how exf settings map to the time axis
---------------------------------------------
File and line numbers refer to ``pkg/exf``, ``pkg/cal`` and ``eesupp/src``
of the MITgcm source this module ships with.

``period = 0`` (constant)
    ``EXF_INIT_FLD`` reads record 1 once (``exf_init_fld.F``, lines 98-126)
    and ``EXF_SET_FLD`` skips the field afterwards (``exf_set_fld.F``, line
    120). In a build with ``ALLOW_GENTIM2D_CONTROL`` the test at line 118 has
    no period condition, and ``EXF_GetFFieldRec`` then returns record 1 with
    weight 1 (``exf_getffieldrec.F``, lines 97-105), which gives the same
    field. Written as one record, ``mitgcm_time_sampling = "constant"``, no
    bounds. Only record 1 of the dense file is used.

``period > 0`` with a start date (``pkg/cal`` in use)
    The start date ``startdate1`` (YYYYMMDD) and ``startdate2`` (HHMMSS) is
    the time of record 1 (``exf_getffield_start.F``, lines 80-104), and
    record ``k`` is at ``start + (k-1)*period``: ``count0 =
    INT((t+0.5)/period) + 1`` with weight ``1 - MOD(t,period)/period``
    (``exf_getffieldrec.F``, lines 116-132 and 149). Written as
    ``"fixed"`` with ``mitgcm_time_period = period`` and ``time`` at the
    record times.

``period > 0`` with a start time in seconds (no ``pkg/cal``)
    Without the calendar the start time defaults to 0 and is given as
    ``runoffStartTime`` in seconds of model time (``exf_getffield_start.F``,
    lines 62-64 and 106), and the records are found by
    ``GET_PERIODIC_INTERVAL`` (``exf_getffieldrec.F``, lines 219-231;
    ``get_periodic_interval.F``, lines 96-133): record ``k`` is at
    ``start_time + (k-1)*period``. Written as ``"fixed"``, with model time 0
    at the reference date of the time units.

Repeat cycle (``repeat_cycle > 0``, the exf ``runoffRepCycle``, which
defaults to ``repeatPeriod``: ``exf_readparms.F``, line 951)
    Time since record 1 is taken modulo the cycle (``exf_getffieldrec.F``,
    lines 134-146; ``get_periodic_interval.F``, lines 83-94 and 115-133), so
    the file holds ``cycle / period`` records. Written with
    ``mitgcm_time_repeat = "annual"``; the cycle is the span of the bounds
    (schema section 7). Schema 1.0 can express only a cycle of exactly one
    year of the file's calendar (rule M05); any other cycle is refused.

``period = -12`` (monthly climatology)
    Twelve records, January to December (``exf_set_fld.F``, lines 133-140),
    placed at the middle of each calendar month and repeated every model year
    by ``cal_GetMonthsRec`` (``cal_getmonthsrec.F``, lines 106-117 and
    134-214). The start date and repeat cycle are not used. Written as
    ``"monthly"`` with ``"annual"`` repeat, on a nominal year
    (``clim_year``).

``period = -1`` (calendar months, not repeated)
    One record per consecutive calendar month, at mid-month
    (``exf_set_fld.F``, lines 142-153). Record 1 is the month of the start
    date: ``count = (year - yy)*12 + month - mm + 1``
    (``exf_getmonthsrec.F``, lines 58-69; start date converted in
    ``exf_getffield_start.F``, lines 80-81). With yearly files each file
    holds the twelve months of its year. Written as ``"monthly"``.

Yearly files (``useExfYearlyFields``)
    The file of year ``YYYY`` is ``<name>_YYYY``, except for ``period = -12``,
    which stays one file (``exf_getyearlyfieldname.F``, lines 45-57). Yearly
    files need ``pkg/cal`` and cannot be combined with a repeat period
    (``exf_check.F``, lines 71-84). Only the start date's offset from
    1 January is kept, and it applies to every year
    (``exf_getffield_start.F``, lines 85-92). Records are counted from that
    offset, and the record after the last one of a year is record 1 of the
    next year's file (``exf_getffieldrec.F``, lines 154-190). Written as
    ``<out>_YYYY.nc`` with ``"fixed"`` sampling, record ``k`` at ``1 January +
    offset + (k-1)*period``. All years share one source and target table.
    Schema 1.0 needs the spacing across the year boundary to equal the period
    (section 7, rule X01), so a year that is not a whole number of periods is
    refused.

Bounds
------
``time_bnds`` is not an exf input; exf interpolates between the record times
above. The converter writes, for ``"monthly"`` files, the calendar months
(``time`` is their midpoint, as ``cal_GetMonthsRec`` uses); for single
``"fixed"`` files, the interval of length ``period`` centred on each record
time; and for yearly ``"fixed"`` files, the intervals of length ``period``
counted from 1 January, because schema section 7 requires every bound of a
``_YYYY`` file to lie in its year and the files to join without a gap. In a
yearly file ``time`` is therefore ``offset`` after the start of its bounds
(at the start for an offset of 0).
"""

import argparse
import datetime
import os
import sys
import warnings

import numpy as np

from . import schema as S

__all__ = ["dense_to_sparse", "sparse_to_dense", "time_axis", "read_dense",
           "read_source_groups", "ConvertError", "main"]

#: MITgcm ``TheCalendar`` names accepted as aliases of CF calendar names.
_CAL_ALIASES = {"noleapyear": "noleap", "model": "360_day"}
#: Reference year of the time units when no date is given (model time 0).
DEFAULT_REF_YEAR = 1
#: Nominal year of a monthly climatology when no start date is given.
DEFAULT_CLIM_YEAR = 2001
#: Largest change of a cell's share of its source's flux between records
#: that still counts as constant.
DEFAULT_SHARE_TOL = 1e-9
_TOL = S.TIME_EQUAL_TOL_SECONDS
_DAY = 86400.0
_CHUNK_BYTES = 4 * 1024 * 1024


class ConvertError(ValueError):
    """Invalid input, or a dense file that schema 1.0 cannot express."""


# ---------------------------------------------------------------------------
# Dense files, grid and grouping


def read_dense(path, nx, ny, prec=32):
    """Records of a dense MITgcm 2D file as a ``(nrec, ny, nx)`` array.

    ``prec`` is 32 or 64 (``exf_iprec``); the file is big-endian, the MITgcm
    default byte order.
    """
    if prec not in (32, 64):
        raise ConvertError("prec must be 32 or 64, not {0!r}".format(prec))
    dtype = np.dtype(">f4" if prec == 32 else ">f8")
    a = np.fromfile(path, dtype=dtype)
    rec = nx * ny
    if a.size == 0 or a.size % rec:
        raise ConvertError(
            "{0}: {1} values of {2} bits is not a whole number of {3} x {4} "
            "records".format(path, a.size, prec, nx, ny))
    return a.reshape(a.size // rec, ny, nx)


def _grid(grid_dir, rac, hfac, xc, yc):
    """``(rac, hfac, xc, yc)`` as float64 ``(ny, nx)`` arrays; xc, yc may be None."""
    if grid_dir is not None:
        from .check import CheckIOError, _read_grid_field

        def get(name, required=True):
            try:
                return _read_grid_field(str(grid_dir), name)
            except CheckIOError as e:
                if required:
                    raise ConvertError(str(e))
                return None
        rac = get("RAC") if rac is None else rac
        hfac = get("hFacC") if hfac is None else hfac
        xc = get("XC", False) if xc is None else xc
        yc = get("YC", False) if yc is None else yc
    if rac is None or hfac is None:
        raise ConvertError("the grid is required: give grid_dir, or rac and hfac")
    rac = np.asarray(rac, dtype=np.float64)
    hfac = np.asarray(hfac, dtype=np.float64)
    if hfac.ndim == 3:
        hfac = hfac[0]
    if rac.ndim != 2 or hfac.shape != rac.shape:
        raise ConvertError("rac and hfac must be 2D (ny, nx) arrays of one shape; got "
                           "{0} and {1}".format(rac.shape, hfac.shape))
    out = [rac, hfac]
    for name, a in (("xc", xc), ("yc", yc)):
        if a is not None:
            a = np.asarray(a, dtype=np.float64)
            if a.shape != rac.shape:
                raise ConvertError("{0} has shape {1}, not {2}".format(
                    name, a.shape, rac.shape))
        out.append(a)
    return tuple(out)


def read_source_groups(path, nx, ny):
    """Read a source-to-cell grouping file; return ``[(source_id, cell), ...]``.

    One line per cell: ``source_id i j [cell ...]`` with 0-based ``i`` and
    ``j`` in the global 2D layout. A fourth column, when present, must equal
    ``i + nx*j``. Further columns, blank lines and ``#`` comments are
    ignored. This is the format of the lab_sea ``runoff_sources.txt``.
    """
    rows, seen = [], {}
    with open(path) as f:
        for n, line in enumerate(f, 1):
            text = line.split("#", 1)[0].split()
            if not text:
                continue
            where = "{0}, line {1}".format(path, n)
            if len(text) < 3:
                raise ConvertError(where + ": expected 'source_id i j [cell ...]'")
            sid = text[0]
            try:
                i, j = int(text[1]), int(text[2])
            except ValueError:
                raise ConvertError(where + ": i and j must be integers")
            if not (0 <= i < nx and 0 <= j < ny):
                raise ConvertError("{0}: (i, j) = ({1}, {2}) is outside the {3} x {4} "
                                   "grid".format(where, i, j, nx, ny))
            cell = i + nx * j
            if len(text) > 3:
                try:
                    given = int(text[3])
                except ValueError:
                    raise ConvertError(where + ": the cell column must be an integer")
                if given != cell:
                    raise ConvertError("{0}: cell {1} is not i + nx*j = {2}".format(
                        where, given, cell))
            if not S.SOURCE_ID_RE.match(sid):
                raise ConvertError("{0}: {1!r} is not a valid source_id".format(where, sid))
            if cell in seen:
                raise ConvertError("{0}: cell {1} is already in source {2!r}".format(
                    where, cell, seen[cell]))
            seen[cell] = sid
            rows.append((sid, cell))
    if not rows:
        raise ConvertError("{0}: no source cells".format(path))
    return rows


# ---------------------------------------------------------------------------
# Time axis


def _calendar(calendar):
    cal = str(calendar).strip().lower()
    cal = _CAL_ALIASES.get(cal, cal)
    if cal not in S.CALENDARS:
        raise ConvertError("calendar {0!r} is not one of {1}, noLeapYear, model".format(
            calendar, ", ".join(S.CALENDARS)))
    return cal


def _start_date(startdate1, startdate2, cal):
    """cftime date of exf's ``startdate1`` (YYYYMMDD) and ``startdate2`` (HHMMSS)."""
    import cftime
    d1, d2 = int(startdate1), int(startdate2)
    try:
        return cftime.datetime(d1 // 10000, d1 // 100 % 100, d1 % 100,
                               d2 // 10000, d2 // 100 % 100, d2 % 100, calendar=cal)
    except ValueError as e:
        raise ConvertError("startdate1 = {0}, startdate2 = {1} is not a date of the "
                           "{2} calendar ({3})".format(startdate1, startdate2, cal, e))


def _units(year):
    return "days since {0:04d}-01-01 00:00:00".format(year)


def _month_edges(year, month, count, ref, cal):
    """Days since ``ref`` of the starts of ``count + 1`` consecutive months."""
    import cftime
    edges = []
    for _ in range(count + 1):
        d = cftime.datetime(year, month, 1, calendar=cal)
        edges.append((d - ref).total_seconds() / _DAY)
        year, month = year + month // 12, month % 12 + 1
    return np.array(edges, dtype=np.float64)


def time_axis(period, nrec, *, startdate1=0, startdate2=0, start_time=None,
              repeat_cycle=0.0, calendar="standard", year=None, clim_year=None,
              ref_year=None):
    """Time axis of ``nrec`` dense records under the given exf timing.

    Returns a dict with ``time`` (days), ``bounds`` (``(n, 2)`` days, or None),
    ``units``, ``calendar``, ``sampling``, ``period`` (seconds, or None),
    ``repeat`` and ``nrec``, the number of dense records exf uses, which may
    be smaller than the number in the file. ``year`` selects the file of that
    year of a yearly set (``useExfYearlyFields``); ``ref_year`` is then the
    year of the time units, shared by all the files. See the module notes
    for the mapping and its source references.
    """
    import cftime
    cal = _calendar(calendar)
    period = float(period)
    cycle = float(repeat_cycle or 0.0)
    has_date = int(startdate1) != 0
    out = {"calendar": cal, "period": None, "repeat": "none"}

    if period == 0.0:
        if year is not None:
            raise ConvertError("a constant field (period 0) has no yearly files")
        ry = ref_year or (int(startdate1) // 10000 if has_date else DEFAULT_REF_YEAR)
        out.update(time=np.array([0.0]), bounds=None, units=_units(ry),
                   sampling="constant", nrec=1)
        return out

    if period in (-12.0, -1.0):
        if start_time is not None:
            raise ConvertError("period {0:g} needs pkg/cal dates, not a start time in "
                               "seconds (exf_set_fld.F stops without useCAL)".format(period))
        if period == -12.0:
            if year is not None:
                raise ConvertError("period -12 is one repeating file; yearly files of "
                                   "monthly records use period -1")
            y0 = clim_year or (int(startdate1) // 10000 if has_date else DEFAULT_CLIM_YEAR)
            m0, n = 1, 12
            out["repeat"] = "annual"
        elif year is not None:
            y0, m0, n = int(year), 1, 12
        else:
            if not has_date:
                raise ConvertError("period -1 needs startdate1: record 1 is the month "
                                   "of the start date")
            d = _start_date(startdate1, startdate2, cal)
            y0, m0, n = d.year, d.month, nrec
        if nrec < n:
            raise ConvertError("period {0:g} needs {1} records but the dense file has "
                               "{2}".format(period, n, nrec))
        ry = ref_year or y0
        ref = cftime.datetime(ry, 1, 1, calendar=cal)
        edges = _month_edges(y0, m0, n, ref, cal)
        out.update(time=0.5 * (edges[:-1] + edges[1:]),
                   bounds=np.stack([edges[:-1], edges[1:]], axis=1),
                   units=_units(ry), sampling="monthly", nrec=n)
        return out

    if period < 0.0:
        raise ConvertError("period {0:g} is not valid: exf accepts 0, a positive "
                           "number of seconds, -12 or -1".format(period))

    out.update(sampling="fixed", period=period)
    if year is not None:
        # exf_getffield_start.F 85-92 and exf_getffieldrec.F 154-190.
        if start_time is not None or not has_date:
            raise ConvertError("yearly files need pkg/cal: give startdate1, not a start "
                               "time in seconds")
        if cycle:
            raise ConvertError("yearly files cannot be combined with a repeat cycle "
                               "(exf_check.F, lines 71-84)")
        d = _start_date(startdate1, startdate2, cal)
        offset = (d - cftime.datetime(d.year, 1, 1, calendar=cal)).total_seconds()
        year = int(year)
        jan1 = cftime.datetime(year, 1, 1, calendar=cal)
        length = (cftime.datetime(year + 1, 1, 1, calendar=cal) - jan1).total_seconds()
        if offset >= period:
            raise ConvertError(
                "yearly files: the start date is {0:.10g} s after 1 January, which is "
                "not less than the period {1:.10g} s, so the first record of a year "
                "is not in its first interval".format(offset, period))
        n = int(np.ceil((length - offset) / period - 1e-9))
        if abs(n * period - length) > _TOL:
            raise ConvertError(
                "yearly files: year {0} has {1:.10g} s, which is not a whole number of "
                "periods of {2:.10g} s. Schema 1.0 (section 7, rule X01) needs the "
                "spacing across the year boundary to equal the period".format(
                    year, length, period))
        if nrec < n:
            raise ConvertError("year {0} needs {1} records but the dense file has "
                               "{2}".format(year, n, nrec))
        ry = ref_year or d.year
        base = (jan1 - cftime.datetime(ry, 1, 1, calendar=cal)).total_seconds()
        k = np.arange(n + 1, dtype=np.float64)
        edges = (base + k * period) / _DAY
        out.update(time=(base + offset + k[:-1] * period) / _DAY,
                   bounds=np.stack([edges[:-1], edges[1:]], axis=1),
                   units=_units(ry), nrec=n)
        return out

    if start_time is not None:
        if has_date:
            raise ConvertError("give either a start date (pkg/cal) or a start time in "
                               "seconds, not both")
        ry = ref_year or DEFAULT_REF_YEAR
        first = float(start_time)
    elif has_date:
        d = _start_date(startdate1, startdate2, cal)
        ry = ref_year or d.year
        first = (d - cftime.datetime(ry, 1, 1, calendar=cal)).total_seconds()
    else:
        raise ConvertError("a positive period needs startdate1 (with pkg/cal) or "
                           "start_time in seconds (without)")
    n = nrec
    if cycle:
        n = int(round(cycle / period))
        if n < 1 or abs(n * period - cycle) > _TOL:
            raise ConvertError("the repeat cycle {0:.10g} s is not a whole number of "
                               "periods of {1:.10g} s".format(cycle, period))
        if nrec < n:
            raise ConvertError("a repeat cycle of {0:.10g} s needs {1} records but the "
                               "dense file has {2}".format(cycle, n, nrec))
        out["repeat"] = "annual"
    k = np.arange(n + 1, dtype=np.float64)
    edges = (first - 0.5 * period + k * period) / _DAY
    if cycle:
        b0 = cftime.datetime(ry, 1, 1, calendar=cal) + datetime.timedelta(
            seconds=first - 0.5 * period)
        try:
            b1 = cftime.datetime(b0.year + 1, b0.month, b0.day, b0.hour, b0.minute,
                                 b0.second, b0.microsecond, calendar=cal)
            span = (b1 - b0).total_seconds()
        except ValueError:
            span = None
        if span is None or abs(span - cycle) > _TOL:
            raise ConvertError(
                "the repeat cycle {0:.10g} s is not one year of the {1} calendar. "
                "The year is counted from {2}, the start of the first record's "
                "interval, half a period before record 1. Schema 1.0 expresses a "
                "repeat only as mitgcm_time_repeat = 'annual', whose records cover "
                "exactly one year (section 7, rule M05)".format(cycle, cal, b0))
    out.update(time=(first + k[:-1] * period) / _DAY,
               bounds=np.stack([edges[:-1], edges[1:]], axis=1),
               units=_units(ry), nrec=n)
    return out


# ---------------------------------------------------------------------------
# Dense -> sparse


def _records(dense, nx, ny, prec, what):
    if isinstance(dense, (str, bytes, os.PathLike)):
        return read_dense(os.fspath(dense), nx, ny, prec)
    a = np.asarray(dense)
    if a.ndim == 2:
        a = a[None]
    if a.ndim != 3 or a.shape[1:] != (ny, nx):
        raise ConvertError("{0} has shape {1}, not (nrec, {2}, {3})".format(
            what, a.shape, ny, nx))
    return a


def _cell_label(cell, nx):
    return "cell {0} (i = {1}, j = {2}, 0-based)".format(int(cell), int(cell) % nx,
                                                         int(cell) // nx)


def _tables(sets, temps, rac, hfac, nx, groups, split_varying, share_tol):
    """Source and target tables, and the per-file flux and temperature series."""
    ncell = rac.size
    area, wet = rac.ravel(), hfac.ravel() > 0
    flat = [np.asarray(a, dtype=np.float64).reshape(a.shape[0], ncell) for a in sets]
    for a in flat:
        bad = ~np.isfinite(a)
        if bad.any():
            r, c = np.argwhere(bad)[0]
            raise ConvertError("dense runoff is {0} at {1}, record {2}: a missing or "
                               "non-finite flux is not allowed".format(
                                   a[r, c], _cell_label(c, nx), r + 1))
        land = (a != 0.0) & ~wet[None, :]
        if land.any():
            r, c = np.argwhere(land)[0]
            raise ConvertError(
                "dense runoff is {0:g} m/s on land at {1}, record {2}: hFacC is 0 at "
                "level 1 there ({3} land cells have nonzero runoff)".format(
                    a[r, c], _cell_label(c, nx), r + 1, int(land.any(axis=0).sum())))
        if (a < 0.0).any():
            r, c = np.argwhere(a < 0.0)[0]
            warnings.warn("dense runoff is negative in {0} values, first {1:g} m/s at "
                          "{2}, record {3}".format(int((a < 0.0).sum()), a[r, c],
                                                   _cell_label(c, nx), r + 1))
    active = np.zeros(ncell, dtype=bool)
    for a in flat:
        active |= (a != 0.0).any(axis=0)
    if (active & (area <= 0.0)).any():
        c = np.nonzero(active & (area <= 0.0))[0][0]
        raise ConvertError("rA is {0:g} at {1}, which has runoff".format(
            area[c], _cell_label(c, nx)))

    width = len(str(ncell - 1))
    members = []                       # (source_id, [cells], [fractions])
    grouped = np.zeros(ncell, dtype=bool)
    if groups:
        by_id = {}
        for sid, cell in groups:
            by_id.setdefault(sid, []).append(cell)
        for sid, cells in by_id.items():
            cells = np.array(sorted(cells))
            if (~wet[cells]).any():
                raise ConvertError("source {0!r}: {1} is on land (hFacC is 0 at level "
                                   "1)".format(sid, _cell_label(cells[~wet[cells]][0], nx)))
            grouped[cells] = True
            cells = cells[active[cells]]
            if cells.size == 0:
                warnings.warn("source {0!r} has no runoff in any record and is "
                              "left out".format(sid))
                continue
            if cells.size == 1:
                members.append((sid, list(cells), [1.0]))
                continue
            vol = np.concatenate([a[:, cells] * area[cells] for a in flat], axis=0)
            total = vol.sum(axis=1)
            net = total.sum()
            undefined = net == 0.0 or bool((vol.sum(axis=0) / net < 0.0).any())
            if not undefined:
                frac = vol.sum(axis=0) / net
                live = total != 0.0
                share = vol[live] / total[live, None]
                dev = np.abs(share - frac[None, :]).max() if live.any() else 0.0
            if undefined or dev > share_tol:
                if undefined and not split_varying:
                    raise ConvertError(
                        "source {0!r}: its cells have no fractions in [0, 1], because "
                        "{1}. Give each cell its own source, or set split_varying "
                        "(--split-varying)".format(
                            sid, "its flux sums to zero over all records" if net == 0.0
                            else "the runoff of some cell has the opposite sign of the "
                                 "source's total"))
                if not split_varying:
                    k = int(np.abs(share - frac[None, :]).max(axis=0).argmax())
                    raise ConvertError(
                        "source {0!r}: the cells' shares of the flux change in time "
                        "(the share of {1} ranges from {2:.9g} to {3:.9g}), and a "
                        "sparse source has one fixed fraction per cell. Give each cell "
                        "its own source, or set split_varying (--split-varying)".format(
                            sid, _cell_label(cells[k], nx), share[:, k].min(),
                            share[:, k].max()))
                for c in cells:
                    members.append(("{0}_{1:0{2}d}".format(sid, c, width), [c], [1.0]))
            else:
                members.append((sid, list(cells), list(frac / frac.sum())))
    for c in np.nonzero(active & ~grouped)[0]:
        members.append(("cell_{0:0{1}d}".format(c, width), [c], [1.0]))
    if not members:
        raise ConvertError("the dense runoff is zero everywhere: there is no source "
                           "to write")
    members.sort(key=lambda m: m[1][0])
    ids = [m[0] for m in members]
    if len(set(ids)) != len(ids):
        raise ConvertError("source ids are not unique after grouping")

    t_src = np.concatenate([[s] * len(m[1]) for s, m in enumerate(members)]).astype(np.int64)
    t_cell = np.concatenate([m[1] for m in members]).astype(np.int64)
    t_frac = np.concatenate([m[2] for m in members]).astype(np.float64)
    single = np.array([len(m[1]) == 1 for m in members])
    first = np.array([m[1][0] for m in members])
    fluxes, temperatures = [], []
    for n, a in enumerate(flat):
        flux = np.zeros((a.shape[0], len(members)))
        flux[:, single] = a[:, first[single]] * area[first[single]][None, :]
        temp = None
        if temps is not None:
            tt = np.asarray(temps[n], dtype=np.float64).reshape(temps[n].shape[0], ncell)
            if tt.shape[0] < a.shape[0]:
                raise ConvertError("the temperature file has {0} records, fewer than "
                                   "the {1} runoff records".format(tt.shape[0], a.shape[0]))
            tt = tt[:a.shape[0]]
            if not np.isfinite(tt[:, t_cell]).all():
                raise ConvertError("the dense temperature is not finite on a source cell")
            temp = np.zeros_like(flux)
            temp[:, single] = tt[:, first[single]]
        for s in np.nonzero(~single)[0]:
            cells, frac = np.array(members[s][1]), np.array(members[s][2])
            flux[:, s] = (a[:, cells] * area[cells]).sum(axis=1)
            if temp is not None:
                temp[:, s] = (tt[:, cells] * frac).sum(axis=1)
        fluxes.append(flux)
        temperatures.append(temp)
    return ids, t_src, t_cell, t_frac, fluxes, temperatures


def _chars(ids):
    n = max(len(s) for s in ids)
    return np.array([s.encode("ascii") for s in ids],
                    dtype="S{0}".format(n)).view("S1").reshape(len(ids), n), n


def _write(path, axis, ids, t_src, t_cell, t_frac, flux, temp, rac, xc, yc, nx, ny,
           flux_dtype, zlib, attrs):
    import netCDF4
    ns = len(ids)
    ds = netCDF4.Dataset(path, "w", format="NETCDF4")
    try:
        ds.set_ncstring_attrs(False)        # model-read text attributes are NC_CHAR
        g = {"Conventions": S.CONVENTIONS,
             "mitgcm_runoff_schema_version": S.SCHEMA_VERSION,
             "mitgcm_grid_nx": np.int32(nx), "mitgcm_grid_ny": np.int32(ny),
             "mitgcm_time_sampling": axis["sampling"]}
        if axis["period"] is not None:
            g["mitgcm_time_period"] = np.float64(axis["period"])
        g["mitgcm_time_repeat"] = axis["repeat"]
        g.update(attrs)
        for k, v in g.items():
            if v is not None:
                ds.setncattr(k, v)
        nt = axis["time"].size
        ds.createDimension("time", None)
        ds.createDimension("source", ns)
        ds.createDimension("target", t_cell.size)
        chars, strlen = _chars(ids)
        ds.createDimension("id_strlen", strlen)
        v = ds.createVariable("time", "f8", ("time",))
        v.setncatts({"long_name": "time", "standard_name": "time", "axis": "T",
                     "units": axis["units"], "calendar": axis["calendar"]})
        if axis["bounds"] is not None:
            v.setncattr("bounds", "time_bnds")
        v[:] = axis["time"]
        if axis["bounds"] is not None:
            ds.createDimension("nv", 2)
            v = ds.createVariable("time_bnds", "f8", ("time", "nv"))
            v.setncattr("long_name", "start and end of the interval each record covers")
            v[:] = axis["bounds"]
        v = ds.createVariable("source_id", "S1", ("source", "id_strlen"))
        v.set_auto_chartostring(False)
        v.setncatts({"long_name": "source identifier", "cf_role": "timeseries_id"})
        v[:] = chars

        def put(name, dtype, dim, values, **att):
            var = ds.createVariable(name, dtype, (dim,))
            var.setncatts(att)
            var[:] = values
        big = t_cell.max() > np.iinfo(np.int32).max
        put("target_source", "i4", "target", t_src,
            long_name="index of the source feeding this target",
            instance_dimension="source")
        put("target_cell", "i8" if big else "i4", "target", t_cell,
            long_name="0-based global cell index",
            comment="cell = i + mitgcm_grid_nx * j, 0-based (i, j) in the global 2D "
                    "layout of a dense runoffFile")
        put("target_fraction", "f8", "target", t_frac,
            long_name="share of the source flux sent to this cell", units="1")
        put("target_cell_area", "f8", "target", rac.ravel()[t_cell],
            long_name="horizontal cell area rA", units="m2")
        if xc is not None and yc is not None:
            put("target_lon", "f8", "target", xc.ravel()[t_cell],
                long_name="cell-center longitude", standard_name="longitude",
                units="degrees_east")
            put("target_lat", "f8", "target", yc.ravel()[t_cell],
                long_name="cell-center latitude", standard_name="latitude",
                units="degrees_north")
        item = np.dtype(flux_dtype).itemsize
        pieces = max(1, int(np.ceil(ns * item / _CHUNK_BYTES)))
        kw = {"chunksizes": (1, int(np.ceil(ns / pieces)))}
        if zlib:
            kw.update(zlib=True, complevel=2, shuffle=True)
        v = ds.createVariable(S.FLUX_VAR, flux_dtype, S.TIMESERIES_DIMS, **kw)
        v.setncatts({"long_name": "runoff volume flux of the source", "units": "m3 s-1",
                     "comment": "dense runoff (m/s) times cell area, summed over the "
                                "source's cells; split among targets by target_fraction"})
        v[:nt, :] = flux
        if temp is not None:
            v = ds.createVariable(S.TEMPERATURE_VAR, flux_dtype, S.TIMESERIES_DIMS, **kw)
            v.setncatts({"long_name": "runoff temperature", "units": "degC"})
            v[:nt, :] = temp
    finally:
        ds.close()


def dense_to_sparse(dense, out, *, grid_dir=None, rac=None, hfac=None, xc=None, yc=None,
                    prec=32, period=0.0, startdate1=0, startdate2=0, start_time=None,
                    repeat_cycle=0.0, years=None, calendar="standard", clim_year=None,
                    temperature=None, sources=None, split_varying=False,
                    share_tol=DEFAULT_SHARE_TOL, flux_dtype="f8", grid_name=None,
                    attrs=None, zlib=True):
    """Convert a dense exf runoff file to a sparse schema-1.0 NetCDF file.

    Parameters
    ----------
    dense : path or array
        The dense ``runoffFile`` (m/s), or its records as ``(nrec, ny, nx)``.
        With ``years`` it is the exf base name, and ``<dense>_YYYY`` is read
        for each year.
    out : path
        Output file. With ``years`` it is a base name, and ``<out>_YYYY.nc``
        is written for each year (a trailing ``.nc`` of ``out`` is dropped).
    grid_dir, rac, hfac, xc, yc
        The grid: a directory of MITgcm grid output (``RAC``, ``hFacC`` and
        optionally ``XC``, ``YC``, read with ``MITgcmutils.mds.rdmds``), or
        the arrays, each ``(ny, nx)`` (``hfac`` may be 3D; level 1 is used).
        ``(ny, nx)`` is the shape of the dense record.
    prec : 32 or 64
        ``exf_iprec`` of the dense files.
    period, startdate1, startdate2, start_time, repeat_cycle, years, calendar,
    clim_year
        exf timing, as in :func:`time_axis`. ``period`` is 0, a positive
        number of seconds, -12 or -1. ``years`` (a list) selects yearly
        files (``useExfYearlyFields``).
    temperature : path or array, optional
        Dense ``runoftempFile`` (degC), same layout and precision.
    sources : path or list, optional
        Grouping file (:func:`read_source_groups`) or ``[(source_id, cell)]``.
        Cells with runoff that are in no group are one-cell sources.
    split_varying : bool
        Write a group whose shares change in time, or that has no fractions
        in [0, 1], as one source per cell instead of raising
        :class:`ConvertError`.
    share_tol : float
        Largest change of a share between records that counts as constant.
    flux_dtype : ``"f8"`` or ``"f4"``
        Storage type of the time series. float64 reproduces the dense values
        to round-off; float32 does not (about 6e-8 relative).
    grid_name, attrs
        ``mitgcm_grid_name`` and further global attributes.

    Returns
    -------
    list of str
        The files written, in time order.

    Raises
    ------
    ConvertError
        Runoff on land, a non-finite value, a group whose shares vary in time
        or that has no fractions in [0, 1], or a timing that schema 1.0
        cannot express. Negative runoff only warns.
    """
    rac, hfac, xc, yc = _grid(grid_dir, rac, hfac, xc, yc)
    ny, nx = rac.shape
    if flux_dtype not in ("f8", "f4"):
        raise ConvertError("flux_dtype must be 'f8' or 'f4'")
    yearly = years is not None
    labels = [int(y) for y in years] if yearly else [None]
    if yearly and (not labels or labels != sorted(set(labels))):
        raise ConvertError("years must be a non-empty increasing list")
    if yearly and not isinstance(dense, (str, bytes, os.PathLike)):
        raise ConvertError("with years, dense is the base name of the yearly files")

    def name(base, y):
        # exf_getyearlyfieldname.F: <name>_YYYY
        return base if y is None else "{0}_{1:04d}".format(os.fspath(base), y)

    ref_year = None
    if yearly:
        ref_year = int(startdate1) // 10000 if int(startdate1) else labels[0]
    sets, temps, axes = [], ([] if temperature is not None else None), []
    for y in labels:
        a = _records(name(dense, y), nx, ny, prec, "dense")
        axis = time_axis(period, a.shape[0], startdate1=startdate1, startdate2=startdate2,
                         start_time=start_time, repeat_cycle=repeat_cycle,
                         calendar=calendar, year=y, clim_year=clim_year,
                         ref_year=ref_year)
        if a.shape[0] > axis["nrec"]:
            warnings.warn("{0}: exf uses {1} of the {2} records with this timing; the "
                          "others are not converted".format(
                              name(dense, y) if isinstance(dense, (str, os.PathLike))
                              else "dense", axis["nrec"], a.shape[0]))
        sets.append(a[:axis["nrec"]])
        axes.append(axis)
        if temps is not None:
            temps.append(_records(name(temperature, y), nx, ny, prec, "temperature"))
    groups = sources
    if isinstance(sources, (str, bytes, os.PathLike)):
        groups = read_source_groups(os.fspath(sources), nx, ny)
    ids, t_src, t_cell, t_frac, fluxes, temperatures = _tables(
        sets, temps, rac, hfac, nx, groups, split_varying, share_tol)

    base = os.fspath(out)
    if yearly and base.endswith(".nc"):
        base = base[:-3]
    paths = []
    for y, axis, flux, temp in zip(labels, axes, fluxes, temperatures):
        path = base if y is None else "{0}_{1:04d}.nc".format(base, y)
        g = {"title": "Sparse runoff converted from a dense MITgcm runoff file",
             "source": "MITgcmutils.runoff.convert.dense_to_sparse",
             "history": "converted from {0} by MITgcmutils.runoff.convert".format(
                 os.path.basename(name(dense, y)) if isinstance(
                     dense, (str, os.PathLike)) else "an array"),
             "references": "docs/runoff_schema.md (sparse runoff schema 1.0)",
             "exf_period": np.float64(period),
             "exf_startdate1": np.int32(startdate1), "exf_startdate2": np.int32(startdate2),
             "exf_repeat_cycle": np.float64(repeat_cycle or 0.0),
             "exf_use_yearly_fields": np.int32(yearly),
             "dense_precision": np.int32(prec)}
        if start_time is not None:
            g["exf_start_time"] = np.float64(start_time)
        if grid_name:
            g["mitgcm_grid_name"] = grid_name
        g.update(attrs or {})
        _write(path, axis, ids, t_src, t_cell, t_frac, flux, temp, rac, xc, yc, nx, ny,
               flux_dtype, zlib, g)
        paths.append(path)
    return paths


# ---------------------------------------------------------------------------
# Sparse -> dense


def sparse_to_dense(path, out=None, *, variable=S.FLUX_VAR, grid_dir=None, rac=None,
                    prec=None, fill=0.0):
    """Dense ``(nrec, ny, nx)`` records of a sparse runoff file.

    For ``runoff_flux`` the result is the runoff in m/s, ``sum_s flux_s *
    fraction / rA`` on each target cell and 0 elsewhere. ``rA`` is the file's
    ``target_cell_area``, or ``rac`` (an ``(ny, nx)`` array), or ``RAC`` of
    ``grid_dir``. For ``runoff_temperature`` (or another concentration) it is
    the mean over the sources feeding a cell, weighted by ``flux *
    fraction`` (by ``fraction`` in a record where that flux is zero), the
    source's own value where one source feeds the cell, and ``fill``
    elsewhere.

    ``prec`` (32 or 64) is the precision of the result; the default is the
    file's ``dense_precision`` attribute, or 64. With ``out`` the records are
    also written as a big-endian MITgcm binary file.
    """
    import netCDF4
    ds = netCDF4.Dataset(path, "r")
    try:
        ds.set_auto_maskandscale(False)
        nx, ny = int(ds.getncattr("mitgcm_grid_nx")), int(ds.getncattr("mitgcm_grid_ny"))
        src = np.asarray(ds.variables["target_source"][:]).astype(np.int64)
        cell = np.asarray(ds.variables["target_cell"][:]).astype(np.int64)
        frac = np.asarray(ds.variables["target_fraction"][:]).astype(np.float64)
        if variable not in ds.variables:
            raise ConvertError("{0} has no variable {1}".format(path, variable))
        flux = np.asarray(ds.variables[S.FLUX_VAR][:]).astype(np.float64)
        if prec is None:
            prec = int(ds.getncattr("dense_precision")) \
                if "dense_precision" in ds.ncattrs() else 64
        if grid_dir is not None or rac is not None:
            area = _grid(grid_dir, rac, np.ones((ny, nx)) if rac is not None else None,
                         None, None)[0]
            if area.shape != (ny, nx):
                raise ConvertError("the grid has shape {0}, not ({1}, {2})".format(
                    area.shape, ny, nx))
            area = area.ravel()[cell]
        elif "target_cell_area" in ds.variables:
            area = np.asarray(ds.variables["target_cell_area"][:]).astype(np.float64)
        else:
            raise ConvertError("{0} has no target_cell_area: give rac or "
                               "grid_dir".format(path))
        values = None
        if variable != S.FLUX_VAR:
            values = np.asarray(ds.variables[variable][:]).astype(np.float64)
    finally:
        ds.close()
    if prec not in (32, 64):
        raise ConvertError("prec must be 32 or 64, not {0!r}".format(prec))
    nt = flux.shape[0]
    dense = np.zeros((nt, nx * ny))
    count = np.bincount(cell, minlength=nx * ny)
    alone = count[cell] == 1
    if variable == S.FLUX_VAR:
        part = flux[:, src] * frac[None, :] / area[None, :]
        dense[:, cell[alone]] = part[:, alone]
        for r in range(nt):
            np.add.at(dense[r], cell[~alone], part[r, ~alone])
    else:
        dense[:] = fill
        dense[:, cell[alone]] = values[:, src[alone]]
        if (~alone).any():
            c, s, f = cell[~alone], src[~alone], frac[~alone]
            for r in range(nt):
                w = flux[r, s] * f
                den = np.bincount(c, weights=w, minlength=nx * ny)
                zero = den[c] == 0.0
                w = np.where(zero, f, w)
                den = np.bincount(c, weights=w, minlength=nx * ny)
                num = np.bincount(c, weights=w * values[r, s], minlength=nx * ny)
                dense[r, np.unique(c)] = (num / np.where(den == 0.0, 1.0, den))[np.unique(c)]
    dense = dense.reshape(nt, ny, nx).astype(">f4" if prec == 32 else ">f8")
    if out is not None:
        dense.tofile(out)
    return dense


# ---------------------------------------------------------------------------
# Command line


def main(argv=None):
    """Command-line entry point; returns the exit status.

    0: converted (and, unless ``--no-check``, the output has no checker
    errors); 1: invalid input or checker errors; 2: usage or I/O problem.
    """
    p = argparse.ArgumentParser(
        prog="python -m MITgcmutils.runoff.convert",
        description="Convert a dense MITgcm exf runoff file (m/s) to a sparse-runoff "
                    "NetCDF file (schema {0}), or back with --to-dense.".format(
                        S.SCHEMA_VERSION))
    p.add_argument("input", metavar="INPUT",
                   help="dense runoffFile (base name with --years), or with "
                        "--to-dense a sparse NetCDF file")
    p.add_argument("-o", "--output", required=True,
                   help="output file (base name with --years: OUT_YYYY.nc)")
    p.add_argument("--to-dense", action="store_true",
                   help="convert a sparse file back to a dense big-endian binary")
    p.add_argument("--grid-dir", help="directory with MITgcm grid output (RAC, hFacC, "
                                      "and XC, YC if present)")
    p.add_argument("--prec", type=int, choices=(32, 64), default=None,
                   help="exf_iprec of the dense files (default 32; with --to-dense, "
                        "the precision recorded in the sparse file)")
    p.add_argument("--period", type=float, default=0.0,
                   help="runoffperiod: 0 (constant), seconds, -12 (monthly "
                        "climatology) or -1 (calendar months)")
    p.add_argument("--startdate1", type=int, default=0, help="runoffstartdate1, YYYYMMDD")
    p.add_argument("--startdate2", type=int, default=0, help="runoffstartdate2, HHMMSS")
    p.add_argument("--start-time", type=float, default=None,
                   help="runoffStartTime in seconds of model time, for a run "
                        "without pkg/cal")
    p.add_argument("--repeat-cycle", type=float, default=0.0,
                   help="runoffRepCycle in seconds (exf default: repeatPeriod)")
    p.add_argument("--years", type=int, nargs="+", metavar="YYYY",
                   help="useExfYearlyFields: convert INPUT_YYYY for each year")
    p.add_argument("--calendar", default="standard",
                   help="CF calendar, or the MITgcm names gregorian, noLeapYear, model")
    p.add_argument("--clim-year", type=int, default=None,
                   help="nominal year of a monthly climatology (period -12)")
    p.add_argument("--temperature", help="dense runoftempFile (degC)")
    p.add_argument("--sources", help="grouping file: 'source_id i j [cell ...]' per cell")
    p.add_argument("--split-varying", action="store_true",
                   help="write a group whose shares change in time, or that has no "
                        "fractions in [0, 1], as one source per cell instead of "
                        "stopping")
    p.add_argument("--flux-dtype", choices=("f8", "f4"), default="f8")
    p.add_argument("--variable", default=S.FLUX_VAR,
                   help="with --to-dense: the time series to convert")
    p.add_argument("--grid-name", help="value of mitgcm_grid_name")
    p.add_argument("--no-check", action="store_true",
                   help="do not run the integrity checker on the output")
    args = p.parse_args(argv)
    try:
        if args.to_dense:
            a = sparse_to_dense(args.input, args.output, variable=args.variable,
                                grid_dir=args.grid_dir, prec=args.prec)
            print("{0}: {1} records of {2} x {3}".format(args.output, a.shape[0],
                                                         a.shape[2], a.shape[1]))
            return 0
        if args.grid_dir is None:
            p.error("--grid-dir is required")
        paths = dense_to_sparse(
            args.input, args.output, grid_dir=args.grid_dir, prec=args.prec or 32,
            period=args.period, startdate1=args.startdate1, startdate2=args.startdate2,
            start_time=args.start_time, repeat_cycle=args.repeat_cycle, years=args.years,
            calendar=args.calendar, clim_year=args.clim_year,
            temperature=args.temperature, sources=args.sources,
            split_varying=args.split_varying, flux_dtype=args.flux_dtype,
            grid_name=args.grid_name)
    except ConvertError as e:
        print("error: {0}".format(e), file=sys.stderr)
        return 1
    except OSError as e:
        print("error: {0}".format(e), file=sys.stderr)
        return 2
    for path in paths:
        print(path)
    if args.no_check:
        return 0
    from .check import check_files
    report = check_files(paths, grid_dir=args.grid_dir)
    print(report.format_text())
    return report.exit_code()


if __name__ == "__main__":
    sys.exit(main())
