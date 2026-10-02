"""Build the source, alias and target tables of a sparse-runoff file (schema 1.0).

:func:`build_targets` places each source (a river mouth, a glacier terminus,
...) on an MITgcm grid and decides which wet surface cells receive its water,
and with what fractions. :func:`write_targets` writes the resulting tables to
a new NetCDF file, or into an existing runoff file, replacing its tables. The
time series (``runoff_flux`` and so on) are added separately.

Library use::

    from MITgcmutils.runoff import build_targets, write_targets
    tables = build_targets("sources.csv", "run/", emission="spread",
                           spread_type="gaussian", spread_scale="25km")
    write_targets("targets.nc", tables)

Command line::

    python -m MITgcmutils.runoff.targets SOURCES.csv|SOURCES.nc --grid-dir DIR
        -o OUT.nc [--into RUNOFF.nc] [--emission {pointwise,spread}]
        [--spread-type {gaussian,exponential,linear}] [--spread-scale X]
        [--cutoff C] [--max-snap-distance D] [--earth-radius R]
        [--connectivity {auto,index,corners}] [--grid-name NAME] [--no-check]

Exit status: 0 success, 1 invalid input or checker errors in the output, 2
usage or I/O problem. Unless ``--no-check`` is given, the output is checked
with :func:`MITgcmutils.runoff.check_files` (``tables_only`` for a new file).

Inputs
------
* **Sources:** a CSV file with columns ``source_id``, ``lon``, ``lat``
  (degrees) and optional ``name``, ``type``, ``notes``, ``reference``,
  ``alt_names`` (semicolon-separated aliases) and per-source options
  ``emission``, ``spread_type``, ``spread_scale``, ``cutoff``,
  ``max_snap_distance``; an empty cell means "use the default". Or a
  schema-1.0 NetCDF file (``source_id``, ``source_lon``, ``source_lat``, the
  other ``source_*`` metadata and the alias table); its sources take every
  option from the defaults. Or a list of dicts with the CSV column names.
* **Grid:** MITgcm grid output read with :func:`MITgcmutils.mds.rdmds`
  (global or tiled files): ``hFacC`` (level 1 is the wet mask), ``XC``,
  ``YC``, ``XG``, ``YG`` (degrees) and ``RAC`` (m2), or a dict of these
  arrays. Every field is a global 2D ``(ny, nx)`` array and cell
  ``c = i + nx*j`` is its C-order ravel index, as in the checker.
* **Distances** are meters; strings such as ``"25km"`` or ``"25000 m"`` are
  accepted everywhere.

Algorithm (``docs/runoff_schema.md`` section 13)
-------------------------------------------------
1. **Snap:** each source goes to the nearest wet surface cell (``hFacC`` > 0
   at level 1, ``RAC`` > 0) by great-circle distance from its point to the
   cell center (haversine, radius ``earth_radius``). A nearest cell farther
   than ``max_snap_distance`` is an error naming the source, the distance and
   that cell. Exact ties go to the lowest cell index.
2. **Pointwise:** one target, fraction 1, at the snapped cell.
3. **Spread:** the distance ``r`` of a wet cell is the shortest path from the
   snapped cell through edge-connected wet cells, each step being the
   great-circle distance between the two cell centers (Dijkstra, stopped at
   the cutoff). ``r = 0`` at the snapped cell: the distance from the source
   point to that cell is limited by ``max_snap_distance`` and recorded in
   ``source_snap_distance``, but is not part of ``r``. Every cell reached
   with ``r <= cutoff`` is a candidate. Its weight is ``W(r) * RAC``, with
   ``X = spread_scale``:

   * exponential ``W = exp(-r/X)``, default cutoff ``3X``;
   * gaussian ``W = exp(-(r/X)**2)``, default cutoff ``3X``;
   * linear ``W = max(0, 1 - r/Rcut)``, ``Rcut = X/(1 - 1/e)``; the cutoff is
     ``Rcut``, or a smaller ``cutoff`` if one is given.

   ``W(X) = W(0)/e`` for all three. Zero weights are dropped and the rest
   normalized to sum to 1 (``math.fsum``; the largest share absorbs the
   rounding residual). The snapped cell has ``W(0) = 1``, the kernel peak,
   so every spread source has at least that target.
4. **Tables:** schema 1.0 sections 3.2-3.4, targets sorted by source and cell,
   with the provenance described in :func:`write_targets`.

Implementation notes
--------------------
* **Neighbour graph** (:func:`wet_graph`): two wet cells are neighbours when
  they share a cell edge, i.e. two corners. ``connectivity="auto"`` picks one
  of two methods from the grid geometry:

  - ``index``, for a grid stored as one logically rectangular block
    (lat-lon, and curvilinear grids without exch2): neighbours are
    ``(i+-1, j)`` and ``(i, j+-1)``, with a zonal wrap between ``i = nx-1``
    and ``i = 0`` when, in every row, ``SW(0, j)`` is as far from the center
    of ``(nx-1, j)`` as that cell's own SW corner (within
    :data:`WRAP_TOL`), i.e. it is that cell's SE corner.
  - ``corners``, for exch2 global I/O layouts (cubed sphere, LLC), whose
    array neighbours at face or tile boundaries are not grid neighbours.
    A cell's four corners are found as vertex ids: its SW corner
    ``(XG, YG)`` is its own vertex; the SW corners of its east, north and
    north-east array neighbours are accepted as its SE, NW and NE corners
    when the four points are a quadrilateral centered on the cell center
    (rotating the SW corner by 180 degrees about the center gives the NE
    corner, and the SE corner gives the NW, both within :data:`CORNER_TOL`
    of the shortest side). Otherwise, at face edges, the corners are found
    geometrically among all SW corners: the NE corner is the vertex nearest
    the 180-degree rotation of the SW corner; one of the other two is the
    nearest vertex to the center that is neither of those, and its
    180-degree rotation locates the last. A corner with no vertex within
    :data:`CORNER_TOL` of its prediction is not the SW corner of any cell
    (a cube corner, or the far side of an open boundary) and gets a
    synthesized position. Edges between two vertex ids are matched exactly;
    an edge from a vertex to a synthesized corner is matched with the other
    cell having the same vertex and a synthesized corner within
    :data:`SYNTH_MERGE_TOL` of the shortest side. An edge shared by more
    than two cells is an error. ``auto`` chooses ``corners`` when any array
    neighbour's SW corner is more than :data:`MOSAIC_RATIO` times as far from
    a cell's center as its own SW corner.

  On the ``global_ocean.cs32x15`` grid with every cell treated as wet, the
  ``corners`` method gives every cell 4 neighbours: 12288 edges, 384 of them
  across faces (12 cube edges x 32), 6144 stored and 2 synthesized vertices
  (Euler characteristic 2). The largest in-face rotation error there is 0.33
  of the shortest side; the synthesized cube corners are 0.73 of the
  shortest side from the nearest stored vertex, and the three predictions of
  one cube corner are 0.61 apart. Applied to all 5766 cells whose array
  neighbours lie on the same face, the geometric search
  (:func:`_geometric_corners`) finds the same four corners as the array. The
  ``corners`` method assumes cells that are close to parallelograms around
  their centers; no LLC grid was available to test it.
* **Search:** nearest-cell and nearest-vertex searches use
  ``scipy.spatial.cKDTree`` on unit vectors (chord distance, which orders
  points as great-circle distance does) when scipy is installed, and
  otherwise an exact pure-numpy search in bands of the z coordinate, which
  is slower for large grids. Both give the same results.
* **Dijkstra** is a bounded search with :mod:`heapq`, run separately for each
  spread source; its cost grows with the number of cells within the cutoff.
* SW corners whose unit vectors round to the same multiples of
  :data:`VERTEX_EPS` (about 6 mm on the Earth) are one vertex. Cells with
  ``RAC <= 0`` (e.g. exch2 blank tiles) are neither wet nor vertices.
* Only spherical grids in degrees are supported (Cartesian grids are not).
"""

import argparse
import csv
import datetime
import heapq
import math
import os
import re
import shlex
import sys
import tempfile
import warnings
from dataclasses import dataclass, field

import numpy as np

from . import schema as S

__all__ = ["build_targets", "write_targets", "read_sources", "read_grid", "wet_graph",
           "kernel_weight", "parse_distance", "TargetTables", "WetGraph", "BuildError",
           "main", "EARTH_RADIUS", "DEFAULT_MAX_SNAP_DISTANCE", "EMISSIONS",
           "SPREAD_TYPES", "CUTOFF_FACTOR", "LINEAR_CUTOFF_FACTOR"]

#: Default Earth radius for great-circle distances, m.
EARTH_RADIUS = 6371000.0
#: Default largest distance from a source to its snapped cell center, m.
DEFAULT_MAX_SNAP_DISTANCE = 50000.0
DEFAULT_EMISSION = "pointwise"
EMISSIONS = ("pointwise", "spread")
SPREAD_TYPES = ("gaussian", "exponential", "linear")
#: Default cutoff of the gaussian and exponential kernels, in units of X.
CUTOFF_FACTOR = 3.0
#: Linear kernel: ``Rcut = LINEAR_CUTOFF_FACTOR * X``, so that ``W(X) = 1/e``.
LINEAR_CUTOFF_FACTOR = 1.0 / (1.0 - math.exp(-1.0))
CONNECTIVITIES = ("auto", "index", "corners")

#: Neighbour-graph tolerances (module notes). Corner predictions must be
#: within CORNER_TOL of the cell's shortest side of a vertex.
CORNER_TOL = 0.5
#: Two synthesized predictions of one corner are at most this far apart, in
#: units of the shorter of the two cells' shortest sides.
SYNTH_MERGE_TOL = 1.0
#: ``auto``: an array neighbour's SW corner this many times farther from a
#: cell's center than the cell's own SW corner means an exch2 mosaic.
MOSAIC_RATIO = 2.5
#: ``index``: zonal wrap when every row's distance ratio is within 1 +- this.
WRAP_TOL = 0.1
#: SW corners whose unit vectors round to the same multiples of this (about
#: 6 mm on the Earth) are one vertex.
VERTEX_EPS = 1e-9
#: Use scipy's k-d tree when it is installed (tests set this False).
USE_SCIPY = True

BUILDER = "MITgcmutils.runoff.targets"
BUILDER_VERSION = "1.0"

#: CSV columns (module notes).
CSV_REQUIRED = ("source_id", "lon", "lat")
CSV_METADATA = {"name": "source_name", "type": "source_type", "notes": "source_notes",
                "reference": "source_reference"}
CSV_ALIASES = "alt_names"
OPTIONS = ("emission", "spread_type", "spread_scale", "cutoff", "max_snap_distance")

#: Errors listed in a BuildError message before summarizing.
MAX_LISTED = 20

_KERNEL_COMMENT = (
    "Targets built by {builder} {version}. Each source is snapped to the nearest wet "
    "surface cell (great-circle distance, radius {radius:.10g} m). A pointwise source puts "
    "fraction 1 in that cell. A spread source weights every wet cell reached within "
    "source_cutoff through edge-connected wet cells by W(r)*rA and normalizes the "
    "weights to sum to 1, where r is the shortest path between cell centers from the "
    "snapped cell, r = 0 there (target_distance; the source-to-cell distance is "
    "source_snap_distance and is not part of r), X is "
    "source_spread_scale, exponential W = exp(-r/X), gaussian W = exp(-(r/X)^2), "
    "linear W = max(0, 1 - r/Rcut) with Rcut = X/(1 - 1/e); W(X) = W(0)/e for all "
    "three. The default cutoff is 3X for exponential and gaussian and Rcut for linear.")


class BuildError(ValueError):
    """Invalid sources, options or grid; ``errors`` lists every problem found."""

    def __init__(self, errors):
        if isinstance(errors, str):
            errors = [errors]
        self.errors = list(errors)
        shown = self.errors[:MAX_LISTED]
        more = len(self.errors) - len(shown)
        text = "\n".join(shown)
        if more > 0:
            text += "\n... and {0} more error(s)".format(more)
        super().__init__(text)


# ---------------------------------------------------------------------------
# Small helpers

_DISTANCE_RE = re.compile(
    r"^\s*([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)\s*(m|km)?\s*$", re.IGNORECASE)


def parse_distance(value):
    """Distance in meters from a number, or a string like ``"25000"``,
    ``"25000 m"`` or ``"25km"`` (units ``m`` or ``km``, any letter case)."""
    if isinstance(value, (int, float, np.integer, np.floating)) and not isinstance(value, bool):
        return float(value)
    if isinstance(value, str):
        m = _DISTANCE_RE.match(value)
        if m:
            x = float(m.group(1))
            return x * 1000.0 if (m.group(2) or "").lower() == "km" else x
    raise ValueError("{0!r} is not a distance (a number of meters, or a number followed "
                     "by m or km)".format(value))


def _blank(value):
    """True for an unset option or metadata value (None, NaN or blank text)."""
    if value is None:
        return True
    if isinstance(value, str):
        return value.strip() == ""
    if isinstance(value, float):
        return math.isnan(value)
    return False


def _unit(lon, lat):
    """Unit vectors ``(..., 3)`` of points given in degrees."""
    lo = np.deg2rad(np.asarray(lon, dtype=np.float64))
    la = np.deg2rad(np.asarray(lat, dtype=np.float64))
    c = np.cos(la)
    return np.stack([c * np.cos(lo), c * np.sin(lo), np.sin(la)], axis=-1)


def _reflect(p, v):
    """Rotate unit vectors ``v`` by 180 degrees about the axes ``p``."""
    return 2.0 * np.sum(p * v, axis=-1)[..., None] * p - v


def _norm(v):
    return np.sqrt(np.sum(v * v, axis=-1))


def _haversine(lon1, lat1, lon2, lat2, radius):
    """Great-circle distance (m) between points in degrees (haversine formula)."""
    lon1, lat1, lon2, lat2 = (np.deg2rad(np.asarray(a, dtype=np.float64))
                              for a in (lon1, lat1, lon2, lat2))
    a = (np.sin(0.5 * (lat2 - lat1)) ** 2
         + np.cos(lat1) * np.cos(lat2) * np.sin(0.5 * (lon2 - lon1)) ** 2)
    return 2.0 * radius * np.arcsin(np.sqrt(np.clip(a, 0.0, 1.0)))


def kernel_weight(r, spread_type, spread_scale):
    """Kernel weight ``W(r)`` of a spread source (module notes); ``r`` and
    ``spread_scale`` in the same units. ``W(0) = 1`` and ``W(X) = 1/e``."""
    r = np.asarray(r, dtype=np.float64)
    x = float(spread_scale)
    if spread_type == "exponential":
        return np.exp(-r / x)
    if spread_type == "gaussian":
        return np.exp(-(r / x) ** 2)
    if spread_type == "linear":
        return np.maximum(0.0, 1.0 - r / (x * LINEAR_CUTOFF_FACTOR))
    raise ValueError("spread_type {0!r} is not one of {1}".format(
        spread_type, ", ".join(SPREAD_TYPES)))


# ---------------------------------------------------------------------------
# Nearest-point search


def _scipy_kdtree():
    if not USE_SCIPY:
        return None
    try:
        from scipy.spatial import cKDTree
    except ImportError:
        return None
    return cKDTree


class _PointIndex:
    """Nearest points among unit vectors, by chord distance.

    scipy's ``cKDTree`` when available; otherwise an exact search over the
    points whose z coordinate lies within the search radius of the query's,
    doubling the radius until ``k`` points lie within it.
    """

    def __init__(self, xyz):
        self.xyz = np.ascontiguousarray(xyz, dtype=np.float64)
        self.n = len(self.xyz)
        kdtree = _scipy_kdtree()
        self.tree = kdtree(self.xyz) if kdtree is not None and self.n else None
        if self.tree is None:
            self.order = np.argsort(self.xyz[:, 2], kind="stable")
            self.zs = self.xyz[self.order, 2]
            self.r0 = 4.0 * math.sqrt(4.0 * math.pi / max(self.n, 1))

    def query(self, q, k=1):
        """``(chord distances, indices)``, each ``(len(q), k)``, nearest first;
        equal distances are ordered by index. ``k`` is capped at the number of points."""
        q = np.atleast_2d(np.asarray(q, dtype=np.float64))
        k = min(k, self.n)
        if self.tree is not None:
            d, i = self.tree.query(q, k=k)
            d, i = d.reshape(len(q), k), i.reshape(len(q), k).astype(np.int64)
            order = np.lexsort((i, d), axis=1)   # ties by index, as the fallback
            return np.take_along_axis(d, order, 1), np.take_along_axis(i, order, 1)
        dist = np.empty((len(q), k))
        idx = np.empty((len(q), k), dtype=np.int64)
        for row, p in enumerate(q):
            r = self.r0
            while True:
                lo = np.searchsorted(self.zs, p[2] - r, side="left")
                hi = np.searchsorted(self.zs, p[2] + r, side="right")
                cand = self.order[lo:hi]
                d = _norm(self.xyz[cand] - p)
                if np.count_nonzero(d <= r) >= k or r >= 2.0:
                    sel = np.lexsort((cand, d))[:k]
                    dist[row], idx[row] = d[sel], cand[sel]
                    break
                r *= 2.0
        return dist, idx


# ---------------------------------------------------------------------------
# Sources


@dataclass
class _Sources:
    ids: list
    lon: np.ndarray
    lat: np.ndarray
    meta: dict            # schema variable name -> list of str, one per source
    aliases: list         # (source index, name, scheme or None)
    options: list         # per source: {option: raw value}
    where: list           # per source: where it was read, for messages


def read_sources(path):
    """Read a source table: a CSV file, or a schema-1.0 NetCDF file.

    Returns a list of dicts with the CSV column names (``source_id``,
    ``lon``, ``lat``, and those of ``name``, ``type``, ``notes``,
    ``reference``, ``alt_names``, and the option columns that are present);
    ``alt_names`` is a list. A NetCDF alias table also gives
    ``alt_name_schemes``. Each dict has ``_where`` naming its row or index.
    """
    path = os.fspath(path)
    with open(path, "rb") as f:
        magic = f.read(8)
    if magic.startswith(b"CDF") or magic.startswith(b"\x89HDF"):
        return _read_sources_nc(path)
    return _read_sources_csv(path)


def _read_sources_csv(path):
    known = set(CSV_REQUIRED) | set(CSV_METADATA) | {CSV_ALIASES} | set(OPTIONS)
    records, errors = [], []
    with open(path, newline="", encoding="utf-8-sig") as f:
        reader = csv.reader(f)
        try:
            header = next(reader)
        except StopIteration:
            raise BuildError("{0}: the file is empty (no header row)".format(path))
        names = [h.strip() for h in header]
        dup = sorted({n for n in names if names.count(n) > 1})
        if dup:
            raise BuildError("{0}: duplicate column(s) {1}".format(path, ", ".join(dup)))
        missing = [c for c in CSV_REQUIRED if c not in names]
        if missing:
            raise BuildError("{0}: required column(s) {1} missing (header: {2})".format(
                path, ", ".join(missing), ", ".join(names)))
        unknown = [n for n in names if n not in known]
        if unknown:
            warnings.warn("{0}: ignoring unknown column(s) {1}".format(
                path, ", ".join(unknown)), stacklevel=3)
        for row in reader:
            line = reader.line_num
            if not any(cell.strip() for cell in row):
                continue
            if len(row) > len(names):
                errors.append("{0}, line {1}: {2} fields, but the header has {3}".format(
                    path, line, len(row), len(names)))
                continue
            row = list(row) + [""] * (len(names) - len(row))
            rec = {n: v for n, v in zip(names, row) if n in known}
            if CSV_ALIASES in rec:
                rec[CSV_ALIASES] = rec[CSV_ALIASES].split(";")
            rec["_where"] = "{0}, line {1}".format(path, line)
            records.append(rec)
    if errors:
        raise BuildError(errors)
    return records


def _read_sources_nc(path):
    import netCDF4
    from .check import _read_strings, _unpacked
    records = []
    with netCDF4.Dataset(path, "r") as ds:
        ds.set_auto_maskandscale(False)
        ds.set_auto_chartostring(False)
        missing = [v for v in ("source_id", "source_lon", "source_lat") if v not in ds.variables]
        if missing:
            raise BuildError("{0}: variable(s) {1} missing; building targets needs "
                             "source_id, source_lon and source_lat".format(
                                 path, ", ".join(missing)))
        raw = np.asarray(ds.variables["source_id"][:])
        ids = [row.tobytes().rstrip(b"\x00").rstrip(b" ").decode("utf-8", errors="replace")
               for row in raw.reshape(raw.shape[0], -1)]
        lon = _unpacked(ds.variables["source_lon"])
        lat = _unpacked(ds.variables["source_lat"])
        meta = {}
        for col, var in CSV_METADATA.items():
            if var in ds.variables:
                meta[col] = _read_strings(ds.variables[var])[0]
        aliases = [[] for _ in ids]
        schemes = [[] for _ in ids]
        if "alias_source" in ds.variables and "alias_name" in ds.variables:
            asrc = np.asarray(ds.variables["alias_source"][:]).astype(np.int64)
            names = _read_strings(ds.variables["alias_name"])[0]
            sch = (_read_strings(ds.variables["alias_scheme"])[0]
                   if "alias_scheme" in ds.variables else [None] * len(names))
            for s, n, c in zip(asrc.tolist(), names, sch):
                if 0 <= s < len(ids):
                    aliases[s].append(n)
                    schemes[s].append(c if c else None)
    for k, sid in enumerate(ids):
        rec = {"source_id": sid, "lon": float(lon[k]), "lat": float(lat[k]),
               "alt_names": aliases[k], "alt_name_schemes": schemes[k],
               "_where": "{0}, source index {1}".format(path, k)}
        for col, values in meta.items():
            rec[col] = values[k]
        records.append(rec)
    return records


def _load_sources(sources):
    """Validate a source table (path or records) into :class:`_Sources`."""
    if isinstance(sources, (str, bytes, os.PathLike)):
        records = read_sources(sources)
    else:
        records = [dict(r) for r in sources]
        for k, r in enumerate(records):
            r.setdefault("_where", "source record {0}".format(k))
    if not records:
        raise BuildError("the source table has no sources")
    errors = []
    ids, lons, lats, where, options = [], [], [], [], []
    present = [col for col in CSV_METADATA if any(col in r for r in records)]
    meta = {CSV_METADATA[col]: [] for col in present}
    aliases = []
    seen = {}
    for k, r in enumerate(records):
        w = r["_where"]
        sid = r.get("source_id")
        sid = "" if sid is None else str(sid).strip()
        label = "source {0!r} ({1})".format(sid, w)
        if sid == "":
            errors.append("{0}: source_id is empty".format(w))
        elif len(sid) > S.MAX_ID_LEN:
            errors.append("{0}: source_id has {1} characters (> {2})".format(
                label, len(sid), S.MAX_ID_LEN))
        elif not S.SOURCE_ID_RE.match(sid):
            errors.append("{0}: source_id must use only ASCII letters, digits, '_', '-' "
                          "and '.', and start with a letter or digit".format(label))
        elif sid in seen:
            errors.append("{0}: source_id is also used by {1}".format(label, seen[sid]))
        else:
            seen[sid] = w
        coords = []
        for key, lim in (("lon", None), ("lat", 90.0)):
            v = r.get(key)
            try:
                x = float(v.strip() if isinstance(v, str) else v)
            except (TypeError, ValueError):
                x = float("nan")
            if not math.isfinite(x) or (lim is not None and abs(x) > lim):
                errors.append("{0}: {1} = {2!r} is not a finite number{3}".format(
                    label, key, v, "" if lim is None else " in [-90, 90]"))
            coords.append(x)
        ids.append(sid)
        lons.append(coords[0])
        lats.append(coords[1])
        where.append(w)
        for col in present:
            v = r.get(col)
            meta[CSV_METADATA[col]].append("" if _blank(v) else str(v).strip())
        names = r.get(CSV_ALIASES) or []
        if isinstance(names, str):
            names = names.split(";")
        schemes = r.get("alt_name_schemes") or [None] * len(names)
        done = set()
        for n, sc in zip(names, schemes):
            n = str(n).strip()
            if n and n not in done:
                done.add(n)
                aliases.append((k, n, sc))
        options.append({o: r.get(o) for o in OPTIONS})
    if errors:
        raise BuildError(errors)
    return _Sources(ids, np.array(lons), np.array(lats), meta, aliases, options, where)


def _resolve_options(src, defaults):
    """Per-source settings: the source's own value where set, else the default.

    Returns a list of dicts with ``emission``, ``spread_type``,
    ``spread_scale``, ``cutoff`` (the effective one) and ``max_snap``.
    """
    errors = []

    def value(name, raw, label):
        if _blank(raw):
            return None
        if name in ("emission", "spread_type"):
            v = str(raw).strip().lower()
            allowed = EMISSIONS if name == "emission" else SPREAD_TYPES
            if v not in allowed:
                errors.append("{0}: {1} {2!r} is not one of {3}".format(
                    label, name, raw, ", ".join(allowed)))
                return None
            return v
        try:
            x = parse_distance(raw)
        except ValueError as e:
            errors.append("{0}: {1}: {2}".format(label, name, e))
            return None
        ok = math.isfinite(x) and (x >= 0 if name == "max_snap_distance" else x > 0)
        if not ok:
            errors.append("{0}: {1} = {2!r} must be a finite distance {3} 0".format(
                label, name, raw, ">=" if name == "max_snap_distance" else ">"))
            return None
        return x

    base = {name: value(name, defaults.get(name), "default") for name in OPTIONS}
    if errors:
        raise BuildError(errors)
    out = []
    for k, opts in enumerate(src.options):
        label = "source {0!r} ({1})".format(src.ids[k], src.where[k])
        own = {name: value(name, opts.get(name), label) for name in OPTIONS}
        eff = {name: own[name] if own[name] is not None else base[name] for name in OPTIONS}
        em = eff["emission"] or DEFAULT_EMISSION
        res = {"emission": em, "spread_type": None, "spread_scale": None, "cutoff": None,
               "max_snap": (eff["max_snap_distance"] if eff["max_snap_distance"] is not None
                            else DEFAULT_MAX_SNAP_DISTANCE)}
        if em == "spread":
            st, x, cut = eff["spread_type"], eff["spread_scale"], eff["cutoff"]
            if st is None:
                errors.append("{0}: spread emission needs a spread_type ({1})".format(
                    label, ", ".join(SPREAD_TYPES)))
            if x is None:
                errors.append("{0}: spread emission needs a spread_scale".format(label))
            if st is not None and x is not None:
                if st == "linear":
                    rcut = x * LINEAR_CUTOFF_FACTOR
                    cut = rcut if cut is None else min(cut, rcut)
                elif cut is None:
                    cut = CUTOFF_FACTOR * x
                res.update(spread_type=st, spread_scale=x, cutoff=cut)
        out.append(res)
    if errors:
        raise BuildError(errors)
    return out


# ---------------------------------------------------------------------------
# Grid

_GRID_FIELDS = ("hFacC", "XC", "YC", "XG", "YG", "RAC")


def read_grid(grid_dir):
    """Global 2D ``(ny, nx)`` arrays of ``hFacC`` (level 1), ``XC``, ``YC``,
    ``XG``, ``YG`` and ``RAC`` from MITgcm grid output (global or tiled
    ``.meta/.data``). Raises :class:`OSError` when a field can't be read."""
    from .check import _read_grid_field
    return {name: _read_grid_field(os.fspath(grid_dir), name) for name in _GRID_FIELDS}


class _Geom:
    """Flattened grid geometry: cell ``c = i + nx*j``."""

    def __init__(self, grid, radius):
        if isinstance(grid, (str, bytes, os.PathLike)):
            grid = read_grid(grid)
        arrays = {}
        for name in _GRID_FIELDS:
            if name not in grid:
                raise BuildError("grid: field {0} is missing (needed: {1})".format(
                    name, ", ".join(_GRID_FIELDS)))
            a = np.asarray(grid[name], dtype=np.float64)
            if name == "hFacC" and a.ndim == 3:
                a = a[0]
            if a.ndim != 2:
                raise BuildError("grid: {0} has shape {1}; expected a 2D (ny, nx) "
                                 "array".format(name, a.shape))
            arrays[name] = a
        shapes = {name: a.shape for name, a in arrays.items()}
        if len(set(shapes.values())) != 1:
            raise BuildError("grid: fields have different shapes {0}".format(shapes))
        self.ny, self.nx = arrays["XC"].shape
        self.n = self.nx * self.ny
        self.radius = float(radius)
        self.lon, self.lat = arrays["XC"].ravel(), arrays["YC"].ravel()
        self.glon, self.glat = arrays["XG"].ravel(), arrays["YG"].ravel()
        self.area = arrays["RAC"].ravel()
        with np.errstate(invalid="ignore"):
            self.valid = (np.isfinite(self.lon) & np.isfinite(self.lat)
                          & np.isfinite(self.glon) & np.isfinite(self.glat)
                          & np.isfinite(self.area) & (self.area > 0))
        if (np.abs(self.lat[self.valid]) > 90.0 + 1e-9).any() or \
                (np.abs(self.glat[self.valid]) > 90.0 + 1e-9).any():
            raise BuildError("grid: YC or YG has values outside [-90, 90]; the builder "
                             "needs a spherical grid in degrees (Cartesian grids are "
                             "not supported)")
        hf = arrays["hFacC"].ravel()
        with np.errstate(invalid="ignore"):
            self.wet = self.valid & (hf > 0)
        self.wet_cells = np.nonzero(self.wet)[0]
        self.wet_pos = np.full(self.n, -1, dtype=np.int64)
        self.wet_pos[self.wet_cells] = np.arange(self.wet_cells.size)
        self._centers = self._corners = None

    @property
    def centers(self):
        """Unit vectors of every cell center, ``(n, 3)``."""
        if self._centers is None:
            self._centers = _unit(self.lon, self.lat)
        return self._centers

    @property
    def corners(self):
        """Unit vectors of every SW corner, ``(n, 3)``."""
        if self._corners is None:
            self._corners = _unit(self.glon, self.glat)
        return self._corners

    def describe(self, c):
        c = int(c)
        return "cell {0} (i={1}, j={2}; lon {3:.6g}, lat {4:.6g})".format(
            c, c % self.nx, c // self.nx, self.lon[c], self.lat[c])


# ---------------------------------------------------------------------------
# Neighbour graph


@dataclass
class WetGraph:
    """Edge-neighbour graph of the wet surface cells (module notes).

    ``cells[k]`` is the global cell index of node ``k`` (ascending);
    ``neighbours[k]`` holds up to 4 node indices (-1 for none) and
    ``distances[k]`` the great-circle distances (m) between the centers.
    ``diagnostics`` counts how the corners were found (``corners`` method).
    """
    cells: np.ndarray
    neighbours: np.ndarray
    distances: np.ndarray
    connectivity: str
    periodic: bool = False
    diagnostics: dict = field(default_factory=dict)

    def edges(self):
        """Node pairs ``(a, b)`` with ``a < b``, one per edge, as two arrays."""
        a = np.repeat(np.arange(len(self.cells)), 4)
        b = self.neighbours.ravel()
        keep = (b >= 0) & (a < b)
        return a[keep], b[keep]

    def degree(self):
        return np.count_nonzero(self.neighbours >= 0, axis=1)


def wet_graph(grid, *, connectivity="auto", earth_radius=EARTH_RADIUS):
    """Edge-neighbour graph of the wet surface cells of ``grid`` (a grid
    directory or a dict of arrays, as for :func:`build_targets`)."""
    geom = grid if isinstance(grid, _Geom) else _Geom(grid, earth_radius)
    return _build_graph(geom, connectivity)


def _is_mosaic(geom):
    """True when an array neighbour's SW corner is far from a cell (exch2 layout)."""
    nx, ny = geom.nx, geom.ny
    P, G, v = geom.centers, geom.corners, geom.valid
    h = _norm(G - P)
    c = np.arange(geom.n)
    for step, ok in ((1, (c % nx) + 1 < nx), (nx, (c // nx) + 1 < ny)):
        a = c[ok]
        b = a + step
        use = v[a] & v[b] & (h[a] > 0)
        a, b = a[use], b[use]
        if a.size and (_norm(G[b] - P[a]) > MOSAIC_RATIO * h[a]).any():
            return True
    return False


def _build_graph(geom, connectivity):
    """:class:`WetGraph` of ``geom``: node pairs from :func:`_index_pairs` or
    :func:`_corner_pairs` (``auto``: :func:`_is_mosaic` decides), without
    duplicates or self-pairs, stored as up to 4 neighbours per node with the
    great-circle distance between the two centers. ``index`` on an exch2
    layout is refused."""
    if connectivity not in CONNECTIVITIES:
        raise BuildError("connectivity {0!r} is not one of {1}".format(
            connectivity, ", ".join(CONNECTIVITIES)))
    mosaic = _is_mosaic(geom)
    if connectivity == "index" and mosaic:
        raise BuildError("connectivity 'index': array neighbours are not grid neighbours "
                         "on this grid (an exch2 layout); use 'corners' or 'auto'")
    method = connectivity if connectivity != "auto" else ("corners" if mosaic else "index")
    diag = {}
    if method == "index":
        pairs, periodic = _index_pairs(geom)
    else:
        pairs, diag = _corner_pairs(geom)
        periodic = False
    nw = geom.wet_cells.size
    nbr = np.full((nw, 4), -1, dtype=np.int64)
    dist = np.full((nw, 4), np.nan)
    if pairs.size:
        pairs = np.unique(np.sort(pairs, axis=1), axis=0)
        pairs = pairs[pairs[:, 0] != pairs[:, 1]]
        a = np.concatenate([pairs[:, 0], pairs[:, 1]])
        b = np.concatenate([pairs[:, 1], pairs[:, 0]])
        order = np.lexsort((b, a))
        a, b = a[order], b[order]
        start = np.searchsorted(a, np.arange(nw))
        slot = np.arange(a.size) - start[a]
        if (slot >= 4).any():
            raise BuildError("grid: a wet cell has more than 4 edge neighbours; the grid "
                             "geometry is inconsistent")
        ca, cb = geom.wet_cells[a], geom.wet_cells[b]
        nbr[a, slot] = b
        dist[a, slot] = _haversine(geom.lon[ca], geom.lat[ca], geom.lon[cb], geom.lat[cb],
                                   geom.radius)
    return WetGraph(geom.wet_cells, nbr, dist, method, periodic, diag)


def _index_pairs(geom):
    """Wet node pairs of (i+1, j) and (i, j+1) neighbours, with the zonal wrap."""
    nx, ny = geom.nx, geom.ny
    rows = np.arange(ny)
    first, last = rows * nx, rows * nx + nx - 1
    periodic = False
    if nx > 2 and geom.valid[first].all() and geom.valid[last].all():
        P, G = geom.centers, geom.corners
        ratio = _norm(G[first] - P[last]) / _norm(G[last] - P[last])
        periodic = bool(np.all(np.abs(ratio - 1.0) <= WRAP_TOL))
    w = geom.wet_cells
    i, j = w % nx, w // nx
    east = np.where(i + 1 < nx, w + 1, (w - (nx - 1)) if periodic else -1)
    north = np.where(j + 1 < ny, w + nx, -1)
    out = []
    for nb in (east, north):
        ok = nb >= 0
        ok[ok] = geom.wet[nb[ok]]
        out.append(np.stack([geom.wet_pos[w[ok]], geom.wet_pos[nb[ok]]], axis=1))
    return np.concatenate(out), periodic


def _vertex_ids(geom):
    """Vertex id of each cell's SW corner: the lowest cell index among valid
    cells whose SW-corner unit vectors round to the same multiples of
    VERTEX_EPS (the same point); -1 for invalid cells."""
    vid = np.full(geom.n, -1, dtype=np.int64)
    cells = np.nonzero(geom.valid)[0]
    q = np.round(geom.corners[cells] / VERTEX_EPS).astype(np.int64)
    _, first, inverse = np.unique(q, axis=0, return_index=True, return_inverse=True)
    vid[cells] = cells[first][np.ravel(inverse)]
    return vid


def _geometric_corners(geom, vid, cells):
    """Corner vertex ids of ``cells`` found among all SW corners (module notes).

    Returns ``(ids, synth, scale)``. ``ids`` is ``(len(cells), 4)`` in cyclic
    order: SW, one of SE/NW, NE, the other of SE/NW. A corner that is no
    cell's SW corner has id -1 and its predicted unit vector in
    ``synth[(row, slot)]``. ``scale`` is each cell's estimated shortest side
    (unit sphere), from its SW-NE diagonal and area as for a rectangle.
    """
    P, G = geom.centers, geom.corners
    valid_cells = np.nonzero(geom.valid)[0]
    index = _PointIndex(G[valid_cells])
    rows = np.arange(len(cells))
    pr, sr = P[cells], G[cells]
    ne_pred = _reflect(pr, sr)
    diag2 = np.sum((sr - ne_pred) ** 2, axis=1)
    area = geom.area[cells] / geom.radius ** 2
    m_est = np.minimum(np.sqrt(area), 0.5 * (np.sqrt(diag2 + 2 * area)
                                             - np.sqrt(np.maximum(diag2 - 2 * area, 0))))
    own = vid[cells]
    d_ne, k_ne = index.query(ne_pred, 1)
    ne_ok = d_ne[:, 0] <= CORNER_TOL * m_est
    ne_id = np.where(ne_ok, vid[valid_cells[k_ne[:, 0]]], -1)
    d_c, k_c = index.query(pr, 8)
    cand = vid[valid_cells[k_c]]
    usable = (cand != own[:, None]) & (cand != ne_id[:, None])
    if not usable.any(axis=1).all():
        bad = cells[~usable.any(axis=1)][0]
        raise BuildError("grid: can't find the corners of {0}".format(geom.describe(bad)))
    first = np.argmax(usable, axis=1)
    u_cell = valid_cells[k_c[rows, first]]
    u_id = cand[rows, first]
    w_pred = _reflect(pr, G[u_cell])
    d_w, k_w = index.query(w_pred, 1)
    w_id = vid[valid_cells[k_w[:, 0]]]
    w_ok = ((d_w[:, 0] <= CORNER_TOL * m_est) & (w_id != own) & (w_id != ne_id)
            & (w_id != u_id))
    ids = np.stack([own, u_id, np.where(ne_ok, ne_id, -1), np.where(w_ok, w_id, -1)], axis=1)
    synth = {}
    for r in np.nonzero(~ne_ok)[0].tolist():
        synth[(r, 2)] = ne_pred[r]
    for r in np.nonzero(~w_ok)[0].tolist():
        synth[(r, 3)] = w_pred[r]
    return ids, synth, m_est


def _corner_pairs(geom):
    """Wet node pairs sharing a cell edge, from corner vertices (module notes).

    Cells whose array neighbours pass the quadrilateral test take their
    corners from the array; the others from :func:`_geometric_corners`.
    Returns the pairs and the counts reported in ``WetGraph.diagnostics``.
    """
    nx, ny = geom.nx, geom.ny
    P, G = geom.centers, geom.corners
    vid = _vertex_ids(geom)
    w = geom.wet_cells
    nw = w.size
    i, j = w % nx, w // nx
    ok = (i + 1 < nx) & (j + 1 < ny)
    e = np.where(ok, w + 1, w)
    n = np.where(ok, w + nx, w)
    ne = np.where(ok, w + nx + 1, w)
    ok &= geom.valid[e] & geom.valid[n] & geom.valid[ne]
    p, s = P[w], G[w]
    side = np.min(np.stack([_norm(G[e] - s), _norm(G[n] - s), _norm(G[ne] - G[e]),
                            _norm(G[ne] - G[n])], axis=1), axis=1)
    with np.errstate(invalid="ignore"):
        quad = (ok & (side > 0)
                & (_norm(_reflect(p, s) - G[ne]) <= CORNER_TOL * side)
                & (_norm(_reflect(p, G[e]) - G[n]) <= CORNER_TOL * side))
    ids = np.stack([vid[w], vid[e], vid[ne], vid[n]], axis=1)   # cyclic: SW, SE, NE, NW
    synth = {}            # (node, slot) -> synthesized unit vector
    scale = {}            # node -> estimated shortest side (unit sphere)
    rest = np.nonzero(~quad)[0]
    if rest.size:
        gids, gsynth, m_est = _geometric_corners(geom, vid, w[rest])
        ids[rest] = gids
        for (r, slot), pos in gsynth.items():
            synth[(int(rest[r]), slot)] = pos
        scale = dict(zip(rest.tolist(), m_est.tolist()))
    # Edges: consecutive corners. Vertex-vertex edges match exactly.
    nodes = np.arange(nw)
    keys, owners = [], []
    half = {}             # vertex id -> [(node, synthesized position)]
    for s0, s1 in ((0, 1), (1, 2), (2, 3), (3, 0)):
        a, b = ids[:, s0], ids[:, s1]
        both = (a >= 0) & (b >= 0) & (a != b)
        keys.append(np.minimum(a, b)[both] * geom.n + np.maximum(a, b)[both])
        owners.append(nodes[both])
        for node in np.nonzero((a >= 0) != (b >= 0))[0].tolist():
            slot = s1 if a[node] >= 0 else s0
            vert = int(a[node] if a[node] >= 0 else b[node])
            half.setdefault(vert, []).append((node, synth[(node, slot)]))
    keys, owners = np.concatenate(keys), np.concatenate(owners)
    order = np.argsort(keys, kind="stable")
    keys, owners = keys[order], owners[order]
    starts = np.flatnonzero(np.r_[True, keys[1:] != keys[:-1]])
    counts = np.diff(np.r_[starts, keys.size])
    if (counts > 2).any():
        k = starts[np.argmax(counts > 2)]
        raise BuildError("grid: the edge between vertices of cells {0} and {1} is shared "
                         "by more than two wet cells; the grid geometry is "
                         "inconsistent".format(keys[k] // geom.n, keys[k] % geom.n))
    two = starts[counts == 2]
    pairs = [np.stack([owners[two], owners[two + 1]], axis=1)]
    # Vertex-to-synthesized edges: same vertex, synthesized corners close together.
    matched = unmatched = 0
    extra = []
    for vert, items in half.items():
        used = set()
        cand = sorted(
            (_norm(pa - pb) / min(scale[na], scale[nb]), na, nb)
            for x, (na, pa) in enumerate(items) for (nb, pb) in items[x + 1:])
        for sep, na, nb in cand:
            if sep <= SYNTH_MERGE_TOL and na not in used and nb not in used:
                used.update((na, nb))
                extra.append((na, nb))
                matched += 1
        unmatched += len(items) - len(used)
    if extra:
        pairs.append(np.array(extra, dtype=np.int64))
    diag = {"array_corner_cells": int(np.count_nonzero(quad)),
            "geometric_corner_cells": int(rest.size),
            "synthesized_corners": len(synth),
            "synthesized_edges": matched,
            "unmatched_synthesized_edges": unmatched}
    return np.concatenate(pairs), diag


# ---------------------------------------------------------------------------
# Snapping and spreading


def _snap(src, geom, opts):
    """Nearest wet cell of each source and its great-circle distance (m)."""
    if geom.wet_cells.size == 0:
        raise BuildError("grid: no wet surface cell (hFacC > 0 at level 1 and RAC > 0)")
    index = _PointIndex(geom.centers[geom.wet_cells])
    q = _unit(src.lon, src.lat)
    d, k = index.query(q, 4)
    tie = d <= d[:, :1] * (1.0 + 1e-12) + 1e-15
    cells = np.where(tie, geom.wet_cells[k], np.iinfo(np.int64).max).min(axis=1)
    dist = _haversine(src.lon, src.lat, geom.lon[cells], geom.lat[cells], geom.radius)
    errors = []
    for s in range(len(src.ids)):
        if dist[s] > opts[s]["max_snap"]:
            errors.append("source {0!r} ({1}): the nearest wet surface cell is {2}, "
                          "{3:.6g} km away, farther than max_snap_distance = {4:.6g} "
                          "km".format(src.ids[s], src.where[s], geom.describe(cells[s]),
                                      dist[s] / 1000.0, opts[s]["max_snap"] / 1000.0))
    if errors:
        raise BuildError(errors)
    return cells, dist


def _dijkstra(graph, start, cutoff):
    """``{node: r}`` of every node reached with ``r <= cutoff``, starting at
    ``start`` with ``r = 0`` (bounded Dijkstra)."""
    nbr, dst = graph.neighbours, graph.distances
    done, best = {}, {start: 0.0}
    heap = [(0.0, start)]
    pop, push = heapq.heappop, heapq.heappush
    while heap:
        d, v = pop(heap)
        if v in done:
            continue
        done[v] = d
        for u, step in zip(nbr[v].tolist(), dst[v].tolist()):
            if u < 0 or u in done:
                continue
            nd = d + step
            if nd <= cutoff and nd < best.get(u, math.inf):
                best[u] = nd
                push(heap, (nd, u))
    return done


def _normalize(weights):
    """Fractions summing to 1: divide by the exact sum, then give the rounding
    residual to the largest fraction."""
    f = weights / math.fsum(weights)
    k = int(np.argmax(f))
    f[k] += 1.0 - math.fsum(f)
    return f


@dataclass
class TargetTables:
    """Source, alias and target tables built by :func:`build_targets`.

    Arrays are per source (``source_*``), per alias (``alias_*``) or per
    target (``target_*``, sorted by source then cell). ``source_meta`` maps
    ``source_name``/``source_type``/``source_notes``/``source_reference`` to
    lists of strings (only those given). ``source_spread_type`` is ``"none"``
    and ``source_spread_scale``/``source_cutoff`` are NaN for pointwise
    sources. ``source_snap_cell`` is the snapped cell of each source.
    """
    nx: int
    ny: int
    source_id: list
    source_lon: np.ndarray
    source_lat: np.ndarray
    source_meta: dict
    alias_source: np.ndarray
    alias_name: list
    alias_scheme: list
    source_snap_cell: np.ndarray
    source_snap_distance: np.ndarray
    source_emission: list
    source_spread_type: list
    source_spread_scale: np.ndarray
    source_cutoff: np.ndarray
    target_source: np.ndarray
    target_cell: np.ndarray
    target_fraction: np.ndarray
    target_cell_area: np.ndarray
    target_lon: np.ndarray
    target_lat: np.ndarray
    target_distance: np.ndarray
    connectivity: str = None
    periodic: bool = False
    earth_radius: float = EARTH_RADIUS

    @property
    def comment(self):
        return _KERNEL_COMMENT.format(builder=BUILDER, version=BUILDER_VERSION,
                                      radius=self.earth_radius)


def build_targets(sources, grid, *, emission=DEFAULT_EMISSION, spread_type=None,
                  spread_scale=None, cutoff=None, max_snap_distance=DEFAULT_MAX_SNAP_DISTANCE,
                  earth_radius=EARTH_RADIUS, connectivity="auto"):
    """Build the source, alias and target tables (module notes, algorithm).

    Parameters
    ----------
    sources : str, os.PathLike or list of dict
        CSV or schema-1.0 NetCDF source table, or records with the CSV column
        names (see :func:`read_sources`).
    grid : str, os.PathLike or dict
        MITgcm grid output directory, or a dict of global 2D arrays
        ``hFacC`` (2D, or 3D whose level 1 is used), ``XC``, ``YC``, ``XG``,
        ``YG`` and ``RAC``.
    emission, spread_type, spread_scale, cutoff, max_snap_distance
        Defaults for sources that don't set their own. Distances are meters or
        strings such as ``"25km"``. ``cutoff`` defaults to ``3 * spread_scale``
        (exponential, gaussian); for linear the cutoff is ``Rcut``, or
        ``cutoff`` if smaller.
    earth_radius : float or str
        Sphere radius for great-circle distances, m.
    connectivity : {"auto", "index", "corners"}
        Neighbour-graph method (module notes).

    Returns
    -------
    TargetTables

    Raises
    ------
    BuildError
        Invalid sources, options or grid; a source farther than its
        ``max_snap_distance`` from every wet cell. ``errors`` lists all.
    OSError
        The source table or a grid field can't be read.
    """
    try:
        radius = parse_distance(earth_radius)
    except ValueError as e:
        raise BuildError("earth_radius: {0}".format(e))
    if not (math.isfinite(radius) and radius > 0):
        raise BuildError("earth_radius = {0!r} must be a positive distance".format(
            earth_radius))
    if connectivity not in CONNECTIVITIES:
        raise BuildError("connectivity {0!r} is not one of {1}".format(
            connectivity, ", ".join(CONNECTIVITIES)))
    src = _load_sources(sources)
    opts = _resolve_options(src, {"emission": emission, "spread_type": spread_type,
                                  "spread_scale": spread_scale, "cutoff": cutoff,
                                  "max_snap_distance": max_snap_distance})
    geom = _Geom(grid, radius)
    snap_cell, snap_dist = _snap(src, geom, opts)
    graph = None
    if any(o["emission"] == "spread" for o in opts):
        graph = _build_graph(geom, connectivity)
    ns = len(src.ids)
    t_src, t_cell, t_frac, t_dist = [], [], [], []
    for s in range(ns):
        o = opts[s]
        if o["emission"] == "pointwise":
            cells = np.array([snap_cell[s]])
            frac = np.array([1.0])
            r = np.array([0.0])          # the snapped cell itself
        else:
            reached = _dijkstra(graph, int(geom.wet_pos[snap_cell[s]]), o["cutoff"])
            nodes = np.fromiter(reached.keys(), dtype=np.int64, count=len(reached))
            r = np.fromiter(reached.values(), dtype=np.float64, count=len(reached))
            cells = geom.wet_cells[nodes]
            wts = kernel_weight(r, o["spread_type"], o["spread_scale"]) * geom.area[cells]
            # The snapped cell has r = 0, W = 1 and RAC > 0, so at least it is kept.
            keep = wts > 0
            order = np.argsort(cells[keep], kind="stable")
            cells, r = cells[keep][order], r[keep][order]
            frac = _normalize(wts[keep][order])
        t_src.append(np.full(cells.size, s, dtype=np.int64))
        t_cell.append(cells)
        t_frac.append(frac)
        t_dist.append(r)
    t_cell = np.concatenate(t_cell)
    cell_type = np.int32 if geom.n - 1 <= np.iinfo(np.int32).max else np.int64
    spread = [o["emission"] == "spread" for o in opts]
    return TargetTables(
        nx=geom.nx, ny=geom.ny, source_id=list(src.ids),
        source_lon=src.lon, source_lat=src.lat, source_meta=dict(src.meta),
        alias_source=np.array([a[0] for a in src.aliases], dtype=np.int32),
        alias_name=[a[1] for a in src.aliases], alias_scheme=[a[2] for a in src.aliases],
        source_snap_cell=np.asarray(snap_cell, dtype=np.int64),
        source_snap_distance=np.asarray(snap_dist, dtype=np.float64),
        source_emission=[o["emission"] for o in opts],
        source_spread_type=[o["spread_type"] if sp else "none" for o, sp in zip(opts, spread)],
        source_spread_scale=np.array([o["spread_scale"] if sp else np.nan
                                      for o, sp in zip(opts, spread)], dtype=np.float64),
        source_cutoff=np.array([o["cutoff"] if sp else np.nan for o, sp in zip(opts, spread)],
                               dtype=np.float64),
        target_source=np.concatenate(t_src).astype(np.int32),
        target_cell=t_cell.astype(cell_type),
        target_fraction=np.concatenate(t_frac),
        target_cell_area=geom.area[t_cell].copy(),
        target_lon=geom.lon[t_cell].copy(), target_lat=geom.lat[t_cell].copy(),
        target_distance=np.concatenate(t_dist),
        connectivity=None if graph is None else graph.connectivity,
        periodic=False if graph is None else graph.periodic,
        earth_radius=radius)


# ---------------------------------------------------------------------------
# Writing

#: Provenance variables (dimension ``source``) written by :func:`write_targets`.
PROVENANCE_VARIABLES = ("source_snap_distance", "source_emission", "source_spread_type",
                        "source_spread_scale", "source_cutoff")
_COPY_BLOCK_BYTES = 50 * 2**20


def _chars(strings, strlen):
    """``(n, strlen)`` ``S1`` array of NUL-padded UTF-8 strings."""
    enc = [s.encode("utf-8") for s in strings]
    return np.array(enc, dtype="S{0}".format(max(strlen, 1))).view("S1").reshape(
        len(enc), max(strlen, 1))


def _put(ds, name, dims, values, dtype, attrs):
    v = ds.createVariable(name, dtype, dims)
    v.set_auto_maskandscale(False)
    v.setncatts({k: a for k, a in attrs.items() if a is not None})
    if np.size(values):
        v[:] = np.asarray(values).astype(dtype)
    return v


def _put_strings(ds, name, dim, values, attrs):
    """A string variable: variable-length in NETCDF4, else a char array."""
    if ds.data_model == "NETCDF4":
        v = ds.createVariable(name, str, (dim,))
        if len(values):
            v[:] = np.array(list(values), dtype=object)
    else:
        n = max([len(s.encode("utf-8")) for s in values] + [1])
        sdim = name + "_strlen"
        ds.createDimension(sdim, n)
        v = ds.createVariable(name, "S1", (dim, sdim))
        v.set_auto_chartostring(False)
        if len(values):
            v[:] = _chars(values, n)
    v.setncatts(attrs)


def _write_tables(ds, t, write_ids):
    """Define and write the builder's variables (``source_id`` if ``write_ids``)."""
    if write_ids:
        strlen = max([len(s) for s in t.source_id] + [1])
        ds.createDimension(S.DIM_SOURCE, len(t.source_id))
        ds.createDimension("id_strlen", strlen)
        v = ds.createVariable("source_id", "S1", (S.DIM_SOURCE, "id_strlen"))
        v.set_auto_chartostring(False)
        v[:] = _chars(t.source_id, strlen)
        v.setncatts({"long_name": "source identifier", "cf_role": "timeseries_id"})
    meta_attrs = {
        "source_name": {"long_name": "primary name of the source"},
        "source_type": {"long_name": "kind of source",
                        "comment": "river, glacier, ice_sheet_basin, iceberg_melt, "
                                   "groundwater or other"},
        "source_notes": {"long_name": "notes on provenance and processing"},
        "source_reference": {"long_name": "citation, DOI or URL of the source data"},
    }
    for name, values in t.source_meta.items():
        _put_strings(ds, name, S.DIM_SOURCE, values, meta_attrs[name])
    _put(ds, "source_lon", (S.DIM_SOURCE,), t.source_lon, "f8",
         {"long_name": "longitude of the mouth or terminus", "standard_name": "longitude",
          "units": "degrees_east"})
    _put(ds, "source_lat", (S.DIM_SOURCE,), t.source_lat, "f8",
         {"long_name": "latitude of the mouth or terminus", "standard_name": "latitude",
          "units": "degrees_north"})
    _put(ds, "source_snap_distance", (S.DIM_SOURCE,), t.source_snap_distance, "f8",
         {"long_name": "great-circle distance from the source to the center of its "
                       "snapped (nearest wet) cell", "units": "m"})
    _put_strings(ds, "source_emission", S.DIM_SOURCE, t.source_emission,
                 {"long_name": "emission type used to build the targets",
                  "comment": "pointwise or spread"})
    _put_strings(ds, "source_spread_type", S.DIM_SOURCE, t.source_spread_type,
                 {"long_name": "spread kernel",
                  "comment": "gaussian, exponential or linear; none for pointwise sources"})
    _put(ds, "source_spread_scale", (S.DIM_SOURCE,), t.source_spread_scale, "f8",
         {"long_name": "spread scale X, the distance at which the kernel falls to 1/e of "
                       "its peak", "units": "m", "comment": "NaN for pointwise sources"})
    _put(ds, "source_cutoff", (S.DIM_SOURCE,), t.source_cutoff, "f8",
         {"long_name": "largest distance r of a target cell", "units": "m",
          "comment": "NaN for pointwise sources"})
    if len(t.alias_name):
        ds.createDimension(S.DIM_ALIAS, len(t.alias_name))
        _put(ds, "alias_source", (S.DIM_ALIAS,), t.alias_source, "i4",
             {"long_name": "index of the source this alias names",
              "instance_dimension": S.DIM_SOURCE})
        _put_strings(ds, "alias_name", S.DIM_ALIAS, t.alias_name,
                     {"long_name": "alternative name or catalogue id"})
        if any(t.alias_scheme):
            _put_strings(ds, "alias_scheme", S.DIM_ALIAS,
                         [s or "" for s in t.alias_scheme],
                         {"long_name": "naming system or authority of the alias"})
    ds.createDimension(S.DIM_TARGET, len(t.target_cell))
    _put(ds, "target_source", (S.DIM_TARGET,), t.target_source, "i4",
         {"long_name": "index of the source feeding this target",
          "instance_dimension": S.DIM_SOURCE})
    _put(ds, "target_cell", (S.DIM_TARGET,), t.target_cell, t.target_cell.dtype.str[1:],
         {"long_name": "0-based global cell index",
          "comment": "cell = i + mitgcm_grid_nx * j, 0-based (i, j) in the global 2D "
                     "layout of a dense runoffFile"})
    _put(ds, "target_fraction", (S.DIM_TARGET,), t.target_fraction, "f8",
         {"long_name": "share of the source flux sent to this cell", "units": "1",
          "comment": t.comment})
    _put(ds, "target_cell_area", (S.DIM_TARGET,), t.target_cell_area, "f8",
         {"long_name": "horizontal cell area rA", "units": "m2"})
    _put(ds, "target_lon", (S.DIM_TARGET,), t.target_lon, "f8",
         {"long_name": "cell-center longitude", "standard_name": "longitude",
          "units": "degrees_east"})
    _put(ds, "target_lat", (S.DIM_TARGET,), t.target_lat, "f8",
         {"long_name": "cell-center latitude", "standard_name": "latitude",
          "units": "degrees_north"})
    _put(ds, "target_distance", (S.DIM_TARGET,), t.target_distance, "f8",
         {"long_name": "distance r from the snapped cell along edge-connected wet cells "
                       "(0 at the snapped cell)", "units": "m"})


def _replaced(t):
    """Source variables a write into an existing file replaces."""
    return {"source_lon", "source_lat", *PROVENANCE_VARIABLES, *t.source_meta}


def _history_line(history):
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    return "{0}: {1}".format(stamp, history or "created by {0}.write_targets".format(BUILDER))


def write_targets(path, tables, *, into=None, history=None, attrs=None):
    """Write ``tables`` (from :func:`build_targets`) to ``path``; returns ``path``.

    Without ``into``, a new NETCDF4 file holds only the tables: ``source_id``,
    ``source_lon``/``source_lat``, the ``source_*`` metadata that was given,
    the alias table (from ``alt_names``), the target table with
    ``target_cell_area``, ``target_lon``, ``target_lat`` and the user
    variable ``target_distance`` (m), and the per-source provenance
    ``source_snap_distance`` (m), ``source_emission``, ``source_spread_type``,
    ``source_spread_scale`` (m) and ``source_cutoff`` (m). Global attributes:
    ``Conventions``, ``mitgcm_runoff_schema_version``, ``mitgcm_grid_nx``/``ny``,
    ``title``, ``source`` (builder and version), ``history`` (``history``, or
    a default, with a UTC time stamp), ``comment`` (the kernel convention,
    also on ``target_fraction``) and ``date_created``; ``attrs`` adds more
    (e.g. ``mitgcm_grid_name``). Check it with ``check_files(path,
    tables_only=True)``.

    With ``into``, an existing runoff file whose ``source_id`` list equals
    ``tables.source_id`` (same ids in the same order; otherwise
    :class:`BuildError`), ``path`` becomes a copy of it in which the target
    table, the alias table and the source variables written by the builder are
    replaced; every other dimension, variable (time series included, copied in
    bounded blocks with their chunking, deflate, shuffle and fletcher32
    settings) and attribute is kept. ``mitgcm_grid_nx``/``ny`` are set from the
    tables, the history line is prepended, and ``comment`` is set only if the
    file had none. ``path`` may equal ``into``: the file is written to a
    temporary file in the same directory and then renamed over it.
    """
    import netCDF4
    path = os.fspath(path)
    t = tables
    line = _history_line(history)
    extra = {k: v for k, v in (attrs or {}).items() if v is not None}
    directory = os.path.dirname(os.path.abspath(path))
    fd, tmp = tempfile.mkstemp(prefix=".targets-", suffix=".nc", dir=directory)
    os.close(fd)
    try:
        if into is None:
            with netCDF4.Dataset(tmp, "w", format="NETCDF4") as ds:
                ds.set_ncstring_attrs(False)
                ds.setncatts(_new_attrs(t, line, extra))
                _write_tables(ds, t, write_ids=True)
        else:
            _write_into(os.fspath(into), tmp, t, line, extra)
        os.replace(tmp, path)
    except BaseException:
        if os.path.exists(tmp):
            os.remove(tmp)
        raise
    return path


def _new_attrs(t, line, extra):
    g = {"Conventions": S.CONVENTIONS,
         "mitgcm_runoff_schema_version": S.SCHEMA_VERSION,
         "mitgcm_grid_nx": np.int32(t.nx), "mitgcm_grid_ny": np.int32(t.ny),
         "title": "Sparse runoff source and target tables",
         "source": "{0} {1}".format(BUILDER, BUILDER_VERSION),
         "history": line,
         "comment": t.comment,
         "references": "docs/runoff_schema.md (sparse runoff schema 1.0), section 13",
         "date_created": line.split(": ", 1)[0]}
    g.update(extra)
    return g


def _write_into(into, tmp, t, line, extra):
    """Write ``into`` to ``tmp`` with its tables replaced by ``t`` (see
    :func:`write_targets`): check the ``source_id`` list, copy every kept
    dimension, variable and attribute, then write the builder's tables.
    Dimensions used only by replaced variables (``target``, ``alias``, a
    string length) are dropped."""
    import netCDF4
    with netCDF4.Dataset(into, "r") as src:
        src.set_auto_maskandscale(False)
        src.set_auto_chartostring(False)
        if src.groups:
            raise BuildError("{0}: groups are not supported".format(into))
        if "source_id" not in src.variables:
            raise BuildError("{0}: no source_id variable".format(into))
        raw = np.asarray(src.variables["source_id"][:])
        ids = [row.tobytes().rstrip(b"\x00").rstrip(b" ").decode("utf-8", errors="replace")
               for row in raw.reshape(raw.shape[0], -1)]
        if ids != list(t.source_id):
            diff = next((k for k, (a, b) in enumerate(zip(ids, t.source_id)) if a != b),
                        min(len(ids), len(t.source_id)))
            raise BuildError(
                "{0}: its source_id list differs from the tables' ({1} vs {2} sources; "
                "first difference at index {3}: {4!r} there, {5!r} here); write into a "
                "file only with the same sources in the same order".format(
                    into, len(ids), len(t.source_id), diff,
                    ids[diff] if diff < len(ids) else None,
                    t.source_id[diff] if diff < len(t.source_id) else None))
        drop_dims = {S.DIM_TARGET, S.DIM_ALIAS}
        replaced = _replaced(t)
        drop_vars = {name for name, v in src.variables.items()
                     if drop_dims & set(v.dimensions) or name in replaced}
        with netCDF4.Dataset(tmp, "w", format=src.data_model) as dst:
            dst.set_ncstring_attrs(False)
            g = {a: src.getncattr(a) for a in src.ncattrs()}
            old = g.get("history")
            g["history"] = line if not old else line + "\n" + str(old)
            g["mitgcm_grid_nx"], g["mitgcm_grid_ny"] = np.int32(t.nx), np.int32(t.ny)
            g.setdefault("comment", t.comment)
            g.update(extra)
            dst.setncatts(g)
            # dimensions used only by dropped variables (e.g. a string length) go too
            kept_dims, dropped_dims = set(), set(drop_dims)
            for name, v in src.variables.items():
                (dropped_dims if name in drop_vars else kept_dims).update(v.dimensions)
            for name, d in src.dimensions.items():
                if name in drop_dims or (name in dropped_dims and name not in kept_dims):
                    continue
                dst.createDimension(name, None if d.isunlimited() else len(d))
            for name, v in src.variables.items():
                if name not in drop_vars:
                    _copy_variable(v, dst, name)
            _write_tables(dst, t, write_ids=False)


def _copy_variable(v, dst, name):
    """Copy one variable with its storage settings and attributes, in blocks."""
    kw = {}
    is_str = v.dtype is str
    try:
        chunks = v.chunking()
    except Exception:  # noqa: BLE001 - NETCDF3 files have no chunking information
        chunks = None
    if isinstance(chunks, (list, tuple)):
        kw["chunksizes"] = tuple(chunks)
    try:
        filters = v.filters() or {}
    except Exception:  # noqa: BLE001 - NETCDF3
        filters = {}
    if filters.get("zlib"):
        kw.update(zlib=True, complevel=filters.get("complevel") or 4)
    if filters.get("shuffle"):
        kw["shuffle"] = True
    if filters.get("fletcher32"):
        kw["fletcher32"] = True
    if not is_str:
        if "_FillValue" in v.ncattrs():
            kw["fill_value"] = v.getncattr("_FillValue")
        kw["endian"] = v.endian()
    out = dst.createVariable(name, str if is_str else v.datatype, v.dimensions, **kw)
    out.set_auto_maskandscale(False)
    out.set_auto_chartostring(False)
    out.setncatts({a: v.getncattr(a) for a in v.ncattrs() if a != "_FillValue"})
    if v.ndim == 0:
        out.assignValue(v.getValue())
        return
    n0 = v.shape[0]
    if n0 == 0:
        return
    if is_str:
        out[:] = v[:]
        return
    row = np.dtype(v.dtype).itemsize * int(np.prod(v.shape[1:], dtype=np.int64))
    step = max(1, _COPY_BLOCK_BYTES // max(row, 1))
    for a in range(0, n0, step):
        b = min(n0, a + step)     # an unlimited dimension needs the exact slice
        out[a:b] = v[a:b]


# ---------------------------------------------------------------------------
# Command line


def main(argv=None):
    """Command-line entry point; returns the exit status."""
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = argparse.ArgumentParser(
        prog="python -m MITgcmutils.runoff.targets",
        description="Build the source and target tables of a sparse-runoff file "
                    "(schema {0}) from source locations and MITgcm grid output. Exit 0: "
                    "success; 1: invalid input or checker errors; 2: usage or I/O "
                    "problem. Distances are meters or e.g. 25km.".format(S.SCHEMA_VERSION))
    parser.add_argument("sources", metavar="SOURCES",
                        help="CSV source table, or a schema-1.0 NetCDF file")
    parser.add_argument("--grid-dir", required=True, metavar="DIR",
                        help="MITgcm grid output (hFacC, XC, YC, XG, YG, RAC)")
    parser.add_argument("-o", "--output", required=True, metavar="OUT",
                        help="output NetCDF file")
    parser.add_argument("--into", metavar="RUNOFF",
                        help="copy this runoff file (same source_id list) to OUT with its "
                             "tables replaced; OUT may be the same file")
    parser.add_argument("--emission", choices=EMISSIONS, default=DEFAULT_EMISSION)
    parser.add_argument("--spread-type", choices=SPREAD_TYPES)
    parser.add_argument("--spread-scale", metavar="X",
                        help="distance at which the kernel falls to 1/e")
    parser.add_argument("--cutoff", metavar="C",
                        help="largest distance r of a target (default 3X; linear: "
                             "X/(1-1/e), or C if smaller)")
    parser.add_argument("--max-snap-distance", metavar="D", default="50km",
                        help="largest source-to-cell distance (default 50km)")
    parser.add_argument("--earth-radius", metavar="R", default=repr(EARTH_RADIUS))
    parser.add_argument("--connectivity", choices=CONNECTIVITIES, default="auto",
                        help="neighbour-graph method (default auto)")
    parser.add_argument("--grid-name", metavar="NAME", help="mitgcm_grid_name attribute")
    parser.add_argument("--no-check", action="store_true",
                        help="don't run the checker on the output")
    try:
        args = parser.parse_args(argv)
    except SystemExit as e:  # usage error (2) or --help (0)
        return int(e.code or 0)
    history = " ".join(shlex.quote(a) for a in ["python", "-m", BUILDER] + argv)
    try:
        tables = build_targets(args.sources, args.grid_dir, emission=args.emission,
                               spread_type=args.spread_type, spread_scale=args.spread_scale,
                               cutoff=args.cutoff, max_snap_distance=args.max_snap_distance,
                               earth_radius=args.earth_radius,
                               connectivity=args.connectivity)
        attrs = {"mitgcm_grid_name": args.grid_name} if args.grid_name else None
        write_targets(args.output, tables, into=args.into, history=history, attrs=attrs)
    except BuildError as e:
        print("error: {0}".format(e), file=sys.stderr)
        return 1
    except (OSError, ValueError) as e:
        print("error: {0}".format(e), file=sys.stderr)
        return 2
    nsp = sum(e == "spread" for e in tables.source_emission)
    print("wrote {0}: {1} source(s) ({2} pointwise, {3} spread), {4} target(s); largest "
          "snap distance {5:.6g} km{6}".format(
              args.output, len(tables.source_id), len(tables.source_id) - nsp, nsp,
              len(tables.target_cell), float(np.max(tables.source_snap_distance)) / 1e3,
              "" if tables.connectivity is None else "; connectivity {0}{1}".format(
                  tables.connectivity, " (zonally periodic)" if tables.periodic else "")))
    if args.no_check:
        return 0
    from .check import check_files
    try:
        report = check_files(args.output, grid_dir=args.grid_dir,
                             tables_only=args.into is None)
    except OSError as e:
        print("error: {0}".format(e), file=sys.stderr)
        return 2
    print(report.format_text())
    return report.exit_code()


if __name__ == "__main__":
    sys.exit(main())
