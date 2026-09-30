"""Write a small example sparse-runoff file (schema 1.0).

:func:`write_example` writes 3 sources, 5 targets, 2 aliases and 4 daily
records for the 20 x 16 lab_sea verification layout (2-degree cells from
280E, 46N). It includes every optional variable, one passive tracer
(``runoff_ptracer_dye``), the full ACDD/CF discovery attributes, one-record
chunks and zlib compression, and passes :func:`MITgcmutils.runoff.check_files`
with no errors and no warnings. The glacier source's temperature is the fill
value, which means "enters at the surface water temperature".

Keyword overrides replace parts of the file, so tests and users can build
variants (see :func:`write_example`). Cell areas and centers come from
:func:`lab_sea_grid`, a spherical 2-degree grid with the lab_sea layout; its
areas are computed analytically, not read from lab_sea model output.
"""

import numpy as np

from . import schema as S

__all__ = ["write_example", "lab_sea_grid", "NX", "NY"]

#: lab_sea global layout.
NX, NY = 20, 16
X0, Y0, DXY = 280.0, 46.0, 2.0
#: Earth radius used for the example cell areas (MITgcm default ``rSphere``).
R_SPHERE = 6370.0e3

_AUTO = "auto"


def lab_sea_grid():
    """``(XC, YC, RAC)`` of the example grid, each shaped ``(NY, NX)``.

    Flattening in C order of ``(NY, NX)`` gives the schema's global cell index
    ``cell = i + NX * j``.
    """
    xc1 = X0 + DXY * (np.arange(NX) + 0.5)
    yc1 = Y0 + DXY * (np.arange(NY) + 0.5)
    xc, yc = np.meshgrid(xc1, yc1)
    lat_s = np.deg2rad(Y0 + DXY * np.arange(NY))
    lat_n = lat_s + np.deg2rad(DXY)
    rac1 = R_SPHERE ** 2 * np.deg2rad(DXY) * (np.sin(lat_n) - np.sin(lat_s))
    rac = np.repeat(rac1[:, None], NX, axis=1)
    return xc, yc, rac


def _default_spec():
    return {
        "format": "NETCDF4",
        "nx": NX,
        "ny": NY,
        "source_id": ["churchill_labrador", "koksoak", "jakobshavn"],
        "source_name": ["Churchill River (Labrador)", "Koksoak River",
                        "Jakobshavn Isbræ"],
        "source_type": ["river", "river", "glacier"],
        "source_lon": [299.6, 291.7, 309.9],
        "source_lat": [53.3, 58.5, 69.2],
        "source_notes": ["Synthetic example values.\nEnters Lake Melville.",
                         "Synthetic example values.",
                         "Synthetic example values; temperature is the fill "
                         "value (surface temperature)."],
        "source_reference": ["https://mitgcm.org", "https://mitgcm.org",
                             "https://mitgcm.org"],
        "alias_source": [0, 2],
        "alias_name": ["Grand River", "Sermeq Kujalleq"],
        "alias_scheme": ["historical", "Greenlandic"],
        # cell = i + NX*j with i = floor((lon-280)/2), j = floor((lat-46)/2)
        "target_source": [0, 0, 1, 2, 2],
        "target_cell": [69, 70, 125, 233, 234],
        "target_fraction": [0.7, 0.3, 1.0, 0.5, 0.5],
        "target_level": [1, 1, 1, 1, 1],
        "target_cell_area": _AUTO,
        "target_lon": _AUTO,
        "target_lat": _AUTO,
        "time_units": "days since 2000-01-01 00:00:00",
        "calendar": "standard",
        "time": [0.5, 1.5, 2.5, 3.5],
        "time_bnds": _AUTO,
        "runoff_flux": _AUTO,
        "runoff_temperature": _AUTO,
        "runoff_salinity": _AUTO,
        "ptracers": {"dye": (_AUTO, "mol m-3")},
        "time_sampling": "fixed",
        "time_period": 86400.0,
        "time_repeat": "none",
        "id_strlen": None,         # None: longest id
        "unlimited_time": True,
        "chunk_time": 1,           # records per chunk along time
        "zlib": True,
        "dtypes": {},              # variable -> numpy dtype string
        "fill_values": {},         # variable -> explicit _FillValue
        "var_attrs": {},           # variable -> {attr: value or None (delete)}
        "global_attrs": {},        # attr -> value or None (delete)
        "extra_vars": {},          # name -> (dims, values, attrs)
        "var_create": {},          # variable -> extra createVariable keywords
    }


_DEFAULT_DTYPES = {
    "time": "f8", "time_bnds": "f8",
    "source_lon": "f8", "source_lat": "f8",
    "alias_source": "i4",
    "target_source": "i4", "target_cell": "i4", "target_fraction": "f8",
    "target_level": "i4", "target_cell_area": "f8",
    "target_lon": "f8", "target_lat": "f8",
    "runoff_flux": "f4", "runoff_temperature": "f4", "runoff_salinity": "f4",
}


def _auto_bounds(time):
    t = np.asarray(time, dtype=np.float64)
    if t.size == 1:
        return None
    mid = 0.5 * (t[1:] + t[:-1])
    edges = np.concatenate([[t[0] - (mid[0] - t[0])], mid, [t[-1] + (t[-1] - mid[-1])]])
    return np.stack([edges[:-1], edges[1:]], axis=1)


def _auto_series(spec, nt, ns):
    """Default time series for ``nt`` records of ``ns`` sources."""
    base = np.array([1500.0, 800.0, 950.0]) if ns == 3 else 100.0 * (np.arange(ns) + 1)
    phase = np.sin(2 * np.pi * np.arange(nt) / max(nt, 2))
    flux = base[None, :] * (1.0 + 0.05 * phase[:, None])
    temp = np.tile(1.0 + 0.5 * np.arange(ns), (nt, 1))
    types = spec.get("source_type") or []
    glacier = [k for k, v in enumerate(types) if v == "glacier" and k < ns]
    salt = np.zeros((nt, ns))
    dye = np.zeros((nt, ns))
    dye[:, 0] = 1.0
    return flux, temp, glacier, salt, dye


def _chars(byte_strings, strlen):
    """``(n, strlen)`` ``S1`` array of NUL-padded byte strings.

    Built with numpy rather than ``netCDF4.stringtochar``, which fails on
    ``bytes_`` input with some netCDF4/numpy builds.
    """
    n = len(byte_strings)
    return np.array(byte_strings, dtype="S{0}".format(strlen)).view("S1").reshape(
        n, strlen)


def _merge(base, over):
    out = dict(base)
    for k, v in over.items():
        if v is None:
            out.pop(k, None)
        else:
            out[k] = v
    return out


def write_example(path, **overrides):
    """Write the example runoff file to ``path`` and return ``path``.

    Every key of the default specification can be overridden by keyword:
    table columns (``source_id``, ``target_cell``, ...), ``time``,
    ``time_bnds``, ``time_units``, ``calendar``, ``runoff_flux`` etc. (arrays
    shaped ``(time, source)``), ``ptracers`` (``{NAME: (values or "auto",
    units)}``), ``nx``/``ny``, ``time_sampling``/``time_period``/``time_repeat``,
    ``format``, ``id_strlen``, ``unlimited_time``, ``chunk_time``, ``zlib``.
    A value of ``None`` omits that variable or attribute. ``"auto"`` values are
    derived: time bounds from ``time`` (none for one record), cell area and
    center from :func:`lab_sea_grid`, time series from a smooth pattern.
    Dictionaries ``dtypes``, ``fill_values``, ``var_attrs`` and
    ``global_attrs`` are merged into the defaults (a ``None`` value deletes an
    attribute), and ``extra_vars`` adds ``{name: (dims, values, attrs)}``.
    ``var_create`` passes extra ``netCDF4.Dataset.createVariable`` keywords per
    variable, e.g. ``{"runoff_flux": {"compression": "zstd"}}`` (an explicit
    ``compression`` replaces the default zlib).
    """
    import netCDF4
    import cftime

    spec = _default_spec()
    unknown = set(overrides) - set(spec)
    if unknown:
        raise TypeError("unknown write_example override(s): {0}".format(sorted(unknown)))
    spec.update(overrides)
    dtypes = dict(_DEFAULT_DTYPES, **spec["dtypes"])

    fmt = spec["format"]
    netcdf4 = fmt.startswith("NETCDF4")
    vlen_ok = fmt == "NETCDF4"
    ids = list(spec["source_id"])
    ns = len(ids)
    time = None if spec["time"] is None else np.asarray(spec["time"], dtype=np.float64)
    nt = 0 if time is None else time.size

    # --- derived values ----------------------------------------------------
    bnds = spec["time_bnds"]
    if isinstance(bnds, str) and bnds == _AUTO:
        bnds = None if time is None else _auto_bounds(time)
    flux, temp, glacier, salt, dye = _auto_series(spec, nt, ns)
    series = {}
    for name, auto in ((S.FLUX_VAR, flux), (S.TEMPERATURE_VAR, temp),
                       (S.SALINITY_VAR, salt)):
        v = spec[name]
        if isinstance(v, str) and v == _AUTO:
            v = auto
            if name == S.TEMPERATURE_VAR and glacier:
                fill = netCDF4.default_fillvals[np.dtype(dtypes[name]).str[1:]]
                spec["fill_values"] = dict(spec["fill_values"])
                spec["fill_values"].setdefault(name, fill)
                v = np.array(v, dtype=np.float64)
                v[:, glacier] = spec["fill_values"][name]
        series[name] = v
    tracers = {}
    for tname, (vals, units) in (spec["ptracers"] or {}).items():
        if isinstance(vals, str) and vals == _AUTO:
            vals = dye
        tracers[S.PTRACER_PREFIX + tname] = (vals, units)

    xc, yc, rac = lab_sea_grid()
    cells = None if spec["target_cell"] is None else np.asarray(spec["target_cell"])

    def from_grid(key, field):
        v = spec[key]
        if isinstance(v, str) and v == _AUTO:
            if cells is None:
                return None
            flat = field.ravel()
            ok = (cells >= 0) & (cells < flat.size)
            return np.where(ok, flat[np.clip(cells, 0, flat.size - 1)], 0.0)
        return v

    derived = {"target_cell_area": from_grid("target_cell_area", rac),
               "target_lon": from_grid("target_lon", xc),
               "target_lat": from_grid("target_lat", yc)}

    # --- file ---------------------------------------------------------------
    ds = netCDF4.Dataset(path, "w", format=fmt)
    try:
        # Text attributes as NC_CHAR (the Fortran reader can't read NC_STRING,
        # rule S08). This is netCDF4's default; set it explicitly.
        ds.set_ncstring_attrs(False)
        # global attributes: model-read first, then discovery
        g = {"Conventions": S.CONVENTIONS,
             "mitgcm_runoff_schema_version": S.SCHEMA_VERSION,
             "mitgcm_grid_nx": np.int32(spec["nx"]) if spec["nx"] is not None else None,
             "mitgcm_grid_ny": np.int32(spec["ny"]) if spec["ny"] is not None else None,
             "mitgcm_grid_name": "lab_sea",
             "mitgcm_grid_description": "lat-lon 2-degree grid, 20 x 16 cells from "
                                        "280E, 46N; cell = i + 20*j (0-based)",
             "mitgcm_time_sampling": spec["time_sampling"],
             "mitgcm_time_period": (None if spec["time_period"] is None
                                    else np.float64(spec["time_period"])),
             "mitgcm_time_repeat": spec["time_repeat"],
             "title": "Example sparse runoff file for the lab_sea layout",
             "summary": "Three synthetic sources (two rivers, one glacier) feeding "
                        "five ocean cells, four daily records.",
             "institution": "MITgcm",
             "source": "MITgcmutils.runoff.example.write_example",
             "history": "created by MITgcmutils.runoff.example.write_example",
             "references": "docs/runoff_schema.md (sparse runoff schema 1.0)",
             "comment": "Synthetic values for testing and documentation only.",
             "creator_name": "MITgcm developers",
             "creator_email": "mitgcm-support@mitgcm.org",
             "creator_url": "https://mitgcm.org",
             "contributor_name": "MITgcm developers",
             "contributor_role": "author",
             "project": "MITgcm sparse runoff",
             "license": "MIT",
             "date_created": "2026-09-29T00:00:00Z",
             "date_modified": "2026-09-29T00:00:00Z",
             "product_version": "1.0",
             "keywords": "runoff, river discharge, glacier discharge, MITgcm"}
        # coverage attributes, when computable
        try:
            cov = bnds if bnds is not None else time
            if cov is not None and spec["time_units"] and spec["calendar"]:
                lo = cftime.num2date(float(np.min(cov)), spec["time_units"],
                                     calendar=spec["calendar"])
                hi = cftime.num2date(float(np.max(cov)), spec["time_units"],
                                     calendar=spec["calendar"])
                g["time_coverage_start"] = lo.strftime("%Y-%m-%dT%H:%M:%SZ")
                g["time_coverage_end"] = hi.strftime("%Y-%m-%dT%H:%M:%SZ")
        except Exception:  # noqa: BLE001 - invalid override: leave coverage out
            pass
        lat = derived["target_lat"]
        lon = derived["target_lon"]
        if lat is not None and lon is not None and len(lat):
            g["geospatial_lat_min"] = float(np.min(lat))
            g["geospatial_lat_max"] = float(np.max(lat))
            g["geospatial_lon_min"] = float(np.min(lon))
            g["geospatial_lon_max"] = float(np.max(lon))
        g = _merge(g, spec["global_attrs"])
        for k, v in g.items():
            if v is None:
                continue
            if isinstance(v, bool):
                v = int(v)
            if isinstance(v, int):
                v = np.int32(v)
            ds.setncattr(k, v)

        # dimensions
        ds.createDimension("time", None if spec["unlimited_time"] else nt)
        ds.createDimension("source", ns)
        ntarget = 0 if spec["target_source"] is None else len(spec["target_source"])
        ds.createDimension("target", ntarget)
        aliases = spec["alias_name"]
        if aliases is not None and len(aliases):
            ds.createDimension("alias", len(aliases))
        if bnds is not None:
            ds.createDimension("nv", 2)
        id_bytes = [s.encode("utf-8") for s in ids]
        strlen = spec["id_strlen"] or max([len(b) for b in id_bytes] + [1])
        ds.createDimension("id_strlen", strlen)

        def attrs_for(name, default):
            return _merge(default, spec["var_attrs"].get(name, {}))

        def put(name, dims, values, attrs, dtype=None, series_var=False):
            dtype = dtype or dtypes.get(name, "f8")
            kw = {}
            if name in spec["fill_values"]:
                kw["fill_value"] = np.asarray(spec["fill_values"][name]).astype(dtype)
            if series_var and netcdf4:
                kw.update(chunksizes=(max(1, spec["chunk_time"]), max(1, ns)))
                if spec["zlib"]:
                    kw.update(zlib=True, complevel=2, shuffle=True)
            extra = spec["var_create"].get(name, {})
            if "compression" in extra:
                kw.pop("zlib", None)       # an explicit compression replaces zlib
            kw.update(extra)
            v = ds.createVariable(name, dtype, dims, **kw)
            for k, a in attrs_for(name, attrs).items():
                if a is not None:
                    v.setncattr(k, a)
            if values is not None and np.size(values):
                v[:] = np.asarray(values).astype(dtype)
            return v

        def put_strings(name, dim, values, attrs):
            if values is None:
                return
            if vlen_ok:
                v = ds.createVariable(name, str, (dim,))
                v[:] = np.array(list(values), dtype=object)
            else:
                enc = [s.encode("utf-8") for s in values]
                n = max([len(b) for b in enc] + [1])
                sdim = name + "_strlen"
                ds.createDimension(sdim, n)
                v = ds.createVariable(name, "S1", (dim, sdim))
                v.set_auto_chartostring(False)
                v[:] = _chars(enc, n)
            for k, a in attrs_for(name, attrs).items():
                if a is not None:
                    v.setncattr(k, a)

        # time
        if time is not None:
            tattrs = {"long_name": "time", "standard_name": "time", "axis": "T",
                      "units": spec["time_units"], "calendar": spec["calendar"],
                      "bounds": "time_bnds" if bnds is not None else None}
            put("time", ("time",), time, tattrs)
        if bnds is not None:
            put("time_bnds", ("time", "nv"), bnds,
                {"long_name": "start and end of the interval each record covers"})

        # source table
        if ids is not None:
            v = ds.createVariable("source_id", "S1", ("source", "id_strlen"))
            v.set_auto_chartostring(False)
            v[:] = _chars(id_bytes, strlen)
            for k, a in attrs_for("source_id", {
                    "long_name": "source identifier", "cf_role": "timeseries_id"}).items():
                if a is not None:
                    v.setncattr(k, a)
        put_strings("source_name", "source", spec["source_name"],
                    {"long_name": "primary name of the source"})
        put_strings("source_type", "source", spec["source_type"],
                    {"long_name": "kind of source",
                     "comment": "river, glacier, ice_sheet_basin, iceberg_melt, "
                                "groundwater or other"})
        if spec["source_lon"] is not None:
            put("source_lon", ("source",), spec["source_lon"],
                {"long_name": "longitude of the mouth or terminus",
                 "standard_name": "longitude", "units": "degrees_east"})
        if spec["source_lat"] is not None:
            put("source_lat", ("source",), spec["source_lat"],
                {"long_name": "latitude of the mouth or terminus",
                 "standard_name": "latitude", "units": "degrees_north"})
        put_strings("source_notes", "source", spec["source_notes"],
                    {"long_name": "notes on provenance and processing"})
        put_strings("source_reference", "source", spec["source_reference"],
                    {"long_name": "citation, DOI or URL of the source data"})

        # alias table
        if aliases is not None and len(aliases):
            if spec["alias_source"] is not None:
                put("alias_source", ("alias",), spec["alias_source"],
                    {"long_name": "index of the source this alias names",
                     "instance_dimension": "source"})
            put_strings("alias_name", "alias", aliases,
                        {"long_name": "alternative name or catalogue id"})
            put_strings("alias_scheme", "alias", spec["alias_scheme"],
                        {"long_name": "naming system or authority of the alias"})

        # target table
        if spec["target_source"] is not None:
            put("target_source", ("target",), spec["target_source"],
                {"long_name": "index of the source feeding this target",
                 "instance_dimension": "source"})
        if cells is not None:
            put("target_cell", ("target",), cells,
                {"long_name": "0-based global cell index",
                 "comment": "cell = i + mitgcm_grid_nx * j, 0-based (i, j) in the "
                            "global 2D layout of a dense runoffFile"})
        if spec["target_fraction"] is not None:
            put("target_fraction", ("target",), spec["target_fraction"],
                {"long_name": "share of the source flux sent to this cell",
                 "units": "1"})
        if spec["target_level"] is not None:
            put("target_level", ("target",), spec["target_level"],
                {"long_name": "1-based model level k (schema 1.0: always 1)"})
        if derived["target_cell_area"] is not None:
            put("target_cell_area", ("target",), derived["target_cell_area"],
                {"long_name": "horizontal cell area rA", "units": "m2"})
        if derived["target_lon"] is not None:
            put("target_lon", ("target",), derived["target_lon"],
                {"long_name": "cell-center longitude", "standard_name": "longitude",
                 "units": "degrees_east"})
        if derived["target_lat"] is not None:
            put("target_lat", ("target",), derived["target_lat"],
                {"long_name": "cell-center latitude", "standard_name": "latitude",
                 "units": "degrees_north"})

        # time series
        ts_attrs = {
            S.FLUX_VAR: {"long_name": "runoff volume flux of the source",
                         "units": "m3 s-1",
                         "comment": "split among targets by target_fraction"},
            S.TEMPERATURE_VAR: {"long_name": "runoff temperature", "units": "degC",
                                "comment": "fill value: enters at the surface water "
                                           "temperature"},
            S.SALINITY_VAR: {"long_name": "runoff salinity", "units": "g kg-1"},
        }
        for name in S.TIMESERIES_VARIABLES:
            if series[name] is not None:
                put(name, ("time", "source"), series[name], ts_attrs[name],
                    series_var=True)
        for name, (vals, units) in tracers.items():
            put(name, ("time", "source"), vals,
                {"long_name": "concentration of ptracer {0} in runoff".format(
                    name[len(S.PTRACER_PREFIX):]),
                 "units": units}, dtype=dtypes.get(name, "f4"), series_var=True)

        # additions
        for name, (dims, values, attrs) in spec["extra_vars"].items():
            arr = np.asarray(values)
            dtype = dtypes.get(name, arr.dtype.str[1:] if arr.dtype.kind in "fiu" else "f8")
            put(name, tuple(dims), arr, attrs or {}, dtype=dtype,
                series_var=tuple(dims) == S.TIMESERIES_DIMS)
    finally:
        ds.close()
    return path
