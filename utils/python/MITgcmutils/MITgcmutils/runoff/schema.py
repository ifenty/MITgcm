"""Constants of the sparse-runoff NetCDF schema, version 1.0.

This module is the single Python statement of the fixed vocabulary in
``docs/runoff_schema.md`` (project repository): the schema version, reserved
names, the variable table, allowed units, calendars, recommended
``source_type`` values, the ``source_id`` pattern and the numerical
tolerances. :mod:`MITgcmutils.runoff.check` enforces these values and
:mod:`MITgcmutils.runoff.example` writes a file that satisfies them. Section
numbers in comments refer to that schema document.
"""

import re

#: Version written by this package, and the major version the checker accepts
#: (schema section 12: any minor version of a known major version is readable).
SCHEMA_VERSION = "1.0"
SUPPORTED_MAJOR_VERSIONS = (1,)

#: Accepted NetCDF data models (section 1).
NETCDF_FORMATS = ("NETCDF4", "NETCDF4_CLASSIC")

# ---------------------------------------------------------------------------
# Names (sections 1-5)

#: Prefixes of schema variable names; ``time`` and ``time_bnds`` are reserved too.
RESERVED_VAR_PREFIXES = ("runoff_", "source_", "target_", "alias_")
RESERVED_VAR_NAMES = ("time", "time_bnds")
#: Prefix of model-control global attributes.
RESERVED_ATTR_PREFIX = "mitgcm_"

#: Every ``mitgcm_*`` global attribute defined by schema 1.0 (sections 4.1, 4.2).
MITGCM_ATTRS = (
    "mitgcm_runoff_schema_version",
    "mitgcm_grid_nx",
    "mitgcm_grid_ny",
    "mitgcm_time_sampling",
    "mitgcm_time_period",
    "mitgcm_time_repeat",
    "mitgcm_grid_name",
    "mitgcm_grid_description",
)

#: Dimension names (section 2).
DIM_TIME = "time"
DIM_SOURCE = "source"
DIM_TARGET = "target"
DIM_ALIAS = "alias"
DIM_NV = "nv"
REQUIRED_DIMS = (DIM_TIME, DIM_SOURCE, DIM_TARGET)
#: Maximum length of a ``source_id`` and of its character dimension.
MAX_ID_LEN = 64

#: Static-table variables: name -> (dimensions, type class, required).
#: ``"*"`` stands for the string-length dimension of a char array. Type classes
#: are ``"double"`` (float64 only), ``"float"`` (float32 or float64), ``"int"``,
#: ``"char"`` (fixed-length char array) and ``"string"`` (either a char array or
#: a variable-length string). ``required`` is True, False, or the name of a
#: dimension whose presence makes it required. S04 checks every one present.
TABLE_VARIABLES = {
    "time": (("time",), "double", True),
    "time_bnds": (("time", "nv"), "double", False),  # required if > 1 record (M04)
    "source_id": (("source", "*"), "char", True),
    "source_name": (("source",), "string", False),
    "source_type": (("source",), "string", False),
    "source_lon": (("source",), "float", False),
    "source_lat": (("source",), "float", False),
    "source_notes": (("source",), "string", False),
    "source_reference": (("source",), "string", False),
    "alias_source": (("alias",), "int", "alias"),
    "alias_name": (("alias",), "string", "alias"),
    "alias_scheme": (("alias",), "string", False),
    "target_source": (("target",), "int", True),
    "target_cell": (("target",), "int", True),
    "target_fraction": (("target",), "float", True),
    "target_level": (("target",), "int", False),
    "target_cell_area": (("target",), "float", False),
    "target_lon": (("target",), "float", False),
    "target_lat": (("target",), "float", False),
}

#: Time series (section 3.5): all have dims ``(time, source)`` and a float type.
TIMESERIES_DIMS = ("time", "source")
FLUX_VAR = "runoff_flux"
TEMPERATURE_VAR = "runoff_temperature"
SALINITY_VAR = "runoff_salinity"
PTRACER_PREFIX = "runoff_ptracer_"
TIMESERIES_VARIABLES = (FLUX_VAR, TEMPERATURE_VAR, SALINITY_VAR)
REQUIRED_TIMESERIES = (FLUX_VAR,)
#: ``<NAME>`` of ``runoff_ptracer_<NAME>``: letters, digits and ``_``.
PTRACER_NAME_RE = re.compile(r"^[A-Za-z0-9_]+$")

#: Variables the model reads (section 1): these names, and every variable whose
#: name starts with ``runoff_``. They must not be packed (S07). ``target_lon``,
#: ``target_lat`` and user-added ``target_*`` variables are not model-read.
MODEL_READ_VARIABLES = ("time", "time_bnds", "source_id", "target_source",
                        "target_cell", "target_fraction", "target_level",
                        "target_cell_area")
MODEL_READ_PREFIXES = ("runoff_",)
PACKING_ATTRS = ("scale_factor", "add_offset")

#: S08: text attributes the model reads are ASCII NC_CHAR, not NC_STRING
#: (section 1): these global attributes, these attributes of ``time`` and
#: ``time_bnds``, and ``units`` of every ``runoff_*`` variable. Descriptive
#: ``mitgcm_*`` attributes the model doesn't read (``mitgcm_grid_name``,
#: ``mitgcm_grid_description``) may be any type and any UTF-8 text.
MODEL_READ_GLOBAL_TEXT_ATTRS = ("mitgcm_runoff_schema_version", "mitgcm_time_sampling",
                                "mitgcm_time_repeat")
MODEL_READ_TIME_ATTRS = ("units", "calendar")
NC_CHAR = 2      # netCDF external type codes (netcdf.h)
NC_STRING = 12

#: S09: on a numeric model-read variable these are numbers of the variable's
#: own type (section 1); the reader reads them with ``NF_GET_ATT_DOUBLE``.
MISSING_VALUE_ATTRS = ("_FillValue", "missing_value")

#: P02: keys of ``netCDF4.Variable.filters()`` allowed on model-read variables
#: (section 8): deflate (``zlib``, with its ``complevel`` parameter), ``shuffle``
#: and ``fletcher32``. Any other filter that is set (zstd, bzip2, szip, blosc)
#: is an error, because the model's netCDF build may lack its plugin.
ALLOWED_FILTER_KEYS = ("zlib", "complevel", "shuffle", "fletcher32")

#: Index variables carry no ``units`` attribute (section 6.2).
INDEX_VARIABLES = ("target_source", "target_cell", "target_level", "alias_source")

# ---------------------------------------------------------------------------
# Units (section 6). Compared after trimming surrounding spaces, case-sensitive
# (``celsius`` and ``Celsius`` are both listed explicitly).

UNITS = {
    FLUX_VAR: ("m3 s-1", "m3/s", "m^3/s", "m3.s-1", "m^3 s^-1", "m3 s^-1"),
    TEMPERATURE_VAR: ("degC", "degree_Celsius", "degrees_Celsius", "degree_C",
                      "degrees_C", "celsius", "Celsius"),
    SALINITY_VAR: ("g kg-1", "g/kg", "1e-3", "0.001", "psu", "PSU", "PSS-78", "1"),
    "target_fraction": ("1",),
    "target_cell_area": ("m2", "m^2"),
}
LON_UNITS = ("degrees_east", "degree_east", "degree_E", "degrees_E")
LAT_UNITS = ("degrees_north", "degree_north", "degree_N", "degrees_N")
#: U01: table/coordinate variables that must carry a ``units`` attribute from
#: this list whenever they are present. User-added ``*_lon``/``*_lat``
#: variables are not checked. Time units are M01; time series D01/D07.
TABLE_UNITS = {
    "target_fraction": UNITS["target_fraction"],
    "target_cell_area": UNITS["target_cell_area"],
    "source_lon": LON_UNITS,
    "source_lat": LAT_UNITS,
    "target_lon": LON_UNITS,
    "target_lat": LAT_UNITS,
}

#: ``<days|hours|minutes|seconds> since <YYYY-MM-DD>[ hh:mm[:ss]][Z]``;
#: ``T`` is also accepted as the date/time separator.
TIME_UNITS_RE = re.compile(
    r"^(days|hours|minutes|seconds) since "
    r"(\d{4})-(\d{2})-(\d{2})"
    r"(?:[ T](\d{2}):(\d{2})(?::(\d{2}))?)?"
    r"(Z)?$"
)
#: Seconds per time unit.
TIME_UNIT_SECONDS = {"days": 86400.0, "hours": 3600.0, "minutes": 60.0, "seconds": 1.0}

# ---------------------------------------------------------------------------
# Calendars (section 6.3): CF calendar -> MITgcm cal package ``TheCalendar``.

CALENDARS = {
    "standard": "gregorian",
    "gregorian": "gregorian",
    "proleptic_gregorian": "gregorian",
    "noleap": "noLeapYear",
    "365_day": "noLeapYear",
    "360_day": "model",
}
#: Calendar attribute values are matched ignoring letter case (CF); keys are
#: lower case. Files compare equal for X01 when their calendars map to the
#: same MITgcm calendar.
#: A missing ``calendar`` attribute means this (CF), with warning M02.
DEFAULT_CALENDAR = "standard"

# ---------------------------------------------------------------------------
# Timing attributes (sections 4.1, 7).

TIME_SAMPLINGS = ("constant", "fixed", "monthly", "yearly")
TIME_REPEATS = ("none", "annual")
DEFAULT_TIME_REPEAT = "none"
#: Yearly file names: ``<anything>_YYYY.nc`` (section 7).
YEARLY_FILE_RE = re.compile(r"_(\d{4})\.nc$")

# ---------------------------------------------------------------------------
# Sources (section 3.2).

SOURCE_TYPES = ("river", "glacier", "ice_sheet_basin", "iceberg_melt",
                "groundwater", "other")
#: ASCII letters, digits, ``_``, ``-``, ``.``; starts with a letter or digit;
#: 1 to 64 characters.
SOURCE_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.\-]{0,63}$")

# ---------------------------------------------------------------------------
# Tolerances and ranges (sections 3, 9).

#: T05: each source's fractions, summed in float64, satisfy |sum - 1| <= this.
FRACTION_SUM_TOL = 1e-6
#: R02 (and the model): relative tolerance of ``target_cell_area`` against ``rA``.
AREA_RTOL = 1e-4
#: R03: tolerance of ``target_lon``/``target_lat`` against ``XC``/``YC``, degrees.
LONLAT_TOL = 1e-3
#: D04: warning range of ``runoff_temperature``, degrees Celsius.
TEMPERATURE_RANGE = (-2.5, 40.0)
#: D06: warning threshold of ``runoff_salinity``.
SALINITY_MAX = 45.0
#: The single time-equality tolerance (section 7): contiguous bounds (M04),
#: fixed spacing, month/year edges, midpoints and annual coverage (M05),
#: yearly-file limits (M06), file-to-file continuity, fixed spacing across
#: files and the common offset of each yearly file's first time value (X01).
#: It absorbs only floating-point
#: representation error; a bound within it of a month or year edge counts as
#: that edge.
TIME_EQUAL_TOL_SECONDS = 1e-3

# ---------------------------------------------------------------------------
# Discovery attributes (section 4.3), plus the CF ``Conventions`` attribute
# (section 1). Missing ones are reported as information (S06).

CONVENTIONS = "CF-1.11, ACDD-1.3"
RECOMMENDED_GLOBAL_ATTRS = (
    "Conventions",
    "title", "summary", "institution", "source", "history", "references",
    "comment", "creator_name", "creator_email", "creator_url",
    "contributor_name", "contributor_role", "project", "license",
    "date_created", "date_modified", "product_version", "keywords",
    "time_coverage_start", "time_coverage_end",
    "geospatial_lat_min", "geospatial_lat_max",
    "geospatial_lon_min", "geospatial_lon_max",
)

#: Tables-only checking (``check_files(..., tables_only=True)``, section 9): a
#: file holding only the source, alias and target tables, such as the output of
#: :mod:`MITgcmutils.runoff.targets`, has no time axis or time series. These
#: rules, or parts of rules, are then not checked, and ``S10`` (I) says so.
TABLES_ONLY_SKIPPED = (
    "S03 (the time dimension requirement)",
    "S04 (the time and runoff_flux requirements)",
    "M01-M06", "D01-D09", "P01",
)
#: Also skipped when several files are checked together in tables-only mode.
TABLES_ONLY_SKIPPED_MULTI = ("X01 (time order and continuity across files)",)

#: Rule id -> level, for every rule of schema section 9.
RULES = {
    "S01": "E", "S02": "E", "S03": "E", "S04": "E", "S05": "E", "S06": "I",
    "S07": "E", "S08": "E", "S09": "E", "S10": "I",
    "G01": "E",
    "I01": "E", "I02": "E", "I03": "W", "I04": "W",
    "A01": "E", "A02": "W",
    "T01": "E", "T02": "E", "T03": "E", "T04": "E", "T05": "E", "T06": "W",
    "T07": "E", "T08": "I",
    "M01": "E", "M02": "W", "M03": "E", "M04": "E", "M05": "E", "M06": "E",
    "D01": "E", "D02": "E", "D03": "W", "D04": "W", "D05": "E", "D06": "W",
    "D07": "E", "D08": "W", "D09": "E",
    "U01": "E", "P01": "W", "P02": "E", "X01": "E",
    "R01": "E", "R02": "E", "R03": "W",
}
