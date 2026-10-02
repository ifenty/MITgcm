#!/usr/bin/env python3
"""Generate the dense exf runoff forcing of the lab_sea ``input.rnof_*`` cases.

The same file is kept in every ``input.rnof_<X>`` directory; the case is taken
from the name of the directory holding the script (or from the first argument).
Run it inside that directory::

    python3 gendata.py

It reads ``../input/bathy.labsea1979`` to confirm that every source cell is a
wet cell next to land, then writes, into its own directory:

* the runoff binary or binaries: big-endian float32 (``exf_iprec = 32``),
  shape (nrec, 16, 20), in m/s, zero away from the source cells;
* ``runoff_sources.txt``: one line per source cell.

The output is deterministic (no random numbers), so rerunning the script
reproduces the committed files byte for byte.

Rates stay below 1e-6 m/s because ``input/data.exf`` sets
``useExfCheckRange = .TRUE.`` and ``pkg/exf/exf_check_range.F`` stops the run
at the first time step when runoff exceeds 1e-6 m/s on a wet cell.
"""
import os
import sys

import numpy as np

NX, NY = 20, 16            # global grid
SNX, SNY = 10, 8           # tile size (code/SIZE.h and code/SIZE.h_mpi)
XG0, YG0, DXY = 280.0, 46.0, 2.0   # xgOrigin, ygOrigin, delX = delY (input/data)

# Source cells: (group, i, j, peak rate in m/s), with 0-based global indices.
#   baffin   : three adjacent cells on the north coast. They straddle i = 9|10,
#              which is a tile boundary in code/SIZE.h and the boundary between
#              the two processes in code/SIZE.h_mpi (nPx = 2).
#   labrador : two adjacent cells that straddle j = 7|8, a tile boundary that
#              lies inside one process.
#   greenland: one cell on the open west-Greenland coast, away from tile edges.
#   newfound : one cell on the southern, ice-free coast.
SOURCES = [
    ("baffin",    8, 14, 6.0e-7),
    ("baffin",    9, 14, 9.0e-7),
    ("baffin",   10, 14, 7.5e-7),
    ("labrador",  7,  7, 8.0e-7),
    ("labrador",  7,  8, 5.0e-7),
    ("greenland", 13, 11, 9.5e-7),
    ("newfound", 12,  2, 7.0e-7),
]

RATE_LIMIT = 1.0e-6        # exf_check_range.F upper bound for runoff, m/s

# Seasonal factor of the 12-record files, January to December. December and
# January differ by a factor of four, so a wrong wrap changes the answer.
MONTH_FACTOR = [0.22, 0.30, 0.42, 0.60, 0.85, 1.00,
                0.92, 0.78, 0.66, 0.55, 0.47, 0.88]

# Factor of the six monthly records of the non-repeating monthly file,
# December 1978 to May 1979. The values differ from MONTH_FACTOR, so reading
# this file as a repeating climatology gives a different answer.
MONTH1_FACTOR = [0.95, 0.35, 0.70, 0.50, 1.00, 0.25]

# Scale of each yearly file: the last days of 1978 and the first days of 1979
# differ by a factor of about two.
YEAR_FACTOR = {1978: 1.00, 1979: 0.45}


def daily_factor(nrec, phase):
    """Factor in [0.3, 1.0] that differs for every record and source cell."""
    n = np.arange(nrec, dtype=np.float64)
    return 0.65 + 0.35 * np.sin(0.7 * n + phase)


def month_factor(phase):
    """Twelve factors in (0, 1], with a small source-dependent modulation."""
    m = np.arange(12, dtype=np.float64)
    return np.array(MONTH_FACTOR) * (0.95 + 0.05 * np.sin(m + phase))


def month1_factor(phase):
    """Six factors in (0, 1] for the non-repeating monthly records."""
    m = np.arange(len(MONTH1_FACTOR), dtype=np.float64)
    return np.array(MONTH1_FACTOR) * (0.95 + 0.05 * np.sin(1.3 * m + phase))


def field(factors):
    """Dense (nrec, NY, NX) float32 runoff from per-source factor series."""
    nrec = len(factors[0])
    out = np.zeros((nrec, NY, NX), dtype=np.float64)
    for (_, i, j, peak), fac in zip(SOURCES, factors):
        out[:, j, i] = peak * fac
    out = out.astype(">f4")
    if out.min() < 0.0 or out.max() >= RATE_LIMIT:
        raise SystemExit("runoff outside [0, %g) m/s" % RATE_LIMIT)
    return out


def check_sources(here):
    """Stop unless each source cell is wet and has a land neighbour."""
    bathy = np.fromfile(os.path.join(here, "..", "input", "bathy.labsea1979"),
                        dtype=">f4").reshape(NY, NX)
    wet = bathy < 0.0
    for name, i, j, _ in SOURCES:
        if not wet[j, i]:
            raise SystemExit("source %s (%d,%d) is on land" % (name, i, j))
        near = [(i - 1, j), (i + 1, j), (i, j - 1), (i, j + 1)]
        if all(wet[jj, ii] for ii, jj in near
               if 0 <= ii < NX and 0 <= jj < NY):
            raise SystemExit("source %s (%d,%d) is not coastal" % (name, i, j))


def write_sources(here):
    """Write runoff_sources.txt, the table later converted to sparse form."""
    lines = [
        "# lab_sea runoff source cells (written by gendata.py)",
        "# i, j     : 0-based global cell indices (x, y) of the 20x16 grid",
        "# cell     : 0-based index in the flattened dense record, j*20 + i",
        "# lon, lat : cell centre, degrees east and north",
        "# tile     : tile number (bj-1)*2 + bi of the 2x2 layout of 10x8",
        "#            tiles, with bi and bj as in code/SIZE.h",
        "# proc     : 0-based MPI process with code/SIZE.h_mpi (nPx=2, nSy=2)",
        "# peak     : largest runoff rate the generator can give, m/s",
        "# group         i   j  cell     lon    lat  bi  bj  tile  proc"
        "       peak",
    ]
    for name, i, j, peak in SOURCES:
        bi, bj = i // SNX + 1, j // SNY + 1
        lines.append("%-12s %4d %3d %5d %7.1f %6.1f %3d %3d %5d %5d  %9.3e" % (
            name, i, j, j * NX + i, XG0 + DXY * (i + 0.5),
            YG0 + DXY * (j + 0.5), bi, bj, (bj - 1) * 2 + bi, bi - 1, peak))
    with open(os.path.join(here, "runoff_sources.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    case = sys.argv[1] if len(sys.argv) > 1 else os.path.basename(here)
    case = case.replace("input.rnof_", "")
    check_sources(here)
    phases = [0.9 * k for k in range(len(SOURCES))]
    files = {}
    if case == "const":
        # One record; exf reads record 1 once (runoffperiod = 0).
        files["runoff_const.bin"] = field([np.array([0.8])] * len(SOURCES))
    elif case == "daily":
        # 40 daily records starting 1 January 1979, 00:00.
        files["runoff_daily.bin"] = field(
            [daily_factor(40, p) for p in phases])
    elif case in ("month", "clim"):
        # 12 records, January to December.
        files["runoff_%s.bin" % case] = field(
            [month_factor(p) for p in phases])
    elif case == "month1":
        # 6 calendar-month records, December 1978 to May 1979, not repeated.
        files["runoff_month1.bin"] = field(
            [month1_factor(p) for p in phases])
    elif case == "yearly":
        # One file per year, 365 daily records each, first record at
        # 1 January 00:00. Neither 1978 nor 1979 is a leap year.
        for year, scale in sorted(YEAR_FACTOR.items()):
            files["runoff_yearly_%4d" % year] = field(
                [scale * daily_factor(365, p + 0.3 * (year - 1978))
                 for p in phases])
    else:
        raise SystemExit("unknown case: %s" % case)
    for name, data in files.items():
        data.tofile(os.path.join(here, name))
        print("%s: %d records, max %.3e m/s" % (name, data.shape[0],
                                                data.max()))
    write_sources(here)


if __name__ == "__main__":
    main()
