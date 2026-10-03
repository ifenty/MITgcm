#!/usr/bin/env python3
"""Write the sparse-runoff NetCDF files of the dense runoff test cases.

One script converts every dense case, so the sparse files can be regenerated
with one command from the committed inputs::

    python3 gen_sparse.py            # all cases
    python3 gen_sparse.py clim cs32  # only these

It needs ``MITgcmutils`` (``utils/python/MITgcmutils`` of this source tree is
used when the package is not installed), ``numpy`` and ``netCDF4``.

lab_sea, ``input.rnof_<X>/`` (float32 dense files, 20 x 16 lat-lon grid). The
timing of each case is its ``runoff*`` settings in ``data.exf`` and the
calendar of ``data.cal``; the mapping to a time axis is described in
``MITgcmutils.runoff.convert``:

=======  ========================  ===========================================
case     output                    dense timing (data.exf)
=======  ========================  ===========================================
const    runoff_sparse.nc,         runoffperiod = 0.
         runoff_sparse_cells.nc
daily    runoff_sparse.nc          startdate 19790101 000000, period 86400.,
                                   runoffRepCycle = 0.
month    runoff_sparse.nc          runoffperiod = -12.
month1   runoff_sparse.nc          startdate 19781201 000000, period -1.
clim     runoff_sparse.nc          startdate 19780116 120000, period
                                   2628000., runoffRepCycle = 31536000.
yearly   runoff_sparse_1978.nc,    useExfYearlyFields, startdate 19780101
         runoff_sparse_1979.nc     000000, period 86400.
=======  ========================  ===========================================

Every case is converted with the grouping of its ``runoff_sources.txt``. In
``const`` the single record gives every group fixed shares, so its
``runoff_sparse.nc`` has the four sources ``baffin`` (3 cells), ``labrador``
(2 cells), ``greenland`` and ``newfound``; ``runoff_sparse_cells.nc`` holds
the same runoff as one source per cell. In the other cases the cells of
``baffin`` and ``labrador`` vary differently in time, so their shares are not
constant and each of their cells is its own source (``baffin_288`` ...), which
gives seven sources.

The lab_sea model writes its grid with ``pkg/mnc``, not as ``.meta/.data``
files, so the grid is computed here as the model does: cell area from
``model/src/ini_spherical_polar_grid.F`` (lines 187-188) with ``rSphere =
6371.D3``, ``delX = delY = 2`` degrees, ``xgOrigin = 280``, ``ygOrigin = 46``
(``input/data``), and the wet mask from ``input/bathy.labsea1979`` (wet where
the depth is negative). ``tests/runoff/test_convert.py`` compares both with
the model's own grid output when a run directory exists.

cs32, ``global_ocean.cs32x15/input.rnof_sparse/runoff_sparse.nc`` (float64
dense files, 192 x 32 exch2 global layout): ``core_rnof_1_cs32.bin`` of
``input.icedyn`` and ``runoff_temperature.bin`` of ``input.seaice``, with
``runoffStartTime = 1296000.``, ``runoffperiod = 2592000.`` and
``repeatPeriod = 31104000.`` (their ``data.exf``). This experiment does not
compile ``pkg/cal`` (``-cal`` in ``code/packages.conf``), so its times are
seconds of model time, and model time 0 is the reference date of the time
units. The file is written on the ``360_day`` calendar, the only calendar of
schema 1.0 in which the 360-day repeat cycle is one year. The cubed-sphere
cell areas and mask are read from the grid output of a model run
(``output_esx_input.seaice`` or ``output_esx_input.icedyn``), so this case is
skipped when no run directory exists.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
VERIF = os.path.dirname(os.path.dirname(HERE))
try:
    from MITgcmutils.runoff import convert
except ImportError:
    sys.path.insert(0, os.path.join(os.path.dirname(VERIF), "utils", "python",
                                    "MITgcmutils"))
    from MITgcmutils.runoff import convert

LAB = os.path.join(VERIF, "lab_sea")
CS32 = os.path.join(VERIF, "global_ocean.cs32x15")
NX, NY = 20, 16

#: exf timing of each lab_sea case, from its data.exf and data.cal.
LAB_CASES = {
    "const": dict(dense="runoff_const.bin", period=0.0),
    "daily": dict(dense="runoff_daily.bin", period=86400.0, startdate1=19790101,
                  startdate2=0, repeat_cycle=0.0),
    "month": dict(dense="runoff_month.bin", period=-12.0, clim_year=1979),
    "month1": dict(dense="runoff_month1.bin", period=-1.0, startdate1=19781201,
                   startdate2=0),
    "clim": dict(dense="runoff_clim.bin", period=2628000.0, startdate1=19780116,
                 startdate2=120000, repeat_cycle=31536000.0),
    "yearly": dict(dense="runoff_yearly", period=86400.0, startdate1=19780101,
                   startdate2=0, years=[1978, 1979]),
}
#: Run directories that hold the cs32 grid output, in order of preference.
CS32_GRID_DIRS = ("output_esx_input.seaice", "output_esx_input.icedyn", "run")


def lab_sea_grid():
    """``(rac, hfac, xc, yc)`` of the lab_sea grid, each ``(NY, NX)``."""
    deg2rad = 2.0 * np.pi / 360.0
    r_sphere, dxy, xg0, yg0 = 6371.0e3, 2.0, 280.0, 46.0
    lat = yg0 + dxy * np.arange(NY)
    rac1 = r_sphere * r_sphere * dxy * deg2rad * np.abs(
        np.sin((lat + dxy) * deg2rad) - np.sin(lat * deg2rad))
    rac = np.repeat(rac1[:, None], NX, axis=1)
    xc, yc = np.meshgrid(xg0 + dxy * (np.arange(NX) + 0.5),
                         yg0 + dxy * (np.arange(NY) + 0.5))
    bathy = np.fromfile(os.path.join(LAB, "input", "bathy.labsea1979"),
                        dtype=">f4").reshape(NY, NX)
    return rac, (bathy < 0.0).astype(np.float64), xc, yc


def lab_sea(case, out_dir=None):
    """Convert one lab_sea case; return the files written."""
    kw = dict(LAB_CASES[case])
    src = os.path.join(LAB, "input.rnof_" + case)
    out_dir = out_dir or src
    rac, hfac, xc, yc = lab_sea_grid()
    common = dict(rac=rac, hfac=hfac, xc=xc, yc=yc, prec=32, calendar="gregorian",
                  grid_name="lab_sea")
    dense = os.path.join(src, kw.pop("dense"))
    paths = convert.dense_to_sparse(
        dense, os.path.join(out_dir, "runoff_sparse.nc"),
        sources=os.path.join(src, "runoff_sources.txt"), split_varying=True,
        **common, **kw)
    if case == "const":
        paths += convert.dense_to_sparse(
            dense, os.path.join(out_dir, "runoff_sparse_cells.nc"), **common, **kw)
    return paths


def cs32_grid_dir():
    """The first run directory holding the cs32 grid output, or None."""
    for name in CS32_GRID_DIRS:
        d = os.path.join(CS32, name)
        if os.path.exists(os.path.join(d, "RAC.meta")) and \
                os.path.exists(os.path.join(d, "hFacC.meta")):
            return d
    return None


def cs32(out_dir=None, grid_dir=None):
    """Convert the cs32 runoff and runoff temperature; return the files written."""
    grid_dir = grid_dir or cs32_grid_dir()
    if grid_dir is None:
        print("cs32: skipped, no grid output (RAC, hFacC) in " + ", ".join(CS32_GRID_DIRS))
        return []
    out_dir = out_dir or os.path.join(CS32, "input.rnof_sparse")
    os.makedirs(out_dir, exist_ok=True)
    return convert.dense_to_sparse(
        os.path.join(CS32, "input.icedyn", "core_rnof_1_cs32.bin"),
        os.path.join(out_dir, "runoff_sparse.nc"), grid_dir=grid_dir, prec=64,
        period=2592000.0, start_time=1296000.0, repeat_cycle=31104000.0,
        calendar="360_day", grid_name="global_ocean.cs32x15",
        temperature=os.path.join(CS32, "input.seaice", "runoff_temperature.bin"))


def main(argv=None):
    cases = (argv if argv is not None else sys.argv[1:]) or list(LAB_CASES) + ["cs32"]
    for case in cases:
        if case == "cs32":
            paths = cs32()
        elif case in LAB_CASES:
            paths = lab_sea(case)
        else:
            raise SystemExit("unknown case: {0}".format(case))
        for path in paths:
            print(os.path.relpath(path, VERIF))


if __name__ == "__main__":
    main()
