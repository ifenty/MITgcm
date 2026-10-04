# pkg/rnf: sparse runoff

## Purpose

`pkg/rnf` reads runoff from a NetCDF file organized by source (a river, a
glacier outlet) and hands it to `pkg/exf` as the `runoff` field. A dense
`runofffile` stores one value per surface cell per record, and almost all of
them are zero at high resolution. Here each source lists the ocean cells it
feeds and the fraction of its flux that goes to each of them, and the time
series are stored once per source.

- **Volume:** the flux in m³/s is split by fraction, divided by the cell area
  and assigned to the exf `runoff` field in m/s. Everything downstream in exf
  and in the model is the code a dense file uses.
- **Temperature, salinity and passive tracers:** optional. They enter as
  tendency terms, in the way `pkg/icefront` and `pkg/shelfice` add theirs.

The package needs `pkg/exf` (`pkg_depend`: `rnf +exf`) with `ALLOW_RUNOFF`
defined, and a build with NetCDF.

## Status

The volume flux of a file with **one constant record** works. The package
reads the file, places its target cells on the tiles of each process, checks
the file and the targets, and assigns the exf `runoff` field at every step.
A file with any other time sampling is refused, and temperature, salinity and
tracers are not applied yet.

- With `useRNF=.FALSE.` (the default) the package changes no result.
- With `useRNF=.TRUE.` and a valid file, a sparse run reproduces the dense
  `runofffile` run of the same runoff to round-off.

| File | Content | State |
|---|---|---|
| `RNF_OPTIONS.h` | CPP options | no option yet |
| `RNF_SIZE.h` | array bounds | `RNF_nSrcTile` 2000, `RNF_nTgtTile` 10000, `RNF_nBuf` 1000; a count above a bound stops the run and prints the value needed |
| `RNF.h` | parameters, per-tile lists and dense fields in common blocks | done for the volume flux |
| `rnf_readparms.F` | reads `data.rnf`; refuses a blank `RNF_file` and a build without NetCDF | done |
| `rnf_check.F` | the other configuration refusals | done |
| `rnf_summary.F` | prints the parameters, the bounds and what the read found | done |
| `rnf_nc_utils.F` | NetCDF helpers, all inside `#ifdef HAVE_NETCDF` | done for the static read and the flux record |
| `rnf_init_fixed.F` | static read of the file, placement of targets on tiles, file and target checks, fraction sums | done |
| `rnf_init_varia.F` | zeroes the fields, reads the one constant record, reports the flux | done for a constant record |
| `rnf_fields_load.F` | builds the dense fields at each step | done for a constant record; time records to come |
| `rnf_exf_runoff.F` | assigns the exf `runoff` field | done |
| `rnf_tendency_apply.F` | `RNF_TENDENCY_APPLY_T`, `_S` and `_PTR` | do nothing yet |

Time records and interpolation, temperature, salinity, passive tracers,
diagnostics, the monitor and the adjoint (TAF) list files are not there yet.

## Parameters

`data.pkg` switches the package on with `useRNF=.TRUE.`. `data.rnf` holds the
namelist `RNF_PARM01`:

| Parameter | Default | Meaning |
|---|---|---|
| `RNF_file` | `' '` | NetCDF file, or its base name with `RNF_useYearlyFiles`; a blank name is fatal |
| `RNF_holdRecord` | `.FALSE.` | false: linear interpolation in time, as exf; true: hold each record over its interval |
| `RNF_startDate1`, `RNF_startDate2` | 0 | start date override (`YYYYMMDD`, `HHMMSS`) |
| `RNF_startTime` | unset | start time override (s) |
| `RNF_period` | unset | period override: 0 constant, > 0 seconds, -12 monthly climatology, -1 monthly |
| `RNF_repCycle` | unset | repeat cycle override (s) |
| `RNF_useYearlyFiles` | `.FALSE.` | append `_YYYY` to `RNF_file` by model year |
| `RNF_useTemp` | `.TRUE.` | apply the runoff temperature if the file has it |
| `RNF_useSalt` | `.TRUE.` | apply the runoff salinity if the file has it |
| `RNF_usePtracers` | `.TRUE.` | apply the runoff tracers if the file has them |
| `RNF_monFreq` | `monitorFreq` | monitor interval (s) |
| `RNF_debugLev` | `debugLevel` | message level |

An unset timing override means the value in the file is used. Of these,
`RNF_file` and `RNF_debugLev` have an effect so far; `RNF_useYearlyFiles` and
an `RNF_period` other than 0 stop the run, because only a file with one
constant record can be read.

The constants fixed by the model contract are in `RNF.h`: `RNF_fracTol`
(1e-6) on the fraction sum of a source, `RNF_areaTol` (1e-4) between
`target_cell_area` and `rA`, `RNF_fluxMax` (1e30) above which a flux counts
as a missing value, and `RNF_maxErrMsg` (20) messages per error counter and
process.

## Configurations that are refused

Each one stops the run with an error that names the parameter.

| Refused | Where | Reason |
|---|---|---|
| blank `RNF_file` | `RNF_READPARMS` | there is nothing to read |
| build without NetCDF (`HAVE_NETCDF` undefined) | `RNF_READPARMS` | the file cannot be read |
| `useEXF=.FALSE.`, or `ALLOW_RUNOFF` undefined | `RNF_CHECK` | the package fills the exf `runoff` field |
| `runofffile` not blank | `RNF_CHECK` | dense and sparse runoff are mutually exclusive |
| `runoftempfile` not blank | `RNF_CHECK`; in a build without `ALLOW_RUNOFTEMP`, `EXF_CHECK` stops the run first, because exf cannot read that file | exf would apply the dense runoff temperature to the sparse volume |
| `runoffconst` not 0 | `RNF_CHECK` | the package would overwrite it without notice |
| `exf_outscal_sflux` not 1 | `RNF_CHECK` | it scales runoff in `EmPmR` but not in the other places that use the field |
| `SHI_update_kTopC` with `useShelfIce` | `RNF_CHECK` | a moving ice-shelf edge can cover a target cell, which then loses its volume |
| `USE_OLD_EXTERNAL_FORCING` defined | `RNF_CHECK` | the package has no hooks in the old forcing routines |

The first two are in `RNF_READPARMS` because they guard the read of the file
in `RNF_INIT_FIXED`, which the model calls before `RNF_CHECK`.

The file itself is refused in `RNF_INIT_FIXED` (and in `RNF_NC_READ_FLUX` or
`RNF_NC_ATT_REAL` for the flux record), with a message that names the source
id, or the table entry where no source can be named:

| Refused | Reason |
|---|---|
| a schema major version other than 1 | the reader knows schema 1 |
| `mitgcm_grid_nx`, `mitgcm_grid_ny` that differ from the model's global I/O layout, are missing, or are not one number | `target_cell` counts in that layout |
| `mitgcm_time_sampling` other than `constant`, more than one record, `RNF_useYearlyFiles`, `RNF_period` ≠ 0 | time records are not implemented yet (RUNOFF-005) |
| a missing dimension (`time`, `source`, `target`) or variable (`source_id`, `target_source`, `target_cell`, `target_fraction`, `runoff_flux`), or `runoff_flux` not dimensioned `(time, source)` | required input |
| `target_source` or `target_cell` out of range | an out-of-range index would be placed on a tile by the integer division, or silently dropped |
| `target_level` ≠ 1 | schema 1.0 allows only the surface cell |
| `target_fraction` outside [0, 1] | a negative fraction can hide in a sum that is still 1 |
| a source whose fractions do not sum to 1 within `RNF_fracTol` | mass would be lost; this is also how a target that no tile owns is caught (a blank exch2 tile, or a cell no facet uses) |
| a target on land, beyond an open boundary (`maskInC` = 0), under an ice shelf (`kTopC` ≠ 0), or whose `target_cell_area` differs from `rA` by more than `RNF_areaTol` | the volume would be lost, or the file was built for another grid |
| a missing flux: not a number, above `RNF_fluxMax`, or equal to the `_FillValue` or `missing_value` of `runoff_flux` | the model contract does not allow a missing flux |
| a count above `RNF_nSrcTile` or `RNF_nTgtTile` on one tile | the message prints the value the bound needs |

Every process reads the whole file, so the refusals of the header, of a table
entry and of the flux are reached by every process by itself, and the
fraction sums are the same everywhere. A refused target and a count above an
array bound are seen only by the process that owns the tile: they are
counted, the counts are summed over all processes with `GLOBAL_SUM_INT`, and
then every process stops. `ALL_PROC_DIE` ends MPI on the calling process
only, so a stop on one process alone would hang the others.

An entry that the table checks refuse is not placed, so its fraction would be
missing from the sum of its source: the fraction sums are computed only when
every entry was accepted, and the log says so.

## Tests

- `verification/lab_sea/input.rnof_sp_const`: the runoff of the dense case
  `input.rnof_const` in sparse form (4 grouped sources, 7 target cells, one
  of them spanning a tile and process boundary). It must reproduce
  `results/output.rnof_const.txt`.
- `verification/global_ocean.cs32x15/input.rnof_sp_icedyn`: the runoff of
  `input.icedyn` in sparse form on the cubed sphere (1189 one-cell sources on
  the 192 x 32 exch2 I/O layout, 12 tiles). It must reproduce
  `results/output.icedyn.txt`.

Both files are written by
`verification/lab_sea/input.rnof_const/gen_sparse.py`. The development
repository holds the refusal and placement checks.

## Design

The design, with the reasons for each choice and the source lines it relies
on, is in `docs/package_design.md` of the development repository
<https://github.com/ifenty/MITgcm_new_runoff>. The file format is in
`docs/runoff_schema.md` there.
