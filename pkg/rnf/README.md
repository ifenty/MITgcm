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

The volume flux works, in every time mode, and so do the temperature, the
salinity and the passive tracers. The package
reads the file, places its target cells on the tiles of each process, checks
the file and the targets, selects the time records that bracket the model
time, assigns the exf `runoff` field at every step and adds the heat, the salt
and the tracer of the runoff as tendency terms at the target cells.
What is not there yet: the diagnostics of the input fields and the monitor
(`RNF_MONITOR`), the budget checks over time and over the domain, the adjoint
(TAF) list files, and the `addMass` path that interior and under-shelf
targets need.

**Time modes.** `constant` (one record), a fixed period with or without a
repeat cycle, a monthly climatology (12 records, January to December,
repeated every model year), consecutive calendar months, and `_YYYY` yearly
files of a fixed-period series. Record selection is delegated to the pkg/exf
routine of each mode (`EXF_GetFFieldRec`, `cal_GetMonthsRec`,
`EXF_GetMonthsRec`), so a sparse run picks records with the code a dense
`runofffile` run runs. `RNF_holdRecord` is the package's own addition: it
holds the record whose interval contains the model time instead of
interpolating. `yearly` *sampling* (one record per calendar year) is refused,
because exf has no such mode and schema 1.0 puts its time at the midpoint of
the year.

- With `useRNF=.FALSE.` (the default) the package changes no result.
- With `useRNF=.TRUE.` and a valid file, a sparse run reproduces the dense
  `runofffile` run of the same runoff to round-off.

| File | Content | State |
|---|---|---|
| `RNF_OPTIONS.h` | CPP options | no option yet |
| `RNF_SIZE.h` | array bounds | `RNF_nSrcTile` 2000, `RNF_nTgtTile` 10000, `RNF_nBuf` 1000, `RNF_nTr` 5; a count above a bound stops the run and prints the value needed |
| `RNF.h` | parameters, per-tile lists, record buffers and dense fields in common blocks | done |
| `rnf_readparms.F` | reads `data.rnf`; refuses a blank `RNF_file` and a build without NetCDF | done |
| `rnf_check.F` | the other configuration refusals | done |
| `rnf_summary.F` | prints the parameters, the bounds and what the read found | done |
| `rnf_nc_utils.F` | NetCDF helpers, all inside `#ifdef HAVE_NETCDF`; `RNF_NC_SERIES` finds the optional time series and matches the tracer names, `RNF_NC_READ_FLUX` reads one record of every series and `RNF_NC_READ_ONE` one series of it | done |
| `rnf_init_fixed.F` | static read of the file, placement of targets on tiles, file and target checks, fraction sums | done |
| `rnf_time_setup.F` | resolves the file's time axis and the `data.rnf` overrides into exf's period, start time and repeat cycle | done |
| `rnf_getrec.F` | `RNF_GETREC` (the two records bracketing the model time and the weight, from the pkg/exf routine of the mode, plus hold-exact) and `RNF_FILE_NAME` (`<base>_YYYY.nc`) | done |
| `rnf_init_varia.F` | zeroes the fields, empties the record buffers, loads the records of the start time and reports the flux | done |
| `rnf_fields_load.F` | `RNF_FIELDS_LOAD` decides the two time levels; `RNF_LOAD_AT` selects the records, reads the ones no buffer holds, interpolates or holds every series and builds the dense fields, tracing the selection at `RNF_debugLev` ≥ 3; `RNF_COPY_APPLY`/`RNF_ZERO_APPLY` keep the set the tendency terms read | done |
| `rnf_exf_runoff.F` | assigns the exf `runoff` field | done |
| `rnf_tendency_apply.F` | `RNF_TENDENCY_APPLY_T`, `_S` and `_PTR`: the heat, salt and tracer terms at the target level | done |
| `rnf_diagnostics_init.F` | registers `RNFgT`, `RNFgS` and `RNFtrNN`, the diagnostics of those terms | done; the input-only diagnostics and the monitor are not |

The monitor (`RNF_MONITOR`), the diagnostics of the input fields and the
adjoint (TAF) list files are not there yet.

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
| `RNF_useTemp` | `.TRUE.` | apply `runoff_temperature` if the file has it; false treats it as absent |
| `RNF_useSalt` | `.TRUE.` | apply `runoff_salinity` if the file has it; false treats it as absent |
| `RNF_usePtracers` | `.TRUE.` | apply `runoff_ptracer_*` if the file has them; false treats them as absent, and their names are then not matched to the ptracers |
| `RNF_monFreq` | `monitorFreq` | monitor interval (s) |
| `RNF_debugLev` | `debugLevel` | message level |

An unset timing override means the value in the file is used; what the two
together resolved to is reported by `RNF_SUMMARY` as `RNF_recPeriod`,
`RNF_recStart`, `RNF_recCycle` and `RNF_recDate1`, `RNF_recDate2`.
`RNF_useYearlyFiles` is the one setting a file cannot carry itself; with it,
`RNF_file` is the base name of a `<base>_YYYY.nc` set and a trailing `.nc` of
the name given is replaced, so `runoff.nc` and `runoff` both select
`runoff_1979.nc`. `RNF_startTime` is for a run without `pkg/cal`, where the
file's time is model time in seconds from the reference date of its units;
with `pkg/cal` it is refused, and the start date comes from the file or from
`RNF_startDate1`/`RNF_startDate2`, as it does for a dense exf field.
`RNF_monFreq` waits for the monitor.

**What the three switches do.** `RNF_INIT_FIXED` reports which of the three
optional series the file carries and which ptracer each runoff tracer feeds,
and the summary prints `RNF_hasTemp`, `RNF_hasSalt`, `RNF_nTrUse`,
`RNF_applyT`, `RNF_applyS` and one line per tracer:

```
RNF_SUMMARY: runoff tracer  1 is runoff_ptracer_dye, applied to ptracer  1
```

A switch set to false makes the reader treat those variables as absent and say
so per variable, which is the way to read a file whose tracers this run does
not carry. `RNF_applyT` and `RNF_applyS` say whether each term can be
non-zero at all: without `runoff_temperature` the temperature term is
identically zero, and without `runoff_salinity` the salinity term is zero as
well as long as `salt_EvPrRn` is 0, which is its default. A file that carries
only a flux therefore leaves a run bit-for-bit where it was.

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

The file itself is refused in `RNF_INIT_FIXED` (and in `RNF_NC_SERIES`,
`RNF_NC_READ_FLUX`, `RNF_NC_READ_ONE` or `RNF_NC_ATT_REAL` for the series and
their records), with a message that names the source
id, or the table entry where no source can be named:

| Refused | Reason |
|---|---|
| a schema major version other than 1 | the reader knows schema 1 |
| `mitgcm_grid_nx`, `mitgcm_grid_ny` that differ from the model's global I/O layout, are missing, or are not one number | `target_cell` counts in that layout |
| `mitgcm_time_sampling` other than `constant`, `fixed` or `monthly` | `yearly` sampling has no pkg/exf mode and schema 1.0 puts its time at the midpoint of the year |
| `mitgcm_time_repeat` other than `none` or `annual` | schema 1.0 defines only those two |
| `fixed` sampling without a positive `mitgcm_time_period`, or `fixed` + `annual` without `time_bnds` | the period and the repeat cycle (the span of the bounds) come from them |
| a repeat cycle that is not the record count times the period | the wrap would select a record the file has not got, or never reach the last records |
| a monthly climatology that has not 12 records, or whose record 1 is not in January | exf's period −12 is calendar month k in record k, January to December |
| a `time:units` that is not `<days\|hours\|minutes\|seconds> since <date>`, or a reference date that is not `YYYY-MM-DD[ hh:mm[:ss]]` | the reader converts that date through `pkg/cal` |
| a reference date before 1583 when `pkg/cal` is used | `pkg/cal` counts from 15 October 1582; give the date of record 1 as `RNF_startDate1`/`RNF_startDate2` instead. A `constant` file carries `0001-01-01` and its time axis is never read |
| a `calendar` that `pkg/cal` has not got, or that is not this run's | the record times were computed on another calendar |
| `RNF_startTime` with `pkg/cal`, or `RNF_startDate*` without it, or a monthly period without it, or `RNF_useYearlyFiles` without it | the same rules a dense exf field follows |
| a record outside the time dimension of the file being read | the series does not cover the model time: it needs more records, a repeat cycle, or another period |
| a missing dimension (`time`, `source`, `target`) or variable (`source_id`, `target_source`, `target_cell`, `target_fraction`, `runoff_flux`), or `runoff_flux` not dimensioned `(time, source)` | required input |
| `target_source` or `target_cell` out of range | an out-of-range index would be placed on a tile by the integer division, or silently dropped |
| `target_level` ≠ 1 | schema 1.0 allows only the surface cell |
| `target_fraction` outside [0, 1] | a negative fraction can hide in a sum that is still 1 |
| a source whose fractions do not sum to 1 within `RNF_fracTol` | mass would be lost; this is also how a target that no tile owns is caught (a blank exch2 tile, or a cell no facet uses) |
| a target on land, beyond an open boundary (`maskInC` = 0), under an ice shelf (`kTopC` ≠ 0), or whose `target_cell_area` differs from `rA` by more than `RNF_areaTol` | the volume would be lost, or the file was built for another grid |
| a missing value of `runoff_flux`, of `runoff_salinity` or of a tracer: not a number, above `RNF_fluxMax`, or equal to the `_FillValue` or `missing_value` of that variable | the model contract allows no missing flux, and a missing salinity or tracer has no defined meaning. A missing `runoff_temperature` **is** allowed: that source enters at the reference temperature, as a source without the variable does |
| a `runoff_ptracer_<NAME>` whose `<NAME>` matches no `PTRACERS_names` entry, or more than one, or any such variable in a run that does not use pkg/ptracers | there is no tendency array to add it to, so the water would arrive without the tracer and nothing would say so. `RNF_usePtracers=.FALSE.` is the way to read such a file on purpose |
| more `runoff_ptracer_*` variables than `RNF_nTr`, or a name that is empty or longer than `RNF_idLen` | the message prints the number the file has |
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
repository holds the refusal and placement checks, the cell-by-cell
applied-field and timing checks, and the two checks of the tendency terms:
one against the analytic value of every case of the design's reference
tables, and one against the heat that `pkg/exf` applies through its own
`runoftempfile` path, cell by cell.

## Design

The design, with the reasons for each choice and the source lines it relies
on, is in `docs/package_design.md` of the development repository
<https://github.com/ifenty/MITgcm_new_runoff>. The file format is in
`docs/runoff_schema.md` there.
