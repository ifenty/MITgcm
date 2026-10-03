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

This is the package skeleton. It compiles, registers with the model, reads
`data.rnf`, prints its parameters and checks the configuration. It does not
read the runoff file yet and applies no runoff.

- With `useRNF=.FALSE.` (the default) the package changes no result.
- With `useRNF=.TRUE.` the run stops at the end of `RNF_CHECK` with the message
  `RNF: reader not implemented`, after the configuration checks below.

| File | Content | State |
|---|---|---|
| `RNF_OPTIONS.h` | CPP options | no option yet |
| `RNF_SIZE.h` | array bounds | placeholder values, no array uses them yet |
| `RNF.h` | parameters in common blocks | parameters of `data.rnf` |
| `rnf_readparms.F` | reads `data.rnf`; refuses a blank `RNF_file` and a build without NetCDF | done |
| `rnf_check.F` | the other configuration refusals | done; ends with the "not implemented" stop |
| `rnf_summary.F` | prints the parameters and bounds | done; source and target counts to come |
| `rnf_init_fixed.F` | static read of the file, placement of targets on tiles | prints the summary only |
| `rnf_init_varia.F` | zeroes the fields, loads the first records | does nothing |
| `rnf_fields_load.F` | reads records at each step, builds the dense fields | does nothing |
| `rnf_exf_runoff.F` | assigns the exf `runoff` field | does nothing |
| `rnf_tendency_apply.F` | `RNF_TENDENCY_APPLY_T`, `_S` and `_PTR` | do nothing |

Diagnostics, the monitor and the adjoint (TAF) list files are not there yet.

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

An unset timing override means the value in the file is used. Only `RNF_file`
has an effect in the skeleton.

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

## Design

The design, with the reasons for each choice and the source lines it relies
on, is in `docs/package_design.md` of the development repository
<https://github.com/ifenty/MITgcm_new_runoff>. The file format is in
`docs/runoff_schema.md` there.
