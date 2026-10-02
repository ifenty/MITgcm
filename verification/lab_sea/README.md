Labrador Sea Region with Sea-Ice
=========================================

### Primary test Overview:
This example sets up a small (20x16x23) Labrador Sea experiment
coupled to a dynamic thermodynamic sea-ice model (MITgcm Documentation 8.6.2).

The domain of integration spans $`[280, 320]^\circ`$E and $`[46, 78]^\circ`$N.
Horizontal grid spacing is 2 degrees.
The 23 vertical levels and the bathymetry file

```
  bathyFile      = 'bathy.labsea1979'
```
are obtained from the the 2$`^\circ`$ ECCO configuration.

Integration is initialized from annual-mean Levitus climatology

```
 hydrogThetaFile = 'LevCli_temp.labsea1979'
 hydrogSaltFile  = 'LevCli_salt.labsea1979'
```

Surface salinity relaxation is to the monthly mean Levitus climatology

```
 saltClimFile    = 'SSS.labsea1979'
```

Forcing files are a 1979-1999 monthly climatology computed from the
NCEP reanalysis (see [`SEAICE_PARAMS.h`](https://github.com/MITgcm/MITgcm/blob/master/pkg/seaice/SEAICE_PARAMS.h) for units and signs)

```
  uwindFile      = 'u10m.labsea1979'  # 10-m zonal wind
  vwindFile      = 'v10m.labsea1979'  # 10-m meridional wind
  atempFile      = 'tair.labsea1979'  # 2-m air temperature
  aqhFile        = 'qa.labsea1979'    # 2-m specific humidity
  lwdownFile     = 'flo.labsea1979'   # downward longwave radiation
  swdownFile     = 'fsh.labsea1979'   # downward shortwave radiation
  precipFile     = 'prate.labsea1979' # precipitation
```

The experiment uses `pkg/gmredi`, `pkg/kpp`, `pkg/seaice`, and `pkg/exf`.
The test is a 1-cpu, 10-hour integration. Both the atmospheric
state and the open-water surface fluxes are provided by `pkg/exf`.

More `pkg/seaice` test experiments, configured for low and
high-resolution global cube-sphere domains are described
in `MITgcm_contrib/high_res_cube/README_ice`.

### Lab Sea adjoint
The `code_ad` directory provides files required to compile the adjoint
version of this verification experiment.  This verification
experiment uses the 'divided adjoint'.

To compile the adjoint, one must enable the divided adjoint with the
compile-time flag `USE_DIVA`, the location of which is specified in
the file `build/genmake_local`.
To wit,

```
  USE_DIVA=1
```

To compile the adjoint without the divided adjoint, the compile-time
flag `ALLOW_DIVIDED_ADJOINT` in `code_ad/AUTODIFF_OPTIONS.h` should
be changed from

```
  #define ALLOW_DIVIDED_ADJOINT
```
to

```
  #undef ALLOW_DIVIDED_ADJOINT
```

Note: `testreport` builds in the `lab_sea/build` directory which contains
this `genmake_local` file and so it knows to use the divided adjoint.

## Instructions
Navigate to experiment directory

```
  cd MITgcm/verification/lab_sea
```

### 1-CPU forward experiment
Configure and compile the code:
```
  cd build
  ../../../tools/genmake2 -mods ../code [-of my_platform_optionFile]
 [make Clean]
  make depend
  make
  cd ..
```

To run:
```
  cd run
  ln -s ../input/* .
  ln -s ../build/mitgcmuv .
  ./mitgcmuv > output.txt
  cd ..
```

There is comparison output in the directory:
```
  results/output.txt
```

Use matlab script `lookat_ice.m` to compare the output
 with that from `checkpoint51f` sea-ice code:
```
  cd ../../../verification/lab_sea/matlab
  matlab
  lookat_ice
```

### 2-CPU forward experiment
Configure and compile the code:
```
  cd build
  ../../../tools/genmake2 -mpi -mods ../code [-of my_platform_optionFile]
  ln -s ../code/SIZE.h_mpi SIZE.h
 [make Clean]
  make depend
  make
  cd ..
```

To run:
```
  cd run
  ln -s ../input/* .
  mpirun -np 2 ../build/mitgcmuv
  cd ..
```

### 1-CPU adjoint experiment
Configure and compile the code:
```
  cd build
  ../../../tools/genmake2 -mods ../code_ad [-of my_platform_optionFile]
  make adall
  cd ..
```

To run:
```
  cd run
  ln -s ../input_ad/* .
  ./prepare_run
  ln -s ../build/mitgcmuv_ad .
  ./do_run.sh
  cd ..
```

**Note:** `prepare_run` shell script is also used when running `testreport` (see below)
and could be replaced by these 2 commands:
```
  ln -s ../input/* .
  ln -s ../../isomip/input_ad/ones_64b.bin .
```
And the overly simple shell script "do_run.sh" just executes four times
(as specified in file "run_ADM_DIVA", `add_DIVA_runs = 4`) `mitgcmuv_ad`, saving
output in intermediate files "output_adm.txt.diva_0,1,2,3", plus a final time:
```
  mitgcmuv_ad > output_adm.txt
```
where output file `output_adm.txt` can be compared with reference output:
```
  results/output_adm.txt
```

## Secondary tests
In addition to the primary tests described above, 5 secondary forward tests and
and 2 secondary adjoint tests can be run using the same executable as the corresponding
primary tests but with specific input parameter files (in `input.$st\` and `input_ad.$st\`).
The secondary forward tests include alternative seaice model formulations:
free-drift in `input.fd/`; using EVP and `useHB87stressCoupling` in `input.hb87/` ;
with `pkg/salt_plume` in `input.salt_plume/`;
and two other ice-free "North-Altlantic box" set-up (formerly in `verification/natl_box/`)
in `input.natl_box\` and with `pkg/longstep` in `input.longstep/`.
The secondary adjoint tests are simpler version of the primary adjoint test,
without seaice in `input_ad.noseaice/` and without seaice dynamics in `input_ad.noseaicedyn/`.

### Runoff forcing tests
Six more secondary forward tests add `pkg/exf` runoff to the primary test, one for
each way `pkg/exf` can time a dense `runoffFile`. Each `input.rnof_<X>/` holds only the
files that differ from `input/`: `data.exf`, `data` (run length), `data.cal` when the
start date changes, the runoff file or files, and the script `gendata.py` that wrote them.

| Test | Runoff timing | `data.exf` settings | Run (first step to last step) |
| --- | --- | --- | --- |
| `input.rnof_const` | one constant field | `runoffperiod = 0.` | 48 steps, 1979-01-01 01:00 to 1979-01-03 00:00 |
| `input.rnof_daily` | daily records, not repeated | `runoffstartdate1 = 19790101`, `runoffstartdate2 = 000000`, `runoffperiod = 86400.`, `runoffRepCycle = 0.` | 768 steps (32 days), 1979-01-01 to 1979-02-02; reads records 1 to 34, of which 1 to 33 get non-zero weight (record 34 is held with weight 0 at the last forcing time, 1979-02-02 00:00) |
| `input.rnof_month` | repeating monthly climatology: 12 calendar-month records, January to December, repeated every year | `runoffperiod = -12.` | 1464 steps (61 days), 1979-01-01 to 1979-03-03; reads records 12, 1, 2 and 3, all with non-zero weight |
| `input.rnof_month1` | calendar-month records, not repeated; record 1 is December 1978 | `runoffstartdate1 = 19781201`, `runoffstartdate2 = 000000`, `runoffperiod = -1.` | 1464 steps (61 days), 1979-01-01 to 1979-03-03; reads records 1 (December 1978) to 4 (March 1979), all with non-zero weight |
| `input.rnof_clim` | 12 equally spaced records with a repeat cycle | `runoffstartdate1 = 19780116`, `runoffstartdate2 = 120000`, `runoffperiod = 2628000.`, `runoffRepCycle = 31536000.` | 1200 steps (50 days), 1978-12-01 to 1979-01-20; reads records 11, 12, 1 and 2, all with non-zero weight |
| `input.rnof_yearly` | daily records in yearly files `runoff_yearly_YYYY` | `useExfYearlyFields = .TRUE.`, `runoffstartdate1 = 19780101`, `runoffstartdate2 = 000000`, `runoffperiod = 86400.` | 624 steps (26 days), 1978-12-20 to 1979-01-15; reads records 354 to 365 of `runoff_yearly_1978` and 1 to 16 of `runoff_yearly_1979`, of which 354 to 365 and 1 to 15 get non-zero weight (record 16 of 1979 is held with weight 0 at the last forcing time, 1979-01-15 00:00) |

The time step is 3600 s. `startDate_1` in `data.cal` is the date at model time 0
(`pkg/cal/cal_set.F`, lines 235-237), and the runs start at `startTime = 3600.`, so the
first step is at 01:00 on the start date. Forcing is evaluated at the start of each
step, so the last forcing time is `endTime - deltaT` (00:00 on the last day above).
`pkg/exf` always holds the record after the current time, and reads it even when the
time falls exactly on a record and that later record has weight 0. The five long
tests write no pickups (`pChkptFreq = 0.`) and print the monitor every 12 h in
`input.rnof_daily` and `input.rnof_yearly` (`monitorFreq = 43200.`, so the monitor
samples times where the interpolation weight is 0.5, including the year-wrap
interval between 31 December 1978 and 1 January 1979), and once a day in the other
three (`monitorFreq = 86400.`); nothing else in `data` differs from `input/data` apart
from `endTime`.

**Runoff sources.** Runoff is in m/s and is zero except at seven coastal cells (wet cells
with a land neighbour in `bathy.labsea1979`), listed in `runoff_sources.txt` with their
0-based global indices, position, tile and process:

- three adjacent cells on the northern coast, (i, j) = (8, 14), (9, 14) and (10, 14).
  They straddle the boundary between columns 9 and 10, which separates two tiles in
  `code/SIZE.h` and the two processes in `code/SIZE.h_mpi` (`nPx = 2`);
- two adjacent cells, (7, 7) and (7, 8), that straddle the tile boundary between rows 7
  and 8, inside one process;
- one cell on the west Greenland coast, (13, 11), away from tile edges;
- one cell on the ice-free southern coast, (12, 2).

Every record differs from the others. In the 12-record files December is four times
January, and the 1979 yearly file is scaled by 0.45 relative to 1978, so a record or
file chosen wrongly at the turn of the year changes the result. The six records of
`input.rnof_month1` (December 1978 to May 1979) differ from the climatology of
`input.rnof_month`, so the two monthly modes cannot be mistaken for each other. Rates stay below
1e-6 m/s: `input/data.exf` sets `useExfCheckRange = .TRUE.`, and at the first time step
`EXF_CHECK_RANGE` stops the run if runoff on a wet cell is negative or above 1e-6 m/s
(`pkg/exf/exf_check_range.F`, lines 177-191 and 211-216, called from
`pkg/exf/exf_getforcing.F`, lines 346-349).

To regenerate the runoff files (needs python3 with numpy; the output is deterministic):

```
  cd MITgcm/verification/lab_sea/input.rnof_daily
  python3 gendata.py
```

**Timing conventions of `pkg/exf` used by these tests.** Line numbers refer to the
source files in `pkg/exf/` unless another package is named.

- *Constant field* (`runoffperiod = 0.`, the default set in `exf_readparms.F`, line 440).
  `EXF_INIT_FLD` reads record 1 once at initialisation (`exf_init_fld.F`, lines 98-126),
  and `EXF_SET_FLD` skips the field afterwards (`exf_set_fld.F`, line 120).
- *Start time* (`exf_getffield_start.F`, lines 80-104). For a positive period,
  `runoffstartdate1` and `runoffstartdate2` give the date of record 1; for
  `runoffperiod = -1.` without yearly files they give the month of record 1. Without yearly
  files it is converted to model time. With yearly files only its offset from 1 January
  is kept, and that offset applies to every yearly file.
- *Records without a repeat cycle* (`exf_getffieldrec.F`, lines 119-132). With
  `t` the time since record 1, the record before is `count0 = INT((t+0.5)/period) + 1`,
  the record after is `count0 + 1`, and the weight of `count0` is
  `1 - MOD(t,period)/period` (line 149). The field is interpolated linearly between the
  two records (`exf_set_fld.F`, lines 299-314). A time before record 1 stops the run.
- *Repeat cycle* (`exf_getffieldrec.F`, lines 134-146). When `runoffRepCycle` is
  positive, `t` is taken modulo the cycle, so the record after the last one is
  record 1. `runoffRepCycle` defaults to `repeatPeriod` (`exf_readparms.F`, line 951),
  which `input/data.exf` sets to 31622400 s. `input.rnof_daily` therefore sets
  `runoffRepCycle = 0.` to read its records without repeating, and `input.rnof_clim`
  sets a cycle of 365 days with a period of one twelfth of it.
- *Repeating monthly climatology* (`exf_set_fld.F`, lines 133-140).
  `runoffperiod = -12.` means the file holds 12 records, January to December, which
  `cal_GetMonthsRec` places at the middle of each calendar month and repeats every year
  (`pkg/cal/cal_getmonthsrec.F`, lines 106-117 and 134-213). Between two mid-month
  times the field is interpolated linearly. The start date and the repeat cycle are
  not used.
- *Calendar months without repetition* (`exf_set_fld.F`, lines 142-153).
  `runoffperiod = -1.` uses the same mid-month times and interpolation, but the file
  holds one record per consecutive calendar month. Without yearly files, record 1 is
  the month in which the runoff start date falls, and the record of year `y`, month `m`
  is `(y - yy)*12 + m - mm + 1`, with `yy` and `mm` the year and month of the start
  date (`exf_getmonthsrec.F`, lines 58-69). A run that starts before the middle of a
  month needs the record of the month before, which is why `input.rnof_month1` starts
  its records in December 1978. Any other negative period stops the run
  (`exf_set_fld.F`, lines 154-160).
- *Yearly files* (`exf_getffieldrec.F`, lines 152-190). `useExfYearlyFields = .TRUE.`
  applies to every `pkg/exf` field. Records are counted within each year from the
  start date's offset, and when the record after `count0` falls beyond the end of the year it is
  record 1 of the next year's file. The file name is the given name followed by `_YYYY`
  (`exf_getyearlyfieldname.F`, lines 46-55). Yearly files cannot be combined with a
  non-zero `repeatPeriod` (`exf_check.F`, lines 71-84). `input.rnof_yearly/data.exf`
  therefore follows `input/data.exf_YearlyFields` for the other fields (record 1 at
  16 January 18:00, no repeat period), and its `prepare_run` links the climatological
  forcing files `<name>.labsea1979` as `<name>.labsea_1978` and `<name>.labsea_1979`.
  The atmospheric forcing of this test thus differs slightly from that of `input/`.

The reference outputs are `results/output.rnof_<X>.txt`. In each test the runoff
statistics printed by the monitor (`exf_runoff_max`, `_min`, `_mean` and `_sd`) equal,
to 1e-12 relative, those of the field interpolated from the input records under the
conventions above, at every monitor time.

### Instruction to run secondary tests
Run the testscript _forward_ experiments:

```
  cd MITgcm/verification
  ./testreport -t lab_sea [-of my_platform_optionFile]
```

Standard testreport output, with all secondary tests:
```
default 10  ----T-----  ----S-----  ----U-----  ----V-----  --PTR 01--  --PTR 02--  --PTR 03--  --PTR 04--  --PTR 05--
G D M    c        m  s        m  s        m  s        m  s        m  s        m  s        m  s        m  s        m  s
e p a R  g  m  m  e  .  m  m  e  .  m  m  e  .  m  m  e  .  m  m  e  .  m  m  e  .  m  m  e  .  m  m  e  .  m  m  e  .
n n k u  2  i  a  a  d  i  a  a  d  i  a  a  d  i  a  a  d  i  a  a  d  i  a  a  d  i  a  a  d  i  a  a  d  i  a  a  d
2 d e n  d  n  x  n  .  n  x  n  .  n  x  n  .  n  x  n  .  n  x  n  .  n  x  n  .  n  x  n  .  n  x  n  .  n  x  n  .

Y Y Y Y>11<13 16 16 16 16 16 16 14 13 13 13 16 16 16 12 14 22 16 16 16 22 16 16 16 pass  lab_sea
Y Y Y Y>12<16 16 16 16 16 16 16 16 16 13 13 16 14 16 13 16 22 16 16 16 22 16 16 16 pass  lab_sea.fd
Y Y Y Y> 4< 9 10  9  9 16 13 11  8  6  8  4  5  7  8  4  6 22  7  7  7 22  6  7  7 FAIL  lab_sea.hb87
Y Y Y Y 11 16 16 16 14 16 16 16 16 16 13 12 14 16 13 12 14 16 16 16>16<pass  lab_sea.longstep
Y Y Y Y>11<16 16 16 16 16 16 16 16 13 12 12 14 14 13 12 14 pass  lab_sea.natl_box
Y Y Y Y>13<16 16 16 16 16 16 16 14 16 14 13 14 16 14 13 16 22 16 16 16 22 16 16 16  pass  lab_sea.salt_plume
```

**Note:** Some differences in accuracy occur across different platforms as seen
here for secondary test "lab_sea.hb87".

Run the testscript _adjoint_ experiments:

```
  cd MITgcm/verification
  ./testreport -t lab_sea -ad [-of my_platform_optionFile]
```

Standard adjoint testreport output, with all secondary tests:
```
Adjoint generated by TAF Version 6.5.1

default    10     ----T-----  ----S-----  ----U-----  ----V-----
G D M    C  A  F        m  s        m  s        m  s        m  s
e p a R  o  d  D  m  m  e  .  m  m  e  .  m  m  e  .  m  m  e  .
n n k u  s  G  G  i  a  a  d  i  a  a  d  i  a  a  d  i  a  a  d
2 d e n  t  r  r  n  x  n  .  n  x  n  .  n  x  n  .  n  x  n  .

Y Y Y Y 16>16<16 16 14 16 16 16 16 16 16 16 16 16 14 14 16 16 16 pass  lab_sea  (e=0, w=0, lfd=1, dop=1, sm=1)
Y Y Y Y 14>15< 6 16 16 16 13 16 13 13 11 13 13 12 13 11 11 12 12 pass  lab_sea.noseaice
Y Y Y Y 16>13<16 16 16 14 16 16 13 13 16 16 16 16 14 14 16 14 16 pass  lab_sea.noseaicedyn
```
