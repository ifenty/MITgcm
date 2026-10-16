.. _sub_phys_pkg_rnf:

RNF: Sparse runoff forcing by source
------------------------------------

.. _ssub_phys_pkg_rnf_intro:

Introduction
~~~~~~~~~~~~

:filelink:`pkg/rnf` reads river, glacier and other runoff from a NetCDF
file that is organized by *source*. A dense :code:`runoffFile` of
:filelink:`pkg/exf` stores one value per surface grid cell and record,
and at high resolution almost all of them are zero. In the sparse file
each source has an identifier, a list of the ocean cells it feeds, the fraction of
its flux that goes into each of them (the fractions of one source sum to 1) and
time series on one shared time axis. The time series are stored once per
source and not once per cell.

The runoff is a volume flux (m\ :sup:`3`/s). It can carry a temperature, a
salinity (0 if the file has none) and any number of passive tracers, up to the
compile-time limit :code:`RNF_nTr` (:numref:`ssub_phys_pkg_rnf_config`). The
package works on every MITgcm grid type, including the exch2 cubed sphere and
LLC layouts: the file addresses cells by their index in the global 2D array
that MITgcm reads from a dense input file on the same grid, and each
tile keeps only the target cells it owns. Every process reads the whole file.

:filelink:`pkg/rnf` needs :filelink:`pkg/exf`, whose record-selection routines it
uses for the time handling (:numref:`ssub_phys_pkg_rnf_time`) and whose
:code:`ALLOW_RUNOFF` option it requires. It is described in the :ref:`exf
section <sub_phys_pkg_exf>` only as far as :code:`useExfCheckRange` is concerned
(:numref:`ssub_phys_pkg_rnf_exfrange`).

CPP options and package dependencies are in
:numref:`ssub_phys_pkg_rnf_config`, the input file in
:numref:`ssub_phys_pkg_rnf_input`, run-time parameters in
:numref:`ssub_phys_pkg_rnf_runtime`, the time handling in
:numref:`ssub_phys_pkg_rnf_time`, diagnostics and the monitor in
:numref:`ssub_phys_pkg_rnf_diagnostics`, example configurations in
:numref:`ssub_phys_pkg_rnf_examples`, the Python tools that write and check the
input file in :numref:`ssub_phys_pkg_rnf_tools`, and how the package is tested,
and what that testing does not cover, in
:numref:`ssub_phys_pkg_rnf_verification`.

Two routes
++++++++++

The volume, and the heat, salt and tracers it carries, enter the model by one
of two routes per run, selected by :code:`RNF_applyAsRain` in
:code:`data.rnf`:

**Generalized route** (:code:`RNF_applyAsRain = .FALSE.`, the default).
Runoff is treated as a lateral input and not as rain. For each target
entry the volume flux :code:`rhoConstFresh*flux*fraction` (kg/s) is added to the
model's mass source :code:`addMass` at the level of that entry
(:code:`RNF_ADD_MASS`). The temperature, salinity and tracers are tendencies at
the same level, :math:`[(mX) - m_X X(k)]\,\mu/(\Delta r_F\, h_C)`, relative
to the ambient value :math:`X(k)` of the target cell. Here :math:`(mX)` is the
sum over the sources feeding the cell of mass flux times the property of the
source, :math:`m_X` is the mass flux of the sources that have the property,
and :math:`\mu` is :code:`mass2rUnit`. A source without a temperature
therefore enters at the ambient temperature, and a salinity that the file lacks
is 0 and dilutes. The model's own generic :code:`temp_addMass` /
:code:`salt_addMass` term acts on the runoff share of :code:`addMass` as well,
and the package cancels it for that share. With a schema 1.1 file a
target can be any wet level of its column, and one source can feed several
levels of one cell. The route needs :code:`#define ALLOW_ADDFLUID`,
:code:`selectAddFluid = 1` or :code:`2` and a blank :code:`addMassFile`; with
:code:`staggerTimeStep` the model's warning about a one-step lag in the effect
of :code:`addMass` on temperature and salinity does not appear. The exf
field :code:`runoff` stays zero on this route.

**Legacy surface route** (:code:`RNF_applyAsRain = .TRUE.`). The package fills
the exf field :code:`runoff` (m/s) with the sparse runoff, and everything
downstream is the code that a dense :code:`runoffFile` uses. The temperature
goes through the exf field :code:`runoftemp` and the :code:`Qnet` term of
:code:`EXF_MAPFIELDS` (this needs :code:`#define ALLOW_RUNOFTEMP`), the salt
is added to :code:`surfaceForcingS` and the tracers to
:code:`surfaceForcingPTr`. The file may only have targets that are the
surface cell of their column.

What the choice means for other packages: :filelink:`pkg/seaice`,
:filelink:`pkg/thsice` and :filelink:`pkg/bling` read the exf array
:code:`runoff`, so on the legacy route they see the sparse runoff and on the
generalized route they do not. On the generalized route the runoff is also not
part of :code:`EmPmR`, so the virtual salt flux, :code:`selectBalanceEmPmR` and
the surface buoyancy flux and non-local transport of :filelink:`pkg/kpp` do not
include it. The package design records this as accepted and, for
:filelink:`pkg/kpp`, as not exercised by any test (:numref:`ssub_phys_pkg_rnf_verification`).

What the package does not do
++++++++++++++++++++++++++++

- Runoff is read from one file per run (yearly sets of one series,
  :code:`RNF_useYearlyFiles`, count as one).
- Targets under an ice shelf, or beyond an open boundary, are refused, and so
  is :code:`SHI_update_kTopC` with :filelink:`pkg/shelfice`, because a moving
  ice-shelf edge could cover a target.
- A dense :code:`runoffFile`, :code:`runoftempfile` or :code:`runoffconst`
  together with the sparse file is refused: they are mutually exclusive.
- The exf input scaling :code:`exf_inscal_runoff` and the :code:`runoff_exfremo_*`
  parameters apply to the dense file only. The sparse flux is used as stored,
  in m\ :sup:`3`/s, and :code:`exf_outscal_sflux` must be 1.
- With :filelink:`pkg/longstep` and :code:`LS_nIter` other than 1, a file that
  feeds a passive tracer is refused.
- Yearly *sampling* (one record per calendar year) is refused, because exf has
  no such mode. Yearly *files* (one file per year of a fixed-period series) are
  supported.
- The adjoint (TAF) list files of the package do not exist, and
  adjoint use of the package is untested.

.. _ssub_phys_pkg_rnf_config:

RNF configuration and CPP options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As for any package, :filelink:`pkg/rnf` is compiled in by adding :code:`rnf` to
``packages.conf`` (or with ``genmake2 -enable=rnf``). The file
:filelink:`pkg/pkg_depend` contains the line ``rnf +exf``, so that
:code:`exf` is compiled in as well. The package needs NetCDF: the file is read
with the Fortran-77 NetCDF interface, and a build in which ``genmake2`` did not
find NetCDF (:code:`HAVE_NETCDF` undefined) stops at run time when
:code:`useRNF = .TRUE.`.

:filelink:`pkg/rnf/RNF_OPTIONS.h` has no package-specific CPP option today.
The options that matter are those of other packages:

.. list-table::
   :header-rows: 1
   :widths: 22 20 58

   * - Option
     - File
     - What it enables for pkg/rnf
   * - :code:`ALLOW_RUNOFF`
     - :filelink:`pkg/exf/EXF_OPTIONS.h`
     - The exf field :code:`runoff` and its record selection. Required on both
       routes (``RNF_CHECK`` stops without it). It is defined in the
       committed :filelink:`pkg/exf/EXF_OPTIONS.h`.
   * - :code:`ALLOW_RUNOFTEMP`
     - :filelink:`pkg/exf/EXF_OPTIONS.h`
     - The exf field :code:`runoftemp` and its heat term in :code:`EXF_MAPFIELDS`.
       Needed on the legacy route for a file with :code:`runoff_temperature`
       (stopped otherwise; :code:`RNF_useTemp = .FALSE.` reads the file
       without it). Not needed on the generalized route. Not defined in the
       committed :filelink:`pkg/exf/EXF_OPTIONS.h`.
   * - :code:`ALLOW_ATM_TEMP`
     - :filelink:`pkg/exf/EXF_OPTIONS.h`
     - On the legacy route, a file with :code:`runoff_temperature` together
       with :code:`temp_EvPrRn` set needs this option, because exf cancels
       the model's :code:`temp_EvPrRn` term for runoff only with it. The
       combination without it is refused. :code:`temp_EvPrRn` unset (the default)
       or :code:`RNF_useTemp = .FALSE.` avoids the refusal.
   * - :code:`ALLOW_ADDFLUID`
     - :filelink:`model/inc/CPP_OPTIONS.h`
     - The mass source :code:`addMass` of the generalized route. Required for
       :code:`RNF_applyAsRain = .FALSE.`. It is undefined in the committed
       :filelink:`model/inc/CPP_OPTIONS.h`, so an experiment that uses the
       default route needs a ``CPP_OPTIONS.h`` in its code directory that defines it.
   * - :code:`USE_OLD_EXTERNAL_FORCING`
     - :filelink:`model/inc/CPP_OPTIONS.h`
     - The legacy route has no hooks in the old forcing routines and is
       refused if the option is defined.

The array bounds are in :filelink:`pkg/rnf/RNF_SIZE.h`. They are
compile-time maxima with run-time counts inside them, and a count above its
bound stops the run and prints the value that is needed:

.. list-table::
   :header-rows: 1
   :widths: 22 12 66

   * - Parameter
     - Value
     - Meaning
   * - :code:`RNF_nSrcTile`
     - 2000
     - maximum number of sources on one tile
   * - :code:`RNF_nTgtTile`
     - 10000
     - maximum number of target entries on one tile
   * - :code:`RNF_nBuf`
     - 1000
     - chunk length of the NetCDF reads and of the global fraction sum
   * - :code:`RNF_nTr`
     - 5
     - maximum number of runoff tracers (``runoff_ptracer_<NAME>`` variables
       in the file)

The constants of the model contract are in :filelink:`pkg/rnf/RNF.h`:
:code:`RNF_fracTol` (1e-6, tolerance of the fraction sum of a source),
:code:`RNF_areaTol` (1e-4, relative tolerance between
:code:`target_cell_area` and the model's :code:`rA`), :code:`RNF_lonLatTol`
(0.5, the largest distance of a stored :code:`target_lon`,
:code:`target_lat` from the cell centre, as a fraction of the grid spacing of the cell),
:code:`RNF_fluxMax` (1e30, a flux above it counts as a missing value),
:code:`RNF_srcFluxMax` (1e7 m\ :sup:`3`/s) and :code:`RNF_cellVolMax` (0.2),
which are described in :numref:`ssub_phys_pkg_rnf_exfrange`.

.. _ssub_phys_pkg_rnf_input:

The sparse runoff file
~~~~~~~~~~~~~~~~~~~~~~

The format is specified in ``docs/runoff_schema.md`` (version 1.1) of the
development repository for this package,
https://github.com/ifenty/MITgcm_new_runoff, which also holds the model
contract (``docs/model_contract.md``). This section only summarizes it; where
they differ, the schema document and the code govern.

The file is NetCDF-4. It has the dimensions :code:`time`, :code:`source` and
:code:`target`, and the model reads these variables:

.. list-table::
   :header-rows: 1
   :widths: 28 22 50

   * - Variable
     - Dimensions
     - Content
   * - :code:`source_id`
     - (source, id_strlen)
     - fixed-length character ids (letters, digits, ``_``, ``-``, ``.``;
       at most 64 characters), used in error messages
   * - :code:`target_source`
     - (target)
     - 0-based index of the source of each (source, cell) pair
   * - :code:`target_cell`
     - (target)
     - 0-based global cell index ``i + mitgcm_grid_nx * j`` in the global 2D
       layout that MITgcm reads for a dense input file (for exch2 grids the
       global I/O layout, which depends on :code:`W2_mapIO`)
   * - :code:`target_fraction`
     - (target)
     - share of the source flux sent to the cell; each in [0, 1], and the
       fractions of one source sum to 1 within 1e-6
   * - :code:`target_level`
     - (target), optional
     - schema 1.0: only 1, the surface cell. Schema 1.1: absolute 1-based
       level index, which must be a wet level of the column. Absent: surface
   * - :code:`target_cell_area`
     - (target), optional
     - area of the cell on the grid the file was built for (m\ :sup:`2`),
       compared with the model's :code:`rA`
   * - :code:`target_lon`, :code:`target_lat`
     - (target), optional
     - cell centre, only used to check the placement against :code:`XC`,
       :code:`YC`
   * - :code:`time`, :code:`time_bnds`
     - (time), (time, nv)
     - record times and intervals, in CF units (``<unit> since <date>``) and
       calendar
   * - :code:`runoff_flux`
     - (time, source)
     - volume flux in m\ :sup:`3` s\ :sup:`-1`; required; a missing value is not
       allowed
   * - :code:`runoff_temperature`
     - (time, source), optional
     - degC; a missing value means that source enters at the ambient
       (reference) temperature
   * - :code:`runoff_salinity`
     - (time, source), optional
     - model salinity units; if absent the salinity is 0
   * - :code:`runoff_ptracer_<NAME>`
     - (time, source), optional
     - concentration per unit volume of runoff of the passive tracer whose
       :code:`PTRACERS_names` entry is ``<NAME>``

The global attributes the model reads are :code:`mitgcm_runoff_schema_version`
(``"1.0"`` or ``"1.1"``), :code:`mitgcm_grid_nx`, :code:`mitgcm_grid_ny`
(these must equal the model's global I/O layout), :code:`mitgcm_time_sampling`
(``constant``, ``fixed`` or ``monthly``; ``yearly`` is refused by the model) and,
where applicable, :code:`mitgcm_time_period` and :code:`mitgcm_time_repeat`.
Different sources may feed the same cell: the model adds their volumes and
weights their temperature, salinity and tracers by flux. Several
(source, cell) pairs with different levels are allowed in a 1.1 file.
The model never converts units and does not unpack ``scale_factor`` /
``add_offset``.

The package checks the file in :code:`RNF_INIT_FIXED` and when it reads each record,
and a violation stops the run with a message that names the source id or the
table entry. The refusals include a major version other than 1, a grid
size that differs from the model's, an index out of range, a fraction outside
[0, 1] or a fraction sum other than 1, a target on land, beyond an open boundary,
under an ice shelf or with a wrong :code:`target_cell_area`, a missing flux,
and a ``runoff_ptracer_<NAME>`` that matches no entry of
:code:`PTRACERS_names`. The list is in ``pkg/rnf/README.md`` of the source tree. The
same rules are applied before a run by the checker described in
:numref:`ssub_phys_pkg_rnf_tools`.

.. _ssub_phys_pkg_rnf_runtime:

RNF run-time parameters
~~~~~~~~~~~~~~~~~~~~~~~

:filelink:`pkg/rnf` is switched on at run time by :code:`useRNF = .TRUE.` in the
namelist :code:`PACKAGES` of ``data.pkg``; with :code:`useRNF` false the
package changes no result. The parameters are in the namelist
:code:`RNF_PARM01` of ``data.rnf``, which is read by
:filelink:`pkg/rnf/rnf_readparms.F`; a missing file ``data.rnf`` is an error when
:code:`useRNF = .TRUE.`. Their declarations are in :filelink:`pkg/rnf/RNF.h`.

.. list-table:: Run-time parameters and default values; all parameters are in namelist group ``RNF_PARM01``
   :header-rows: 1
   :widths: 22 17 61
   :name: tab_phys_pkg_rnf_runtimeparms

   * - Parameter
     - Default
     - Description
   * - :code:`RNF_file`
     - :kbd:`' '`
     - CHARACTER. Sparse runoff NetCDF file, or its base name with
       :code:`RNF_useYearlyFiles`. A blank name stops the run.
   * - :code:`RNF_holdRecord`
     - FALSE
     - LOGICAL. FALSE: linear interpolation in time between records, as
       :code:`pkg/exf`. TRUE: hold each record over its interval.
   * - :code:`RNF_startDate1`, :code:`RNF_startDate2`
     - 0
     - INTEGER. Start date of the series (YYYYMMDD, HHMMSS); overrides the file; 0
       means not set. Needs :code:`pkg/cal`.
   * - :code:`RNF_startTime`
     - unset
     - REAL. Start time of the series (s); overrides the file. For a run
       without :code:`pkg/cal`; refused with it.
   * - :code:`RNF_period`
     - unset
     - REAL. Record period; overrides the file: 0 constant, > 0 seconds, -12
       monthly climatology, -1 calendar months.
   * - :code:`RNF_repCycle`
     - unset
     - REAL. Repeat cycle (s); overrides the file.
   * - :code:`RNF_useYearlyFiles`
     - FALSE
     - LOGICAL. Append ``_YYYY`` to :code:`RNF_file` by model year: the file read
       is ``<base>_YYYY.nc`` (a trailing ``.nc`` of the name given is replaced).
       Needs :code:`pkg/cal`.
   * - :code:`RNF_useTemp`
     - TRUE
     - LOGICAL. Apply :code:`runoff_temperature` if the file has it; FALSE treats
       it as absent.
   * - :code:`RNF_useSalt`
     - TRUE
     - LOGICAL. Apply :code:`runoff_salinity` if the file has it; FALSE treats it
       as absent.
   * - :code:`RNF_usePtracers`
     - TRUE
     - LOGICAL. Apply the :code:`runoff_ptracer_*` variables if the file has them;
       FALSE treats them as absent, and their names are then not matched to the
       passive tracers.
   * - :code:`RNF_monFreq`
     - :code:`monitorFreq`
     - REAL. Interval (s) of the monitor output (:numref:`ssub_phys_pkg_rnf_diagnostics`).
   * - :code:`RNF_debugLev`
     - :code:`debugLevel`
     - INTEGER. Message level. At 3 (:code:`debLevC`) the package prints, at every
       step, the two records it selected and their weight.
   * - :code:`RNF_applyAsRain`
     - FALSE
     - LOGICAL. FALSE: the generalized route; TRUE: the legacy surface route
       (see the introduction).

An unset timing override means that the value in the file is used. The
resolved values are printed by the package summary as :code:`RNF_recPeriod`,
:code:`RNF_recStart`, :code:`RNF_recCycle`, :code:`RNF_recDate1` and
:code:`RNF_recDate2`, together with which of the optional series the file has
(:code:`RNF_hasTemp`, :code:`RNF_hasSalt`, :code:`RNF_nTrUse`) and the ptracer
each runoff tracer feeds.

Settings elsewhere that matter:

- ``data.pkg``: :code:`useEXF = .TRUE.` and :code:`useRNF = .TRUE.`;
  :code:`usePTRACERS = .TRUE.` if the file has tracers; :code:`useCAL` for
  calendar dates (needed by :code:`RNF_useYearlyFiles` and the calendar-month
  modes).
- ``data.exf``: :code:`runoffFile` and :code:`runoftempfile` must be blank (the
  default), :code:`runoffconst` must be 0 and :code:`exf_outscal_sflux` must be
  1. The dense timing parameters (:code:`runoffperiod` and the other
  :code:`runoff*` ones) are not used by the sparse file.
- ``data``, for the generalized route: :code:`selectAddFluid = 1` or :code:`2`
  and :code:`addMassFile` blank; :code:`rigidLid` is an error with
  :code:`selectAddFluid`.

.. _ssub_phys_pkg_rnf_exfrange:

Range checks: interaction with pkg/exf
++++++++++++++++++++++++++++++++++++++

A sparse point source puts a whole river into one cell, and a rate bound in
m/s is a statement about the grid and not about the water. Two tests of
:code:`EXF_CHECK_RANGE` (:code:`useExfCheckRange`, default true) are
therefore conditioned on :code:`useRNF`, as described in
:numref:`ssub_phys_pkg_exf_runtime`: the upper bound of 10\ :sup:`--6` m/s on
:code:`runoff` is skipped, and the :code:`sflux` bound is applied to
:code:`sflux + runoff`, so that it still tests :code:`evap - precip`. The
negative-:code:`runoff` test and all other range checks are unchanged. In their
place :filelink:`pkg/rnf` has two bounds of its own, both printed by the package summary:

- :code:`RNF_srcFluxMax` = 1e7 m\ :sup:`3`/s on the absolute flux of each
  source, checked on **every record as it is read**. It is a filter for unit
  errors (a flux per year instead of per second, or mm instead of m), not a
  bound on the applied field: several sources collapsed onto one cell can add up
  to far more.
- :code:`RNF_cellVolMax` = 0.2, a dimensionless bound checked at every step: one
  time step of runoff may add at most 0.2 of the volume of the target cell (of
  its top layer on the legacy route, of the cell at the target's own level on
  the generalized route). It stops the run and names the cell
  (``i,j,bi,bj`` and ``XC,YC``), the applied value and the limit. The thickness
  used is the live one, so on a grid with a nonlinear free surface the limit moves
  with the state, and a file can pass at the first step and be stopped later.
  It does not see the :code:`xx_runoff` control increment and it bounds neither
  the sign nor the properties of the water.

.. _ssub_phys_pkg_rnf_time:

Time handling
~~~~~~~~~~~~~

The package picks the records that bracket the model time with the record
routines of :filelink:`pkg/exf` (:code:`EXF_GetFFieldRec`,
:code:`cal_GetMonthsRec`, :code:`EXF_GetMonthsRec`), so that a sparse run uses
the selection code of a dense :code:`runoffFile` run. The file states
its timing in :code:`mitgcm_time_sampling` and :code:`mitgcm_time_repeat`;
:code:`RNF_TIME_SETUP` maps them to the exf convention, and the
:code:`RNF_period`, :code:`RNF_repCycle`, :code:`RNF_startDate1/2` and
:code:`RNF_startTime` of ``data.rnf`` override the file:

.. list-table::
   :header-rows: 1
   :widths: 24 16 60

   * - Mode
     - exf period
     - File
   * - Constant
     - 0
     - ``constant`` sampling, one record, read once.
   * - Fixed period, not repeated
     - seconds
     - ``fixed`` sampling with :code:`mitgcm_time_period`; ``repeat = none``.
       Interpolation between records, as exf.
   * - Fixed period with a repeat cycle (a climatology)
     - seconds
     - ``fixed`` with ``repeat = annual``; the cycle is the span of
       :code:`time_bnds`, and is anchored on the real date of record 1.
   * - Monthly climatology
     - -12
     - ``monthly`` with ``repeat = annual``: 12 records, January to December,
       repeated every model year.
   * - Calendar months
     - -1
     - ``monthly`` with ``repeat = none``: consecutive calendar months.
   * - Yearly files
     - as the series
     - :code:`RNF_useYearlyFiles = .TRUE.`: ``<base>_YYYY.nc`` for the model year;
       each file holds a fixed-period series.

All series share **one set of record weights** per step. The flux, the
temperature, the salinity and every tracer are each interpolated with it
(or held, see below), and the products are formed afterwards: a cell receives
:math:`m(t)=\sum_s m_s(t)` and :math:`(mX)(t)=\sum_s m_s(t)X_s(t)`, with the
interpolated values of the sources feeding it. Consequently the temperature, salinity and tracer
of a source follow the same time axis as its flux and cannot be on another.
A source whose temperature is missing in the records used enters at the ambient
temperature, and the diagnostic :code:`RNFtemp` reports 0 for a cell where no
source has one.

With :code:`RNF_holdRecord = .TRUE.` the package does not interpolate: it holds
the record whose interval contains the model time, for the flux and for every
property alike. For a fixed period the interval runs from the time of the record to
the time of the next one; for the monthly modes it is the calendar month, which
is not the nearer of two mid-month times. exf has no such mode, so there is no
dense counterpart.

If the file's reference date precedes the start of the :code:`pkg/cal`
calendar (15 October 1582), give the date of record 1 as
:code:`RNF_startDate1` / :code:`RNF_startDate2`. A run without
:code:`pkg/cal` reads the file's times as model time in seconds from the
reference date of the units.

.. _ssub_phys_pkg_rnf_diagnostics:

RNF diagnostics and monitor
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Diagnostics output is available via the diagnostics package (see
:numref:`outp_pack`) and is registered in
:filelink:`pkg/rnf/rnf_diagnostics_init.F`. Five fields describe the input and
exist on both routes:

::

    ---------+----+----+----------------+-----------------
     <-Name->|Levs|grid|<--  Units   -->|<- Tile (max=80c)
    ---------+----+----+----------------+-----------------
     RNFvflx |  1 |SM  |m/s             |runoff volume flux per area, before controls (= EXFroff)
     RNFmflx |  1 |SM  |kg/m^2/s        |runoff mass flux, rhoConstFresh*RNFvflx
     RNFtemp |  1 |SM  |degC            |flux-weighted runoff temperature of the sources with one; 0: none has one, or no runoff
     RNFsaln |  1 |SM  |g/kg            |flux-weighted runoff salinity; 0: no runoff
     RNFnsrc |  1 |SM  |1               |number of target entries feeding the cell

:code:`RNFnsrc` counts the target entries, which is the number of sources
for a schema 1.0 file; in a 1.1 file a source spread over several levels of
the cell counts once per level.

The generalized route (:code:`RNF_applyAsRain = .FALSE.`) has four
three-dimensional fields. Each is the physical term of the route,
including the share that the model's own :code:`addMass` term provides:

::

    ---------+----+----+----------------+-----------------
     <-Name->|Levs|grid|<--  Units   -->|<- Tile (max=80c)
    ---------+----+----+----------------+-----------------
     RNFaddM | Nr |SM  |kg/s            |runoff share of addMass, >0 adds water
     RNFgT   | Nr |SM  |degC/s          |runoff temperature tendency
     RNFgS   | Nr |SM  |g/kg/s          |runoff salinity tendency
     RNFtrNN | Nr |SM  |tracer units/s  |tendency of runoff tracer NN (NN = 01 to 05)

The legacy route (:code:`RNF_applyAsRain = .TRUE.`) has these surface fields,
each the contribution of the runoff to the field it feeds:

::

    ---------+----+----+----------------+-----------------
     <-Name->|Levs|grid|<--  Units   -->|<- Tile (max=80c)
    ---------+----+----+----------------+-----------------
     RNFqnet |  1 |SM  |W/m^2           |runoff heat added to Qnet, >0 decreases theta (only with ALLOW_RUNOFTEMP)
     RNFsflx |  1 |SM  |g/m^2/s         |runoff salt flux in surfaceForcingS, >0 increases salinity
     RNFtfNN |  1 |SM  |kg/m^2/s.trc    |surface flux of runoff tracer NN, >0 increases it (NN = 01 to 05)

:code:`NN` numbers the runoff tracers in the order of the
``runoff_ptracer_<NAME>`` variables of the file, not the ptracers. The package
summary prints the name and ptracer number of each, in lines such as
``RNF_SUMMARY: runoff tracer  1 is runoff_ptracer_dye, applied to ptracer  1``.
The diagnostics for all :code:`RNF_nTr` tracers are registered whether the file
has them or not.

**Monitor.** :code:`RNF_MONITOR` (:filelink:`pkg/rnf/rnf_monitor.F`) prints a
block ``Begin MONITOR RNF statistics`` to standard output every
:code:`RNF_monFreq` seconds, in the style of the exf monitor. The names are
prefixed ``rnf`` and printed as ``%MON rnf_...``:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Name
     - Meaning
   * - :code:`rnf_tsnumber`, :code:`rnf_time_sec`
     - iteration number and model time (s)
   * - :code:`rnf_vflx_max`, :code:`_min`, :code:`_mean`, :code:`_sd`, :code:`_del2`
     - statistics of the applied volume flux per area (m/s) over the ocean
   * - :code:`rnf_volsum`
     - global sum of the volume flux, :math:`\sum RNF\_vflx \cdot rA` (m\ :sup:`3`/s)
   * - :code:`rnf_heatsum`
     - global sum of :math:`C_p \sum m T\, rA`, the heat the temperatures
       of the sources carry (W)
   * - :code:`rnf_saltsum`
     - global sum of :math:`\sum m S\, rA` (g/s)
   * - :code:`rnf_trc01sum`, ...
     - one per runoff tracer in use, global sum of :math:`\sum m C_n\, rA`
       (tracer units kg/s)
   * - :code:`rnf_nsrc`, :code:`rnf_ntgt`
     - number of sources in the file, and number of target entries owned over
       all processes

The sums are exact global sums, so the output text of a one-process and a
two-process run is the same. The sums refer to the applied input fields, not
to what the model does with them; they are the same on both routes.

.. _ssub_phys_pkg_rnf_examples:

Examples
~~~~~~~~

Minimal configuration of the legacy route, as in the verification experiment
:filelink:`lab_sea <verification/lab_sea>`, case
``input.rnof_sp_const`` (the runoff of the dense case ``input.rnof_const`` in
sparse form: 4 sources, 7 target cells).

``packages.conf`` (``code/``) contains the lines ``exf`` and ``rnf``:

::

    oceanic
    cd_code
    exf
    rnf
    seaice
    salt_plume
    ptracers
    longstep
    diagnostics
    mnc

``data.pkg``:

::

     &PACKAGES
      useGMRedi = .TRUE.,
      useKPP    = .TRUE.,
      useEXF    = .TRUE.,
      useCAL    = .TRUE.,
      useSEAICE = .TRUE.,
      useDiagnostics = .TRUE.,
      useMNC    = .TRUE.,
      useRNF    = .TRUE.,
     &

``data.exf`` has ``runoffFile = ' '`` (no dense runoff), and ``data.rnf``:

::

     &RNF_PARM01
      RNF_file = '../input.rnof_const/runoff_sparse.nc',
      RNF_debugLev = 3,
      RNF_applyAsRain = .TRUE.,
     &

The timing parameters are left unset because the file says
``mitgcm_time_sampling = "constant"``. The file is built from the dense
file by ``verification/lab_sea/input.rnof_const/gen_sparse.py const``, which uses the
converter of :numref:`ssub_phys_pkg_rnf_tools`.

The cubed-sphere case
:filelink:`global_ocean.cs32x15 <verification/global_ocean.cs32x15>`,
``input.rnof_sp_icedyn``, has 1189 one-cell sources on the 192 x 32 exch2 I/O
layout. Its ``packages.conf`` lists ``exf`` and ``rnf``, its ``data.pkg`` has
:code:`useEXF`, :code:`useSEAICE`, :code:`useThSIce` and :code:`useRNF` set to
true, and its ``data.rnf`` is

::

     &RNF_PARM01
      RNF_file = 'runoff_sparse_const.nc',
      RNF_applyAsRain = .TRUE.,
     &

with ``runoffFile = ' '`` in ``data.exf``; the file is made by
``gen_sparse.py cs32const`` in the same directory as above.

Both committed cases select the legacy route, because they are surface-only
files that reproduce a dense reference. For the default, generalized
route, the build needs ``#define ALLOW_ADDFLUID`` in ``CPP_OPTIONS.h`` (the
``code/`` directory of ``global_ocean.cs32x15`` defines it), ``data`` needs
:code:`selectAddFluid = 1` (or 2) with a blank :code:`addMassFile`, and
``data.rnf`` needs only :code:`RNF_file`. No committed verification input
runs the generalized route; it is exercised by the scratch inputs of the
checks listed in :numref:`ssub_phys_pkg_rnf_verification`.

.. _ssub_phys_pkg_rnf_tools:

Python tools for the input file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The tools are the module :code:`MITgcmutils.runoff` of
:filelink:`utils/python/MITgcmutils` (see :numref:`MITgcmutils`). They need
``netCDF4`` (install with ``pip install MITgcmutils[runoff]``) and are run as
modules, for example from the directory ``utils/python/MITgcmutils``; each
has ``--help``. A grid directory (``--grid-dir``) is the output of a run
of the model with the grid files :code:`hFacC`, :code:`RAC`, :code:`XC`,
:code:`YC` (and, for the target builder, :code:`XG`, :code:`YG`).

**Converter** (``python -m MITgcmutils.runoff.convert DENSE -o OUT --grid-dir DIR``).
Converts a dense exf :code:`runoffFile` (m/s) into a sparse file, and back
with ``--to-dense``. Grid-cell areas turn m/s into m\ :sup:`3`/s, each wet
cell with runoff becomes a source of its own unless a grouping file
(``--sources``) puts cells together, the file written is schema 1.0, and the exf timing settings are given with
``--period``, ``--startdate1``, ``--startdate2``, ``--start-time``,
``--repeat-cycle``, ``--calendar``, ``--clim-year`` and ``--years``. A dense
:code:`runoftempfile` is passed with ``--temperature``. The result is checked
with the checker unless ``--no-check`` is given.

**Checker** (``python -m MITgcmutils.runoff.check FILE [FILE ...]``). Validates one file, or
a set of yearly files given in time order, against the schema and returns
exit code 0 (no errors), 1 (errors, or warnings with ``--strict``) or 2 (usage
or I/O problem). With ``--grid-dir`` it also compares the target cells with the
grid (wet cell, area, centre); ``--surface-level {top,bottom}`` says which
level is the surface cell, ``--tables-only`` checks a file without time series
and ``--json OUT`` writes the report. The Python functions are
:code:`check_files` and, for a small valid example file, :code:`write_example`.

**Target-table builder** (``python -m MITgcmutils.runoff.targets SOURCES --grid-dir DIR -o OUT``).
Builds the source, alias and target tables from a table of source locations
(CSV, or a NetCDF file) and the model grid: each source is put in the nearest wet cell
(``--emission pointwise``) or spread over cells with a gaussian, exponential
or linear kernel (``--emission spread``, ``--spread-type``, ``--spread-scale``,
``--cutoff``). ``--connectivity {latlon,exch2}`` says how cells are connected
for spreading, and ``--max-snap-distance`` bounds the distance from a source to
its cell. ``--depth-range TOP:BOTTOM`` spreads each source over the levels in
that depth range and writes schema 1.1 ``target_level`` values; ``--into RUNOFF``
replaces the tables of an existing runoff file. The output of the builder has
the tables only: the time series
(:code:`runoff_flux` and the optional ones) are written by the user, for example
with netCDF4 or xarray as in the schema document.

.. _ssub_phys_pkg_rnf_verification:

Verification
~~~~~~~~~~~~

The list of tests, the oracles and their tolerances are in
``docs/verification_matrix.md`` of the development repository, and the design
and the reasons for each choice in ``docs/package_design.md``. In summary:

- **No change.** With :code:`pkg/rnf` compiled in and :code:`useRNF` false, and
  with it compiled out, the results of the affected verification experiments
  are unchanged.
- **Sparse equals dense.** The runoff of dense :code:`runoffFile` cases is
  converted to sparse form and must reproduce the dense result. The cases are
  ``lab_sea/input.rnof_sp_const`` (constant record, with a source spanning a tile
  and a process boundary) and ``global_ocean.cs32x15/input.rnof_sp_icedyn`` (cubed
  sphere, 12 tiles). Both run on the legacy route, since a dense reference
  exists only for surface runoff. The five timed ``lab_sea`` cases
  (``input.rnof_sp_daily``, ``_month``, ``_month1``, ``_clim``, ``_yearly``) have
  their own reference outputs. As the digit threshold of the comparison is not sensitive
  enough to timing, a separate check compares the applied runoff field of the sparse run with
  that of the dense run, cell by cell, at every step.
- **Direct oracles.** Python checks in ``tests/rnf`` of the development
  repository compare the model with values computed from the file: the
  applied field and the volume conservation at every record, the heat, salt and
  tracer terms against analytic values, budget closure of volume, heat, salt
  and tracers over the domain, the diagnostics and monitor against a
  reconstruction, and, for the generalized route, the placement at each level
  (``tests/rnf/level_check.py``), including the level independence of the
  first-step free surface.
- **Refusals.** Each configuration and file refusal described above is a case
  of ``tests/rnf/refusal_check.py``, on one and on two processes, with controls
  that must run.
- **Python tools.** The converter, checker and target builder have their own unit tests
  (``tests/runoff``), including the round trip dense to sparse to dense.

Limits that are stated in the matrix and the design:

- No LLC experiment exists in ``verification/``; the placement code for a compact
  layout wider than the global array is compiled and has not been run, and
  no test has a blank exch2 tile.
- The generalized route has no dense counterpart. Its oracles are analytic and
  budget-based, on :filelink:`lab_sea <verification/lab_sea>` with scratch inputs
  (``ALLOW_ADDFLUID``, linear free surface). The loss of runoff for
  :filelink:`pkg/kpp`, :filelink:`pkg/seaice` and :filelink:`pkg/thsice` on that
  route is a documented consequence and has no test.
- The adjoint is not tested. The scale of a global 2 km daily case
  (10\ :sup:`5` to 10\ :sup:`6` sources) has not been measured.
- Restart behavior of the package is not separately tested, and the refusals that
  need an open boundary or an ice shelf are not run, since no experiment that
  compiles the package has either.
- :code:`RNF_cellVolMax` and :code:`RNF_srcFluxMax` refuse implausible runoff
  and do not certify a plausible one.
