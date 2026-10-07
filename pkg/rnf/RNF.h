#ifdef ALLOW_RNF

CBOP
C !ROUTINE: RNF.h

C !DESCRIPTION: \bv
C     *==========================================================*
C     | RNF.h
C     | o Header of the sparse runoff (RNF) package.
C     |   Holds the run-time parameters read from data.rnf.
C     *==========================================================*
C \ev

C-----------------------------------------------------------------------
C--   Parameters that can be set in data.rnf, namelist /RNF_PARM01/
C     RNF_file           :: name of the sparse runoff NetCDF file, or
C                           its base name with RNF_useYearlyFiles, in
C                           which case the reader reads
C                           <base>_YYYY.nc (a trailing ".nc" of the
C                           name given is replaced, so both
C                           "runoff.nc" and "runoff" select
C                           runoff_YYYY.nc); a blank name stops the
C                           run when useRNF=T
C     RNF_holdRecord     :: F: interpolate linearly in time between
C                           records, as pkg/exf does;
C                           T: hold each record over its interval.
C                           The interval is the record's own time up
C                           to the next record's time for a fixed
C                           period (so a yearly-file record at
C                           1 January 00:00, which is the start of
C                           its bounds, is held over that day), and
C                           the calendar month for the monthly modes,
C                           which is NOT the nearer of two mid-month
C                           times (package design, decision 7)
C     RNF_startDate1     :: start date of the series (YYYYMMDD);
C                           overrides the file, 0 = not set
C     RNF_startDate2     :: start date of the series (HHMMSS);
C                           overrides the file, 0 = not set
C     RNF_startTime      :: start time of the series (s); overrides
C                           the file, UNSET_RL = use the file
C     RNF_period         :: record period; overrides the file,
C                           UNSET_RL = use the file;
C                           0 = constant, > 0 = seconds,
C                           -12 = monthly climatology, -1 = monthly
C     RNF_repCycle       :: repeat cycle (s); overrides the file,
C                           UNSET_RL = use the file
C     RNF_useYearlyFiles :: append _YYYY to RNF_file by model year
C     RNF_useTemp        :: apply runoff_temperature if the file has it
C     RNF_useSalt        :: apply runoff_salinity if the file has it
C     RNF_usePtracers    :: apply runoff_ptracer_* if the file has them
C     RNF_monFreq        :: monitor interval (s), default monitorFreq
C     RNF_debugLev       :: message level, default debugLevel
C
C--   Internal parameters (not in the namelist)
C     RNFisON            :: package activation status
C
C--   Constants fixed by the model contract (not in the namelist)
C     RNF_fracTol        :: tolerance on the sum over the domain of
C                           the fractions of one source (must be 1)
C     RNF_areaTol        :: relative tolerance between the cell area
C                           stored in the file and rA
C     RNF_lonLatTol      :: tolerance between the cell centre stored
C                           in the file (target_lon, target_lat) and
C                           the owning cell's XC,YC, as a fraction
C                           of that cell's own spacing
C                           MIN(dxF,dyF). 0.5 means "the stored
C                           centre has to lie inside the owning
C                           cell".
C                           The margin of a corrupted target_cell is
C                           the distance between the two centres
C                           over this threshold. It is NOT 2 in
C                           general: 1 + spacing(from)/spacing(to)
C                           holds only where the spacing is locally
C                           uniform, and at the cs32 facet corner
C                           dxF jumps from 120208 to 156359 m, where
C                           the margin is 1.78. What does hold is
C                           that the margin is strictly above 1 on
C                           any grid with positive cell sizes: the
C                           centres of two distinct cells are about
C                           0.5*(s_from + s_to) apart along the
C                           move, and MIN(dxF,dyF) of the
C                           destination is at most s_to, so the
C                           ratio is about 1 + s_from/s_to.
C                           "About", not "at least": for a zonal
C                           move the great-circle distance is a
C                           little SHORTER than the along-parallel
C                           spacing (0.99995 of it at 77N with a
C                           2 degree step), which is exactly why
C                           lab_sea below measures 1.999904 rather
C                           than 2. The
C                           measured floor over EVERY ordered pair
C                           of distinct cells is 1.779673 on cs32
C                           (37,742,592 pairs) and 1.999904 on
C                           lab_sea, so no corruption of a single
C                           target_cell - by one cell or by any
C                           other amount - can evade the check on
C                           either grid. Those two numbers do not
C                           depend on rSphere: the margin is a
C                           ratio of two lengths that both scale
C                           with it, so a consistent change cancels
C                           (6370 and 6371 km both give 1.779673 on
C                           cs32). Only a MISMATCH moves it, by the
C                           ratio of the two radii, so whoever
C                           re-measures this has to take the
C                           coordinates and the spacing from the
C                           same grid dump.
C                           MIN, not MAX, of the two spacings is
C                           necessary rather than merely cautious:
C                           under MAX the margin of a one-cell zonal
C                           move on a 1-degree global lat-lon grid
C                           is 2*cos(lat), i.e. 0.347 at 80N and
C                           0.035 at 89N, so MAX would defeat the
C                           guard on every high-latitude row.
C                           Because both sides scale with the local
C                           cell, the one constant needs no
C                           retuning between grids or resolutions
C     RNF_fluxMax        :: a flux above this value in absolute terms
C                           counts as a missing value (not allowed)
C     RNF_srcFluxMax     :: largest volume flux one source may carry
C                           [m^3/s]. A record with a larger value
C                           stops the run, naming the source
C                           (RNF_NC_READ_ONE). This is the package's
C                           replacement for the pkg/exf runoff upper
C                           bound of 1.E-6 m/s, which EXF_CHECK_RANGE
C                           skips when useRNF is true because a sparse
C                           point source exceeds it by construction:
C                           1000 m^3/s into one 2 km cell is
C                           2.5E-4 m/s.
C                           The value, 1.E7 m^3/s, is set from the
C                           physics and not from taste: the Amazon,
C                           the largest river on Earth, carries about
C                           2.1E5 m^3/s and all rivers of the world
C                           together about 1.2E6 m^3/s, so one source
C                           id may hold 48 Amazons, or every river on
C                           Earth with a factor of 8 to spare. What it
C                           refuses is a unit mistake: a flux given
C                           per year rather than per second is
C                           3.2E7 times too large, so anything above
C                           0.32 m^3/s is caught; a factor of 1000
C                           (mm, or kg/s read as m^3/s) is caught
C                           above 1E4 m^3/s, i.e. for any source the
C                           size of a real river.
C                           What the applied field may then reach
C                           depends on the grid, because the package
C                           applies flux*frac/rA with frac in [0,1]:
C                           1.E7 m^3/s is 3.2E-4 m/s into the lab_sea
C                           target cell (rA = 3.112287E10 m^2) and
C                           2.5 m/s into a 2 km cell (4.0E6 m^2). A
C                           per-cell value of 1E9 m/s, the kind of
C                           number a unit error produces, needs a cell
C                           smaller than 1E-2 m^2 to get past this
C                           bound.
C                           READ THIS AS A FILE-SCALE UNIT-ERROR
C                           FILTER, NOT AS A PER-CELL SAFETY BOUND.
C                           It does not see the cell: four sources
C                           each carrying exactly this value, with
C                           every target on one lab_sea cell, apply
C                           1.285228E-3 m/s - 1285 times the pkg/exf
C                           bound that was relaxed - and the run ends
C                           normally with no EXF WARNING at all
C                           (measured). N is 10^5 to 10^6 sources in
C                           the intended global 2 km case, so the
C                           aggregate headroom is five to six orders
C                           of magnitude, and a converter index bug
C                           that collapses sources onto one cell
C                           reaches it. At 2 km this bound admits
C                           2.5 m/s per cell while a real Amazon is
C                           5.25E-2 m/s, so it constrains the file and
C                           not the applied field. A grid of smaller
C                           cells admits proportionally more. Closing
C                           that gap needs rA, the top-layer
C                           thickness and deltaT, i.e. a different
C                           check, not a different number here.
C                           The per-cell SIGN is still guarded, by the
C                           negative-runoff test of EXF_CHECK_RANGE,
C                           which useRNF does not skip - but that test
C                           runs at nIter0 only, as the bound it sits
C                           beside did before this change, so the
C                           sign asymmetry is pre-existing and not
C                           introduced here; what gained coverage was
C                           the magnitude - per record here, and per
C                           cell per step in RNF_cellVolMax below,
C                           which is what bounds the per-CELL
C                           aggregate this bound cannot see
C     RNF_cellVolMax     :: largest share of a target cell's top-layer
C                           volume that ONE time step of runoff may
C                           add [1]. RNF_EXF_RUNOFF refuses a cell at
C                           which
C                             |RNF_vflx|*deltaTFreeSurf
C                               > RNF_cellVolMax*drF(ks)*hFacC(..,ks)
C                           naming the cell (i,j,bi,bj), the applied
C                           value and the limit in m/s, on every step.
C                           This is the per-cell companion of
C                           RNF_srcFluxMax, which is per source and
C                           per file and so cannot see an aggregate at
C                           all: four sources each carrying exactly
C                           RNF_srcFluxMax with every target on one
C                           lab_sea cell apply 1.285228E-3 m/s and the
C                           run ended normally, with no message of any
C                           kind, before this check existed (measured
C                           on RUNOFF-030 by review B, re-measured on
C                           RUNOFF-040 against the committed build:
C                           exit 0, "Execution ended Normally", 0 EXF
C                           WARNING lines, one non-zero cell, applied
C                           1.2852284E-03 m/s at (i,j) = (12,2)).
C                           THE QUANTITY. The bound is on the
C                           dimensionless ratio
C                             f(c) = |runoff(c)|*deltaTFreeSurf
C                                    / ( drF(ks)*hFacC(c,ks) )
C                           i.e. the depth of water one step of runoff
C                           puts on the cell over the thickness of the
C                           cell it goes into. With a real freshwater
C                           flux f is the fractional change of the
C                           top-cell volume in that step: the runoff
C                           reaches etaN through EmPmR and dEtaHdt and
C                           is integrated with deltaTFreeSurf
C                           (model/src/integr_continuity.F:221). With
C                           a linear free surface no volume moves and
C                           f is instead the fractional freshwater
C                           dilution the surface tracer forcing
C                           applies in that step; the model linearises
C                           the exact dilution h/(h+dh) = 1/(1+f) to
C                           1-f, whose relative error is exactly f^2
C                           (the linear term is EmPmR*(salt -
C                           salt_EvPrRn)*mass2rUnit,
C                           model/src/external_forcing_surf.F:310-316).
C                           THAT SECOND READING HOLDS ONLY WHERE
C                           dTtracerLev(ks) = deltaTFreeSurf. The two
C                           are equal in both test experiments, but
C                           they need not be: deltaTFreeSurf defaults
C                           to deltaTMom and NOT to deltaTtracer
C                           (model/src/ini_parms.F:1068, whose own
C                           comment calls that default "inappropriate"
C                           and advises deltaTFreeSurf = deltaTtracer
C                           under asynchronous stepping). cs32 has
C                           deltaTMom = 1200 against
C                           deltaTtracer = 86400 and escapes the trap
C                           only by setting deltaTFreeSurf = 86400
C                           explicitly; on that ratio an
C                           asynchronously stepped set-up that left
C                           the default would make the dilution
C                           reading wrong by 72x while the volume
C                           reading stayed right. No enrolled case can
C                           see a mismatch. Whether to bound with
C                           MAX(deltaTFreeSurf,dTtracerLev(ks)) is an
C                           open design question, deliberately not
C                           decided here.
C                           Being dimensionless is the point: ONE
C                           number serves every grid, resolution and
C                           time step, which is what the pkg/exf
C                           runoff bound of 1.E-6 m/s could not do and
C                           why it had to be skipped.
C                           WHERE 0.2 COMES FROM. Two legs, both
C                           properties of f itself. READ THIS BEFORE
C                           RAISING THE CONSTANT.
C                           o CFL. The water injected into the cell
C                             has to leave it: f is exactly the
C                             Courant number of the top-layer outflow
C                             the injection requires, since
C                             f <= 0.2 is U*dt/dx <= 0.2 for the
C                             horizontal outflow U that carries the
C                             added volume away. 0.2 is a standard
C                             advective-CFL safety factor, i.e. a
C                             fifth of the stability limit.
C                           o Linearisation. The surface tracer
C                             forcing is first order in f and its
C                             relative error is exactly f^2 (above),
C                             so 0.2 is the share at which that error
C                             is 4%. Beyond it the model's own
C                             dilution term is no longer a
C                             linearisation of anything.
C                           Both legs say "a deliberate share, not an
C                           edge". For scale, MITgcm's own band on the
C                           surface-cell fraction is
C                           hFacInf = 0.2 to hFacSup = 2.0 (defaults
C                           at model/src/set_defaults.F:258-259,
C                           described at model/inc/PARAMS.h:762 as
C                           "Threshold (inf and sup) for fraction size
C                           of surface cell"). NOTE WHAT THAT BAND IS
C                           AND IS NOT: it bounds the FRACTION, not
C                           its per-step change, and runoff THICKENS
C                           the surface cell, so from a full cell
C                           (hFacC = 1.0, measured at every target of
C                           both test grids) the band is first crossed
C                           at +1.0, at hFacSup. f = 0.2 crosses
C                           nothing; it is 4 to 5 times inside the
C                           band, which is the margin this constant
C                           buys and the reason it is a share rather
C                           than a threshold. (An earlier version of
C                           this comment read hFacInf as a bound on
C                           the per-step change and called 0.2 the
C                           edge of the band. That was wrong in both
C                           parts and both reviewers caught it; the
C                           value 0.2 is unchanged, only its
C                           justification is.)
C                           Outside the band the model does NOT merely
C                           warn, and the two routines differ:
C                           CALC_R_STAR warns and then STOPS on the
C                           thin side (model/src/calc_r_star.F:
C                           201-242), while CALC_SURF_DR's thin-side
C                           STOP is commented out
C                           (model/src/calc_surf_dr.F:105-108) and it
C                           clamps the surface to Rmin_surf instead.
C                           The value is a fixed constant and not a
C                           data.rnf parameter, as RNF_srcFluxMax is:
C                           the dimensionless form needs no retuning
C                           per grid, so a run-time scalar would be
C                           cost without benefit.
C                           WHAT IT IS ON EACH GRID (measured). THE
C                           LIMIT IS NOT A CONSTANT OF THE GRID: it
C                           tracks the live hFacC, so on an r* grid it
C                           moves with the state (see "the thickness
C                           is the live one" below). Each figure below
C                           says which basis it is on.
C                           o lab_sea (rA = 3.112287E10 m^2,
C                             drF(1) = 10 m, deltaTFreeSurf = 3600 s):
C                             5.5556E-4 m/s, i.e. 1.729E7 m^3/s into
C                             that cell, 82 Amazons. Here reference
C                             and live agree exactly and for all time:
C                             lab_sea is a LINEAR free surface
C                             (nonlinFreeSurf = 0, select_rStar = 0,
C                             reported by its own run), so nothing
C                             updates hFacC after initialisation;
C                           o cs32 (drF(1) = 50 m, deltaTFreeSurf =
C                             86400 s) is an r* grid
C                             (nonlinFreeSurf = 4, select_rStar = 2),
C                             so two figures are needed.
C                             REFERENCE basis: 1.1574E-4 m/s at every
C                             one of the 1189 target cells of
C                             input.rnof_sp_icedyn, h0FacC being 1.0 at
C                             all of them in the init dump (its
C                             hFacMinDr of 20 m would allow a thinner
C                             surface cell, but no target has one); in
C                             m^3/s, 1.62E6 (7.7 Amazons) on the
C                             smallest target cell, rA = 1.4019E10 m^2,
C                             and 1.03E7 (49 Amazons) on the median,
C                             rA = 8.8743E10 m^2.
C                             LIVE basis, which is what is actually
C                             enforced and applies from nIter0 because
C                             INITIALISE_VARIA updates r* before the
C                             first step (initialise_varia.F:302,307):
C                             rStarFacC over those targets spans
C                             0.8787 to 0.99956, the thickness 43.94 to
C                             49.98 m and the enforced limit 1.0171E-4
C                             to 1.1569E-4 m/s. 309 of the 1189
C                             targets (26%) sit more than 1% below the
C                             reference figure and 18 (1.5%) more than
C                             10% below. Measured from the committed
C                             run's own Depth.data and
C                             Eta.0000036010.data with CALC_R_STAR's
C                             own formula (calc_r_star.F:103-105), and
C                             agreeing with review B's independent
C                             measurement and with review A's executed
C                             cs32 refusal, which printed
C                             thickness 4.69852227E+01 m and
C                             limit 1.08762090E-04;
C                           o a 2 km cell (rA = 4E6 m^2, drF(1) =
C                             10 m, deltaTFreeSurf = 1200 s):
C                             1.6667E-3 m/s, i.e. 6.67E3 m^3/s, 0.032
C                             Amazons.
C                           A PHYSICALLY CORRECT LARGE RIVER, for
C                           comparison: the Amazon's 2.1E5 m^3/s is
C                           f = 2.43E-3 on the lab_sea cell (82 times
C                           under the bound) and f = 4.09E-3 on the
C                           median cs32 cell (49 times under), but
C                           f = 6.3 in ONE
C                           2 km cell at deltaTFreeSurf = 1200 s - 31
C                           times OVER. That refusal is correct and
C                           not a false positive: 6.3 top-layer
C                           volumes in one step is past hFacSup in the
C                           first step and is not a configuration the
C                           model can integrate. What it says is that
C                           a 2 km grid must spread the Amazon over at
C                           least 32 cells, which its ~200 km mouth is
C                           (about 100 cells) - and spreading a source
C                           over its real cells is what target_fraction
C                           exists for. Every committed sparse oracle
C                           is far below: the largest per-cell f over
C                           every record of every one of them is
C                           6.72E-4, on cs32, a margin of 297 - but
C                           that pair is on the REFERENCE basis. On
C                           the live thickness the same cs32 cell
C                           (5903) is f = 7.16E-4 and the margin 279.5
C                           (measured, and the figure review A
C                           measured at the decisive cell at nIter0).
C                           The lab_sea files reach 3.42E-4, a margin
C                           of 585, on both bases at once, that grid
C                           being a linear free surface.
C                           WHAT IT DOES NOT COVER.
C                           o It is per cell and per STEP, so it does
C                             not keep a run physical: a flux just
C                             under the bound, sustained, still adds
C                             0.2 of the surface layer every step. It
C                             refuses the absurd, it does not certify
C                             the plausible.
C                           o It sees the sparse field only. Under
C                             ALLOW_CTRL with ALLOW_GENTIM2D_CONTROL,
C                             xx_runoff is added to the exf runoff
C                             array at pkg/exf/exf_getffields.F:531-534,
C                             AFTER the RNF_EXF_RUNOFF call at :456, so
C                             neither this bound nor RNF_srcFluxMax
C                             sees the controlled field. EXF_CHECK_RANGE
C                             does run after that addition
C                             (exf_getforcing.F:199 then :348), but of
C                             its tests on the runoff array the upper
C                             bound is skipped with useRNF and the
C                             sflux one has the runoff added back. That
C                             leaves the negative test - the SIGN, at
C                             nIter0 - and, only where ALLOW_RUNOFTEMP
C                             is compiled (cs32 defines it, lab_sea
C                             does not), a 36 m/s ceiling that the
C                             runoff-TEMPERATURE test reads from the
C                             runoff array where it means runoftemp
C                             (exf_check_range.F, the ALLOW_RUNOFTEMP
C                             block): an upstream misnaming, not
C                             conditioned on useRNF, and 6.5E4 times
C                             above this bound on the lab_sea cell, so
C                             it constrains nothing in practice.
C                             Nothing therefore bounds the magnitude of
C                             xx_runoff; closing that is a pkg/ctrl
C                             question and not this bound's.
C                           o It bounds the magnitude only: the
C                             temperature, the salinity and the tracer
C                             concentrations the water carries are not
C                             bounded by it, and neither is the sign
C                             (see RNF_srcFluxMax above).
C                           o It cannot tell one wrong source from N
C                             collapsed ones. It names the CELL, which
C                             is what locates a collapsed target; the
C                             sources feeding it are in the file's
C                             target table.
C                           o It is blind to a collapse onto a cell
C                             whose top layer is thick and whose time
C                             step is short, exactly in proportion to
C                             drF(ks)*hFacC/deltaTFreeSurf.
C                           o The NaN arm of the test (it is written
C                             .NOT.(x .LE. lim) so a value that is not
C                             a number is refused rather than compared)
C                             is correct but UNREACHABLE on a supported
C                             input: RNF_NC_READ_FLUX refuses a
C                             non-finite flux first and rA is positive,
C                             so nothing can deliver a NaN here. It is
C                             defence in depth, and its validity rests
C                             on the optfile: -ffinite-math-only would
C                             silently void it AND the pre-existing
C                             rnf_init_fixed.F:728-737 tests of the
C                             same shape. The build measured for this
C                             issue is -O0 with no fast-math.
C                           o Granularity: neither enrolled case
C                             separates 0.2 from any value in
C                             (0.99*limit, limit], so a mutant that
C                             tightened the bound by less than 1% would
C                             pass both.
C                           o The GLOBAL_SUM_INT is unconditional for
C                             every useRNF run, including one with no
C                             target on any tile, and is unmeasured
C                             beyond 2 processes.
C                           o RNF_tgtK is fixed at init while this
C                             guard re-evaluates kSurfC every step, so
C                             the two could diverge under pkg/shelfice
C                             remeshing. Not reachable now: a shelfice
C                             target is refused at init.
C                           THE THICKNESS IS THE LIVE ONE, which is a
C                           strength and not a caveat, but it has two
C                           consequences worth stating. _hFacC resolves
C                           to hFacC (model/inc/HFACC_MACROS.h:37-39,
C                           the macro also adapting to the reduced-
C                           memory HFACC_* options), and the only
C                           run-time writer of hFacC is
C                           model/src/update_r_star.F:55-57, which sets
C                           hFacC = h0FacC*rStarFacC. So:
C                           o With select_rStar the guard enforces the
C                             r*-stretched thickness of the current
C                             state - state-consistent, with no
C                             rStarFacC factor left to apply. But that
C                             limit is NOT clipped: calc_r_star.F:
C                             185-198 only COUNTS cells outside
C                             [hFacInf,hFacSup], so the limit can drift
C                             with the state and ONE FILE CAN PASS AT
C                             nIter0 AND BE REFUSED LATER. That
C                             mid-run abort is deliberate: an
C                             init-only check would be unsound in the
C                             numerator (the flux series is not
C                             static) and in the denominator (the
C                             thickness is not either), and bounding
C                             by h0FacC*hFacInf instead would be 5
C                             times stricter than the physics above and
C                             would refuse legitimate configurations
C                             at init.
C                           o With nonlinFreeSurf and select_rStar = 0
C                             the opposite holds: CALC_SURF_DR writes
C                             hFac_surfC and NOT hFacC
C                             (model/src/calc_surf_dr.F:120-122), so
C                             there the guard uses the REFERENCE
C                             thickness and UNDER-states the actual
C                             departure. No test experiment runs that
C                             combination, so it is unmeasured here.
C                           o Under a linear free surface nothing
C                             updates hFacC at all and live equals
C                             reference for the whole run. That is
C                             lab_sea, and it is the premise the
C                             0.99-of-the-bound control case relies on.
C     RNF_idLen          :: length of a source id in the model
C     RNF_maxErrMsg      :: messages one process prints per error
C                           counter, so that a large file cannot fill
C                           the log
C
C--   What the static read of RNF_file found (RNF_INIT_FIXED)
C     RNF_nSrcFile       :: number of sources in the file
C     RNF_nTgtFile       :: number of target entries in the file
C     RNF_nRecFile       :: number of time records in the file
C     RNF_nSrcProc       :: tile-source entries of this process
C                           (a source on two tiles counts twice)
C     RNF_nTgtProc       :: target entries owned by this process
C     RNF_nTgtOwned      :: target entries owned over all processes
C
C--   Per-tile lists built by RNF_INIT_FIXED. Index k runs over the
C     entries of tile (bi,bj); only what the tile owns is stored.
C     RNF_nSrc           :: number of sources on the tile
C     RNF_nTgt           :: number of target entries on the tile
C     RNF_srcGlob        :: 0-based source index in the file of local
C                           source k
C     RNF_tgtI, RNF_tgtJ :: local cell indices of target entry k
C     RNF_tgtK           :: model level of target entry k (the
C                           surface level of its column)
C     RNF_tgtSrc         :: local source index of target entry k
C     RNF_tgtFrac        :: fraction of that source sent to the cell
C     RNF_srcFracSum     :: fractions of local source k summed over
C                           the entries of this tile (fraction check)
C     RNF_srcFlux        :: volume flux of local source k at the
C                           current model time [m^3/s], combined from
C                           the two record buffers by RNF_FIELDS_LOAD
C     RNF_srcTemp        :: its temperature [degC], idem
C     RNF_srcTvld        :: 1 when that temperature is present in
C                           every record used, 0 otherwise; the
C                           source then contributes to neither (mT)
C                           nor m_T and enters at the reference
C                           temperature (package design, decision 3,
C                           "Missing temperature")
C     RNF_srcSalt        :: its salinity [g/kg], idem
C     RNF_srcTrc         :: its concentration of runoff tracer n, idem
C
C--   Dense per-tile fields, rebuilt at every time step
C     RNF_vflx           :: runoff volume flux per unit area [m/s],
C                           sum_s flux_s*frac_s,c / rA; this is what
C                           RNF_EXF_RUNOFF puts in the exf field
C     RNF_mflx           :: runoff mass flux [kg/m^2/s],
C                           rhoConstFresh*RNF_vflx
C     RNF_mflxT          :: the same sum restricted to the sources
C                           whose temperature is present, m_T
C                           [kg/m^2/s] (package design, decision 3)
C     RNF_mXT            :: (mT) = sum_s m_s*T_s over the sources with
C                           a temperature [kg/m^2/s * degC]
C     RNF_mXS            :: (mS) = sum_s m_s*S_s over all sources
C                           [kg/m^2/s * g/kg]
C     RNF_mXTr           :: (mC_n) = sum_s m_s*C_s,n for runoff tracer
C                           n [kg/m^2/s * tracer units]
C
C--   The same five fields at the time level the tendency terms use
C     (RNF_ap*). They are a copy of the fields above, except in the
C     one case where the model's own freshwater flux lags: see
C     RNF_lagFlds below.
C
C--   Work space of the global sums (filled by the master thread,
C     read by every thread inside GLOBAL_SUM_VECTOR_RL)
C     RNF_sumWrk         :: per-tile fraction sums of one chunk of
C                           sources
C     RNF_sumVec         :: their sum over tiles and processes
C     RNF_fluxBuf        :: NetCDF read buffer of one chunk of fluxes
C     RNF_tileSum        :: one value per tile, for GLOBAL_SUM_TILE_RL
C     RNF_vldBuf         :: "the value was present" flag of each
C                           entry of RNF_fluxBuf
C     RNF_vldWrk         :: per-source "the value was present" flags
C                           of the series RNF_NC_READ_ONE just read.
C                           Only the temperature keeps them (in
C                           RNF_bufTvld); every other series refuses a
C                           missing value, so one work array serves
C                           all of them
C-----------------------------------------------------------------------
CEOP

      LOGICAL RNFisON
      LOGICAL RNF_holdRecord
      LOGICAL RNF_useYearlyFiles
      LOGICAL RNF_useTemp
      LOGICAL RNF_useSalt
      LOGICAL RNF_usePtracers
      COMMON /RNF_PARM_L/
     &     RNFisON,
     &     RNF_holdRecord, RNF_useYearlyFiles,
     &     RNF_useTemp, RNF_useSalt, RNF_usePtracers

      INTEGER RNF_startDate1
      INTEGER RNF_startDate2
      INTEGER RNF_debugLev
      COMMON /RNF_PARM_I/
     &     RNF_startDate1, RNF_startDate2,
     &     RNF_debugLev

      _RL RNF_startTime
      _RL RNF_period
      _RL RNF_repCycle
      _RL RNF_monFreq
      COMMON /RNF_PARM_R/
     &     RNF_startTime, RNF_period, RNF_repCycle,
     &     RNF_monFreq

      CHARACTER*(MAX_LEN_FNAM) RNF_file
      COMMON /RNF_PARM_C/
     &     RNF_file

      _RL RNF_fracTol
      PARAMETER ( RNF_fracTol = 1. _d -6 )
      _RL RNF_areaTol
      PARAMETER ( RNF_areaTol = 1. _d -4 )
      _RL RNF_lonLatTol
      PARAMETER ( RNF_lonLatTol = 0.5 _d 0 )
      _RL RNF_fluxMax
      PARAMETER ( RNF_fluxMax = 1. _d 30 )
      _RL RNF_srcFluxMax
      PARAMETER ( RNF_srcFluxMax = 1. _d 7 )
      _RL RNF_cellVolMax
      PARAMETER ( RNF_cellVolMax = 0.2 _d 0 )
      INTEGER RNF_idLen
      PARAMETER ( RNF_idLen = 64 )
      INTEGER RNF_maxErrMsg
      PARAMETER ( RNF_maxErrMsg = 20 )

C--   RNF_tableRead :: the target table was read and every entry of
C     it was valid, so the fraction sums are worth computing
C     RNF_lonLatChk :: the cell-centre check of RNF_INIT_FIXED ran.
C     It is optional in both directions: target_lon/target_lat are
C     optional in schema 1.0, and XC,YC are degrees only on a
C     spherical-polar or curvilinear grid. RNF_SUMMARY reports it so
C     that a run never looks guarded when it was not.
      LOGICAL RNF_tableRead
      LOGICAL RNF_lonLatChk
      COMMON /RNF_COUNT_L/ RNF_tableRead, RNF_lonLatChk

C--   Effective time handling (RNF_TIME_SETUP). These are the values
C     RNF_GETREC hands to the pkg/exf record-selection routines, so
C     the sparse path picks records with the code the dense path runs.
C     They come from the file's own time axis (mitgcm_time_sampling,
C     mitgcm_time_period, mitgcm_time_repeat, time:units,
C     time:calendar, time, time_bnds) and are overridden by data.rnf.
C     RNF_recPeriod :: period in exf's convention: 0 constant, > 0
C                      seconds, -12 monthly climatology (12 records,
C                      January to December, repeated every model
C                      year), -1 consecutive calendar months
C     RNF_recStart  :: time of record 1 [s], exf's fldStartTime: in
C                      model time, or the offset from 1 January of
C                      its own year when RNF_useYearlyFiles
C     RNF_recCycle  :: repeat cycle [s], 0 = no repeat. For a fixed
C                      period with mitgcm_time_repeat = "annual" it
C                      is the span of time_bnds, which is what makes
C                      the cycle repeat on the file's real dates
C                      rather than on a nominal calendar year
C     RNF_recDate1  :: date of record 1 (YYYYMMDD) that the reader
C     RNF_recDate2  :: derived from the file or was given (HHMMSS);
C                      0 when no date was needed (constant sampling,
C                      or a run without pkg/cal)
      _RL RNF_recPeriod
      _RL RNF_recStart
      _RL RNF_recCycle
      COMMON /RNF_TIME_R/
     &     RNF_recPeriod, RNF_recStart, RNF_recCycle

      INTEGER RNF_recDate1
      INTEGER RNF_recDate2
      COMMON /RNF_TIME_I/
     &     RNF_recDate1, RNF_recDate2

C--   The two record buffers of RNF_FIELDS_LOAD. Each holds the
C     per-tile source fluxes of one time record, tagged with the
C     record number and the file year it came from, so that a buffer
C     is re-read only when the bracket moves off it. A tag of -1 means
C     "nothing read yet"; RNF_INIT_VARIA sets both, which is why a
C     restart needs no pickup (package design, decision 8).
C     RNF_bufRec    :: record number held in each buffer (-1: none)
C     RNF_bufYr     :: file year it came from (0 without yearly files)
C     RNF_bufFlux   :: its per-tile source fluxes [m^3/s]
C     RNF_bufSum    :: its flux summed over the sources of the file
C     RNF_fluxFile  :: flux summed over the sources of the file at the
C                      current model time, combined from RNF_bufSum
C                      with the same weights as RNF_srcFlux
      INTEGER RNF_bufRec(2)
      INTEGER RNF_bufYr(2)
      COMMON /RNF_BUF_I/
     &     RNF_bufRec, RNF_bufYr

C     RNF_bufTemp   :: its per-tile source temperatures [degC]
C     RNF_bufTvld   :: 1 where that temperature is present, 0 where it
C                      is missing (a missing temperature is allowed:
C                      runoff schema 1.0, section 3.5, and the source
C                      then enters at the reference temperature). The
C                      stored temperature of a missing value is 0, so
C                      that a fill value cannot propagate as a NaN
C     RNF_bufSalt   :: its per-tile source salinities [g/kg]
C     RNF_bufTrc    :: its per-tile source tracer concentrations
      _RL RNF_bufFlux(RNF_nSrcTile,nSx,nSy,2)
      _RL RNF_bufTemp(RNF_nSrcTile,nSx,nSy,2)
      _RL RNF_bufTvld(RNF_nSrcTile,nSx,nSy,2)
      _RL RNF_bufSalt(RNF_nSrcTile,nSx,nSy,2)
      _RL RNF_bufTrc (RNF_nSrcTile,nSx,nSy,2,RNF_nTr)
      _RL RNF_bufSum(2)
      _RL RNF_fluxFile
      COMMON /RNF_BUF_R/
     &     RNF_bufFlux, RNF_bufTemp, RNF_bufTvld,
     &     RNF_bufSalt, RNF_bufTrc,
     &     RNF_bufSum, RNF_fluxFile

C--   What time series the file carries, found by RNF_NC_SERIES from
C     the variables of the file (package design, decisions 3 and 4).
C     A variable the user switched off (RNF_useTemp, RNF_useSalt,
C     RNF_usePtracers) is reported and then treated as absent.
C     RNF_hasTemp  :: the file has runoff_temperature and it is used
C     RNF_hasSalt  :: the file has runoff_salinity and it is used
C     RNF_nTrUse   :: number of runoff_ptracer_<NAME> variables used
C     RNF_trPtr    :: ptracer number (1..PTRACERS_numInUse) that
C                     runoff tracer n feeds. A <NAME> with no matching
C                     PTRACERS_names entry stops the run
C     RNF_trNam    :: the <NAME> of runoff tracer n, i.e. the variable
C                     name without the runoff_ptracer_ prefix
C     RNF_applyT   :: the temperature term can be non-zero, so
C                     RNF_TENDENCY_APPLY_T has work to do
C     RNF_applyS   :: idem for salinity. True also without
C                     runoff_salinity when salt_EvPrRn is not 0,
C                     because the water then still arrives at
C                     salinity 0 while the model gave it salt_EvPrRn
C                     (package design, decision 3, "S = 0")
      LOGICAL RNF_hasTemp
      LOGICAL RNF_hasSalt
      LOGICAL RNF_applyT
      LOGICAL RNF_applyS
      COMMON /RNF_SERIES_L/
     &     RNF_hasTemp, RNF_hasSalt, RNF_applyT, RNF_applyS

      INTEGER RNF_nTrUse
      INTEGER RNF_trPtr(RNF_nTr)
      COMMON /RNF_SERIES_I/
     &     RNF_nTrUse, RNF_trPtr

      CHARACTER*(RNF_idLen) RNF_trNam(RNF_nTr)
      COMMON /RNF_SERIES_C/
     &     RNF_trNam

C--   Time level of the tendency terms (package design, decision 3,
C     "Time level"). The package fields must belong to the same step
C     as the freshwater flux the model uses for its own temperature,
C     salinity and tracer terms. That flux is EmPmR of the current
C     step in branches L and U and PmEpR of the current step in
C     branch N with staggerTimeStep, but PmEpR of the PREVIOUS step
C     in branch N without it (model/src/integr_continuity.F:172-179),
C     and zero at the first step of a run that starts at iteration 0
C     (model/src/ini_nlfs_vars.F:59).
C     The two iteration numbers below are what the time levels are
C     tracked by, rather than the model times themselves, because an
C     iteration is exact: myTime - deltaTClock of one step need not be
C     bitwise the myTime of the step before it.
C     RNF_lagFlds :: the terms use the fields of the previous step,
C                    i.e. branch N without staggerTimeStep
C     RNF_curIter :: iteration of RNF_vflx and the other dense fields
C                    (RNF_noIter: nothing built yet)
C     RNF_apIter  :: iteration of the RNF_ap* fields, the ones the
C                    tendency routines read (RNF_noIter: none yet)
      LOGICAL RNF_lagFlds
      COMMON /RNF_TLEV_L/ RNF_lagFlds

      INTEGER RNF_curIter
      INTEGER RNF_apIter
      COMMON /RNF_TLEV_I/ RNF_curIter, RNF_apIter

      INTEGER RNF_noIter
      PARAMETER ( RNF_noIter = -999999999 )

C--   Counts of the static read
      INTEGER RNF_nSrcFile
      INTEGER RNF_nTgtFile
      INTEGER RNF_nRecFile
      INTEGER RNF_nSrcProc
      INTEGER RNF_nTgtProc
      INTEGER RNF_nTgtOwned
      COMMON /RNF_COUNT_I/
     &     RNF_nSrcFile, RNF_nTgtFile, RNF_nRecFile,
     &     RNF_nSrcProc, RNF_nTgtProc, RNF_nTgtOwned

C--   Per-tile source and target lists
      INTEGER RNF_nSrc   (nSx,nSy)
      INTEGER RNF_nTgt   (nSx,nSy)
      INTEGER RNF_srcGlob(RNF_nSrcTile,nSx,nSy)
      INTEGER RNF_tgtI   (RNF_nTgtTile,nSx,nSy)
      INTEGER RNF_tgtJ   (RNF_nTgtTile,nSx,nSy)
      INTEGER RNF_tgtK   (RNF_nTgtTile,nSx,nSy)
      INTEGER RNF_tgtSrc (RNF_nTgtTile,nSx,nSy)
      COMMON /RNF_LIST_I/
     &     RNF_nSrc, RNF_nTgt, RNF_srcGlob,
     &     RNF_tgtI, RNF_tgtJ, RNF_tgtK, RNF_tgtSrc

      _RL RNF_tgtFrac   (RNF_nTgtTile,nSx,nSy)
      _RL RNF_srcFracSum(RNF_nSrcTile,nSx,nSy)
      _RL RNF_srcFlux   (RNF_nSrcTile,nSx,nSy)
      _RL RNF_srcTemp   (RNF_nSrcTile,nSx,nSy)
      _RL RNF_srcTvld   (RNF_nSrcTile,nSx,nSy)
      _RL RNF_srcSalt   (RNF_nSrcTile,nSx,nSy)
      _RL RNF_srcTrc    (RNF_nSrcTile,nSx,nSy,RNF_nTr)
      COMMON /RNF_LIST_R/
     &     RNF_tgtFrac, RNF_srcFracSum, RNF_srcFlux,
     &     RNF_srcTemp, RNF_srcTvld, RNF_srcSalt, RNF_srcTrc

C--   Dense per-tile fields
      _RL RNF_vflx (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL RNF_mflx (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL RNF_mflxT(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL RNF_mXT  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL RNF_mXS  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL RNF_mXTr (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy,RNF_nTr)
      COMMON /RNF_FIELDS_R/
     &     RNF_vflx, RNF_mflx, RNF_mflxT,
     &     RNF_mXT, RNF_mXS, RNF_mXTr

C--   The same, at the time level of the tendency terms
      _RL RNF_apMflx (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL RNF_apMflxT(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL RNF_apXT   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL RNF_apXS   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL RNF_apXTr  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy,RNF_nTr)
      COMMON /RNF_APPLY_R/
     &     RNF_apMflx, RNF_apMflxT,
     &     RNF_apXT, RNF_apXS, RNF_apXTr

C--   Work space of the global sums and of the NetCDF reads
      _RL RNF_sumWrk(nSx,nSy,RNF_nBuf)
      _RL RNF_sumVec(RNF_nBuf)
      _RL RNF_fluxBuf(RNF_nBuf)
      _RL RNF_tileSum(nSx,nSy)
      _RL RNF_vldBuf(RNF_nBuf)
      _RL RNF_vldWrk(RNF_nSrcTile,nSx,nSy)
      COMMON /RNF_WORK_R/
     &     RNF_sumWrk, RNF_sumVec, RNF_fluxBuf, RNF_tileSum,
     &     RNF_vldBuf, RNF_vldWrk

#endif /* ALLOW_RNF */
