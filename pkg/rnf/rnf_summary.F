#include "RNF_OPTIONS.h"

CBOP
C !ROUTINE: RNF_SUMMARY

C !INTERFACE: ==========================================================
      SUBROUTINE RNF_SUMMARY( myThid )

C !DESCRIPTION:
C     Print the RNF parameters (data.rnf), the constants fixed by the
C     model contract, the array bounds of RNF_SIZE.h and what the
C     static read of RNF_file found, to the standard output of every
C     process, in the layout of EXF_SUMMARY.
C     Called at the end of RNF_INIT_FIXED, so the counts are the
C     final ones. RNF_nSrcProc and RNF_nTgtProc are the counts of
C     this process; RNF_nTgtOwned is summed over all processes.
C     RNF_lonLatChk reports whether the optional check of the target
C     cell centres against XC,YC ran on this file and grid.
C     The RNF_rec* values are the time handling RNF_TIME_SETUP
C     resolved from the file and the data.rnf overrides, i.e. what
C     the pkg/exf record selection is actually called with.
C     RNF_hasTemp, RNF_hasSalt, RNF_nTrUse, RNF_applyT, RNF_applyS
C     and one line per runoff tracer report which time series of the
C     file are applied and which ptracer each tracer feeds, so that a
C     run never looks as if it carried a property it did not; and
C     RNF_lagFlds reports which time level the salt and tracer terms
C     use.
C     RNF_srcFluxMax, RNF_cellVolMax and the five lines after them
C     report which runoff bounds the run was held to: EXF_CHECK_RANGE
C     conditions two of its tests on useRNF, so a log that named
C     neither relaxation nor the two package bounds that replace them
C     would leave a reader unable to tell which bounds applied at all.
C     A sixth, conditional line (RUNOFF-041) follows them only when
C     dTtracerLev(1) differs from deltaTFreeSurf: it names both steps
C     and their ratio and says that RNF_cellVolMax's dilution reading
C     (not its primary volume reading) is looser by that factor. It is
C     silent, not merely absent of content, when the two are equal.

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "RNF_SIZE.h"
#include "RNF.h"

C !INPUT PARAMETERS: ===================================================
C  myThid               :: my thread Id number
      INTEGER myThid

C !OUTPUT PARAMETERS: ==================================================
C  none

#ifdef ALLOW_RNF

C !FUNCTIONS: ==========================================================
      INTEGER  ILNBLNK
      EXTERNAL ILNBLNK

C !LOCAL VARIABLES: ====================================================
C  blkLin               :: blank line
C  msgBuf               :: message buffer
C  dtRatio              :: dTtracerLev(1)/deltaTFreeSurf, when they
C                          differ (RUNOFF-041)
      CHARACTER*1 blkLin
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      CHARACTER*(RNF_idLen) trNam
      _RL dtRatio
      INTEGER k, iL
CEOP

      _BEGIN_MASTER(myThid)

      blkLin = ' '
      CALL PRINT_MESSAGE( blkLin, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)')
     &'// ======================================================='
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)')
     &'// Sparse runoff (RNF) configuration >>> START <<<'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)')
     &'// ======================================================='
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )
      CALL PRINT_MESSAGE( blkLin, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )

C--   Parameters of data.rnf
      CALL WRITE_0D_C( RNF_file, -1, INDEX_NONE, 'RNF_file =',
     &  ' /* sparse runoff NetCDF file */')
      CALL WRITE_0D_L( RNF_holdRecord, INDEX_NONE,
     &  'RNF_holdRecord =',
     &  ' /* hold each record (F: interpolate in time) */')
      CALL WRITE_0D_L( RNF_useYearlyFiles, INDEX_NONE,
     &  'RNF_useYearlyFiles =',
     &  ' /* append _YYYY to RNF_file by model year */')
      CALL WRITE_0D_I( RNF_startDate1, INDEX_NONE,
     &  'RNF_startDate1 =',
     &  ' /* start date override YYYYMMDD (0: not set) */')
      CALL WRITE_0D_I( RNF_startDate2, INDEX_NONE,
     &  'RNF_startDate2 =',
     &  ' /* start date override HHMMSS (0: not set) */')
      CALL WRITE_0D_RL( RNF_startTime, INDEX_NONE,
     &  'RNF_startTime =',
     &  ' /* start time override [s] (UNSET: from file) */')
      CALL WRITE_0D_RL( RNF_period, INDEX_NONE,
     &  'RNF_period =',
     &  ' /* record period override (UNSET: from file) */')
      CALL WRITE_0D_RL( RNF_repCycle, INDEX_NONE,
     &  'RNF_repCycle =',
     &  ' /* repeat cycle override [s] (UNSET: from file) */')
      CALL WRITE_0D_L( RNF_useTemp, INDEX_NONE,
     &  'RNF_useTemp =',
     &  ' /* apply runoff_temperature if present */')
      CALL WRITE_0D_L( RNF_useSalt, INDEX_NONE,
     &  'RNF_useSalt =',
     &  ' /* apply runoff_salinity if present */')
      CALL WRITE_0D_L( RNF_usePtracers, INDEX_NONE,
     &  'RNF_usePtracers =',
     &  ' /* apply runoff_ptracer_* if present */')
      CALL WRITE_0D_RL( RNF_monFreq, INDEX_NONE,
     &  'RNF_monFreq =',
     &  ' /* RNF monitor frequency [s] */')
      CALL WRITE_0D_I( RNF_debugLev, INDEX_NONE,
     &  'RNF_debugLev =',
     &  ' /* select RNF-debug printing level */')

C--   The time handling RNF_TIME_SETUP resolved from the file's time
C     axis and these overrides. RNF_recPeriod, RNF_recStart and
C     RNF_recCycle are what the pkg/exf record selection is called
C     with, so a log says which records a run reads and from when.
      CALL WRITE_0D_RL( RNF_recPeriod, INDEX_NONE,
     &  'RNF_recPeriod =',
     &  ' /* period in use: 0, s, -12 or -1 */')
      CALL WRITE_0D_RL( RNF_recStart, INDEX_NONE,
     &  'RNF_recStart =',
     &  ' /* time of record 1 in use [s] */')
      CALL WRITE_0D_RL( RNF_recCycle, INDEX_NONE,
     &  'RNF_recCycle =',
     &  ' /* repeat cycle in use [s], 0: none */')
      CALL WRITE_0D_I( RNF_recDate1, INDEX_NONE,
     &  'RNF_recDate1 =',
     &  ' /* date of record 1 YYYYMMDD (0: no date) */')
      CALL WRITE_0D_I( RNF_recDate2, INDEX_NONE,
     &  'RNF_recDate2 =',
     &  ' /* date of record 1 HHMMSS */')

C--   Constants fixed by the model contract
      CALL WRITE_0D_RL( RNF_fracTol, INDEX_NONE,
     &  'RNF_fracTol =',
     &  ' /* tolerance on the fraction sum of a source */')
      CALL WRITE_0D_RL( RNF_areaTol, INDEX_NONE,
     &  'RNF_areaTol =',
     &  ' /* relative tolerance on the target cell area */')
      CALL WRITE_0D_RL( RNF_lonLatTol, INDEX_NONE,
     &  'RNF_lonLatTol =',
     &  ' /* cell centre tolerance, in cell spacings */')
      CALL WRITE_0D_RL( RNF_srcFluxMax, INDEX_NONE,
     &  'RNF_srcFluxMax =',
     &  ' /* largest flux one source may carry [m^3/s] */')
      CALL WRITE_0D_RL( RNF_cellVolMax, INDEX_NONE,
     &  'RNF_cellVolMax =',
     &  ' /* largest share of a cell top-layer volume one step */')

C--   Which runoff bounds this run was actually held to. EXF_CHECK_RANGE
C     conditions TWO of its tests on useRNF - it skips the runoff upper
C     bound of 1.E-6 m/s and tests sflux+runoff in place of sflux -
C     because one sparse point source exceeds both by construction. A
C     reader of this log would otherwise have no way to tell which
C     bounds applied, so all five statements are printed: the two
C     relaxations, the two package bounds that replace them - per
C     source in RNF_NC_READ_ONE and per cell in RNF_EXF_RUNOFF - and
C     the negative-runoff test that is NOT skipped. Review A of
C     RUNOFF-030 correction round 2 found the freshwater one missing
C     here; the per-cell line is RUNOFF-040.
      WRITE(msgBuf,'(2A)') 'RNF_SUMMARY: pkg/exf skips its runoff',
     &  ' upper bound (1.E-6 m/s) with useRNF=.TRUE.;'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(2A,1PE16.8,A)') 'RNF_SUMMARY: RNF_NC_READ_ONE',
     &  ' refuses a source flux above RNF_srcFluxMax =',
     &  RNF_srcFluxMax, ' m^3/s instead'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(2A)') 'RNF_SUMMARY: its sflux bound tests',
     &  ' sflux+runoff, i.e. evap-precip, with useRNF=.TRUE.'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(2A)') 'RNF_SUMMARY: the negative-runoff test of',
     &  ' EXF_CHECK_RANGE still applies per cell (at nIter0)'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(2A,1PE16.8,A)') 'RNF_SUMMARY: RNF_EXF_RUNOFF',
     &  ' refuses a cell above RNF_cellVolMax =', RNF_cellVolMax,
     &  ' of its top-layer volume per step'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )

C--   RUNOFF-041: RNF_cellVolMax's "fractional dilution" reading (as
C     opposed to its primary "fractional volume change" reading, RNF.h)
C     holds only where dTtracerLev(1) = deltaTFreeSurf. Report rather
C     than bound (both reviewers recommended against bounding with
C     MAX(deltaTFreeSurf,dTtracerLev(ks)), and Arch accepted): print
C     the mismatch and its ratio when it occurs, so a run in that
C     configuration says so in its own STDOUT instead of leaving the
C     looser dilution reading unannounced. Silent (no WRITE at all)
C     when the two are equal, which is every committed experiment
C     today (lab_sea by default, cs32 by setting deltaTFreeSurf
C     explicitly).
      IF ( dTtracerLev(1) .NE. deltaTFreeSurf ) THEN
        dtRatio = dTtracerLev(1)/deltaTFreeSurf
        WRITE(msgBuf,'(A,1PE16.8,A,E16.8,A,E16.8,A)')
     &    'RNF_SUMMARY: dTtracerLev(1) =', dTtracerLev(1),
     &    ' differs from deltaTFreeSurf =', deltaTFreeSurf,
     &    ': RNF_cellVolMax dilution reading is looser by their ratio',
     &    dtRatio, '; the volume reading is unaffected'
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                      SQUEEZE_RIGHT, myThid )
      ENDIF

C--   Array bounds of RNF_SIZE.h
      CALL WRITE_0D_I( RNF_nSrcTile, INDEX_NONE,
     &  'RNF_nSrcTile =',
     &  ' /* max number of sources on one tile */')
      CALL WRITE_0D_I( RNF_nTgtTile, INDEX_NONE,
     &  'RNF_nTgtTile =',
     &  ' /* max number of target entries on one tile */')
      CALL WRITE_0D_I( RNF_nBuf, INDEX_NONE,
     &  'RNF_nBuf =',
     &  ' /* chunk length for NetCDF reads */')
      CALL WRITE_0D_I( RNF_nTr, INDEX_NONE,
     &  'RNF_nTr =',
     &  ' /* max number of runoff tracers */')
      CALL WRITE_0D_I( RNF_nFile, INDEX_NONE,
     &  'RNF_nFile =',
     &  ' /* number of runoff files */')

C--   What the static read of RNF_file found
      CALL WRITE_0D_I( RNF_nSrcFile, INDEX_NONE,
     &  'RNF_nSrcFile =',
     &  ' /* sources in the file */')
      CALL WRITE_0D_I( RNF_nTgtFile, INDEX_NONE,
     &  'RNF_nTgtFile =',
     &  ' /* target entries in the file */')
      CALL WRITE_0D_I( RNF_nRecFile, INDEX_NONE,
     &  'RNF_nRecFile =',
     &  ' /* time records in the file */')
      CALL WRITE_0D_I( RNF_nSrcProc, INDEX_NONE,
     &  'RNF_nSrcProc =',
     &  ' /* tile-source entries of this process */')
      CALL WRITE_0D_I( RNF_nTgtProc, INDEX_NONE,
     &  'RNF_nTgtProc =',
     &  ' /* target entries owned by this process */')
      CALL WRITE_0D_I( RNF_nTgtOwned, INDEX_NONE,
     &  'RNF_nTgtOwned =',
     &  ' /* target entries owned over all processes */')
C--   Whether the optional cell-centre check ran, so that a run is
C     never taken to be guarded by it when the file did not carry
C     target_lon/target_lat or the grid's XC,YC are not degrees.
      CALL WRITE_0D_L( RNF_lonLatChk, INDEX_NONE,
     &  'RNF_lonLatChk =',
     &  ' /* target_lon,lat were checked against XC,YC */')

C--   Which time series of the file are applied, and what each term
C     will do (package design, decisions 3 and 4). RNF_applyT says
C     the package sets the exf runoftemp field, whose heat
C     EXF_MAPFIELDS adds to Qnet; RNF_applyS is false when the salt
C     term is identically zero, in which case surfaceForcingS is not
C     touched. RNF_lagFlds says which time level the salt and tracer
C     terms use: the current step, or the step before it in the one
C     freshwater formulation where the model's own PmEpR lags.
      CALL WRITE_0D_L( RNF_hasTemp, INDEX_NONE,
     &  'RNF_hasTemp =',
     &  ' /* runoff_temperature is read from the file */')
      CALL WRITE_0D_L( RNF_hasSalt, INDEX_NONE,
     &  'RNF_hasSalt =',
     &  ' /* runoff_salinity is read from the file */')
      CALL WRITE_0D_I( RNF_nTrUse, INDEX_NONE,
     &  'RNF_nTrUse =',
     &  ' /* runoff tracers read from the file */')
      CALL WRITE_0D_L( RNF_applyT, INDEX_NONE,
     &  'RNF_applyT =',
     &  ' /* pkg/rnf sets the exf runoftemp field */')
      CALL WRITE_0D_L( RNF_applyS, INDEX_NONE,
     &  'RNF_applyS =',
     &  ' /* the salinity term can be non-zero */')
      CALL WRITE_0D_L( RNF_lagFlds, INDEX_NONE,
     &  'RNF_lagFlds =',
     &  ' /* S, tracer terms use the previous step */')
      DO k = 1, RNF_nTrUse
        trNam = RNF_trNam(k)
        iL = MAX( 1, ILNBLNK(trNam) )
        WRITE(msgBuf,'(A,I3,3A,I3)') 'RNF_SUMMARY: runoff tracer',
     &    k, ' is runoff_ptracer_', trNam(1:iL),
     &    ', applied to ptracer', RNF_trPtr(k)
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                      SQUEEZE_RIGHT, myThid )
      ENDDO

      CALL PRINT_MESSAGE( blkLin, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)')
     &'// ======================================================='
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)')
     &'// Sparse runoff (RNF) configuration  >>> END <<<'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)')
     &'// ======================================================='
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )
      CALL PRINT_MESSAGE( blkLin, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )

      _END_MASTER(myThid)

#endif /* ALLOW_RNF */

      RETURN
      END
