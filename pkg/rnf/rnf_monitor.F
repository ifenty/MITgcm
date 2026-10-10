#include "RNF_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: RNF_MONITOR

C !INTERFACE: ==========================================================
      SUBROUTINE RNF_MONITOR( myTime, myIter, myThid )

C !DESCRIPTION:
C     Write the runoff monitor block at RNF_monFreq, in the style of
C     pkg/exf/exf_monitor.F:189-198 (package design, decision 8). Called
C     from RNF_FIELDS_LOAD at every step, as EXF_GETFORCING calls
C     EXF_MONITOR, so the two monitors are written at the steps the
C     model's own forcing monitor is.
C
C     Three groups:
C
C       o MON_WRITESTATS_RL statistics of the volume-flux field
C         RNF_vflx, as exf's own "_runoff" entry (exf_monitor.F:189-193)
C         reports when a dense runofffile is in use; the package
C         monitor is the only runoff entry of the monitor block when
C         pkg/rnf is, because exf's own entry is gated on runofffile
C         being set, which it is not under pkg/rnf;
C       o the global sums of decision 8: Sum m*rA/rhoConstFresh (m^3/s,
C         the applied volume, which RNF_INIT_VARIA's start-time pair
C         also measures), Sum Cp*(mT)*rA (W) and Sum (mS)*rA (g/s),
C         plus one per runoff tracer, Sum (mC_n)*rA (kg/s times tracer
C         units). Each is cheap -- one more multiply-accumulate over
C         the tile interior, already visited by RNF_DIAGNOSTICS_FILL --
C         and consistent with the volume and heat sums, both of which
C         are the same per-cell quantities (RNF_vflx*rA,
C         RNF_mXT*rA) weighted the same way, so the tracer sums are
C         included rather than left out as an extra cost;
C       o the number of sources and target entries in use,
C         RNF_nSrcFile and RNF_nTgtOwned. Both are set once, by the
C         static read of the file (RNF_INIT_FIXED) and are therefore
C         the same number on every process and at every monitor time:
C         printing them needs no reduction and is tiling-independent
C         by construction, unlike RNF_nSrcProc/RNF_nTgtProc, which
C         are per-process and are not what is printed here.
C
C     **The global sums are tiling- and process-independent.** Each is
C     a per-tile sum over the tile interior (RNF_vflx, RNF_mXT, RNF_mXS
C     and RNF_mXTr are built identically whichever process owns a tile:
C     RNF_FIELDS_BUILD depends only on the input and the tile's own
C     target list), reduced with GLOBAL_SUM_TILE_RL, which sums over
C     every tile of every process in the fixed global tile order
C     (eesupp/src/global_sum_tile.F), independent of how many
C     processes own them. So a 1-process and a 2-process run reduce
C     the same per-tile values in the same order and print the same
C     line, as pkg/monitor's own statistics already do for dynstat_*.

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "RNF_SIZE.h"
#include "RNF.h"
#ifdef ALLOW_MONITOR
# include "MONITOR.h"
#endif /* ALLOW_MONITOR */

C !INPUT PARAMETERS: ===================================================
C  myTime               :: current time in simulation
C  myIter               :: current iteration number
C  myThid               :: my thread Id number
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
CEOP

#ifdef ALLOW_RNF
#ifdef ALLOW_MONITOR

C     !FUNCTIONS:
      LOGICAL  DIFFERENT_MULTIPLE
      EXTERNAL DIFFERENT_MULTIPLE
      LOGICAL  MASTER_CPU_IO
      EXTERNAL MASTER_CPU_IO

C !LOCAL VARIABLES: ====================================================
C  bi, bj, i, j, iTr     :: tile and loop indices
C  dummyRL               :: unused MON_WRITESTATS_RL output
C  tileVol, tileHeat,
C  tileSalt, tileTrc     :: per-tile sums, before GLOBAL_SUM_TILE_RL
C  rnfVol, rnfHeat,
C  rnfSalt, rnfTrc       :: the reduced global sums
C  msgBuf                :: message buffer
      INTEGER bi, bj, i, j, iTr
      _RL     dummyRL(6)
      _RL     tileVol (nSx,nSy)
      _RL     tileHeat(nSx,nSy)
      _RL     tileSalt(nSx,nSy)
      _RL     tileTrc (nSx,nSy,RNF_nTr)
      _RL     rnfVol, rnfHeat, rnfSalt
      _RL     rnfTrc(RNF_nTr)
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      CHARACTER*9 footTrc

      IF ( DIFFERENT_MULTIPLE(RNF_monFreq,myTime,deltaTClock) ) THEN

        IF ( MASTER_CPU_IO(myThid) ) THEN
          IF ( monitor_stdio ) THEN
            mon_write_stdout = .TRUE.
          ELSE
            mon_write_stdout = .FALSE.
          ENDIF
C--       No MNC output of this monitor today: mon_write_mnc stays
C         false, and nothing below is written through MNC_CW_*.
          mon_write_mnc = .FALSE.
          IF ( mon_write_stdout ) THEN
            WRITE(msgBuf,'(2A)') '// ===========================',
     &           '============================'
            CALL PRINT_MESSAGE(msgBuf, mon_ioUnit, SQUEEZE_RIGHT, 1)
            WRITE(msgBuf,'(A)') '// Begin MONITOR RNF statistics'
            CALL PRINT_MESSAGE(msgBuf, mon_ioUnit, SQUEEZE_RIGHT, 1)
            WRITE(msgBuf,'(2A)') '// ===========================',
     &           '============================'
            CALL PRINT_MESSAGE(msgBuf, mon_ioUnit, SQUEEZE_RIGHT, 1)
          ENDIF
        ENDIF

        CALL MON_SET_PREF('rnf',myThid)
        CALL MON_OUT_I ('_tsnumber', myIter,mon_string_none,myThid)
        CALL MON_OUT_RL('_time_sec', myTime,mon_string_none,myThid)

C--     Statistics of the volume-flux field, as exf's own "_runoff"
C       entry does for the dense path.
        CALL MON_WRITESTATS_RL( 1, RNF_vflx, '_vflx',
     &           maskInC, maskInC, rA, drF, dummyRL, myThid )

C--     The global sums: one multiply-accumulate per cell over what
C       RNF_DIAGNOSTICS_FILL has already visited this step.
        DO bj = 1, nSy
         DO bi = 1, nSx
          tileVol (bi,bj) = 0. _d 0
          tileHeat(bi,bj) = 0. _d 0
          tileSalt(bi,bj) = 0. _d 0
          DO iTr = 1, RNF_nTr
            tileTrc(bi,bj,iTr) = 0. _d 0
          ENDDO
         ENDDO
        ENDDO
        DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
          DO j = 1, sNy
           DO i = 1, sNx
            tileVol(bi,bj)  = tileVol(bi,bj)
     &        + RNF_vflx(i,j,bi,bj)*rA(i,j,bi,bj)
            tileHeat(bi,bj) = tileHeat(bi,bj)
     &        + RNF_mXT(i,j,bi,bj)*rA(i,j,bi,bj)
            tileSalt(bi,bj) = tileSalt(bi,bj)
     &        + RNF_mXS(i,j,bi,bj)*rA(i,j,bi,bj)
           ENDDO
          ENDDO
          DO iTr = 1, RNF_nTrUse
           DO j = 1, sNy
            DO i = 1, sNx
             tileTrc(bi,bj,iTr) = tileTrc(bi,bj,iTr)
     &         + RNF_mXTr(i,j,bi,bj,iTr)*rA(i,j,bi,bj)
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO
        CALL GLOBAL_SUM_TILE_RL( tileVol,  rnfVol,  myThid )
        CALL GLOBAL_SUM_TILE_RL( tileHeat, rnfHeat, myThid )
        CALL GLOBAL_SUM_TILE_RL( tileSalt, rnfSalt, myThid )
        DO iTr = 1, RNF_nTrUse
          CALL GLOBAL_SUM_TILE_RL( tileTrc(1,1,iTr), rnfTrc(iTr),
     &                             myThid )
        ENDDO

C--     rnfVol is already Sum RNF_vflx*rA, which is Sum m*rA/rhoConstFresh
C       with m = rhoConstFresh*RNF_vflx: the rhoConstFresh cancels, so no
C       further division is made here (an earlier draft of this routine
C       divided again, which tests/rnf/diagnostics_check.py's monitor
C       comparison against the file's own source sum caught at a factor
C       of 1/rhoConstFresh).
        rnfHeat = HeatCapacity_Cp*rnfHeat

        CALL MON_OUT_RL('_volsum ', rnfVol,  mon_string_none, myThid)
        CALL MON_OUT_RL('_heatsum', rnfHeat, mon_string_none, myThid)
        CALL MON_OUT_RL('_saltsum', rnfSalt, mon_string_none, myThid)
        DO iTr = 1, RNF_nTrUse
          WRITE(footTrc,'(A,I2.2,A)') '_trc', iTr, 'sum'
          CALL MON_OUT_RL(footTrc, rnfTrc(iTr), mon_string_none,
     &                    myThid)
        ENDDO

C--     The number of sources and target entries in use: both are set
C       once by the static read of the file and are the same on every
C       process, so no reduction is needed here.
        CALL MON_OUT_I('_nsrc', RNF_nSrcFile, mon_string_none, myThid)
        CALL MON_OUT_I('_ntgt', RNF_nTgtOwned, mon_string_none,myThid)

        IF ( MASTER_CPU_IO(myThid) ) THEN
          IF ( mon_write_stdout ) THEN
            WRITE(msgBuf,'(2A)') '// ===========================',
     &           '============================'
            CALL PRINT_MESSAGE(msgBuf, mon_ioUnit, SQUEEZE_RIGHT, 1)
            WRITE(msgBuf,'(A)') '// End MONITOR RNF statistics'
            CALL PRINT_MESSAGE(msgBuf, mon_ioUnit, SQUEEZE_RIGHT, 1)
            WRITE(msgBuf,'(2A)') '// ===========================',
     &           '============================'
            CALL PRINT_MESSAGE(msgBuf, mon_ioUnit, SQUEEZE_RIGHT, 1)
          ENDIF
          mon_write_stdout = .FALSE.
          mon_write_mnc    = .FALSE.
        ENDIF

      ENDIF

#endif /* ALLOW_MONITOR */
#endif /* ALLOW_RNF */

      RETURN
      END
