#include "RNF_OPTIONS.h"
#ifdef ALLOW_EXF
# include "EXF_OPTIONS.h"
#endif /* ALLOW_EXF */

CBOP
C !ROUTINE: RNF_EXF_RUNOFF

C !INTERFACE: ==========================================================
      SUBROUTINE RNF_EXF_RUNOFF( myTime, myIter, myThid )

C !DESCRIPTION:
C     Set the exf field "runoff" (m/s) from the sparse runoff volume
C     flux of the current step (package design, decision 2), and refuse
C     a cell at which one step of that runoff would add more than
C     RNF_cellVolMax of the cell's own top-layer volume.
C     Called from EXF_GETFFIELDS after the reads of the dense runoff
C     and runoff temperature files and before the control variables
C     are added, so that exf and everything downstream (EmPmR, the
C     free surface, the sea-ice packages) use the field as they use a
C     dense one.
C
C     The field is assigned, not added to: EXF_SET_FLD does nothing
C     when runofffile is blank, so nothing else resets the array, and
C     the control block adds to it at every step. The assignment
C     covers the whole tile interior, the cells without runoff
C     included. The halo is left alone, as the dense path leaves it:
C     exf does not exchange runoff, and RNF_CHECK refuses a non-zero
C     runoffconst, so the halo holds the zero of EXF_INIT_FLD.
C
C     **The per-cell magnitude bound.** The bound is enforced here
C     because this routine is where the three things it needs are all
C     available together -- rA (through RNF_vflx, which is already a
C     flux per unit area), the top-layer thickness and the time step
C     -- and where the field that is actually applied exists. That is
C     a statement about this routine and not an exclusivity claim
C     about the package: RNF_LOAD_AT also has rA, and the tendency
C     routines also have drF and hFacC. What they do not have is the
C     assembled per-cell field together with the step it is applied
C     over. The bound is on the dimensionless ratio
C
C       f(c) = |RNF_vflx(c)|*deltaTFreeSurf / ( drF(ks)*hFacC(c,ks) )
C
C     so one constant serves every grid; RNF.h derives the value and
C     states what the bound does and does not cover. The level ks is
C     the surface level of the column, by the rule RNF_INIT_FIXED uses
C     for RNF_tgtK (rnf_init_fixed.F:217-218 and :709-712), evaluated
C     here so that it stays right if pkg/shelfice reassigns kSurfC
C     during the run.
C
C     The test is written as a negated .LE. rather than as a .GT. so
C     that a value that is not a number is refused instead of
C     compared: every comparison with a NaN is false, so the .GT. form
C     would accept it. That is the same shape, and for the same
C     reason, as the area and cell-centre tests of RNF_INIT_FIXED
C     (rnf_init_fixed.F:730). It is multiplied out rather than
C     divided, so a column with no thickness at all needs no special
C     case: the limit is then zero and any non-zero runoff on it is
C     refused.
C
C     **Every step, not only nIter0.** The target table is static, so
C     a collapse of targets onto one cell is already visible at
C     nIter0; the flux time series is not, so a file whose record 1 is
C     innocent and whose record 500 is not would pass a check made
C     once. RNF_srcFluxMax is checked on every record for exactly that
C     reason (RNF_NC_READ_ONE), and a per-cell bound that ran only at
C     nIter0 would leave the aggregate of a later record unguarded
C     while the per-source value of that same record was checked.
C     The cost is a second pass over the same tile interior: one
C     compare per cell, and the thickness lookup and the two
C     multiplies only on the cells that carry runoff, which are the
C     target cells. The assignment loop above is deliberately left as
C     a plain copy rather than having the test folded into it, so that
C     it stays a straight vectorisable copy. On top of that there is
C     one GLOBAL_SUM_INT per step, which is what lets a refusal seen
C     on one tile stop every process (the same reason RNF_INIT_FIXED
C     reduces its counters at rnf_init_fixed.F:933-935), and it is the
C     only part of the cost that does not scale down with the tile.
C     Measured on lab_sea/input.rnof_sp_const, 48 steps, serial: the
C     EXF_GETFORCING timer section that contains this routine is 1.1%
C     of MAIN LOOP (4.52E-2 s of 3.999 s before, 4.56-4.80E-2 s of
C     3.962-4.047 s after, four samples), so the added cost is inside
C     the run-to-run scatter of the host, below about 60 microseconds
C     per step on a 20x16 grid of four tiles.

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "RNF_SIZE.h"
#include "RNF.h"
#ifdef ALLOW_EXF
# include "EXF_PARAM.h"
# include "EXF_FIELDS.h"
#endif /* ALLOW_EXF */

C !INPUT PARAMETERS: ===================================================
C  myTime               :: current time in simulation
C  myIter               :: current iteration number
C  myThid               :: my thread Id number
      _RL     myTime
      INTEGER myIter
      INTEGER myThid

C !OUTPUT PARAMETERS: ==================================================
C  none

#ifdef ALLOW_RNF
#ifdef ALLOW_EXF
#ifdef ALLOW_RUNOFF

C !LOCAL VARIABLES: ====================================================
C  bi, bj               :: tile indices
C  i, j                 :: loop indices
C  kSurf                :: surface level of a column where it is the
C                          same for every column of the grid
C  kLoc                 :: surface level of this column
C  hTop                 :: its top-layer thickness [m]
C  vLim                 :: largest runoff this cell may take [m/s]
C  nErr                 :: cells of this thread above the bound
C  nMsg                 :: messages this thread has printed
      INTEGER bi, bj, i, j
      INTEGER kSurf, kLoc
      INTEGER nErr, nMsg
      _RL     hTop, vLim
      CHARACTER*(MAX_LEN_MBUF) msgBuf
CEOP

      kSurf = 1
      IF ( usingPCoords ) kSurf = Nr
      nErr = 0
      nMsg = 0

      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)
        DO j = 1, sNy
         DO i = 1, sNx
          runoff(i,j,bi,bj) = RNF_vflx(i,j,bi,bj)
         ENDDO
        ENDDO

C--     The per-cell magnitude of what was just assigned. A cell with
C       no runoff passes the bound whatever its thickness, so it is
C       skipped before the thickness is looked up.
        DO j = 1, sNy
         DO i = 1, sNx
          IF ( RNF_vflx(i,j,bi,bj) .NE. 0. _d 0 ) THEN
           kLoc = kSurf
#ifdef ALLOW_SHELFICE
           IF ( usingZCoords .AND. useShelfIce )
     &          kLoc = kSurfC(i,j,bi,bj)
#endif /* ALLOW_SHELFICE */
           hTop = 0. _d 0
           IF ( kLoc.GE.1 .AND. kLoc.LE.Nr )
     &          hTop = drF(kLoc)*_hFacC(i,j,kLoc,bi,bj)
           IF ( .NOT.( ABS(RNF_vflx(i,j,bi,bj))*deltaTFreeSurf
     &                 .LE. RNF_cellVolMax*hTop ) ) THEN
            nErr = nErr + 1
            IF ( nMsg .LT. RNF_maxErrMsg ) THEN
             nMsg = nMsg + 1
             vLim = RNF_cellVolMax*hTop/deltaTFreeSurf
             WRITE(msgBuf,'(2A,4I5,A,1PE16.8,A,E16.8)')
     &         'RNF_EXF_RUNOFF: RNF: runoff out of range:',
     &         ' cell (i,j,bi,bj) =', i, j, bi, bj,
     &         ', value', RNF_vflx(i,j,bi,bj), ', limit', vLim
             CALL PRINT_ERROR( msgBuf, myThid )
             WRITE(msgBuf,'(A,1P2E16.8,A,E16.8,A,E16.8)')
     &         'RNF_EXF_RUNOFF: its XC,YC =',
     &         XC(i,j,bi,bj), YC(i,j,bi,bj),
     &         ', top-layer thickness', hTop,
     &         ' m, deltaTFreeSurf', deltaTFreeSurf
             CALL PRINT_ERROR( msgBuf, myThid )
            ENDIF
           ENDIF
          ENDIF
         ENDDO
        ENDDO

       ENDDO
      ENDDO

C--   A cell above the bound is on one tile of one process, and every
C     process has to stop: ALL_PROC_DIE hangs unless all of them call
C     it. So the count is reduced over the threads and the processes,
C     as RNF_INIT_FIXED reduces the counters of its own one-tile
C     refusals (rnf_init_fixed.F:933-935).
      CALL GLOBAL_SUM_INT( nErr, myThid )

      IF ( nErr .GE. 1 ) THEN
       _BEGIN_MASTER(myThid)
       WRITE(msgBuf,'(A,I8,A,I10)') 'RNF_EXF_RUNOFF: RNF:', nErr,
     &   ' cell(s) with too much runoff at iteration', myIter
       CALL PRINT_ERROR( msgBuf, myThid )
       WRITE(msgBuf,'(2A,1PE16.8,A)') 'RNF_EXF_RUNOFF: one step adds',
     &   ' more than RNF_cellVolMax =', RNF_cellVolMax,
     &   ' of the cell top-layer volume'
       CALL PRINT_ERROR( msgBuf, myThid )
       WRITE(msgBuf,'(2A)') 'RNF_EXF_RUNOFF: several sources on one',
     &   ' cell add up here: check target_cell of the file'
       CALL PRINT_ERROR( msgBuf, myThid )
       _END_MASTER(myThid)
       CALL ALL_PROC_DIE( myThid )
       STOP 'ABNORMAL END: S/R RNF_EXF_RUNOFF'
      ENDIF

#endif /* ALLOW_RUNOFF */
#endif /* ALLOW_EXF */
#endif /* ALLOW_RNF */

      RETURN
      END
