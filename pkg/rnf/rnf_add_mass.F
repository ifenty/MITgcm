#include "RNF_OPTIONS.h"

CBOP
C !ROUTINE: RNF_ADD_MASS

C !INTERFACE: ==========================================================
      SUBROUTINE RNF_ADD_MASS( myTime, myIter, myThid )

C !DESCRIPTION:
C     The volume of the generalized runoff route (RNF_applyAsRain = F,
C     the default; package design, decision 2): add the
C     runoff of the current step to the model's interior mass source
C     addMass [kg/s], at the level of each target entry,
C
C       addMass(i,j,k) += rhoConstFresh * flux_s * frac_{s,c}
C
C     for every entry (s,c) of the tile with RNF_tgtK = k. The model
C     then carries the volume: CALC_DIV_GHAT puts the depth sum of
C     addMass into the free-surface right-hand side and
C     INTEGR_CONTINUITY into dEtaHdt and the vertical velocity
C     (model/src/calc_div_ghat.F:127-135,
C     model/src/integr_continuity.F:108-133), so the free surface
C     rises by the same amount whatever the level, and only the
C     vertical velocity under the free surface depends on it. The
C     factor rhoConstFresh, with the model's own mass2rUnit =
C     1/rhoConst applied to addMass, puts in the volume the legacy
C     route puts in through EmPmR = -rhoConstFresh*runoff, so the two
C     routes carry the same volume.
C
C     Called from LOAD_FIELDS_DRIVER right after RNF_FIELDS_LOAD,
C     whose RNF_srcFlux this routine reads, i.e. after the reset of
C     addMass at the start of that routine
C     (model/src/load_fields_driver.F:116-128) and before
C     CALC_DIV_GHAT and INTEGR_CONTINUITY read it. The reset happens
C     only with selectAddFluid /= 0 and a blank addMassFile, which
C     RNF_CHECK requires of this route, so the runoff is not added to
C     the array again at every step. Other packages may add to addMass
C     later in the step (pkg/steep_icecavity, from DO_OCEANIC_PHYS);
C     nothing else writes it between the reset and this routine, which
C     is what makes the RNFaddM diagnostic below the runoff share
C     alone.
C
C     A no-op on the legacy route (RNF_applyAsRain = T), whose volume
C     is the exf runoff field (RNF_EXF_RUNOFF).
C
C     **The per-cell magnitude bound.** The legacy route refuses a
C     cell where one step of runoff exceeds RNF_cellVolMax of the
C     surface cell (RNF_EXF_RUNOFF, RNF.h). The same bound is applied
C     here to the cell each entry is actually in, at its own level:
C
C       f(i,j,k) = |sum_e flux*frac/rA|*deltaTFreeSurf
C                  / ( drF(k)*hFacC(i,j,k) )
C
C     summed over the entries of that cell and level, so that several
C     sources on one cell add up as they do on the legacy route. Same
C     constant, same step, same negated .LE. so that a NaN is refused.
C     What changes is only the thickness, which is the target level's
C     and not the surface level's. Every step, as on the legacy route.

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "FFIELDS.h"
#include "RNF_SIZE.h"
#include "RNF.h"

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
#ifdef ALLOW_ADDFLUID

C !LOCAL VARIABLES: ====================================================
C  bi, bj               :: tile indices
C  i, j, k, kLev        :: cell, level and entry indices
C  kSrc                 :: local source of an entry
C  kMin, kMax           :: range of the levels of the entries of a tile
C  vSum                 :: volume flux per unit area of one level [m/s]
C  hLev                 :: thickness of a target cell [m]
C  vLim                 :: largest runoff the cell may take [m/s]
C  nErr                 :: cells of this thread above the bound
C  nMsg                 :: messages this thread has printed
      INTEGER bi, bj, i, j, k, kLev, kSrc
      INTEGER kMin, kMax
      INTEGER nErr, nMsg
      _RL     vSum(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     hLev, vLim
      CHARACTER*(MAX_LEN_MBUF) msgBuf
CEOP

      IF ( RNF_applyAsRain ) RETURN

      nErr = 0
      nMsg = 0

      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)

C--     The mass source, one entry at a time, in the order of the
C       target table
        kMin = Nr + 1
        kMax = 0
        DO k = 1, RNF_nTgt(bi,bj)
          i = RNF_tgtI(k,bi,bj)
          j = RNF_tgtJ(k,bi,bj)
          kLev = RNF_tgtK(k,bi,bj)
          kSrc = RNF_tgtSrc(k,bi,bj)
          addMass(i,j,kLev,bi,bj) = addMass(i,j,kLev,bi,bj)
     &      + rhoConstFresh*RNF_srcFlux(kSrc,bi,bj)*RNF_tgtFrac(k,bi,bj)
          kMin = MIN( kMin, kLev )
          kMax = MAX( kMax, kLev )
        ENDDO

C--     The per-cell bound, one level at a time, so that the entries
C       of one cell and level are summed before the test. vSum is
C       zero everywhere between levels: each cell is set back to zero
C       once it has been tested.
        DO j = 1-OLy, sNy+OLy
         DO i = 1-OLx, sNx+OLx
          vSum(i,j) = 0. _d 0
         ENDDO
        ENDDO
        DO kLev = kMin, kMax
         DO k = 1, RNF_nTgt(bi,bj)
          IF ( RNF_tgtK(k,bi,bj) .EQ. kLev ) THEN
           i = RNF_tgtI(k,bi,bj)
           j = RNF_tgtJ(k,bi,bj)
           kSrc = RNF_tgtSrc(k,bi,bj)
           vSum(i,j) = vSum(i,j)
     &       + RNF_srcFlux(kSrc,bi,bj)*RNF_tgtFrac(k,bi,bj)
     &         /rA(i,j,bi,bj)
          ENDIF
         ENDDO
         DO k = 1, RNF_nTgt(bi,bj)
          IF ( RNF_tgtK(k,bi,bj) .EQ. kLev ) THEN
           i = RNF_tgtI(k,bi,bj)
           j = RNF_tgtJ(k,bi,bj)
           IF ( vSum(i,j) .NE. 0. _d 0 ) THEN
            hLev = drF(kLev)*_hFacC(i,j,kLev,bi,bj)
            IF ( .NOT.( ABS(vSum(i,j))*deltaTFreeSurf
     &                  .LE. RNF_cellVolMax*hLev ) ) THEN
             nErr = nErr + 1
             IF ( nMsg .LT. RNF_maxErrMsg ) THEN
              nMsg = nMsg + 1
              vLim = RNF_cellVolMax*hLev/deltaTFreeSurf
              WRITE(msgBuf,'(2A,5I5,A,1PE16.8,A,E16.8)')
     &          'RNF_ADD_MASS: RNF: runoff out of range:',
     &          ' cell (i,j,k,bi,bj) =', i, j, kLev, bi, bj,
     &          ', value', vSum(i,j), ', limit', vLim
              CALL PRINT_ERROR( msgBuf, myThid )
              WRITE(msgBuf,'(A,1P2E16.8,A,E16.8,A,E16.8)')
     &          'RNF_ADD_MASS: its XC,YC =',
     &          XC(i,j,bi,bj), YC(i,j,bi,bj),
     &          ', cell thickness', hLev,
     &          ' m, deltaTFreeSurf', deltaTFreeSurf
              CALL PRINT_ERROR( msgBuf, myThid )
             ENDIF
            ENDIF
            vSum(i,j) = 0. _d 0
           ENDIF
          ENDIF
         ENDDO
        ENDDO

#ifdef ALLOW_DIAGNOSTICS
C--     The runoff share of addMass: nothing else has written the
C       array since its reset (see above)
        IF ( useDiagnostics ) THEN
         CALL DIAGNOSTICS_FILL( addMass(1-OLx,1-OLy,1,bi,bj),
     &        'RNFaddM ', 0, Nr, 2, bi, bj, myThid )
        ENDIF
#endif /* ALLOW_DIAGNOSTICS */

       ENDDO
      ENDDO

C--   A cell above the bound is on one tile of one process, and every
C     process has to stop, as in RNF_EXF_RUNOFF
      CALL GLOBAL_SUM_INT( nErr, myThid )

      IF ( nErr .GE. 1 ) THEN
       _BEGIN_MASTER(myThid)
       WRITE(msgBuf,'(A,I8,A,I10)') 'RNF_ADD_MASS: RNF:', nErr,
     &   ' cell(s) with too much runoff at iteration', myIter
       CALL PRINT_ERROR( msgBuf, myThid )
       WRITE(msgBuf,'(2A,1PE16.8,A)') 'RNF_ADD_MASS: one step adds',
     &   ' more than RNF_cellVolMax =', RNF_cellVolMax,
     &   ' of the target cell volume'
       CALL PRINT_ERROR( msgBuf, myThid )
       WRITE(msgBuf,'(2A)') 'RNF_ADD_MASS: several sources on one',
     &   ' cell add up here: check target_cell and target_level'
       CALL PRINT_ERROR( msgBuf, myThid )
       _END_MASTER(myThid)
       CALL ALL_PROC_DIE( myThid )
       STOP 'ABNORMAL END: S/R RNF_ADD_MASS'
      ENDIF

#endif /* ALLOW_ADDFLUID */
#endif /* ALLOW_RNF */

      RETURN
      END
