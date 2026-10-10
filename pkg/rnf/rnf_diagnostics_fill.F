#include "RNF_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: RNF_DIAGNOSTICS_FILL

C !INTERFACE: ==========================================================
      SUBROUTINE RNF_DIAGNOSTICS_FILL( myTime, myIter, myThid )

C !DESCRIPTION:
C     Fill the five input-field diagnostics of decision 8
C     (RNFvflx, RNFmflx, RNFtemp, RNFsaln, RNFnsrc), registered by
C     RNF_DIAGNOSTICS_INIT. They depend only on the input, i.e. on the
C     dense fields RNF_FIELDS_BUILD has just built, so this is called
C     once per step from RNF_FIELDS_LOAD, at the end of
C     RNF_FIELDS_BUILD -- unlike RNFqnet, RNFsflx and RNFtfNN
C     (RUNOFF-031), which are filled where their own term is computed
C     because that term also depends on the model state.
C
C     RNFvflx = RNF_vflx and RNFmflx = RNF_mflx are the dense fields
C     themselves: RNF_EXF_RUNOFF assigns the exf runoff array from
C     RNF_vflx with no arithmetic (rnf_exf_runoff.F:191), so RNFvflx is
C     bitwise EXFroff before any control is added, as the decision-8
C     table requires; RNF_mflx is already rhoConstFresh*RNF_vflx,
C     computed by RNF_LOAD_AT in the same pass that builds RNF_vflx.
C
C     RNFtemp = (mT)/m_T where m_T > 0, and RNFsaln = (mS)/m where
C     m > 0: both divisions RNF_LOAD_AT deliberately leaves undone
C     (its header comment), because the terms that read (mT), (mS) and
C     m directly need the sums and not the ratio. Two kinds of cell
C     have no defined ratio to report:
C
C       o no runoff at all (m = 0, so m_T = 0 too);
C       o runoff present, but every source feeding the cell is
C         missing its temperature in the record used (m_T = 0 with
C         m > 0; RNFsaln is unaffected, since a missing temperature
C         does not make the salinity sum empty).
C
C     This package's own convention for "undefined" (the issue asks
C     for one, deliberate and stated): 0. It is the value RNF_mXT and
C     RNF_mXS themselves already hold at such a cell, so leaving the
C     ratio at 0 rather than writing a special fill value needs no
C     second convention downstream and matches the sign-carrying zero
C     every other RNF diagnostic reports at a cell with no runoff.
C     Nothing reads RNFtemp or RNFsaln as a term (RNF_EXF_RUNOFF and
C     RNF_FORCING_SURF read the sums directly), so this convention
C     reaches no model field.
C
C     RNFnsrc counts the target entries that own cell (i,j): schema
C     rule T03 forbids a duplicate (source, cell) pair in the file
C     (docs/runoff_schema.md), so that count is also the number of
C     distinct sources feeding the cell, with no separate bookkeeping
C     of which source ids have been seen.

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
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
#ifdef ALLOW_DIAGNOSTICS

C !LOCAL VARIABLES: ====================================================
C  bi, bj               :: tile indices
C  i, j, k               :: loop indices
C  tLoc                  :: RNFtemp of this tile  [degC]
C  sLoc                  :: RNFsaln of this tile   [g/kg]
C  nLoc                  :: RNFnsrc of this tile   [1], real for
C                          DIAGNOSTICS_FILL, which takes an _RL array
      INTEGER bi, bj, i, j, k
      _RL     tLoc(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     sLoc(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     nLoc(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
CEOP

      IF ( .NOT.useDiagnostics ) RETURN

      CALL DIAGNOSTICS_FILL( RNF_vflx, 'RNFvflx ', 0, 1, 0,1,1,myThid )
      CALL DIAGNOSTICS_FILL( RNF_mflx, 'RNFmflx ', 0, 1, 0,1,1,myThid )

      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)

        DO j = 1-OLy, sNy+OLy
         DO i = 1-OLx, sNx+OLx
          tLoc(i,j) = 0. _d 0
          sLoc(i,j) = 0. _d 0
          nLoc(i,j) = 0. _d 0
         ENDDO
        ENDDO

        DO j = 1, sNy
         DO i = 1, sNx
          IF ( RNF_mflxT(i,j,bi,bj) .GT. 0. _d 0 ) THEN
            tLoc(i,j) = RNF_mXT(i,j,bi,bj)/RNF_mflxT(i,j,bi,bj)
          ENDIF
          IF ( RNF_mflx(i,j,bi,bj) .GT. 0. _d 0 ) THEN
            sLoc(i,j) = RNF_mXS(i,j,bi,bj)/RNF_mflx(i,j,bi,bj)
          ENDIF
         ENDDO
        ENDDO

C--     One target entry per (source, cell) pair (schema rule T03), so
C       this count of entries owning (i,j) is the number of distinct
C       sources feeding it.
        DO k = 1, RNF_nTgt(bi,bj)
          i = RNF_tgtI(k,bi,bj)
          j = RNF_tgtJ(k,bi,bj)
          nLoc(i,j) = nLoc(i,j) + 1. _d 0
        ENDDO

        CALL DIAGNOSTICS_FILL( tLoc, 'RNFtemp ', 0, 1, 2, bi,bj,myThid)
        CALL DIAGNOSTICS_FILL( sLoc, 'RNFsaln ', 0, 1, 2, bi,bj,myThid)
        CALL DIAGNOSTICS_FILL( nLoc, 'RNFnsrc ', 0, 1, 2, bi,bj,myThid)

       ENDDO
      ENDDO

#endif /* ALLOW_DIAGNOSTICS */
#endif /* ALLOW_RNF */

      RETURN
      END
