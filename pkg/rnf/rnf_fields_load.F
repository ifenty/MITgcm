#include "RNF_OPTIONS.h"

CBOP
C !ROUTINE: RNF_FIELDS_LOAD

C !INTERFACE: ==========================================================
      SUBROUTINE RNF_FIELDS_LOAD( myTime, myIter, myThid )

C !DESCRIPTION:
C     Build the dense per-tile runoff fields of the current time step
C     from the per-tile target lists and the source fluxes.
C     Called from LOAD_FIELDS_DRIVER before EXF_GETFORCING, with the
C     same myTime that pkg/exf uses.
C
C     For every target cell c of the tile
C       RNF_vflx(c) = sum_s flux_s * frac_{s,c} / rA(c)     [m/s]
C       RNF_mflx(c) = rhoConstFresh * RNF_vflx(c)    [kg/m^2/s]
C     Each target adds (flux*frac)/rA in real*8, the order in which
C     the Python converter builds a dense record, so that one target
C     on a cell gives back the dense value to within one unit in the
C     last place. The volume flux is kept in m/s and not derived from
C     the mass flux: multiplying and dividing by rhoConstFresh does
C     not give back the same number.
C
C     The flux of a source is read once, in RNF_INIT_VARIA, because
C     only a file with one constant record can be read. Reading the
C     records that bracket myTime and interpolating between them
C     comes with the time handling (RUNOFF-005); the fields are
C     rebuilt here at every step so that nothing else has to change
C     when they start to depend on time.

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
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

C !LOCAL VARIABLES: ====================================================
C  bi, bj               :: tile indices
C  i, j, k              :: loop indices
C  kSrc                 :: local source index of a target entry
      INTEGER bi, bj, i, j, k, kSrc
CEOP

      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)

        DO j = 1-OLy, sNy+OLy
         DO i = 1-OLx, sNx+OLx
          RNF_vflx(i,j,bi,bj) = 0. _d 0
          RNF_mflx(i,j,bi,bj) = 0. _d 0
         ENDDO
        ENDDO

        DO k = 1, RNF_nTgt(bi,bj)
          i = RNF_tgtI(k,bi,bj)
          j = RNF_tgtJ(k,bi,bj)
          kSrc = RNF_tgtSrc(k,bi,bj)
          RNF_vflx(i,j,bi,bj) = RNF_vflx(i,j,bi,bj)
     &      + RNF_srcFlux(kSrc,bi,bj)*RNF_tgtFrac(k,bi,bj)
     &        /rA(i,j,bi,bj)
        ENDDO

        DO j = 1, sNy
         DO i = 1, sNx
          RNF_mflx(i,j,bi,bj) = rhoConstFresh*RNF_vflx(i,j,bi,bj)
         ENDDO
        ENDDO

       ENDDO
      ENDDO

#endif /* ALLOW_RNF */

      RETURN
      END
