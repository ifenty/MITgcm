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
C     flux of the current step (package design, decision 2).
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

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
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
      INTEGER bi, bj, i, j
CEOP

      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)
        DO j = 1, sNy
         DO i = 1, sNx
          runoff(i,j,bi,bj) = RNF_vflx(i,j,bi,bj)
         ENDDO
        ENDDO
       ENDDO
      ENDDO

#endif /* ALLOW_RUNOFF */
#endif /* ALLOW_EXF */
#endif /* ALLOW_RNF */

      RETURN
      END
