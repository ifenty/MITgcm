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
C     flux of the current step.
C     Called from EXF_GETFFIELDS after the reads of the dense runoff
C     and runoff temperature files and before the control variables
C     are added, so that exf uses the field as it does a dense one.
C     Package skeleton: this routine does nothing and leaves "runoff"
C     untouched. It will assign the field over the tile interior
C     (RUNOFF-004).

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
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
CEOP

#ifdef ALLOW_RNF

#endif /* ALLOW_RNF */

      RETURN
      END
