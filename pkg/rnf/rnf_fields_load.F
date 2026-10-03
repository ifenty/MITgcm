#include "RNF_OPTIONS.h"

CBOP
C !ROUTINE: RNF_FIELDS_LOAD

C !INTERFACE: ==========================================================
      SUBROUTINE RNF_FIELDS_LOAD( myTime, myIter, myThid )

C !DESCRIPTION:
C     Load the sparse runoff records for the current time and build
C     the dense per-tile runoff fields.
C     Called from LOAD_FIELDS_DRIVER before EXF_GETFORCING, with the
C     same myTime that pkg/exf uses.
C     Package skeleton: this routine does nothing. It will read the
C     records, interpolate in time and build the dense fields
C     (RUNOFF-004, RUNOFF-005).

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
