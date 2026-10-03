#include "RNF_OPTIONS.h"

CBOP
C !ROUTINE: RNF_INIT_VARIA

C !INTERFACE: ==========================================================
      SUBROUTINE RNF_INIT_VARIA( myThid )

C !DESCRIPTION:
C     Initialize the variable part of the RNF package.
C     Package skeleton: the package has no field yet and this routine
C     does nothing. It will set the dense runoff fields to zero and
C     load the first records (RUNOFF-004, RUNOFF-005).

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "RNF.h"

C !INPUT PARAMETERS: ===================================================
C  myThid               :: my thread Id number
      INTEGER myThid

C !OUTPUT PARAMETERS: ==================================================
C  none
CEOP

#ifdef ALLOW_RNF

#endif /* ALLOW_RNF */

      RETURN
      END
