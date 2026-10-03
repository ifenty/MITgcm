#include "RNF_OPTIONS.h"

CBOP
C !ROUTINE: RNF_INIT_FIXED

C !INTERFACE: ==========================================================
      SUBROUTINE RNF_INIT_FIXED( myThid )

C !DESCRIPTION:
C     Initialize the fixed part of the RNF package.
C     Package skeleton: it only prints the parameter summary. The
C     static read of RNF_file (source table, target cells, fractions),
C     the placement of the targets on tiles and the global check of
C     the fraction sums come with the reader (RUNOFF-004).
C     This routine is called before RNF_CHECK (S/R INITIALISE_FIXED
C     calls PACKAGES_INIT_FIXED before PACKAGES_CHECK), so what the
C     static read relies on is checked in RNF_READPARMS.

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

C--   Print the package parameters
      CALL RNF_SUMMARY( myThid )

#endif /* ALLOW_RNF */

      RETURN
      END
