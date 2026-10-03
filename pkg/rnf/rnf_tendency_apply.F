#include "RNF_OPTIONS.h"
C--  File rnf_tendency_apply.F: Routines to apply RNF tendencies
C--   Contents
C--   o RNF_TENDENCY_APPLY_T
C--   o RNF_TENDENCY_APPLY_S
C--   o RNF_TENDENCY_APPLY_PTR

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP 0
C !ROUTINE: RNF_TENDENCY_APPLY_T

C !INTERFACE:
      SUBROUTINE RNF_TENDENCY_APPLY_T(
     U                    gT_arr,
     I                    iMin,iMax,jMin,jMax, k, bi, bj,
     I                    myTime, myIter, myThid )

C     !DESCRIPTION:
C     Add the runoff temperature term to the T (theta) tendency.
C     Called from APPLY_FORCING_T; works for one level at a time.
C     Package skeleton: this routine does nothing and leaves gT_arr
C     untouched. The term comes with RUNOFF-013.

C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "RNF.h"

C     !INPUT/OUTPUT PARAMETERS:
C     gT_arr    :: the tendency array
C     iMin,iMax :: working range of tile for applying forcing
C     jMin,jMax :: working range of tile for applying forcing
C     k         :: Current vertical level index
C     bi,bj     :: Current tile indices
C     myTime    :: Current time in simulation
C     myIter    :: Current iteration number
C     myThid    :: my Thread Id number
      _RL     gT_arr(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      INTEGER iMin, iMax, jMin, jMax
      INTEGER k, bi, bj
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
CEOP

#ifdef ALLOW_RNF

#endif /* ALLOW_RNF */

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP 0
C !ROUTINE: RNF_TENDENCY_APPLY_S

C !INTERFACE:
      SUBROUTINE RNF_TENDENCY_APPLY_S(
     U                    gS_arr,
     I                    iMin,iMax,jMin,jMax, k, bi, bj,
     I                    myTime, myIter, myThid )

C     !DESCRIPTION:
C     Add the runoff salinity term to the S tendency.
C     Called from APPLY_FORCING_S; works for one level at a time.
C     Package skeleton: this routine does nothing and leaves gS_arr
C     untouched. The term comes with RUNOFF-013.

C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "RNF.h"

C     !INPUT/OUTPUT PARAMETERS:
C     gS_arr    :: the tendency array
C     iMin,iMax :: working range of tile for applying forcing
C     jMin,jMax :: working range of tile for applying forcing
C     k         :: Current vertical level index
C     bi,bj     :: Current tile indices
C     myTime    :: Current time in simulation
C     myIter    :: Current iteration number
C     myThid    :: my Thread Id number
      _RL     gS_arr(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      INTEGER iMin, iMax, jMin, jMax
      INTEGER k, bi, bj
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
CEOP

#ifdef ALLOW_RNF

#endif /* ALLOW_RNF */

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP 0
C !ROUTINE: RNF_TENDENCY_APPLY_PTR

C !INTERFACE:
      SUBROUTINE RNF_TENDENCY_APPLY_PTR(
     U                    gPtr_arr,
     I                    iMin,iMax,jMin,jMax, k, bi, bj,
     I                    iTracer, myTime, myIter, myThid )

C     !DESCRIPTION:
C     Add the runoff tracer term to the tendency of passive tracer
C     iTracer.
C     Called from PTRACERS_APPLY_FORCING; works for one level and one
C     tracer at a time.
C     Package skeleton: this routine does nothing and leaves gPtr_arr
C     untouched. The term comes with RUNOFF-013.

C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "RNF.h"

C     !INPUT/OUTPUT PARAMETERS:
C     gPtr_arr  :: the tendency array
C     iMin,iMax :: working range of tile for applying forcing
C     jMin,jMax :: working range of tile for applying forcing
C     k         :: Current vertical level index
C     bi,bj     :: Current tile indices
C     iTracer   :: passive tracer number
C     myTime    :: Current time in simulation
C     myIter    :: Current iteration number
C     myThid    :: my Thread Id number
      _RL     gPtr_arr(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      INTEGER iMin, iMax, jMin, jMax
      INTEGER k, bi, bj
      INTEGER iTracer
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
CEOP

#ifdef ALLOW_RNF

#endif /* ALLOW_RNF */

      RETURN
      END
