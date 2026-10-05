#include "RNF_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP 0
C !ROUTINE: RNF_DIAGNOSTICS_INIT

C !INTERFACE:
      SUBROUTINE RNF_DIAGNOSTICS_INIT( myThid )

C     !DESCRIPTION:
C     Register the diagnostics of the sparse runoff package, called
C     from RNF_INIT_FIXED (package design, decision 1).
C
C     **What is registered here, and what is not.** This is the subset
C     RUNOFF-013 needs: the three tendency terms as they are applied,
C     filled inside RNF_TENDENCY_APPLY_T, _S and _PTR, at the level
C     and with the reference value each term actually used. They are
C     the package's own record of what it did to the model state, and
C     the instrument the cell-by-cell comparison with the exf
C     runoftemp term uses, since neither TFLUX nor SFLUX includes
C     these terms (model/src/diags_oceanic_surf_flux.F builds those
C     from surfaceForcingT/S and PmEpR).
C
C     The input-only fields of decision 8 (RNFvflx, RNFmflx, RNFtemp,
C     RNFsaln, RNFnsrc), the heat and salt terms in flux units
C     (RNFheat, RNFsalt) and the monitor belong to RUNOFF-015, which
C     fills them at the end of RNF_FIELDS_LOAD. Decision 8 lists
C     RNFheat and RNFsalt as two-dimensional; that has to be settled
C     with the fill, because a 2-D diagnostic may be filled only once
C     per step and tile, while the terms below are computed at one
C     level per column and that level is kSurfC(i,j) under an ice
C     shelf. The three registered here are three-dimensional, which
C     is how decision 8 lists RNFgT and RNFgS, so the question does
C     not arise for them.
C
C     The tracer diagnostics are numbered by runoff tracer, not by
C     ptracer: RNFtrNN is the term of the NNth runoff_ptracer_<NAME>
C     variable of the file, whose name and ptracer number
C     RNF_SUMMARY prints. RNF_nTr of them are registered, because the
C     list is built before the file is read.

C     !USES:
      IMPLICIT NONE
#include "EEPARAMS.h"
#include "SIZE.h"
#include "RNF_SIZE.h"
#include "RNF.h"

C     !INPUT/OUTPUT PARAMETERS:
C     myThid ::  my Thread Id number
      INTEGER myThid
CEOP

#ifdef ALLOW_DIAGNOSTICS
C     !LOCAL VARIABLES:
C     diagNum   :: diagnostics number in the (long) list of available diag.
C     diagName  :: local short name (8c) of a diagnostics
C     diagCode  :: local parser field with characteristics of the diagnostics
C     diagUnits :: local string (16c): physical units of a diagnostic field
C     diagTitle :: local string (80c): description of field in diagnostic
      INTEGER       diagNum
      INTEGER       iTr
      CHARACTER*8   diagName
      CHARACTER*16  diagCode
      CHARACTER*16  diagUnits
      CHARACTER*(80) diagTitle

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      diagName  = 'RNFgT   '
      diagTitle = 'Runoff temperature tendency term of pkg/rnf'
      diagUnits = 'degC/s          '
      diagCode  = 'SMR     MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'RNFgS   '
      diagTitle = 'Runoff salinity tendency term of pkg/rnf'
      diagUnits = 'g/kg/s          '
      diagCode  = 'SMR     MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      DO iTr = 1, RNF_nTr
        WRITE(diagName,'(A,I2.2,A)') 'RNFtr', iTr, ' '
        WRITE(diagTitle,'(A,I3,A)')
     &    'Runoff tracer tendency term of pkg/rnf, runoff tracer',
     &    iTr, ' (see RNF_SUMMARY for its name)'
        diagUnits = 'tracer units/s  '
        diagCode  = 'SMR     MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      ENDDO

#endif /* ALLOW_DIAGNOSTICS */

      RETURN
      END
