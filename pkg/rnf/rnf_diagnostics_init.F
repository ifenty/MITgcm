#include "RNF_OPTIONS.h"
#ifdef ALLOW_EXF
# include "EXF_OPTIONS.h"
#endif /* ALLOW_EXF */

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP 0
C !ROUTINE: RNF_DIAGNOSTICS_INIT

C !INTERFACE:
      SUBROUTINE RNF_DIAGNOSTICS_INIT( myThid )

C     !DESCRIPTION:
C     Register the diagnostics of the sparse runoff package, called
C     from RNF_INIT_FIXED (package design, decision 1).
C
C     **What is registered here, and what is not.** The runoff's own
C     contributions to the surface flux fields it feeds, one per
C     property, each in the units and sign of the field it is part of
C     (RUNOFF-031), so that a budget can separate the runoff from the
C     rest of TFLUX, SFLUX and the ptracers surface forcing:
C
C       RNFqnet  :: the heat EXF_MAPFIELDS adds to Qnet for the runoff
C                   (RNF_EXF_RUNOFF), W/m^2, > 0 decreases theta, as
C                   EXFqnet and Qnet; registered and filled only with
C                   ALLOW_RUNOFTEMP, which the heat needs
C       RNFsflx  :: the salt added to surfaceForcingS
C                   (RNF_FORCING_SURF), g/m^2/s, > 0 increases
C                   salinity, as SFLUX
C       RNFtfNN  :: the tracer added to surfaceForcingPTr
C                   (RNF_FORCING_SURF_PTR), tracer units * kg/m^2/s,
C                   > 0 increases the tracer
C
C     Each is filled where the term is computed, from the values that
C     term used, once per step. All three are two-dimensional: the
C     surface fields are, and the model applies them at the surface
C     level of the column. They replace the three-dimensional tendency
C     diagnostics RNFgT, RNFgS and RNFtrNN of RUNOFF-013, which ended
C     with the tendency terms themselves.
C
C     The input-only fields of decision 8 (RNFvflx, RNFmflx, RNFtemp,
C     RNFsaln, RNFnsrc) and the monitor belong to RUNOFF-015.
C
C     The tracer diagnostics are numbered by runoff tracer, not by
C     ptracer: RNFtfNN is the term of the NNth runoff_ptracer_<NAME>
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

#ifdef ALLOW_RUNOFTEMP
      diagName  = 'RNFqnet '
      diagTitle = 'Runoff heat added to Qnet by pkg/rnf,'
     &          //' >0 decreases theta'
      diagUnits = 'W/m^2           '
      diagCode  = 'SM      U1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
#endif /* ALLOW_RUNOFTEMP */

      diagName  = 'RNFsflx '
      diagTitle = 'Runoff salt flux of pkg/rnf in surfaceForcingS,'
     &          //' >0 increases salinity'
      diagUnits = 'g/m^2/s         '
      diagCode  = 'SM      U1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      DO iTr = 1, RNF_nTr
        WRITE(diagName,'(A,I2.2,A)') 'RNFtf', iTr, ' '
        WRITE(diagTitle,'(A,I3,A)')
     &    'pkg/rnf surface flux of runoff tracer',
     &    iTr, ', named by RNF_SUMMARY, >0 increases it'
        diagUnits = 'kg/m^2/s.trc    '
        diagCode  = 'SM      U1      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      ENDDO

#endif /* ALLOW_DIAGNOSTICS */

      RETURN
      END
