#include "RNF_OPTIONS.h"

CBOP
C !ROUTINE: RNF_INIT_VARIA

C !INTERFACE: ==========================================================
      SUBROUTINE RNF_INIT_VARIA( myThid )

C !DESCRIPTION:
C     Initialize the variable part of the RNF package: read the one
C     constant record of runoff_flux, report the flux it holds and
C     build the dense fields.
C
C     Two sums are printed: the flux summed over the sources of the
C     file, and the flux summed over the target entries of all tiles
C     and processes, sum_s sum_c flux_s*frac_{s,c}. The second is the
C     volume the model will apply, so the pair is the initialization
C     check of volume conservation. They agree to the fraction
C     tolerance whenever the fraction sums of RNF_INIT_FIXED did,
C     which is what makes the pair a check of the placement: a target
C     that no tile owns, or one owned twice, shows up here.
C     The routine prints their relative difference as well and warns
C     when it exceeds RNF_fracTol; it does not stop, because
C     RNF_INIT_FIXED has already refused a fraction sum that is not 1
C     and a difference here would mean an internal inconsistency
C     rather than bad input. The value is asserted by the control
C     cases of tests/rnf/refusal_check.py, which require the printed
C     file sum to equal the sum computed from the file within 1e-12
C     relative and the target sum within 1e-6.
C
C     Only a file with one constant record can be read. Reading the
C     records that bracket the model time comes with the time
C     handling (RUNOFF-005); this routine then loads the first
C     records instead. The file is opened and closed by
C     RNF_NC_READ_FLUX, so this routine may be called more than once.

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "RNF_SIZE.h"
#include "RNF.h"

C !INPUT PARAMETERS: ===================================================
C  myThid               :: my thread Id number
      INTEGER myThid

C !OUTPUT PARAMETERS: ==================================================
C  none

#ifdef ALLOW_RNF

C !LOCAL VARIABLES: ====================================================
C  fluxSrc              :: flux summed over the sources of the file
C  fluxTgt              :: flux summed over the targets of all tiles
      INTEGER bi, bj, i, j, k, kSrc
      _RL     fluxSrc, fluxTgt, fluxRef
      CHARACTER*(MAX_LEN_MBUF) msgBuf
CEOP

      fluxSrc = 0. _d 0
      fluxTgt = 0. _d 0

      _BEGIN_MASTER(myThid)

C--   Dense fields start at zero, also in the halo region, which
C     RNF_EXF_RUNOFF does not touch
      DO bj = 1, nSy
       DO bi = 1, nSx
        DO j = 1-OLy, sNy+OLy
         DO i = 1-OLx, sNx+OLx
          RNF_vflx(i,j,bi,bj) = 0. _d 0
          RNF_mflx(i,j,bi,bj) = 0. _d 0
         ENDDO
        ENDDO
        RNF_tileSum(bi,bj) = 0. _d 0
       ENDDO
      ENDDO

C--   The one record of the file, by the master thread
      CALL RNF_NC_READ_FLUX( 1,
     O     fluxSrc,
     I     myThid )

      DO bj = 1, nSy
       DO bi = 1, nSx
        DO k = 1, RNF_nTgt(bi,bj)
          kSrc = RNF_tgtSrc(k,bi,bj)
          RNF_tileSum(bi,bj) = RNF_tileSum(bi,bj)
     &      + RNF_srcFlux(kSrc,bi,bj)*RNF_tgtFrac(k,bi,bj)
        ENDDO
       ENDDO
      ENDDO

      _END_MASTER(myThid)
      _BARRIER

C--   The target sum runs over the tiles of every process
      CALL GLOBAL_SUM_TILE_RL( RNF_tileSum, fluxTgt, myThid )

      _BEGIN_MASTER(myThid)
      WRITE(msgBuf,'(2A,1PE20.12)') 'RNF_INIT_VARIA: runoff flux,',
     &  ' sum over the sources of the file [m^3/s] =', fluxSrc
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(2A,1PE20.12)') 'RNF_INIT_VARIA: runoff flux,',
     &  ' sum over the targets of all tiles [m^3/s] =', fluxTgt
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )

C--   The difference of the two, which is the volume the placement
C     would lose or double. Both sums are zero in a run whose fluxes
C     are all zero, which is no error.
      fluxRef = MAX( ABS(fluxSrc), ABS(fluxTgt) )
      IF ( fluxRef .GT. 0. _d 0 ) THEN
        WRITE(msgBuf,'(2A,1PE16.8)') 'RNF_INIT_VARIA: runoff flux,',
     &    ' relative difference of the two sums =',
     &    ( fluxTgt - fluxSrc )/fluxRef
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                      SQUEEZE_RIGHT, myThid )
        IF ( ABS( fluxTgt - fluxSrc ) .GT. RNF_fracTol*fluxRef ) THEN
          WRITE(msgBuf,'(2A)') '** WARNING ** RNF_INIT_VARIA: the',
     &      ' applied volume differs from the flux of the file by'
          CALL PRINT_ERROR( msgBuf, myThid )
          WRITE(msgBuf,'(2A)') '** WARNING ** RNF_INIT_VARIA: more',
     &      ' than RNF_fracTol; the target table loses volume'
          CALL PRINT_ERROR( msgBuf, myThid )
        ENDIF
      ENDIF
      _END_MASTER(myThid)
      _BARRIER

C--   Build the dense fields, so that they are set whichever
C     routine of the model asks for them first
      CALL RNF_FIELDS_LOAD( startTime, nIter0, myThid )

#endif /* ALLOW_RNF */

      RETURN
      END
