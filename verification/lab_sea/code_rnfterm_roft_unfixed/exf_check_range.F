#include "EXF_OPTIONS.h"

      SUBROUTINE EXF_CHECK_RANGE( myTime, myIter, myThid )

C     ==================================================================
C     SUBROUTINE EXF_CHECK_RANGE
C     ==================================================================
C     Warn about every input field that is outside the range it can
C     physically take on a wet cell, then stop the run. Called from
C     EXF_GETFORCING (exf_getforcing.F:346-349) when useExfCheckRange
C     (default .TRUE., exf_readparms.F:307) and either the iteration is
C     nIter0 or exf_debugLev is debLevC or above: i.e. at the first
C     step of a run, and at every step only under that debug level.
C     Two tests are treated conditionally, and both only when useRNF is
C     true, because pkg/rnf feeds the runoff field from sparse point
C     sources whose per-cell volume flux is far above what either bound
C     was written for, and bounds those sources itself:
C      o the runoff UPPER bound is skipped;
C      o the sflux bound is applied to sflux + runoff, i.e. to what
C        sflux is without the runoff EXF_GETFORCING subtracted into it,
C        so it still refuses an out-of-range evap - precip.
C     Everything else is unconditional, the negative-runoff test and
C     every other field's range included, and nothing changes at all
C     with useRNF false, the dense runoffFile path included. See the
C     freshwater-flux and runoff blocks below.
C     ==================================================================

      IMPLICIT NONE

C     == global variables ==
#include "EEPARAMS.h"
#include "SIZE.h"
#include "PARAMS.h"
#include "GRID.h"

#include "EXF_PARAM.h"
c#include "EXF_CONSTANTS.h"
#include "EXF_FIELDS.h"

C     == routine arguments ==
C     myThid - thread number for this instance of the routine.
      _RL myTime
      INTEGER myIter, myThid

C     == local variables ==
C     sfluxLoc - the freshwater flux this routine bounds: sflux itself,
C                or sflux with the runoff added back when useRNF
      INTEGER i, j, ks, bi, bj
      INTEGER exferr
      _RL sfluxLoc

C     == end of interface ==

      exferr = 0

C--   Only master thread can safely write directly to standard output:
      _BARRIER
      _BEGIN_MASTER( myThid )

      ks = 1
      IF ( usingPCoords ) ks = Nr
c     DO bj = myByLo(myThid), myByHi(myThid)
c      DO bi = myBxLo(myThid), myBxHi(myThid)
      DO bj = 1, nSy
       DO bi = 1, nSx

C Change checking range because some atmospheric fields will
C not always have valid values in the tile edges.
c       DO j = 1-OLy, sNy+OLy
c        DO i = 1-OLx, sNx+OLx
        DO j = 1, sNy
         DO i = 1, sNx

C     Heat flux.
          IF ( ( hflux(i,j,bi,bj) .GT. 1600. .OR.
     &         hflux(i,j,bi,bj) .LT. -500. ) .AND.
     &         maskC(i,j,ks,bi,bj) .NE. 0. ) THEN
           WRITE(standardMessageUnit,'(A,5(1X,I6),2X,D22.15)')
     &          'EXF WARNING: hflux out of range for bi,bj,i,j,it= ',
     &          bi, bj, i, j, myIter, hflux(i,j,bi,bj)
           exferr = 1
          ENDIF

C     Freshwater flux. With useRNF the runoff is added back before the
C     test: EXF_GETFORCING sets sflux = sflux - runoff
C     (exf_getforcing.F:313, inside #ifdef ALLOW_RUNOFF) just before
C     calling this routine, and one sparse point source is far above
C     1.E-6 m/s on a fine grid, so testing sflux itself would refuse
C     the run for its runoff volume rather than for the freshwater
C     flux this bound is about - and would do so however the runoff
C     test below is written. What is tested with useRNF is therefore
C     evap - precip, which is the part of sflux this check exists for
C     and which stays guarded: a run with no runoff and an
C     out-of-range evap - precip still stops.
C     Without useRNF sfluxLoc IS sflux and nothing changes, the dense
C     runoffFile path included. The restore sits inside ALLOW_RUNOFF
C     because the runoff array exists only under that option, while
C     this test does not.
          sfluxLoc = sflux(i,j,bi,bj)
#ifdef ALLOW_RUNOFF
          IF ( useRNF ) sfluxLoc = sfluxLoc + runoff(i,j,bi,bj)
#endif
          IF ( ABS(sfluxLoc) .GT. 1.E-6 .AND.
     &         maskC(i,j,ks,bi,bj) .NE. 0. ) THEN
           WRITE(standardMessageUnit,'(A,5(1X,I6),2X,D22.15)')
     &          'EXF WARNING: sflux out of range for bi,bj,i,j,it= ',
     &          bi, bj, i, j, myIter, sfluxLoc
           exferr = 1
          ENDIF

C     Zonal wind stress.
          IF ( ABS(ustress(i,j,bi,bj)) .GT. 2.7 .AND.
     &         maskW(i,j,ks,bi,bj) .NE. 0. ) THEN
           WRITE(standardMessageUnit,'(A,5(1X,I6),2X,D22.15)')
     &          'EXF WARNING: ustress out of range for bi,bj,i,j,it= ',
     &          bi, bj, i, j, myIter, ustress(i,j,bi,bj)
           exferr = 1
          ENDIF

C     Meridional wind stress.
          IF ( ABS(vstress(i,j,bi,bj)) .GT. 2.3 .AND.
     &         maskS(i,j,ks,bi,bj) .NE. 0. ) THEN
           WRITE(standardMessageUnit,'(A,5(1X,I6),2X,D22.15)')
     &          'EXF WARNING: vstress out of range for bi,bj,i,j,it= ',
     &          bi, bj, i, j, myIter, vstress(i,j,bi,bj)
           exferr = 1
          ENDIF

          IF ( useAtmWind ) THEN
C     zonal wind speed
           IF ( ABS(uwind(i,j,bi,bj)) .GT. 100. .AND.
     &         maskW(i,j,ks,bi,bj) .NE. 0. ) THEN
            WRITE(standardMessageUnit,'(A,5(1X,I6),2X,D22.15)')
     &          'EXF WARNING: uwind out of range for bi,bj,i,j,it= ',
     &          bi, bj, i, j, myIter, uwind(i,j,bi,bj)
            exferr = 1
           ENDIF

C     zonal wind speed
           IF ( ABS(vwind(i,j,bi,bj)) .GT. 100. .AND.
     &         maskS(i,j,ks,bi,bj) .NE. 0. ) THEN
            WRITE(standardMessageUnit,'(A,5(1X,I6),2X,D22.15)')
     &          'EXF WARNING: vwind out of range for bi,bj,i,j,it= ',
     &          bi, bj, i, j, myIter, vwind(i,j,bi,bj)
            exferr = 1
           ENDIF
          ENDIF

C     wind speed modulus
          IF ( ( wspeed(i,j,bi,bj) .LT. 0. .OR.
     &         wspeed(i,j,bi,bj) .GT. 100. ) .AND.
     &         maskS(i,j,ks,bi,bj) .NE. 0. ) THEN
           WRITE(standardMessageUnit,'(A,5(1X,I6),2X,D22.15)')
     &          'EXF WARNING: wspeed out of range for bi,bj,i,j,it= ',
     &          bi, bj, i, j, myIter, wspeed(i,j,bi,bj)
           exferr = 1
          ENDIF

#ifdef ALLOW_ATM_TEMP
C     2-m air temperature
          IF ( (atemp(i,j,bi,bj) .LT. 183 .OR.
     &         atemp(i,j,bi,bj) .GT. 343 ) .AND.
     &         maskC(i,j,ks,bi,bj) .NE. 0. ) THEN
           WRITE(standardMessageUnit,'(2A,5(1X,I6),2X,D22.15)')
     &          'EXF WARNING: atemp + exf_offset_atemp ',
     &          'out of range for bi,bj,i,j,it= ',
     &          bi, bj, i, j, myIter, atemp(i,j,bi,bj)
           exferr = 1
          ENDIF

C     2-m specific humidity
          IF ( (aqh(i,j,bi,bj) .LT. 0. .OR.
     &         aqh(i,j,bi,bj) .GT. 0.1 ) .AND.
     &         maskC(i,j,ks,bi,bj) .NE. 0. ) THEN
           WRITE(standardMessageUnit,'(A,5(1X,I6),2X,D22.15)')
     &          'EXF WARNING: aqh out of range for bi,bj,i,j,it= ',
     &          bi, bj, i, j, myIter, aqh(i,j,bi,bj)
           exferr = 1
          ENDIF

C     precipitation rate
          IF ( (precip(i,j,bi,bj) .LT. 0. .OR.
     &         precip(i,j,bi,bj) .GT. 2.E-6 ) .AND.
     &         maskC(i,j,ks,bi,bj) .NE. 0. ) THEN
           WRITE(standardMessageUnit,'(A,5(1X,I6),2X,D22.15)')
     &          'EXF WARNING: precip out of range for bi,bj,i,j,it= ',
     &          bi, bj, i, j, myIter, precip(i,j,bi,bj)
           exferr = 1
          ENDIF

C     snow
          IF ( (snowprecip(i,j,bi,bj) .LT. 0. .OR.
     &         snowprecip(i,j,bi,bj) .GT. 2.E-6 ) .AND.
     &         maskC(i,j,ks,bi,bj) .NE. 0. ) THEN
           WRITE(standardMessageUnit,'(2A,5(1X,I6),2X,D22.15)')
     &          'EXF WARNING: snowprecip out of range ',
     &          'for bi,bj,i,j,it= ',
     &          bi, bj, i, j, myIter, snowprecip(i,j,bi,bj)
           exferr = 1
          ENDIF
#endif

#ifdef SHORTWAVE_HEATING
C     Short wave radiative flux.
          IF ( (swflux(i,j,bi,bj) .GT. 1. .OR.
     &         swflux(i,j,bi,bj) .LT. -1000. ) .AND.
     &         maskC(i,j,ks,bi,bj) .NE. 0. ) THEN
           WRITE(standardMessageUnit,'(A,5(1X,I6),2X,D22.15)')
     &          'EXF WARNING: swflux out of range for bi,bj,i,j,it= ',
     &          bi, bj, i, j, myIter, swflux(i,j,bi,bj)
           exferr = 1
          ENDIF
#endif

#ifdef ALLOW_RUNOFF
C     Runoff. The UPPER bound is skipped when pkg/rnf is in use: a
C     sparse point source puts one whole river into one cell, which on
C     a fine grid is far above 1.E-6 m/s (1000 m3/s into a single 2 km
C     cell is 2.5E-4 m/s), so this bound stops every realistic
C     point-source run at its first step. pkg/rnf bounds the volume
C     flux of each source as it reads it instead (RNF_srcFluxMax of
C     pkg/rnf/RNF.h, enforced in RNF_NC_READ_ONE) and says so in its
C     summary. The negative test below is NOT skipped: it holds with
C     useRNF true as well as false.
C     This skip is not sufficient on its own, which is why the sflux
C     test above restores the runoff when useRNF: EXF_GETFORCING sets
C     sflux = sflux - runoff (exf_getforcing.F:313) before calling this
C     routine (:346-349), so without that restore any wet cell whose
C     runoff is over 1.E-6 m/s would still be refused there, whatever
C     this test does. Both relaxations are needed, and both are
C     guarded by useRNF alone.
C     A DENSE runoffFile above 1.E-6 m/s is still refused twice over,
C     here and by the sflux test, and still needs
C     useExfCheckRange=.FALSE. That is deliberate: the dense path is
C     untouched by design, and the wider upstream question belongs to
C     an upstream discussion.
          IF ( (runoff(i,j,bi,bj) .LT. 0. .OR.
     &         ( runoff(i,j,bi,bj) .GT. 1.E-6
     &           .AND. .NOT.useRNF ) ) .AND.
     &         maskC(i,j,ks,bi,bj) .NE. 0. ) THEN
           WRITE(standardMessageUnit,'(A,5(1X,I6),2X,D22.15)')
     &          'EXF WARNING: runoff out of range for bi,bj,i,j,it= ',
     &          bi, bj, i, j, myIter, runoff(i,j,bi,bj)
C     The m/yr advice belongs to the dense runoffFile, which
C     exf_inscal_runoff scales when EXF_SET_FLD reads it. With useRNF
C     the only trigger left here is a NEGATIVE value, and that
C     parameter does not touch the sparse field at all, so the advice
C     would name an unrelated setting. The text itself is unchanged
C     for the dense path.
           IF ( .NOT.useRNF ) THEN
            WRITE(standardMessageUnit,'(A)')
     &          'Please note that input units for runoff are'
            WRITE(standardMessageUnit,'(A)')
     &          'm/s not m/yr.  If input file is in m/yr, set'
            WRITE(standardMessageUnit,'(A)')
     &          'exf_inscal_runoff=3.170979198E-8'
            WRITE(standardMessageUnit,'(A)')
     &          'in the data.exf input file.'
           ENDIF
           exferr = 1
          ENDIF
# ifdef ALLOW_RUNOFTEMP
C     Runoff temperature.
          IF ( (runoftemp(i,j,bi,bj) .LT. -2. .OR.
     &         runoff(i,j,bi,bj) .GT. 36 ) .AND.
     &         maskC(i,j,ks,bi,bj) .NE. 0. ) THEN
           WRITE(standardMessageUnit,'(A,5(1X,I6),2X,D22.15)')
     &          'EXF WARNING: runoftemp out of range at bi,bj,i,j,it= ',
     &          bi, bj, i, j, myIter, runoff(i,j,bi,bj)
           exferr = 1
          ENDIF
# endif /* ALLOW_RUNOFTEMP */
#endif /* ALLOW_RUNOFF */
         ENDDO
        ENDDO

       ENDDO
      ENDDO

      IF ( exferr .NE. 0 ) THEN
       WRITE(standardMessageUnit,'(A)')
     &      'EXF WARNING: If you think these values are OK '
       WRITE(standardMessageUnit,'(A)')
     &      'EXF WARNING: then set useExfCheckRange=.FALSE.'
       STOP 'ABNORMAL END: S/R EXF_CHECK_RANGE'
      ENDIF
      _END_MASTER( myThid )
      _BARRIER

      RETURN
      END
