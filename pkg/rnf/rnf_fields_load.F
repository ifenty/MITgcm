#include "RNF_OPTIONS.h"

CBOP
C !ROUTINE: RNF_FIELDS_LOAD

C !INTERFACE: ==========================================================
      SUBROUTINE RNF_FIELDS_LOAD( myTime, myIter, myThid )

C !DESCRIPTION:
C     Put the sparse runoff of the current time step into the dense
C     per-tile fields: select the time records, read the ones that are
C     not already in a buffer, combine them into the per-tile source
C     fluxes, and spread those over the target cells.
C     Called from LOAD_FIELDS_DRIVER before EXF_GETFORCING, with the
C     same myTime that pkg/exf uses.
C
C     **Records.** RNF_GETREC returns the two records that bracket
C     myTime, the weight fac of the earlier one and the file year each
C     comes from. The two record buffers of RNF.h are tagged with what
C     they hold, so a record is read when, and only when, the bracket
C     moves onto a record no buffer has; a buffer whose record the
C     other one already holds is copied instead of re-read. Only a
C     record with a non-zero weight is needed, so
C
C       fac = 1 (on a record time, or hold-exact): only the earlier,
C       fac = 0 (hold-exact): only the later,
C       otherwise:            both.
C
C     Skipping the unused record is what pkg/exf cannot do (it keeps
C     two tiled arrays and swaps them), and it changes no value: the
C     weight of the record not read is zero. It also keeps a run from
C     being refused for a record beyond the end of a series that it
C     never needed.
C
C     **Interpolation.** The per-tile source flux is
C
C       flux_s = fac*flux_s(count0) + (1-fac)*flux_s(count1)  [m^3/s]
C
C     which is the form pkg/exf interpolates its dense fields with
C     (pkg/exf/exf_set_fld.F), applied to the source vector instead of
C     to a 2-D field. With RNF_holdRecord the weight is exactly 1 or
C     exactly 0 and the selected record is copied, with no arithmetic,
C     so the applied flux is a value of the file and not a combination
C     of two.
C
C     **Dense fields.** For every target cell c of the tile
C       RNF_vflx(c) = sum_s flux_s * frac_{s,c} / rA(c)     [m/s]
C       RNF_mflx(c) = rhoConstFresh * RNF_vflx(c)    [kg/m^2/s]
C     Each target adds (flux*frac)/rA in real*8, in the order the
C     target table of the file lists them, which is the order the
C     Python converter builds a dense record in, so that one target on
C     a cell gives back the dense value to within one unit in the last
C     place. **That order is a contract:** the applied-field oracle
C     (tests/rnf/applied_field_check.py) rebuilds this sum in the
C     file's table order and compares it bitwise, and
C     tests/rnf/refusal_check.py's order_sensitive_sum case is built
C     so that reversing the order changes the result by one unit in
C     the last place. Reordering the per-tile lists of RNF_INIT_FIXED
C     would void both. The volume flux is kept in m/s and not derived
C     from the mass flux: multiplying and dividing by rhoConstFresh
C     does not give back the same number.

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "RNF_SIZE.h"
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

#ifdef ALLOW_RNF

C !LOCAL VARIABLES: ====================================================
C  bi, bj               :: tile indices
C  i, j, k              :: loop indices
C  kSrc                 :: local source index of a target entry
C  fac                  :: weight of record count0
C  count0, count1       :: the two records that bracket myTime
C  year0, year1         :: the file year each one comes from
C  need0, need1         :: whether that record has a non-zero weight
      INTEGER bi, bj, i, j, k, kSrc
      _RL     fac
      LOGICAL first, changed
      INTEGER count0, count1, year0, year1
      LOGICAL need0, need1
      CHARACTER*(MAX_LEN_MBUF) msgBuf
CEOP

C--   1. Which records, and with which weights
      CALL RNF_GETREC(
     I     myTime, myIter,
     O     fac, first, changed,
     O     count0, count1, year0, year1,
     I     myThid )

      need0 = fac .GT. 0. _d 0
      need1 = fac .LT. 1. _d 0

C--   2. The records themselves, and the source flux of this time.
C     The master thread does the NetCDF reads (model contract) and
C     fills the per-tile source vector of every tile of the process.
      _BEGIN_MASTER(myThid)

      IF ( need0 ) CALL RNF_LOAD_REC( 1, count0, year0, myThid )
      IF ( need1 ) CALL RNF_LOAD_REC( 2, count1, year1, myThid )

      DO bj = 1, nSy
       DO bi = 1, nSx
        DO k = 1, RNF_nSrc(bi,bj)
          IF     ( .NOT.need1 ) THEN
            RNF_srcFlux(k,bi,bj) = RNF_bufFlux(k,bi,bj,1)
          ELSEIF ( .NOT.need0 ) THEN
            RNF_srcFlux(k,bi,bj) = RNF_bufFlux(k,bi,bj,2)
          ELSE
            RNF_srcFlux(k,bi,bj) =
     &             fac             *RNF_bufFlux(k,bi,bj,1)
     &         + ( 1. _d 0 - fac )*RNF_bufFlux(k,bi,bj,2)
          ENDIF
        ENDDO
       ENDDO
      ENDDO
      IF     ( .NOT.need1 ) THEN
        RNF_fluxFile = RNF_bufSum(1)
      ELSEIF ( .NOT.need0 ) THEN
        RNF_fluxFile = RNF_bufSum(2)
      ELSE
        RNF_fluxFile = fac*RNF_bufSum(1)
     &               + ( 1. _d 0 - fac )*RNF_bufSum(2)
      ENDIF

C--   Trace of the record selection, which is what the sparse timing
C     check reads (tests/runoff/lab_sea_runoff_timing_check.py): at
C     debLevC and above every step, below that only when the bracket
C     moves, so that a long run does not fill its log.
C     The two reals are printed with 17 significant digits, which is
C     what a real*8 needs to be read back as the same number: the
C     applied-field oracle rebuilds fac*flux0 + (1-fac)*flux1 from this
C     line and compares the result bitwise, so a shorter field would
C     make it disagree by round-off (measured 2e-15 with 14 digits).
      IF ( RNF_debugLev.GE.debLevC .OR.
     &     ( RNF_debugLev.GE.debLevB .AND.
     &       ( first .OR. changed ) ) ) THEN
        WRITE(msgBuf,'(A,I10,A,I6,A,I6,A,I6,A,I6,A,1PE25.16,A,E25.16)')
     &    ' RNF_FIELDS_LOAD: it=', myIter,
     &    ', rec0=', count0, ', yr0=', year0,
     &    ', rec1=', count1, ', yr1=', year1,
     &    ', t=', myTime, ', fac=', fac
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                      SQUEEZE_RIGHT, myThid )
      ENDIF

      _END_MASTER(myThid)
C--   Every thread waits for the source fluxes of this step
      _BARRIER

C--   3. The dense fields of the tiles of this thread
      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)

        DO j = 1-OLy, sNy+OLy
         DO i = 1-OLx, sNx+OLx
          RNF_vflx(i,j,bi,bj) = 0. _d 0
          RNF_mflx(i,j,bi,bj) = 0. _d 0
         ENDDO
        ENDDO

        DO k = 1, RNF_nTgt(bi,bj)
          i = RNF_tgtI(k,bi,bj)
          j = RNF_tgtJ(k,bi,bj)
          kSrc = RNF_tgtSrc(k,bi,bj)
          RNF_vflx(i,j,bi,bj) = RNF_vflx(i,j,bi,bj)
     &      + RNF_srcFlux(kSrc,bi,bj)*RNF_tgtFrac(k,bi,bj)
     &        /rA(i,j,bi,bj)
        ENDDO

        DO j = 1, sNy
         DO i = 1, sNx
          RNF_mflx(i,j,bi,bj) = rhoConstFresh*RNF_vflx(i,j,bi,bj)
         ENDDO
        ENDDO

       ENDDO
      ENDDO

#endif /* ALLOW_RNF */

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: RNF_LOAD_REC

C !INTERFACE: ==========================================================
      SUBROUTINE RNF_LOAD_REC( iBuf, iRec, iYear, myThid )

C !DESCRIPTION:
C     Make record buffer iBuf (1 or 2) hold record iRec of the file of
C     year iYear, reading it only if neither buffer has it already.
C     Each buffer is tagged with the record and year it holds
C     (RNF_bufRec, RNF_bufYr of RNF.h; -1 means empty), so this
C     routine replaces the record swap of pkg/exf
C     (pkg/exf/exf_SwapFFields) without relying on the new bracket
C     overlapping the old one. That matters at the two places where it
C     does not overlap in the obvious way: the wrap of a repeat cycle,
C     and the first step of a run or a restart, where both records are
C     new and nothing is carried over (package design, decision 8: no
C     pickup).
C
C     Called by the master thread, which does the NetCDF reads.

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "RNF_SIZE.h"
#include "RNF.h"

C !INPUT PARAMETERS: ===================================================
C  iBuf                 :: buffer to fill (1 or 2)
C  iRec                 :: record wanted, counted from 1
C  iYear                :: file year it belongs to (0: a single file)
C  myThid               :: my thread Id number
      INTEGER iBuf
      INTEGER iRec
      INTEGER iYear
      INTEGER myThid

C !OUTPUT PARAMETERS: ==================================================
C  none

#ifdef ALLOW_RNF

C !LOCAL VARIABLES: ====================================================
C  jBuf                 :: the other buffer
      INTEGER jBuf, bi, bj, k
      CHARACTER*(MAX_LEN_FNAM) fileName
CEOP

      IF ( RNF_bufRec(iBuf).EQ.iRec .AND.
     &     RNF_bufYr(iBuf) .EQ.iYear ) RETURN

      jBuf = 3 - iBuf
      IF ( RNF_bufRec(jBuf).EQ.iRec .AND.
     &     RNF_bufYr(jBuf) .EQ.iYear ) THEN
C--     the other buffer has it: copy instead of reading again
        DO bj = 1, nSy
         DO bi = 1, nSx
          DO k = 1, RNF_nSrcTile
            RNF_bufFlux(k,bi,bj,iBuf) = RNF_bufFlux(k,bi,bj,jBuf)
          ENDDO
         ENDDO
        ENDDO
        RNF_bufSum(iBuf) = RNF_bufSum(jBuf)
      ELSE
        CALL RNF_FILE_NAME( iYear,
     O       fileName,
     I       myThid )
        CALL RNF_NC_READ_FLUX( fileName, iRec,
     O       RNF_bufFlux(1,1,1,iBuf), RNF_bufSum(iBuf),
     I       myThid )
      ENDIF
      RNF_bufRec(iBuf) = iRec
      RNF_bufYr(iBuf)  = iYear

#endif /* ALLOW_RNF */

      RETURN
      END
