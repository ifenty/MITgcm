#include "RNF_OPTIONS.h"

CBOP
C !ROUTINE: RNF_FIELDS_LOAD

C !INTERFACE: ==========================================================
      SUBROUTINE RNF_FIELDS_LOAD( myTime, myIter, myThid )

C !DESCRIPTION:
C     Put the sparse runoff of the current time step into the dense
C     per-tile fields, and make the set the salt and tracer terms read
C     belong to the right time level.
C     Called from LOAD_FIELDS_DRIVER before EXF_GETFORCING, with the
C     same myTime that pkg/exf uses.
C
C     RNF_LOAD_AT, below, does the work for one time: select the
C     records, read what no buffer holds, combine the two records into
C     the per-tile source series and spread those over the target
C     cells. This routine only decides which time levels are needed.
C
C     **Two time levels.** The volume flux goes to the exf runoff
C     array of the current step, and the heat with it: RNF_EXF_RUNOFF
C     makes the exf runoftemp field from the fields of the current
C     step, which is the step the dense runoftemp belongs to. The salt
C     and tracer terms (RNF_FORCING_SURF) must instead belong to the
C     same step as the freshwater flux that the model uses for its own
C     salinity and tracer terms (package design, decision 3, "Time
C     level"):
C
C       o branches L and U, and branch N with staggerTimeStep: the
C         current step, so the RNF_ap* fields are a copy of the ones
C         just built;
C       o branch N without staggerTimeStep (RNF_lagFlds): PmEpR is
C         the flux of the PREVIOUS step
C         (model/src/integr_continuity.F:172-179), so the RNF_ap*
C         fields are the previous step's. They are zero at the first
C         step of a run that starts at iteration 0, where PmEpR is
C         zero too (model/src/ini_nlfs_vars.F:59), and at the first
C         step of a restart they are evaluated from the records at
C         myTime - deltaTClock, which is the step the pickup's PmEpR
C         came from. Nothing is stored in a pickup for this: the
C         fields are a function of model time.
C
C     Which level each set holds is tracked by iteration number
C     (RNF_curIter, RNF_apIter), not by model time, because
C     myTime - deltaTClock of one step need not be bitwise the myTime
C     of the step before it. RNF_INIT_VARIA calls this routine for
C     the start time and LOAD_FIELDS_DRIVER calls it again at the
C     first step with the same time and iteration; the second call
C     then finds both sets already current and rebuilds the dense
C     fields only.

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
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
C  tApply               :: time level the salt and tracer terms need
C  iApply               :: its iteration number
C  needAp               :: the RNF_ap* fields are not that level yet
C  zeroAp               :: that level is before the start of the run
      _RL     tApply
      INTEGER iApply
      LOGICAL needAp, zeroAp, reload
      CHARACTER*(MAX_LEN_MBUF) msgBuf
CEOP

      tApply = myTime
      iApply = myIter
      IF ( RNF_lagFlds ) THEN
        tApply = myTime - deltaTClock
        iApply = myIter - 1
      ENDIF

C--   Every thread decides from state that no thread has changed yet
      needAp = RNF_apIter .NE. iApply
      zeroAp = myIter.EQ.nIter0 .AND. nIter0.EQ.0
      reload = RNF_curIter .NE. iApply

      IF ( RNF_lagFlds .AND. needAp ) THEN
        IF ( zeroAp ) THEN
C--      first step of a run that starts at iteration 0: the model
C        applies no freshwater flux to its tracers at this step, so
C        the package applies no salt or tracer either. One step of
C        [(mS) - m*S_ref]*mass2rUnit would otherwise be added with
C        no volume to carry it.
          CALL RNF_ZERO_APPLY( myThid )
        ELSE
          IF ( reload ) THEN
C--        first step of a restart: evaluate the fields of the step
C          before it from the records. The series has to reach that
C          far back; taking the fields as zero was rejected, because
C          the PmEpR rebuilt from the pickup may carry runoff that an
C          earlier segment applied and whose temperature and salinity
C          this file cannot know (package design, decision 3, "Time
C          before the first record"). No case is enrolled for this
C          refusal yet: no verification experiment runs branch N
C          without staggerTimeStep (RUNOFF-017).
            IF ( RNF_recPeriod.GT.0. _d 0 .AND.
     &           RNF_recCycle.EQ.0. _d 0 .AND.
     &           tApply.LT.RNF_recStart ) THEN
              WRITE(msgBuf,'(2A,1PE20.12)') 'RNF_FIELDS_LOAD: RNF:',
     &          ' the runoff series does not reach t =', tApply
              CALL PRINT_ERROR( msgBuf, myThid )
              WRITE(msgBuf,'(2A,1PE20.12)') 'RNF_FIELDS_LOAD: its',
     &          ' first record is at t =', RNF_recStart
              CALL PRINT_ERROR( msgBuf, myThid )
              WRITE(msgBuf,'(2A)') 'RNF_FIELDS_LOAD: a restart in',
     &          ' this freshwater formulation (nonlinear free'
              CALL PRINT_ERROR( msgBuf, myThid )
              WRITE(msgBuf,'(2A)') 'RNF_FIELDS_LOAD: surface, real',
     &          ' freshwater flux, no staggerTimeStep) needs the'
              CALL PRINT_ERROR( msgBuf, myThid )
              WRITE(msgBuf,'(2A)') 'RNF_FIELDS_LOAD: fields one time',
     &          ' step before the restart: give the file a record'
              CALL PRINT_ERROR( msgBuf, myThid )
              WRITE(msgBuf,'(2A)') 'RNF_FIELDS_LOAD: at least one',
     &          ' step earlier, or a repeat cycle'
              CALL PRINT_ERROR( msgBuf, myThid )
              CALL ALL_PROC_DIE( 0 )
              STOP 'ABNORMAL END: S/R RNF_FIELDS_LOAD'
            ENDIF
            CALL RNF_LOAD_AT( tApply, iApply, .FALSE., myThid )
          ENDIF
C--        the dense fields now hold the level the terms need
          CALL RNF_COPY_APPLY( myThid )
        ENDIF
        _BEGIN_MASTER(myThid)
        RNF_apIter = iApply
        _END_MASTER(myThid)
        _BARRIER
      ENDIF

C--   The fields of this step: the volume flux pkg/exf applies, and,
C     without the lag, the properties the terms apply as well
      CALL RNF_LOAD_AT( myTime, myIter, .TRUE., myThid )

      IF ( .NOT.RNF_lagFlds ) THEN
        CALL RNF_COPY_APPLY( myThid )
        _BEGIN_MASTER(myThid)
        RNF_apIter = myIter
        _END_MASTER(myThid)
        _BARRIER
      ENDIF

#endif /* ALLOW_RNF */

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: RNF_LOAD_AT

C !INTERFACE: ==========================================================
      SUBROUTINE RNF_LOAD_AT( tLoad, iLoad, doTrace, myThid )

C !DESCRIPTION:
C     Build the dense per-tile fields of one time: select the time
C     records, read the ones that are not already in a buffer, combine
C     them into the per-tile source series and spread those over the
C     target cells.
C
C     **Records.** RNF_GETREC returns the two records that bracket
C     tLoad, the weight fac of the earlier one and the file year each
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
C     **Interpolation.** Every series is interpolated on its own, with
C     the one set of weights (package design, decision 7):
C
C       flux_s = fac*flux_s(count0) + (1-fac)*flux_s(count1) [m^3/s]
C       T_s    = fac*T_s(count0)    + (1-fac)*T_s(count1)    [degC]
C
C     and the same for the salinity and each tracer. That is the form
C     pkg/exf interpolates its dense fields with
C     (pkg/exf/exf_set_fld.F), applied to the source vectors instead
C     of to a 2-D field. With RNF_holdRecord the weight is exactly 1
C     or exactly 0 and the selected record is copied, with no
C     arithmetic, so the applied values are values of the file and not
C     combinations of two.
C
C     A source whose temperature is missing in any record used is left
C     out of the temperature sums altogether (RNF_srcTvld = 0): it
C     enters at the reference temperature, as a source without a
C     temperature does (package design, decision 3, "Missing
C     temperature"). Every other series refuses a missing value when
C     it is read.
C
C     **Dense fields.** For every target cell c of the tile, with
C     w_{s,c} = flux_s*frac_{s,c}/rA(c) the volume flux per unit area
C     that source s sends to it,
C
C       RNF_vflx(c)  = sum_s w_{s,c}                            [m/s]
C       RNF_mflx(c)  = rhoConstFresh*RNF_vflx(c)          [kg/m^2/s]
C       RNF_mflxT(c) = rhoConstFresh*sum_{s: T present} w_{s,c}
C       RNF_mXT(c)   = rhoConstFresh*sum_{s: T present} w_{s,c}*T_s
C       RNF_mXS(c)   = rhoConstFresh*sum_s w_{s,c}*S_s
C       RNF_mXTr(c,n)= rhoConstFresh*sum_s w_{s,c}*C_{s,n}
C
C     i.e. volumes ADD while the properties are FLUX-WEIGHTED, which
C     is what conserves heat, salt and tracer content when several
C     sources feed one cell (model contract). No division by the
C     total flux is made here, so a cell with zero flux needs no
C     special case: the salt and tracer terms use the sums directly.
C     The one division, (mT)/m for the exf runoftemp field, is made in
C     RNF_EXF_RUNOFF, where m = 0 is a case of its own.
C     The property sums are accumulated in volume-flux units and
C     scaled to mass flux in one pass at the end, as RNF_mflx is.
C
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
C     would void both. The RNF_vflx statement is kept exactly as it
C     was for the same reason; the property sums use a separate
C     expression of the same value. The volume flux is kept in m/s and
C     not derived from the mass flux: multiplying and dividing by
C     rhoConstFresh does not give back the same number.

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "RNF_SIZE.h"
#include "RNF.h"

C !INPUT PARAMETERS: ===================================================
C  tLoad                :: time to build the fields for
C  iLoad                :: iteration number that time belongs to
C  doTrace              :: print the record-selection trace line.
C                          False for the extra load of the lagged
C                          level at a restart, so that the trace
C                          stays one line per step, which is the
C                          interface the timing checks read
C  myThid               :: my thread Id number
      _RL     tLoad
      INTEGER iLoad
      LOGICAL doTrace
      INTEGER myThid

C !OUTPUT PARAMETERS: ==================================================
C  none

#ifdef ALLOW_RNF

C !LOCAL VARIABLES: ====================================================
C  bi, bj               :: tile indices
C  i, j, k              :: loop indices
C  kSrc                 :: local source index of a target entry
C  iTr                  :: runoff tracer index
C  fac                  :: weight of record count0
C  count0, count1       :: the two records that bracket tLoad
C  year0, year1         :: the file year each one comes from
C  need0, need1         :: whether that record has a non-zero weight
C  wVol                 :: volume flux per unit area of one target
C  tVal                 :: interpolated source temperature
C  useProp              :: any property series to spread at all
      INTEGER bi, bj, i, j, k, kSrc, iTr
      _RL     fac
      LOGICAL first, changed
      INTEGER count0, count1, year0, year1
      LOGICAL need0, need1, useProp, tVld
      _RL     wVol, tVal
      CHARACTER*(MAX_LEN_MBUF) msgBuf
CEOP

      useProp = RNF_applyT .OR. RNF_applyS .OR. RNF_nTrUse.GE.1

C--   1. Which records, and with which weights
      CALL RNF_GETREC(
     I     tLoad, iLoad,
     O     fac, first, changed,
     O     count0, count1, year0, year1,
     I     myThid )

      need0 = fac .GT. 0. _d 0
      need1 = fac .LT. 1. _d 0

C--   2. The records themselves, and the source series of this time.
C     The master thread does the NetCDF reads (model contract) and
C     fills the per-tile source vectors of every tile of the process.
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

        IF ( RNF_applyT ) THEN
         DO k = 1, RNF_nSrc(bi,bj)
C--       the temperature counts only where every record used has one
          tVld = .TRUE.
          IF ( need0 ) tVld = RNF_bufTvld(k,bi,bj,1) .GT. 0. _d 0
          IF ( need1 .AND. tVld )
     &         tVld = RNF_bufTvld(k,bi,bj,2) .GT. 0. _d 0
          IF ( tVld ) THEN
            RNF_srcTvld(k,bi,bj) = 1. _d 0
            IF     ( .NOT.need1 ) THEN
              tVal = RNF_bufTemp(k,bi,bj,1)
            ELSEIF ( .NOT.need0 ) THEN
              tVal = RNF_bufTemp(k,bi,bj,2)
            ELSE
              tVal = fac             *RNF_bufTemp(k,bi,bj,1)
     &           + ( 1. _d 0 - fac )*RNF_bufTemp(k,bi,bj,2)
            ENDIF
          ELSE
            RNF_srcTvld(k,bi,bj) = 0. _d 0
            tVal = 0. _d 0
          ENDIF
          RNF_srcTemp(k,bi,bj) = tVal
         ENDDO
        ENDIF

        IF ( RNF_hasSalt ) THEN
         DO k = 1, RNF_nSrc(bi,bj)
          IF     ( .NOT.need1 ) THEN
            RNF_srcSalt(k,bi,bj) = RNF_bufSalt(k,bi,bj,1)
          ELSEIF ( .NOT.need0 ) THEN
            RNF_srcSalt(k,bi,bj) = RNF_bufSalt(k,bi,bj,2)
          ELSE
            RNF_srcSalt(k,bi,bj) =
     &             fac             *RNF_bufSalt(k,bi,bj,1)
     &         + ( 1. _d 0 - fac )*RNF_bufSalt(k,bi,bj,2)
          ENDIF
         ENDDO
        ENDIF

        DO iTr = 1, RNF_nTrUse
         DO k = 1, RNF_nSrc(bi,bj)
          IF     ( .NOT.need1 ) THEN
            RNF_srcTrc(k,bi,bj,iTr) = RNF_bufTrc(k,bi,bj,1,iTr)
          ELSEIF ( .NOT.need0 ) THEN
            RNF_srcTrc(k,bi,bj,iTr) = RNF_bufTrc(k,bi,bj,2,iTr)
          ELSE
            RNF_srcTrc(k,bi,bj,iTr) =
     &             fac             *RNF_bufTrc(k,bi,bj,1,iTr)
     &         + ( 1. _d 0 - fac )*RNF_bufTrc(k,bi,bj,2,iTr)
          ENDIF
         ENDDO
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
      IF ( doTrace .AND. ( RNF_debugLev.GE.debLevC .OR.
     &     ( RNF_debugLev.GE.debLevB .AND.
     &       ( first .OR. changed ) ) ) ) THEN
        WRITE(msgBuf,'(A,I10,A,I6,A,I6,A,I6,A,I6,A,1PE25.16,A,E25.16)')
     &    ' RNF_FIELDS_LOAD: it=', iLoad,
     &    ', rec0=', count0, ', yr0=', year0,
     &    ', rec1=', count1, ', yr1=', year1,
     &    ', t=', tLoad, ', fac=', fac
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                      SQUEEZE_RIGHT, myThid )
      ENDIF

      _END_MASTER(myThid)
C--   Every thread waits for the source series of this step
      _BARRIER

C--   3. The dense fields of the tiles of this thread
      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)

        DO j = 1-OLy, sNy+OLy
         DO i = 1-OLx, sNx+OLx
          RNF_vflx(i,j,bi,bj) = 0. _d 0
          RNF_mflx(i,j,bi,bj) = 0. _d 0
          RNF_mflxT(i,j,bi,bj) = 0. _d 0
          RNF_mXT(i,j,bi,bj) = 0. _d 0
          RNF_mXS(i,j,bi,bj) = 0. _d 0
         ENDDO
        ENDDO
        DO iTr = 1, RNF_nTr
         DO j = 1-OLy, sNy+OLy
          DO i = 1-OLx, sNx+OLx
           RNF_mXTr(i,j,bi,bj,iTr) = 0. _d 0
          ENDDO
         ENDDO
        ENDDO

        DO k = 1, RNF_nTgt(bi,bj)
          i = RNF_tgtI(k,bi,bj)
          j = RNF_tgtJ(k,bi,bj)
          kSrc = RNF_tgtSrc(k,bi,bj)
          RNF_vflx(i,j,bi,bj) = RNF_vflx(i,j,bi,bj)
     &      + RNF_srcFlux(kSrc,bi,bj)*RNF_tgtFrac(k,bi,bj)
     &        /rA(i,j,bi,bj)
          IF ( useProp ) THEN
            wVol = RNF_srcFlux(kSrc,bi,bj)*RNF_tgtFrac(k,bi,bj)
     &             /rA(i,j,bi,bj)
            IF ( RNF_applyT ) THEN
              RNF_mflxT(i,j,bi,bj) = RNF_mflxT(i,j,bi,bj)
     &          + wVol*RNF_srcTvld(kSrc,bi,bj)
              RNF_mXT(i,j,bi,bj) = RNF_mXT(i,j,bi,bj)
     &          + wVol*RNF_srcTvld(kSrc,bi,bj)*RNF_srcTemp(kSrc,bi,bj)
            ENDIF
            IF ( RNF_hasSalt ) THEN
              RNF_mXS(i,j,bi,bj) = RNF_mXS(i,j,bi,bj)
     &          + wVol*RNF_srcSalt(kSrc,bi,bj)
            ENDIF
            DO iTr = 1, RNF_nTrUse
              RNF_mXTr(i,j,bi,bj,iTr) = RNF_mXTr(i,j,bi,bj,iTr)
     &          + wVol*RNF_srcTrc(kSrc,bi,bj,iTr)
            ENDDO
          ENDIF
        ENDDO

        DO j = 1, sNy
         DO i = 1, sNx
          RNF_mflx(i,j,bi,bj) = rhoConstFresh*RNF_vflx(i,j,bi,bj)
          RNF_mflxT(i,j,bi,bj) = rhoConstFresh*RNF_mflxT(i,j,bi,bj)
          RNF_mXT(i,j,bi,bj) = rhoConstFresh*RNF_mXT(i,j,bi,bj)
          RNF_mXS(i,j,bi,bj) = rhoConstFresh*RNF_mXS(i,j,bi,bj)
         ENDDO
        ENDDO
        DO iTr = 1, RNF_nTrUse
         DO j = 1, sNy
          DO i = 1, sNx
           RNF_mXTr(i,j,bi,bj,iTr) =
     &          rhoConstFresh*RNF_mXTr(i,j,bi,bj,iTr)
          ENDDO
         ENDDO
        ENDDO

       ENDDO
      ENDDO

      _BEGIN_MASTER(myThid)
      RNF_curIter = iLoad
      _END_MASTER(myThid)
      _BARRIER

#endif /* ALLOW_RNF */

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: RNF_COPY_APPLY

C !INTERFACE: ==========================================================
      SUBROUTINE RNF_COPY_APPLY( myThid )

C !DESCRIPTION:
C     Copy the dense mass-flux and property fields that RNF_LOAD_AT
C     just built into the set the salt and tracer terms read (RNF_ap*
C     of RNF.h), and exchange it. Called by every thread for its own
C     tiles.
C
C     The two sets exist so that RNF_FORCING_SURF need not know which
C     time level it is applying: the time-level rule of decision 3 is
C     decided once, in RNF_FIELDS_LOAD.
C
C     **The exchange.** RNF_LOAD_AT fills the tile interior only, and
C     the surface forcing fields the terms go into are used over the
C     halo as well: EXTERNAL_FORCING_SURF fills them there and KPP_CALC
C     reads them there (pkg/kpp/kpp_calc.F, iMin = 2-OLx). The fields
C     those routines take from other packages are exchanged before
C     they are used (EmPmR and Qnet in EXF_MAPFIELDS, saltFlux in
C     SEAICE_MODEL), so these are too. Nothing is exchanged when
C     neither term can be non-zero, which leaves a run whose file
C     carries no salinity or tracer exactly as it was.

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "RNF_SIZE.h"
#include "RNF.h"

C !INPUT PARAMETERS: ===================================================
C  myThid               :: my thread Id number
      INTEGER myThid

#ifdef ALLOW_RNF

C !LOCAL VARIABLES: ====================================================
      INTEGER bi, bj, i, j, iTr
CEOP

      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)
        DO j = 1-OLy, sNy+OLy
         DO i = 1-OLx, sNx+OLx
          RNF_apMflx(i,j,bi,bj)  = RNF_mflx(i,j,bi,bj)
          RNF_apXS(i,j,bi,bj)    = RNF_mXS(i,j,bi,bj)
         ENDDO
        ENDDO
        DO iTr = 1, RNF_nTr
         DO j = 1-OLy, sNy+OLy
          DO i = 1-OLx, sNx+OLx
           RNF_apXTr(i,j,bi,bj,iTr) = RNF_mXTr(i,j,bi,bj,iTr)
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      _BARRIER

      IF ( RNF_applyS .OR. RNF_nTrUse.GE.1 ) THEN
        _EXCH_XY_RL( RNF_apMflx, myThid )
      ENDIF
      IF ( RNF_applyS ) THEN
        _EXCH_XY_RL( RNF_apXS, myThid )
      ENDIF
      DO iTr = 1, RNF_nTrUse
        CALL EXCH_XY_RL( RNF_apXTr(1-OLx,1-OLy,1,1,iTr), myThid )
      ENDDO

#endif /* ALLOW_RNF */

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: RNF_ZERO_APPLY

C !INTERFACE: ==========================================================
      SUBROUTINE RNF_ZERO_APPLY( myThid )

C !DESCRIPTION:
C     Set the fields the salt and tracer terms read to zero, so that
C     no salt or tracer is applied. Used at the first step of a run
C     that starts at iteration 0 when the terms lag the step
C     (RNF_lagFlds): the model's own PmEpR is zero there
C     (model/src/ini_nlfs_vars.F:59), so there is no volume to carry
C     the properties. The heat is not lagged: it follows the exf
C     runoff field of the current step, as the dense runoftemp does.

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "RNF_SIZE.h"
#include "RNF.h"

C !INPUT PARAMETERS: ===================================================
C  myThid               :: my thread Id number
      INTEGER myThid

#ifdef ALLOW_RNF

C !LOCAL VARIABLES: ====================================================
      INTEGER bi, bj, i, j, iTr
CEOP

      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)
        DO j = 1-OLy, sNy+OLy
         DO i = 1-OLx, sNx+OLx
          RNF_apMflx(i,j,bi,bj)  = 0. _d 0
          RNF_apXS(i,j,bi,bj)    = 0. _d 0
         ENDDO
        ENDDO
        DO iTr = 1, RNF_nTr
         DO j = 1-OLy, sNy+OLy
          DO i = 1-OLx, sNx+OLx
           RNF_apXTr(i,j,bi,bj,iTr) = 0. _d 0
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      _BARRIER

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
C     A buffer holds every series of its record, so the copy and the
C     read carry the temperature, salinity and tracers with the flux.
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
      INTEGER jBuf, bi, bj, k, iTr
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
            RNF_bufTemp(k,bi,bj,iBuf) = RNF_bufTemp(k,bi,bj,jBuf)
            RNF_bufTvld(k,bi,bj,iBuf) = RNF_bufTvld(k,bi,bj,jBuf)
            RNF_bufSalt(k,bi,bj,iBuf) = RNF_bufSalt(k,bi,bj,jBuf)
          ENDDO
          DO iTr = 1, RNF_nTr
           DO k = 1, RNF_nSrcTile
            RNF_bufTrc(k,bi,bj,iBuf,iTr) =
     &           RNF_bufTrc(k,bi,bj,jBuf,iTr)
           ENDDO
          ENDDO
         ENDDO
        ENDDO
        RNF_bufSum(iBuf) = RNF_bufSum(jBuf)
      ELSE
        CALL RNF_FILE_NAME( iYear,
     O       fileName,
     I       myThid )
        CALL RNF_NC_READ_FLUX( fileName, iRec, iBuf,
     I       myThid )
      ENDIF
      RNF_bufRec(iBuf) = iRec
      RNF_bufYr(iBuf)  = iYear

#endif /* ALLOW_RNF */

      RETURN
      END
