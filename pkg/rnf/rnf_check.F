#include "RNF_OPTIONS.h"
#ifdef ALLOW_EXF
# include "EXF_OPTIONS.h"
#endif /* ALLOW_EXF */

CBOP
C !ROUTINE: RNF_CHECK

C !INTERFACE: ==========================================================
      SUBROUTINE RNF_CHECK( myThid )

C !DESCRIPTION:
C     Check the RNF configuration and its dependencies on other
C     packages. Each refusal is fatal: all of them are printed, then
C     the run stops. The refusals of a blank RNF_file and of a build
C     without NetCDF are in RNF_READPARMS, because they guard the
C     static read in RNF_INIT_FIXED, which is called before this
C     routine.
C
C     Refused here, with useRNF true:
C     o pkg/exf not in use, or ALLOW_RUNOFF undefined: the legacy
C       route fills the exf field "runoff", which exists only with
C       that option (EXF_FIELDS.h). The generalized route keeps the
C       requirement: the record selection of both routes is pkg/exf's
C       own code (RNF_GETREC), and the dense-file exclusions below
C       read EXF_PARAM.h.
C     o runofffile not blank: a dense runoff file and the sparse
C       RNF_file are mutually exclusive.
C     o runoftempfile not blank: EXF_MAPFIELDS would apply the dense
C       runoff temperature to the sparse runoff volume.
C     o runoffconst not zero: EXF_INIT_VARIA sets "runoff" to this
C       constant and the package would overwrite it without notice.
C     o exf_outscal_sflux not 1: EXF_MAPFIELDS scales the whole
C       sflux, runoff included, when it builds EmPmR, while the exf
C       runoff heat terms and the sea-ice packages use the unscaled
C       field.
C     o SHI_update_kTopC with useShelfIce: an ice-shelf edge that
C       moves during the run can cover a target cell, which then
C       loses its volume, because SHELFICE_FORCING_SURF sets EmPmR to
C       zero under the shelf. SHI_update_kTopC is tested only when
C       useShelfIce is true: SHELFICE_READPARMS returns before it
C       sets the default otherwise.
C     o on the generalized route (RNF_applyAsRain = F, the default):
C       a build without ALLOW_ADDFLUID, selectAddFluid
C       other than 1 or 2, or a non-blank addMassFile (see the test
C       below for why each one breaks the addMass contribution of
C       RNF_ADD_MASS). The model's own temp_addMass and salt_addMass
C       are NOT refused: they are one value for the whole addMass
C       array and so cannot carry per-source properties, and
C       RNF_TENDENCY_APPLY_T/S cancel the model's term on the runoff
C       share of the array instead (package design, decision 3).
C       salt_addMass is 0, not unset, in a default run
C       (model/src/ini_parms.F:1580 copies salt_EvPrRn, whose default
C       is 0), so refusing it would refuse every default set-up;
C     o USE_OLD_EXTERNAL_FORCING defined with RNF_applyAsRain: the
C       legacy route has no hooks in EXTERNAL_FORCING_T/S (the
C       generalized route has, beside ICEFRONT's).
C     o a file with runoff_temperature in a build without
C       ALLOW_RUNOFTEMP, on the legacy route (the generalized route
C       applies the heat as a tendency): the heat of the runoff goes
C       through the exf
C       field runoftemp and the EXF_MAPFIELDS term that uses it
C       (package design, decision 3), which exist only with that
C       option, as for a dense runoftempfile (EXF_CHECK). Without it
C       the temperature would be read and silently not applied.
C       RNF_useTemp = .FALSE. reads such a file on purpose.
C     o on the legacy route, a file with runoff_temperature
C       (RNF_applyT) together with
C       temp_EvPrRn set, in a build without ALLOW_ATM_TEMP (owner
C       decision 2026-10-09, option b). The model gives runoff
C       water the temperature temp_EvPrRn whenever it is set
C       (model/src/external_forcing_surf.F:299-305), exf cancels that
C       assumption for runoff only inside "#ifdef ALLOW_ATM_TEMP"
C       (pkg/exf/exf_mapfields.F:136-197), and the exf runoff-heat
C       block that follows adds the runoff heat term either way
C       (exf_mapfields.F:199-210). Without this refusal the runoff
C       would silently enter at temp_EvPrRn + T_r - theta (branches N
C       and L) or T_r - theta + temp_EvPrRn - tRef (branch U) instead
C       of at the source temperature T_r. This is a pre-existing
C       defect of the dense runoftempfile path, which the sparse path
C       inherited when its heat moved onto the same exf
C       code; it is reported upstream and not fixed here (that is
C       option c, declined). temp_EvPrRn = UNSET_RL (the default), or
C       RNF_useTemp = .FALSE., avoids it. RNF_applyT is set by
C       RNF_INIT_FIXED, which has already run when this routine is
C       called (see above).
C     o pkg/longstep stepping the ptracers with LS_nIter other than 1
C       while the file feeds at least one of them (owner decision
C       2026-10-09, option C): the tracer term is that of the
C       current step while pkg/longstep applies the long-step average
C       of the freshwater flux (LS_fwFlux), so the two would belong to
C       different times. LS_nIter = 1, or a file without tracers, is
C       exact and is not refused. RNF_usePtracers = .FALSE. reads the
C       file without its tracers. RNF_nTrUse is set by RNF_NC_SERIES in
C       RNF_INIT_FIXED, which has run by now.
C
C     The refusals that need the file itself (its header, its table
C     entries and its targets) are in RNF_INIT_FIXED, which has
C     already run when this routine is called.
C
C     NOT refused here: the magnitude of the runoff itself. pkg/exf
C     skips its runoff upper bound of 1.E-6 m/s when useRNF is true
C     (EXF_CHECK_RANGE), and the package's replacement is
C     RNF_srcFluxMax (RNF.h), enforced in RNF_NC_READ_ONE as each
C     record is read and reported by RNF_SUMMARY. It cannot live in
C     this routine: the file's time series are first read by
C     RNF_INIT_VARIA, which runs after PACKAGES_CHECK, so no flux
C     value exists yet when this routine is called. Checking at every
C     record also covers a file whose later records are the bad ones,
C     which an init-time check could not.

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "RNF_SIZE.h"
#include "RNF.h"
#ifdef ALLOW_EXF
# include "EXF_PARAM.h"
#endif /* ALLOW_EXF */
#ifdef ALLOW_SHELFICE
# include "SHELFICE.h"
#endif /* ALLOW_SHELFICE */
#ifdef ALLOW_LONGSTEP
# include "LONGSTEP_PARAMS.h"
#endif /* ALLOW_LONGSTEP */

C !INPUT PARAMETERS: ===================================================
C  myThid               :: my thread Id number
      INTEGER myThid

C !OUTPUT PARAMETERS: ==================================================
C  none

#ifdef ALLOW_RNF

C !LOCAL VARIABLES: ====================================================
C  errCount             :: number of fatal errors detected
C  msgBuf               :: message buffer
      INTEGER errCount
      CHARACTER*(MAX_LEN_MBUF) msgBuf
#ifdef ALLOW_LONGSTEP
C  trNam                :: name of one runoff tracer, for the message
      INTEGER k, iL
      CHARACTER*(RNF_idLen) trNam
      INTEGER  ILNBLNK
      EXTERNAL ILNBLNK
#endif /* ALLOW_LONGSTEP */
CEOP

      _BEGIN_MASTER(myThid)

      errCount = 0

      WRITE(msgBuf,'(A)') 'RNF_CHECK: #define ALLOW_RNF'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )

C--   pkg/exf has to be compiled, with its runoff field, and in use
#ifdef ALLOW_EXF
      IF ( .NOT.useEXF ) THEN
        WRITE(msgBuf,'(2A)') 'RNF_CHECK: useRNF=.TRUE. needs',
     &    ' useEXF=.TRUE. (pkg/rnf fills the exf runoff field)'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
      ENDIF
# ifndef ALLOW_RUNOFF
      WRITE(msgBuf,'(2A)') 'RNF_CHECK: useRNF=.TRUE. needs',
     &  ' "#define ALLOW_RUNOFF" in EXF_OPTIONS.h'
      CALL PRINT_ERROR( msgBuf, myThid )
      errCount = errCount + 1
# endif /* ndef ALLOW_RUNOFF */
#else /* ALLOW_EXF */
      WRITE(msgBuf,'(2A)') 'RNF_CHECK: useRNF=.TRUE. needs',
     &  ' pkg/exf (useEXF), which is not compiled'
      CALL PRINT_ERROR( msgBuf, myThid )
      errCount = errCount + 1
#endif /* ALLOW_EXF */

#ifdef ALLOW_EXF
C--   Dense exf runoff settings that conflict with the sparse file.
C     They are tested only when exf is in use: EXF_READPARMS returns
C     before it sets its defaults otherwise.
      IF ( useEXF ) THEN
       IF ( runofffile .NE. ' ' ) THEN
        WRITE(msgBuf,'(2A)') 'RNF_CHECK: runofffile must be blank',
     &    ' with useRNF=.TRUE. (dense or sparse runoff, not both)'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
       ENDIF
       IF ( runoftempfile .NE. ' ' ) THEN
        WRITE(msgBuf,'(2A)') 'RNF_CHECK: runoftempfile must be blank',
     &    ' with useRNF=.TRUE. (it would apply to the sparse runoff)'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
       ENDIF
       IF ( runoffconst .NE. 0. _d 0 ) THEN
        WRITE(msgBuf,'(2A)') 'RNF_CHECK: runoffconst must be 0',
     &    ' with useRNF=.TRUE. (pkg/rnf overwrites the exf runoff)'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
       ENDIF
       IF ( exf_outscal_sflux .NE. 1. _d 0 ) THEN
        WRITE(msgBuf,'(2A)') 'RNF_CHECK: exf_outscal_sflux must be 1',
     &    ' with useRNF=.TRUE. (it scales runoff in EmPmR only)'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
       ENDIF
      ENDIF
#endif /* ALLOW_EXF */

#ifdef ALLOW_SHELFICE
C--   An ice-shelf edge that moves during the run
      IF ( useShelfIce ) THEN
       IF ( SHI_update_kTopC ) THEN
        WRITE(msgBuf,'(2A)') 'RNF_CHECK: SHI_update_kTopC must be',
     &    ' .FALSE. with useRNF=.TRUE. (a moving ice-shelf edge'
        CALL PRINT_ERROR( msgBuf, myThid )
        WRITE(msgBuf,'(2A)') 'RNF_CHECK: can cover a runoff target,',
     &    ' which then loses its volume)'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
       ENDIF
      ENDIF
#endif /* ALLOW_SHELFICE */

C--   The generalized route (RNF_applyAsRain = F, the default) puts
C     the volume in the model's addMass array (RNF_ADD_MASS), which
C     exists only with ALLOW_ADDFLUID, is used by the free surface
C     only with selectAddFluid = 1 or 2 (CALC_DIV_GHAT,
C     INTEGR_CONTINUITY; -1 is a virtual flux with no volume), and is
C     reset to zero at every step only when addMassFile is blank
C     (model/src/load_fields_driver.F:116). With a file the package
C     would add its runoff to the file's field again at every step.
      IF ( .NOT.RNF_applyAsRain ) THEN
#ifndef ALLOW_ADDFLUID
        WRITE(msgBuf,'(2A)') 'RNF_CHECK: RNF_applyAsRain=.FALSE.',
     &    ' (the generalized runoff route) needs'
        CALL PRINT_ERROR( msgBuf, myThid )
        WRITE(msgBuf,'(2A)') 'RNF_CHECK: "#define ALLOW_ADDFLUID" in',
     &    ' CPP_OPTIONS.h, or RNF_applyAsRain=.TRUE. in data.rnf'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
#endif /* ndef ALLOW_ADDFLUID */
        IF ( selectAddFluid.NE.1 .AND. selectAddFluid.NE.2 ) THEN
          WRITE(msgBuf,'(2A,I4,A)') 'RNF_CHECK: RNF_applyAsRain=',
     &      '.FALSE. needs selectAddFluid = 1 or 2 (it is',
     &      selectAddFluid, ')'
          CALL PRINT_ERROR( msgBuf, myThid )
          errCount = errCount + 1
        ENDIF
        IF ( addMassFile .NE. ' ' ) THEN
          WRITE(msgBuf,'(2A)') 'RNF_CHECK: RNF_applyAsRain=.FALSE.',
     &      ' needs addMassFile blank (addMass is then not reset at'
          CALL PRINT_ERROR( msgBuf, myThid )
          WRITE(msgBuf,'(2A)') 'RNF_CHECK: every step, and the',
     &      ' runoff would accumulate in it)'
          CALL PRINT_ERROR( msgBuf, myThid )
          errCount = errCount + 1
        ENDIF
      ENDIF

#if defined(ALLOW_EXF) && defined(ALLOW_RUNOFF)
# ifndef ALLOW_RUNOFTEMP
C--   The runoff heat needs the exf runoftemp field on the legacy
C     route; the generalized route applies it as a tendency
      IF ( RNF_hasTemp .AND. RNF_applyAsRain ) THEN
        WRITE(msgBuf,'(2A)') 'RNF_CHECK: the file has',
     &    ' runoff_temperature, which needs "#define ALLOW_RUNOFTEMP"'
        CALL PRINT_ERROR( msgBuf, myThid )
        WRITE(msgBuf,'(2A)') 'RNF_CHECK: in EXF_OPTIONS.h',
     &    ' (pkg/rnf applies it through the exf runoftemp field),'
        CALL PRINT_ERROR( msgBuf, myThid )
        WRITE(msgBuf,'(2A)') 'RNF_CHECK: or RNF_useTemp=.FALSE.',
     &    ' in data.rnf to read the file without it'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
      ENDIF
# endif /* ndef ALLOW_RUNOFTEMP */
#endif /* ALLOW_EXF and ALLOW_RUNOFF */

#if defined(ALLOW_EXF) && defined(ALLOW_RUNOFF)
# ifndef ALLOW_ATM_TEMP
C--   The runoff heat without exf's ALLOW_ATM_TEMP cancellation of its
C     own temp_EvPrRn term: a pre-existing defect of the
C     dense runoftempfile path, reported upstream and not fixed here.
C     temp_EvPrRn and UNSET_RL are the model's own (PARAMS.h,
C     EEPARAMS.h), already included above with no extra header.
C     Legacy route only: the generalized route does not use exf's
C     runoff heat, so neither the defect nor its refusal applies.
      IF ( RNF_applyT .AND. RNF_applyAsRain
     &     .AND. temp_EvPrRn.NE.UNSET_RL ) THEN
        WRITE(msgBuf,'(2A)') 'RNF_CHECK: the file has',
     &    ' runoff_temperature and temp_EvPrRn is set, which needs'
        CALL PRINT_ERROR( msgBuf, myThid )
        WRITE(msgBuf,'(2A)') 'RNF_CHECK: "#define ALLOW_ATM_TEMP" in',
     &    ' EXF_OPTIONS.h (without it exf does not cancel its own'
        CALL PRINT_ERROR( msgBuf, myThid )
        WRITE(msgBuf,'(2A)') 'RNF_CHECK: temp_EvPrRn term for runoff,',
     &    ' which then enters at the wrong temperature),'
        CALL PRINT_ERROR( msgBuf, myThid )
        WRITE(msgBuf,'(2A)') 'RNF_CHECK: or temp_EvPrRn = UNSET_RL',
     &    ' (the default) in data, or RNF_useTemp=.FALSE. in data.rnf'
        CALL PRINT_ERROR( msgBuf, myThid )
        WRITE(msgBuf,'(A)') 'RNF_CHECK: to read the file without it'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
      ENDIF
# endif /* ndef ALLOW_ATM_TEMP */
#endif /* ALLOW_EXF and ALLOW_RUNOFF */

#ifdef ALLOW_LONGSTEP
C--   pkg/longstep with a long step, and runoff tracers. LS_nIter is
C     set only when usePTRACERS is (LONGSTEP_READPARMS returns before
C     it otherwise), hence the nested test.
      IF ( usePTRACERS .AND. RNF_nTrUse.GE.1 ) THEN
       IF ( LS_nIter .NE. 1 ) THEN
        WRITE(msgBuf,'(2A,I6,A)') 'RNF_CHECK: pkg/longstep with',
     &    ' LS_nIter =', LS_nIter, ' cannot take runoff tracers:'
        CALL PRINT_ERROR( msgBuf, myThid )
        DO k = 1, RNF_nTrUse
          trNam = RNF_trNam(k)
          iL = MAX( 1, ILNBLNK(trNam) )
          WRITE(msgBuf,'(3A)') 'RNF_CHECK:   runoff_ptracer_',
     &      trNam(1:iL), ' of the file'
          CALL PRINT_ERROR( msgBuf, myThid )
        ENDDO
        WRITE(msgBuf,'(2A)') 'RNF_CHECK: use LS_nIter = 1, or',
     &    ' RNF_usePtracers=.FALSE. to read the file without them'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
       ENDIF
      ENDIF
#endif /* ALLOW_LONGSTEP */

#ifdef USE_OLD_EXTERNAL_FORCING
C--   No legacy-route hooks in the old forcing routines. The
C     generalized route has its RNF_TENDENCY_APPLY_T/S calls there
C     (model/src/external_forcing.F, beside ICEFRONT's), so it is not
C     refused; no experiment compiles that path, so they are unrun.
      IF ( RNF_applyAsRain ) THEN
       WRITE(msgBuf,'(2A)') 'RNF_CHECK: RNF_applyAsRain=.TRUE. needs',
     &  ' "#undef USE_OLD_EXTERNAL_FORCING" (no pkg/rnf hooks'
       CALL PRINT_ERROR( msgBuf, myThid )
       WRITE(msgBuf,'(2A)') 'RNF_CHECK: in S/R EXTERNAL_FORCING_T',
     &  ' and EXTERNAL_FORCING_S for the legacy route)'
       CALL PRINT_ERROR( msgBuf, myThid )
       errCount = errCount + 1
      ENDIF
#endif /* USE_OLD_EXTERNAL_FORCING */

      IF ( errCount.GE.1 ) THEN
        WRITE(msgBuf,'(A,I3,A)')
     &     'RNF_CHECK: detected', errCount,' fatal error(s)'
        CALL PRINT_ERROR( msgBuf, myThid )
        CALL ALL_PROC_DIE( 0 )
        STOP 'ABNORMAL END: S/R RNF_CHECK'
      ENDIF

      WRITE(msgBuf,'(A)') 'RNF_CHECK: configuration checks passed'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )

      _END_MASTER(myThid)

#endif /* ALLOW_RNF */

      RETURN
      END
