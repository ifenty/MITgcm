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
C     o pkg/exf not in use, or ALLOW_RUNOFF undefined: the package
C       fills the exf field "runoff", which exists only with that
C       option (EXF_FIELDS.h).
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
C     o USE_OLD_EXTERNAL_FORCING defined: the package has no hooks in
C       EXTERNAL_FORCING_T/S.
C     o a file with runoff_temperature in a build without
C       ALLOW_RUNOFTEMP: the heat of the runoff goes through the exf
C       field runoftemp and the EXF_MAPFIELDS term that uses it
C       (package design, decision 3), which exist only with that
C       option, as for a dense runoftempfile (EXF_CHECK). Without it
C       the temperature would be read and silently not applied.
C       RNF_useTemp = .FALSE. reads such a file on purpose.
C     o pkg/longstep stepping the ptracers with LS_nIter other than 1
C       while the file feeds at least one of them (RUNOFF-043, owner
C       decision 2026-10-09, option C): the tracer term is that of the
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

#if defined(ALLOW_EXF) && defined(ALLOW_RUNOFF)
# ifndef ALLOW_RUNOFTEMP
C--   The runoff heat needs the exf runoftemp field
      IF ( RNF_hasTemp ) THEN
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
C--   No pkg/rnf hooks in the old forcing routines
      WRITE(msgBuf,'(2A)') 'RNF_CHECK: useRNF=.TRUE. needs',
     &  ' "#undef USE_OLD_EXTERNAL_FORCING" (no pkg/rnf hooks'
      CALL PRINT_ERROR( msgBuf, myThid )
      WRITE(msgBuf,'(2A)') 'RNF_CHECK: in S/R EXTERNAL_FORCING_T',
     &  ' and EXTERNAL_FORCING_S)'
      CALL PRINT_ERROR( msgBuf, myThid )
      errCount = errCount + 1
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
