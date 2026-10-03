#include "RNF_OPTIONS.h"

CBOP
C !ROUTINE: RNF_READPARMS

C !INTERFACE: ==========================================================
      SUBROUTINE RNF_READPARMS( myThid )

C !DESCRIPTION:
C     Initialize RNF parameters, read in data.rnf (namelist RNF_PARM01).
C     Returns at once when useRNF is false.
C     When useRNF is true, stops the run for a blank RNF_file and for
C     a build without NetCDF. These two refusals guard the static read
C     of RNF_file in RNF_INIT_FIXED. They cannot wait for RNF_CHECK,
C     because S/R INITIALISE_FIXED calls PACKAGES_INIT_FIXED before
C     PACKAGES_CHECK.

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

#ifdef ALLOW_RNF

C !LOCAL VARIABLES: ====================================================
C  iUnit                :: unit number for I/O
C  errCount             :: number of fatal errors detected
C  msgBuf               :: message buffer
      INTEGER iUnit
      INTEGER errCount
      CHARACTER*(MAX_LEN_MBUF) msgBuf
CEOP

      NAMELIST /RNF_PARM01/
     &     RNF_file, RNF_holdRecord,
     &     RNF_startDate1, RNF_startDate2,
     &     RNF_startTime, RNF_period, RNF_repCycle,
     &     RNF_useYearlyFiles,
     &     RNF_useTemp, RNF_useSalt, RNF_usePtracers,
     &     RNF_monFreq, RNF_debugLev

      IF ( .NOT.useRNF ) THEN
C-    pkg RNF is not used
        _BEGIN_MASTER(myThid)
C-    Track pkg activation status:
         RNFisON = .FALSE.
C     print a (weak) warning if data.rnf is found
         CALL PACKAGES_UNUSED_MSG( 'useRNF', ' ', ' ' )
        _END_MASTER(myThid)
        RETURN
      ENDIF

      _BEGIN_MASTER(myThid)

      errCount = 0

C This routine has been called by the main model so we set our
C internal flag to indicate we are in business
      RNFisON = .TRUE.

C Set defaults values for parameters in RNF.h
C     An unset timing override means: use the value from the file
      RNF_file           = ' '
      RNF_holdRecord     = .FALSE.
      RNF_startDate1     = 0
      RNF_startDate2     = 0
      RNF_startTime      = UNSET_RL
      RNF_period         = UNSET_RL
      RNF_repCycle       = UNSET_RL
      RNF_useYearlyFiles = .FALSE.
      RNF_useTemp        = .TRUE.
      RNF_useSalt        = .TRUE.
      RNF_usePtracers    = .TRUE.
      RNF_monFreq        = monitorFreq
      RNF_debugLev       = debugLevel

C Open and read the data.rnf file
      WRITE(msgBuf,'(A)') ' RNF_READPARMS: opening data.rnf'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )
      CALL OPEN_COPY_DATA_FILE(
     I                   'data.rnf', 'RNF_READPARMS',
     O                   iUnit,
     I                   myThid )

      READ(UNIT=iUnit,NML=RNF_PARM01)

      WRITE(msgBuf,'(A)')
     &  ' RNF_READPARMS: finished reading data.rnf'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )

C Close the open data file
#ifdef SINGLE_DISK_IO
      CLOSE(iUnit)
#else
      CLOSE(iUnit,STATUS='DELETE')
#endif /* SINGLE_DISK_IO */

C--   Refusals that guard the static read in RNF_INIT_FIXED
      IF ( RNF_file .EQ. ' ' ) THEN
        WRITE(msgBuf,'(2A)') 'RNF_READPARMS: RNF_file is blank;',
     &    ' useRNF=.TRUE. needs a sparse runoff file (data.rnf)'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
      ENDIF

#ifndef HAVE_NETCDF
      WRITE(msgBuf,'(2A)') 'RNF_READPARMS: useRNF=.TRUE. needs',
     &  ' NetCDF to read RNF_file, but HAVE_NETCDF is undefined'
      CALL PRINT_ERROR( msgBuf, myThid )
      WRITE(msgBuf,'(2A)') 'RNF_READPARMS: in this build',
     &  ' (the genmake2 NetCDF test failed)'
      CALL PRINT_ERROR( msgBuf, myThid )
      errCount = errCount + 1
#endif /* ndef HAVE_NETCDF */

      IF ( errCount.GE.1 ) THEN
        WRITE(msgBuf,'(A,I3,A)')
     &     'RNF_READPARMS: detected', errCount,' fatal error(s)'
        CALL PRINT_ERROR( msgBuf, myThid )
        CALL ALL_PROC_DIE( 0 )
        STOP 'ABNORMAL END: S/R RNF_READPARMS'
      ENDIF

      _END_MASTER(myThid)
C Everyone else must wait for the parameters to be loaded
      _BARRIER

#endif /* ALLOW_RNF */

      RETURN
      END
