#include "RNF_OPTIONS.h"

C--  File rnf_nc_utils.F: NetCDF helpers of the sparse runoff package
C--   Contents
C--   o RNF_NC_ERROR     :: stop the run on a NetCDF error
C--   o RNF_NC_ATT_TEXT  :: read a text attribute
C--   o RNF_NC_ATT_REAL  :: read an attribute that holds one number
C--   o RNF_NC_SOURCE_ID :: read the id of one source
C--   o RNF_NC_TGT_ERROR :: print one error line about a target entry
C--   o RNF_NC_SERIES    :: which time series the file carries, and
C--                         the ptracer each runoff tracer feeds
C--   o RNF_NC_READ_FLUX :: read one record of every time series from
C--                         one file of the set
C--   o RNF_NC_READ_ONE  :: read one record of one time series
C--
C--  Every NF_* call and the netcdf.inc include are inside
C--  #ifdef HAVE_NETCDF (package design, decision 6). Without NetCDF
C--  the package still compiles; RNF_READPARMS stops the run when
C--  useRNF is true, before any of these routines can be reached.

CBOP
C !ROUTINE: RNF_NC_ERROR

C !INTERFACE: ==========================================================
      SUBROUTINE RNF_NC_ERROR(
     I                   caller, operation, errNC, myThid )

C !DESCRIPTION:
C     Stop the run when a NetCDF call failed. Does nothing when
C     errNC reports success. The message names the routine that
C     called, the operation that failed, the NetCDF error text and
C     RNF_file. Every process reads the same file, so every process
C     reaches the same error by itself.

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "RNF_SIZE.h"
#include "RNF.h"
#ifdef HAVE_NETCDF
# include "netcdf.inc"
#endif

C !INPUT PARAMETERS: ===================================================
C  caller               :: name of the calling routine
C  operation            :: what was attempted, e.g. "opening the file"
C  errNC                :: status returned by the NetCDF call
C  myThid               :: my thread Id number
      CHARACTER*(*) caller
      CHARACTER*(*) operation
      INTEGER errNC
      INTEGER myThid

C !OUTPUT PARAMETERS: ==================================================
C  none

#ifdef ALLOW_RNF
#ifdef HAVE_NETCDF

C !FUNCTIONS: ==========================================================
      INTEGER  ILNBLNK
      EXTERNAL ILNBLNK

C !LOCAL VARIABLES: ====================================================
C  iL, jL, kL           :: length of the trimmed strings
C  msgBuf               :: message buffer
      INTEGER iL, jL, kL
      CHARACTER*(MAX_LEN_MBUF) msgBuf
CEOP

      IF ( errNC .EQ. NF_NOERR ) RETURN

      iL = MAX( 1, ILNBLNK(caller) )
      jL = MAX( 1, ILNBLNK(operation) )
      WRITE(msgBuf,'(4A)') caller(1:iL), ': RNF: ', operation(1:jL),
     &     ' failed'
      CALL PRINT_ERROR( msgBuf, myThid )
      WRITE(msgBuf,'(2A)') 'RNF_NC_ERROR: NetCDF says: ',
     &     NF_STRERROR(errNC)
      CALL PRINT_ERROR( msgBuf, myThid )
      kL = MAX( 1, ILNBLNK(RNF_file) )
      WRITE(msgBuf,'(2A)') 'RNF_NC_ERROR: RNF_file = ',
     &     RNF_file(1:kL)
      CALL PRINT_ERROR( msgBuf, myThid )
      CALL ALL_PROC_DIE( 0 )
      STOP 'ABNORMAL END: S/R RNF_NC_ERROR'

#endif /* HAVE_NETCDF */
#endif /* ALLOW_RNF */

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: RNF_NC_ATT_TEXT

C !INTERFACE: ==========================================================
      SUBROUTINE RNF_NC_ATT_TEXT(
     I                   fid, varId, attName, caller,
     O                   attVal, attLen, found,
     I                   myThid )

C !DESCRIPTION:
C     Read a text attribute of an open sparse runoff file. Only "no
C     such attribute" returns found = .FALSE.; an attribute that is
C     not of type char, or longer than the output string, stops the
C     run. The caller decides whether a missing attribute is fatal.

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "RNF_SIZE.h"
#include "RNF.h"
#ifdef HAVE_NETCDF
# include "netcdf.inc"
#endif

C !INPUT PARAMETERS: ===================================================
C  fid                  :: NetCDF id of the open file
C  varId                :: variable id, or NF_GLOBAL
C  attName              :: name of the attribute
C  caller               :: name of the calling routine
C  myThid               :: my thread Id number
      INTEGER fid
      INTEGER varId
      CHARACTER*(*) attName
      CHARACTER*(*) caller
      INTEGER myThid

C !OUTPUT PARAMETERS: ==================================================
C  attVal               :: the text, blank-padded
C  attLen               :: its length
C  found                :: .FALSE. if the file has no such attribute
      CHARACTER*(*) attVal
      INTEGER attLen
      LOGICAL found

#ifdef ALLOW_RNF

C !FUNCTIONS: ==========================================================
      INTEGER  ILNBLNK
      EXTERNAL ILNBLNK

C !LOCAL VARIABLES: ====================================================
      INTEGER i, iL, jL
      CHARACTER*(MAX_LEN_MBUF) msgBuf
#ifdef HAVE_NETCDF
      INTEGER errNC, attType, nVals
#endif
CEOP

      DO i = 1, LEN(attVal)
        attVal(i:i) = ' '
      ENDDO
      attLen = 0
      found  = .FALSE.
      iL = MAX( 1, ILNBLNK(caller) )
      jL = MAX( 1, ILNBLNK(attName) )

#ifdef HAVE_NETCDF
      errNC = NF_INQ_ATT( fid, varId, attName, attType, nVals )
      IF ( errNC .EQ. NF_ENOTATT ) RETURN
      WRITE(msgBuf,'(2A)') 'reading the attribute ', attName(1:jL)
      CALL RNF_NC_ERROR( caller, msgBuf, errNC, myThid )

      IF ( attType.NE.NF_CHAR .OR. nVals.GT.LEN(attVal) ) THEN
        WRITE(msgBuf,'(4A,I6,A)') caller(1:iL),
     &    ': RNF: the attribute ', attName(1:jL),
     &    ' must be char text of at most', LEN(attVal),
     &    ' characters'
        CALL PRINT_ERROR( msgBuf, myThid )
        CALL ALL_PROC_DIE( 0 )
        STOP 'ABNORMAL END: S/R RNF_NC_ATT_TEXT'
      ENDIF

      errNC = NF_GET_ATT_TEXT( fid, varId, attName, attVal )
      WRITE(msgBuf,'(2A)') 'reading the attribute ', attName(1:jL)
      CALL RNF_NC_ERROR( caller, msgBuf, errNC, myThid )
C-    a char attribute may be NUL-terminated rather than blank-padded
      DO i = 1, LEN(attVal)
        IF ( attVal(i:i) .EQ. CHAR(0) ) attVal(i:i) = ' '
      ENDDO
      DO i = nVals+1, LEN(attVal)
        attVal(i:i) = ' '
      ENDDO
      attLen = ILNBLNK( attVal )
      found  = .TRUE.
#endif /* HAVE_NETCDF */

#endif /* ALLOW_RNF */

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: RNF_NC_ATT_REAL

C !INTERFACE: ==========================================================
      SUBROUTINE RNF_NC_ATT_REAL(
     I                   fid, varId, varName, attName, caller,
     O                   attVal, found,
     I                   myThid )

C !DESCRIPTION:
C     Read an attribute that must hold exactly one number, such as
C     mitgcm_grid_nx or the _FillValue of runoff_flux. Only "no such
C     attribute" returns found = .FALSE.; an attribute that exists
C     but is text, or is a list, or cannot be read as a number, stops
C     the run. Treating such a marker as absent would apply the
C     marker value as runoff.

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "RNF_SIZE.h"
#include "RNF.h"
#ifdef HAVE_NETCDF
# include "netcdf.inc"
#endif

C !INPUT PARAMETERS: ===================================================
C  fid                  :: NetCDF id of the open file
C  varId                :: variable id, or NF_GLOBAL
C  varName              :: name of the variable, for the message
C                          ("the file" for a global attribute)
C  attName              :: name of the attribute
C  caller               :: name of the calling routine
C  myThid               :: my thread Id number
      INTEGER fid
      INTEGER varId
      CHARACTER*(*) varName
      CHARACTER*(*) attName
      CHARACTER*(*) caller
      INTEGER myThid

C !OUTPUT PARAMETERS: ==================================================
C  attVal               :: the value
C  found                :: .FALSE. if the file has no such attribute
      _RL attVal
      LOGICAL found

#ifdef ALLOW_RNF

C !FUNCTIONS: ==========================================================
      INTEGER  ILNBLNK
      EXTERNAL ILNBLNK

C !LOCAL VARIABLES: ====================================================
      INTEGER iL, jL, kL
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      LOGICAL oneNumber
#ifdef HAVE_NETCDF
      INTEGER errNC, attType, nVals
      Real*8  tmpVal(1)
#endif
CEOP

      attVal = 0. _d 0
      found  = .FALSE.
      iL = MAX( 1, ILNBLNK(caller) )
      jL = MAX( 1, ILNBLNK(attName) )
      kL = MAX( 1, ILNBLNK(varName) )

#ifdef HAVE_NETCDF
      errNC = NF_INQ_ATT( fid, varId, attName, attType, nVals )
      IF ( errNC .EQ. NF_ENOTATT ) RETURN
      WRITE(msgBuf,'(2A)') 'reading the attribute ', attName(1:jL)
      CALL RNF_NC_ERROR( caller, msgBuf, errNC, myThid )

      oneNumber = nVals.EQ.1 .AND. attType.NE.NF_CHAR
      IF ( oneNumber ) THEN
        errNC = NF_GET_ATT_DOUBLE( fid, varId, attName, tmpVal )
        oneNumber = errNC .EQ. NF_NOERR
        attVal = tmpVal(1)
      ENDIF
      IF ( .NOT.oneNumber ) THEN
        WRITE(msgBuf,'(6A)') caller(1:iL),
     &    ': RNF: the attribute ', attName(1:jL), ' of ',
     &    varName(1:kL), ' must be one number'
        CALL PRINT_ERROR( msgBuf, myThid )
        WRITE(msgBuf,'(A,I6,A,I6,A)') 'RNF_NC_ATT_REAL: it has type',
     &    attType, ' and', nVals, ' value(s)'
        CALL PRINT_ERROR( msgBuf, myThid )
        CALL ALL_PROC_DIE( 0 )
        STOP 'ABNORMAL END: S/R RNF_NC_ATT_REAL'
      ENDIF
      found = .TRUE.
#endif /* HAVE_NETCDF */

#endif /* ALLOW_RNF */

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: RNF_NC_SOURCE_ID

C !INTERFACE: ==========================================================
      SUBROUTINE RNF_NC_SOURCE_ID(
     I                   fid, iSrc, caller,
     O                   srcId, idLen,
     I                   myThid )

C !DESCRIPTION:
C     Read the id of source iSrc (0-based, as the file counts) from
C     the open sparse runoff file. Ids are read only when a message
C     needs one, so no array has the length of the source dimension.
C     NUL padding of the char array is turned into blanks.

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "RNF_SIZE.h"
#include "RNF.h"
#ifdef HAVE_NETCDF
# include "netcdf.inc"
#endif

C !INPUT PARAMETERS: ===================================================
C  fid                  :: NetCDF id of the open file
C  iSrc                 :: 0-based source index in the file
C  caller               :: name of the calling routine
C  myThid               :: my thread Id number
      INTEGER fid
      INTEGER iSrc
      CHARACTER*(*) caller
      INTEGER myThid

C !OUTPUT PARAMETERS: ==================================================
C  srcId                :: the id, blank-padded
C  idLen                :: its length (1 if the id is empty)
      CHARACTER*(*) srcId
      INTEGER idLen

#ifdef ALLOW_RNF

C !FUNCTIONS: ==========================================================
      INTEGER  ILNBLNK
      EXTERNAL ILNBLNK

C !LOCAL VARIABLES: ====================================================
      INTEGER i
#ifdef HAVE_NETCDF
      INTEGER errNC, varId, nDims, nChar
      INTEGER dimIds(NF_MAX_VAR_DIMS)
      INTEGER iStart(2), iCount(2)
#endif
CEOP

      DO i = 1, LEN(srcId)
        srcId(i:i) = ' '
      ENDDO
      idLen = 1

#ifdef HAVE_NETCDF
      errNC = NF_INQ_VARID( fid, 'source_id', varId )
      CALL RNF_NC_ERROR( caller,
     &     'looking for the required variable source_id',
     &     errNC, myThid )
      errNC = NF_INQ_VARNDIMS( fid, varId, nDims )
      CALL RNF_NC_ERROR( caller, 'reading the shape of source_id',
     &     errNC, myThid )
      nChar = 1
      IF ( nDims .GE. 2 ) THEN
        errNC = NF_INQ_VARDIMID( fid, varId, dimIds )
        CALL RNF_NC_ERROR( caller, 'reading the shape of source_id',
     &       errNC, myThid )
C-      in Fortran the character dimension varies fastest
        errNC = NF_INQ_DIMLEN( fid, dimIds(1), nChar )
        CALL RNF_NC_ERROR( caller, 'reading the length of the id',
     &       errNC, myThid )
      ENDIF
      nChar = MIN( nChar, LEN(srcId) )
      iStart(1) = 1
      iStart(2) = iSrc + 1
      iCount(1) = nChar
      iCount(2) = 1
      errNC = NF_GET_VARA_TEXT( fid, varId, iStart, iCount, srcId )
      CALL RNF_NC_ERROR( caller, 'reading source_id', errNC, myThid )
      DO i = 1, LEN(srcId)
        IF ( srcId(i:i) .EQ. CHAR(0) ) srcId(i:i) = ' '
      ENDDO
      DO i = nChar+1, LEN(srcId)
        srcId(i:i) = ' '
      ENDDO
      idLen = MAX( 1, ILNBLNK(srcId) )
#endif /* HAVE_NETCDF */

#endif /* ALLOW_RNF */

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: RNF_NC_TGT_ERROR

C !INTERFACE: ==========================================================
      SUBROUTINE RNF_NC_TGT_ERROR(
     I                   caller, text, srcId, idLen, cell,
     I                   i, j, bi, bj, myThid )

C !DESCRIPTION:
C     Print one error line about a target entry that the model
C     refuses, naming the source id, the cell index of the file and
C     the local indices where the entry was placed. The run is
C     stopped by the caller, after the counts of all processes have
C     been summed, so that every process stops.

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"

C !INPUT PARAMETERS: ===================================================
C  caller               :: name of the calling routine
C  text                 :: what is wrong with the entry
C  srcId, idLen         :: id of the source and its length
C  cell                 :: 0-based global cell index from the file
C  i, j, bi, bj         :: where the entry was placed
C  myThid               :: my thread Id number
      CHARACTER*(*) caller
      CHARACTER*(*) text
      CHARACTER*(*) srcId
      INTEGER idLen
      INTEGER cell
      INTEGER i, j, bi, bj
      INTEGER myThid

C !OUTPUT PARAMETERS: ==================================================
C  none

#ifdef ALLOW_RNF

C !FUNCTIONS: ==========================================================
      INTEGER  ILNBLNK
      EXTERNAL ILNBLNK

C !LOCAL VARIABLES: ====================================================
      INTEGER iL, jL
      CHARACTER*(MAX_LEN_MBUF) msgBuf
CEOP

      iL = MAX( 1, ILNBLNK(caller) )
      jL = MAX( 1, ILNBLNK(text) )
      WRITE(msgBuf,'(6A,I12,A,4I6)') caller(1:iL), ': RNF: ',
     &     text(1:jL), ': source ', srcId(1:idLen),
     &     ', target_cell', cell, ', i,j,bi,bj =', i, j, bi, bj
      CALL PRINT_ERROR( msgBuf, myThid )

#endif /* ALLOW_RNF */

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: RNF_NC_SERIES

C !INTERFACE: ==========================================================
      SUBROUTINE RNF_NC_SERIES(
     I                   fid,
     U                   errCount,
     I                   myThid )

C !DESCRIPTION:
C     Find out which of the optional time series the open file fid
C     carries, and which ptracer each runoff tracer feeds. Called once
C     by the master thread of RNF_INIT_FIXED, which has the file open.
C
C     Sets RNF_hasTemp, RNF_hasSalt, RNF_nTrUse, RNF_trNam and
C     RNF_trPtr of RNF.h. A series that the file has and the user
C     switched off (RNF_useTemp, RNF_useSalt, RNF_usePtracers of
C     data.rnf) is reported and then treated as absent.
C
C     **Name matching** (package design, decision 4). Each variable
C     runoff_ptracer_<NAME> is compared with PTRACERS_names(iTr) for
C     iTr = 1 ... PTRACERS_numInUse, exactly and with trailing blanks
C     removed. Each of these raises errCount, so the run stops in
C     RNF_INIT_FIXED after every one of them has been reported:
C
C       o no PTRACERS_names entry matches <NAME>;
C       o more than one does;
C       o the file has runoff tracers while pkg/ptracers is not in use
C         (usePTRACERS false, or the model was built without it);
C       o the file has more of them than RNF_nTr (RNF_SIZE.h);
C       o the name after the prefix is empty;
C       o it is longer than RNF_idLen (RNF.h), so it would not fit
C         RNF_trNam.
C
C     The two name-length checks come before the matching, so they are
C     reached whatever usePTRACERS is. Five of these seven refusals
C     have a case in tests/rnf/refusal_check.py; the two that do not
C     need another build (two ptracers of the same name, and a model
C     compiled without pkg/ptracers).
C
C     Refusing an unmatched name is what keeps a tracer from being
C     silently dropped: with no match there is no tendency array to
C     add it to, so the water would arrive without it and nothing in
C     the log would say so. Setting RNF_usePtracers=.FALSE. is the
C     way to read such a file on purpose.

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#ifdef ALLOW_PTRACERS
# include "PTRACERS_SIZE.h"
# include "PTRACERS_PARAMS.h"
#endif
#include "RNF_SIZE.h"
#include "RNF.h"
#ifdef HAVE_NETCDF
# include "netcdf.inc"
#endif

C !INPUT PARAMETERS: ===================================================
C  fid                  :: NetCDF id of the open file
C  myThid               :: my thread Id number
      INTEGER fid
      INTEGER myThid

C !INPUT/OUTPUT PARAMETERS: ============================================
C  errCount             :: fatal errors of RNF_INIT_FIXED
      INTEGER errCount

#ifdef ALLOW_RNF

C !FUNCTIONS: ==========================================================
      INTEGER  ILNBLNK
      EXTERNAL ILNBLNK

C !LOCAL VARIABLES: ====================================================
C  nTrFile              :: runoff_ptracer_* variables in the file
C  nMatch               :: ptracers whose name matches this one
      INTEGER k, iL, nTrFile
      CHARACTER*(MAX_LEN_MBUF) msgBuf
#ifdef HAVE_NETCDF
      INTEGER errNC, nVars
      CHARACTER*(NF_MAX_NAME) varNam
      CHARACTER*(RNF_idLen) trNam
      INTEGER nTrLen
# ifdef ALLOW_PTRACERS
      INTEGER iTr, nMatch, iHit, nmLen, nPtr
      CHARACTER*(MAX_LEN_FNAM) pNam
# endif
#endif
CEOP

      RNF_hasTemp = .FALSE.
      RNF_hasSalt = .FALSE.
      RNF_nTrUse  = 0
      nTrFile     = 0
      DO k = 1, RNF_nTr
        RNF_trNam(k) = ' '
        RNF_trPtr(k) = 0
      ENDDO

#ifdef HAVE_NETCDF
      errNC = NF_INQ_NVARS( fid, nVars )
      CALL RNF_NC_ERROR( 'RNF_NC_SERIES',
     &     'counting the variables of the file', errNC, myThid )

      DO k = 1, nVars
        varNam = ' '
        errNC = NF_INQ_VARNAME( fid, k, varNam )
        CALL RNF_NC_ERROR( 'RNF_NC_SERIES',
     &       'reading the name of a variable', errNC, myThid )
        iL = MAX( 1, ILNBLNK(varNam) )

        IF     ( varNam(1:iL) .EQ. 'runoff_temperature' ) THEN
          IF ( RNF_useTemp ) THEN
            RNF_hasTemp = .TRUE.
          ELSE
            WRITE(msgBuf,'(2A)') 'RNF_NC_SERIES: runoff_temperature',
     &        ' is in the file and RNF_useTemp is false: not applied'
            CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                          SQUEEZE_RIGHT, myThid )
          ENDIF

        ELSEIF ( varNam(1:iL) .EQ. 'runoff_salinity' ) THEN
          IF ( RNF_useSalt ) THEN
            RNF_hasSalt = .TRUE.
          ELSE
            WRITE(msgBuf,'(2A)') 'RNF_NC_SERIES: runoff_salinity',
     &        ' is in the file and RNF_useSalt is false: not applied'
            CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                          SQUEEZE_RIGHT, myThid )
          ENDIF

        ELSEIF ( iL.GE.15 .AND.
     &           varNam(1:15) .EQ. 'runoff_ptracer_' ) THEN
          IF ( .NOT.RNF_usePtracers ) THEN
            WRITE(msgBuf,'(3A)') 'RNF_NC_SERIES: ', varNam(1:iL),
     &        ' is in the file and RNF_usePtracers is false:'
            CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                          SQUEEZE_RIGHT, myThid )
            WRITE(msgBuf,'(2A)') 'RNF_NC_SERIES: not applied, and',
     &        ' its name is not matched to the ptracers'
            CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                          SQUEEZE_RIGHT, myThid )
            GOTO 10
          ENDIF
          nTrFile = nTrFile + 1
          nTrLen = iL - 15
          IF ( nTrLen .LT. 1 ) THEN
            WRITE(msgBuf,'(3A)') 'RNF_NC_SERIES: RNF: ',
     &        varNam(1:iL), ': the tracer name after the prefix'
     &        //' runoff_ptracer_ is empty'
            CALL PRINT_ERROR( msgBuf, myThid )
            errCount = errCount + 1
            GOTO 10
          ENDIF
          IF ( nTrLen .GT. RNF_idLen ) THEN
            WRITE(msgBuf,'(3A,I6)') 'RNF_NC_SERIES: RNF: ',
     &        varNam(1:iL), ': the tracer name is longer than',
     &        RNF_idLen
            CALL PRINT_ERROR( msgBuf, myThid )
            errCount = errCount + 1
            GOTO 10
          ENDIF
          trNam = ' '
          trNam(1:nTrLen) = varNam(16:iL)
#ifdef ALLOW_PTRACERS
          IF ( .NOT.usePTRACERS ) THEN
            WRITE(msgBuf,'(3A)') 'RNF_NC_SERIES: RNF: ',
     &        varNam(1:iL), ': the file has runoff tracers but'
     &        //' pkg/ptracers is not in use'
            CALL PRINT_ERROR( msgBuf, myThid )
            WRITE(msgBuf,'(2A)') 'RNF_NC_SERIES: set usePTRACERS in',
     &        ' data.pkg, or RNF_usePtracers=.FALSE. in data.rnf'
            CALL PRINT_ERROR( msgBuf, myThid )
            errCount = errCount + 1
            GOTO 10
          ENDIF
          nMatch = 0
          iHit = 0
          DO iTr = 1, PTRACERS_numInUse
            nmLen = ILNBLNK( PTRACERS_names(iTr) )
            IF ( nmLen.EQ.nTrLen ) THEN
             IF ( PTRACERS_names(iTr)(1:nmLen).EQ.trNam(1:nTrLen) )
     &       THEN
              nMatch = nMatch + 1
              iHit = iTr
             ENDIF
            ENDIF
          ENDDO
          IF ( nMatch .EQ. 0 ) THEN
            nPtr = PTRACERS_numInUse
            WRITE(msgBuf,'(3A)') 'RNF_NC_SERIES: RNF: ',
     &        varNam(1:iL), ': no PTRACERS_names entry matches'
            CALL PRINT_ERROR( msgBuf, myThid )
            WRITE(msgBuf,'(2A,I4,A)') 'RNF_NC_SERIES: the tracer',
     &        ' name must be one of the', nPtr,
     &        ' PTRACERS_names of data.ptracers:'
            CALL PRINT_ERROR( msgBuf, myThid )
            DO iTr = 1, PTRACERS_numInUse
              pNam = PTRACERS_names(iTr)
              nmLen = MAX( 1, ILNBLNK( pNam ) )
              WRITE(msgBuf,'(A,I4,2A)') 'RNF_NC_SERIES: ptracer',
     &          iTr, ' is called ', pNam(1:nmLen)
              CALL PRINT_ERROR( msgBuf, myThid )
            ENDDO
            errCount = errCount + 1
            GOTO 10
          ELSEIF ( nMatch .GE. 2 ) THEN
            WRITE(msgBuf,'(3A,I4,A)') 'RNF_NC_SERIES: RNF: ',
     &        varNam(1:iL), ': the name matches', nMatch,
     &        ' PTRACERS_names entries'
            CALL PRINT_ERROR( msgBuf, myThid )
            errCount = errCount + 1
            GOTO 10
          ENDIF
#else /* ALLOW_PTRACERS */
          iL = MAX( 1, ILNBLNK(varNam) )
          WRITE(msgBuf,'(3A)') 'RNF_NC_SERIES: RNF: ', varNam(1:iL),
     &      ': the model was built without pkg/ptracers'
          CALL PRINT_ERROR( msgBuf, myThid )
          errCount = errCount + 1
          GOTO 10
#endif /* ALLOW_PTRACERS */
          IF ( RNF_nTrUse .LT. RNF_nTr ) THEN
            RNF_nTrUse = RNF_nTrUse + 1
            RNF_trNam(RNF_nTrUse) = trNam
#ifdef ALLOW_PTRACERS
            RNF_trPtr(RNF_nTrUse) = iHit
#endif
          ENDIF
        ENDIF
   10   CONTINUE
      ENDDO

      IF ( nTrFile .GT. RNF_nTr ) THEN
        WRITE(msgBuf,'(A,I6,A,I6,A)') 'RNF_NC_SERIES: RNF: RNF_nTr =',
     &    RNF_nTr, ' is too small: the file has', nTrFile,
     &    ' runoff tracer(s)'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
      ENDIF
#endif /* HAVE_NETCDF */

#endif /* ALLOW_RNF */

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: RNF_NC_READ_FLUX

C !INTERFACE: ==========================================================
      SUBROUTINE RNF_NC_READ_FLUX(
     I                   fileName, iRec, iBuf,
     I                   myThid )

C !DESCRIPTION:
C     Read record iRec of every time series the file carries from
C     fileName into record buffer iBuf (1 or 2) of RNF.h, keeping for
C     every tile of this process the value of each source that tile
C     needs:
C
C       runoff_flux            -> RNF_bufFlux  [m^3/s]     (required)
C       runoff_temperature     -> RNF_bufTemp  [degC]      (optional)
C       runoff_salinity        -> RNF_bufSalt  [g/kg]      (optional)
C       runoff_ptracer_<NAME>  -> RNF_bufTrc               (optional)
C
C     Which of the optional ones exist and are used was settled at
C     init by RNF_NC_SERIES (RNF_hasTemp, RNF_hasSalt, RNF_nTrUse,
C     RNF_trNam of RNF.h). RNF_bufSum(iBuf) returns the flux summed
C     over all the sources of the file. Every series is read
C     RNF_nBuf sources at a time, so no array has the length of the
C     source dimension.
C
C     fileName is RNF_file, or the file of one year of a yearly set,
C     as RNF_FILE_NAME builds it.
C
C     A record outside the time dimension of this file stops the run.
C     That is the file being too short for the run: a non-repeating
C     series ends, and the step after its last record asks for a
C     record that does not exist. NetCDF would refuse the read with
C     "index exceeds dimension bound", which says nothing about
C     runoff, so the record, the count and the file are named here.
C
C     A missing value stops the run for every series except the
C     temperature (model contract; runoff schema 1.0, section 3.5):
C     not a number, above RNF_fluxMax in absolute value, or equal to
C     the _FillValue or missing_value of the variable. Such an
C     attribute that exists but is not one number also stops the run,
C     in RNF_NC_ATT_REAL. A missing temperature is allowed and means
C     "this source enters at the reference temperature"; it is
C     recorded in RNF_bufTvld and leaves the source out of (mT) and
C     m_T (package design, decision 3). Each series reports its own
C     missing values and how many of them in RNF_NC_READ_ONE; this
C     routine only ends the run afterwards, so that one read reports
C     every series rather than stopping at the first.
C
C     Called by the master thread, which opens and closes the file.

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "RNF_SIZE.h"
#include "RNF.h"
#ifdef HAVE_NETCDF
# include "netcdf.inc"
#endif

C !INPUT PARAMETERS: ===================================================
C  fileName             :: file of the set to read from
C  iRec                 :: time record to read (from 1)
C  iBuf                 :: record buffer of RNF.h to fill (1 or 2)
C  myThid               :: my thread Id number
      CHARACTER*(*) fileName
      INTEGER iRec
      INTEGER iBuf
      INTEGER myThid

C !OUTPUT PARAMETERS: ==================================================
C  none (the record buffers iBuf of RNF.h)

#ifdef ALLOW_RNF

C !FUNCTIONS: ==========================================================
      INTEGER  ILNBLNK
      EXTERNAL ILNBLNK

C !LOCAL VARIABLES: ====================================================
C  nBad                 :: missing values found (not allowed)
C  nMsg                 :: messages printed about them
C  varNam               :: name of the series being read
C  dumSum               :: sum of a series that has no use for one
      INTEGER bi, bj, k, n, iL, iT
      INTEGER nBad, nMsg
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      CHARACTER*(80) varNam
      _RL dumSum
#ifdef HAVE_NETCDF
      INTEGER errNC, fid
      INTEGER timDimId, nRecLoc
#endif
CEOP

      nBad = 0
      nMsg = 0
      RNF_bufSum(iBuf) = 0. _d 0
      DO bj = 1, nSy
       DO bi = 1, nSx
        DO k = 1, RNF_nSrcTile
          RNF_bufFlux(k,bi,bj,iBuf) = 0. _d 0
          RNF_bufTemp(k,bi,bj,iBuf) = 0. _d 0
          RNF_bufTvld(k,bi,bj,iBuf) = 0. _d 0
          RNF_bufSalt(k,bi,bj,iBuf) = 0. _d 0
        ENDDO
        DO iT = 1, RNF_nTr
         DO k = 1, RNF_nSrcTile
          RNF_bufTrc(k,bi,bj,iBuf,iT) = 0. _d 0
         ENDDO
        ENDDO
       ENDDO
      ENDDO

#ifdef HAVE_NETCDF
      iL = MAX( 1, ILNBLNK(fileName) )
      WRITE(msgBuf,'(2A)') 'opening ', fileName(1:iL)
      errNC = NF_OPEN( fileName(1:iL), NF_NOWRITE, fid )
      CALL RNF_NC_ERROR( 'RNF_NC_READ_FLUX', msgBuf, errNC, myThid )

C--   The record has to be in this file. A non-repeating series that
C     ends before the run does would otherwise be read past its end.
      errNC = NF_INQ_DIMID( fid, 'time', timDimId )
      CALL RNF_NC_ERROR( 'RNF_NC_READ_FLUX',
     &     'looking for the dimension time', errNC, myThid )
      errNC = NF_INQ_DIMLEN( fid, timDimId, nRecLoc )
      CALL RNF_NC_ERROR( 'RNF_NC_READ_FLUX',
     &     'reading the length of the dimension time',
     &     errNC, myThid )
      IF ( iRec.LT.1 .OR. iRec.GT.nRecLoc ) THEN
        WRITE(msgBuf,'(A,I9,A,I9,2A)') 'RNF_NC_READ_FLUX: RNF:'
     &    //' record', iRec, ' is outside the', nRecLoc,
     &    ' record(s) of ', fileName(1:iL)
        CALL PRINT_ERROR( msgBuf, myThid )
        WRITE(msgBuf,'(2A)') 'RNF_NC_READ_FLUX: the series does not',
     &    ' cover the model time; it needs more records,'
        CALL PRINT_ERROR( msgBuf, myThid )
        WRITE(msgBuf,'(2A)') 'RNF_NC_READ_FLUX: a repeat cycle',
     &    ' (mitgcm_time_repeat) or another RNF_period'
        CALL PRINT_ERROR( msgBuf, myThid )
        CALL ALL_PROC_DIE( 0 )
        STOP 'ABNORMAL END: S/R RNF_NC_READ_FLUX'
      ENDIF

C--   The volume flux, which every file has and no file may leave out
      CALL RNF_NC_READ_ONE( fid, iRec, 'runoff_flux', .FALSE.,
     I     'RNF_NC_READ_FLUX',
     O     RNF_bufFlux(1,1,1,iBuf), RNF_vldWrk, RNF_bufSum(iBuf),
     U     nBad, nMsg,
     I     myThid )

C--   The temperature, the only series that may have missing values.
C     RNF_vldWrk is work space shared by every series, so the flags
C     are copied out right after the read.
      IF ( RNF_hasTemp ) THEN
        CALL RNF_NC_READ_ONE( fid, iRec, 'runoff_temperature', .TRUE.,
     I       'RNF_NC_READ_FLUX',
     O       RNF_bufTemp(1,1,1,iBuf), RNF_vldWrk, dumSum,
     U       nBad, nMsg,
     I       myThid )
        DO bj = 1, nSy
         DO bi = 1, nSx
          DO k = 1, RNF_nSrc(bi,bj)
            RNF_bufTvld(k,bi,bj,iBuf) = RNF_vldWrk(k,bi,bj)
          ENDDO
         ENDDO
        ENDDO
      ENDIF

      IF ( RNF_hasSalt ) THEN
        CALL RNF_NC_READ_ONE( fid, iRec, 'runoff_salinity', .FALSE.,
     I       'RNF_NC_READ_FLUX',
     O       RNF_bufSalt(1,1,1,iBuf), RNF_vldWrk, dumSum,
     U       nBad, nMsg,
     I       myThid )
      ENDIF

      DO iT = 1, RNF_nTrUse
        n = MAX( 1, ILNBLNK(RNF_trNam(iT)) )
        WRITE(varNam,'(2A)') 'runoff_ptracer_', RNF_trNam(iT)(1:n)
        CALL RNF_NC_READ_ONE( fid, iRec, varNam, .FALSE.,
     I       'RNF_NC_READ_FLUX',
     O       RNF_bufTrc(1,1,1,iBuf,iT), RNF_vldWrk, dumSum,
     U       nBad, nMsg,
     I       myThid )
      ENDDO

      errNC = NF_CLOSE( fid )
      CALL RNF_NC_ERROR( 'RNF_NC_READ_FLUX', 'closing the file',
     &     errNC, myThid )

C--   Each series has already said how many values of its own were
C     missing, and which, so this only ends the run.
      IF ( nBad .GE. 1 ) THEN
        CALL ALL_PROC_DIE( 0 )
        STOP 'ABNORMAL END: S/R RNF_NC_READ_FLUX'
      ENDIF
#endif /* HAVE_NETCDF */

#endif /* ALLOW_RNF */

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: RNF_NC_READ_ONE

C !INTERFACE: ==========================================================
      SUBROUTINE RNF_NC_READ_ONE(
     I                   fid, iRec, varNam, allowMiss, caller,
     O                   srcVal, srcVld, valSum,
     U                   nBad, nMsg,
     I                   myThid )

C !DESCRIPTION:
C     Read record iRec of one (time, source) series of the open file
C     fid and store, for every tile of this process, the value of each
C     source that tile needs. Returns the value summed over all the
C     sources of the file, which only the flux has a use for.
C
C     The record is read RNF_nBuf sources at a time, so no array has
C     the length of the source dimension. The variable must have the
C     dimensions (time, source), i.e. (source, time) in Fortran order,
C     so that one record is one contiguous hyperslab.
C
C     A value is missing when it is not a number, when it is above
C     RNF_fluxMax in absolute value, or when it equals the _FillValue
C     or missing_value of the variable. srcVld is set to 1 for a
C     value that is present and to 0 for one that is missing, and a
C     missing value is stored as 0 so that a fill value cannot
C     propagate as a NaN. With allowMiss false, each missing value
C     also raises nBad, and the caller stops the run; the flux and
C     everything but the temperature is in that class (model
C     contract; runoff schema 1.0, section 3.5).
C
C     nBad and nMsg are carried across the series of one record, so
C     that one read reports every missing value it found and one
C     process still prints at most RNF_maxErrMsg of them. The count
C     line this routine prints at the end is per series and names it:
C     one read covers every series of a record, so a count of the
C     record as a whole would say how many values were lost and not of
C     what. The caller stops the run once every series has reported.
C
C     Called by the master thread, with the file already open.

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "RNF_SIZE.h"
#include "RNF.h"
#ifdef HAVE_NETCDF
# include "netcdf.inc"
#endif

C !INPUT PARAMETERS: ===================================================
C  fid                  :: NetCDF id of the open file
C  iRec                 :: time record to read (from 1)
C  varNam               :: name of the series
C  allowMiss            :: a missing value is legal for this series
C  caller               :: routine named in the messages
C  myThid               :: my thread Id number
      INTEGER fid
      INTEGER iRec
      CHARACTER*(*) varNam
      LOGICAL allowMiss
      CHARACTER*(*) caller
      INTEGER myThid

C !OUTPUT PARAMETERS: ==================================================
C  srcVal               :: per-tile value of each source the tile needs
C  srcVld               :: 1 where that value was present, 0 where it
C                          was missing
C  valSum               :: value summed over the sources of the file
C  nBad                 :: raised by each missing value that is not
C                          allowed
C  nMsg                 :: raised by each message printed
      _RL srcVal(RNF_nSrcTile,nSx,nSy)
      _RL srcVld(RNF_nSrcTile,nSx,nSy)
      _RL valSum
      INTEGER nBad
      INTEGER nMsg

#ifdef ALLOW_RNF

C !FUNCTIONS: ==========================================================
      INTEGER  ILNBLNK
      EXTERNAL ILNBLNK

C !LOCAL VARIABLES: ====================================================
      INTEGER bi, bj, k, n, iV
      INTEGER iSrc0, nLoc, nChunk, iChunk, nBad0
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      CHARACTER*(RNF_idLen) srcId
      INTEGER idLen
      _RL val
      LOGICAL isBad
#ifdef HAVE_NETCDF
      INTEGER errNC, varId, nDims
      INTEGER dimIds(NF_MAX_VAR_DIMS)
      INTEGER srcDimId, timDimId
      INTEGER iStart(2), iCount(2)
      LOGICAL hasFill, hasMiss
      _RL fillVal, missVal
#endif
CEOP

      valSum = 0. _d 0
      nBad0 = nBad
      iV = MAX( 1, ILNBLNK(varNam) )
      DO bj = 1, nSy
       DO bi = 1, nSx
        DO k = 1, RNF_nSrcTile
          srcVal(k,bi,bj) = 0. _d 0
          srcVld(k,bi,bj) = 0. _d 0
        ENDDO
       ENDDO
      ENDDO

#ifdef HAVE_NETCDF
      WRITE(msgBuf,'(2A)') 'looking for the required variable ',
     &  varNam(1:iV)
      errNC = NF_INQ_VARID( fid, varNam(1:iV), varId )
      CALL RNF_NC_ERROR( caller, msgBuf, errNC, myThid )

C--   (time, source), i.e. (source, time) in Fortran order
      errNC = NF_INQ_DIMID( fid, 'source', srcDimId )
      CALL RNF_NC_ERROR( caller,
     &     'looking for the dimension source', errNC, myThid )
      errNC = NF_INQ_DIMID( fid, 'time', timDimId )
      CALL RNF_NC_ERROR( caller,
     &     'looking for the dimension time', errNC, myThid )
      WRITE(msgBuf,'(2A)') 'reading the shape of ', varNam(1:iV)
      errNC = NF_INQ_VARNDIMS( fid, varId, nDims )
      CALL RNF_NC_ERROR( caller, msgBuf, errNC, myThid )
      IF ( nDims .EQ. 2 ) THEN
        errNC = NF_INQ_VARDIMID( fid, varId, dimIds )
        CALL RNF_NC_ERROR( caller, msgBuf, errNC, myThid )
      ELSE
        dimIds(1) = -1
        dimIds(2) = -1
      ENDIF
      IF ( nDims.NE.2 .OR. dimIds(1).NE.srcDimId
     &     .OR. dimIds(2).NE.timDimId ) THEN
        WRITE(msgBuf,'(4A)') caller, ': RNF: ', varNam(1:iV),
     &    ' must have the dimensions (time, source)'
        CALL PRINT_ERROR( msgBuf, myThid )
        CALL ALL_PROC_DIE( 0 )
        STOP 'ABNORMAL END: S/R RNF_NC_READ_ONE'
      ENDIF

C--   markers of a missing value; absent is the normal case
      CALL RNF_NC_ATT_REAL( fid, varId, varNam(1:iV), '_FillValue',
     I     caller,
     O     fillVal, hasFill,
     I     myThid )
      CALL RNF_NC_ATT_REAL( fid, varId, varNam(1:iV),
     I     'missing_value', caller,
     O     missVal, hasMiss,
     I     myThid )

C--   one chunk of sources at a time
      nChunk = 1 + (RNF_nSrcFile-1)/RNF_nBuf
      DO iChunk = 1, nChunk
        iSrc0 = (iChunk-1)*RNF_nBuf
        nLoc = MIN( RNF_nBuf, RNF_nSrcFile - iSrc0 )
        iStart(1) = iSrc0 + 1
        iStart(2) = iRec
        iCount(1) = nLoc
        iCount(2) = 1
        errNC = NF_GET_VARA_DOUBLE( fid, varId, iStart, iCount,
     &                              RNF_fluxBuf )
        WRITE(msgBuf,'(2A)') 'reading a record of ', varNam(1:iV)
        CALL RNF_NC_ERROR( caller, msgBuf, errNC, myThid )

        DO n = 1, nLoc
          val = RNF_fluxBuf(n)
C-        not a number, or too large to be a physical value
          isBad = val .NE. val
          IF ( .NOT.isBad ) isBad = ABS(val) .GT. RNF_fluxMax
          IF ( .NOT.isBad .AND. hasFill ) isBad = val .EQ. fillVal
          IF ( .NOT.isBad .AND. hasMiss ) isBad = val .EQ. missVal
          RNF_vldBuf(n) = 1. _d 0
          IF ( isBad ) THEN
            RNF_fluxBuf(n) = 0. _d 0
            RNF_vldBuf(n) = 0. _d 0
            IF ( .NOT.allowMiss ) THEN
              nBad = nBad + 1
              IF ( nMsg .LT. RNF_maxErrMsg ) THEN
                nMsg = nMsg + 1
                CALL RNF_NC_SOURCE_ID( fid, iSrc0+n-1, caller,
     O               srcId, idLen,
     I               myThid )
                WRITE(msgBuf,'(6A,I8,A,1PE16.8)')
     &            caller, ': RNF: missing ', varNam(1:iV),
     &            ': source ', srcId(1:idLen), ', record', iRec,
     &            ', value', val
                CALL PRINT_ERROR( msgBuf, myThid )
              ENDIF
            ENDIF
          ELSE
            valSum = valSum + val
          ENDIF
        ENDDO

C-      hand the chunk to every tile that needs part of it
        DO bj = 1, nSy
         DO bi = 1, nSx
          DO k = 1, RNF_nSrc(bi,bj)
           n = RNF_srcGlob(k,bi,bj) - iSrc0 + 1
           IF ( n.GE.1 .AND. n.LE.nLoc ) THEN
             srcVal(k,bi,bj) = RNF_fluxBuf(n)
             srcVld(k,bi,bj) = RNF_vldBuf(n)
           ENDIF
          ENDDO
         ENDDO
        ENDDO
      ENDDO

C--   How many values of THIS series were missing when they may not be.
C     The count is per series, and names it, because one read covers
C     every series of a record: a count of the record as a whole would
C     say how many values were lost and not of what.
      IF ( nBad .GT. nBad0 ) THEN
        WRITE(msgBuf,'(2A,I8,3A)') caller, ': RNF:', nBad-nBad0,
     &    ' missing value(s) of ', varNam(1:iV), ' (not allowed)'
        CALL PRINT_ERROR( msgBuf, myThid )
      ENDIF
#endif /* HAVE_NETCDF */

#endif /* ALLOW_RNF */

      RETURN
      END
