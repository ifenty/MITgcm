#include "RNF_OPTIONS.h"

C--  File rnf_nc_utils.F: NetCDF helpers of the sparse runoff package
C--   Contents
C--   o RNF_NC_ERROR     :: stop the run on a NetCDF error
C--   o RNF_NC_ATT_TEXT  :: read a text attribute
C--   o RNF_NC_ATT_REAL  :: read an attribute that holds one number
C--   o RNF_NC_SOURCE_ID :: read the id of one source
C--   o RNF_NC_TGT_ERROR :: print one error line about a target entry
C--   o RNF_NC_READ_FLUX :: read one record of runoff_flux
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
        WRITE(msgBuf,'(A,2(A,I6))') 'RNF_NC_ATT_REAL: it has type',
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
      WRITE(msgBuf,'(5A,I12,A,4I6)') caller(1:iL), ': RNF: ',
     &     text(1:jL), ': source ', srcId(1:idLen),
     &     ', target_cell', cell, ', i,j,bi,bj =', i, j, bi, bj
      CALL PRINT_ERROR( msgBuf, myThid )

#endif /* ALLOW_RNF */

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: RNF_NC_READ_FLUX

C !INTERFACE: ==========================================================
      SUBROUTINE RNF_NC_READ_FLUX(
     I                   iRec,
     O                   fluxSum,
     I                   myThid )

C !DESCRIPTION:
C     Read record iRec of runoff_flux (m^3/s) and store, for every
C     tile of this process, the flux of each source that tile needs
C     (RNF_srcFlux). Returns the flux summed over all the sources of
C     the file. The record is read RNF_nBuf sources at a time, so no
C     array has the length of the source dimension.
C
C     A missing flux stops the run (model contract): not a number,
C     above RNF_fluxMax in absolute value, or equal to the _FillValue
C     or missing_value of the variable. Such an attribute that exists
C     but is not one number also stops the run, in RNF_NC_ATT_REAL.
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
C  iRec                 :: time record to read (from 1)
C  myThid               :: my thread Id number
      INTEGER iRec
      INTEGER myThid

C !OUTPUT PARAMETERS: ==================================================
C  fluxSum              :: flux summed over the sources of the file
      _RL fluxSum

#ifdef ALLOW_RNF

C !FUNCTIONS: ==========================================================
      INTEGER  ILNBLNK
      EXTERNAL ILNBLNK

C !LOCAL VARIABLES: ====================================================
C  nBad                 :: missing values found
C  nMsg                 :: messages printed about them
      INTEGER bi, bj, k, n, iL
      INTEGER iSrc0, nLoc, nBad, nMsg, nChunk, iChunk
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      CHARACTER*(RNF_idLen) srcId
      INTEGER idLen
      _RL fluxVal
      LOGICAL isBad
#ifdef HAVE_NETCDF
      INTEGER errNC, fid, varId, nDims
      INTEGER dimIds(NF_MAX_VAR_DIMS)
      INTEGER srcDimId, timDimId
      INTEGER iStart(2), iCount(2)
      LOGICAL hasFill, hasMiss
      _RL fillVal, missVal
#endif
CEOP

      fluxSum = 0. _d 0
      nBad = 0
      nMsg = 0
      DO bj = 1, nSy
       DO bi = 1, nSx
        DO k = 1, RNF_nSrcTile
          RNF_srcFlux(k,bi,bj) = 0. _d 0
        ENDDO
       ENDDO
      ENDDO

#ifdef HAVE_NETCDF
      iL = MAX( 1, ILNBLNK(RNF_file) )
      errNC = NF_OPEN( RNF_file(1:iL), NF_NOWRITE, fid )
      CALL RNF_NC_ERROR( 'RNF_NC_READ_FLUX', 'opening the file',
     &     errNC, myThid )

      errNC = NF_INQ_VARID( fid, 'runoff_flux', varId )
      CALL RNF_NC_ERROR( 'RNF_NC_READ_FLUX',
     &     'looking for the required variable runoff_flux',
     &     errNC, myThid )

C--   runoff_flux must be (time, source), i.e. (source, time) in
C     Fortran order, so that one record is one contiguous hyperslab
      errNC = NF_INQ_DIMID( fid, 'source', srcDimId )
      CALL RNF_NC_ERROR( 'RNF_NC_READ_FLUX',
     &     'looking for the dimension source', errNC, myThid )
      errNC = NF_INQ_DIMID( fid, 'time', timDimId )
      CALL RNF_NC_ERROR( 'RNF_NC_READ_FLUX',
     &     'looking for the dimension time', errNC, myThid )
      errNC = NF_INQ_VARNDIMS( fid, varId, nDims )
      CALL RNF_NC_ERROR( 'RNF_NC_READ_FLUX',
     &     'reading the shape of runoff_flux', errNC, myThid )
      IF ( nDims .EQ. 2 ) THEN
        errNC = NF_INQ_VARDIMID( fid, varId, dimIds )
        CALL RNF_NC_ERROR( 'RNF_NC_READ_FLUX',
     &       'reading the shape of runoff_flux', errNC, myThid )
      ELSE
        dimIds(1) = -1
        dimIds(2) = -1
      ENDIF
      IF ( nDims.NE.2 .OR. dimIds(1).NE.srcDimId
     &     .OR. dimIds(2).NE.timDimId ) THEN
        WRITE(msgBuf,'(2A)') 'RNF_NC_READ_FLUX: RNF: runoff_flux',
     &    ' must have the dimensions (time, source)'
        CALL PRINT_ERROR( msgBuf, myThid )
        CALL ALL_PROC_DIE( 0 )
        STOP 'ABNORMAL END: S/R RNF_NC_READ_FLUX'
      ENDIF

C--   markers of a missing value; absent is the normal case
      CALL RNF_NC_ATT_REAL( fid, varId, 'runoff_flux', '_FillValue',
     I     'RNF_NC_READ_FLUX',
     O     fillVal, hasFill,
     I     myThid )
      CALL RNF_NC_ATT_REAL( fid, varId, 'runoff_flux',
     I     'missing_value', 'RNF_NC_READ_FLUX',
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
        CALL RNF_NC_ERROR( 'RNF_NC_READ_FLUX',
     &       'reading a record of runoff_flux', errNC, myThid )

        DO n = 1, nLoc
          fluxVal = RNF_fluxBuf(n)
C-        not a number, or too large to be a flux
          isBad = fluxVal .NE. fluxVal
          IF ( .NOT.isBad ) isBad = ABS(fluxVal) .GT. RNF_fluxMax
          IF ( .NOT.isBad .AND. hasFill )
     &         isBad = fluxVal .EQ. fillVal
          IF ( .NOT.isBad .AND. hasMiss )
     &         isBad = fluxVal .EQ. missVal
          IF ( isBad ) THEN
            nBad = nBad + 1
            IF ( nMsg .LT. RNF_maxErrMsg ) THEN
              nMsg = nMsg + 1
              CALL RNF_NC_SOURCE_ID( fid, iSrc0+n-1,
     I             'RNF_NC_READ_FLUX',
     O             srcId, idLen,
     I             myThid )
              WRITE(msgBuf,'(3A,I8,A,1PE16.8)')
     &          'RNF_NC_READ_FLUX: RNF: missing runoff_flux: source ',
     &          srcId(1:idLen), ', record', iRec, ', value',
     &          fluxVal
              CALL PRINT_ERROR( msgBuf, myThid )
            ENDIF
          ELSE
            fluxSum = fluxSum + fluxVal
          ENDIF
        ENDDO

C-      hand the chunk to every tile that needs part of it
        DO bj = 1, nSy
         DO bi = 1, nSx
          DO k = 1, RNF_nSrc(bi,bj)
           n = RNF_srcGlob(k,bi,bj) - iSrc0 + 1
           IF ( n.GE.1 .AND. n.LE.nLoc )
     &          RNF_srcFlux(k,bi,bj) = RNF_fluxBuf(n)
          ENDDO
         ENDDO
        ENDDO
      ENDDO

      errNC = NF_CLOSE( fid )
      CALL RNF_NC_ERROR( 'RNF_NC_READ_FLUX', 'closing the file',
     &     errNC, myThid )

      IF ( nBad .GE. 1 ) THEN
        WRITE(msgBuf,'(A,I8,A)') 'RNF_NC_READ_FLUX: RNF:', nBad,
     &    ' missing value(s) of runoff_flux (not allowed)'
        CALL PRINT_ERROR( msgBuf, myThid )
        CALL ALL_PROC_DIE( 0 )
        STOP 'ABNORMAL END: S/R RNF_NC_READ_FLUX'
      ENDIF
#endif /* HAVE_NETCDF */

#endif /* ALLOW_RNF */

      RETURN
      END
