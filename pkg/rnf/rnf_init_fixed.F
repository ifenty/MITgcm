#include "RNF_OPTIONS.h"

CBOP
C !ROUTINE: RNF_INIT_FIXED

C !INTERFACE: ==========================================================
      SUBROUTINE RNF_INIT_FIXED( myThid )

C !DESCRIPTION:
C     Static read of the sparse runoff file, placement of its target
C     cells on the tiles of this process and the file checks
C     (package design, decision 6). The master thread of every
C     process reads the file, RNF_nBuf entries at a time, and keeps
C     only what its own tiles own.
C
C     Steps, in the order they are made:
C     1. the global attributes: schema major version, time sampling,
C        mitgcm_grid_nx and mitgcm_grid_ny against the model's own
C        global I/O layout;
C     2. the dimensions time, source and target, and the required
C        variables;
C     3. the target table in chunks. Each entry is checked before it
C        is placed: target_source and target_cell in range,
C        target_level = 1 and target_fraction in [0,1];
C     4. placement with the arithmetic that pkg/mdsio uses to read a
C        global file (mdsio_read_field.F), never from the grid
C        geometry: row j of tile (bi,bj) is the global record
C        (tBx + (j-1)*iGjLoc)/sNx + (tBy + (j-1)*jGjLoc)*(xSize/sNx),
C        with tBx, tBy, iGjLoc and jGjLoc from myXGlobalLo and
C        myYGlobalLo, or from the exch2 I/O variables including the
C        fold and long-line cases. A target whose record is a record
C        of the tile belongs to it;
C     5. the checks of an owned target: on land, beyond an open
C        boundary (maskInC), under an ice shelf (kTopC) and, when the
C        file has target_cell_area, against rA;
C     6. the fraction sum of every source over all tiles and
C        processes (GLOBAL_SUM_VECTOR_RL, one chunk of sources at a
C        time). This also catches a target that no tile owns, when
C        its fraction exceeds the tolerance: on a blank exch2 tile,
C        or on a cell of the global layout that no facet uses.
C
C     Errors that only the process owning the tile can see (a refused
C     target, an array bound of RNF_SIZE.h that is too small) are
C     counted, the counts are summed over all processes with
C     GLOBAL_SUM_INT, and then every process stops. ALL_PROC_DIE ends
C     MPI on the calling process only, so a stop on one process alone
C     would hang the others.
C
C     This routine is called before RNF_CHECK (S/R INITIALISE_FIXED
C     calls PACKAGES_INIT_FIXED before PACKAGES_CHECK), so what the
C     static read relies on (a non-blank RNF_file, a build with
C     NetCDF) is checked in RNF_READPARMS.

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#ifdef ALLOW_EXCH2
# include "W2_EXCH2_SIZE.h"
# include "W2_EXCH2_TOPOLOGY.h"
# include "W2_EXCH2_PARAMS.h"
#endif /* ALLOW_EXCH2 */
#ifdef ALLOW_SHELFICE
# include "SHELFICE.h"
#endif /* ALLOW_SHELFICE */
#include "RNF_SIZE.h"
#include "RNF.h"
#ifdef HAVE_NETCDF
# include "netcdf.inc"
#endif

C !INPUT PARAMETERS: ===================================================
C  myThid               :: my thread Id number
      INTEGER myThid

C !OUTPUT PARAMETERS: ==================================================
C  none

#ifdef ALLOW_RNF

C !FUNCTIONS: ==========================================================
      INTEGER  ILNBLNK
      EXTERNAL ILNBLNK

C !LOCAL VARIABLES: ====================================================
C  errCount             :: fatal errors this process reports
C  nErrEnt              :: invalid entries of the target table
C  nErrFrac             :: sources whose fractions do not sum to 1
C  nErrTgt              :: refused targets (summed over processes)
C  nErrBnd              :: tiles needing a larger array bound (idem)
C  nTgtOwn              :: target entries owned (summed over processes)
C  rowRec               :: global record of each row of each tile
C  nSrcNeed, nTgtNeed   :: entries a tile would need per array bound
C  lastOvf              :: last source that did not fit, so that a
C                          file sorted by source counts it once
      INTEGER errCount
      INTEGER nErrEnt, nErrFrac, nErrTgt, nErrBnd, nTgtOwn
      INTEGER nMsgEnt, nMsgTgt
      INTEGER bi, bj, i, j, k, n, iL
      INTEGER kSurf
      INTEGER xSize, ySize, global_nTx
      INTEGER tBx, tBy, iGjLoc, jGjLoc
      LOGICAL useExch2ioLayOut
      INTEGER rowRec(sNy,nSx,nSy)
      INTEGER nSrcNeed(nSx,nSy), nTgtNeed(nSx,nSy)
      INTEGER lastOvf(nSx,nSy)
      INTEGER fileNx, fileNy, majVer, iDot, ios
      INTEGER iChunk, nChunk, iEnt0, nLoc, kEnt
      INTEGER iSrcF, cellF, levF, iRec, iLoc, kSrc
      INTEGER refuse
      _RL     fracF, areaF
      LOGICAL badEnt, placed, found, foundNx, foundNy
      LOGICAL fileOpen
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      CHARACTER*(80) attTxt
      INTEGER attLen
      _RL     attVal
      CHARACTER*(RNF_idLen) srcId
      INTEGER idLen
#ifdef ALLOW_EXCH2
      INTEGER tN
#endif /* ALLOW_EXCH2 */
#ifdef HAVE_NETCDF
      INTEGER errNC, fid, dimId
      INTEGER idSrc, idCell, idFrac, idLev, idArea
      INTEGER iStart(1), iCount(1)
      INTEGER tgtSrcBuf(RNF_nBuf)
      INTEGER tgtCellBuf(RNF_nBuf)
      INTEGER tgtLevBuf(RNF_nBuf)
      _RL     tgtFracBuf(RNF_nBuf)
      _RL     tgtAreaBuf(RNF_nBuf)
      LOGICAL hasLev, hasArea
#endif /* HAVE_NETCDF */
CEOP

C--   Counters that have to be summed over all processes: every
C     thread contributes its own value to GLOBAL_SUM_INT, so they
C     start at zero on every thread and only the master thread of
C     each process raises them.
      nErrTgt = 0
      nErrBnd = 0
      nTgtOwn = 0
      errCount = 0
      nErrEnt = 0
      nErrFrac = 0
      nMsgEnt = 0
      nMsgTgt = 0
      fileOpen = .FALSE.

      _BEGIN_MASTER(myThid)

C--   Empty lists
      DO bj = 1, nSy
       DO bi = 1, nSx
        RNF_nSrc(bi,bj) = 0
        RNF_nTgt(bi,bj) = 0
        nSrcNeed(bi,bj) = 0
        nTgtNeed(bi,bj) = 0
        lastOvf(bi,bj) = -1
        DO k = 1, RNF_nSrcTile
          RNF_srcGlob(k,bi,bj) = -1
          RNF_srcFracSum(k,bi,bj) = 0. _d 0
          RNF_srcFlux(k,bi,bj) = 0. _d 0
        ENDDO
       ENDDO
      ENDDO
      RNF_nSrcFile = 0
      RNF_nTgtFile = 0
      RNF_nRecFile = 0
      RNF_nSrcProc = 0
      RNF_nTgtProc = 0
      RNF_nTgtOwned = 0
      RNF_tableRead = .FALSE.

C--   The surface level of a column, as S/R APPLY_FORCING_T selects
C     it (1 in z coordinates, Nr in pressure coordinates, kSurfC with
C     pkg/shelfice). Schema 1.0 allows only the surface.
      kSurf = 1
      IF ( usingPCoords ) kSurf = Nr

C--   Global I/O layout, as S/R MDSIO_READ_FIELD sets it
      xSize = Nx
      ySize = Ny
      useExch2ioLayOut = .FALSE.
#ifdef ALLOW_EXCH2
      IF ( W2_useE2ioLayOut ) THEN
        xSize = exch2_global_Nx
        ySize = exch2_global_Ny
        useExch2ioLayOut = .TRUE.
      ENDIF
#endif /* ALLOW_EXCH2 */
      global_nTx = xSize/sNx

C--   Time handling that this reader does not have yet
      IF ( RNF_useYearlyFiles ) THEN
        WRITE(msgBuf,'(2A)') 'RNF_INIT_FIXED: RNF:',
     &    ' RNF_useYearlyFiles=.TRUE. cannot be read yet'
        CALL PRINT_ERROR( msgBuf, myThid )
        WRITE(msgBuf,'(2A)') 'RNF_INIT_FIXED: yearly files come',
     &    ' with the time handling (RUNOFF-005)'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
      ENDIF
      IF ( RNF_period .NE. UNSET_RL .AND.
     &     RNF_period .NE. 0. _d 0 ) THEN
        WRITE(msgBuf,'(2A,1PE16.8,A)') 'RNF_INIT_FIXED: RNF:',
     &    ' RNF_period =', RNF_period, ' cannot be read yet'
        CALL PRINT_ERROR( msgBuf, myThid )
        WRITE(msgBuf,'(2A)') 'RNF_INIT_FIXED: only a constant',
     &    ' record is read; time records come with RUNOFF-005'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
      ENDIF
      IF ( errCount .GE. 1 ) GOTO 900

#ifdef HAVE_NETCDF

C--   1. Open the file and read the global attributes
      iL = MAX( 1, ILNBLNK(RNF_file) )
      errNC = NF_OPEN( RNF_file(1:iL), NF_NOWRITE, fid )
      CALL RNF_NC_ERROR( 'RNF_INIT_FIXED', 'opening the file',
     &     errNC, myThid )
      fileOpen = .TRUE.

      CALL RNF_NC_ATT_TEXT( fid, NF_GLOBAL,
     I     'mitgcm_runoff_schema_version', 'RNF_INIT_FIXED',
     O     attTxt, attLen, found,
     I     myThid )
      IF ( .NOT.found ) THEN
        WRITE(msgBuf,'(2A)') 'RNF_INIT_FIXED: RNF: the file has no',
     &    ' global attribute mitgcm_runoff_schema_version'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
      ELSE
        iDot = INDEX( attTxt(1:attLen), '.' )
        IF ( iDot .EQ. 0 ) iDot = attLen + 1
        majVer = -1
        ios = 1
        IF ( iDot .GE. 2 ) READ(attTxt(1:iDot-1),*,IOSTAT=ios) majVer
        IF ( ios.NE.0 .OR. majVer.NE.1 ) THEN
          WRITE(msgBuf,'(4A)') 'RNF_INIT_FIXED: RNF: the schema',
     &      ' major version of the file (', attTxt(1:attLen),
     &      ') is not 1'
          CALL PRINT_ERROR( msgBuf, myThid )
          errCount = errCount + 1
        ENDIF
      ENDIF

      CALL RNF_NC_ATT_TEXT( fid, NF_GLOBAL, 'mitgcm_time_sampling',
     I     'RNF_INIT_FIXED',
     O     attTxt, attLen, found,
     I     myThid )
      IF ( .NOT.found ) THEN
        WRITE(msgBuf,'(2A)') 'RNF_INIT_FIXED: RNF: the file has no',
     &    ' global attribute mitgcm_time_sampling'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
      ELSEIF ( attTxt(1:attLen) .NE. 'constant' ) THEN
        WRITE(msgBuf,'(4A)') 'RNF_INIT_FIXED: RNF:',
     &    ' mitgcm_time_sampling = "', attTxt(1:attLen),
     &    '" cannot be read yet'
        CALL PRINT_ERROR( msgBuf, myThid )
        WRITE(msgBuf,'(2A)') 'RNF_INIT_FIXED: only "constant" is',
     &    ' implemented; time records come with RUNOFF-005'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
      ENDIF

      CALL RNF_NC_ATT_REAL( fid, NF_GLOBAL, 'the file',
     I     'mitgcm_grid_nx', 'RNF_INIT_FIXED',
     O     attVal, foundNx,
     I     myThid )
      fileNx = NINT( attVal )
      IF ( .NOT.foundNx ) THEN
        WRITE(msgBuf,'(2A)') 'RNF_INIT_FIXED: RNF: the file has no',
     &    ' global attribute mitgcm_grid_nx'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
      ENDIF
      CALL RNF_NC_ATT_REAL( fid, NF_GLOBAL, 'the file',
     I     'mitgcm_grid_ny', 'RNF_INIT_FIXED',
     O     attVal, foundNy,
     I     myThid )
      fileNy = NINT( attVal )
      IF ( .NOT.foundNy ) THEN
        WRITE(msgBuf,'(2A)') 'RNF_INIT_FIXED: RNF: the file has no',
     &    ' global attribute mitgcm_grid_ny'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
      ENDIF
      IF ( foundNx .AND. foundNy .AND.
     &     ( fileNx.NE.xSize .OR. fileNy.NE.ySize ) ) THEN
        WRITE(msgBuf,'(A,2I8,A,2I8)') 'RNF_INIT_FIXED: RNF: grid'
     &    //' size mismatch: the file has mitgcm_grid_nx,ny =',
     &    fileNx, fileNy, ' and the model', xSize, ySize
        CALL PRINT_ERROR( msgBuf, myThid )
        WRITE(msgBuf,'(2A)') 'RNF_INIT_FIXED: target_cell counts in',
     &    ' the global 2-D layout of a dense runoff file'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
      ENDIF

C--   2. Dimensions and required variables
      errNC = NF_INQ_DIMID( fid, 'time', dimId )
      CALL RNF_NC_ERROR( 'RNF_INIT_FIXED',
     &     'looking for the dimension time', errNC, myThid )
      errNC = NF_INQ_DIMLEN( fid, dimId, RNF_nRecFile )
      CALL RNF_NC_ERROR( 'RNF_INIT_FIXED',
     &     'reading the length of the dimension time',
     &     errNC, myThid )
      errNC = NF_INQ_DIMID( fid, 'source', dimId )
      CALL RNF_NC_ERROR( 'RNF_INIT_FIXED',
     &     'looking for the dimension source', errNC, myThid )
      errNC = NF_INQ_DIMLEN( fid, dimId, RNF_nSrcFile )
      CALL RNF_NC_ERROR( 'RNF_INIT_FIXED',
     &     'reading the length of the dimension source',
     &     errNC, myThid )
      errNC = NF_INQ_DIMID( fid, 'target', dimId )
      CALL RNF_NC_ERROR( 'RNF_INIT_FIXED',
     &     'looking for the dimension target', errNC, myThid )
      errNC = NF_INQ_DIMLEN( fid, dimId, RNF_nTgtFile )
      CALL RNF_NC_ERROR( 'RNF_INIT_FIXED',
     &     'reading the length of the dimension target',
     &     errNC, myThid )

      IF ( RNF_nRecFile .NE. 1 ) THEN
        WRITE(msgBuf,'(2A,I8,A)') 'RNF_INIT_FIXED: RNF: "constant"',
     &    ' sampling but the file has', RNF_nRecFile, ' record(s)'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
      ENDIF
      IF ( RNF_nSrcFile.LT.1 .OR. RNF_nTgtFile.LT.1 ) THEN
        WRITE(msgBuf,'(A,2(A,I8))') 'RNF_INIT_FIXED: RNF: the file',
     &    ' has', RNF_nSrcFile, ' source(s) and', RNF_nTgtFile,
     &    ' target entry(ies)'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
      ENDIF

C--   Nothing can be placed when the header is wrong
      IF ( errCount .GE. 1 ) GOTO 900

      errNC = NF_INQ_VARID( fid, 'source_id', idSrc )
      CALL RNF_NC_ERROR( 'RNF_INIT_FIXED',
     &     'looking for the required variable source_id',
     &     errNC, myThid )
      errNC = NF_INQ_VARID( fid, 'target_source', idSrc )
      CALL RNF_NC_ERROR( 'RNF_INIT_FIXED',
     &     'looking for the required variable target_source',
     &     errNC, myThid )
      errNC = NF_INQ_VARID( fid, 'target_cell', idCell )
      CALL RNF_NC_ERROR( 'RNF_INIT_FIXED',
     &     'looking for the required variable target_cell',
     &     errNC, myThid )
      errNC = NF_INQ_VARID( fid, 'target_fraction', idFrac )
      CALL RNF_NC_ERROR( 'RNF_INIT_FIXED',
     &     'looking for the required variable target_fraction',
     &     errNC, myThid )
      hasLev = .FALSE.
      errNC = NF_INQ_VARID( fid, 'target_level', idLev )
      IF ( errNC .EQ. NF_NOERR ) THEN
        hasLev = .TRUE.
      ELSEIF ( errNC .NE. NF_ENOTVAR ) THEN
        CALL RNF_NC_ERROR( 'RNF_INIT_FIXED',
     &       'looking for the variable target_level',
     &       errNC, myThid )
      ENDIF
      hasArea = .FALSE.
      errNC = NF_INQ_VARID( fid, 'target_cell_area', idArea )
      IF ( errNC .EQ. NF_NOERR ) THEN
        hasArea = .TRUE.
      ELSEIF ( errNC .NE. NF_ENOTVAR ) THEN
        CALL RNF_NC_ERROR( 'RNF_INIT_FIXED',
     &       'looking for the variable target_cell_area',
     &       errNC, myThid )
      ENDIF

C--   4a. The global record of every row of every tile, with the
C     arithmetic of S/R MDSIO_READ_FIELD. The grid kind is never
C     inferred from the geometry of the array.
      DO bj = 1, nSy
       DO bi = 1, nSx
        tBx = myXGlobalLo - 1 + (bi-1)*sNx
        tBy = myYGlobalLo - 1 + (bj-1)*sNy
        iGjLoc = 0
        jGjLoc = 1
#ifdef ALLOW_EXCH2
        IF ( useExch2ioLayOut ) THEN
          tN = W2_myTileList(bi,bj)
          tBx = exch2_txGlobalo(tN) - 1
          tBy = exch2_tyGlobalo(tN) - 1
          IF     ( exch2_mydNx(tN) .GT. xSize ) THEN
C-          face x-size larger than glob-size: fold it
            iGjLoc = 0
            jGjLoc = exch2_mydNx(tN) / xSize
          ELSEIF ( exch2_tNy(tN) .GT. ySize ) THEN
C-          tile y-size larger than glob-size: make a long line
            iGjLoc = exch2_mydNx(tN)
            jGjLoc = 0
          ELSE
C-          default (face fits into the global I/O array)
            iGjLoc = 0
            jGjLoc = 1
          ENDIF
        ENDIF
#endif /* ALLOW_EXCH2 */
        DO j = 1, sNy
          rowRec(j,bi,bj) = ( tBx + (j-1)*iGjLoc )/sNx
     &                    + ( tBy + (j-1)*jGjLoc )*global_nTx
        ENDDO
       ENDDO
      ENDDO

C--   3. and 4b. The target table, one chunk of entries at a time
      nChunk = 1 + (RNF_nTgtFile-1)/RNF_nBuf
      DO iChunk = 1, nChunk
        iEnt0 = (iChunk-1)*RNF_nBuf
        nLoc = MIN( RNF_nBuf, RNF_nTgtFile - iEnt0 )
        iStart(1) = iEnt0 + 1
        iCount(1) = nLoc
        errNC = NF_GET_VARA_INT( fid, idSrc, iStart, iCount,
     &                           tgtSrcBuf )
        CALL RNF_NC_ERROR( 'RNF_INIT_FIXED', 'reading target_source',
     &       errNC, myThid )
        errNC = NF_GET_VARA_INT( fid, idCell, iStart, iCount,
     &                           tgtCellBuf )
        CALL RNF_NC_ERROR( 'RNF_INIT_FIXED', 'reading target_cell',
     &       errNC, myThid )
        errNC = NF_GET_VARA_DOUBLE( fid, idFrac, iStart, iCount,
     &                              tgtFracBuf )
        CALL RNF_NC_ERROR( 'RNF_INIT_FIXED',
     &       'reading target_fraction', errNC, myThid )
        IF ( hasLev ) THEN
          errNC = NF_GET_VARA_INT( fid, idLev, iStart, iCount,
     &                             tgtLevBuf )
          CALL RNF_NC_ERROR( 'RNF_INIT_FIXED',
     &         'reading target_level', errNC, myThid )
        ENDIF
        IF ( hasArea ) THEN
          errNC = NF_GET_VARA_DOUBLE( fid, idArea, iStart, iCount,
     &                                tgtAreaBuf )
          CALL RNF_NC_ERROR( 'RNF_INIT_FIXED',
     &         'reading target_cell_area', errNC, myThid )
        ENDIF

        DO n = 1, nLoc
          kEnt = iEnt0 + n - 1
          iSrcF = tgtSrcBuf(n)
          cellF = tgtCellBuf(n)
          fracF = tgtFracBuf(n)
          levF = 1
          IF ( hasLev ) levF = tgtLevBuf(n)
          areaF = 0. _d 0
          IF ( hasArea ) areaF = tgtAreaBuf(n)
          badEnt = .FALSE.

C--       3. The entry itself, before any placement. A target_source
C         out of range has no source to name, so its message gives
C         the value and the entry.
          IF ( iSrcF.LT.0 .OR. iSrcF.GE.RNF_nSrcFile ) THEN
            badEnt = .TRUE.
            IF ( nMsgEnt .LT. RNF_maxErrMsg ) THEN
              nMsgEnt = nMsgEnt + 1
              WRITE(msgBuf,'(A,I12,A,I12)') 'RNF_INIT_FIXED: RNF:'
     &          //' target_source out of range:', iSrcF,
     &          ' at target entry (from 0)', kEnt
              CALL PRINT_ERROR( msgBuf, myThid )
            ENDIF
          ELSE
            IF ( cellF.LT.0 .OR. cellF.GE.xSize*ySize ) THEN
              badEnt = .TRUE.
              IF ( nMsgEnt .LT. RNF_maxErrMsg ) THEN
                nMsgEnt = nMsgEnt + 1
                CALL RNF_NC_SOURCE_ID( fid, iSrcF, 'RNF_INIT_FIXED',
     O               srcId, idLen,
     I               myThid )
                WRITE(msgBuf,'(3A,I12,A,I12)') 'RNF_INIT_FIXED:'
     &            //' RNF: target_cell out of range (valid: 0 to'
     &            //' nx*ny-1): source ', srcId(1:idLen),
     &            ', target_cell', cellF,
     &            ', target entry (from 0)', kEnt
                CALL PRINT_ERROR( msgBuf, myThid )
              ENDIF
            ENDIF
            IF ( levF .NE. 1 ) THEN
              badEnt = .TRUE.
              IF ( nMsgEnt .LT. RNF_maxErrMsg ) THEN
                nMsgEnt = nMsgEnt + 1
                CALL RNF_NC_SOURCE_ID( fid, iSrcF, 'RNF_INIT_FIXED',
     O               srcId, idLen,
     I               myThid )
                WRITE(msgBuf,'(3A,I12,A,I12)') 'RNF_INIT_FIXED:'
     &            //' RNF: target_level is not 1 (only the surface'
     &            //' cell is allowed): source ', srcId(1:idLen),
     &            ', target_level', levF,
     &            ', target entry (from 0)', kEnt
                CALL PRINT_ERROR( msgBuf, myThid )
              ENDIF
            ENDIF
            IF ( .NOT.( fracF.GE.0. _d 0 .AND.
     &                  fracF.LE.1. _d 0 ) ) THEN
C-            catches a negative fraction, which a sum of exactly 1
C             would hide, and one that is not a number
              badEnt = .TRUE.
              IF ( nMsgEnt .LT. RNF_maxErrMsg ) THEN
                nMsgEnt = nMsgEnt + 1
                CALL RNF_NC_SOURCE_ID( fid, iSrcF, 'RNF_INIT_FIXED',
     O               srcId, idLen,
     I               myThid )
                WRITE(msgBuf,'(3A,1PE16.8,A,I12)') 'RNF_INIT_FIXED:'
     &            //' RNF: target_fraction is not in [0,1]: source ',
     &            srcId(1:idLen), ', target_fraction', fracF,
     &            ', target entry (from 0)', kEnt
                CALL PRINT_ERROR( msgBuf, myThid )
              ENDIF
            ENDIF
          ENDIF
          IF ( badEnt ) THEN
            nErrEnt = nErrEnt + 1
            GOTO 100
          ENDIF

C--       4b. Place the entry: the tile and row whose global record
C         is the record of this cell
          iRec = cellF/sNx
          iLoc = cellF - iRec*sNx + 1
          placed = .FALSE.
          DO bj = 1, nSy
           DO bi = 1, nSx
            DO j = 1, sNy
             IF ( .NOT.placed .AND. rowRec(j,bi,bj).EQ.iRec ) THEN
              placed = .TRUE.
              nTgtOwn = nTgtOwn + 1
              nTgtNeed(bi,bj) = nTgtNeed(bi,bj) + 1

C--           the local source list of this tile
              kSrc = 0
              IF ( RNF_nSrc(bi,bj) .GE. 1 ) THEN
               IF ( RNF_srcGlob(RNF_nSrc(bi,bj),bi,bj) .EQ. iSrcF )
     &              kSrc = RNF_nSrc(bi,bj)
              ENDIF
              IF ( kSrc .EQ. 0 ) THEN
               DO k = 1, RNF_nSrc(bi,bj)
                IF ( RNF_srcGlob(k,bi,bj) .EQ. iSrcF ) kSrc = k
               ENDDO
              ENDIF
              IF ( kSrc .EQ. 0 ) THEN
               IF ( RNF_nSrc(bi,bj) .LT. RNF_nSrcTile ) THEN
                RNF_nSrc(bi,bj) = RNF_nSrc(bi,bj) + 1
                kSrc = RNF_nSrc(bi,bj)
                RNF_srcGlob(kSrc,bi,bj) = iSrcF
                nSrcNeed(bi,bj) = nSrcNeed(bi,bj) + 1
               ELSE
C-              no room: count what would be needed. Exact for a
C               file sorted by source, an upper bound otherwise.
                IF ( iSrcF .NE. lastOvf(bi,bj) )
     &               nSrcNeed(bi,bj) = nSrcNeed(bi,bj) + 1
                lastOvf(bi,bj) = iSrcF
               ENDIF
              ENDIF

              IF ( kSrc.GE.1 .AND.
     &             RNF_nTgt(bi,bj).LT.RNF_nTgtTile ) THEN
               RNF_nTgt(bi,bj) = RNF_nTgt(bi,bj) + 1
               k = RNF_nTgt(bi,bj)
               RNF_tgtI(k,bi,bj) = iLoc
               RNF_tgtJ(k,bi,bj) = j
               RNF_tgtSrc(k,bi,bj) = kSrc
               RNF_tgtFrac(k,bi,bj) = fracF
               RNF_srcFracSum(kSrc,bi,bj) =
     &              RNF_srcFracSum(kSrc,bi,bj) + fracF

C--            5. The checks of an owned target. The level is the
C              surface level of this column. The source id is read
C              from the file only when a message needs it.
               RNF_tgtK(k,bi,bj) = kSurf
#ifdef ALLOW_SHELFICE
               IF ( usingZCoords .AND. useShelfIce )
     &              RNF_tgtK(k,bi,bj) = kSurfC(iLoc,j,bi,bj)
#endif /* ALLOW_SHELFICE */
               refuse = 0
               IF ( RNF_tgtK(k,bi,bj).LT.1 .OR.
     &              RNF_tgtK(k,bi,bj).GT.Nr ) THEN
                refuse = 1
               ELSEIF ( maskC(iLoc,j,RNF_tgtK(k,bi,bj),bi,bj)
     &                  .EQ. 0. _d 0 ) THEN
                refuse = 1
               ELSEIF ( maskInC(iLoc,j,bi,bj) .EQ. 0. _d 0 ) THEN
                refuse = 2
#ifdef ALLOW_SHELFICE
               ELSEIF ( useShelfIce .AND.
     &                  kTopC(iLoc,j,bi,bj).NE.0 ) THEN
                refuse = 3
#endif /* ALLOW_SHELFICE */
               ELSEIF ( hasArea .AND. ABS( areaF
     &                  - rA(iLoc,j,bi,bj) ) .GT.
     &                  RNF_areaTol*rA(iLoc,j,bi,bj) ) THEN
                refuse = 4
               ENDIF
               IF ( refuse .NE. 0 ) THEN
                nErrTgt = nErrTgt + 1
                IF ( nMsgTgt .LT. RNF_maxErrMsg ) THEN
                 nMsgTgt = nMsgTgt + 1
                 CALL RNF_NC_SOURCE_ID( fid, iSrcF, 'RNF_INIT_FIXED',
     O                srcId, idLen,
     I                myThid )
                 IF     ( refuse .EQ. 1 ) THEN
                  CALL RNF_NC_TGT_ERROR( 'RNF_INIT_FIXED',
     I                 'target on land', srcId, idLen, cellF,
     I                 iLoc, j, bi, bj, myThid )
                 ELSEIF ( refuse .EQ. 2 ) THEN
                  CALL RNF_NC_TGT_ERROR( 'RNF_INIT_FIXED',
     I                 'target beyond an open boundary (maskInC = 0)',
     I                 srcId, idLen, cellF, iLoc, j, bi, bj, myThid )
                 ELSEIF ( refuse .EQ. 3 ) THEN
                  CALL RNF_NC_TGT_ERROR( 'RNF_INIT_FIXED',
     I                 'target under an ice shelf (kTopC /= 0)',
     I                 srcId, idLen, cellF, iLoc, j, bi, bj, myThid )
                 ELSE
                  CALL RNF_NC_TGT_ERROR( 'RNF_INIT_FIXED',
     I                 'target_cell_area differs from the cell area'
     I                 //' rA', srcId, idLen, cellF,
     I                 iLoc, j, bi, bj, myThid )
                 ENDIF
                ENDIF
               ENDIF
              ENDIF
             ENDIF
            ENDDO
           ENDDO
          ENDDO
  100     CONTINUE
        ENDDO
      ENDDO

C--   Report an entry count that exceeds an array bound of
C     RNF_SIZE.h. The counts are final here, so the message can say
C     what the bound would have to be.
      DO bj = 1, nSy
       DO bi = 1, nSx
        IF ( nSrcNeed(bi,bj) .GT. RNF_nSrcTile ) THEN
          nErrBnd = nErrBnd + 1
          WRITE(msgBuf,'(A,I8,A,I8)') 'RNF_INIT_FIXED: RNF:'
     &      //' RNF_nSrcTile =', RNF_nSrcTile,
     &      ' is too small: this tile needs', nSrcNeed(bi,bj)
          CALL PRINT_ERROR( msgBuf, myThid )
        ENDIF
        IF ( nTgtNeed(bi,bj) .GT. RNF_nTgtTile ) THEN
          nErrBnd = nErrBnd + 1
          WRITE(msgBuf,'(A,I8,A,I8)') 'RNF_INIT_FIXED: RNF:'
     &      //' RNF_nTgtTile =', RNF_nTgtTile,
     &      ' is too small: this tile needs', nTgtNeed(bi,bj)
          CALL PRINT_ERROR( msgBuf, myThid )
        ENDIF
        RNF_nSrcProc = RNF_nSrcProc + RNF_nSrc(bi,bj)
        RNF_nTgtProc = RNF_nTgtProc + RNF_nTgt(bi,bj)
       ENDDO
      ENDDO

      IF ( nErrEnt .GE. 1 ) THEN
        WRITE(msgBuf,'(A,I8,A)') 'RNF_INIT_FIXED: RNF:', nErrEnt,
     &    ' invalid entry(ies) in the target table'
        CALL PRINT_ERROR( msgBuf, myThid )
        WRITE(msgBuf,'(2A)') 'RNF_INIT_FIXED: the fraction sums are',
     &    ' not checked when the target table is invalid'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
      ELSE
C-      an entry that was dropped would come up as a missing
C       fraction, so the sums are only worth checking when every
C       entry of the table was read and accepted
        RNF_tableRead = .TRUE.
      ENDIF

#endif /* HAVE_NETCDF */

  900 CONTINUE
      _END_MASTER(myThid)
C--   Everyone waits for the master thread to have read the file
      _BARRIER

C--   6. The fraction sum of every source, over all tiles and
C     processes, one chunk of sources at a time. Every thread has to
C     take part in the global sum.
      nChunk = 0
      IF ( RNF_tableRead .AND. RNF_nSrcFile.GE.1 )
     &     nChunk = 1 + (RNF_nSrcFile-1)/RNF_nBuf
      DO iChunk = 1, nChunk
        iEnt0 = (iChunk-1)*RNF_nBuf
        nLoc = MIN( RNF_nBuf, RNF_nSrcFile - iEnt0 )
        _BEGIN_MASTER(myThid)
        DO n = 1, nLoc
         DO bj = 1, nSy
          DO bi = 1, nSx
            RNF_sumWrk(bi,bj,n) = 0. _d 0
          ENDDO
         ENDDO
        ENDDO
        DO bj = 1, nSy
         DO bi = 1, nSx
          DO k = 1, RNF_nSrc(bi,bj)
            n = RNF_srcGlob(k,bi,bj) - iEnt0 + 1
            IF ( n.GE.1 .AND. n.LE.nLoc )
     &           RNF_sumWrk(bi,bj,n) = RNF_sumWrk(bi,bj,n)
     &                               + RNF_srcFracSum(k,bi,bj)
          ENDDO
         ENDDO
        ENDDO
        _END_MASTER(myThid)
        _BARRIER
        CALL GLOBAL_SUM_VECTOR_RL( nLoc, RNF_sumWrk, RNF_sumVec,
     &                             myThid )
        _BEGIN_MASTER(myThid)
        DO n = 1, nLoc
          IF ( ABS( RNF_sumVec(n) - 1. _d 0 ) .GT. RNF_fracTol ) THEN
            nErrFrac = nErrFrac + 1
#ifdef HAVE_NETCDF
            IF ( nErrFrac .LE. RNF_maxErrMsg ) THEN
              CALL RNF_NC_SOURCE_ID( fid, iEnt0+n-1,
     I             'RNF_INIT_FIXED',
     O             srcId, idLen,
     I             myThid )
              WRITE(msgBuf,'(3A,1PE20.12)') 'RNF_INIT_FIXED: RNF:'
     &          //' fraction sum is not 1: source ', srcId(1:idLen),
     &          ', sum =', RNF_sumVec(n)
              CALL PRINT_ERROR( msgBuf, myThid )
            ENDIF
#endif /* HAVE_NETCDF */
          ENDIF
        ENDDO
        _END_MASTER(myThid)
        _BARRIER
      ENDDO

C--   Errors that only one process can see have to reach every
C     process before the stop
      CALL GLOBAL_SUM_INT( nErrTgt, myThid )
      CALL GLOBAL_SUM_INT( nErrBnd, myThid )
      CALL GLOBAL_SUM_INT( nTgtOwn, myThid )

      _BEGIN_MASTER(myThid)
      RNF_nTgtOwned = nTgtOwn

      IF ( nErrFrac .GE. 1 ) THEN
        WRITE(msgBuf,'(A,I8,A)') 'RNF_INIT_FIXED: RNF:', nErrFrac,
     &    ' source(s) whose fractions do not sum to 1'
        CALL PRINT_ERROR( msgBuf, myThid )
        WRITE(msgBuf,'(2A)') 'RNF_INIT_FIXED: a target that no tile',
     &    ' owns shows up here as a missing fraction'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
      ENDIF
      IF ( nErrTgt .GE. 1 ) THEN
        WRITE(msgBuf,'(A,I8,A)') 'RNF_INIT_FIXED: RNF:', nErrTgt,
     &    ' refused target(s) on all processes'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
      ENDIF
      IF ( nErrBnd .GE. 1 ) THEN
        WRITE(msgBuf,'(A,I8,A)') 'RNF_INIT_FIXED: RNF:', nErrBnd,
     &    ' array bound(s) of RNF_SIZE.h that are too small'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
      ENDIF

#ifdef HAVE_NETCDF
C--   The file is closed here and opened again by the routines that
C     read a record, so that RNF_INIT_VARIA can be called twice.
      IF ( fileOpen ) THEN
        errNC = NF_CLOSE( fid )
        CALL RNF_NC_ERROR( 'RNF_INIT_FIXED', 'closing the file',
     &       errNC, myThid )
      ENDIF
#endif /* HAVE_NETCDF */

      IF ( errCount .GE. 1 ) THEN
        WRITE(msgBuf,'(A,I8,A)') 'RNF_INIT_FIXED: detected',
     &    errCount, ' fatal error(s)'
        CALL PRINT_ERROR( msgBuf, myThid )
        CALL ALL_PROC_DIE( 0 )
        STOP 'ABNORMAL END: S/R RNF_INIT_FIXED'
      ENDIF

      WRITE(msgBuf,'(A,3(A,I8))') 'RNF_INIT_FIXED: read',
     &  ' sources:', RNF_nSrcFile, ', target entries:',
     &  RNF_nTgtFile, ', owned here:', RNF_nTgtProc
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )
      _END_MASTER(myThid)
      _BARRIER

C--   Print the package parameters and what the read found
      CALL RNF_SUMMARY( myThid )

#endif /* ALLOW_RNF */

      RETURN
      END
