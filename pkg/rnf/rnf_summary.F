#include "RNF_OPTIONS.h"

CBOP
C !ROUTINE: RNF_SUMMARY

C !INTERFACE: ==========================================================
      SUBROUTINE RNF_SUMMARY( myThid )

C !DESCRIPTION:
C     Print the RNF parameters (data.rnf), the constants fixed by the
C     model contract, the array bounds of RNF_SIZE.h and what the
C     static read of RNF_file found, to the standard output of every
C     process, in the layout of EXF_SUMMARY.
C     Called at the end of RNF_INIT_FIXED, so the counts are the
C     final ones. RNF_nSrcProc and RNF_nTgtProc are the counts of
C     this process; RNF_nTgtOwned is summed over all processes.

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "RNF_SIZE.h"
#include "RNF.h"

C !INPUT PARAMETERS: ===================================================
C  myThid               :: my thread Id number
      INTEGER myThid

C !OUTPUT PARAMETERS: ==================================================
C  none

#ifdef ALLOW_RNF

C !LOCAL VARIABLES: ====================================================
C  blkLin               :: blank line
C  msgBuf               :: message buffer
      CHARACTER*1 blkLin
      CHARACTER*(MAX_LEN_MBUF) msgBuf
CEOP

      _BEGIN_MASTER(myThid)

      blkLin = ' '
      CALL PRINT_MESSAGE( blkLin, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)')
     &'// ======================================================='
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)')
     &'// Sparse runoff (RNF) configuration >>> START <<<'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)')
     &'// ======================================================='
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )
      CALL PRINT_MESSAGE( blkLin, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )

C--   Parameters of data.rnf
      CALL WRITE_0D_C( RNF_file, -1, INDEX_NONE, 'RNF_file =',
     &  ' /* sparse runoff NetCDF file */')
      CALL WRITE_0D_L( RNF_holdRecord, INDEX_NONE,
     &  'RNF_holdRecord =',
     &  ' /* hold each record (F: interpolate in time) */')
      CALL WRITE_0D_L( RNF_useYearlyFiles, INDEX_NONE,
     &  'RNF_useYearlyFiles =',
     &  ' /* append _YYYY to RNF_file by model year */')
      CALL WRITE_0D_I( RNF_startDate1, INDEX_NONE,
     &  'RNF_startDate1 =',
     &  ' /* start date override YYYYMMDD (0: not set) */')
      CALL WRITE_0D_I( RNF_startDate2, INDEX_NONE,
     &  'RNF_startDate2 =',
     &  ' /* start date override HHMMSS (0: not set) */')
      CALL WRITE_0D_RL( RNF_startTime, INDEX_NONE,
     &  'RNF_startTime =',
     &  ' /* start time override [s] (UNSET: from file) */')
      CALL WRITE_0D_RL( RNF_period, INDEX_NONE,
     &  'RNF_period =',
     &  ' /* record period override (UNSET: from file) */')
      CALL WRITE_0D_RL( RNF_repCycle, INDEX_NONE,
     &  'RNF_repCycle =',
     &  ' /* repeat cycle override [s] (UNSET: from file) */')
      CALL WRITE_0D_L( RNF_useTemp, INDEX_NONE,
     &  'RNF_useTemp =',
     &  ' /* apply runoff_temperature if present */')
      CALL WRITE_0D_L( RNF_useSalt, INDEX_NONE,
     &  'RNF_useSalt =',
     &  ' /* apply runoff_salinity if present */')
      CALL WRITE_0D_L( RNF_usePtracers, INDEX_NONE,
     &  'RNF_usePtracers =',
     &  ' /* apply runoff_ptracer_* if present */')
      CALL WRITE_0D_RL( RNF_monFreq, INDEX_NONE,
     &  'RNF_monFreq =',
     &  ' /* RNF monitor frequency [s] */')
      CALL WRITE_0D_I( RNF_debugLev, INDEX_NONE,
     &  'RNF_debugLev =',
     &  ' /* select RNF-debug printing level */')

C--   Constants fixed by the model contract
      CALL WRITE_0D_RL( RNF_fracTol, INDEX_NONE,
     &  'RNF_fracTol =',
     &  ' /* tolerance on the fraction sum of a source */')
      CALL WRITE_0D_RL( RNF_areaTol, INDEX_NONE,
     &  'RNF_areaTol =',
     &  ' /* relative tolerance on the target cell area */')

C--   Array bounds of RNF_SIZE.h
      CALL WRITE_0D_I( RNF_nSrcTile, INDEX_NONE,
     &  'RNF_nSrcTile =',
     &  ' /* max number of sources on one tile */')
      CALL WRITE_0D_I( RNF_nTgtTile, INDEX_NONE,
     &  'RNF_nTgtTile =',
     &  ' /* max number of target entries on one tile */')
      CALL WRITE_0D_I( RNF_nBuf, INDEX_NONE,
     &  'RNF_nBuf =',
     &  ' /* chunk length for NetCDF reads */')
      CALL WRITE_0D_I( RNF_nTr, INDEX_NONE,
     &  'RNF_nTr =',
     &  ' /* max number of runoff tracers */')
      CALL WRITE_0D_I( RNF_nFile, INDEX_NONE,
     &  'RNF_nFile =',
     &  ' /* number of runoff files */')

C--   What the static read of RNF_file found
      CALL WRITE_0D_I( RNF_nSrcFile, INDEX_NONE,
     &  'RNF_nSrcFile =',
     &  ' /* sources in the file */')
      CALL WRITE_0D_I( RNF_nTgtFile, INDEX_NONE,
     &  'RNF_nTgtFile =',
     &  ' /* target entries in the file */')
      CALL WRITE_0D_I( RNF_nRecFile, INDEX_NONE,
     &  'RNF_nRecFile =',
     &  ' /* time records in the file */')
      CALL WRITE_0D_I( RNF_nSrcProc, INDEX_NONE,
     &  'RNF_nSrcProc =',
     &  ' /* tile-source entries of this process */')
      CALL WRITE_0D_I( RNF_nTgtProc, INDEX_NONE,
     &  'RNF_nTgtProc =',
     &  ' /* target entries owned by this process */')
      CALL WRITE_0D_I( RNF_nTgtOwned, INDEX_NONE,
     &  'RNF_nTgtOwned =',
     &  ' /* target entries owned over all processes */')

      CALL PRINT_MESSAGE( blkLin, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)')
     &'// ======================================================='
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)')
     &'// Sparse runoff (RNF) configuration  >>> END <<<'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)')
     &'// ======================================================='
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )
      CALL PRINT_MESSAGE( blkLin, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )

      _END_MASTER(myThid)

#endif /* ALLOW_RNF */

      RETURN
      END
