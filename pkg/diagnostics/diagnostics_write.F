#include "DIAG_OPTIONS.h"

      SUBROUTINE DIAGNOSTICS_WRITE(
     I                              modelEnd,
     I                              myTime, myIter, myThid )
C***********************************************************************
C  Purpose
C  -------
C    Output sequence for the (multiple) diagnostics output files
C
C  Arguments  Description
C  ----------------------
C     modelEnd :: true if call at end of model run.
C     myTime   :: Current time of simulation ( s )
C     myIter   :: Current Iteration Number
C     myThid   :: my Thread Id number
C***********************************************************************
       IMPLICIT NONE
#include "EEPARAMS.h"
#include "SIZE.h"
#include "PARAMS.h"
#include "DIAGNOSTICS_SIZE.h"
#include "DIAGNOSTICS_P2SHARE.h"
#include "DIAGNOSTICS.h"

C     !INPUT PARAMETERS:
      LOGICAL modelEnd
      _RL     myTime
      INTEGER myIter, myThid

C     !FUNCTIONS:
      LOGICAL  DIFF_PHASE_MULTIPLE
      EXTERNAL DIFF_PHASE_MULTIPLE
#ifdef ALLOW_FIZHI
      LOGICAL  ALARM2
      EXTERNAL ALARM2
#endif

c Local variables
c ===============
      INTEGER   n, nd
      INTEGER   myItM1, wrIter
      LOGICAL   dump2fileNow, write2file
      LOGICAL   writeDiags(numLists), writeStats(numLists)
      _RL       phiSec, freqSec, wrTime
#ifdef ALLOW_FIZHI
      CHARACTER *9 tagname
#endif

      IF ( myIter.NE.nIter0 ) THEN
        myItM1 = myIter - 1

C***********************************************************************
C***   Check to see if its time for Diagnostic Output                ***
C***********************************************************************

        write2file = .FALSE.
        DO n = 1,nlists
          freqSec = freq(n)
          phiSec = phase(n)

          IF ( freqSec.LT.0. ) THEN
C--     write snap-shot with suffix = myIter to be consistent with
C       time-average diagnostics (e.g., freq=-1 & freq=1):
c           wrIter = myIter
c           wrTime = myTime
C--     write snap-shot with suffix = myIter-1 to be consistent with
C       state-variable time-step:
            wrIter = myItM1
            wrTime = myTime - deltaTClock
          ELSE
            wrIter = myIter
            wrTime = myTime
          ENDIF
          dump2fileNow = DIFF_PHASE_MULTIPLE( phiSec, freqSec,
     &                                        wrTime, deltaTClock )
#ifdef ALLOW_FIZHI
          IF ( useFIZHI ) THEN
            WRITE(tagname,'(A,I2.2)')'diagtag',n
            dump2fileNow = ALARM2(tagname)
          ENDIF
#endif
#ifdef ALLOW_CAL
          IF ( useCAL ) THEN
            CALL CAL_TIME2DUMP( phiSec, freqSec, deltaTClock,
     U                          dump2fileNow,
     I                          wrTime, myIter, myThid )
          ENDIF
#endif /* ALLOW_CAL */
          IF ( dumpAtLast .AND. modelEnd
     &                    .AND. freqSec.GE.0. ) dump2fileNow = .TRUE.
          IF ( useDiag4AdjOutp ) THEN
            nd = ABS(jdiag(1,n))
            IF ( gdiag(nd)(4:4).EQ.'A' ) dump2fileNow = .FALSE.
          ENDIF
          IF ( dump2fileNow ) THEN
            write2file = .TRUE.
            CALL DIAGNOSTICS_OUT(n,wrTime,wrIter,myThid)
          ENDIF
          writeDiags(n) = dump2fileNow
C-      end loop on list id number n
        ENDDO

C---   Check to see if its time for Statistics Diag. Output

        DO n = 1,diagSt_nbLists
          freqSec = diagSt_freq(n)
          phiSec = diagSt_phase(n)

          IF ( freqSec.LT.0. ) THEN
C--     write snap-shot with suffix = myIter to be consistent with
C       time-average diagnostics (e.g., freq=-1 & freq=1):
c           wrIter = myIter
c           wrTime = myTime
C--     write snap-shot with suffix = myIter-1 to be consistent with
C       state-variable time-step:
            wrIter = myItM1
            wrTime = myTime - deltaTClock
          ELSE
            wrIter = myIter
            wrTime = myTime
          ENDIF
          dump2fileNow = DIFF_PHASE_MULTIPLE( phiSec, freqSec,
     &                                        wrTime, deltaTClock )
#ifdef ALLOW_FIZHI
          IF ( useFIZHI ) THEN
           WRITE(tagname,'(A,I2.2)')'diagStg',n
           dump2fileNow = ALARM2(tagname)
          ENDIF
#endif
#ifdef ALLOW_CAL
          IF ( useCAL ) THEN
            CALL CAL_TIME2DUMP( phiSec, freqSec, deltaTClock,
     U                          dump2fileNow,
     I                          wrTime, myIter, myThid )
          ENDIF
#endif /* ALLOW_CAL */
          IF ( dumpAtLast .AND. modelEnd
     &                    .AND. freqSec.GE.0. ) dump2fileNow = .TRUE.
          IF ( dump2fileNow ) THEN
            write2file = .TRUE.
            CALL DIAGSTATS_OUTPUT(n,wrTime,wrIter,myThid)
          ENDIF
          writeStats(n) = dump2fileNow
C-      end loop on list id number n
        ENDDO

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

        IF ( write2file ) THEN
          IF ( diag_dBugLevel.GE.debLevC ) THEN
            CALL DIAGNOSTICS_SUMMARY( 1, myTime, myIter, myThid )
          ENDIF
C-      wait for everyone before setting arrays to zero:
          _BARRIER
        ENDIF
        IF ( modelEnd ) THEN
C-      Track diagnostics pkg activation status:
c         IF ( diag_pkgStatus.NE.ready2fillDiags ) STOP
          _BARRIER
          _BEGIN_MASTER(myThid)
          diag_pkgStatus = 99
          _END_MASTER(myThid)
          _BARRIER
C       Close all Stat-diags output files
          CALL DIAGSTATS_CLOSE_IO( myThid )
        ENDIF

C--     Clear storage space:
        DO n = 1,nlists
          IF ( writeDiags(n) ) CALL DIAGNOSTICS_CLEAR(n,myThid)
        ENDDO
        DO n = 1,diagSt_nbLists
          IF ( writeStats(n) ) CALL DIAGSTATS_CLEAR( n, myThid )
        ENDDO

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      ENDIF

      RETURN
      END
