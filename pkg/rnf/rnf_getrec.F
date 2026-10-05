#include "RNF_OPTIONS.h"

C--  File rnf_getrec.F: record selection of the sparse runoff package
C--   Contents
C--   o RNF_GETREC    :: the two records that bracket the model time,
C--                      and the weight of the earlier one
C--   o RNF_FILE_NAME :: name of the file of one model year
C--
C--  Both are thin: RNF_GETREC delegates every mode pkg/exf has to the
C--  pkg/exf routine that implements it, so a sparse run selects
C--  records with the same code a dense runoffFile run does and the
C--  dense oracle stays an oracle (package design, decision 7).

CBOP
C !ROUTINE: RNF_GETREC

C !INTERFACE: ==========================================================
      SUBROUTINE RNF_GETREC(
     I                   myTime, myIter,
     O                   fac, first, changed,
     O                   count0, count1, year0, year1,
     I                   myThid )

C !DESCRIPTION:
C     Return the two time records of the sparse runoff file that
C     bracket myTime and the weight fac of the earlier one, so that
C
C       flux(myTime) = fac*flux(count0) + (1-fac)*flux(count1)
C
C     with count0 from the file of year0 and count1 from that of
C     year1 (both 0 without RNF_useYearlyFiles).
C
C     The selection is pkg/exf's own, with the timing RNF_TIME_SETUP
C     resolved from the file (RNF_recPeriod, RNF_recStart,
C     RNF_recCycle of RNF.h):
C
C     | RNF_recPeriod | routine                                     |
C     | 0 (constant)  | EXF_GetFFieldRec, which returns record 1    |
C     |               | with weight 1 (exf_getffieldrec.F)          |
C     | > 0 seconds   | EXF_GetFFieldRec: a fixed period, with the  |
C     |               | repeat cycle and the yearly-file wrap       |
C     | -12           | cal_GetMonthsRec: 12 records, January to    |
C     |               | December, repeated every model year         |
C     | -1            | EXF_GetMonthsRec: consecutive calendar      |
C     |               | months from the month of RNF_recStart       |
C
C     This is why the repeat cycle of a fixed-period climatology
C     follows the file's real dates: EXF_GetFFieldRec wraps
C     myTime - RNF_recStart modulo RNF_recCycle, and RNF_recStart is
C     the model time of the real date record 1 carries. A reader that
C     mapped a record to a day of the model's own year instead would
C     agree only until a leap day intervenes.
C
C     **Hold-exact** (RNF_holdRecord) is this package's own: it
C     returns the same bracketing pair but a weight of exactly 1 or
C     exactly 0, so the caller applies one record of the file and no
C     arithmetic at all. The record held is the one whose interval
C     contains myTime:
C
C     - a fixed period: count0, whose interval runs from its own time
C       to the next record's time. For the yearly files of
C       exf (and of the schema, section 14) the records sit at the
C       start of their bounds - 1 January 00:00 for a daily file - so
C       holding count0 holds each record over its own bounds;
C     - monthly (-12 and -1): the record of the calendar month
C       containing myTime, which is NOT the nearer of the two
C       mid-month times the interpolation brackets with. In the first
C       half of a month the current month is count1, in the second
C       half it is count0, and the two halves are not equal because
C       the bracketing midpoints belong to months of different
C       lengths;
C     - constant: record 1, as always.
C
C     first and changed are passed through from the exf routine.
C     Nothing depends on them: the caller tags each record buffer
C     with what it holds, so a buffer is re-read exactly when the
C     bracket moves off it. They are reported in the debug trace.

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
C  fac                  :: weight of record count0
C  first                :: first step of the run (both records new)
C  changed              :: the bracket moved since the previous step
C  count0, count1       :: the two records, counted from 1
C  year0, year1         :: the file year each comes from (0: one file)
      _RL     fac
      LOGICAL first, changed
      INTEGER count0, count1, year0, year1

#ifdef ALLOW_RNF

C !LOCAL VARIABLES: ====================================================
C  mon0, mon1           :: calendar month of the two brackets
C  mYr0, mYr1           :: calendar year of the two brackets
C  nowDate              :: model date of myTime
      CHARACTER*(MAX_LEN_MBUF) msgBuf
#ifdef ALLOW_CAL
      INTEGER mon0, mon1, mYr0, mYr1
      INTEGER nowDate(4)
      INTEGER yy, mm, dd, ss, lp, wd
      _RL     facLoc
      LOGICAL firstLoc, chgLoc
#endif /* ALLOW_CAL */
CEOP

      fac     = 1. _d 0
      first   = myIter .EQ. nIter0
      changed = .FALSE.
      count0  = 1
      count1  = 1
      year0   = 0
      year1   = 0
#ifdef ALLOW_CAL
      mon0 = 0
      mon1 = 0
      mYr0 = 0
      mYr1 = 0
#endif /* ALLOW_CAL */

      IF ( RNF_recPeriod .EQ. -12. _d 0 ) THEN
#ifdef ALLOW_CAL
C--     12 monthly records, January to December, repeated every model
C       calendar year. The record number is the calendar month, so
C       the year cal_GetMonthsRec returns names no file: pkg/exf
C       ignores it too (exf_getyearlyfieldname.F returns the base
C       name for period -12).
        CALL CAL_GETMONTHSREC(
     O       fac, first, changed,
     O       count0, count1, year0, year1,
     I       myTime, myIter, myThid )
        mon0  = count0
        mon1  = count1
        mYr0  = year0
        mYr1  = year1
        year0 = 0
        year1 = 0
#endif /* ALLOW_CAL */
      ELSEIF ( RNF_recPeriod .EQ. -1. _d 0 ) THEN
#ifdef ALLOW_CAL
C--     Consecutive calendar months. EXF_GetMonthsRec turns the
C       month/year pair of cal_GetMonthsRec into a record number
C       counted from the month RNF_recStart falls in; the pair itself
C       is what hold-exact needs, so it is taken from the same
C       routine EXF_GetMonthsRec calls.
        CALL CAL_GETMONTHSREC(
     O       facLoc, firstLoc, chgLoc,
     O       mon0, mon1, mYr0, mYr1,
     I       myTime, myIter, myThid )
        CALL EXF_GETMONTHSREC(
     I       RNF_recStart, RNF_useYearlyFiles,
     O       fac, first, changed,
     O       count0, count1, year0, year1,
     I       myTime, myIter, myThid )
#endif /* ALLOW_CAL */
      ELSE
C--     Constant (period 0) and any fixed period, including the
C       repeat cycle and the yearly-file wrap
        CALL EXF_GETFFIELDREC(
     I       RNF_recStart, RNF_recPeriod, RNF_recCycle,
     I       'rnf_flux', RNF_useYearlyFiles,
     O       fac, first, changed,
     O       count0, count1, year0, year1,
     I       myTime, myIter, myThid )
      ENDIF

C--   Hold-exact: one record of the file, with no arithmetic
      IF ( RNF_holdRecord .AND. RNF_recPeriod.NE.0. _d 0 ) THEN
        IF ( RNF_recPeriod .GT. 0. _d 0 ) THEN
C-        the record at or before myTime, held to the next one
          fac = 1. _d 0
        ELSE
#ifdef ALLOW_CAL
C-        the record of the calendar month holding myTime
          CALL CAL_GETDATE( myIter, myTime, nowDate, myThid )
          CALL CAL_CONVDATE( nowDate, yy, mm, dd, ss, lp, wd, myThid )
          IF     ( mon0.EQ.mm .AND. mYr0.EQ.yy ) THEN
            fac = 1. _d 0
          ELSEIF ( mon1.EQ.mm .AND. mYr1.EQ.yy ) THEN
            fac = 0. _d 0
          ELSE
C-          Cannot happen: the two brackets of cal_GetMonthsRec are
C           the current month and one of its neighbours. Refuse
C           rather than apply a record that was not asked for.
            WRITE(msgBuf,'(2A,I9,I3,A,4I6)') 'RNF_GETREC: RNF:',
     &        ' hold-exact: the model month', yy, mm,
     &        ' is neither bracket', mYr0, mon0, mYr1, mon1
            CALL PRINT_ERROR( msgBuf, myThid )
            CALL ALL_PROC_DIE( 0 )
            STOP 'ABNORMAL END: S/R RNF_GETREC'
          ENDIF
#endif /* ALLOW_CAL */
        ENDIF
      ENDIF

#endif /* ALLOW_RNF */

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: RNF_FILE_NAME

C !INTERFACE: ==========================================================
      SUBROUTINE RNF_FILE_NAME(
     I                   year,
     O                   fileName,
     I                   myThid )

C !DESCRIPTION:
C     Name of the sparse runoff file that holds the records of model
C     year "year". Without RNF_useYearlyFiles it is RNF_file itself
C     and "year" is ignored.
C
C     With RNF_useYearlyFiles, RNF_file is the base name of a set and
C     the file of a year is <base>_YYYY.nc, which is what the schema
C     calls a yearly file (section 7) and what the converter writes.
C     A trailing ".nc" of RNF_file is replaced rather than kept, so
C     RNF_file = 'runoff.nc' and RNF_file = 'runoff' both select
C     runoff_1979.nc. This differs from exf, which appends "_YYYY" to
C     a dense file name that has no extension
C     (pkg/exf/exf_getyearlyfieldname.F).

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "RNF_SIZE.h"
#include "RNF.h"

C !INPUT PARAMETERS: ===================================================
C  year                 :: model year of the records wanted
C  myThid               :: my thread Id number
      INTEGER year
      INTEGER myThid

C !OUTPUT PARAMETERS: ==================================================
C  fileName             :: the name to open
      CHARACTER*(MAX_LEN_FNAM) fileName

#ifdef ALLOW_RNF

C !FUNCTIONS: ==========================================================
      INTEGER  ILNBLNK
      EXTERNAL ILNBLNK

C !LOCAL VARIABLES: ====================================================
      INTEGER iL
      CHARACTER*(MAX_LEN_MBUF) msgBuf
CEOP

      fileName = RNF_file
      IF ( .NOT.RNF_useYearlyFiles ) RETURN

      iL = MAX( 1, ILNBLNK(RNF_file) )
      IF ( iL.GT.3 .AND. RNF_file(iL-2:iL).EQ.'.nc' ) iL = iL - 3
      IF ( iL+8 .GT. MAX_LEN_FNAM ) THEN
        WRITE(msgBuf,'(2A)') 'RNF_FILE_NAME: RNF: RNF_file is too',
     &    ' long to append _YYYY.nc to it'
        CALL PRINT_ERROR( msgBuf, myThid )
        CALL ALL_PROC_DIE( 0 )
        STOP 'ABNORMAL END: S/R RNF_FILE_NAME'
      ENDIF
      fileName = ' '
      WRITE(fileName(1:iL+8),'(2A,I4.4,A)') RNF_file(1:iL), '_',
     &     year, '.nc'

#endif /* ALLOW_RNF */

      RETURN
      END
