#include "RNF_OPTIONS.h"

CBOP
C !ROUTINE: RNF_TIME_SETUP

C !INTERFACE: ==========================================================
      SUBROUTINE RNF_TIME_SETUP( fid, errCount, myThid )

C !DESCRIPTION:
C     Resolve the time handling of the sparse runoff file into the
C     three values pkg/exf selects records with: a period in exf's
C     convention, the time of record 1 and a repeat cycle
C     (RNF_recPeriod, RNF_recStart, RNF_recCycle of RNF.h).
C     RNF_GETREC then hands them to the exf routines themselves, so
C     the sparse path picks records with the code the dense path
C     runs (package design, decision 7).
C
C     Called by the master thread of RNF_INIT_FIXED with the file
C     open, after the dimensions have been read (RNF_nRecFile is
C     used). Every problem raises errCount; the caller stops the run.
C
C     What the file states (runoff schema 1.0, section 7) and how it
C     maps:
C
C     | mitgcm_time_sampling | mitgcm_time_repeat | RNF_recPeriod |
C     | "constant"           | -                  | 0             |
C     | "fixed"              | "none"             | the period    |
C     | "fixed"              | "annual"           | the period,   |
C     |                      |                    | cycle = span  |
C     |                      |                    | of time_bnds  |
C     | "monthly"            | "none"             | -1            |
C     | "monthly"            | "annual"           | -12           |
C
C     Anything else stops the run. "yearly" sampling (one record per
C     calendar year) has no exf mode and schema 1.0 puts its time at
C     the midpoint of the year, so mapping it needs its own code; it
C     is refused rather than approximated (RUNOFF-001 review 1).
C
C     **The repeat cycle of a fixed-period climatology is anchored on
C     the file's real dates, not on a nominal calendar year.** The
C     cycle is the span of time_bnds in seconds, and record 1 sits at
C     the real date its own time value names, which is what
C     EXF_GETFFIELD_START converts to a model time through pkg/cal.
C     EXF_GetFFieldRec then wraps the elapsed time modulo that cycle,
C     exactly as it does for a dense runoffFile with
C     runoffRepCycle = the same number. Anchoring on a nominal year
C     instead - mapping a record to a day of the model's year - agrees
C     with this only while no leap day has intervened: a 365-day cycle
C     started in 1978 and the day-of-year mapping part company on
C     1980-02-29 and stay one day apart, which on the lab_sea
C     climatology is up to 2.2% of the peak runoff
C     (tests/rnf/clim_anchor_check.py measures it).
C
C     Three things the time axis is NOT read for:
C     1. "constant" sampling. The reader returns before touching
C        time:units, because a constant file carries no meaningful
C        reference date: the converter writes 0001-01-01, which
C        precedes pkg/cal's own reference date (15 October 1582,
C        pkg/cal/cal_set.F) and must not be handed to it;
C     2. a run without pkg/cal. The file's time is then model time in
C        seconds counted from the reference date of its units (schema
C        section 14), no date is converted, and the file's calendar is
C        not compared with the model's, because the model has none.
C        The cs32 runoff is that case: its build has -cal and its file
C        is written on the 360_day calendar;
C     3. a start date or start time given in data.rnf. data.rnf
C        overrides the file (model contract), so an override is used
C        as given and the file's own date is not needed - which is
C        also the way out when a file's reference date precedes
C        pkg/cal's.

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "RNF_SIZE.h"
#include "RNF.h"
#ifdef ALLOW_CAL
# include "cal.h"
#endif
#ifdef HAVE_NETCDF
# include "netcdf.inc"
#endif

C !INPUT PARAMETERS: ===================================================
C  fid                  :: NetCDF id of the open sparse runoff file
C  myThid               :: my thread Id number
      INTEGER fid
      INTEGER myThid

C !INPUT/OUTPUT PARAMETERS: ============================================
C  errCount             :: fatal errors, raised by this routine
      INTEGER errCount

#ifdef ALLOW_RNF

C !LOCAL VARIABLES: ====================================================
C  sampTxt, repTxt      :: mitgcm_time_sampling, mitgcm_time_repeat
C  annual               :: the records are a climatology
C  filePer              :: period the file states [s]
C  fileCyc              :: repeat cycle the file states [s]
C  recSecs              :: time of record 1 in the file's units [s
C                          after the reference date]
C  unitSec              :: seconds per unit of time:units
C  refDat1, refDat2     :: reference date of time:units
C  dat1, dat2           :: date of record 1 (YYYYMMDD, HHMMSS)
C  nNum, num            :: integers scanned out of the reference date
      CHARACTER*(80) sampTxt, repTxt, unitTxt
      INTEGER sampLen, repLen
      LOGICAL found, annual, haveDate
      _RL     filePer, fileCyc, recSecs, unitSec
      INTEGER refDat1, refDat2, dat1, dat2
      INTEGER num(7), nNum, iS, iL, i, n
      LOGICAL inNum
      CHARACTER*1 chr
      CHARACTER*(MAX_LEN_MBUF) msgBuf
#ifdef ALLOW_CAL
      INTEGER refArr(4), intv(4), datArr(4)
      INTEGER yy, mm, dd, ss, lp, wd
      LOGICAL calOk
#endif /* ALLOW_CAL */
#ifdef HAVE_NETCDF
      CHARACTER*(256) unitAll
      INTEGER errNC, varId
      INTEGER iStart(2), iCount(2)
      _RL     tmpVal(2)
      _RL     bnd0, bnd1
      LOGICAL hasBnds
#endif /* HAVE_NETCDF */
CEOP

C--   Nothing resolved yet: a constant single record, which is what a
C     file without any timing attribute would be.
      RNF_recPeriod = 0. _d 0
      RNF_recStart  = 0. _d 0
      RNF_recCycle  = 0. _d 0
      RNF_recDate1  = 0
      RNF_recDate2  = 0
      filePer  = 0. _d 0
      fileCyc  = 0. _d 0
      recSecs  = 0. _d 0
      haveDate = .FALSE.
      annual   = .FALSE.
      dat1     = 0
      dat2     = 0
      refDat1  = 0
      refDat2  = 0
      unitSec  = 0. _d 0

#ifdef HAVE_NETCDF

C--   1. The sampling kind, and the repeat flag that completes it
      CALL RNF_NC_ATT_TEXT( fid, NF_GLOBAL, 'mitgcm_time_sampling',
     I     'RNF_TIME_SETUP',
     O     sampTxt, sampLen, found,
     I     myThid )
      IF ( .NOT.found ) THEN
        WRITE(msgBuf,'(2A)') 'RNF_TIME_SETUP: RNF: the file has no',
     &    ' global attribute mitgcm_time_sampling'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
        RETURN
      ENDIF

      CALL RNF_NC_ATT_TEXT( fid, NF_GLOBAL, 'mitgcm_time_repeat',
     I     'RNF_TIME_SETUP',
     O     repTxt, repLen, found,
     I     myThid )
      IF ( .NOT.found ) THEN
C-      optional; "none" is the default (schema 1.0, section 7)
        repTxt = 'none'
        repLen = 4
      ENDIF
      IF     ( repTxt(1:repLen) .EQ. 'none' ) THEN
        annual = .FALSE.
      ELSEIF ( repTxt(1:repLen) .EQ. 'annual' ) THEN
        annual = .TRUE.
      ELSE
        WRITE(msgBuf,'(4A)') 'RNF_TIME_SETUP: RNF:',
     &    ' mitgcm_time_repeat = "', repTxt(1:repLen),
     &    '" is not "none" or "annual"'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
        RETURN
      ENDIF

      IF ( sampTxt(1:sampLen) .EQ. 'constant' ) THEN
C--     One record, read once. The time axis is not read at all: see
C       the description (a constant file's reference date precedes
C       the pkg/cal reference date).
        filePer = 0. _d 0
        IF ( RNF_nRecFile .NE. 1 ) THEN
          WRITE(msgBuf,'(2A,I8,A)') 'RNF_TIME_SETUP: RNF: "constant"',
     &      ' sampling but the file has', RNF_nRecFile, ' record(s)'
          CALL PRINT_ERROR( msgBuf, myThid )
          errCount = errCount + 1
        ENDIF
      ELSEIF ( sampTxt(1:sampLen) .EQ. 'fixed' ) THEN
        CALL RNF_NC_ATT_REAL( fid, NF_GLOBAL, 'the file',
     I       'mitgcm_time_period', 'RNF_TIME_SETUP',
     O       filePer, found,
     I       myThid )
        IF ( .NOT.found .OR. filePer.LE.0. _d 0 ) THEN
          WRITE(msgBuf,'(2A)') 'RNF_TIME_SETUP: RNF: "fixed"',
     &      ' sampling needs mitgcm_time_period > 0 [s]'
          CALL PRINT_ERROR( msgBuf, myThid )
          errCount = errCount + 1
          RETURN
        ENDIF
      ELSEIF ( sampTxt(1:sampLen) .EQ. 'monthly' ) THEN
        filePer = -1. _d 0
        IF ( annual ) filePer = -12. _d 0
      ELSE
        WRITE(msgBuf,'(4A)') 'RNF_TIME_SETUP: RNF:',
     &    ' mitgcm_time_sampling = "', sampTxt(1:sampLen),
     &    '" is not supported'
        CALL PRINT_ERROR( msgBuf, myThid )
        WRITE(msgBuf,'(2A)') 'RNF_TIME_SETUP: use "constant",',
     &    ' "fixed" or "monthly"; "yearly" has no pkg/exf mode'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
        RETURN
      ENDIF

C--   2. The repeat cycle of a fixed-period climatology: the span of
C     time_bnds, i.e. the first start to the last end. This is the
C     number a dense run would give as runoffRepCycle, and it is
C     anchored on the file's real dates (see the description).
      IF ( annual .AND. filePer.GT.0. _d 0 ) THEN
        errNC = NF_INQ_VARID( fid, 'time_bnds', varId )
        hasBnds = errNC .EQ. NF_NOERR
        IF ( .NOT.hasBnds ) THEN
          WRITE(msgBuf,'(2A)') 'RNF_TIME_SETUP: RNF: a "fixed"',
     &      ' climatology (mitgcm_time_repeat = "annual") needs'
          CALL PRINT_ERROR( msgBuf, myThid )
          WRITE(msgBuf,'(2A)') 'RNF_TIME_SETUP: time_bnds: its span',
     &      ' is the repeat cycle (schema 1.0, section 7)'
          CALL PRINT_ERROR( msgBuf, myThid )
          errCount = errCount + 1
          RETURN
        ENDIF
      ENDIF

C--   3. The units of the time axis, which every non-constant mode
C     needs: the seconds per unit, and the reference date
      IF ( filePer .NE. 0. _d 0 ) THEN
        errNC = NF_INQ_VARID( fid, 'time', varId )
        CALL RNF_NC_ERROR( 'RNF_TIME_SETUP',
     &       'looking for the required variable time',
     &       errNC, myThid )
        CALL RNF_NC_ATT_TEXT( fid, varId, 'units', 'RNF_TIME_SETUP',
     O       unitAll, iL, found,
     I       myThid )
        IF ( .NOT.found ) THEN
          WRITE(msgBuf,'(2A)') 'RNF_TIME_SETUP: RNF: the variable',
     &      ' time has no attribute units'
          CALL PRINT_ERROR( msgBuf, myThid )
          errCount = errCount + 1
          RETURN
        ENDIF
        iS = INDEX( unitAll(1:iL), ' since ' )
        unitTxt = ' '
        IF ( iS .GE. 2 ) unitTxt = unitAll(1:iS-1)
        unitSec = 0. _d 0
        IF     ( unitTxt(1:5).EQ.'days ' .OR.
     &           unitTxt(1:4).EQ.'day ' ) THEN
          unitSec = 86400. _d 0
        ELSEIF ( unitTxt(1:6).EQ.'hours ' .OR.
     &           unitTxt(1:5).EQ.'hour ' ) THEN
          unitSec = 3600. _d 0
        ELSEIF ( unitTxt(1:8).EQ.'minutes ' .OR.
     &           unitTxt(1:7).EQ.'minute ' ) THEN
          unitSec = 60. _d 0
        ELSEIF ( unitTxt(1:8).EQ.'seconds ' .OR.
     &           unitTxt(1:7).EQ.'second ' ) THEN
          unitSec = 1. _d 0
        ENDIF
        IF ( iS.LT.2 .OR. unitSec.EQ.0. _d 0 ) THEN
          WRITE(msgBuf,'(3A)') 'RNF_TIME_SETUP: RNF: time:units = "',
     &      unitAll(1:iL), '" is not understood'
          CALL PRINT_ERROR( msgBuf, myThid )
          WRITE(msgBuf,'(2A)') 'RNF_TIME_SETUP: it must be',
     &      ' "<days|hours|minutes|seconds> since <date>"'
          CALL PRINT_ERROR( msgBuf, myThid )
          errCount = errCount + 1
          RETURN
        ENDIF

C--     Scan the reference date into at most 6 integers:
C       YYYY-MM-DD[ hh:mm[:ss]], with "T" and a trailing "Z" allowed
C       as separators (schema 1.0, section 6.2). A seventh number
C       means a fraction of a second or a time-zone offset, neither
C       of which this reader applies, so it is refused rather than
C       dropped.
        DO n = 1, 7
          num(n) = 0
        ENDDO
        nNum = 0
        inNum = .FALSE.
        DO i = iS+7, iL
          chr = unitAll(i:i)
          IF ( chr.GE.'0' .AND. chr.LE.'9' ) THEN
            IF ( .NOT.inNum ) THEN
              nNum = nNum + 1
              inNum = .TRUE.
              IF ( nNum .LE. 7 ) num(nNum) = 0
            ENDIF
            IF ( nNum .LE. 7 ) num(nNum) = num(nNum)*10
     &                                   + ( ICHAR(chr)-ICHAR('0') )
          ELSE
            inNum = .FALSE.
          ENDIF
        ENDDO
        IF ( nNum.LT.3 .OR. nNum.GT.6 ) THEN
          WRITE(msgBuf,'(3A)') 'RNF_TIME_SETUP: RNF: the reference',
     &      ' date of time:units = "', unitAll(1:iL)
          CALL PRINT_ERROR( msgBuf, myThid )
          WRITE(msgBuf,'(2A,I4,A)') 'RNF_TIME_SETUP: " must be',
     &      ' YYYY-MM-DD[ hh:mm[:ss]]; it has', nNum, ' number(s)'
          CALL PRINT_ERROR( msgBuf, myThid )
          errCount = errCount + 1
          RETURN
        ENDIF
        refDat1 = num(1)*10000 + num(2)*100 + num(3)
        refDat2 = num(4)*10000 + num(5)*100 + num(6)

C--     4. The calendar. Only a run with pkg/cal has one to compare
C       with, so a run without it does not demand a match
C       (schema 1.0, section 6.3).
        CALL RNF_NC_ATT_TEXT( fid, varId, 'calendar',
     I       'RNF_TIME_SETUP',
     O       unitTxt, iL, found,
     I       myThid )
        IF ( .NOT.found ) THEN
C-        a missing calendar means "standard", as in CF
          unitTxt = 'standard'
          iL = 8
        ENDIF
        IF ( useCAL ) THEN
#ifdef ALLOW_CAL
          calOk = .FALSE.
          IF ( unitTxt(1:iL).EQ.'standard' .OR.
     &         unitTxt(1:iL).EQ.'gregorian' .OR.
     &         unitTxt(1:iL).EQ.'proleptic_gregorian' ) THEN
            calOk = usingGregorianCalendar
          ELSEIF ( unitTxt(1:iL).EQ.'noleap' .OR.
     &             unitTxt(1:iL).EQ.'365_day' ) THEN
            calOk = usingNoLeapYearCal
          ELSEIF ( unitTxt(1:iL).EQ.'360_day' ) THEN
            calOk = usingModelCalendar
          ELSE
            WRITE(msgBuf,'(4A)') 'RNF_TIME_SETUP: RNF: calendar "',
     &        unitTxt(1:iL), '" of the time axis is not one of',
     &        ' standard/gregorian, noleap/365_day, 360_day'
            CALL PRINT_ERROR( msgBuf, myThid )
            errCount = errCount + 1
            RETURN
          ENDIF
          IF ( .NOT.calOk ) THEN
            WRITE(msgBuf,'(4A)') 'RNF_TIME_SETUP: RNF: calendar "',
     &        unitTxt(1:iL), '" of the time axis is not the calendar',
     &        ' of this run (data.cal, TheCalendar)'
            CALL PRINT_ERROR( msgBuf, myThid )
            errCount = errCount + 1
            RETURN
          ENDIF
#endif /* ALLOW_CAL */
        ELSE
          WRITE(msgBuf,'(2A)') 'RNF_TIME_SETUP: no pkg/cal: the time',
     &      ' axis is model time in seconds from its reference date'
          CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                        SQUEEZE_RIGHT, myThid )
          WRITE(msgBuf,'(3A)') 'RNF_TIME_SETUP: the calendar of the',
     &      ' file, "', unitTxt(1:iL)
          CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                        SQUEEZE_RIGHT, myThid )
          WRITE(msgBuf,'(2A)') 'RNF_TIME_SETUP: ", is not compared:',
     &      ' this run has no calendar to compare it with'
          CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                        SQUEEZE_RIGHT, myThid )
        ENDIF

C--     5. The time of record 1, and the span of the bounds when the
C       records are a fixed-period climatology
        iStart(1) = 1
        iCount(1) = 1
        errNC = NF_GET_VARA_DOUBLE( fid, varId, iStart, iCount,
     &                              tmpVal )
        CALL RNF_NC_ERROR( 'RNF_TIME_SETUP',
     &       'reading the first value of time', errNC, myThid )
        recSecs = tmpVal(1)*unitSec

        IF ( annual .AND. filePer.GT.0. _d 0 ) THEN
          errNC = NF_INQ_VARID( fid, 'time_bnds', varId )
          CALL RNF_NC_ERROR( 'RNF_TIME_SETUP',
     &         'looking for the variable time_bnds', errNC, myThid )
          errNC = NF_INQ_VARNDIMS( fid, varId, n )
          CALL RNF_NC_ERROR( 'RNF_TIME_SETUP',
     &         'reading the shape of time_bnds', errNC, myThid )
          IF ( n .NE. 2 ) THEN
            WRITE(msgBuf,'(2A)') 'RNF_TIME_SETUP: RNF: time_bnds',
     &        ' must have the dimensions (time, nv)'
            CALL PRINT_ERROR( msgBuf, myThid )
            errCount = errCount + 1
            RETURN
          ENDIF
C-        in Fortran the bound pair varies fastest
          iStart(1) = 1
          iStart(2) = 1
          iCount(1) = 2
          iCount(2) = 1
          errNC = NF_GET_VARA_DOUBLE( fid, varId, iStart, iCount,
     &                                tmpVal )
          CALL RNF_NC_ERROR( 'RNF_TIME_SETUP',
     &         'reading the first bounds of time_bnds',
     &         errNC, myThid )
          bnd0 = tmpVal(1)*unitSec
          iStart(2) = RNF_nRecFile
          errNC = NF_GET_VARA_DOUBLE( fid, varId, iStart, iCount,
     &                                tmpVal )
          CALL RNF_NC_ERROR( 'RNF_TIME_SETUP',
     &         'reading the last bounds of time_bnds',
     &         errNC, myThid )
          bnd1 = tmpVal(2)*unitSec
          fileCyc = bnd1 - bnd0
        ENDIF
      ENDIF

#endif /* HAVE_NETCDF */

C--   6. data.rnf overrides the file (model contract)
      RNF_recPeriod = filePer
      IF ( RNF_period .NE. UNSET_RL ) RNF_recPeriod = RNF_period
      RNF_recCycle = fileCyc
      IF ( RNF_repCycle .NE. UNSET_RL ) RNF_recCycle = RNF_repCycle

C--   7. Consistency of the effective values with what the file holds
      IF ( RNF_recPeriod .EQ. -12. _d 0 .AND.
     &     RNF_nRecFile .NE. 12 ) THEN
        WRITE(msgBuf,'(2A,I8,A)') 'RNF_TIME_SETUP: RNF: a monthly',
     &    ' climatology needs 12 records; the file has',
     &    RNF_nRecFile, ' record(s)'
        CALL PRINT_ERROR( msgBuf, myThid )
        WRITE(msgBuf,'(2A)') 'RNF_TIME_SETUP: record k is calendar',
     &    ' month k, January to December (pkg/cal/cal_getmonthsrec)'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
      ENDIF
      IF ( RNF_recPeriod.GT.0. _d 0 .AND.
     &     RNF_recCycle .GT.0. _d 0 ) THEN
C-      One record per period of the cycle: the wrap of
C       EXF_GetFFieldRec counts records from the start of the cycle,
C       so a cycle longer than the records would select a record the
C       file does not have, and a shorter one would never reach the
C       last records. The tolerance is the 1e-3 s of the schema.
        IF ( ABS( RNF_recCycle
     &            - RNF_nRecFile*RNF_recPeriod ) .GT. 1. _d -3 ) THEN
          WRITE(msgBuf,'(2A,1PE16.8,A,1PE16.8)') 'RNF_TIME_SETUP:',
     &      ' RNF: repeat cycle', RNF_recCycle,
     &      ' s is not the record count times the period,',
     &      RNF_nRecFile*RNF_recPeriod
          CALL PRINT_ERROR( msgBuf, myThid )
          WRITE(msgBuf,'(2A)') 'RNF_TIME_SETUP: a repeating series',
     &      ' must hold exactly one record per period of its cycle'
          CALL PRINT_ERROR( msgBuf, myThid )
          errCount = errCount + 1
        ENDIF
      ENDIF

C--   8. The time of record 1, in the form pkg/exf wants
      IF ( RNF_recPeriod .EQ. 0. _d 0 ) THEN
C-      constant: never used
        RNF_recStart = 0. _d 0
      ELSEIF ( .NOT.useCAL ) THEN
C-      No calendar: the file's time is model time in seconds counted
C       from the reference date of its units, and data.rnf may say so
C       directly instead (schema 1.0, section 14).
        RNF_recStart = recSecs
        IF ( RNF_startTime .NE. UNSET_RL ) RNF_recStart = RNF_startTime
        IF ( RNF_startDate1.NE.0 .OR. RNF_startDate2.NE.0 ) THEN
          WRITE(msgBuf,'(2A)') 'RNF_TIME_SETUP: RNF: RNF_startDate1',
     &      '/RNF_startDate2 need pkg/cal; use RNF_startTime [s]'
          CALL PRINT_ERROR( msgBuf, myThid )
          errCount = errCount + 1
        ENDIF
        IF ( RNF_recPeriod .LT. 0. _d 0 ) THEN
          WRITE(msgBuf,'(2A,1PE16.8,A)') 'RNF_TIME_SETUP: RNF:',
     &      ' a monthly period (', RNF_recPeriod, ') needs pkg/cal'
          CALL PRINT_ERROR( msgBuf, myThid )
          errCount = errCount + 1
        ENDIF
      ELSE
#ifdef ALLOW_CAL
C-      With a calendar, record 1 sits at a real date, and
C       EXF_GETFFIELD_START converts it to the model time the same way
C       it does for a dense runoffFile. data.rnf overrides the date.
        IF ( RNF_startTime .NE. UNSET_RL ) THEN
          WRITE(msgBuf,'(2A)') 'RNF_TIME_SETUP: RNF: RNF_startTime',
     &      ' cannot be set with pkg/cal (useCAL=T):'
          CALL PRINT_ERROR( msgBuf, myThid )
          WRITE(msgBuf,'(2A)') 'RNF_TIME_SETUP: the start time comes',
     &      ' from the file, or from RNF_startDate1/RNF_startDate2'
          CALL PRINT_ERROR( msgBuf, myThid )
          errCount = errCount + 1
        ENDIF
        IF ( RNF_startDate1 .NE. 0 ) THEN
          dat1 = RNF_startDate1
          dat2 = RNF_startDate2
          haveDate = .TRUE.
        ELSEIF ( refDat1/10000 .LT. 1583 ) THEN
C-        pkg/cal counts from 15 October 1582 (pkg/cal/cal_set.F), so
C         an earlier reference date cannot be converted. A constant
C         file carries 0001-01-01 and never reaches this point; a
C         timed file that does has to say its start date in data.rnf.
          WRITE(msgBuf,'(2A,I9,A)') 'RNF_TIME_SETUP: RNF: the',
     &      ' reference date of time:units,', refDat1,
     &      ', precedes the'
          CALL PRINT_ERROR( msgBuf, myThid )
          WRITE(msgBuf,'(2A)') 'RNF_TIME_SETUP: pkg/cal reference',
     &      ' date 15821015 and cannot be converted;'
          CALL PRINT_ERROR( msgBuf, myThid )
          WRITE(msgBuf,'(2A)') 'RNF_TIME_SETUP: give the date of',
     &      ' record 1 as RNF_startDate1/RNF_startDate2 (data.rnf)'
          CALL PRINT_ERROR( msgBuf, myThid )
          errCount = errCount + 1
        ELSE
C-        the real date of record 1 = reference date + its own time
          CALL CAL_FULLDATE( refDat1, refDat2, refArr, myThid )
          CALL CAL_TIMEINTERVAL( recSecs, 'secs', intv, myThid )
          CALL CAL_ADDTIME( refArr, intv, datArr, myThid )
          CALL CAL_CONVDATE( datArr, yy, mm, dd, ss, lp, wd, myThid )
          dat1 = yy*10000 + mm*100 + dd
          dat2 = (ss/3600)*10000 + (MOD(ss,3600)/60)*100
     &         + MOD(ss,60)
          haveDate = .TRUE.
        ENDIF
        IF ( haveDate ) THEN
          RNF_recDate1 = dat1
          RNF_recDate2 = dat2
C-        A monthly climatology is January to December in pkg/exf
C         (period -12 ignores the start date), so a file whose first
C         record is another month would be read with every record
C         shifted, silently.
          IF ( RNF_recPeriod.EQ.-12. _d 0 .AND.
     &         MOD( dat1/100, 100 ).NE.1 ) THEN
            WRITE(msgBuf,'(A,I9,A)') 'RNF_TIME_SETUP: RNF: record 1'
     &        //' of a monthly climatology is dated', dat1,
     &        ' but has'
            CALL PRINT_ERROR( msgBuf, myThid )
            WRITE(msgBuf,'(2A)') 'RNF_TIME_SETUP: to be in January:',
     &        ' record k is calendar month k (pkg/exf period -12)'
            CALL PRINT_ERROR( msgBuf, myThid )
            errCount = errCount + 1
          ENDIF
          RNF_recStart = UNSET_RL
          CALL EXF_GETFFIELD_START(
     I         RNF_useYearlyFiles, 'rnf', 'rnf_flux',
     I         RNF_recPeriod, dat1, dat2,
     U         RNF_recStart, errCount,
     I         myThid )
        ENDIF
#endif /* ALLOW_CAL */
      ENDIF

C--   9. Yearly files need a calendar: the file of a model year is
C     selected from the model date (pkg/exf/exf_check.F does the same
C     for useExfYearlyFields).
      IF ( RNF_useYearlyFiles .AND. .NOT.useCAL ) THEN
        WRITE(msgBuf,'(2A)') 'RNF_TIME_SETUP: RNF:',
     &    ' RNF_useYearlyFiles=.TRUE. needs pkg/cal (useCAL=T)'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
      ENDIF
      IF ( RNF_useYearlyFiles .AND. RNF_recCycle.NE.0. _d 0 ) THEN
        WRITE(msgBuf,'(2A)') 'RNF_TIME_SETUP: RNF:',
     &    ' RNF_useYearlyFiles=.TRUE. allows no repeat cycle'
        CALL PRINT_ERROR( msgBuf, myThid )
        errCount = errCount + 1
      ENDIF

C--   10. What was resolved, so that a log says which records a run
C     will read and from when
      WRITE(msgBuf,'(2A,1PE16.8,A,1PE16.8)') 'RNF_TIME_SETUP:',
     &  ' period =', RNF_recPeriod, ' s, repeat cycle =',
     &  RNF_recCycle
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(2A,1PE16.8,A,2I9)') 'RNF_TIME_SETUP:',
     &  ' record 1 at', RNF_recStart, ' s, date', RNF_recDate1,
     &  RNF_recDate2
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )

#endif /* ALLOW_RNF */

      RETURN
      END
