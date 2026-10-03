#ifdef ALLOW_RNF

CBOP
C !ROUTINE: RNF.h

C !DESCRIPTION: \bv
C     *==========================================================*
C     | RNF.h
C     | o Header of the sparse runoff (RNF) package.
C     |   Holds the run-time parameters read from data.rnf.
C     *==========================================================*
C \ev

C-----------------------------------------------------------------------
C--   Parameters that can be set in data.rnf, namelist /RNF_PARM01/
C     RNF_file           :: name of the sparse runoff NetCDF file, or
C                           its base name with RNF_useYearlyFiles; a
C                           blank name stops the run when useRNF=T
C     RNF_holdRecord     :: F: interpolate linearly in time between
C                           records, as pkg/exf does;
C                           T: hold each record over its interval
C     RNF_startDate1     :: start date of the series (YYYYMMDD);
C                           overrides the file, 0 = not set
C     RNF_startDate2     :: start date of the series (HHMMSS);
C                           overrides the file, 0 = not set
C     RNF_startTime      :: start time of the series (s); overrides
C                           the file, UNSET_RL = use the file
C     RNF_period         :: record period; overrides the file,
C                           UNSET_RL = use the file;
C                           0 = constant, > 0 = seconds,
C                           -12 = monthly climatology, -1 = monthly
C     RNF_repCycle       :: repeat cycle (s); overrides the file,
C                           UNSET_RL = use the file
C     RNF_useYearlyFiles :: append _YYYY to RNF_file by model year
C     RNF_useTemp        :: apply runoff_temperature if the file has it
C     RNF_useSalt        :: apply runoff_salinity if the file has it
C     RNF_usePtracers    :: apply runoff_ptracer_* if the file has them
C     RNF_monFreq        :: monitor interval (s), default monitorFreq
C     RNF_debugLev       :: message level, default debugLevel
C
C--   Internal parameters (not in the namelist)
C     RNFisON            :: package activation status
C
C--   Constants fixed by the model contract (not in the namelist)
C     RNF_fracTol        :: tolerance on the sum over the domain of
C                           the fractions of one source (must be 1)
C     RNF_areaTol        :: relative tolerance between the cell area
C                           stored in the file and rA
C-----------------------------------------------------------------------
CEOP

      LOGICAL RNFisON
      LOGICAL RNF_holdRecord
      LOGICAL RNF_useYearlyFiles
      LOGICAL RNF_useTemp
      LOGICAL RNF_useSalt
      LOGICAL RNF_usePtracers
      COMMON /RNF_PARM_L/
     &     RNFisON,
     &     RNF_holdRecord, RNF_useYearlyFiles,
     &     RNF_useTemp, RNF_useSalt, RNF_usePtracers

      INTEGER RNF_startDate1
      INTEGER RNF_startDate2
      INTEGER RNF_debugLev
      COMMON /RNF_PARM_I/
     &     RNF_startDate1, RNF_startDate2,
     &     RNF_debugLev

      _RL RNF_startTime
      _RL RNF_period
      _RL RNF_repCycle
      _RL RNF_monFreq
      COMMON /RNF_PARM_R/
     &     RNF_startTime, RNF_period, RNF_repCycle,
     &     RNF_monFreq

      CHARACTER*(MAX_LEN_FNAM) RNF_file
      COMMON /RNF_PARM_C/
     &     RNF_file

      _RL RNF_fracTol
      PARAMETER ( RNF_fracTol = 1. _d -6 )
      _RL RNF_areaTol
      PARAMETER ( RNF_areaTol = 1. _d -4 )

#endif /* ALLOW_RNF */
