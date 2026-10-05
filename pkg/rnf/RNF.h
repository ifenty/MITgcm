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
C     RNF_lonLatTol      :: tolerance between the cell centre stored
C                           in the file (target_lon, target_lat) and
C                           the owning cell's XC,YC, as a fraction
C                           of that cell's own spacing
C                           MIN(dxF,dyF). 0.5 means "the stored
C                           centre has to lie inside the owning
C                           cell".
C                           The margin of a corrupted target_cell is
C                           the distance between the two centres
C                           over this threshold. It is NOT 2 in
C                           general: 1 + spacing(from)/spacing(to)
C                           holds only where the spacing is locally
C                           uniform, and at the cs32 facet corner
C                           dxF jumps from 120208 to 156359 m, where
C                           the margin is 1.78. What does hold is
C                           that the margin is strictly above 1 on
C                           any grid with positive cell sizes: the
C                           centres of two distinct cells are at
C                           least 0.5*(s_from + s_to) apart along
C                           the move, and MIN(dxF,dyF) of the
C                           destination is at most s_to, so the
C                           ratio is at least 1 + s_from/s_to. The
C                           measured floor over EVERY ordered pair
C                           of distinct cells is 1.779673 on cs32
C                           (37,742,592 pairs) and 1.999904 on
C                           lab_sea, so no corruption of a single
C                           target_cell - by one cell or by any
C                           other amount - can evade the check on
C                           either grid.
C                           MIN, not MAX, of the two spacings is
C                           necessary rather than merely cautious:
C                           under MAX the margin of a one-cell zonal
C                           move on a 1-degree global lat-lon grid
C                           is 2*cos(lat), i.e. 0.347 at 80N and
C                           0.035 at 89N, so MAX would defeat the
C                           guard on every high-latitude row.
C                           Because both sides scale with the local
C                           cell, the one constant needs no
C                           retuning between grids or resolutions
C     RNF_fluxMax        :: a flux above this value in absolute terms
C                           counts as a missing value (not allowed)
C     RNF_idLen          :: length of a source id in the model
C     RNF_maxErrMsg      :: messages one process prints per error
C                           counter, so that a large file cannot fill
C                           the log
C
C--   What the static read of RNF_file found (RNF_INIT_FIXED)
C     RNF_nSrcFile       :: number of sources in the file
C     RNF_nTgtFile       :: number of target entries in the file
C     RNF_nRecFile       :: number of time records in the file
C     RNF_nSrcProc       :: tile-source entries of this process
C                           (a source on two tiles counts twice)
C     RNF_nTgtProc       :: target entries owned by this process
C     RNF_nTgtOwned      :: target entries owned over all processes
C
C--   Per-tile lists built by RNF_INIT_FIXED. Index k runs over the
C     entries of tile (bi,bj); only what the tile owns is stored.
C     RNF_nSrc           :: number of sources on the tile
C     RNF_nTgt           :: number of target entries on the tile
C     RNF_srcGlob        :: 0-based source index in the file of local
C                           source k
C     RNF_tgtI, RNF_tgtJ :: local cell indices of target entry k
C     RNF_tgtK           :: model level of target entry k (the
C                           surface level of its column)
C     RNF_tgtSrc         :: local source index of target entry k
C     RNF_tgtFrac        :: fraction of that source sent to the cell
C     RNF_srcFracSum     :: fractions of local source k summed over
C                           the entries of this tile (fraction check)
C     RNF_srcFlux        :: current volume flux of local source k
C                           [m^3/s], set by RNF_NC_READ_FLUX
C
C--   Dense per-tile fields, rebuilt at every time step
C     RNF_vflx           :: runoff volume flux per unit area [m/s],
C                           sum_s flux_s*frac_s,c / rA; this is what
C                           RNF_EXF_RUNOFF puts in the exf field
C     RNF_mflx           :: runoff mass flux [kg/m^2/s],
C                           rhoConstFresh*RNF_vflx
C
C--   Work space of the global sums (filled by the master thread,
C     read by every thread inside GLOBAL_SUM_VECTOR_RL)
C     RNF_sumWrk         :: per-tile fraction sums of one chunk of
C                           sources
C     RNF_sumVec         :: their sum over tiles and processes
C     RNF_fluxBuf        :: NetCDF read buffer of one chunk of fluxes
C     RNF_tileSum        :: one value per tile, for GLOBAL_SUM_TILE_RL
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
      _RL RNF_lonLatTol
      PARAMETER ( RNF_lonLatTol = 0.5 _d 0 )
      _RL RNF_fluxMax
      PARAMETER ( RNF_fluxMax = 1. _d 30 )
      INTEGER RNF_idLen
      PARAMETER ( RNF_idLen = 64 )
      INTEGER RNF_maxErrMsg
      PARAMETER ( RNF_maxErrMsg = 20 )

C--   RNF_tableRead :: the target table was read and every entry of
C     it was valid, so the fraction sums are worth computing
C     RNF_lonLatChk :: the cell-centre check of RNF_INIT_FIXED ran.
C     It is optional in both directions: target_lon/target_lat are
C     optional in schema 1.0, and XC,YC are degrees only on a
C     spherical-polar or curvilinear grid. RNF_SUMMARY reports it so
C     that a run never looks guarded when it was not.
      LOGICAL RNF_tableRead
      LOGICAL RNF_lonLatChk
      COMMON /RNF_COUNT_L/ RNF_tableRead, RNF_lonLatChk

C--   Counts of the static read
      INTEGER RNF_nSrcFile
      INTEGER RNF_nTgtFile
      INTEGER RNF_nRecFile
      INTEGER RNF_nSrcProc
      INTEGER RNF_nTgtProc
      INTEGER RNF_nTgtOwned
      COMMON /RNF_COUNT_I/
     &     RNF_nSrcFile, RNF_nTgtFile, RNF_nRecFile,
     &     RNF_nSrcProc, RNF_nTgtProc, RNF_nTgtOwned

C--   Per-tile source and target lists
      INTEGER RNF_nSrc   (nSx,nSy)
      INTEGER RNF_nTgt   (nSx,nSy)
      INTEGER RNF_srcGlob(RNF_nSrcTile,nSx,nSy)
      INTEGER RNF_tgtI   (RNF_nTgtTile,nSx,nSy)
      INTEGER RNF_tgtJ   (RNF_nTgtTile,nSx,nSy)
      INTEGER RNF_tgtK   (RNF_nTgtTile,nSx,nSy)
      INTEGER RNF_tgtSrc (RNF_nTgtTile,nSx,nSy)
      COMMON /RNF_LIST_I/
     &     RNF_nSrc, RNF_nTgt, RNF_srcGlob,
     &     RNF_tgtI, RNF_tgtJ, RNF_tgtK, RNF_tgtSrc

      _RL RNF_tgtFrac   (RNF_nTgtTile,nSx,nSy)
      _RL RNF_srcFracSum(RNF_nSrcTile,nSx,nSy)
      _RL RNF_srcFlux   (RNF_nSrcTile,nSx,nSy)
      COMMON /RNF_LIST_R/
     &     RNF_tgtFrac, RNF_srcFracSum, RNF_srcFlux

C--   Dense per-tile fields
      _RL RNF_vflx(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL RNF_mflx(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      COMMON /RNF_FIELDS_R/
     &     RNF_vflx, RNF_mflx

C--   Work space of the global sums and of the NetCDF reads
      _RL RNF_sumWrk(nSx,nSy,RNF_nBuf)
      _RL RNF_sumVec(RNF_nBuf)
      _RL RNF_fluxBuf(RNF_nBuf)
      _RL RNF_tileSum(nSx,nSy)
      COMMON /RNF_WORK_R/
     &     RNF_sumWrk, RNF_sumVec, RNF_fluxBuf, RNF_tileSum

#endif /* ALLOW_RNF */
