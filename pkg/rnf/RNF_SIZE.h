C     *==========================================================*
C     | RNF_SIZE.h
C     | o Array bounds of the sparse runoff (RNF) package.
C     *==========================================================*
C     | Compile-time maxima with run-time counts inside them
C     | (package design, decision 9). A count that exceeds its
C     | bound stops the run, and the message prints the value
C     | that would be needed.
C     | This file must be included before RNF.h, which uses these
C     | bounds to dimension the per-tile lists.
C     *==========================================================*
C     RNF_nSrcTile :: maximum number of sources on one tile
C     RNF_nTgtTile :: maximum number of target entries on one tile
C     RNF_nBuf     :: chunk length for NetCDF reads and for the
C                     global fraction sum
C     RNF_nTr      :: maximum number of runoff tracers
C                     (no tracer is read yet: RUNOFF-013)
C     RNF_nFile    :: number of runoff files (1 in phase 1)

      INTEGER RNF_nSrcTile
      PARAMETER ( RNF_nSrcTile = 2000 )
      INTEGER RNF_nTgtTile
      PARAMETER ( RNF_nTgtTile = 10000 )
      INTEGER RNF_nBuf
      PARAMETER ( RNF_nBuf = 1000 )
      INTEGER RNF_nTr
      PARAMETER ( RNF_nTr = 1 )
      INTEGER RNF_nFile
      PARAMETER ( RNF_nFile = 1 )
