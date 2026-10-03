C     *==========================================================*
C     | RNF_SIZE.h
C     | o Array bounds of the sparse runoff (RNF) package.
C     *==========================================================*
C     | Placeholder values: no array is dimensioned with these
C     | bounds yet. The reader (RUNOFF-004) sets the values and
C     | adds the errors that stop the run and print the bound
C     | needed when a count exceeds its bound.
C     *==========================================================*
C     RNF_nSrcTile :: maximum number of sources on one tile
C     RNF_nTgtTile :: maximum number of target entries on one tile
C     RNF_nBuf     :: chunk length for NetCDF reads and for the
C                     global fraction sum
C     RNF_nTr      :: maximum number of runoff tracers
C     RNF_nFile    :: number of runoff files (1 for now)

      INTEGER RNF_nSrcTile
      PARAMETER ( RNF_nSrcTile = 1000 )
      INTEGER RNF_nTgtTile
      PARAMETER ( RNF_nTgtTile = 1000 )
      INTEGER RNF_nBuf
      PARAMETER ( RNF_nBuf = 1000 )
      INTEGER RNF_nTr
      PARAMETER ( RNF_nTr = 1 )
      INTEGER RNF_nFile
      PARAMETER ( RNF_nFile = 1 )
