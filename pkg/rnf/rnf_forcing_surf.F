#include "RNF_OPTIONS.h"
C--  File rnf_forcing_surf.F: the runoff salt and passive-tracer
C--  contributions to the model's surface forcing fields
C--   Contents
C--   o RNF_FORCING_SURF
C--   o RNF_FORCING_SURF_PTR
C--
C--  At a surface target the runoff salt and tracers enter the model
C--  the way every other surface flux does, through surfaceForcingS
C--  and surfaceForcingPTr (package design, decisions 3 and 4, as
C--  revised by RUNOFF-031), so that KPP, the SFLUX diagnostic and the
C--  ptracers surface path see them. The runoff heat takes the exf
C--  path instead (RNF_EXF_RUNOFF fills runoftemp, and EXF_MAPFIELDS
C--  adds it to Qnet), because that is the convention MITgcm already
C--  has for the heat of runoff; salt and tracers have none.
C--
C--  Each routine adds, to the surface field of each target cell,
C--
C--    surfaceForcingX += [ (mX) - m*X_ref ] * mass2rUnit
C--
C--  where (mX) is the mass flux of runoff times its property, summed
C--  over the sources that feed the cell, m that mass flux, and X_ref
C--  the property that the model's own freshwater formulation has
C--  already given this water. The bracket is therefore the mass flux
C--  times the difference between what the source carries and what the
C--  water was taken to carry. surfaceForcingS and surfaceForcingPTr
C--  are applied at the surface level of the column by APPLY_FORCING_S
C--  and PTRACERS_APPLY_FORCING, which is the level of every schema-1.0
C--  target (package design, decision 5).
C--
C--  The fields read here are the RNF_ap* set of RNF.h, which
C--  RNF_FIELDS_LOAD has put at the time level of the model's own
C--  freshwater flux and exchanged, so that the term is valid over the
C--  whole range EXTERNAL_FORCING_SURF fills, halo included, which is
C--  what KPP_CALC reads (package design, decision 3, "Time level").

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP 0
C !ROUTINE: RNF_FORCING_SURF

C !INTERFACE:
      SUBROUTINE RNF_FORCING_SURF(
     I                    bi, bj, iMin, iMax, jMin, jMax,
     I                    myTime, myIter, myThid )

C     !DESCRIPTION:
C     Add the runoff salt and tracer terms to surfaceForcingS and
C     surfaceForcingPTr. Called from EXTERNAL_FORCING_SURF at the end
C     of its tile loop, after every other contribution to those two
C     fields (PTRACERS_FORCING_SURF and SHELFICE_FORCING_SURF
C     included), so that none of them resets this one.
C
C     **S_ref** is the salinity the model has already given this
C     water (package design, decision 3):
C
C       o salt_EvPrRn when it is set, which is the default (0,
C         model/src/set_defaults.F:265): every branch then gives the
C         water that salinity;
C       o salt(i,j,ks) when salt_EvPrRn is unset and the formulation
C         is branch N or L, where the model adds nothing;
C       o convertFW2Salt in branch U, the uniform reference salinity
C         of the virtual salt flux.
C
C     Branch N is the nonlinear-free-surface (or pressure-coordinate)
C     real-freshwater-flux case, branch L the linear free surface with
C     convertFW2Salt = -1 and branch U the virtual salt flux
C     (model/src/external_forcing_surf.F:261-349). The local salinity
C     is read at ks, the level that routine reads it at.
C
C     Unlike the temperature, the salinity of runoff has a default:
C     water without a runoff_salinity value is fresh, (mS) = 0. With
C     the default salt_EvPrRn = 0 the term is then identically zero
C     and nothing is added (RNF_applyS is false); with any other value
C     the term is the correction from the salinity the model assumed
C     to zero, which is a real difference from the dense path and the
C     contract asks for it (package design, decision 3, "S = 0").
C
C     The term is not added to saltFlux, the field a user's saltflxfile
C     reaches surfaceForcingS through: pkg/seaice and pkg/thsice
C     overwrite saltFlux with their own ice-ocean salt flux before
C     EXTERNAL_FORCING_SURF reads it (pkg/seaice/seaice_growth.F:2007,
C     pkg/thsice/thsice_step_fwd.F:276), so a runoff salt added there
C     would be lost in every sea-ice run (package design, decision 3).
C
C     RNFsflx, the diagnostic of this term, is in the units and sign of
C     the model's SFLUX: g/m^2/s, > 0 increases salinity.

C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "DYNVARS.h"
#include "FFIELDS.h"
#include "RNF_SIZE.h"
#include "RNF.h"

C     !INPUT PARAMETERS:
C     bi,bj     :: Current tile indices
C     iMin,iMax :: range of tile for the surface forcing fields
C     jMin,jMax :: range of tile for the surface forcing fields
C     myTime    :: Current time in simulation
C     myIter    :: Current iteration number
C     myThid    :: my Thread Id number
      INTEGER bi, bj
      INTEGER iMin, iMax, jMin, jMax
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
CEOP

#ifdef ALLOW_RNF

C     !LOCAL VARIABLES:
C     ks        :: surface level index, as EXTERNAL_FORCING_SURF uses it
C     useSalt   :: S_ref is the local salinity
C     sRefLoc   :: S_ref where it does not depend on the cell
C     sLoc      :: the salt term [g/m^2/s]
      INTEGER i, j, ks
      LOGICAL useSalt
      _RL     sRefLoc
      _RL     sLoc(1-OLx:sNx+OLx,1-OLy:sNy+OLy)

      IF ( RNF_applyS ) THEN

       ks = 1
       IF ( usingPCoords ) ks = Nr

C--   Which salinity the model's own freshwater term used
       useSalt = .FALSE.
       sRefLoc = 0. _d 0
       IF ( salt_EvPrRn .NE. UNSET_RL ) THEN
        sRefLoc = salt_EvPrRn
       ELSEIF ( ( nonlinFreeSurf.GT.0 .OR. usingPCoords )
     &          .AND. useRealFreshWaterFlux ) THEN
        useSalt = .TRUE.
       ELSEIF ( convertFW2Salt .EQ. -1. _d 0 ) THEN
        useSalt = .TRUE.
       ELSE
        sRefLoc = convertFW2Salt
       ENDIF

       DO j = 1-OLy, sNy+OLy
        DO i = 1-OLx, sNx+OLx
         sLoc(i,j) = 0. _d 0
        ENDDO
       ENDDO
       DO j = jMin, jMax
        DO i = iMin, iMax
         IF ( useSalt ) THEN
          sLoc(i,j) = RNF_apXS(i,j,bi,bj)
     &              - RNF_apMflx(i,j,bi,bj)*salt(i,j,ks,bi,bj)
         ELSE
          sLoc(i,j) = RNF_apXS(i,j,bi,bj)
     &              - RNF_apMflx(i,j,bi,bj)*sRefLoc
         ENDIF
         surfaceForcingS(i,j,bi,bj) = surfaceForcingS(i,j,bi,bj)
     &                              + sLoc(i,j)*mass2rUnit
        ENDDO
       ENDDO

#ifdef ALLOW_DIAGNOSTICS
       IF ( useDiagnostics ) THEN
        CALL DIAGNOSTICS_FILL( sLoc, 'RNFsflx ', 0, 1, 2,
     &                         bi, bj, myThid )
       ENDIF
#endif /* ALLOW_DIAGNOSTICS */

      ENDIF

      IF ( RNF_nTrUse .GE. 1 ) THEN
        CALL RNF_FORCING_SURF_PTR( .TRUE.,
     I                    bi, bj, iMin, iMax, jMin, jMax,
     I                    myTime, myIter, myThid )
      ENDIF

#endif /* ALLOW_RNF */

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP 0
C !ROUTINE: RNF_FORCING_SURF_PTR

C !INTERFACE:
      SUBROUTINE RNF_FORCING_SURF_PTR( doDiags,
     I                    bi, bj, iMin, iMax, jMin, jMax,
     I                    myTime, myIter, myThid )

C     !DESCRIPTION:
C     Add the runoff tracer term to surfaceForcingPTr of every ptracer
C     that some runoff_ptracer_<NAME> variable of the file feeds:
C
C       surfaceForcingPTr(n) += [ (mC_n) - m*C_ref,n ] * mass2rUnit
C
C     **C_ref,n** is PTRACERS_EvPrRn(n) when it is set, and otherwise
C     the local pTracer in branches N and L and PTRACERS_ref(ks,n) in
C     branch U, which is what PTRACERS_FORCING_SURF uses
C     (pkg/ptracers/ptracers_forcing_surf.F:114-189). The branch test
C     is the same one as for salinity, so the three freshwater
C     formulations stay consistent between salt and tracers (package
C     design, decision 4).
C
C     Called from RNF_FORCING_SURF, i.e. after PTRACERS_FORCING_SURF
C     has set the field, and with doDiags false from
C     LONGSTEP_FORCING_SURF, because pkg/longstep zeroes and rebuilds
C     surfaceForcingPTr itself before it steps the ptracers
C     (pkg/longstep/longstep_thermodynamics.F:157), which would
C     otherwise erase the term. The diagnostic is filled on the first
C     call only, so that it is filled once per step. RNF_CHECK refuses
C     a runoff tracer under pkg/longstep with LS_nIter other than 1,
C     where this term (of the current step) and LS_fwFlux (a long-step
C     average) would belong to different times (RUNOFF-043).
C
C     The variable-to-ptracer match is made once, at init, by
C     RNF_NC_SERIES: RNF_trPtr(n) is the ptracer that runoff tracer n
C     feeds, and a name that matches none stops the run there. A
C     ptracer that no variable feeds gets no term and behaves as on
C     the dense path.
C
C     RNFtfNN, the diagnostic of runoff tracer NN (the NN-th
C     runoff_ptracer_<NAME> variable of the file, not ptracer NN), is
C     the bracket above in tracer units times kg/m^2/s, > 0 increases
C     the tracer; it is the surfaceForcingPTr contribution times
C     rUnit2mass, as SFLUX is surfaceForcingS times rUnit2mass.

C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#ifdef ALLOW_PTRACERS
# include "PTRACERS_SIZE.h"
# include "PTRACERS_PARAMS.h"
# include "PTRACERS_FIELDS.h"
#endif
#include "RNF_SIZE.h"
#include "RNF.h"

C     !INPUT PARAMETERS:
C     doDiags   :: fill the RNFtfNN diagnostics
C     bi,bj     :: Current tile indices
C     iMin,iMax :: range of tile for the surface forcing field
C     jMin,jMax :: range of tile for the surface forcing field
C     myTime    :: Current time in simulation
C     myIter    :: Current iteration number
C     myThid    :: my Thread Id number
      LOGICAL doDiags
      INTEGER bi, bj
      INTEGER iMin, iMax, jMin, jMax
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
CEOP

#if ( defined ALLOW_RNF && defined ALLOW_PTRACERS )

C     !LOCAL VARIABLES:
C     iRnf      :: runoff tracer index
C     iTracer   :: the ptracer it feeds
C     ks        :: surface level index, as PTRACERS_FORCING_SURF uses it
C     useLocal  :: C_ref is the local tracer concentration
C     cLoc      :: the tracer term [tracer units * kg/m^2/s]
      INTEGER i, j, ks
      INTEGER iRnf, iTracer
      LOGICAL useLocal
      _RL     cRefLoc
      _RL     cLoc(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      CHARACTER*8 diagNam

      ks = 1
      IF ( usingPCoords ) ks = Nr

      DO iRnf = 1, RNF_nTrUse
       iTracer = RNF_trPtr(iRnf)

       useLocal = .FALSE.
       cRefLoc = 0. _d 0
       IF ( PTRACERS_EvPrRn(iTracer) .NE. UNSET_RL ) THEN
        cRefLoc = PTRACERS_EvPrRn(iTracer)
       ELSEIF ( ( nonlinFreeSurf.GT.0 .OR. usingPCoords )
     &          .AND. useRealFreshWaterFlux ) THEN
        useLocal = .TRUE.
       ELSEIF ( convertFW2Salt .EQ. -1. _d 0 ) THEN
        useLocal = .TRUE.
       ELSE
        cRefLoc = PTRACERS_ref(ks,iTracer)
       ENDIF

       DO j = 1-OLy, sNy+OLy
        DO i = 1-OLx, sNx+OLx
         cLoc(i,j) = 0. _d 0
        ENDDO
       ENDDO
       DO j = jMin, jMax
        DO i = iMin, iMax
         IF ( useLocal ) THEN
          cLoc(i,j) = RNF_apXTr(i,j,bi,bj,iRnf)
     &              - RNF_apMflx(i,j,bi,bj)
     &               *pTracer(i,j,ks,bi,bj,iTracer)
         ELSE
          cLoc(i,j) = RNF_apXTr(i,j,bi,bj,iRnf)
     &              - RNF_apMflx(i,j,bi,bj)*cRefLoc
         ENDIF
         surfaceForcingPTr(i,j,bi,bj,iTracer) =
     &        surfaceForcingPTr(i,j,bi,bj,iTracer)
     &      + cLoc(i,j)*mass2rUnit
        ENDDO
       ENDDO

#ifdef ALLOW_DIAGNOSTICS
       IF ( useDiagnostics .AND. doDiags ) THEN
        WRITE(diagNam,'(A,I2.2,A)') 'RNFtf', iRnf, ' '
        CALL DIAGNOSTICS_FILL( cLoc, diagNam, 0, 1, 2,
     &                         bi, bj, myThid )
       ENDIF
#endif /* ALLOW_DIAGNOSTICS */

      ENDDO

#endif /* ALLOW_RNF and ALLOW_PTRACERS */

      RETURN
      END
