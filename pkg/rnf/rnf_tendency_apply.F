#include "RNF_OPTIONS.h"
#ifdef ALLOW_EXF
# include "EXF_OPTIONS.h"
#endif
C--  File rnf_tendency_apply.F: Routines to apply RNF tendencies
C--   Contents
C--   o RNF_TENDENCY_APPLY_T
C--   o RNF_TENDENCY_APPLY_S
C--   o RNF_TENDENCY_APPLY_PTR
C--
C--  Each routine adds, at the target level of each target cell,
C--
C--    g_X += [ (mX) - m_X * X_ref ] * mass2rUnit
C--                  * recip_drF(k) * recip_hFacC(i,j,k,bi,bj)
C--
C--  where (mX) is the mass flux of runoff times its property,
C--  summed over the sources that feed the cell, m_X is that mass
C--  flux, and X_ref is the property that the model's own freshwater
C--  formulation has already given this water. The term is therefore
C--  the mass flux times the difference between what the source
C--  carries and what the water was taken to carry, which is the
C--  form the model uses for addMass
C--  (model/src/apply_forcing.F:504-531, 874-901).
C--
C--  The reference values, and the total each case adds up to, are
C--  the two tables of decision 3 of docs/package_design.md; the
C--  tracer case is decision 4. EXF_OPTIONS.h is included because
C--  one of the references depends on ALLOW_ATM_TEMP, as pkg/seaice
C--  includes it for the same kind of reason
C--  (pkg/seaice/seaice_growth.F:1-4).
C--
C--  The fields read here are the RNF_ap* set of RNF.h, which
C--  RNF_FIELDS_LOAD has already put at the time level of the
C--  model's own freshwater flux (package design, decision 3,
C--  "Time level"), so these routines need not know which level
C--  that is.

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP 0
C !ROUTINE: RNF_TENDENCY_APPLY_T

C !INTERFACE:
      SUBROUTINE RNF_TENDENCY_APPLY_T(
     U                    gT_arr,
     I                    iMin,iMax,jMin,jMax, k, bi, bj,
     I                    myTime, myIter, myThid )

C     !DESCRIPTION:
C     Add the runoff temperature term to the T (theta) tendency:
C
C       gT += [ (mT) - m_T*T_ref ] * mass2rUnit * D
C
C     with D = recip_drF(k)*recip_hFacC(i,j,k,bi,bj). Called from
C     APPLY_FORCING_T after the ICEFRONT call; works for one level at
C     a time, and does something only at the surface level of the
C     column, which is the level every schema-1.0 target has
C     (package design, decision 5).
C
C     **T_ref** is the temperature the model has already given this
C     water (package design, decision 3):
C
C       o theta(i,j,k), the ambient temperature, when temp_EvPrRn is
C         unset, because then the model adds nothing
C         (model/src/config_summary.F:407-409); and also when it is
C         set in a build with ALLOW_ATM_TEMP, because pkg/exf then
C         cancels the model's term for the runoff part through Qnet
C         (pkg/exf/exf_mapfields.F:175-185), so the water again
C         arrives at theta;
C       o temp_EvPrRn when it is set and ALLOW_ATM_TEMP is undefined,
C         where that cancellation is not compiled
C         (pkg/exf/exf_mapfields.F:132-198) and the dense path also
C         delivers runoff at temp_EvPrRn.
C
C     In both cases the three terms (model, exf, package) add up to
C     [(mT) - m*theta]*mass2rUnit in branches N and L, i.e. the heat
C     of the source water relative to the ambient temperature, and to
C     [(mT) - m*tRef(ks)]*mass2rUnit in branch U with temp_EvPrRn
C     set, which is the mix the dense path has there as well. That
C     cancellation is complete only where there is no sea ice: the
C     exf term is part of Qnet, which pkg/seaice and pkg/thsice scale
C     by the open-water fraction while this term is not scaled
C     (package design, decision 3, "Comparison with exf runoftemp").
C
C     A source whose temperature is missing, and every source in a
C     file without runoff_temperature, is in neither (mT) nor m_T, so
C     it contributes nothing here and enters at T_ref. Without the
C     variable the term is identically zero and the routine returns
C     at once, which is what keeps a file that carries only a volume
C     flux bit-for-bit on the dense path.

C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "RNF_SIZE.h"
#include "RNF.h"

C     !INPUT/OUTPUT PARAMETERS:
C     gT_arr    :: the tendency array
C     iMin,iMax :: working range of tile for applying forcing
C     jMin,jMax :: working range of tile for applying forcing
C     k         :: Current vertical level index
C     bi,bj     :: Current tile indices
C     myTime    :: Current time in simulation
C     myIter    :: Current iteration number
C     myThid    :: my Thread Id number
      _RL     gT_arr(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      INTEGER iMin, iMax, jMin, jMax
      INTEGER k, bi, bj
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
CEOP

#ifdef ALLOW_RNF

C     !LOCAL VARIABLES:
C     kSurface  :: surface level index, as APPLY_FORCING_T selects it
C     useTheta  :: T_ref is the local temperature
C     tLoc      :: the term before the thickness factor
      INTEGER i, j
      INTEGER kSurface
      LOGICAL useTheta
      _RL     tLoc(1-OLx:sNx+OLx,1-OLy:sNy+OLy)

      IF ( .NOT.RNF_applyT ) RETURN

C--   The surface level of this column, with the rule of
C     S/R APPLY_FORCING_T (model/src/apply_forcing.F:466-474),
C     evaluated at call time so that it stays right if pkg/shelfice
C     reassigns kSurfC during the run
      IF ( fluidIsAir ) THEN
       kSurface = 0
      ELSEIF ( usingZCoords .AND. useShelfIce ) THEN
       kSurface = -1
      ELSEIF ( usingPCoords ) THEN
       kSurface = Nr
      ELSE
       kSurface = 1
      ENDIF

#if defined(ALLOW_ATM_TEMP) && defined(ALLOW_RUNOFF)
      useTheta = .TRUE.
#else
      useTheta = temp_EvPrRn .EQ. UNSET_RL
#endif

      DO j = 1-OLy, sNy+OLy
       DO i = 1-OLx, sNx+OLx
        tLoc(i,j) = 0. _d 0
       ENDDO
      ENDDO

C--   Only target cells have a non-zero mass flux, and every target is
C     in the tile interior (RNF_INIT_FIXED places them there), so the
C     loop range is the interior.
      IF ( k.EQ.kSurface .OR. kSurface.EQ.-1 ) THEN
       DO j = 1, sNy
        DO i = 1, sNx
         IF ( kSurface.NE.-1 .OR. kSurfC(i,j,bi,bj).EQ.k ) THEN
          IF ( useTheta ) THEN
           tLoc(i,j) = ( RNF_apXT(i,j,bi,bj)
     &                 - RNF_apMflxT(i,j,bi,bj)*theta(i,j,k,bi,bj)
     &                 )*mass2rUnit
          ELSE
           tLoc(i,j) = ( RNF_apXT(i,j,bi,bj)
     &                 - RNF_apMflxT(i,j,bi,bj)*temp_EvPrRn
     &                 )*mass2rUnit
          ENDIF
          gT_arr(i,j) = gT_arr(i,j)
     &      + tLoc(i,j)*recip_drF(k)*_recip_hFacC(i,j,k,bi,bj)
         ENDIF
        ENDDO
       ENDDO
      ENDIF

#ifdef ALLOW_DIAGNOSTICS
      IF ( useDiagnostics ) THEN
       DO j = 1, sNy
        DO i = 1, sNx
         tLoc(i,j) = tLoc(i,j)
     &     *recip_drF(k)*_recip_hFacC(i,j,k,bi,bj)
        ENDDO
       ENDDO
       CALL DIAGNOSTICS_FILL( tLoc, 'RNFgT   ', k, 1, 2,
     &                        bi, bj, myThid )
      ENDIF
#endif /* ALLOW_DIAGNOSTICS */

#endif /* ALLOW_RNF */

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP 0
C !ROUTINE: RNF_TENDENCY_APPLY_S

C !INTERFACE:
      SUBROUTINE RNF_TENDENCY_APPLY_S(
     U                    gS_arr,
     I                    iMin,iMax,jMin,jMax, k, bi, bj,
     I                    myTime, myIter, myThid )

C     !DESCRIPTION:
C     Add the runoff salinity term to the S tendency:
C
C       gS += [ (mS) - m*S_ref ] * mass2rUnit * D
C
C     with D = recip_drF(k)*recip_hFacC(i,j,k,bi,bj). Called from
C     APPLY_FORCING_S after the ICEFRONT call; works for one level at
C     a time, and does something only at the surface level of the
C     column.
C
C     **S_ref** is the salinity the model has already given this
C     water (package design, decision 3):
C
C       o salt_EvPrRn when it is set, which is the default (0,
C         model/src/set_defaults.F:265): every branch then gives the
C         water that salinity;
C       o salt(i,j,k) when salt_EvPrRn is unset and the formulation
C         is branch N or L, where the model adds nothing;
C       o convertFW2Salt in branch U, the uniform reference salinity
C         of the virtual salt flux.
C
C     Branch N is the nonlinear-free-surface (or pressure-coordinate)
C     real-freshwater-flux case, branch L the linear free surface with
C     convertFW2Salt = -1 and branch U the virtual salt flux
C     (model/src/external_forcing_surf.F:261-349).
C
C     Unlike the temperature, the salinity of runoff has a default:
C     water without a runoff_salinity value is fresh, (mS) = 0. With
C     the default salt_EvPrRn = 0 the term is then identically zero
C     and the routine returns at once; with any other value the term
C     is the correction from the salinity the model assumed to zero,
C     which is a real difference from the dense path and the contract
C     asks for it (package design, decision 3, "S = 0").

C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "RNF_SIZE.h"
#include "RNF.h"

C     !INPUT/OUTPUT PARAMETERS:
C     gS_arr    :: the tendency array
C     iMin,iMax :: working range of tile for applying forcing
C     jMin,jMax :: working range of tile for applying forcing
C     k         :: Current vertical level index
C     bi,bj     :: Current tile indices
C     myTime    :: Current time in simulation
C     myIter    :: Current iteration number
C     myThid    :: my Thread Id number
      _RL     gS_arr(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      INTEGER iMin, iMax, jMin, jMax
      INTEGER k, bi, bj
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
CEOP

#ifdef ALLOW_RNF

C     !LOCAL VARIABLES:
C     kSurface  :: surface level index, as APPLY_FORCING_S selects it
C     useSalt   :: S_ref is the local salinity
C     sRefLoc   :: S_ref where it does not depend on the cell
      INTEGER i, j
      INTEGER kSurface
      LOGICAL useSalt
      _RL     sRefLoc
      _RL     sLoc(1-OLx:sNx+OLx,1-OLy:sNy+OLy)

      IF ( .NOT.RNF_applyS ) RETURN

      IF ( fluidIsAir ) THEN
       kSurface = 0
      ELSEIF ( usingZCoords .AND. useShelfIce ) THEN
       kSurface = -1
      ELSEIF ( usingPCoords ) THEN
       kSurface = Nr
      ELSE
       kSurface = 1
      ENDIF

C--   Which salinity the model's own freshwater term used
      useSalt = .FALSE.
      sRefLoc = 0. _d 0
      IF ( salt_EvPrRn .NE. UNSET_RL ) THEN
       sRefLoc = salt_EvPrRn
      ELSEIF ( ( nonlinFreeSurf.GT.0 .OR. usingPCoords )
     &         .AND. useRealFreshWaterFlux ) THEN
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

      IF ( k.EQ.kSurface .OR. kSurface.EQ.-1 ) THEN
       DO j = 1, sNy
        DO i = 1, sNx
         IF ( kSurface.NE.-1 .OR. kSurfC(i,j,bi,bj).EQ.k ) THEN
          IF ( useSalt ) THEN
           sLoc(i,j) = ( RNF_apXS(i,j,bi,bj)
     &                 - RNF_apMflx(i,j,bi,bj)*salt(i,j,k,bi,bj)
     &                 )*mass2rUnit
          ELSE
           sLoc(i,j) = ( RNF_apXS(i,j,bi,bj)
     &                 - RNF_apMflx(i,j,bi,bj)*sRefLoc
     &                 )*mass2rUnit
          ENDIF
          gS_arr(i,j) = gS_arr(i,j)
     &      + sLoc(i,j)*recip_drF(k)*_recip_hFacC(i,j,k,bi,bj)
         ENDIF
        ENDDO
       ENDDO
      ENDIF

#ifdef ALLOW_DIAGNOSTICS
      IF ( useDiagnostics ) THEN
       DO j = 1, sNy
        DO i = 1, sNx
         sLoc(i,j) = sLoc(i,j)
     &     *recip_drF(k)*_recip_hFacC(i,j,k,bi,bj)
        ENDDO
       ENDDO
       CALL DIAGNOSTICS_FILL( sLoc, 'RNFgS   ', k, 1, 2,
     &                        bi, bj, myThid )
      ENDIF
#endif /* ALLOW_DIAGNOSTICS */

#endif /* ALLOW_RNF */

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP 0
C !ROUTINE: RNF_TENDENCY_APPLY_PTR

C !INTERFACE:
      SUBROUTINE RNF_TENDENCY_APPLY_PTR(
     U                    gPtr_arr,
     I                    iMin,iMax,jMin,jMax, k, bi, bj,
     I                    iTracer, myTime, myIter, myThid )

C     !DESCRIPTION:
C     Add the runoff tracer term to the tendency of passive tracer
C     iTracer:
C
C       gC_n += [ (mC_n) - m*C_ref,n ] * mass2rUnit * D
C
C     Called from PTRACERS_APPLY_FORCING beside the GCHEM call; works
C     for one level and one tracer at a time, and does something only
C     for a tracer that some runoff_ptracer_<NAME> variable of the
C     file feeds, at the surface level of the column.
C
C     **C_ref,n** is PTRACERS_EvPrRn(n) when it is set, and otherwise
C     the local pTracer in branches N and L and PTRACERS_ref(ks,n) in
C     branch U, which is what PTRACERS_FORCING_SURF uses
C     (pkg/ptracers/ptracers_forcing_surf.F:114-189). The branch test
C     is the same one as for salinity, so the three freshwater
C     formulations stay consistent between T, S and tracers
C     (package design, decision 4).
C
C     The variable-to-ptracer match is made once, at init, by
C     RNF_NC_SERIES: RNF_trPtr(n) is the ptracer that runoff tracer n
C     feeds, and a name that matches none stops the run there. A
C     ptracer that no variable feeds gets no term and behaves as on
C     the dense path.

C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#ifdef ALLOW_PTRACERS
# include "PTRACERS_SIZE.h"
# include "PTRACERS_PARAMS.h"
# include "PTRACERS_FIELDS.h"
#endif
#include "RNF_SIZE.h"
#include "RNF.h"

C     !INPUT/OUTPUT PARAMETERS:
C     gPtr_arr  :: the tendency array
C     iMin,iMax :: working range of tile for applying forcing
C     jMin,jMax :: working range of tile for applying forcing
C     k         :: Current vertical level index
C     bi,bj     :: Current tile indices
C     iTracer   :: passive tracer number
C     myTime    :: Current time in simulation
C     myIter    :: Current iteration number
C     myThid    :: my Thread Id number
      _RL     gPtr_arr(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      INTEGER iMin, iMax, jMin, jMax
      INTEGER k, bi, bj
      INTEGER iTracer
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
CEOP

#if ( defined ALLOW_RNF && defined ALLOW_PTRACERS )

C     !LOCAL VARIABLES:
C     iRnf      :: runoff tracer that feeds ptracer iTracer (0: none)
C     kSurface  :: surface level index, as PTRACERS_APPLY_FORCING
C                  selects it
C     useLocal  :: C_ref is the local tracer concentration
      INTEGER i, j, iTr
      INTEGER iRnf, kSurface, ks
      LOGICAL useLocal
      _RL     cRefLoc
      _RL     cLoc(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      CHARACTER*8 diagNam

      iRnf = 0
      DO iTr = 1, RNF_nTrUse
        IF ( RNF_trPtr(iTr) .EQ. iTracer ) iRnf = iTr
      ENDDO
      IF ( iRnf .EQ. 0 ) RETURN

      IF ( fluidIsAir ) THEN
       kSurface = 0
      ELSEIF ( usingZCoords .AND. useShelfIce ) THEN
       kSurface = -1
      ELSEIF ( usingPCoords ) THEN
       kSurface = Nr
      ELSE
       kSurface = 1
      ENDIF
      ks = 1
      IF ( usingPCoords ) ks = Nr

      useLocal = .FALSE.
      cRefLoc = 0. _d 0
      IF ( PTRACERS_EvPrRn(iTracer) .NE. UNSET_RL ) THEN
       cRefLoc = PTRACERS_EvPrRn(iTracer)
      ELSEIF ( ( nonlinFreeSurf.GT.0 .OR. usingPCoords )
     &         .AND. useRealFreshWaterFlux ) THEN
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

      IF ( k.EQ.kSurface .OR. kSurface.EQ.-1 ) THEN
       DO j = 1, sNy
        DO i = 1, sNx
         IF ( kSurface.NE.-1 .OR. kSurfC(i,j,bi,bj).EQ.k ) THEN
          IF ( useLocal ) THEN
           cLoc(i,j) = ( RNF_apXTr(i,j,bi,bj,iRnf)
     &                 - RNF_apMflx(i,j,bi,bj)
     &                  *pTracer(i,j,k,bi,bj,iTracer)
     &                 )*mass2rUnit
          ELSE
           cLoc(i,j) = ( RNF_apXTr(i,j,bi,bj,iRnf)
     &                 - RNF_apMflx(i,j,bi,bj)*cRefLoc
     &                 )*mass2rUnit
          ENDIF
          gPtr_arr(i,j) = gPtr_arr(i,j)
     &      + cLoc(i,j)*recip_drF(k)*_recip_hFacC(i,j,k,bi,bj)
         ENDIF
        ENDDO
       ENDDO
      ENDIF

#ifdef ALLOW_DIAGNOSTICS
      IF ( useDiagnostics ) THEN
       DO j = 1, sNy
        DO i = 1, sNx
         cLoc(i,j) = cLoc(i,j)
     &     *recip_drF(k)*_recip_hFacC(i,j,k,bi,bj)
        ENDDO
       ENDDO
       WRITE(diagNam,'(A,I2.2,A)') 'RNFtr', iRnf, ' '
       CALL DIAGNOSTICS_FILL( cLoc, diagNam, k, 1, 2,
     &                        bi, bj, myThid )
      ENDIF
#endif /* ALLOW_DIAGNOSTICS */

#endif /* ALLOW_RNF and ALLOW_PTRACERS */

      RETURN
      END
