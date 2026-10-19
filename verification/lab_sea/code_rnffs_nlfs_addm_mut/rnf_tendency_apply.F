#include "RNF_OPTIONS.h"
C--  File rnf_tendency_apply.F: the heat, salt and tracer terms of the
C--  generalized runoff route (RNF_applyAsRain = F, the default;
C--  package design, decisions 3 and 4; RUNOFF-046)
C--   Contents
C--   o RNF_TENDENCY_APPLY_T
C--   o RNF_TENDENCY_APPLY_S
C--   o RNF_TENDENCY_APPLY_PTR
C--
C--  Restored from the routines of the same names that RUNOFF-031
C--  removed (fork commit 46684c3af, which moved the then surface-only
C--  terms onto the surface flux fields; that route survives as the
C--  legacy one, RNF_applyAsRain = T) and generalized in two ways:
C--
C--  1. The level. The old routines acted only at the surface level
C--     of the column. These act, for each target entry e of the
C--     tile, at that entry's own level RNF_tgtK(e), so entries of one
C--     cell (i,j) at different levels each land where the file puts
C--     them. Each routine works on one level k at a time, as it is
C--     called from the k loop of its caller, and adds the entries
C--     with RNF_tgtK = k.
C--
C--  2. The reference value. The old routines added
C--     [(mX) - m*X_ref]*mass2rUnit*D, with X_ref the property the
C--     model's own SURFACE freshwater term (EmPmR, temp_EvPrRn,
C--     salt_EvPrRn, convertFW2Salt) had already given that water.
C--     On this route the water is not a surface freshwater flux but an
C--     addMass source (RNF_ADD_MASS), and the model treats an addMass
C--     source as water at the LOCAL value of each property: the
C--     tracer tendency is computed in advective form, which removes
C--     X*div(U) (pkg/generic_advdiff/gad_calc_rhs.F:768-779, and the
C--     same in gad_advection.F), so a mass source changes no
C--     concentration by itself. The physical term of a source of
C--     water at X_s entering ambient water at X(i,j,k) is therefore
C--
C--       g_X += [ (mX) - m_X*X(i,j,k) ]*mass2rUnit
C--              *recip_drF(k)*recip_hFacC(i,j,k)
C--
C--     with (mX) = sum_e m_e*X_s and m_X = sum_e m_e over the entries
C--     of the cell and level that carry the property, and m_e =
C--     rhoConstFresh*flux_s*frac_e/rA. No branch on the freshwater
C--     formulation is needed or wanted: the N/L/U branches of the
C--     legacy route exist only to match what the surface freshwater
C--     term assumes, and an addMass source has no such assumption.
C--
C--     One exception, which these routines cancel: the model's own
C--     generic addMass term (model/src/apply_forcing.F:504-531 for T,
C--     :874-901 for S, and the same in external_forcing.F), which
C--     adds addMass*(temp_addMass - X')*mass2rUnit*... with X' =
C--     theta in the nonlinear-free-surface (or convertFW2Salt = -1)
C--     branch and tRef(k) otherwise, whenever selectAddFluid /= 0 and
C--     temp_addMass (salt_addMass) is set. That term acts on the
C--     whole addMass array, the runoff share included, with one
C--     global property, which cannot be the per-source one; so each
C--     routine subtracts exactly that term for the runoff share
C--     m_e, with the model's own branch test. salt_addMass is set in
C--     every default run (model/src/ini_parms.F:1580 copies
C--     salt_EvPrRn, default 0, into it), so for salinity this
C--     cancellation is the normal case, not a corner. pkg/ptracers
C--     has no such term, so the tracer routine has nothing to cancel.
C--
C--  The diagnostics RNFgT, RNFgS and RNFtrNN are the physical term
C--  above, i.e. what the runoff does to the property in total, not
C--  the share of it these routines add after the cancellation.
C--
C--  The fields read are the per-tile source series (RNF_srcFlux,
C--  RNF_srcTemp, RNF_srcTvld, RNF_srcSalt, RNF_srcTrc) of the current
C--  step, built by RNF_FIELDS_LOAD, which is the step whose addMass
C--  RNF_ADD_MASS filled: these terms and the model's own addMass
C--  terms belong to the same step. Without staggerTimeStep the model
C--  warns that addMass acts on T and S one step out of phase with
C--  the flow (model/src/config_check.F:846-854); these terms then
C--  share that lag with the model's own.
C--
C--  Every routine is a no-op on the legacy route.

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP 0
C !ROUTINE: RNF_TENDENCY_APPLY_T

C !INTERFACE:
      SUBROUTINE RNF_TENDENCY_APPLY_T(
     U                    gT_arr,
     I                    iMin,iMax,jMin,jMax, k, bi, bj,
     I                    myTime, myIter, myThid )

C     !DESCRIPTION:
C     Add the runoff temperature term at level k to the T (theta)
C     tendency:
C
C       gT += [ (mT) - m_T*theta(k) ]*mass2rUnit*D - G
C
C     with D = recip_drF(k)*recip_hFacC(i,j,k,bi,bj) and G the
C     model's generic addMass term on the runoff share (file header).
C     Called from APPLY_FORCING_T (and EXTERNAL_FORCING_T) after the
C     ICEFRONT call.
C
C     A source whose temperature is missing, and every source of a
C     file without runoff_temperature, is in neither (mT) nor m_T: it
C     enters at the ambient temperature of the target cell. Without
C     the variable and without a generic term the routine returns at
C     once.

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
C     n, kSrc   :: target entry and its local source
C     genTerm   :: the model's generic addMass T term is active
C     genLoc    :: ... and uses theta (otherwise tRef(k))
C     mE        :: mass flux of one entry [kg/m^2/s]
C     tLoc      :: the physical term before the thickness factor
C     gLoc      :: the generic term on the runoff share, idem
      INTEGER i, j, n, kSrc
      LOGICAL genTerm, genLoc
      _RL     mE
      _RL     tLoc(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     gLoc(1-OLx:sNx+OLx,1-OLy:sNy+OLy)

      IF ( RNF_applyAsRain ) RETURN

C--   The model's own term, with its own test
C     (model/src/apply_forcing.F:504-507)
      genTerm = .FALSE.
      genLoc  = .FALSE.
#ifdef ALLOW_ADDFLUID
      genTerm = selectAddFluid.NE.0 .AND. temp_addMass.NE.UNSET_RL
      genLoc  = ( selectAddFluid.GE.1 .AND. nonlinFreeSurf.GT.0 )
     &          .OR. convertFW2Salt.EQ.-1. _d 0
#endif /* ALLOW_ADDFLUID */
      IF ( .NOT.RNF_applyT .AND. .NOT.genTerm ) RETURN

      DO j = 1-OLy, sNy+OLy
       DO i = 1-OLx, sNx+OLx
        tLoc(i,j) = 0. _d 0
        gLoc(i,j) = 0. _d 0
       ENDDO
      ENDDO

C--   The entries at this level. Every target is in the tile interior
C     (RNF_INIT_FIXED places them there).
      DO n = 1, RNF_nTgt(bi,bj)
       IF ( RNF_tgtK(n,bi,bj) .EQ. k ) THEN
        i = RNF_tgtI(n,bi,bj)
        j = RNF_tgtJ(n,bi,bj)
        kSrc = RNF_tgtSrc(n,bi,bj)
        mE = rhoConstFresh*RNF_srcFlux(kSrc,bi,bj)
     &      *RNF_tgtFrac(n,bi,bj)*recip_rA(i,j,bi,bj)
        IF ( RNF_applyT ) THEN
          tLoc(i,j) = tLoc(i,j) + mE*RNF_srcTvld(kSrc,bi,bj)
     &      *( RNF_srcTemp(kSrc,bi,bj) - theta(i,j,k,bi,bj) )
        ENDIF
        IF ( genTerm ) THEN
         IF ( genLoc ) THEN
          gLoc(i,j) = gLoc(i,j)
     &      + mE*( temp_addMass - theta(i,j,k,bi,bj) )
         ELSE
          gLoc(i,j) = gLoc(i,j) + mE*( temp_addMass - tRef(k) )
         ENDIF
        ENDIF
       ENDIF
      ENDDO

      IF ( nonlinFreeSurf.GT.0 ) THEN
C       mutant: the generic term is not cancelled
       DO j = 1-OLy, sNy+OLy
        DO i = 1-OLx, sNx+OLx
         gLoc(i,j) = 0. _d 0
        ENDDO
       ENDDO
      ENDIF
      DO j = 1, sNy
       DO i = 1, sNx
        tLoc(i,j) = tLoc(i,j)*mass2rUnit
     &             *recip_drF(k)*_recip_hFacC(i,j,k,bi,bj)
        gT_arr(i,j) = gT_arr(i,j) + tLoc(i,j)
     &    - gLoc(i,j)*mass2rUnit
     &      *recip_drF(k)*_recip_hFacC(i,j,k,bi,bj)
       ENDDO
      ENDDO

#ifdef ALLOW_DIAGNOSTICS
      IF ( useDiagnostics ) THEN
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
C     Add the runoff salinity term at level k to the S tendency:
C
C       gS += [ (mS) - m*salt(k) ]*mass2rUnit*D - G
C
C     with D = recip_drF(k)*recip_hFacC(i,j,k,bi,bj) and G the
C     model's generic addMass term on the runoff share (file header),
C     which is active in a default run (salt_addMass = 0). Called from
C     APPLY_FORCING_S (and EXTERNAL_FORCING_S) after the ICEFRONT
C     call.
C
C     Unlike the temperature, the salinity of runoff has a default:
C     water without a runoff_salinity value is fresh, (mS) = 0, so the
C     term is the dilution -m*salt(k) of the target cell. It is never
C     identically zero while there is runoff on the level.

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
C     n, kSrc   :: target entry and its local source
C     genTerm   :: the model's generic addMass S term is active
C     genLoc    :: ... and uses salt (otherwise sRef(k))
C     mE        :: mass flux of one entry [kg/m^2/s]
C     sSrc      :: salinity of its source
C     sLoc      :: the physical term before the thickness factor
C     gLoc      :: the generic term on the runoff share, idem
      INTEGER i, j, n, kSrc
      LOGICAL genTerm, genLoc
      _RL     mE, sSrc
      _RL     sLoc(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     gLoc(1-OLx:sNx+OLx,1-OLy:sNy+OLy)

      IF ( RNF_applyAsRain ) RETURN

C--   The model's own term, with its own test
C     (model/src/apply_forcing.F:875-877)
      genTerm = .FALSE.
      genLoc  = .FALSE.
#ifdef ALLOW_ADDFLUID
      genTerm = selectAddFluid.NE.0 .AND. salt_addMass.NE.UNSET_RL
      genLoc  = ( selectAddFluid.GE.1 .AND. nonlinFreeSurf.GT.0 )
     &          .OR. convertFW2Salt.EQ.-1. _d 0
#endif /* ALLOW_ADDFLUID */

      DO j = 1-OLy, sNy+OLy
       DO i = 1-OLx, sNx+OLx
        sLoc(i,j) = 0. _d 0
        gLoc(i,j) = 0. _d 0
       ENDDO
      ENDDO

      DO n = 1, RNF_nTgt(bi,bj)
       IF ( RNF_tgtK(n,bi,bj) .EQ. k ) THEN
        i = RNF_tgtI(n,bi,bj)
        j = RNF_tgtJ(n,bi,bj)
        kSrc = RNF_tgtSrc(n,bi,bj)
        mE = rhoConstFresh*RNF_srcFlux(kSrc,bi,bj)
     &      *RNF_tgtFrac(n,bi,bj)*recip_rA(i,j,bi,bj)
        sSrc = 0. _d 0
        IF ( RNF_hasSalt ) sSrc = RNF_srcSalt(kSrc,bi,bj)
        sLoc(i,j) = sLoc(i,j) + mE*( sSrc - salt(i,j,k,bi,bj) )
        IF ( genTerm ) THEN
         IF ( genLoc ) THEN
          gLoc(i,j) = gLoc(i,j)
     &      + mE*( salt_addMass - salt(i,j,k,bi,bj) )
         ELSE
          gLoc(i,j) = gLoc(i,j) + mE*( salt_addMass - sRef(k) )
         ENDIF
        ENDIF
       ENDIF
      ENDDO

      DO j = 1, sNy
       DO i = 1, sNx
        sLoc(i,j) = sLoc(i,j)*mass2rUnit
     &             *recip_drF(k)*_recip_hFacC(i,j,k,bi,bj)
        gS_arr(i,j) = gS_arr(i,j) + sLoc(i,j)
     &    - gLoc(i,j)*mass2rUnit
     &      *recip_drF(k)*_recip_hFacC(i,j,k,bi,bj)
       ENDDO
      ENDDO

#ifdef ALLOW_DIAGNOSTICS
      IF ( useDiagnostics ) THEN
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
C     Add the runoff tracer term at level k to the tendency of passive
C     tracer iTracer:
C
C       gC_n += [ (mC_n) - m*pTracer(k,n) ]*mass2rUnit*D
C
C     Called from PTRACERS_APPLY_FORCING beside the GCHEM call; works
C     for one level and one tracer at a time, and does something only
C     for a tracer that some runoff_ptracer_<NAME> variable of the
C     file feeds (RNF_trPtr, matched once at init by RNF_NC_SERIES). A
C     ptracer that no variable feeds gets no term: its runoff water
C     enters at the ambient concentration. pkg/ptracers has no generic
C     addMass term, so there is nothing to cancel.

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
C     n, kSrc   :: target entry and its local source
C     mE        :: mass flux of one entry [kg/m^2/s]
C     cLoc      :: the term
      INTEGER i, j, n, kSrc, iTr, iRnf
      _RL     mE
      _RL     cLoc(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      CHARACTER*8 diagNam

      IF ( RNF_applyAsRain ) RETURN

      iRnf = 0
      DO iTr = 1, RNF_nTrUse
        IF ( RNF_trPtr(iTr) .EQ. iTracer ) iRnf = iTr
      ENDDO
      IF ( iRnf .EQ. 0 ) RETURN

      DO j = 1-OLy, sNy+OLy
       DO i = 1-OLx, sNx+OLx
        cLoc(i,j) = 0. _d 0
       ENDDO
      ENDDO

      DO n = 1, RNF_nTgt(bi,bj)
       IF ( RNF_tgtK(n,bi,bj) .EQ. k ) THEN
        i = RNF_tgtI(n,bi,bj)
        j = RNF_tgtJ(n,bi,bj)
        kSrc = RNF_tgtSrc(n,bi,bj)
        mE = rhoConstFresh*RNF_srcFlux(kSrc,bi,bj)
     &      *RNF_tgtFrac(n,bi,bj)*recip_rA(i,j,bi,bj)
        cLoc(i,j) = cLoc(i,j) + mE
     &    *( RNF_srcTrc(kSrc,bi,bj,iRnf)
     &     - pTracer(i,j,k,bi,bj,iTracer) )
       ENDIF
      ENDDO

      DO j = 1, sNy
       DO i = 1, sNx
        cLoc(i,j) = cLoc(i,j)*mass2rUnit
     &             *recip_drF(k)*_recip_hFacC(i,j,k,bi,bj)
        gPtr_arr(i,j) = gPtr_arr(i,j) + cLoc(i,j)
       ENDDO
      ENDDO

#ifdef ALLOW_DIAGNOSTICS
      IF ( useDiagnostics ) THEN
       WRITE(diagNam,'(A,I2.2,A)') 'RNFtr', iRnf, ' '
       CALL DIAGNOSTICS_FILL( cLoc, diagNam, k, 1, 2,
     &                        bi, bj, myThid )
      ENDIF
#endif /* ALLOW_DIAGNOSTICS */

#endif /* ALLOW_RNF and ALLOW_PTRACERS */

      RETURN
      END
