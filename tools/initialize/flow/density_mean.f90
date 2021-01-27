#include "types.h"
#include "dns_const.h"
#include "dns_error.h"

!########################################################################
!# DESCRIPTION
!#
!# Setting the density field.
!# The default case is when temperature field is set similar to a passive
!# scalar, and the density is then derived taking pressure constant.
!# The inputs are however density maximum/minimum:
!#
!# T = T_0 + s(T_1-T_0), s.t. 0<=s<=1, =>
!#           => rho = W / [ W_0/rho_0 + s(W_1/rho_1-W_0/rho_0) ]
!#
!# This allow to compare the evolution of the passive scalar with that
!# of the temperature.
!#
!# If volumetric forces are present, then rho is computed from
!# dp/dy = rho g, to ensure hydrostatic equilibrium including numerical
!# errors.
!#
!# Spatial case with mulstispecies is not jet done
!#
!########################################################################
SUBROUTINE DENSITY_MEAN(rho, p,T,s, txc, wrk1d,wrk2d,wrk3d)

  USE DNS_CONSTANTS, ONLY : efile
  USE DNS_GLOBAL,    ONLY : g
  USE DNS_GLOBAL,    ONLY : imode_sim, inb_scal, imax,jmax,kmax
  USE DNS_GLOBAL,    ONLY : rbg, tbg, sbg, qbg
  USE DNS_GLOBAL,    ONLY : buoyancy
  USE THERMO_GLOBAL, ONLY : imixture

  IMPLICIT NONE

#include "integers.h"

  TREAL, DIMENSION(imax,jmax,kmax),   INTENT(IN)    :: p,T
  TREAL, DIMENSION(imax,jmax,kmax),   INTENT(OUT)   :: rho
  TREAL, DIMENSION(imax,jmax,kmax,*), INTENT(OUT)   :: s
  TREAL, DIMENSION(imax,jmax,kmax),   INTENT(INOUT) :: txc, wrk3d
  TREAL, DIMENSION(jmax,*),           INTENT(INOUT) :: wrk1d,wrk2d

  ! -------------------------------------------------------------------
  TREAL ycenter, dummy
  TINTEGER j, k, is, bcs(2,2)
  TREAL FLOW_SHEAR_TEMPORAL, FLOW_JET_TEMPORAL
  EXTERNAL FLOW_SHEAR_TEMPORAL, FLOW_JET_TEMPORAL

  bcs = 0

  ! -------------------------------------------------------------------
  ! Temporal shear layer case without volumetric force:
  ! Calculate density from equation of state
  ! -------------------------------------------------------------------
  IF ( imode_sim .EQ. DNS_MODE_TEMPORAL ) THEN
    IF ( buoyancy%type .EQ. EQNS_NONE ) THEN

#define TEM_MEAN_LOC(i,j,k) wrk3d(i,j,k)
#define RHO_MEAN_LOC(i,j,k) txc(i,j,k)

      ! temperature/mixture profile are given
      IF ( rbg%type .EQ. PROFILE_NONE ) THEN
        ycenter = g(2)%nodes(1) + g(2)%scale *tbg%ymean
        DO j = 1,jmax
          dummy =  FLOW_SHEAR_TEMPORAL&
          (tbg%type, tbg%thick, tbg%delta, tbg%mean, ycenter, tbg%parameters, g(2)%nodes(j))
          TEM_MEAN_LOC(:,j,:) = dummy
        ENDDO

        DO is = 1,inb_scal
          ycenter = g(2)%nodes(1) + g(2)%scale *sbg(is)%ymean
          DO j = 1,jmax
            dummy =  FLOW_SHEAR_TEMPORAL&
            (sbg(is)%type, sbg(is)%thick, sbg(is)%delta, sbg(is)%mean, ycenter, sbg(is)%parameters, g(2)%nodes(j))
            s(:,j,:,is) = dummy
          ENDDO
        ENDDO

        ! define liquid content in AirWater case: (p,T) given
        IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
          CALL THERMO_AIRWATER_PT(imax,jmax,kmax, s, p, TEM_MEAN_LOC(1,1,1))
        ENDIF

        CALL THERMO_THERMAL_DENSITY(imax,jmax,kmax, s, p, TEM_MEAN_LOC(1,1,1), RHO_MEAN_LOC(1,1,1))
        rho(:,:,:) = rho(:,:,:) + RHO_MEAN_LOC(:,:,:)

        ! density profile itself is given
      ELSE
        ycenter = g(2)%nodes(1) + g(2)%scale*rbg%ymean
        DO j = 1,jmax
          dummy =  FLOW_SHEAR_TEMPORAL&
          (rbg%type, rbg%thick, rbg%delta, rbg%mean, ycenter, rbg%parameters, g(2)%nodes(j))
          rho(:,j,:) = rho(:,j,:) + dummy
        ENDDO

      ENDIF

#undef TEM_MEAN_LOC
#undef RHO_MEAN_LOC

      ! -------------------------------------------------------------------
      ! Temporal shear layer case with volumetric force:
      ! Calculate density from hydrostatic equilibrium
      ! assuming a volumetric force along OY
      ! -------------------------------------------------------------------
    ELSE
      ! AIRWATER case. Routine OPR_PARTIAL_Y introduces small errors in equilibrium
      IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
        CALL THERMO_THERMAL_DENSITY(imax, jmax, kmax, s, p, T, rho)

        ! General case
      ELSE
        CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), p, txc, wrk3d, wrk2d,wrk3d)
        dummy = C_1_R /buoyancy%vector(2)
        rho(:,:,:) = rho(:,:,:) + txc(:,:,:) *dummy
      ENDIF

    ENDIF

    ! -------------------------------------------------------------------
    ! Spatial jet
    ! Only if there is a density variation. Constant density is already
    ! initialized in previous routine segment.
    ! -------------------------------------------------------------------
  ELSE IF ( imode_sim .EQ. DNS_MODE_SPATIAL .AND. rbg%delta .NE. C_0_R ) THEN

    ! temperature/mixture profile are given
    IF ( rbg%type .EQ. PROFILE_NONE ) THEN
#define rho_vi(j) wrk1d(j,1)
#define u_vi(j)   wrk1d(j,2)
#define aux1(j)   wrk1d(j,3)
#define aux2(j)   wrk1d(j,4)
#define aux3(j)   wrk1d(j,5)
#define aux4(j)   wrk1d(j,6)
      ! Inflow profile of density
      ! rho_vi(:) = rho(1,:,1) ! I need to update this because rho(1,:,1) is now undefined

      ! Inflow profile of axial velocity
      ycenter = g(2)%nodes(1) + g(2)%scale*qbg(1)%ymean
      DO j = 1,jmax
        u_vi(j) = FLOW_JET_TEMPORAL&
        (qbg(1)%type, qbg(1)%thick, qbg(1)%delta, qbg(1)%mean, qbg(1)%diam, ycenter, qbg(1)%parameters,g(2)%nodes(j))
      ENDDO

      ! 2D distribution of density
      CALL FLOW_JET_SPATIAL_DENSITY(imax, jmax, &
      tbg%type, tbg%thick, tbg%delta, tbg%mean, tbg%ymean, tbg%diam, tbg%parameters, &
      qbg(1)%type, qbg(1)%thick, qbg(1)%delta, qbg(1)%mean, qbg(1)%ymean, qbg(1)%diam, qbg(1)%parameters, &
      g(2)%scale, g(1)%nodes, g(2)%nodes, s,p,rho_vi(1),u_vi(1),aux1(1),rho,aux2(1),aux3(1),aux4(1))

      DO k = 2,kmax
        rho(:,:,k) = rho(:,:,1)
      ENDDO

      ! density profile itself is given
    ELSE
      ycenter = g(2)%nodes(1) + g(2)%scale*rbg%ymean
      DO j = 1,jmax
        dummy =  FLOW_JET_TEMPORAL&
        (rbg%type, rbg%thick, rbg%delta, rbg%mean, rbg%diam, ycenter, rbg%parameters, g(2)%nodes(j))
        rho(:,j,:) = rho(:,j,:) + dummy
      ENDDO

    ENDIF
  ENDIF

RETURN
END SUBROUTINE DENSITY_MEAN
