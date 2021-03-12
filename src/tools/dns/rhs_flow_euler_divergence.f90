#include "types.h"
#include "dns_const.h"

!########################################################################
!# DESCRIPTION
!#
!# 15 first derivative operations.
!#
!########################################################################
SUBROUTINE RHS_FLOW_EULER_DIVERGENCE(rho,u,v,w,p,e, h0,h1,h2,h3,h4, tmp1,tmp2,tmp3,tmp4,tmp5, wrk2d,wrk3d)

  USE DNS_GLOBAL,    ONLY : imax,jmax,kmax, isize_field, imode_eqns
  USE DNS_GLOBAL,    ONLY : g, buoyancy
  USE DNS_GLOBAL,    ONLY : mach
  USE THERMO_GLOBAL, ONLY : gama0

  IMPLICIT NONE

  TREAL, DIMENSION(isize_field), INTENT(IN)    :: rho,u,v,w,p,e
  TREAL, DIMENSION(isize_field), INTENT(OUT)   :: h0,h1,h2,h3,h4
  TREAL, DIMENSION(isize_field), INTENT(INOUT) :: tmp1,tmp2,tmp3,tmp4,tmp5
  TREAL, DIMENSION(*),           INTENT(INOUT) :: wrk2d,wrk3d

! -------------------------------------------------------------------
  TINTEGER bcs(2,1), i
  TREAL g1, g2, g3, prefactor, dummy

! ###################################################################
  bcs = 0
  
  g1 = buoyancy%vector(1)
  g2 = buoyancy%vector(2)
  g3 = buoyancy%vector(3)

! ###################################################################
! Terms \rho u in mass and u-momentum equations
! ###################################################################
  DO i = 1,imax*jmax*kmax
     tmp4(i) = rho(i)*u(i)
     tmp3(i) = tmp4(i)*w(i)
     tmp2(i) = tmp4(i)*v(i)
     tmp1(i) = tmp4(i)*u(i) + p(i)
  ENDDO

! -------------------------------------------------------------------
! mass
! -------------------------------------------------------------------
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), tmp4, tmp5, wrk3d, wrk2d,wrk3d)
  h0 = h0 - tmp5

! -------------------------------------------------------------------
! u-momentum
! -------------------------------------------------------------------
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), tmp3, tmp4, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), tmp2, tmp3, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), tmp1, tmp2, wrk3d, wrk2d,wrk3d)
  h1 = h1 - ( tmp2 + tmp3 + tmp4 ) + g1*rho

! ###################################################################
! Terms \rho v in mass and v-momentum equations
! ###################################################################
  DO i = 1,imax*jmax*kmax
     tmp4(i) = rho(i)*v(i)
     tmp3(i) = tmp4(i)*w(i)
     tmp2(i) = tmp4(i)*v(i) + p(i)
     tmp1(i) = tmp4(i)*u(i)
  ENDDO

! -------------------------------------------------------------------
! mass
! -------------------------------------------------------------------
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), tmp4, tmp5, wrk3d, wrk2d,wrk3d)
  h0 = h0 - tmp5

! -------------------------------------------------------------------
! v-momentum
! -------------------------------------------------------------------
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), tmp3, tmp4, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), tmp2, tmp3, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), tmp1, tmp2, wrk3d, wrk2d,wrk3d)
  h2 = h2 - ( tmp2 + tmp3 + tmp4 ) + g2*rho

! ###################################################################
! Terms \rho w in mass and w-momentum equations
! ###################################################################
  DO i = 1,imax*jmax*kmax
     tmp4(i) = rho(i)*w(i)
     tmp3(i) = tmp4(i)*w(i) + p(i)
     tmp2(i) = tmp4(i)*v(i)
     tmp1(i) = tmp4(i)*u(i)
  ENDDO

! -------------------------------------------------------------------
! mass
! -------------------------------------------------------------------
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), tmp4, tmp5, wrk3d, wrk2d,wrk3d)
  h0 = h0 - tmp5

! -------------------------------------------------------------------
! w-momentum
! -------------------------------------------------------------------
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), tmp3, tmp4, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), tmp2, tmp3, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), tmp1, tmp2, wrk3d, wrk2d,wrk3d)
  h3 = h3 - ( tmp2 + tmp3 + tmp4 ) + g3*rho

! ###################################################################
! Total enery equation
! ###################################################################
! -------------------------------------------------------------------
! Total energy formulation
! -------------------------------------------------------------------
  IF ( imode_eqns .EQ. DNS_EQNS_TOTAL ) THEN
     prefactor = (gama0-C_1_R)*mach*mach
     DO i = 1,imax*jmax*kmax
        dummy   = rho(i)*( e(i) + prefactor*C_05_R*(u(i)*u(i)+v(i)*v(i)+w(i)*w(i)) ) + prefactor*p(i)
        tmp3(i) = dummy*w(i)
        tmp2(i) = dummy*v(i)
        tmp1(i) = dummy*u(i)
     ENDDO
     CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), tmp3, tmp4, wrk3d, wrk2d,wrk3d)
     CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), tmp2, tmp3, wrk3d, wrk2d,wrk3d)
     CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), tmp1, tmp2, wrk3d, wrk2d,wrk3d)
     h4 = h4 - ( tmp2 + tmp3 + tmp4 ) + prefactor*rho*( g1*u+g2*v+g3*w )

! -------------------------------------------------------------------
! Internal energy formulation
! the term p div u is add in RHS_FLOW_VISCOUS to save derivatives
! -------------------------------------------------------------------
  ELSE IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
     DO i = 1,imax*jmax*kmax
        dummy   = rho(i)*e(i)
        tmp3(i) = dummy*w(i)
        tmp2(i) = dummy*v(i)
        tmp1(i) = dummy*u(i)
     ENDDO
     CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), tmp3, tmp4, wrk3d, wrk2d,wrk3d)
     CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), tmp2, tmp3, wrk3d, wrk2d,wrk3d)
     CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), tmp1, tmp2, wrk3d, wrk2d,wrk3d)
     h4 = h4 - ( tmp2 + tmp3 + tmp4 )

  ENDIF

  RETURN
END SUBROUTINE RHS_FLOW_EULER_DIVERGENCE
