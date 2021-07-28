#include "types.h"
#include "dns_const.h"

!########################################################################
!# DESCRIPTION
!#
!# Skewsymmetric formulation according to Erlebacher, 1992.
!# Derived from RHS_FLOW_EULER_DIVERGENCE, 12 additional derivative 
!# operations are added.
!# 27 derivative operations.
!#
!########################################################################
SUBROUTINE RHS_FLOW_EULER_SKEWSYMMETRIC(rho,u,v,w,p,e, z1, h0,h1,h2,h3,h4, zh1,&
     tmp1,tmp2,tmp3,tmp4,tmp5, wrk2d,wrk3d)
#ifdef TRACE_ON 
  USE DNS_CONSTANTS, ONLY : tfile 
#endif 

  USE TLAB_VARS,    ONLY : imax,jmax,kmax, isize_field, inb_scal, imode_eqns
  USE TLAB_VARS,    ONLY : g, buoyancy
  USE TLAB_VARS,    ONLY : mach
  USE THERMO_GLOBAL, ONLY : gama0

  IMPLICIT NONE

  TREAL, DIMENSION(isize_field),   INTENT(IN)    :: rho,u,v,w,p,e
  TREAL, DIMENSION(isize_field,*), INTENT(IN)    :: z1
  TREAL, DIMENSION(isize_field),   INTENT(OUT)   :: h0,h1,h2,h3,h4
  TREAL, DIMENSION(isize_field,*), INTENT(OUT)   :: zh1
  TREAL, DIMENSION(isize_field),   INTENT(INOUT) :: tmp1,tmp2,tmp3,tmp4,tmp5
  TREAL, DIMENSION(*),             INTENT(INOUT) :: wrk2d,wrk3d

! -------------------------------------------------------------------
  TINTEGER bcs(2,1), i, is
  TREAL g1, g2, g3, prefactor, dummy

! ###################################################################
#ifdef TRACE_ON
  CALL TLAB_WRITE_ASCII(tfile, 'ENTERING RHS_FLOW_EULER_SKEWSYMMETRIC')
#endif

  bcs = 0
  
  g1 = buoyancy%vector(1)
  g2 = buoyancy%vector(2)
  g3 = buoyancy%vector(3)
  prefactor = (gama0-C_1_R)*mach*mach

! ###################################################################
! Terms \rho u in mass, u-momentum and energy equations
! ###################################################################
  DO i = 1,imax*jmax*kmax
     tmp4(i) = C_05_R*rho(i)*u(i)
     tmp3(i) = tmp4(i)*w(i)
     tmp2(i) = tmp4(i)*v(i)
     tmp1(i) = tmp4(i)*u(i) + p(i)
  ENDDO

! -------------------------------------------------------------------
! mass
! -------------------------------------------------------------------
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), tmp4, tmp5, wrk3d, wrk2d,wrk3d)
  h0 = h0 - C_2_R*tmp5

! -------------------------------------------------------------------
! u-momentum
! -------------------------------------------------------------------
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), tmp3, tmp4, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), tmp2, tmp3, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), tmp1, tmp2, wrk3d, wrk2d,wrk3d)
  DO i = 1,imax*jmax*kmax
     h1(i) = h1(i) - ( tmp2(i) + tmp3(i) + tmp4(i) + u(i)*tmp5(i) ) + g1*rho(i)
     h2(i) = h2(i) - v(i)*tmp5(i)
     h3(i) = h3(i) - w(i)*tmp5(i)
  ENDDO

! -------------------------------------------------------------------
! energy equation
! -------------------------------------------------------------------
! Total energy
  IF ( imode_eqns .EQ. DNS_EQNS_TOTAL ) THEN
     h4 = h4 - (e + prefactor*C_05_R*(u*u+v*v+w*w))*tmp5
! Internal energy
  ELSE IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
     h4 = h4 - e*tmp5
  ENDIF

! -------------------------------------------------------------------
! scalar equations
! -------------------------------------------------------------------
  DO is = 1, inb_scal
     zh1(:,is) = zh1(:,is) - z1(:,is)*tmp5
  ENDDO

! ###################################################################
! Terms \rho v in mass, v-momentum and energy equations
! ###################################################################
  DO i = 1,imax*jmax*kmax
     tmp4(i) = C_05_R*rho(i)*v(i)
     tmp3(i) = tmp4(i)*w(i)
     tmp2(i) = tmp4(i)*v(i) + p(i)
     tmp1(i) = tmp4(i)*u(i)
  ENDDO

! -------------------------------------------------------------------
! mass
! -------------------------------------------------------------------
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), tmp4, tmp5, wrk3d, wrk2d,wrk3d)
  h0 = h0 - C_2_R*tmp5

! -------------------------------------------------------------------
! v-momentum
! -------------------------------------------------------------------
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), tmp3, tmp4, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), tmp2, tmp3, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), tmp1, tmp2, wrk3d, wrk2d,wrk3d)
  DO i = 1,imax*jmax*kmax
     h1(i) = h1(i) - u(i)*tmp5(i)
     h2(i) = h2(i) - ( tmp2(i) + tmp3(i) + tmp4(i) + v(i)*tmp5(i) ) + g2*rho(i)
     h3(i) = h3(i) - w(i)*tmp5(i)
  ENDDO

! -------------------------------------------------------------------
! energy equation
! -------------------------------------------------------------------
! Total energy
  IF ( imode_eqns .EQ. DNS_EQNS_TOTAL ) THEN
     h4 = h4 - (e + prefactor*C_05_R*(u*u+v*v+w*w))*tmp5
! Internal energy
  ELSE IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
     h4 = h4 - e*tmp5
  ENDIF

! -------------------------------------------------------------------
! scalar equations
! -------------------------------------------------------------------
  DO is = 1, inb_scal
     zh1(:,is) = zh1(:,is) - z1(:,is)*tmp5
  ENDDO

! ###################################################################
! Terms \rho w in mass, w-momentum and energy equations
! ###################################################################
  DO i = 1,imax*jmax*kmax
     tmp4(i) = C_05_R*rho(i)*w(i)
     tmp3(i) = tmp4(i)*w(i) + p(i)
     tmp2(i) = tmp4(i)*v(i)
     tmp1(i) = tmp4(i)*u(i)
  ENDDO

! -------------------------------------------------------------------
! mass
! -------------------------------------------------------------------
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), tmp4, tmp5, wrk3d, wrk2d,wrk3d)
  h0 = h0 - C_2_R*tmp5

! -------------------------------------------------------------------
! w-momentum
! -------------------------------------------------------------------
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), tmp3, tmp4, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), tmp2, tmp3, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), tmp1, tmp2, wrk3d, wrk2d,wrk3d)
  DO i = 1,imax*jmax*kmax
     h1(i) = h1(i) - u(i)*tmp5(i)
     h2(i) = h2(i) - v(i)*tmp5(i)
     h3(i) = h3(i) - ( tmp2(i) + tmp3(i) + tmp4(i) + w(i)*tmp5(i) ) + g3*rho(i)
  ENDDO

! -------------------------------------------------------------------
! energy equation
! -------------------------------------------------------------------
! Total energy
  IF ( imode_eqns .EQ. DNS_EQNS_TOTAL ) THEN
     h4 = h4 - (e + prefactor*C_05_R*(u*u+v*v+w*w))*tmp5
! Internal energy
  ELSE IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
     h4 = h4 - e*tmp5
  ENDIF

! -------------------------------------------------------------------
! scalar equations
! -------------------------------------------------------------------
  DO is = 1, inb_scal
     zh1(:,is) = zh1(:,is) - z1(:,is)*tmp5
  ENDDO

! ###################################################################
! Energy equation
! ###################################################################
! -------------------------------------------------------------------
! Total energy formulation
! -------------------------------------------------------------------
  IF ( imode_eqns .EQ. DNS_EQNS_TOTAL ) THEN
     DO i = 1,imax*jmax*kmax
        dummy   = C_05_R*rho(i)*( e(i) + prefactor*C_05_R*(u(i)*u(i)+v(i)*v(i)+w(i)*w(i)) )              + prefactor*p(i)
        tmp3(i) = dummy*w(i)
        tmp2(i) = dummy*v(i)
        tmp1(i) = dummy*u(i)
        tmp5(i) = (dummy-prefactor*p(i))/rho(i)
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
        dummy   = C_05_R*rho(i)*e(i)
        tmp3(i) = dummy*w(i)
        tmp2(i) = dummy*v(i)
        tmp1(i) = dummy*u(i)
        tmp5(i) = dummy/rho(i)
     ENDDO
     CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), tmp3, tmp4, wrk3d, wrk2d,wrk3d)
     CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), tmp2, tmp3, wrk3d, wrk2d,wrk3d)
     CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), tmp1, tmp2, wrk3d, wrk2d,wrk3d)
     h4 = h4 - ( tmp2 + tmp3 + tmp4 )
  ENDIF

! ###################################################################
! Additional convective part due to skewsymmetric formulation
! ###################################################################
! energy equation (array tmp5)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), tmp5, tmp4, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), tmp5, tmp3, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), tmp5, tmp2, wrk3d, wrk2d,wrk3d)
  h4 = h4 - rho*( u*tmp2 + v*tmp3 + w*tmp4 )

! u-momentum equation
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), u, tmp4, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), u, tmp3, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), u, tmp2, wrk3d, wrk2d,wrk3d)
  h1 = h1 - C_05_R*rho*( u*tmp2 + v*tmp3 + w*tmp4 )

! v-momentum equation
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), v, tmp4, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), v, tmp3, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), v, tmp2, wrk3d, wrk2d,wrk3d)
  h2 = h2 - C_05_R*rho*( u*tmp2 + v*tmp3 + w*tmp4 )

! w-momentum equation
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), w, tmp4, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), w, tmp3, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), w, tmp2, wrk3d, wrk2d,wrk3d)
  h3 = h3 - C_05_R*rho*( u*tmp2 + v*tmp3 + w*tmp4 )

#ifdef TRACE_ON
  CALL TLAB_WRITE_ASCII(tfile, 'LEAVING RHS_FLOW_EULER_SKEWSYMMETRIC')
#endif

  RETURN
END SUBROUTINE RHS_FLOW_EULER_SKEWSYMMETRIC
