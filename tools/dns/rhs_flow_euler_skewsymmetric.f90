!########################################################################
!# Tool/Library DNS
!#
!########################################################################
!# HISTORY
!#
!# 2007/08/30 - J.P. Mellado
!#              Created.
!#
!########################################################################
!# DESCRIPTION
!#
!# Skewsymmetric formulation according to Erlebacher, 1992.
!# Derived from RHS_FLOW_EULER_DIVERGENCE, 12 additional derivative 
!# operations are added.
!# 27 derivative operations.
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
#include "types.h"
#include "dns_const.h"

SUBROUTINE RHS_FLOW_EULER_SKEWSYMMETRIC(dx,dy,dz, rho,u,v,w,p,e, z1, h0,h1,h2,h3,h4, zh1,&
     tmp1,tmp2,tmp3,tmp4,tmp5, wrk1d,wrk2d,wrk3d)

  USE DNS_GLOBAL
  USE THERMO_GLOBAL, ONLY : gama0

  IMPLICIT NONE

#include "integers.h"

  TREAL dx(imax)
  TREAL dy(jmax)
  TREAL dz(kmax_total)

  TREAL rho(*), u(*), v(*), w(*), p(*), e(*), z1(imax*jmax*kmax,*)
  TREAL h0(*), h1(*), h2(*), h3(*), h4(*), zh1(imax*jmax*kmax,*)
  TREAL tmp1(*), tmp2(*), tmp3(*), tmp4(*), tmp5(*)
  TREAL wrk1d(*), wrk2d(*), wrk3d(*)

! -------------------------------------------------------------------
  TINTEGER i, is
  TREAL g1, g2, g3, prefactor, dummy

! ###################################################################
#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'ENTERING RHS_FLOW_EULER_SKEWSYMMETRIC')
#endif

  g1 = body_vector(1)
  g2 = body_vector(2)
  g3 = body_vector(3)
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
  CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
       dx, tmp4, tmp5, i0, i0, wrk1d, wrk2d, wrk3d)
  DO i = 1,imax*jmax*kmax
     h0(i) = h0(i) - C_2_R*tmp5(i)
  ENDDO

! -------------------------------------------------------------------
! u-momentum
! -------------------------------------------------------------------
  CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
       dz, tmp3, tmp4, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
       dy, tmp2, tmp3, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
       dx, tmp1, tmp2, i0, i0, wrk1d, wrk2d, wrk3d)
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
     DO i = 1,imax*jmax*kmax
        h4(i) = h4(i) - (e(i) + prefactor*C_05_R*(u(i)*u(i)+v(i)*v(i)+w(i)*w(i)))*tmp5(i)
     ENDDO
! Internal energy
  ELSE IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
     DO i = 1,imax*jmax*kmax
        h4(i) = h4(i) - e(i)*tmp5(i)
     ENDDO
  ENDIF

! -------------------------------------------------------------------
! scalar equations
! -------------------------------------------------------------------
  DO is = 1, inb_scal
     DO i = 1,imax*jmax*kmax
        zh1(i,is) = zh1(i,is) - z1(i,is)*tmp5(i)
     ENDDO
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
  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
       dy, tmp4, tmp5, i0, i0, wrk1d, wrk2d, wrk3d)
  DO i = 1,imax*jmax*kmax
     h0(i) = h0(i) - C_2_R*tmp5(i)
  ENDDO

! -------------------------------------------------------------------
! v-momentum
! -------------------------------------------------------------------
  CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
       dz, tmp3, tmp4, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
       dy, tmp2, tmp3, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
       dx, tmp1, tmp2, i0, i0, wrk1d, wrk2d, wrk3d)
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
     DO i = 1,imax*jmax*kmax
        h4(i) = h4(i) - (e(i) + prefactor*C_05_R*(u(i)*u(i)+v(i)*v(i)+w(i)*w(i)))*tmp5(i)
     ENDDO
! Internal energy
  ELSE IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
     DO i = 1,imax*jmax*kmax
        h4(i) = h4(i) - e(i)*tmp5(i)
     ENDDO
  ENDIF

! -------------------------------------------------------------------
! scalar equations
! -------------------------------------------------------------------
  DO is = 1, inb_scal
     DO i = 1,imax*jmax*kmax
        zh1(i,is) = zh1(i,is) - z1(i,is)*tmp5(i)
     ENDDO
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
  CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
       dz, tmp4, tmp5, i0, i0, wrk1d, wrk2d, wrk3d)
  DO i = 1,imax*jmax*kmax
     h0(i) = h0(i) - C_2_R*tmp5(i)
  ENDDO

! -------------------------------------------------------------------
! w-momentum
! -------------------------------------------------------------------
  CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
       dz, tmp3, tmp4, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
       dy, tmp2, tmp3, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
       dx, tmp1, tmp2, i0, i0, wrk1d, wrk2d, wrk3d)
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
     DO i = 1,imax*jmax*kmax
        h4(i) = h4(i) - (e(i) + prefactor*C_05_R*(u(i)*u(i)+v(i)*v(i)+w(i)*w(i)))*tmp5(i)
     ENDDO
! Internal energy
  ELSE IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
     DO i = 1,imax*jmax*kmax
        h4(i) = h4(i) - e(i)*tmp5(i)
     ENDDO
  ENDIF

! -------------------------------------------------------------------
! scalar equations
! -------------------------------------------------------------------
  DO is = 1, inb_scal
     DO i = 1,imax*jmax*kmax
        zh1(i,is) = zh1(i,is) - z1(i,is)*tmp5(i)
     ENDDO
  ENDDO

! ###################################################################
! Energy equation
! ###################################################################
! -------------------------------------------------------------------
! Total energy formulation
! -------------------------------------------------------------------
  IF ( imode_eqns .EQ. DNS_EQNS_TOTAL ) THEN
     DO i = 1,imax*jmax*kmax
        dummy   = C_05_R*rho(i)*( e(i) + prefactor*C_05_R*(u(i)*u(i)+v(i)*v(i)+w(i)*w(i)) ) &
             + prefactor*p(i)
        tmp3(i) = dummy*w(i)
        tmp2(i) = dummy*v(i)
        tmp1(i) = dummy*u(i)
        tmp5(i) = (dummy-prefactor*p(i))/rho(i)
     ENDDO
     CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
          dz, tmp3, tmp4, i0, i0, wrk1d, wrk2d, wrk3d)
     CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
          dy, tmp2, tmp3, i0, i0, wrk1d, wrk2d, wrk3d)
     CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
          dx, tmp1, tmp2, i0, i0, wrk1d, wrk2d, wrk3d)
     DO i = 1,imax*jmax*kmax
        h4(i) = h4(i) - ( tmp2(i) + tmp3(i) + tmp4(i) )&
             + prefactor*rho(i)*( g1*u(i)+g2*v(i)+g3*w(i) )
     ENDDO

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
     CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
          dz, tmp3, tmp4, i0, i0, wrk1d, wrk2d, wrk3d)
     CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
          dy, tmp2, tmp3, i0, i0, wrk1d, wrk2d, wrk3d)
     CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
          dx, tmp1, tmp2, i0, i0, wrk1d, wrk2d, wrk3d)
     DO i = 1,imax*jmax*kmax
        h4(i) = h4(i) - ( tmp2(i) + tmp3(i) + tmp4(i) )
     ENDDO
  ENDIF

! ###################################################################
! Additional convective part due to skewsymmetric formulation
! ###################################################################
! energy equation (array tmp5)
  CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
       dz, tmp5, tmp4, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
       dy, tmp5, tmp3, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
       dx, tmp5, tmp2, i0, i0, wrk1d, wrk2d, wrk3d)
  DO i = 1,imax*jmax*kmax
     h4(i) = h4(i) - rho(i)*( u(i)*tmp2(i) + v(i)*tmp3(i) + w(i)*tmp4(i) )
  ENDDO

! u-momentum equation
  CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
       dz, u, tmp4, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
       dy, u, tmp3, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
       dx, u, tmp2, i0, i0, wrk1d, wrk2d, wrk3d)
  DO i = 1,imax*jmax*kmax
     h1(i) = h1(i) - C_05_R*rho(i)*( u(i)*tmp2(i) + v(i)*tmp3(i) + w(i)*tmp4(i) )
  ENDDO

! v-momentum equation
  CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
       dz, v, tmp4, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
       dy, v, tmp3, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
       dx, v, tmp2, i0, i0, wrk1d, wrk2d, wrk3d)
  DO i = 1,imax*jmax*kmax
     h2(i) = h2(i) - C_05_R*rho(i)*( u(i)*tmp2(i) + v(i)*tmp3(i) + w(i)*tmp4(i) )
  ENDDO

! w-momentum equation
  CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
       dz, w, tmp4, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
       dy, w, tmp3, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
       dx, w, tmp2, i0, i0, wrk1d, wrk2d, wrk3d)
  DO i = 1,imax*jmax*kmax
     h3(i) = h3(i) - C_05_R*rho(i)*( u(i)*tmp2(i) + v(i)*tmp3(i) + w(i)*tmp4(i) )
  ENDDO

#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'LEAVING RHS_FLOW_EULER_SKEWSYMMETRIC')
#endif

  RETURN
END SUBROUTINE RHS_FLOW_EULER_SKEWSYMMETRIC
