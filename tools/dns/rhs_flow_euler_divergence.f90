#include "types.h"
#include "dns_const.h"

!########################################################################
!# Tool/Library DNS
!#
!########################################################################
!# HISTORY
!#
!# 2007/06/08 - J.P. Mellado
!#              Created
!# 2007/08/16 - J.P. Mellado
!#              Case of internal energy formulation added
!#
!########################################################################
!# DESCRIPTION
!#
!# 15 first derivative operations.
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
SUBROUTINE RHS_FLOW_EULER_DIVERGENCE(dx,dy,dz, rho,u,v,w,p,e, h0,h1,h2,h3,h4,&
     tmp1,tmp2,tmp3,tmp4,tmp5, wrk1d,wrk2d,wrk3d)

  USE DNS_GLOBAL
  USE THERMO_GLOBAL, ONLY : gama0

  IMPLICIT NONE

#include "integers.h"

  TREAL dx(imax)
  TREAL dy(jmax)
  TREAL dz(kmax_total)

  TREAL rho(*), u(*), v(*), w(*), p(*), e(*)
  TREAL h0(*), h1(*), h2(*), h3(*), h4(*)
  TREAL tmp1(*), tmp2(*), tmp3(*), tmp4(*), tmp5(*)
  TREAL wrk1d(*), wrk2d(*), wrk3d(*)

! -------------------------------------------------------------------
  TINTEGER i
  TREAL g1, g2, g3, prefactor, dummy

! ###################################################################
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
  CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
       dx, tmp4, tmp5, i0, i0, wrk1d, wrk2d, wrk3d)
  DO i = 1,imax*jmax*kmax
     h0(i) = h0(i) - tmp5(i)
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
     h1(i) = h1(i) - ( tmp2(i) + tmp3(i) + tmp4(i) ) + g1*rho(i)
  ENDDO

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
  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
       dy, tmp4, tmp5, i0, i0, wrk1d, wrk2d, wrk3d)
  DO i = 1,imax*jmax*kmax
     h0(i) = h0(i) - tmp5(i)
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
     h2(i) = h2(i) - ( tmp2(i) + tmp3(i) + tmp4(i) ) + g2*rho(i)
  ENDDO

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
  CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
       dz, tmp4, tmp5, i0, i0, wrk1d, wrk2d, wrk3d)
  DO i = 1,imax*jmax*kmax
     h0(i) = h0(i) - tmp5(i)
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
     h3(i) = h3(i) - ( tmp2(i) + tmp3(i) + tmp4(i) ) + g3*rho(i)
  ENDDO

! ###################################################################
! Total enery equation
! ###################################################################
! -------------------------------------------------------------------
! Total energy formulation
! -------------------------------------------------------------------
  IF ( imode_eqns .EQ. DNS_EQNS_TOTAL ) THEN
     prefactor = (gama0-C_1_R)*mach*mach
     DO i = 1,imax*jmax*kmax
        dummy   = rho(i)*( e(i) + prefactor*C_05_R*(u(i)*u(i)+v(i)*v(i)+w(i)*w(i)) ) &
             + prefactor*p(i)
        tmp3(i) = dummy*w(i)
        tmp2(i) = dummy*v(i)
        tmp1(i) = dummy*u(i)
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
        dummy   = rho(i)*e(i)
        tmp3(i) = dummy*w(i)
        tmp2(i) = dummy*v(i)
        tmp1(i) = dummy*u(i)
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

  RETURN
END SUBROUTINE RHS_FLOW_EULER_DIVERGENCE
