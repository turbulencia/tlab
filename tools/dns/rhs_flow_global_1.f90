#include "types.h"
#include "dns_const.h"
#include "dns_error.h"

!########################################################################
!# Tool/Library DNS
!#
!########################################################################
!# HISTORY
!#
!# 2007/09/03 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Implementation of skewsymmetric formulation and viscous explicit
!# in one routine, to avoid duplication of derivatives. Internal energy
!# formulation only.
!# 30 first derivative operations and 9 second derivative operations.
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
SUBROUTINE RHS_FLOW_GLOBAL_1(dx,dy,dz, vis, rho,u,v,w,p,e, z1, h0,h1,h2,h3,h4, zh1,&
     tmp1,tmp2,tmp3,tmp4,tmp5,tmp6, wrk1d,wrk2d,wrk3d)

  USE DNS_CONSTANTS
  USE DNS_GLOBAL
  USE THERMO_GLOBAL, ONLY : gama0
  USE DNS_LOCAL

  IMPLICIT NONE

#include "integers.h"

  TREAL dx(imax)
  TREAL dy(jmax)
  TREAL dz(kmax_total)

  TREAL vis(*), rho(*), u(*), v(*), w(*), p(*), e(*), z1(imax*jmax*kmax,*)
  TREAL h0(*), h1(*), h2(*), h3(*), h4(*), zh1(imax*jmax*kmax,*)
  TREAL tmp1(*), tmp2(*), tmp3(*), tmp4(*), tmp5(*), tmp6(*)
  TREAL wrk1d(*), wrk2d(*), wrk3d(*)

! -------------------------------------------------------------------
  TINTEGER i1vsin, i1vsout, imxvsin, imxvsout
  TINTEGER j1vsin, j1vsout, jmxvsin, jmxvsout
  TINTEGER k1vsin, k1vsout, kmxvsin, kmxvsout
  TINTEGER i, is
  TREAL g1, g2, g3, prefactor, dummy, c13, dum1, dum2, dum3

! ###################################################################
#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'ENTERING RHS_FLOW_GLOBAL_1')
#endif

  IF ( imode_eqns .EQ. DNS_EQNS_TOTAL ) THEN
     CALL IO_WRITE_ASCII(efile,'RHS_FLOW_GLOBAL_1. No total energy formulation.')
     CALL DNS_STOP(DNS_ERROR_UNDEVELOP)
  ENDIF

  g1 = buoyancy%vector(1)
  g2 = buoyancy%vector(2)
  g3 = buoyancy%vector(3)
  prefactor = (gama0-C_1_R)*mach*mach
  c13 = C_1_R/C_3_R

#include "dns_bcs_inf.h"
#include "dns_bcs_out.h"

! ###################################################################
! Terms \rho u in mass, u-momentum equations
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
       dx, tmp4, tmp6, i0, i0, wrk1d, wrk2d, wrk3d)
  DO i = 1,imax*jmax*kmax
     h0(i) = h0(i) - C_2_R*tmp6(i)
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
! Terms \rho v in mass, v-momentum equations
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
     tmp6(i) = tmp6(i) + tmp5(i)
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
! Terms \rho w in mass, w-momentum equations
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
     tmp6(i) = tmp6(i) + tmp5(i)
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
! Term \rho e in energy equation
! ###################################################################
  DO i = 1,imax*jmax*kmax
     dummy   = C_05_R*rho(i)*e(i)
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

! ###################################################################
! Additional convective part due to skewsymmetric formulation
! Terms u_i d(\rho u_k)/dx_k
! ###################################################################
  DO i = 1,imax*jmax*kmax
     h1(i) = h1(i) - u(i)*tmp6(i)
     h2(i) = h2(i) - v(i)*tmp6(i)
     h3(i) = h3(i) - w(i)*tmp6(i)
     h4(i) = h4(i) - e(i)*tmp6(i)
  ENDDO
  DO is = 1, inb_scal
     DO i = 1,imax*jmax*kmax
        zh1(i,is) = zh1(i,is) - z1(i,is)*tmp6(i)
     ENDDO
  ENDDO

! ###################################################################
! Additional convective part due to skewsymmetric formulation
! Terms \rho u_k du_i/dx_k
! First derivatives in internal energy equation from dissipation
! ###################################################################
! -------------------------------------------------------------------
! Cross derivatives 
! -------------------------------------------------------------------
! energy equation
  CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
       dz, e, tmp4, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
       dy, e, tmp3, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
       dx, e, tmp2, i0, i0, wrk1d, wrk2d, wrk3d)
  DO i = 1,imax*jmax*kmax
     h4(i) = h4(i) - C_05_R*rho(i)*( u(i)*tmp2(i) + v(i)*tmp3(i) + w(i)*tmp4(i) )
  ENDDO

! momentum equations
  dummy = prefactor*visc

  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
       dy, u, tmp1, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
       dx, v, tmp2, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
       dz, v, tmp3, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
       dy, w, tmp4, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
       dz, u, tmp5, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
       dx, w, tmp6, i0, i0, wrk1d, wrk2d, wrk3d)
  DO i = 1,imax*jmax*kmax
     h1(i) = h1(i) - C_05_R*rho(i)*(v(i)*tmp1(i)+w(i)*tmp5(i))
     h2(i) = h2(i) - C_05_R*rho(i)*(u(i)*tmp2(i)+w(i)*tmp3(i))
     h3(i) = h3(i) - C_05_R*rho(i)*(u(i)*tmp6(i)+v(i)*tmp4(i))
     dum1  = tmp1(i)+tmp2(i)
     dum2  = tmp3(i)+tmp4(i)
     dum3  = tmp5(i)+tmp6(i)
     h4(i) = h4(i) + dummy*vis(i)*( dum1*dum1 + dum2*dum2 + dum3*dum3 )
  ENDDO

! -------------------------------------------------------------------
! Dilatation part
! -------------------------------------------------------------------
  CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
       dz, w, tmp4, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
       dy, v, tmp3, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
       dx, u, tmp2, i0, i0, wrk1d, wrk2d, wrk3d)
  dummy = C_2_R*visc
  DO i = 1,imax*jmax*kmax
     h1(i) = h1(i) - C_05_R*rho(i)*u(i)*tmp2(i)
     h2(i) = h2(i) - C_05_R*rho(i)*v(i)*tmp3(i)
     h3(i) = h3(i) - C_05_R*rho(i)*w(i)*tmp4(i)
     tmp1(i) = tmp2(i) + tmp3(i) + tmp4(i)
! the second implementation of this term is much faster in the SUN machine
!     h4(i) = h4(i) + prefactor*( 
! $        c23*visc*vis(i)*( 
! $        (tmp2(i)-tmp3(i))*(tmp2(i)-tmp3(i))+
! $        (tmp3(i)-tmp4(i))*(tmp3(i)-tmp4(i))+
! $        (tmp4(i)-tmp2(i))*(tmp4(i)-tmp2(i)) ) - p(i)*tmp1(i) )
     h4(i) = h4(i) + prefactor*( &
          dummy*vis(i)*( &
          tmp4(i)*tmp4(i) + tmp3(i)*tmp3(i) + tmp2(i)*tmp2(i) - c13*tmp1(i)*tmp1(i) )&
          - p(i)*tmp1(i) )
  ENDDO

! ###################################################################
! Second derivative terms in the momentum equation
! ###################################################################
  CALL PARTIAL_ZZ(i0, iunifz, imode_fdm, imax, jmax, kmax, k1bc,&
       dz, u, tmp5, i0,i0, k1vsout, kmxvsout, tmp6, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_YY(i0, iunify, imode_fdm, imax, jmax, kmax, j1bc,&
       dy, u, tmp4, i0,i0, j1vsout, jmxvsout, tmp6, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_XX(i0, iunifx, imode_fdm, imax, jmax, kmax, i1bc,&
       dx, u, tmp3, i0,i0, i1vsin, imxvsin, tmp6, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
       dx, tmp1, tmp2, i1vsin, imxvsin, wrk1d, wrk2d, wrk3d)
  DO i = 1,imax*jmax*kmax
     h1(i) = h1(i) + vis(i)*visc*(tmp5(i) + tmp4(i) + tmp3(i) + c13*tmp2(i))
  ENDDO

  CALL PARTIAL_ZZ(i0, iunifz, imode_fdm, imax, jmax, kmax, k1bc,&
       dz, v, tmp5, i0,i0, k1vsout, kmxvsout, tmp6, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_YY(i0, iunify, imode_fdm, imax, jmax, kmax, j1bc,&
       dy, v, tmp4, i0,i0, j1vsin, jmxvsin, tmp6, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_XX(i0, iunifx, imode_fdm, imax, jmax, kmax, i1bc,&
       dx, v, tmp3, i0,i0, i1vsout, imxvsout, tmp6, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
       dy, tmp1, tmp2, j1vsin, jmxvsin, wrk1d, wrk2d, wrk3d)
  DO i = 1,imax*jmax*kmax
     h2(i) = h2(i) + vis(i)*visc*(tmp5(i) + tmp4(i) + tmp3(i) + c13*tmp2(i))
  ENDDO

  CALL PARTIAL_ZZ(i0, iunifz, imode_fdm, imax, jmax, kmax, k1bc,&
       dz, w, tmp5, i0,i0, k1vsin, kmxvsin, tmp6, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_YY(i0, iunify, imode_fdm, imax, jmax, kmax, j1bc,&
       dy, w, tmp4, i0,i0, j1vsout, jmxvsout, tmp6, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_XX(i0, iunifx, imode_fdm, imax, jmax, kmax, i1bc,&
       dx, w, tmp3, i0,i0, i1vsout, imxvsout, tmp6, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
       dz, tmp1, tmp2, k1vsin, kmxvsin, wrk1d, wrk2d, wrk3d)
  DO i = 1,imax*jmax*kmax
     h3(i) = h3(i) + vis(i)*visc*(tmp5(i) + tmp4(i) + tmp3(i) + c13*tmp2(i))
  ENDDO

#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'LEAVING RHS_FLOW_GLOBAL_1')
#endif

  RETURN
END SUBROUTINE RHS_FLOW_GLOBAL_1
