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
!# 2008/04/11 - J.P. Mellado
!#              Array vis removed
!#
!########################################################################
!# DESCRIPTION
!#
!# Implementation of skewsymmetric formulation and viscous explicit
!# in one routine, to avoid duplication of derivatives. Internal energy
!# formulation only.
!# 30 first derivative operations and 9 second derivative operations.
!# Viscosity needs to be homogeneous, because of the explicit treatment
!# of diffusion terms.
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
SUBROUTINE RHS_FLOW_GLOBAL_2&
     (dx, dy, dz, rho, u, v, w, p, e, T, z1, h0, h1, h2, h3, h4, zh1,&
     tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, wrk1d, wrk2d, wrk3d)

#ifdef USE_OPENMP
  USE OMP_LIB
#endif
  USE DNS_CONSTANTS
  USE DNS_GLOBAL
  USE THERMO_GLOBAL, ONLY : gama0
  USE DNS_LOCAL

  IMPLICIT NONE

#include "integers.h"

  TREAL dx(imax)
  TREAL dy(jmax)
  TREAL dz(kmax_total)

  TREAL rho(*), u(*), v(*), w(*), p(*), e(*), T(*), z1(imax*jmax*kmax,*)
  TREAL h0(*), h1(*), h2(*), h3(*), h4(*), zh1(imax*jmax*kmax,*)
  TREAL tmp1(*), tmp2(*), tmp3(*), tmp4(*), tmp5(*), tmp6(*)
  TREAL wrk1d(*), wrk2d(*), wrk3d(*)

! -------------------------------------------------------------------
  TINTEGER i1vsin, i1vsout, imxvsin, imxvsout
  TINTEGER j1vsin, j1vsout, jmxvsin, jmxvsout
  TINTEGER k1vsin, k1vsout, kmxvsin, kmxvsout
  TINTEGER i, is
  TREAL g1, g2, g3, prefactor, cond, dummy, c13, dum1, dum2, dum3

! ###################################################################
#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'ENTERING RHS_FLOW_GLOBAL_2')
#endif

  IF ( itransport .EQ. EQNS_TRANS_SUTHERLAND .OR. itransport .EQ. EQNS_TRANS_POWERLAW ) THEN
     CALL IO_WRITE_ASCII(efile,'RHS_FLOW_GLOBAL_2. Only constant viscosity.')
     CALL DNS_STOP(DNS_ERROR_UNDEVELOP)
  ENDIF

  IF ( imode_eqns .EQ. DNS_EQNS_TOTAL ) THEN
     CALL IO_WRITE_ASCII(efile,'RHS_FLOW_GLOBAL_2. No total energy formulation.')
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
!$omp parallel default( shared ) private( i )
!$omp do
  DO i = 1,imax*jmax*kmax
     tmp4(i) = C_05_R*rho(i)*u(i)
     tmp3(i) = tmp4(i)*w(i)
     tmp2(i) = tmp4(i)*v(i)
     tmp1(i) = tmp4(i)*u(i) + p(i)
  ENDDO
!$omp end do
!$omp end parallel

! -------------------------------------------------------------------
! mass
! -------------------------------------------------------------------
  CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc, dx, tmp4, tmp6, i0, i0, wrk1d, wrk2d, wrk3d)
!$omp parallel default( shared ) private( i )
!$omp do
  DO i = 1,imax*jmax*kmax
     h0(i) = h0(i) - C_2_R*tmp6(i)
  ENDDO
!$omp end do
!$omp end parallel

! -------------------------------------------------------------------
! u-momentum
! -------------------------------------------------------------------
  CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc, dz, tmp3, tmp4, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc, dy, tmp2, tmp3, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc, dx, tmp1, tmp2, i0, i0, wrk1d, wrk2d, wrk3d)
!$omp parallel default( shared ) private( i )
!$omp do
  DO i = 1,imax*jmax*kmax
     h1(i) = h1(i) - ( tmp2(i) + tmp3(i) + tmp4(i) ) + g1*rho(i)
  ENDDO
!$omp end do
!$omp end parallel

! ###################################################################
! Terms \rho v in mass, v-momentum equations
! ###################################################################
!$omp parallel default( shared ) private( i )
!$omp do
  DO i = 1,imax*jmax*kmax
     tmp4(i) = C_05_R*rho(i)*v(i)
     tmp3(i) = tmp4(i)*w(i)
     tmp2(i) = tmp4(i)*v(i) + p(i)
     tmp1(i) = tmp4(i)*u(i)
  ENDDO
!$omp end do
!$omp end parallel

! -------------------------------------------------------------------
! mass
! -------------------------------------------------------------------
  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc, dy, tmp4, tmp5, i0, i0, wrk1d, wrk2d, wrk3d)
!$omp parallel default( shared ) private( i )
!$omp do
  DO i = 1,imax*jmax*kmax
     h0(i) = h0(i) - C_2_R*tmp5(i)
     tmp6(i) = tmp6(i) + tmp5(i)
  ENDDO
!$omp end do
!$omp end parallel

! -------------------------------------------------------------------
! v-momentum
! -------------------------------------------------------------------
  CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc, dz, tmp3, tmp4, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc, dy, tmp2, tmp3, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc, dx, tmp1, tmp2, i0, i0, wrk1d, wrk2d, wrk3d)
!$omp parallel default( shared ) private( i )
!$omp do
  DO i = 1,imax*jmax*kmax
     h2(i) = h2(i) - ( tmp2(i) + tmp3(i) + tmp4(i) ) + g2*rho(i)
  ENDDO
!$omp end do
!$omp end parallel

! ###################################################################
! Terms \rho w in mass, w-momentum equations
! ###################################################################
!$omp parallel default( shared ) private( i )
!$omp do
  DO i = 1,imax*jmax*kmax
     tmp4(i) = C_05_R*rho(i)*w(i)
     tmp3(i) = tmp4(i)*w(i) + p(i)
     tmp2(i) = tmp4(i)*v(i)
     tmp1(i) = tmp4(i)*u(i)
  ENDDO
!$omp end do
!$omp end parallel

! -------------------------------------------------------------------
! mass
! -------------------------------------------------------------------
  CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc, dz, tmp4, tmp5, i0, i0, wrk1d, wrk2d, wrk3d)
!$omp parallel default( shared ) private( i )
!$omp do
  DO i = 1,imax*jmax*kmax
     h0(i) = h0(i) - C_2_R*tmp5(i)
     tmp6(i) = tmp6(i) + tmp5(i)
  ENDDO
!$omp end do
!$omp end parallel

! -------------------------------------------------------------------
! w-momentum
! -------------------------------------------------------------------
  CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc, dz, tmp3, tmp4, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc, dy, tmp2, tmp3, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc, dx, tmp1, tmp2, i0, i0, wrk1d, wrk2d, wrk3d)
!$omp parallel default( shared ) private( i )
!$omp do
  DO i = 1,imax*jmax*kmax
     h3(i) = h3(i) - ( tmp2(i) + tmp3(i) + tmp4(i) ) + g3*rho(i)
  ENDDO
!$omp end do
!$omp end parallel

! ###################################################################
! Term \rho e in energy equation
! ###################################################################
!$omp parallel default( shared ) private( i, dummy )
!$omp do
  DO i = 1,imax*jmax*kmax
     dummy   = C_05_R*rho(i)*e(i)
     tmp3(i) = dummy*w(i)
     tmp2(i) = dummy*v(i)
     tmp1(i) = dummy*u(i)
  ENDDO
!$omp end do
!$omp end parallel
  CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc, dz, tmp3, tmp4, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc, dy, tmp2, tmp3, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc, dx, tmp1, tmp2, i0, i0, wrk1d, wrk2d, wrk3d)
!$omp parallel default( shared ) private( i )
!$omp do
  DO i = 1,imax*jmax*kmax
     h4(i) = h4(i) - ( tmp2(i) + tmp3(i) + tmp4(i) )
  ENDDO
!$omp end do
!$omp end parallel
 
! ###################################################################
! Additional convective part due to skewsymmetric formulation
! Terms u_i d(\rho u_k)/dx_k
! ###################################################################
!$omp parallel default( shared ) private( i )
!$omp do
  DO i = 1,imax*jmax*kmax
     h1(i) = h1(i) - u(i)*tmp6(i)
     h2(i) = h2(i) - v(i)*tmp6(i)
     h3(i) = h3(i) - w(i)*tmp6(i)
     h4(i) = h4(i) - e(i)*tmp6(i)
  ENDDO
!$omp end do
!$omp end parallel

  DO is = 1, inb_scal
!$omp parallel default( shared ) private( i )
!$omp do
     DO i = 1,imax*jmax*kmax
        zh1(i,is) = zh1(i,is) - z1(i,is)*tmp6(i)
     ENDDO
!$omp end do
!$omp end parallel
  ENDDO

! ###################################################################
! Additional convective part due to skewsymmetric formulation
! Terms \rho u_k du_i/dx_k
! First derivatives in internal energy equation from dissipation
! Second derivative terms in the momentum equation
! Note that second derivative routines give also first derivatives !
! ###################################################################
! -------------------------------------------------------------------
! Cross derivatives 
! -------------------------------------------------------------------
! energy equation
  CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc, dz, e, tmp4, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc, dy, e, tmp3, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc, dx, e, tmp2, i0, i0, wrk1d, wrk2d, wrk3d)
!$omp parallel default( shared ) private( i )
!$omp do
  DO i = 1,imax*jmax*kmax
     h4(i) = h4(i) - C_05_R*rho(i)*( u(i)*tmp2(i) + v(i)*tmp3(i) + w(i)*tmp4(i) )
  ENDDO
!$omp end do
!$omp end parallel

! momentum equations
  dummy = prefactor*visc

  CALL PARTIAL_YY(i1, iunify, imode_fdm, imax, jmax, kmax, j1bc,&
       dy, u, tmp5, i0,i0, j1vsout, jmxvsout, tmp1, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_ZZ(i1, iunifz, imode_fdm, imax, jmax, kmax, k1bc,&
       dz, u, tmp6, i0,i0, k1vsout, kmxvsout, tmp2, wrk1d, wrk2d, wrk3d)
!$omp parallel default( shared ) private( i )
!$omp do
  DO i = 1,imax*jmax*kmax
     h1(i) = h1(i) + visc*(tmp5(i) + tmp6(i)) - C_05_R*rho(i)*(v(i)*tmp1(i)+w(i)*tmp2(i))
  ENDDO
!$omp end do
!$omp end parallel
  CALL PARTIAL_XX(i1, iunifx, imode_fdm, imax, jmax, kmax, i1bc,&
       dx, v, tmp5, i0,i0, i1vsout, imxvsout, tmp3, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_ZZ(i1, iunifz, imode_fdm, imax, jmax, kmax, k1bc,&
       dz, v, tmp6, i0,i0, k1vsout, kmxvsout, tmp4, wrk1d, wrk2d, wrk3d)
!$omp parallel default( shared ) private( i, dum1 )
!$omp do
  DO i = 1,imax*jmax*kmax
     h2(i) = h2(i) + visc*(tmp5(i) + tmp6(i)) - C_05_R*rho(i)*(u(i)*tmp3(i)+w(i)*tmp4(i))
     dum1  = tmp1(i)+tmp3(i)
     h4(i) = h4(i) + dummy*dum1*dum1
  ENDDO
!$omp end do
!$omp end parallel

! arrays tmp1 and tmp3 can be reused
  CALL PARTIAL_XX(i1, iunifx, imode_fdm, imax, jmax, kmax, i1bc,&
       dx, w, tmp5, i0,i0, i1vsout, imxvsout, tmp1, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_YY(i1, iunify, imode_fdm, imax, jmax, kmax, j1bc,&
       dy, w, tmp6, i0,i0, j1vsout, jmxvsout, tmp3, wrk1d, wrk2d, wrk3d)
!$omp parallel default( shared ) private( i, dum2, dum3 )
!$omp do
  DO i = 1,imax*jmax*kmax
     h3(i) = h3(i) + visc*(tmp5(i) + tmp6(i)) - C_05_R*rho(i)*(u(i)*tmp1(i)+v(i)*tmp3(i))
     dum2  = tmp4(i)+tmp3(i)
     dum3  = tmp1(i)+tmp2(i)
     h4(i) = h4(i) + dummy*( dum2*dum2 + dum3*dum3 )
  ENDDO
!$omp end do
!$omp end parallel

! -------------------------------------------------------------------
! Dilatation part
! -------------------------------------------------------------------
  CALL PARTIAL_ZZ(i1, iunifz, imode_fdm, imax, jmax, kmax, k1bc,&
       dz, w, tmp6, i0,i0, k1vsin, kmxvsin, tmp3, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_YY(i1, iunify, imode_fdm, imax, jmax, kmax, j1bc,&
       dy, v, tmp5, i0,i0, j1vsin, jmxvsin, tmp2, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_XX(i1, iunifx, imode_fdm, imax, jmax, kmax, i1bc,&
       dx, u, tmp4, i0,i0, i1vsin, imxvsin, tmp1, wrk1d, wrk2d, wrk3d)
  dummy = C_2_R*visc
!$omp parallel default( shared ) private( i, dum1 )
!$omp do
  DO i = 1,imax*jmax*kmax
     h1(i) = h1(i) - C_05_R*rho(i)*u(i)*tmp1(i)
     h2(i) = h2(i) - C_05_R*rho(i)*v(i)*tmp2(i)
     h3(i) = h3(i) - C_05_R*rho(i)*w(i)*tmp3(i)

     dum1    = tmp1(i) + tmp2(i) + tmp3(i)
     h4(i) = h4(i) + prefactor*( &
          dummy*( tmp3(i)*tmp3(i) + tmp2(i)*tmp2(i) + tmp1(i)*tmp1(i) - c13*dum1*dum1 )& 
          - p(i)*dum1 )
! array tmp1 no longer needed
     tmp1(i) = c13*dum1
  ENDDO
!$omp end do
!$omp end parallel

! Second derivative terms in the momentum equation
  CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
       dx, tmp1, tmp2, i1vsin, imxvsin, wrk1d, wrk2d, wrk3d)
!$omp parallel default( shared ) private( i )
!$omp do
  DO i = 1,imax*jmax*kmax
     h1(i) = h1(i) + visc*(tmp4(i) + tmp2(i))
  ENDDO
!$omp end do
!$omp end parallel

  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
       dy, tmp1, tmp2, j1vsin, jmxvsin, wrk1d, wrk2d, wrk3d)
!$omp parallel default( shared ) private( i )
!$omp do
  DO i = 1,imax*jmax*kmax
     h2(i) = h2(i) + visc*(tmp5(i) + tmp2(i))
  ENDDO
!$omp end do
!$omp end parallel

  CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
       dz, tmp1, tmp2, k1vsin, kmxvsin, wrk1d, wrk2d, wrk3d)
!$omp parallel default( shared ) private( i )
!$omp do
  DO i = 1,imax*jmax*kmax
     h3(i) = h3(i) + visc*(tmp6(i) + tmp2(i))
  ENDDO
!$omp end do
!$omp end parallel

! ###################################################################
! Enthalpy diffusion in energy equation
! ###################################################################
  IF ( idiffusion .EQ. EQNS_NONE ) THEN; cond = C_0_R
  ELSE;                                  cond = visc/prandtl; ENDIF

! calculate the enthalpy
  CALL THERMO_CALORIC_ENTHALPY(imax, jmax, kmax, z1, T, tmp4)

! total flux
  CALL PARTIAL_ZZ(i0, iunifz, imode_fdm, imax, jmax, kmax, k1bc,&
       dz, tmp4, tmp3, i0,i0, k1vsout, kmxvsout, tmp5, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_YY(i0, iunify, imode_fdm, imax, jmax, kmax, j1bc,&
       dy, tmp4, tmp2, i0,i0, j1vsout, jmxvsout, tmp5, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_XX(i0, iunifx, imode_fdm, imax, jmax, kmax, i1bc,&
       dx, tmp4, tmp1, i0,i0, i1vsout, imxvsout, tmp5, wrk1d, wrk2d, wrk3d)
!$omp parallel default( shared ) private( i )
!$omp do
  DO i = 1,imax*jmax*kmax
     h4(i) = h4(i) + cond*( tmp1(i) + tmp2(i) + tmp3(i) )
  ENDDO
!$omp end do
!$omp end parallel

#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'LEAVING RHS_FLOW_GLOBAL_2')
#endif

  RETURN
END SUBROUTINE RHS_FLOW_GLOBAL_2
