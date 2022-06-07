#include "types.h"
#include "dns_const.h"
#include "dns_error.h"

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
SUBROUTINE RHS_FLOW_GLOBAL_2(rho,u,v,w,p,e,T, z1, h0,h1,h2,h3,h4, zh1, tmp1,tmp2,tmp3,tmp4,tmp5,tmp6, wrk2d,wrk3d)

  USE TLAB_CONSTANTS, ONLY : efile
#ifdef TRACE_ON
  USE TLAB_CONSTANTS, ONLY : tfile
  USE TLAB_PROCS,     ONLY : TLAB_WRITE_ASCII 
#endif
  USE TLAB_VARS,    ONLY : imax,jmax,kmax, isize_field, inb_scal
  USE TLAB_VARS,    ONLY : g, buoyancy
  USE TLAB_VARS,    ONLY : idiffusion, visc,prandtl,mach
  USE THERMO_VARS, ONLY : gama0
  USE BOUNDARY_BCS

#ifdef USE_OPENMP
  USE OMP_LIB
#endif

  IMPLICIT NONE

  TREAL, DIMENSION(isize_field),   INTENT(IN)    :: rho,u,v,w,p,e,T
  TREAL, DIMENSION(isize_field,*), INTENT(IN)    :: z1
  TREAL, DIMENSION(isize_field),   INTENT(OUT)   :: h0,h1,h2,h3,h4
  TREAL, DIMENSION(isize_field,*), INTENT(OUT)   :: zh1
  TREAL, DIMENSION(isize_field),   INTENT(INOUT) :: tmp1,tmp2,tmp3,tmp4,tmp5,tmp6
  TREAL, DIMENSION(*),             INTENT(INOUT) :: wrk2d,wrk3d

! -------------------------------------------------------------------
  TINTEGER bcs(2,1), i, is
  TREAL g1, g2, g3, prefactor, cond, dummy, c13, dum1, dum2, dum3

! ###################################################################
#ifdef TRACE_ON
  CALL TLAB_WRITE_ASCII(tfile, 'ENTERING RHS_FLOW_GLOBAL_2')
#endif

  bcs = 0

  g1 = buoyancy%vector(1)
  g2 = buoyancy%vector(2)
  g3 = buoyancy%vector(3)
  prefactor = (gama0-C_1_R)*mach*mach
  c13 = C_1_R/C_3_R

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
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), tmp4, tmp6, wrk3d, wrk2d,wrk3d)
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
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), tmp3, tmp4, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), tmp2, tmp3, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), tmp1, tmp2, wrk3d, wrk2d,wrk3d)
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
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), tmp4, tmp5, wrk3d, wrk2d,wrk3d)
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
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), tmp3, tmp4, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), tmp2, tmp3, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), tmp1, tmp2, wrk3d, wrk2d,wrk3d)
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
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), tmp4, tmp5, wrk3d, wrk2d,wrk3d)
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
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), tmp3, tmp4, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), tmp2, tmp3, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), tmp1, tmp2, wrk3d, wrk2d,wrk3d)
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
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), tmp3, tmp4, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), tmp2, tmp3, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), tmp1, tmp2, wrk3d, wrk2d,wrk3d)
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
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), e, tmp4, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), e, tmp3, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), e, tmp2, wrk3d, wrk2d,wrk3d)
!$omp parallel default( shared ) private( i )
!$omp do
  DO i = 1,imax*jmax*kmax
     h4(i) = h4(i) - C_05_R*rho(i)*( u(i)*tmp2(i) + v(i)*tmp3(i) + w(i)*tmp4(i) )
  ENDDO
!$omp end do
!$omp end parallel

! momentum equations
  dummy = prefactor*visc

  CALL OPR_PARTIAL_Y(OPR_P2_P1, imax,jmax,kmax, bcs_out(1,1,2), g(2), u, tmp5, tmp1, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P2_P1, imax,jmax,kmax, bcs_out(1,1,3), g(3), u, tmp6, tmp2, wrk2d,wrk3d)
!$omp parallel default( shared ) private( i )
!$omp do
  DO i = 1,imax*jmax*kmax
     h1(i) = h1(i) + visc*(tmp5(i) + tmp6(i)) - C_05_R*rho(i)*(v(i)*tmp1(i)+w(i)*tmp2(i))
  ENDDO
!$omp end do
!$omp end parallel
  CALL OPR_PARTIAL_X(OPR_P2_P1, imax,jmax,kmax, bcs_out(1,1,1), g(1), v, tmp5, tmp3, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P2_P1, imax,jmax,kmax, bcs_out(1,1,3), g(3), v, tmp6, tmp4, wrk2d,wrk3d)
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
  CALL OPR_PARTIAL_X(OPR_P2_P1, imax,jmax,kmax, bcs_out(1,1,1), g(1), w, tmp5, tmp1, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P2_P1, imax,jmax,kmax, bcs_out(1,1,2), g(2), w, tmp6, tmp3, wrk2d,wrk3d)
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
  CALL OPR_PARTIAL_Z(OPR_P2_P1, imax,jmax,kmax, bcs_inf(1,1,3), g(3), w, tmp6, tmp3, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P2_P1, imax,jmax,kmax, bcs_inf(1,1,2), g(2), v, tmp5, tmp2, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P2_P1, imax,jmax,kmax, bcs_inf(1,1,1), g(1), u, tmp4, tmp1, wrk2d,wrk3d)
  dummy = C_2_R*visc
!$omp parallel default( shared ) private( i, dum1 )
!$omp do
  DO i = 1,imax*jmax*kmax
     h1(i) = h1(i) - C_05_R*rho(i)*u(i)*tmp1(i)
     h2(i) = h2(i) - C_05_R*rho(i)*v(i)*tmp2(i)
     h3(i) = h3(i) - C_05_R*rho(i)*w(i)*tmp3(i)

     dum1    = tmp1(i) + tmp2(i) + tmp3(i)
     h4(i) = h4(i) + prefactor*( &
          dummy*( tmp3(i)*tmp3(i) + tmp2(i)*tmp2(i) + tmp1(i)*tmp1(i) - c13*dum1*dum1 ) - p(i)*dum1 )
! array tmp1 no longer needed
     tmp1(i) = c13*dum1
  ENDDO
!$omp end do
!$omp end parallel

! Second derivative terms in the momentum equation
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs_inf(1,2,1), g(1), tmp1, tmp2, wrk3d, wrk2d,wrk3d)
!$omp parallel default( shared ) private( i )
!$omp do
  DO i = 1,imax*jmax*kmax
     h1(i) = h1(i) + visc*(tmp4(i) + tmp2(i))
  ENDDO
!$omp end do
!$omp end parallel

  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs_inf(1,2,2), g(2), tmp1, tmp2, wrk3d, wrk2d,wrk3d)
!$omp parallel default( shared ) private( i )
!$omp do
  DO i = 1,imax*jmax*kmax
     h2(i) = h2(i) + visc*(tmp5(i) + tmp2(i))
  ENDDO
!$omp end do
!$omp end parallel

  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs_inf(1,2,3), g(3), tmp1, tmp2, wrk3d, wrk2d,wrk3d)
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
  CALL THERMO_CALORIC_ENTHALPY(imax,jmax,kmax, z1, T, tmp4)

! total flux
  CALL OPR_PARTIAL_Z(OPR_P2, imax,jmax,kmax, bcs_out(1,1,3), g(3), tmp4, tmp3, tmp5, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P2, imax,jmax,kmax, bcs_out(1,1,2), g(2), tmp4, tmp2, tmp5, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P2, imax,jmax,kmax, bcs_out(1,1,1), g(1), tmp4, tmp1, tmp5, wrk2d,wrk3d)
!$omp parallel default( shared ) private( i )
!$omp do
  DO i = 1,imax*jmax*kmax
     h4(i) = h4(i) + cond*( tmp1(i) + tmp2(i) + tmp3(i) )
  ENDDO
!$omp end do
!$omp end parallel

#ifdef TRACE_ON
  CALL TLAB_WRITE_ASCII(tfile, 'LEAVING RHS_FLOW_GLOBAL_2')
#endif

  RETURN
END SUBROUTINE RHS_FLOW_GLOBAL_2
