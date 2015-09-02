#include "types.h"
#include "dns_const.h"
#include "dns_error.h"

!########################################################################
!# Tool/Library DNS
!#
!########################################################################
!# HISTORY
!#
!# 2008/01/17 - J.P. Mellado
!#              Created
!# 2008/04/11 - J.P. Mellado
!#              Array vis removed
!#
!########################################################################
!# DESCRIPTION
!#
!# Modified from RHS_SCAL_EULER_SKEWSYMMETRIC to include diffusion terms
!# from RHS_SCAL_DIFFUSION_EXPLICIT and avoid duplication of derivatives 
!# in routines PARTIAL_XX, PARTIAL_YY, PARTIAL_ZZ.
!# Internal energy formulation only.
!# Additional convective part due to skewsymmetric formulation Y_i d(\rho u_k)/dx_k
!# done in RHS_FLOW_GLOBAL_2
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
SUBROUTINE RHS_SCAL_GLOBAL_2(is, dx,dy,dz, rho,u,v,w, z1, T, zh1, h4,&
     tmp1,tmp2,tmp3,tmp4,tmp5,tmp6, wrk1d,wrk2d,wrk3d)

#ifdef USE_OPENMP
  USE OMP_LIB
#endif
  USE DNS_CONSTANTS
  USE THERMO_GLOBAL
  USE DNS_GLOBAL
  USE DNS_LOCAL

  IMPLICIT NONE

#include "integers.h"

  TINTEGER is

  TREAL, DIMENSION(*)              :: dx, dy, dz
  TREAL, DIMENSION(imax*jmax*kmax) :: T, rho, u, v, w, z1, zh1, h4
  TREAL, DIMENSION(imax*jmax*kmax) :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6
  TREAL, DIMENSION(*)              :: wrk1d, wrk2d, wrk3d

! -------------------------------------------------------------------
  TINTEGER i1vsout, imxvsout
  TINTEGER j1vsout, jmxvsout
  TINTEGER k1vsout, kmxvsout
  TINTEGER i, im, icp
  TREAL diff, cond, dummy

! ###################################################################
#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'ENTERING RHS_SCAL_GLOBAL_2')
#endif

  IF ( itransport .EQ. EQNS_TRANS_POWERLAW ) THEN
     CALL IO_WRITE_ASCII(efile,'RHS_SCAL_GLOBAL_2. Only constant viscosity.')
     CALL DNS_STOP(DNS_ERROR_UNDEVELOP)
  ENDIF

  IF ( imode_eqns .EQ. DNS_EQNS_TOTAL ) THEN
     CALL IO_WRITE_ASCII(efile,'RHS_SCAL_GLOBAL_2. No total energy formulation.')
     CALL DNS_STOP(DNS_ERROR_UNDEVELOP)
  ENDIF

#include "dns_bcs_out.h"

  IF ( idiffusion .EQ. EQNS_NONE ) THEN; diff = C_0_R;            cond = C_0_R
  ELSE;                                  diff = visc/schmidt(is); cond = visc/prandtl; ENDIF

! ###################################################################
! divergence terms 
! ###################################################################
!$omp parallel default( shared ) private( i, dummy )
!$omp do
  DO i = 1,imax*jmax*kmax
     dummy   = C_05_R*rho(i)*z1(i)
     tmp3(i) = dummy*w(i)
     tmp2(i) = dummy*v(i)
     tmp1(i) = dummy*u(i)
  ENDDO
!$omp end do
!$omp end parallel
  CALL PARTIAL_Z(imode_fdm, imax,jmax,kmax, k1bc, dz, tmp3,tmp4, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, tmp2,tmp3, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc, dx, tmp1,tmp2, i0,i0, wrk1d,wrk2d,wrk3d)
!$omp parallel default( shared ) private( i )
!$omp do
  DO i = 1,imax*jmax*kmax
     zh1(i) = zh1(i) - ( tmp2(i) + tmp3(i) + tmp4(i) )
  ENDDO
!$omp end do
!$omp end parallel

! ###################################################################
! convective part + diffusion
! ###################################################################
  CALL PARTIAL_ZZ(i1, iunifz, imode_fdm, imax,jmax,kmax, k1bc,&
       dz, z1,tmp6, i0,i0, k1vsout,kmxvsout, tmp3, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_YY(i1, iunify, imode_fdm, imax,jmax,kmax, j1bc,&
       dy, z1,tmp5, i0,i0, j1vsout,jmxvsout, tmp2, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_XX(i1, iunifx, imode_fdm, imax,jmax,kmax, i1bc,&
       dx, z1,tmp4, i0,i0, i1vsout,imxvsout, tmp1, wrk1d,wrk2d,wrk3d)

!$omp parallel default( shared ) private( i )
!$omp do
  DO i = 1,imax*jmax*kmax
     zh1(i) = zh1(i) - C_05_R*rho(i)*( u(i)*tmp1(i) + v(i)*tmp2(i) + w(i)*tmp3(i) )&
          + diff*( tmp4(i) + tmp5(i) + tmp6(i) )
  ENDDO
!$omp end do
!$omp end parallel

! -------------------------------------------------------------------
! enthalpy transport by diffusion velocities
! -------------------------------------------------------------------
! enthalpy difference (h_is-h_NSP) is calculated in array tmp4
  IF ( imixture .GT. 0 .AND. is .LT. NSP .AND. schmidt(is) .NE. prandtl ) THEN
     tmp4(:) = C_0_R
     DO i = 1, imax*jmax*kmax
        IF ( T(i) .LT. THERMO_TLIM(3,is) ) THEN; im = 2
        ELSE;                                    im = 1; ENDIF
        DO icp = NCP_CHEMKIN,1,-1
           tmp4(i) = tmp4(i)*T(i) &
                + (THERMO_AI(icp,im,is)-THERMO_AI(icp,im,NSP))/M_REAL(icp)
        ENDDO
! factor (diff-cond) added now
        tmp4(i) = (diff-cond)*( tmp4(i)*T(i) + THERMO_AI(6,im,is)-THERMO_AI(6,im,NSP) )
     ENDDO
!$omp parallel default( shared ) private( i )
!$omp do
     DO i = 1,imax*jmax*kmax
        h4(i) = h4(i) + tmp4(i)*( tmp1(i) + tmp2(i) + tmp3(i) )
     ENDDO
!$omp end do
!$omp end parallel

! cross-gradients
     CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc, dx, tmp4,tmp1, i0,i0, wrk1d,wrk2d,wrk3d)
     CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, tmp4,tmp2, i0,i0, wrk1d,wrk2d,wrk3d)
     CALL PARTIAL_Z(imode_fdm, imax,jmax,kmax, k1bc, dz, tmp4,tmp3, i0,i0, wrk1d,wrk2d,wrk3d)
     CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc, dx, z1,  tmp4, i0,i0, wrk1d,wrk2d,wrk3d)
     CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, z1,  tmp5, i0,i0, wrk1d,wrk2d,wrk3d)
     CALL PARTIAL_Z(imode_fdm, imax,jmax,kmax, k1bc, dz, z1,  tmp6, i0,i0, wrk1d,wrk2d,wrk3d)
!$omp parallel default( shared ) private( i )
!$omp do
     DO i = 1,imax*jmax*kmax
        h4(i) = h4(i) + ( tmp1(i)*tmp4(i) + tmp2(i)*tmp5(i) + tmp3(i)*tmp6(i) )
     ENDDO
!$omp end do
!$omp end parallel

  ENDIF

#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'LEAVING RHS_SCAL_GLOBAL_2')
#endif

  RETURN
END SUBROUTINE RHS_SCAL_GLOBAL_2
