#include "types.h"
#include "dns_const.h"
#include "dns_error.h"

!########################################################################
!# DESCRIPTION
!#
!# Modified from RHS_SCAL_EULER_SKEWSYMMETRIC to include diffusion terms
!# from RHS_SCAL_DIFFUSION_EXPLICIT and avoid duplication of derivatives
!# in routines OPR_PARTIAL_XX, OPR_PARTIAL_YY, OPR_PARTIAL_ZZ.
!# Internal energy formulation only.
!# Additional convective part due to skewsymmetric formulation Y_i d(\rho u_k)/dx_k
!# done in RHS_FLOW_GLOBAL_2
!#
!########################################################################
SUBROUTINE RHS_SCAL_GLOBAL_2(is, rho,u,v,w, z1, T, zh1, h4, tmp1,tmp2,tmp3,tmp4,tmp5,tmp6, wrk2d,wrk3d)

  USE TLAB_CONSTANTS, ONLY : efile
#ifdef TRACE_ON
  USE TLAB_CONSTANTS, ONLY : tfile
#endif
  USE TLAB_VARS,    ONLY : imax,jmax,kmax, isize_field
  USE TLAB_VARS,    ONLY : g
  USE TLAB_VARS,    ONLY : idiffusion, visc,prandtl,schmidt
  USE THERMO_VARS, ONLY : imixture, THERMO_AI, THERMO_TLIM, NSP, NCP_CHEMKIN
  USE BOUNDARY_BCS

#ifdef USE_OPENMP
  USE OMP_LIB
#endif

  IMPLICIT NONE

  TINTEGER is
  TREAL, DIMENSION(isize_field),   INTENT(IN)    :: rho,u,v,w,T, z1
  TREAL, DIMENSION(isize_field),   INTENT(OUT)   :: zh1,h4
  TREAL, DIMENSION(isize_field),   INTENT(INOUT) :: tmp1,tmp2,tmp3,tmp4,tmp5,tmp6
  TREAL, DIMENSION(*),             INTENT(INOUT) :: wrk2d,wrk3d

! -------------------------------------------------------------------
  TINTEGER bcs(2,1), i, im, icp
  TREAL diff, cond, dummy

! ###################################################################
#ifdef TRACE_ON
  CALL TLAB_WRITE_ASCII(tfile, 'ENTERING RHS_SCAL_GLOBAL_2')
#endif

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
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), tmp3,tmp4, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), tmp2,tmp3, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), tmp1,tmp2, wrk3d, wrk2d,wrk3d)
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
  CALL OPR_PARTIAL_Z(OPR_P2_P1, imax,jmax,kmax, bcs_out(1,1,3), g(3), z1,tmp6, tmp3, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P2_P1, imax,jmax,kmax, bcs_out(1,1,2), g(2), z1,tmp5, tmp2, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P2_P1, imax,jmax,kmax, bcs_out(1,1,1), g(1), z1,tmp4, tmp1, wrk2d,wrk3d)

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
     CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), tmp4,tmp1, wrk3d, wrk2d,wrk3d)
     CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), tmp4,tmp2, wrk3d, wrk2d,wrk3d)
     CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), tmp4,tmp3, wrk3d, wrk2d,wrk3d)
     CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), z1,  tmp4, wrk3d, wrk2d,wrk3d)
     CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), z1,  tmp5, wrk3d, wrk2d,wrk3d)
     CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), z1,  tmp6, wrk3d, wrk2d,wrk3d)
!$omp parallel default( shared ) private( i )
!$omp do
     DO i = 1,imax*jmax*kmax
        h4(i) = h4(i) + ( tmp1(i)*tmp4(i) + tmp2(i)*tmp5(i) + tmp3(i)*tmp6(i) )
     ENDDO
!$omp end do
!$omp end parallel

  ENDIF

#ifdef TRACE_ON
  CALL TLAB_WRITE_ASCII(tfile, 'LEAVING RHS_SCAL_GLOBAL_2')
#endif

  RETURN
END SUBROUTINE RHS_SCAL_GLOBAL_2
