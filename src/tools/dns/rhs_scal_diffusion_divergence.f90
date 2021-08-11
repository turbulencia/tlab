#include "types.h"
#include "dns_const.h"
#include "dns_error.h"

!########################################################################
!# DESCRIPTION
!#
!# Computing terms depending on diffusion velocities, i.e. diffusion
!# term in the species eqns. and enthapy transport in energy eqn. The
!# latter is stored cumulatively in arrays diff to be employed later
!# in the routine RHS_FLOW_CONDUCTION.
!#
!########################################################################
SUBROUTINE RHS_SCAL_DIFFUSION_DIVERGENCE&
     (is, vis, z1, T, zh1, diff_x,diff_y,diff_z, tmp1,tmp2,tmp3,tmp4, wrk2d,wrk3d)
#ifdef TRACE_ON 
  USE TLAB_CONSTANTS, ONLY : tfile 
#endif 
  USE TLAB_VARS,    ONLY : imax,jmax,kmax, isize_field
  USE TLAB_VARS,    ONLY : g
  USE TLAB_VARS,    ONLY : idiffusion, visc, prandtl, schmidt
  USE THERMO_VARS, ONLY : imixture, THERMO_AI, THERMO_TLIM, NSP, NCP_CHEMKIN
  USE BOUNDARY_BCS

  IMPLICIT NONE

  TINTEGER is
  TREAL, DIMENSION(isize_field),   INTENT(IN)    :: vis, T
  TREAL, DIMENSION(isize_field,*), INTENT(IN)    :: z1
  TREAL, DIMENSION(isize_field),   INTENT(OUT)   :: zh1
  TREAL, DIMENSION(isize_field),   INTENT(OUT)   :: diff_x, diff_y, diff_z
  TREAL, DIMENSION(isize_field),   INTENT(INOUT) :: tmp1,tmp2,tmp3,tmp4
  TREAL, DIMENSION(*),             INTENT(INOUT) :: wrk2d,wrk3d

! -------------------------------------------------------------------
  TINTEGER bcs(2,1)

  TINTEGER i, im, icp
  TREAL diff, cond, dummy
  TREAL ENTHALPY_L, ENTHALPY_G

! ###################################################################
#ifdef TRACE_ON
  CALL TLAB_WRITE_ASCII(tfile, 'ENTERING RHS_SCAL_DIFFUSION_DIVERGENCE')
#endif

  bcs = 0

  IF ( idiffusion .EQ. EQNS_NONE ) THEN; diff = C_0_R;            cond = C_0_R
  ELSE;                                  diff = visc/schmidt(is); cond = visc/prandtl; ENDIF

! -------------------------------------------------------------------
! mass fraction gradients
! -------------------------------------------------------------------
! diffusion velocities in special cases
  IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
     tmp4 = (z1(:,1)-z1(:,2))/(C_1_R-z1(:,2))
     CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), tmp4, tmp1, wrk3d, wrk2d,wrk3d)
     CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), tmp4, tmp2, wrk3d, wrk2d,wrk3d)
     CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), tmp4, tmp3, wrk3d, wrk2d,wrk3d)

! standard diffusion velocities
  ELSE
     CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), z1, tmp1, wrk3d, wrk2d,wrk3d)
     CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), z1, tmp2, wrk3d, wrk2d,wrk3d)
     CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), z1, tmp3, wrk3d, wrk2d,wrk3d)

  ENDIF

! -------------------------------------------------------------------
! enthalpy transport by diffusion velocities
! -------------------------------------------------------------------
! enthalpy difference (h_is-h_NSP) is calculated in array tmp4
  IF ( imixture .GT. 0 .AND. is .LT. NSP .AND. schmidt(is) .NE. prandtl ) THEN
     tmp4 = C_0_R

     DO i = 1, imax*jmax*kmax
        IF ( T(i) .LT. THERMO_TLIM(3,is) ) THEN
           im = 2
        ELSE
           im = 1
        ENDIF
        DO icp = NCP_CHEMKIN,1,-1
           tmp4(i) = tmp4(i)*T(i) + (THERMO_AI(icp,im,is)-THERMO_AI(icp,im,NSP))/M_REAL(icp)
        ENDDO
! factor (diff-cond) added now
        tmp4(i) = (diff-cond)*( tmp4(i)*T(i) + THERMO_AI(6,im,is)-THERMO_AI(6,im,NSP) )
     ENDDO
     diff_x = diff_x + tmp4 *tmp1
     diff_y = diff_y + tmp4 *tmp2
     diff_z = diff_z + tmp4 *tmp3

  ENDIF

! -------------------------------------------------------------------
! mass fraction transport by diffusion velocities
! -------------------------------------------------------------------
  tmp1 = diff *vis *tmp1
  tmp2 = diff *vis *tmp2
  tmp3 = diff *vis *tmp3
  
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs_out(1,2,3), g(3), tmp3, tmp4, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs_out(1,2,2), g(2), tmp2, tmp3, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs_out(1,2,1), g(1), tmp1, tmp2, wrk3d, wrk2d,wrk3d)
  zh1 = zh1 + tmp2 + tmp3 + tmp4

! ###################################################################
! Numerical liquid correction term in the AIRWATER case
! ###################################################################
  IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
     IF ( idiffusion .EQ. EQNS_NONE ) THEN; diff = C_0_R
     ELSE;                                  diff = visc/schmidt(3); ENDIF

! gradient of liquid content
     CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), z1(1,2), tmp1, wrk3d, wrk2d,wrk3d)
     CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), z1(1,2), tmp2, wrk3d, wrk2d,wrk3d)
     CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), z1(1,2), tmp3, wrk3d, wrk2d,wrk3d)

! enthalpy equation
     DO i = 1,imax*jmax*kmax
        ENTHALPY_G = &
             ( (C_1_R  -z1(i,1))*(THERMO_AI(6,1,2)+THERMO_AI(1,1,2)*T(i))&
             + (z1(i,1)-z1(i,2))*(THERMO_AI(6,1,1)+THERMO_AI(1,1,1)*T(i))&
             ) / (C_1_R-z1(i,2))
        ENTHALPY_L = &
             THERMO_AI(6,1,3)+THERMO_AI(1,1,3)*T(i)
        diff_x(i) = diff_x(i) + diff *(ENTHALPY_L-ENTHALPY_G) *tmp1(i)
        diff_y(i) = diff_y(i) + diff *(ENTHALPY_L-ENTHALPY_G) *tmp2(i)
        diff_z(i) = diff_z(i) + diff *(ENTHALPY_L-ENTHALPY_G) *tmp3(i)
     ENDDO

! scalar equation
     DO i = 1,imax*jmax*kmax
        dummy = (C_1_R-z1(i,1))/(C_1_R-z1(i,2))
        tmp1(i) = diff *vis(i) *tmp1(i) *dummy
        tmp2(i) = diff *vis(i) *tmp2(i) *dummy
        tmp3(i) = diff *vis(i) *tmp3(i) *dummy
     ENDDO

     CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs_out(1,2,3), g(3), tmp3, tmp4, wrk3d, wrk2d,wrk3d)
     CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs_out(1,2,2), g(2), tmp2, tmp3, wrk3d, wrk2d,wrk3d)
     CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs_out(1,2,1), g(1), tmp1, tmp2, wrk3d, wrk2d,wrk3d)
     zh1 = zh1 + tmp2 + tmp3 + tmp4

  ENDIF

#ifdef TRACE_ON
  CALL TLAB_WRITE_ASCII(tfile, 'LEAVING RHS_SCAL_DIFFUSION_DIVERGENCE')
#endif

  RETURN
END SUBROUTINE RHS_SCAL_DIFFUSION_DIVERGENCE
