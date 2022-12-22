#include "types.h"
#include "dns_const.h"
#include "dns_error.h"

!########################################################################
!# DESCRIPTION
!#
!# Computing terms depending on diffusion velocities, i.e. diffusion
!# term in the species eqns. and enthalpy transport in energy eqn.
!#
!# Using 2nd order derivative finite difference operators
!#
!########################################################################
SUBROUTINE RHS_SCAL_DIFFUSION_EXPLICIT(is, vis, z1, T, zh1, h4, tmp1,tmp2,tmp3,tmp4,tmp5,tmp6, wrk2d,wrk3d)

  USE TLAB_CONSTANTS, ONLY : efile
#ifdef TRACE_ON
  USE TLAB_CONSTANTS, ONLY : tfile
  USE TLAB_PROCS,     ONLY : TLAB_WRITE_ASCII 
#endif
  USE TLAB_VARS,    ONLY : imax,jmax,kmax, isize_field
  USE TLAB_VARS,    ONLY : g
  USE TLAB_VARS,    ONLY : idiffusion, visc,prandtl,schmidt
  USE THERMO_VARS, ONLY : imixture, THERMO_AI, THERMO_TLIM, NSP, NCP
  use OPR_PARTIAL

  IMPLICIT NONE

  TINTEGER is
  TREAL, DIMENSION(isize_field),   INTENT(IN)    :: vis, T
  TREAL, DIMENSION(isize_field,*), INTENT(IN)    :: z1
  TREAL, DIMENSION(isize_field),   INTENT(OUT)   :: zh1, h4
  TREAL, DIMENSION(isize_field),   INTENT(INOUT) :: tmp1,tmp2,tmp3,tmp4,tmp5,tmp6
  TREAL, DIMENSION(*),             INTENT(INOUT) :: wrk2d,wrk3d

! -------------------------------------------------------------------
  TINTEGER bcs(2,2), i, im, icp
  TREAL diff, cond

! ###################################################################
#ifdef TRACE_ON
  CALL TLAB_WRITE_ASCII(tfile, 'ENTERING RHS_SCAL_DIFFUSION_EXPLICIT')
#endif

  bcs = 0

  IF ( idiffusion .EQ. EQNS_NONE ) THEN; diff = C_0_R;            cond = C_0_R
  ELSE;                                  diff = visc/schmidt(is); cond = visc/prandtl; ENDIF

! ###################################################################
  CALL OPR_PARTIAL_X(OPR_P2, imax,jmax,kmax, bcs, g(1), z1, tmp1, tmp4, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P2, imax,jmax,kmax, bcs, g(2), z1, tmp2, tmp4, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P2, imax,jmax,kmax, bcs, g(3), z1, tmp3, tmp4, wrk2d,wrk3d)
  zh1 = zh1 + diff*vis*( tmp1 + tmp2 + tmp3 )

! -------------------------------------------------------------------
! enthalpy transport by diffusion velocities
! -------------------------------------------------------------------
! enthalpy difference (h_is-h_NSP) is calculated in array tmp4
  IF ( imixture .GT. 0 .AND. is .LT. NSP .AND. schmidt(is) .NE. prandtl ) THEN
     tmp4 = C_0_R

     DO i = 1, imax*jmax*kmax
        IF ( T(i) .LT. THERMO_TLIM(3,is) ) THEN; im = 2
        ELSE;                                    im = 1; ENDIF
        DO icp = NCP,1,-1
           tmp4(i) = tmp4(i)*T(i) +( THERMO_AI(icp,im,is) -THERMO_AI(icp,im,NSP) )/M_REAL(icp)
        ENDDO
! factor (diff-cond) added now
        tmp4(i) = (diff-cond)*( tmp4(i)*T(i) + THERMO_AI(6,im,is)-THERMO_AI(6,im,NSP) )
     ENDDO
     h4 = h4 + vis*tmp4*( tmp1 + tmp2 + tmp3 )

! cross-gradients
     CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), tmp4, tmp1, wrk3d, wrk2d,wrk3d)
     CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), tmp4, tmp2, wrk3d, wrk2d,wrk3d)
     CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), tmp4, tmp3, wrk3d, wrk2d,wrk3d)
     CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), z1,   tmp4, wrk3d, wrk2d,wrk3d)
     CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), z1,   tmp5, wrk3d, wrk2d,wrk3d)
     CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), z1,   tmp6, wrk3d, wrk2d,wrk3d)
     h4 = h4 + vis*( tmp1*tmp4 + tmp2*tmp5 + tmp3*tmp6 )

  ENDIF

#ifdef TRACE_ON
  CALL TLAB_WRITE_ASCII(tfile, 'LEAVING RHS_SCAL_DIFFUSION_EXPLICIT')
#endif

  RETURN
END SUBROUTINE RHS_SCAL_DIFFUSION_EXPLICIT
