#include "types.h"
#include "dns_const.h"
#include "dns_error.h"

!########################################################################
!# Tool/Library DNS
!#
!########################################################################
!# HISTORY
!#
!# 2007/06/08 - J.P. Mellado
!#              Created
!# 2007/06/20 - J.P. Mellado
!#              AIRWATER case included. 
!# 2007/07/03 - J.P. Mellado
!#              AIRWATER correction term included.
!#
!########################################################################
!# DESCRIPTION
!#
!# Computing terms depending on diffusion velocities, i.e. diffusion
!# term in the species eqns. and enthapy transport in energy eqn. The
!# latter is stored cumulatively in arrays diff to be employed later
!# in the routine RHS_FLOW_CONDUCTION.
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
SUBROUTINE RHS_SCAL_DIFFUSION_DIVERGENCE&
     (is, dx, dy, dz, vis, z1, T, zh1, diff_x, diff_y, diff_z,&
     tmp1, tmp2, tmp3, tmp4, wrk1d, wrk2d, wrk3d)

  USE THERMO_GLOBAL
  USE DNS_GLOBAL
  USE DNS_LOCAL

  IMPLICIT NONE

#include "integers.h"

  TINTEGER is

  TREAL, DIMENSION(*)                :: dx, dy, dz
  TREAL, DIMENSION(imax*jmax*kmax)   :: T, vis, zh1
  TREAL, DIMENSION(imax*jmax*kmax,*) :: z1
  TREAL, DIMENSION(imax*jmax*kmax)   :: tmp1, tmp2, tmp3, tmp4
  TREAL, DIMENSION(*)                :: wrk1d, wrk2d, wrk3d

  TREAL diff_x(*), diff_y(*), diff_z(*)

! -------------------------------------------------------------------
  TINTEGER i1vsout, imxvsout
  TINTEGER j1vsout, jmxvsout
  TINTEGER k1vsout, kmxvsout

  TINTEGER i, im, icp
  TREAL diff, cond, dummy
  TREAL ENTHALPY_L, ENTHALPY_G

! ###################################################################
#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'ENTERING RHS_SCAL_DIFFUSION_DIVERGENCE')
#endif

#include "dns_bcs_out.h"

  IF ( idiffusion .EQ. EQNS_NONE ) THEN; diff = C_0_R;            cond = C_0_R
  ELSE;                                  diff = visc/schmidt(is); cond = visc/prandtl; ENDIF

! -------------------------------------------------------------------
! mass fraction gradients
! -------------------------------------------------------------------
! diffusion velocities in special cases
  IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
     DO i = 1,imax*jmax*kmax
        tmp4(i) = (z1(i,1)-z1(i,2))/(C_1_R-z1(i,2))
     ENDDO
     CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
          dx, tmp4, tmp1, i0, i0, wrk1d, wrk2d, wrk3d)
     CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
          dy, tmp4, tmp2, i0, i0, wrk1d, wrk2d, wrk3d)
     CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
          dz, tmp4, tmp3, i0, i0, wrk1d, wrk2d, wrk3d)

! standard diffusion velocities
  ELSE
     CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
          dx, z1, tmp1, i0, i0, wrk1d, wrk2d, wrk3d)
     CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
          dy, z1, tmp2, i0, i0, wrk1d, wrk2d, wrk3d)
     CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
          dz, z1, tmp3, i0, i0, wrk1d, wrk2d, wrk3d)
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
           tmp4(i) = tmp4(i)*T(i) &
                + (THERMO_AI(icp,im,is)-THERMO_AI(icp,im,NSP))/M_REAL(icp)
        ENDDO
! factor (diff-cond) added now
        tmp4(i) = (diff-cond)*( tmp4(i)*T(i) + THERMO_AI(6,im,is)-THERMO_AI(6,im,NSP) )
     ENDDO
     DO i = 1,imax*jmax*kmax
        diff_x(i) = diff_x(i) + tmp4(i)*tmp1(i)
        diff_y(i) = diff_y(i) + tmp4(i)*tmp2(i)
        diff_z(i) = diff_z(i) + tmp4(i)*tmp3(i)
     ENDDO

  ENDIF

! -------------------------------------------------------------------
! mass fraction transport by diffusion velocities
! -------------------------------------------------------------------
  DO i = 1,imax*jmax*kmax
     tmp1(i) = diff*vis(i)*tmp1(i)
     tmp2(i) = diff*vis(i)*tmp2(i)
     tmp3(i) = diff*vis(i)*tmp3(i)
  ENDDO

  CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
       dz, tmp3, tmp4, k1vsout, kmxvsout, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
       dy, tmp2, tmp3, j1vsout, jmxvsout, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
       dx, tmp1, tmp2, i1vsout, imxvsout, wrk1d, wrk2d, wrk3d)
  DO i = 1,imax*jmax*kmax
     zh1(i) = zh1(i) + tmp2(i) + tmp3(i) + tmp4(i)
  ENDDO

! ###################################################################
! Numerical liquid correction term in the AIRWATER case
! ###################################################################
  IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
     IF ( idiffusion .EQ. EQNS_NONE ) THEN; diff = C_0_R
     ELSE;                                  diff = visc/schmidt(3); ENDIF

! gradient of liquid content
     CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
          dx, z1(1,2), tmp1, i0, i0, wrk1d, wrk2d, wrk3d)
     CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
          dy, z1(1,2), tmp2, i0, i0, wrk1d, wrk2d, wrk3d)
     CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
          dz, z1(1,2), tmp3, i0, i0, wrk1d, wrk2d, wrk3d)

! enthalpy equation
     DO i = 1,imax*jmax*kmax
        ENTHALPY_G = &
             ( (C_1_R  -z1(i,1))*(THERMO_AI(6,1,2)+THERMO_AI(1,1,2)*T(i))&
             + (z1(i,1)-z1(i,2))*(THERMO_AI(6,1,1)+THERMO_AI(1,1,1)*T(i))&
             ) / (C_1_R-z1(i,2))
        ENTHALPY_L = &
             THERMO_AI(6,1,3)+THERMO_AI(1,1,3)*T(i)
        diff_x(i) = diff_x(i) + diff*(ENTHALPY_L-ENTHALPY_G)*tmp1(i)
        diff_y(i) = diff_y(i) + diff*(ENTHALPY_L-ENTHALPY_G)*tmp2(i)
        diff_z(i) = diff_z(i) + diff*(ENTHALPY_L-ENTHALPY_G)*tmp3(i)
     ENDDO

! scalar equation
     DO i = 1,imax*jmax*kmax
        dummy = (C_1_R-z1(i,1))/(C_1_R-z1(i,2))
        tmp1(i) = diff*vis(i)*tmp1(i)*dummy
        tmp2(i) = diff*vis(i)*tmp2(i)*dummy
        tmp3(i) = diff*vis(i)*tmp3(i)*dummy
     ENDDO

     CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
          dz, tmp3, tmp4, k1vsout, kmxvsout, wrk1d, wrk2d, wrk3d)
     CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
          dy, tmp2, tmp3, j1vsout, jmxvsout, wrk1d, wrk2d, wrk3d)
     CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
          dx, tmp1, tmp2, i1vsout, imxvsout, wrk1d, wrk2d, wrk3d)
     DO i = 1,imax*jmax*kmax
        zh1(i) = zh1(i) + tmp2(i) + tmp3(i) + tmp4(i)
     ENDDO

  ENDIF

#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'LEAVING RHS_SCAL_DIFFUSION_DIVERGENCE')
#endif

  RETURN
END SUBROUTINE RHS_SCAL_DIFFUSION_DIVERGENCE
