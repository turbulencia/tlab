#include "types.h"
#include "dns_const.h"
#include "dns_error.h"

!########################################################################
!# Tool/Library DNS
!#
!########################################################################
!# HISTORY
!#
!# 2007/09/21 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Computing terms depending on diffusion velocities, i.e. diffusion
!# term in the species eqns. and enthalpy transport in energy eqn. 
!#
!# Using 2nd order derivative finite difference operators
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
SUBROUTINE RHS_SCAL_DIFFUSION_EXPLICIT&
     (is, dx, dy, dz, vis, z1, T, zh1, h4,&
     tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, wrk1d, wrk2d, wrk3d)

  USE DNS_CONSTANTS
  USE THERMO_GLOBAL
  USE DNS_GLOBAL
  USE DNS_LOCAL

  IMPLICIT NONE

#include "integers.h"

  TINTEGER is

  TREAL, DIMENSION(*)                :: dx, dy, dz
  TREAL, DIMENSION(imax*jmax*kmax)   :: T, vis, zh1, h4
  TREAL, DIMENSION(imax*jmax*kmax,*) :: z1
  TREAL, DIMENSION(imax*jmax*kmax)   :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6
  TREAL, DIMENSION(*)                :: wrk1d, wrk2d, wrk3d

! -------------------------------------------------------------------
  TINTEGER i1vsout, imxvsout
  TINTEGER j1vsout, jmxvsout
  TINTEGER k1vsout, kmxvsout

  TINTEGER i, im, icp
  TREAL diff, cond
!  TREAL ENTHALPY_L, ENTHALPY_G

! ###################################################################
#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'ENTERING RHS_SCAL_DIFFUSION_EXPLICIT')
#endif

  IF ( itransport .EQ. EQNS_TRANS_SUTHERLAND .OR. itransport .EQ. EQNS_TRANS_POWERLAW ) THEN
     CALL IO_WRITE_ASCII(efile,'RHS_SCAL_DIFFUSION_EXPLICIT. Only constant viscosity.')
     CALL DNS_STOP(DNS_ERROR_UNDEVELOP)
  ENDIF

#include "dns_bcs_out.h"

  IF ( idiffusion .EQ. EQNS_NONE ) THEN; diff = C_0_R;            cond = C_0_R
  ELSE;                                  diff = visc/schmidt(is); cond = visc/prandtl; ENDIF

! ###################################################################
! -------------------------------------------------------------------
! species equation
! -------------------------------------------------------------------
! diffusion velocities in special cases
!  IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
!     DO i = 1,imax*jmax*kmax
!        tmp4(i) = (z1(i,1)-z1(i,2))/(C_1_R-z1(i,2))
!     ENDDO
!     CALL PARTIAL_XX(imode_fdm, imax, jmax, kmax, i1bc,
! $        dx, tmp4, tmp1, i0, i0, i0, i0, wrk1d, wrk2d, wrk3d)
!     CALL PARTIAL_YY(imode_fdm, imax, jmax, kmax, j1bc,
! $        dy, tmp4, tmp2, i0, i0, i0, i0, wrk1d, wrk2d, wrk3d)
!     CALL PARTIAL_ZZ(imode_fdm, imax, jmax, kmax, k1bc,
! $        dz, tmp4, tmp3, i0, i0, i0, i0, wrk1d, wrk2d, wrk3d)
!
!c     standard diffusion velocities
!  ELSE
  CALL PARTIAL_XX(i0, iunifx, imode_fdm, imax, jmax, kmax, i1bc,&
       dx, z1, tmp1, i0, i0, i0, i0, tmp4, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_YY(i0, iunify, imode_fdm, imax, jmax, kmax, j1bc,&
       dy, z1, tmp2, i0, i0, i0, i0, tmp4, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_ZZ(i0, iunifz, imode_fdm, imax, jmax, kmax, k1bc,&
       dz, z1, tmp3, i0, i0, i0, i0, tmp4, wrk1d, wrk2d, wrk3d)
!  ENDIF

  DO i = 1,imax*jmax*kmax
     zh1(i) = zh1(i) + diff*vis(i)*( tmp1(i) + tmp2(i) + tmp3(i) )
  ENDDO

! -------------------------------------------------------------------
! enthalpy transport by diffusion velocities
! -------------------------------------------------------------------
! enthalpy difference (h_is-h_NSP) is calculated in array tmp4
  IF ( imixture .GT. 0 .AND. is .LT. NSP .AND. schmidt(is) .NE. prandtl ) THEN
     tmp4 = C_0_R

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
     DO i = 1,imax*jmax*kmax
        h4(i) = h4(i) + vis(i)*tmp4(i)*( tmp1(i) + tmp2(i) + tmp3(i) )
     ENDDO

! cross-gradients
     CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
          dx, tmp4, tmp1, i0, i0, wrk1d, wrk2d, wrk3d)
     CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
          dy, tmp4, tmp2, i0, i0, wrk1d, wrk2d, wrk3d)
     CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
          dz, tmp4, tmp3, i0, i0, wrk1d, wrk2d, wrk3d)
     CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
          dx, z1, tmp4, i0, i0, wrk1d, wrk2d, wrk3d)
     CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
          dy, z1, tmp5, i0, i0, wrk1d, wrk2d, wrk3d)
     CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
          dz, z1, tmp6, i0, i0, wrk1d, wrk2d, wrk3d)
     DO i = 1,imax*jmax*kmax
        h4(i) = h4(i) + vis(i)*( tmp1(i)*tmp4(i) + tmp2(i)*tmp5(i) + tmp3(i)*tmp6(i) )
     ENDDO

  ENDIF

! ###################################################################
! Numerical liquid correction term in the AIRWATER case
! The factor q_d/q_g in species equation is droped to avoid 6 more
! derivatives: relative error 10e-2. See report.
! ###################################################################
!  IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
!     diff = diffcoe*visc/schmidt(3)
!
!c     -------------------------------------------------------------------
!c     species equation
!c     -------------------------------------------------------------------
!     CALL PARTIAL_XX(imode_fdm, imax, jmax, kmax, i1bc,
! $        dx, z1(1,2), tmp1, i0, i0, i0, i0, wrk1d, wrk2d, wrk3d)
!     CALL PARTIAL_YY(imode_fdm, imax, jmax, kmax, j1bc,
! $        dy, z1(1,2), tmp2, i0, i0, i0, i0, wrk1d, wrk2d, wrk3d)
!     CALL PARTIAL_ZZ(imode_fdm, imax, jmax, kmax, k1bc,
! $        dz, z1(1,2), tmp3, i0, i0, i0, i0, wrk1d, wrk2d, wrk3d)
!
!     DO i = 1,imax*jmax*kmax
!        zh1(i) = zh1(i) + diff*vis(i)*( tmp1(i) + tmp2(i) + tmp3(i) )
!     ENDDO
!
!c     -------------------------------------------------------------------
!c     enthalpy equation
!c     -------------------------------------------------------------------
!     DO i = 1,imax*jmax*kmax
!        ENTHALPY_G = 
! $           ( (C_1_R  -z1(i,1))*(THERMO_AI(6,1,2)+THERMO_AI(1,1,2)*T(i))
! $           + (z1(i,1)-z1(i,2))*(THERMO_AI(6,1,1)+THERMO_AI(1,1,1)*T(i))
! $           ) / (C_1_R-z1(i,2))
!        ENTHALPY_L = 
! $           THERMO_AI(6,1,3)+THERMO_AI(1,1,3)*T(i)
!
!        tmp4(i) = diff*(ENTHALPY_L-ENTHALPY_G)
!        h4(i) = h4(i) + vis(i)*tmp4(i)*( tmp1(i) + tmp2(i) + tmp3(i) )
!     ENDDO
!
!c     cross terms
!     CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,
! $        dx, tmp4, tmp1, i0, i0, wrk1d, wrk2d, wrk3d)
!     CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,
! $        dy, tmp4, tmp2, i0, i0, wrk1d, wrk2d, wrk3d)
!     CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,
! $        dz, tmp4, tmp3, i0, i0, wrk1d, wrk2d, wrk3d)
!     CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,
! $        dx, z1(1,2), tmp4, i1vsout, imxvsout, wrk1d, wrk2d, wrk3d)
!     CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,
! $        dy, z1(1,2), tmp5, j1vsout, jmxvsout, wrk1d, wrk2d, wrk3d)
!     CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,
! $        dz, z1(1,2), tmp6, k1vsout, kmxvsout, wrk1d, wrk2d, wrk3d)
!     DO i = 1,imax*jmax*kmax
!        h4(i) = h4(i) + vis(i)*( tmp1(i)*tmp4(i) + tmp2(i)*tmp5(i) + tmp3(i)*tmp6(i) )
!     ENDDO
!
!  ENDIF

#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'LEAVING RHS_SCAL_DIFFUSION_EXPLICIT')
#endif

  RETURN
END SUBROUTINE RHS_SCAL_DIFFUSION_EXPLICIT
