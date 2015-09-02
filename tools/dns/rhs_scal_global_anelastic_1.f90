!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2007/01/01 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Scalar equation, nonlinear term in convective form and the 
!# diffusion term explicit. 3 2nd order + 3 1st order derivatives.
!#
!########################################################################
!# ARGUMENTS 
!#
!# wrk1d     In     Reference density profile
!#
!########################################################################
#include "types.h"
#include "dns_const.h"
#include "dns_error.h"

SUBROUTINE  RHS_SCAL_GLOBAL_ANELASTIC_1&
     (is, dx,dy,dz, u,v,w, s,hs, tmp1,tmp2,tmp3,tmp4,tmp5,tmp6, wrk1d,wrk2d,wrk3d)

  USE DNS_GLOBAL

  IMPLICIT NONE

#include "integers.h"

  TINTEGER is
  TREAL, DIMENSION(*)              :: dx,dy,dz
  TREAL, DIMENSION(imax,jmax,kmax) :: u,v,w, s, hs
  TREAL, DIMENSION(imax,jmax,kmax) :: tmp1,tmp2,tmp3,tmp4,tmp5,tmp6, wrk2d,wrk3d
  TREAL wrk1d(isize_wrk1d,*)

! -----------------------------------------------------------------------
  TINTEGER i,j,k
  TREAL diff

! #######################################################################
  IF ( idiffusion .EQ. EQNS_NONE ) THEN; diff = C_0_R
  ELSE;                                  diff = visc/schmidt(is); ENDIF

! #######################################################################
! Diffusion and convection terms in scalar equations
! No need to impose dZ/dy=0 because v is already zero
! #######################################################################
  CALL PARTIAL_ZZ(i1, iunifz, imode_fdm, imax, jmax, kmax, k1bc,&
       dz, s, tmp6, i0, i0, i0, i0, tmp3, wrk1d(1,2), wrk2d, wrk3d)
  CALL PARTIAL_YY(i1, iunify, imode_fdm, imax, jmax, kmax, j1bc,&
       dy, s, tmp5, i0, i0, i0, i0, tmp2, wrk1d(1,2), wrk2d, wrk3d)
  CALL PARTIAL_XX(i1, iunifx, imode_fdm, imax, jmax, kmax, i1bc,&
       dx, s, tmp4, i0, i0, i0, i0, tmp1, wrk1d(1,2), wrk2d, wrk3d)
  DO k = 1,kmax
     DO j = 1,jmax
        DO i = 1,imax
           hs(i,j,k) = hs(i,j,k) + visc*( tmp6(i,j,k)+tmp5(i,j,k)+tmp4(i,j,k) ) &
                - wrk1d(j,1)*( w(i,j,k)*tmp3(i,j,k)+v(i,j,k)*tmp2(i,j,k)+u(i,j,k)*tmp1(i,j,k) )
        ENDDO
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE RHS_SCAL_GLOBAL_ANELASTIC_1
