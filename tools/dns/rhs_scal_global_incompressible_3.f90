#include "types.h"
#include "dns_const.h"
#include "dns_error.h"

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
!# Scalar equation, nonlinear term in divergence form and the 
!# diffusion term explicit. 3 2nd order + 3 1st order derivatives.
!#
!########################################################################
SUBROUTINE  RHS_SCAL_GLOBAL_INCOMPRESSIBLE_3&
     (is, u,v,w,z1,zh1, tmp1,tmp2,tmp3,tmp4,tmp5,tmp6, wrk1d,wrk2d,wrk3d)

  USE DNS_GLOBAL

  IMPLICIT NONE

#include "integers.h"

  TINTEGER is

  TREAL u(*), v(*), w(*), z1(*), zh1(*)
  TREAL tmp1(*), tmp2(*), tmp3(*), tmp4(*), tmp5(*), tmp6(*)
  TREAL wrk1d(*), wrk2d(*), wrk3d(*)

! -----------------------------------------------------------------------
  TINTEGER ij, i, k
  TREAL diff

! Pointers to existing allocated space
  TREAL, DIMENSION(:), POINTER :: dx,dy,dz

! #######################################################################
! Define pointers
  dx => g(1)%jac(:,1)
  dy => g(2)%jac(:,1)
  dz => g(3)%jac(:,1)

  IF ( idiffusion .EQ. EQNS_NONE ) THEN; diff = C_0_R
  ELSE;                                  diff = visc/schmidt(is); ENDIF

! #######################################################################
! Diffusion and convection terms in scalar equations
! #######################################################################
  CALL PARTIAL_ZZ(i0, iunifz, imode_fdm, imax, jmax, kmax,  k1bc,&
       dz, z1, tmp6, i0, i0, i0, i0, tmp3, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_YY(i0, iunify, imode_fdm, imax, jmax, kmax, j1bc,&
       dy, z1, tmp5, i0, i0, i0, i0, tmp2, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_XX(i0, iunifx, imode_fdm, imax, jmax, kmax, i1bc,&
       dx, z1, tmp4, i0, i0, i0, i0, tmp1, wrk1d, wrk2d, wrk3d)
  DO ij = 1,imax*jmax*kmax
     zh1(ij) = zh1(ij) + diff*( tmp6(ij)+tmp5(ij)+tmp4(ij) )
  ENDDO

  DO ij = 1,imax*jmax*kmax
     tmp6(ij)=z1(ij)*w(ij)
     tmp5(ij)=z1(ij)*v(ij)
     tmp4(ij)=z1(ij)*u(ij)
  ENDDO
  CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
       dz, tmp6, tmp3, i0, i0, wrk1d, wrk2d, wrk3d)
! which BCs should I use here ?
  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
       dy, tmp5, tmp2, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
       dx, tmp4, tmp1, i0, i0, wrk1d, wrk2d, wrk3d)
  DO ij = 1,imax*jmax*kmax
     zh1(ij) = zh1(ij) - ( tmp3(ij) + tmp2(ij) + tmp1(ij) )
  ENDDO

! -----------------------------------------------------------------------
! Dilatation term
! -----------------------------------------------------------------------
  CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
       dz, w, tmp3, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
       dy, v, tmp2, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
       dx, u, tmp1, i0, i0, wrk1d, wrk2d, wrk3d)
!  DO ij = 1,imax*jmax*kmax
!     zh1(ij) = zh1(ij) + z1(ij)*( tmp3(ij) + tmp2(ij) + tmp1(ij) )
!  ENDDO
  DO k = 1,kmax
     DO i = 1,imax
        ij = i                 + imax*jmax*(k-1) ! bottom
        zh1(ij) = zh1(ij) + z1(ij)*( tmp3(ij) + tmp2(ij) + tmp1(ij) )
        ij = i + imax*(jmax-1) + imax*jmax*(k-1) ! top
        zh1(ij) = zh1(ij) + z1(ij)*( tmp3(ij) + tmp2(ij) + tmp1(ij) )
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE RHS_SCAL_GLOBAL_INCOMPRESSIBLE_3
