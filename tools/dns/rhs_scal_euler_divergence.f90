!########################################################################
!# Tool/Library DNS
!#
!########################################################################
!# HISTORY
!#
!# 2007/06/08 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
SUBROUTINE RHS_SCAL_EULER_DIVERGENCE&
     (dx, dy, dz, rho, u, v, w, z1, zh1, tmp1, tmp2, tmp3, tmp4, wrk1d, wrk2d, wrk3d)

#include "types.h"

  USE DNS_GLOBAL

  IMPLICIT NONE

#include "integers.h"

  TREAL dx(imax)
  TREAL dy(jmax)
  TREAL dz(kmax_total)

  TREAL rho(*), u(*), v(*), w(*), z1(*)
  TREAL zh1(*)
  TREAL tmp1(*), tmp2(*), tmp3(*), tmp4(*)
  TREAL wrk1d(*), wrk2d(*), wrk3d(*)

! -------------------------------------------------------------------
  TINTEGER i
  TREAL dummy

! ###################################################################
  DO i = 1,imax*jmax*kmax
     dummy   = rho(i)*z1(i)
     tmp3(i) = dummy*w(i)
     tmp2(i) = dummy*v(i)
     tmp1(i) = dummy*u(i)
  ENDDO
  CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
       dz, tmp3, tmp4, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
       dy, tmp2, tmp3, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
       dx, tmp1, tmp2, i0, i0, wrk1d, wrk2d, wrk3d)
  DO i = 1,imax*jmax*kmax
     zh1(i) = zh1(i) - ( tmp2(i) + tmp3(i) + tmp4(i) )
  ENDDO

  RETURN
END SUBROUTINE RHS_SCAL_EULER_DIVERGENCE
