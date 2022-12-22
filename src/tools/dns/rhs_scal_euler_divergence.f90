#include "types.h"
#include "dns_const.h"

SUBROUTINE RHS_SCAL_EULER_DIVERGENCE(rho,u,v,w, z1, zh1, tmp1,tmp2,tmp3,tmp4, wrk2d,wrk3d)

  USE TLAB_VARS, ONLY : imax,jmax,kmax, isize_field
  USE TLAB_VARS, ONLY : g
  use OPR_PARTIAL

  IMPLICIT NONE

  TREAL, DIMENSION(isize_field), INTENT(IN)    :: rho,u,v,w,z1
  TREAL, DIMENSION(isize_field), INTENT(OUT)   :: zh1
  TREAL, DIMENSION(isize_field), INTENT(INOUT) :: tmp1,tmp2,tmp3,tmp4
  TREAL, DIMENSION(*),           INTENT(INOUT) :: wrk2d,wrk3d

! -------------------------------------------------------------------
  TINTEGER bcs(2,1), i
  TREAL dummy

! ###################################################################
  DO i = 1,imax*jmax*kmax
     dummy   = rho(i)*z1(i)
     tmp3(i) = dummy*w(i)
     tmp2(i) = dummy*v(i)
     tmp1(i) = dummy*u(i)
  ENDDO
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), tmp3, tmp4, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), tmp2, tmp3, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), tmp1, tmp2, wrk3d, wrk2d,wrk3d)
  zh1 = zh1 - ( tmp2 + tmp3 + tmp4 )

  RETURN
END SUBROUTINE RHS_SCAL_EULER_DIVERGENCE
