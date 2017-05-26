#include "types.h"
#include "dns_const.h"

!########################################################################
!# DESCRIPTION
!#
!# remove divergence part of a vector field a=(u,v,w)
!#
!# Calculate scalar field phi s.t. lap phi = -div a, with BCs phi = 0
!# at top and bottom (i.e. zero tangential component of vector grad phi) 
!# Then, add grad phi to vector a, where n grad phi = 0 at top and bottom. 
!#
!# The BCs are such that a and a + grad phi are the same at top and bottom 
!#
!########################################################################
SUBROUTINE SOLENOIDAL(iwall, u,v,w, tmp1,tmp2,tmp3,tmp4,tmp5, wrk1d,wrk2d,wrk3d)

  USE DNS_GLOBAL, ONLY : g
  USE DNS_GLOBAL, ONLY : imax,jmax,kmax, isize_wrk1d

  IMPLICIT NONE

#include "integers.h"

  TINTEGER iwall

  TREAL, DIMENSION(imax,jmax,kmax) :: u,v,w
  TREAL, DIMENSION(imax,jmax,kmax) :: tmp1,tmp2,tmp3,tmp4,tmp5, wrk3d
  TREAL, DIMENSION(imax,kmax,*)    :: wrk2d
  TREAL, DIMENSION(isize_wrk1d,*)  :: wrk1d

! -------------------------------------------------------------------
  TINTEGER  ibc, bcs(2,2)

! ###################################################################
  bcs = 0
  
!  IF ( iwall .EQ. 1) THEN; ibc = 4
!  ELSE;                    ibc = 3; ENDIF
  ibc = 3

! ###################################################################
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), u, tmp2, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), v, tmp3, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), w, tmp4, wrk3d, wrk2d,wrk3d)
  tmp2 = tmp2 +tmp3 +tmp4

! -------------------------------------------------------------------
! Solve lap(phi) = - div(u)
! -------------------------------------------------------------------
! wx, wy are aux arrays
  IF ( g(1)%periodic .AND. g(3)%periodic ) THEN ! Doubly periodic in xOz
     wrk2d(:,:,1:2) = C_0_R  ! bcs
     tmp1 = -tmp2            ! change of forcing term sign
     CALL OPR_POISSON_FXZ(.FALSE., imax,jmax,kmax, g, ibc, &
          tmp1,wrk3d, tmp4,tmp5, wrk2d(1,1,1),wrk2d(1,1,2), wrk1d,wrk1d(1,5),wrk3d)

  ELSE                                          ! General treatment
#ifdef USE_CGLOC
! Need to define global variable with ipos,jpos,kpos,ci,cj,ck,
     CALL CGPOISSON(i1, imax,jmax,kmax,kmax_total, tmp1, tmp2,tmp3,tmp4, ipos,jpos,kpos,ci,cj,ck, wrk2d)
#endif
  ENDIF

! -------------------------------------------------------------------
! Eliminate solenoidal part of u by adding grad(phi)
! -------------------------------------------------------------------
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), tmp1, tmp2, wrk3d, wrk2d,wrk3d)
  u = u + tmp2
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), tmp1, tmp2, wrk3d, wrk2d,wrk3d)
  v = v + tmp2
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), tmp1, tmp2, wrk3d, wrk2d,wrk3d)     
  w = w + tmp2
  
  RETURN
END SUBROUTINE SOLENOIDAL

