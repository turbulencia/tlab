#include "types.h"

!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2007/11/12 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!########################################################################
!# ARGUMENTS 
!#
!# itype      In    Flag indicating type of filter:
!#                  1 Compact 4th order
!#                  2 Explicit 6th order
!#                  3 Explicit 4th order
!#                  4 Explicit ADM
!# tmp        In    Auxilar 3D array of size only used in ADM type
!#
!########################################################################
SUBROUTINE OPR_FILTER_X(itype, nx, ny, nz, i1bc, bcs_imin, bcs_imax, u, cx, tmp, wrk1d,wrk2d,wrk3d)

  IMPLICIT NONE

  TINTEGER itype, nx, ny, nz, i1bc, bcs_imin, bcs_imax
  TREAL, DIMENSION(*)    :: u, tmp, cx, wrk2d,wrk3d
  TREAL, DIMENSION(nx,*) :: wrk1d
! -----------------------------------------------------------------------
  TINTEGER nyz

! #######################################################################
  nyz = ny*nz

! Make  x  direction the last one 
  CALL DNS_TRANSPOSE(u, nx, nyz, nx, wrk3d, nyz)

! Filter
  IF ( itype .EQ. 1 ) THEN
     CALL FILT4C_KERNEL(nx, nyz, wrk3d, u, i1bc, bcs_imin, bcs_imax, wrk1d, cx)
     IF ( i1bc .EQ. 0 ) THEN
        CALL TRIDPFS(nx,      wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5))
        CALL TRIDPSS(nx, nyz, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5), u, wrk2d)
     ELSE
        CALL TRIDFS(nx,      wrk1d(1,1),wrk1d(1,2),wrk1d(1,3))
        CALL TRIDSS(nx, nyz, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3), u)
     ENDIF

  ELSE IF ( itype .EQ. 2 ) THEN
     CALL FILT6E_KERNEL(nyz, nx, i1bc, bcs_imin, bcs_imax, wrk3d, u)

  ELSE IF ( itype .EQ. 3 ) THEN
     CALL FILT4E_KERNEL(nx, nyz, i1bc, wrk3d, u, cx)

  ELSE IF ( itype .EQ. 4 ) THEN
     CALL FILTADM_KERNEL(nx, nyz, i1bc, wrk3d, u, tmp, cx)

  ENDIF

! Make  x  direction the first one  
  CALL DNS_TRANSPOSE(u, nyz, nx, nyz, wrk3d, nx)

  u(1:nx*nyz) = wrk3d(1:nx*nyz)

  RETURN
END SUBROUTINE OPR_FILTER_X
