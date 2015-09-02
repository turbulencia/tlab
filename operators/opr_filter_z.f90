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
SUBROUTINE OPR_FILTER_Z(itype, nxy,nz, k1bc,bcs_kmin,bcs_kmax, u, cz, tmp, wrk1d,wrk2d,wrk3d)
  
  IMPLICIT NONE
  
  TINTEGER itype, nxy,nz, k1bc, bcs_kmin, bcs_kmax
  TREAL, DIMENSION(*)    :: u, tmp, cz, wrk2d,wrk3d
  TREAL, DIMENSION(nz,*) :: wrk1d
  
! -----------------------------------------------------------------------
  
! #######################################################################

! Filter
  IF ( itype .EQ. 1 ) THEN
     CALL FILT4C_KERNEL(nz, nxy, u, wrk3d, k1bc, bcs_kmin, bcs_kmax, wrk1d, cz)
     IF ( k1bc .EQ. 0 ) THEN
        CALL TRIDPFS(nz,     wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5))
        CALL TRIDPSS(nz,nxy, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5), wrk3d, wrk2d)
     ELSE
        CALL TRIDFS(nz,      wrk1d(1,1),wrk1d(1,2),wrk1d(1,3))
        CALL TRIDSS(nz, nxy, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3), wrk3d)
     ENDIF

  ELSE IF ( itype .EQ. 2 ) THEN
     CALL FILT6E_KERNEL(nxy, nz, k1bc, bcs_kmin, bcs_kmax, u, wrk3d)

  ELSE IF ( itype .EQ. 3 ) THEN
     CALL FILT4E_KERNEL(nz, nxy, k1bc, u, wrk3d, cz)

  ELSE IF ( itype .EQ. 4 ) THEN
     CALL FILTADM_KERNEL(nz, nxy, k1bc, u, wrk3d, tmp, cz)

  ENDIF
  
  u(1:nxy*nz) = wrk3d(1:nxy*nz)
  
  RETURN
END SUBROUTINE OPR_FILTER_Z
