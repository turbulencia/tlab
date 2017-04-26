#include "types.h"
#include "dns_const.h"

!########################################################################
!# DESCRIPTION
!#
!########################################################################
!# ARGUMENTS 
!#
!# j1bc      In  grid structure: 0  periodic
!#                               1  non-periodic
!# bcs_jmin  In  BC derivative at jmin: 0  biased, non-zero
!#                                      1  forced to zero
!# bcs_jmax  In  BC derivative at jmax: 0  biased, non-zero
!#                                      1  forced to zero
!#
!########################################################################
SUBROUTINE PARTIAL_Y(imode_fdm, nx,ny,nz, j1bc, dy, u,up, bcs_jmin,bcs_jmax, wrk1d,wrk2d,wrk3d)

  USE DNS_GLOBAL, ONLY : g

  IMPLICIT NONE

  TINTEGER                            :: nx,ny,nz, bcs_jmin,bcs_jmax
  TREAL, DIMENSION(nx*ny*nz), TARGET  :: u, up, wrk3d
  TREAL, DIMENSION(nx*nz)             :: wrk2d
!
  TINTEGER            :: imode_fdm, j1bc      ! not used, to be removed
  TREAL, DIMENSION(*) :: dy, wrk1d            ! not used, to be removed

! -------------------------------------------------------------------
  TINTEGER nxz, nxy, bcs(2)
  TREAL, DIMENSION(:), POINTER :: p_org, p_dst

! ###################################################################
  bcs(1) = bcs_jmin; bcs(2) = bcs_jmax

  IF ( g(2)%size .EQ. 1 ) THEN ! Set to zero in 2D case
     up = C_0_R

  ELSE
! ###################################################################
  nxy = nx*ny 
  nxz = nx*nz

! -------------------------------------------------------------------
! Local transposition: Make y-direction the last one
! -------------------------------------------------------------------
  IF ( nz .EQ. 1 ) THEN
     p_org => u
     p_dst => up
  ELSE
#ifdef USE_ESSL
     CALL DGETMO       (u, nxy, nxy, nz, up, nz)
#else
     CALL DNS_TRANSPOSE(u, nxy, nz, nxy, up, nz)
#endif
     p_org => up
     p_dst => wrk3d
  ENDIF

! ###################################################################
  CALL OPR_PARTIAL1(nxz, bcs, g(2), p_org,p_dst, wrk2d)
  
! ###################################################################
! Put arrays back in the order in which they came in
  IF ( nz .GT. 1 ) THEN
#ifdef USE_ESSL
     CALL DGETMO       (p_dst, nz, nz, nxy, up, nxy)
#else
     CALL DNS_TRANSPOSE(p_dst, nz, nxy, nz, up, nxy)
#endif
  ENDIF

  NULLIFY(p_org,p_dst)

  ENDIF

  RETURN
END SUBROUTINE PARTIAL_Y
