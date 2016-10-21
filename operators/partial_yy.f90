#include "types.h"
#include "dns_const.h"

!########################################################################
!# Tool/Library PADE
!#
!########################################################################
!# HISTORY
!#
!# 2007/08/30 - J.P. Mellado
!#              Created
!# 2013/01/20 - J.P. Mellado
!#              Introducing direct formulation of non-uniform grid
!#
!########################################################################
!# DESCRIPTION
!#
!########################################################################
!# ARGUMENTS 
!#
!# ifirst     In  flag indicating to compute first derivative
!# j1bc       In  grid structure: 0  periodic
!#                                1  non-periodic
!# bcs2_jmin  In  BC derivative at imin: 0  biased, non-zero
!#                                       1  forced to zero
!# bcs2_jmax  In  BC derivative at imax: 0  biased, non-zero
!#                                       1  forced to zero
!#
!# up2            Out second derivative
!# up1            Out first derivative, if ifirst=1
!#
!########################################################################
SUBROUTINE PARTIAL_YY(ifirst,iunif,imode_fdm, nx,ny,nz, j1bc, dy, u, up2, &
     bcs1_jmin,bcs1_jmax, bcs2_jmin,bcs2_jmax, up1, wrk1d,wrk2d,wrk3d)

  USE DNS_GLOBAL, ONLY : g

  IMPLICIT NONE

  TINTEGER ifirst,iunif
  TINTEGER imode_fdm, nx,ny,nz, j1bc, bcs1_jmin,bcs1_jmax, bcs2_jmin,bcs2_jmax
  TREAL, DIMENSION(nx*ny*nz),    TARGET :: u, up1, up2, wrk3d
  TREAL, DIMENSION(*)                   :: dy, wrk1d ! not used, to be removed
  TREAL, DIMENSION(nx*nz)               :: wrk2d

! -------------------------------------------------------------------
  TINTEGER nxy, nxz, bcs_min(2), bcs_max(2)!, ifirst_loc
  TREAL, DIMENSION(:), POINTER :: p_org, p_dst1, p_dst2

! ###################################################################
  bcs_min(1) = bcs1_jmin; bcs_max(1) = bcs1_jmax
  bcs_min(2) = bcs2_jmin; bcs_max(2) = bcs2_jmax

  IF ( g(2)%size .EQ. 1 ) THEN ! Set to zero in 2D case
     up2 = C_0_R
     IF ( ifirst .EQ. 1 ) up1 = C_0_R
     
  ELSE
! ###################################################################
  nxy = nx*ny 
  nxz = nx*nz

! -------------------------------------------------------------------
! Local transposition: Make y direction the last one
! -------------------------------------------------------------------
  IF ( nz .EQ. 1 ) THEN
     p_org  => u
     p_dst1 => up1
     p_dst2 => up2
  ELSE
#ifdef USE_ESSL
     CALL DGETMO       (u, nxy, nxy, nz, up2, nz)
#else
     CALL DNS_TRANSPOSE(u, nxy, nz, nxy, up2, nz)
#endif
     p_org  => up2
     p_dst1 => wrk3d
     p_dst2 => up1
  ENDIF

! ###################################################################
  CALL OPR_PARTIAL2(nxz, g(2), p_org,p_dst2, bcs_min,bcs_max, wrk2d,p_dst1)
  
! Check whether we need to calculate the 1. order derivative
  IF ( ifirst .EQ. 1 ) THEN
     IF ( g(2)%uniform .OR. imode_fdm .EQ. FDM_COM6_DIRECT ) THEN
        CALL OPR_PARTIAL1(nxz, g(2), p_org,p_dst1, bcs_min(1),bcs_max(1), wrk2d)
     ENDIF
  ENDIF
  
! ###################################################################
! Put arrays back in the order in which they came in
  IF ( nz .GT. 1 ) THEN
#ifdef USE_ESSL
     CALL DGETMO       (p_dst2, nz, nz, nxy, up2, nxy)
#else
     CALL DNS_TRANSPOSE(p_dst2, nz, nxy, nz, up2, nxy)
#endif
     IF ( ifirst .EQ. 1 ) THEN
#ifdef USE_ESSL
        CALL DGETMO       (p_dst1, nz, nz, nxy, up1, nxy)
#else
        CALL DNS_TRANSPOSE(p_dst1, nz, nxy, nz, up1, nxy)
#endif
     ENDIF
  ENDIF

  ENDIF

  RETURN
END SUBROUTINE PARTIAL_YY
