#include "types.h"
#include "dns_const.h"

!########################################################################
!# Tool/Library PADE
!#
!########################################################################
!# HISTORY
!#
!# 1999/01/01 - C. Pantano
!#              Created
!# 2003/01/01 - J.P. Mellado
!#              Modified
!#
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

  USE DNS_GLOBAL, ONLY : g !jmax_total, inb_grid_1

  IMPLICIT NONE

  TINTEGER imode_fdm, nx, ny, nz, j1bc, bcs_jmin, bcs_jmax
  TREAL, DIMENSION(*)         :: dy
  TREAL, DIMENSION(nx*ny*nz),    TARGET  :: u, up, wrk3d
  TREAL, DIMENSION(*)                    :: wrk1d ! not used, to be removed
  TREAL, DIMENSION(nx*nz)                :: wrk2d

! -------------------------------------------------------------------
  TINTEGER nxz, nxy
  TREAL, DIMENSION(:), POINTER :: p_org, p_dst

! ###################################################################
  IF ( g(2)%size .EQ. 1 ) THEN ! Set to zero in 2D case
  up = C_0_R

  ELSE
! ###################################################################
  nxy = nx*ny 
  nxz = nx*nz

! -------------------------------------------------------------------
! Make y  direction the last one
! -------------------------------------------------------------------
  IF ( nz .EQ. 1 ) THEN
     p_org => u
     p_dst => up
  ELSE
#ifdef USE_ESSL
     CALL DGETMO(u, nxy, nxy, nz, up, nz)
#else
     CALL DNS_TRANSPOSE(u, nxy, nz, nxy, up, nz)
#endif
     p_org => up
     p_dst => wrk3d
  ENDIF

! ###################################################################
  CALL OPR_PARTIAL1(nxz, g(2), p_org,p_dst, bcs_jmin,bcs_jmax, wrk2d)
  
! ! ###################################################################
! ! -------------------------------------------------------------------
! ! Periodic case
! ! -------------------------------------------------------------------
!   IF ( j1bc .EQ. 0 ) THEN
!      IF      ( imode_fdm .EQ. FDM_COM4_JACOBIAN                                 ) THEN; CALL FDM_C1N4P_RHS(ny,nxz, p_org, p_dst)
!      ELSE IF ( imode_fdm .EQ. FDM_COM6_JACOBIAN .OR. imode_fdm .EQ. FDM_COM6_DIRECT ) THEN; CALL FDM_C1N6P_RHS(ny,nxz, p_org, p_dst)
!      ELSE IF ( imode_fdm .EQ. FDM_COM8_JACOBIAN                                 ) THEN; CALL FDM_C1N8P_RHS(ny,nxz, p_org, p_dst)
!      ENDIF
     
!      ip  = inb_grid_1 - 1
!      CALL TRIDPSS(ny,nxz, dy(1,ip+1),dy(1,ip+2),dy(1,ip+3),dy(1,ip+4),dy(1,ip+5), p_dst,wrk2d)

! ! -------------------------------------------------------------------
! ! Nonperiodic case
! ! -------------------------------------------------------------------
!   ELSE
!      IF      ( imode_fdm .eq. FDM_COM4_JACOBIAN ) THEN; CALL FDM_C1N4_RHS(ny,nxz, bcs_jmin,bcs_jmax, p_org, p_dst)
!      ELSE IF ( imode_fdm .eq. FDM_COM6_JACOBIAN ) THEN; CALL FDM_C1N6_RHS(ny,nxz, bcs_jmin,bcs_jmax, p_org, p_dst)
!      ELSE IF ( imode_fdm .eq. FDM_COM8_JACOBIAN ) THEN; CALL FDM_C1N8_RHS(ny,nxz, bcs_jmin,bcs_jmax, p_org, p_dst)
!      ELSE IF ( imode_fdm .eq. FDM_COM6_DIRECT   ) THEN; CALL FDM_C1N6_RHS(ny,nxz, bcs_jmin,bcs_jmax, p_org, p_dst) ! not yet implemented
!      ENDIF
     
!      ip = inb_grid_1 + (bcs_jmin + bcs_jmax*2)*3 - 1
!      CALL TRIDSS(ny,nxz, dy(1,ip+1),dy(1,ip+2),dy(1,ip+3), p_dst)

!   ENDIF
  
! ###################################################################
! -------------------------------------------------------------------
! Put arrays back in the order in which they came in
! -------------------------------------------------------------------
  IF ( nz .GT. 1 ) THEN
#ifdef USE_ESSL
     CALL DGETMO(p_dst, nz, nz, nxy, up, nxy)
#else
     CALL DNS_TRANSPOSE(p_dst, nz, nxy, nz, up, nxy)
#endif
  ENDIF

  NULLIFY(p_org,p_dst)

  ENDIF

  RETURN
END SUBROUTINE PARTIAL_Y
