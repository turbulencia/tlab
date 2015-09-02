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
!# i1bc      In  grid structure: 0  periodic
!#                               1  non-periodic
!# bcs_jmin  In  BC derivative at jmin: 0  biased, non-zero
!#                                      1  forced to zero
!# bcs_jmax  In  BC derivative at jmax: 0  biased, non-zero
!#                                      1  forced to zero
!#
!########################################################################
SUBROUTINE PARTIAL(&
     imode_fdm, nxy,nz, i1bc, dx, u,up, &
     bcs_jmin,bcs_jmax, wrk1d,wrk2d,wrk3d)

  USE DNS_GLOBAL, ONLY : inb_grid_1

  IMPLICIT NONE

  TINTEGER,                   INTENT(IN) :: imode_fdm, nxy, nz, i1bc, bcs_jmin, bcs_jmax
  TREAL, DIMENSION(nz,*),     INTENT(IN) :: dx
  TREAL, DIMENSION(nxy*nz),      TARGET  :: u, up, wrk3d
  TREAL, DIMENSION(*)                    :: wrk1d ! not used, to be removed
  TREAL, DIMENSION(*)                    :: wrk2d

! -------------------------------------------------------------------
  TINTEGER ip
  TREAL, DIMENSION(:), POINTER :: p_org, p_dst

  p_org => u
  p_dst => up

! ###################################################################
! -------------------------------------------------------------------
! Periodic case
! -------------------------------------------------------------------
  IF ( i1bc .EQ. 0 ) THEN
     IF      ( imode_fdm .EQ. FDM_COM4_JACOBIAN                                 ) THEN; CALL FDM_C1N4P_RHS(nz,nxy, p_org, p_dst)
     ELSE IF ( imode_fdm .EQ. FDM_COM6_JACOBIAN .OR. &
          imode_fdm .EQ. FDM_COM6_DIRECT ) THEN; 
        CALL FDM_C1N6P_RHS(nz,nxy, p_org, p_dst)
     ELSE IF ( imode_fdm .EQ. FDM_COM8_JACOBIAN                                 ) THEN; CALL FDM_C1N8P_RHS(nz,nxy, p_org, p_dst)
     ENDIF
     
     ip  = inb_grid_1 - 1
     CALL TRIDPSS(nz,nxy, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3),dx(1,ip+4),dx(1,ip+5), p_dst,wrk2d)

! -------------------------------------------------------------------
! Nonperiodic case
! -------------------------------------------------------------------
  ELSE
     IF      ( imode_fdm .eq. FDM_COM4_JACOBIAN ) THEN; CALL FDM_C1N4_RHS(nz,nxy, bcs_jmin,bcs_jmax, p_org, p_dst)
     ELSE IF ( imode_fdm .eq. FDM_COM6_JACOBIAN ) THEN; CALL FDM_C1N6_RHS(nz,nxy, bcs_jmin,bcs_jmax, p_org, p_dst)
     ELSE IF ( imode_fdm .eq. FDM_COM8_JACOBIAN ) THEN; CALL FDM_C1N8_RHS(nz,nxy, bcs_jmin,bcs_jmax, p_org, p_dst)
     ELSE IF ( imode_fdm .eq. FDM_COM6_DIRECT   ) THEN; CALL FDM_C1N6_RHS(nz,nxy, bcs_jmin,bcs_jmax, p_org, p_dst) ! not yet implemented
     ENDIF
     
     ip = inb_grid_1 + (bcs_jmin + bcs_jmax*2)*3 - 1
     CALL TRIDSS(nz,nxy, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3), p_dst)

  ENDIF
  
  RETURN
END SUBROUTINE PARTIAL
