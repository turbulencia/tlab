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

  USE DNS_GLOBAL, ONLY : jmax_total, inb_grid_1, inb_grid_2

  IMPLICIT NONE

  TINTEGER ifirst,iunif
  TINTEGER imode_fdm, nx,ny,nz, j1bc, bcs1_jmin,bcs1_jmax, bcs2_jmin,bcs2_jmax
  TREAL, DIMENSION(jmax_total,*)        :: dy
  TREAL, DIMENSION(nx*ny*nz),    TARGET :: u, up1, up2, wrk3d
  TREAL, DIMENSION(*)                   :: wrk1d ! not used, to be removed
  TREAL, DIMENSION(nx*nz)               :: wrk2d

! -------------------------------------------------------------------
  TINTEGER ip, nxy, nxz, ifirst_loc
  TREAL, DIMENSION(:), POINTER :: p_org, p_dst1, p_dst2

! ###################################################################
! Check whether we need to calculate the 1. order derivative
  ifirst_loc = ifirst
  IF ( iunif .NE. 0 ) THEN
     IF ( imode_fdm .eq. FDM_COM4_JACOBIAN .OR. &
          imode_fdm .eq. FDM_COM6_JACOBIAN .OR. &
          imode_fdm .eq. FDM_COM8_JACOBIAN      ) THEN; ifirst_loc = MAX(ifirst_loc,1)
     ENDIF
  ENDIF

  IF ( jmax_total .EQ. 1 ) THEN ! Set to zero in 2D case
  up2 = C_0_R
  IF ( ifirst_loc .EQ. 1 ) up1 = C_0_R
     
  ELSE
! ###################################################################
  nxy = nx*ny 
  nxz = nx*nz

! -------------------------------------------------------------------
! Make y direction the last one
! -------------------------------------------------------------------
  IF ( nz .EQ. 1 ) THEN
     p_org  => u
     p_dst1 => up1
     p_dst2 => up2
  ELSE
#ifdef USE_ESSL
     CALL DGETMO(u, nxy, nxy, nz, up2, nz)
#else
     CALL DNS_TRANSPOSE(u, nxy, nz, nxy, up2, nz)
#endif
     p_org  => up2
     p_dst1 => wrk3d
     p_dst2 => up1
  ENDIF

! ###################################################################
! -------------------------------------------------------------------
! Periodic case
! -------------------------------------------------------------------
  IF ( j1bc  .EQ. 0 ) THEN
     IF ( ifirst_loc .EQ. 1 ) THEN ! First derivative
     IF      ( imode_fdm .eq. FDM_COM4_JACOBIAN                                 ) THEN; CALL FDM_C1N4P_RHS(jmax_total,nxz, p_org, p_dst1)
     ELSE IF ( imode_fdm .eq. FDM_COM6_JACOBIAN .OR. imode_fdm .EQ. FDM_COM6_DIRECT ) THEN; CALL FDM_C1N6P_RHS(jmax_total,nxz, p_org, p_dst1)
     ELSE IF ( imode_fdm .eq. FDM_COM8_JACOBIAN                                 ) THEN; CALL FDM_C1N8P_RHS(jmax_total,nxz, p_org, p_dst1)
     ENDIF
     ip = inb_grid_1 - 1
     CALL TRIDPSS(jmax_total,nxz, dy(1,ip+1),dy(1,ip+2),dy(1,ip+3),dy(1,ip+4),dy(1,ip+5), p_dst1,wrk2d)
     ENDIF

! Calculate second derivative
     IF      ( imode_fdm .eq. FDM_COM4_JACOBIAN                                 ) THEN; CALL FDM_C2N4P_RHS(jmax_total,nxz, p_org, p_dst2)
     ELSE IF ( imode_fdm .eq. FDM_COM6_JACOBIAN .OR. imode_fdm .EQ. FDM_COM6_DIRECT ) THEN; CALL FDM_C2N6P_RHS(jmax_total,nxz, p_org, p_dst2)
     ELSE IF ( imode_fdm .eq. FDM_COM8_JACOBIAN                                 ) THEN; CALL FDM_C2N6P_RHS(jmax_total,nxz, p_org, p_dst2) ! 8th not yet developed
     ENDIF
     ip = inb_grid_2 - 1
     CALL TRIDPSS(jmax_total,nxz, dy(1,ip+1),dy(1,ip+2),dy(1,ip+3),dy(1,ip+4),dy(1,ip+5), p_dst2,wrk2d)

! -------------------------------------------------------------------
! Nonperiodic case
! -------------------------------------------------------------------
  ELSE
     IF ( ifirst_loc .EQ. 1 ) THEN ! First derivative
     IF      ( imode_fdm .eq. FDM_COM4_JACOBIAN ) THEN; CALL FDM_C1N4_RHS(jmax_total,nxz, bcs1_jmin,bcs1_jmax, p_org, p_dst1)
     ELSE IF ( imode_fdm .eq. FDM_COM6_JACOBIAN ) THEN; CALL FDM_C1N6_RHS(jmax_total,nxz, bcs1_jmin,bcs1_jmax, p_org, p_dst1)
     ELSE IF ( imode_fdm .eq. FDM_COM8_JACOBIAN ) THEN; CALL FDM_C1N8_RHS(jmax_total,nxz, bcs1_jmin,bcs1_jmax, p_org, p_dst1)
     ELSE IF ( imode_fdm .eq. FDM_COM6_DIRECT   ) THEN; CALL FDM_C1N6_RHS(jmax_total,nxz, bcs1_jmin,bcs1_jmax, p_org, p_dst1)  ! not yet implemented
     ENDIF
     ip = inb_grid_1 + (bcs1_jmin + bcs1_jmax*2)*3 - 1
     CALL TRIDSS(jmax_total,nxz, dy(1,ip+1),dy(1,ip+2),dy(1,ip+3), p_dst1)
     ENDIF
     
! Calculate second derivative
     IF      ( imode_fdm .eq. FDM_COM4_JACOBIAN ) THEN; CALL FDM_C2N4_RHS(iunif, jmax_total,nxz, bcs2_jmin,bcs2_jmax, dy, p_org, p_dst1, p_dst2)
     ELSE IF ( imode_fdm .eq. FDM_COM6_JACOBIAN ) THEN; CALL FDM_C2N6_RHS(iunif, jmax_total,nxz, bcs2_jmin,bcs2_jmax, dy, p_org, p_dst1, p_dst2)
     ELSE IF ( imode_fdm .eq. FDM_COM8_JACOBIAN ) THEN; CALL FDM_C2N6_RHS(iunif, jmax_total,nxz, bcs2_jmin,bcs2_jmax, dy, p_org, p_dst1, p_dst2)
     ELSE IF ( imode_fdm .eq. FDM_COM6_DIRECT   ) THEN; CALL FDM_C2N6N_RHS(jmax_total,nxz, dy(1,inb_grid_2+3), p_org, p_dst2)
     ENDIF
     ip = inb_grid_2 + (bcs2_jmin + bcs2_jmax*2)*3 - 1
     CALL TRIDSS(jmax_total,nxz, dy(1,ip+1),dy(1,ip+2),dy(1,ip+3), p_dst2)

  ENDIF
     
! ###################################################################
! Put arrays back in the order in which they came in
  IF ( nz .GT. 1 ) THEN
#ifdef USE_ESSL
     CALL DGETMO(p_dst2, nz, nz, nxy, up2, nxy)
#else
     CALL DNS_TRANSPOSE(p_dst2, nz, nxy, nz, up2, nxy)
#endif
     IF ( ifirst .EQ. 1 ) THEN
#ifdef USE_ESSL
        CALL DGETMO(p_dst1, nz, nz, nxy, up1, nxy)
#else
        CALL DNS_TRANSPOSE(p_dst1, nz, nxy, nz, up1, nxy)
#endif
     ENDIF
  ENDIF

  ENDIF

  RETURN
END SUBROUTINE PARTIAL_YY
