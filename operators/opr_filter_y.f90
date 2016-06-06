#include "types.h"
#include "dns_const.h"

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
SUBROUTINE OPR_FILTER_Y(imode_filter, nx,ny,nz, j1bc,bcs_jmin,bcs_jmax, u, cy, tmp, wrk1d,wrk2d,wrk3d)
  
  USE DNS_GLOBAL, ONLY : jmax_total

  IMPLICIT NONE

  TINTEGER,                       INTENT(IN)            :: imode_filter, nx,ny,nz, j1bc, bcs_jmin,bcs_jmax
  TREAL, DIMENSION(nx*ny*nz),     INTENT(INOUT), TARGET :: u, wrk3d    ! in-place operation
  TREAL, DIMENSION(jmax_total,6), INTENT(IN)            :: cy          ! Filter kernel information
  TREAL, DIMENSION(*),            INTENT(INOUT)         :: wrk2d       ! Aux arrays
  TREAL, DIMENSION(nx*ny*nz),     INTENT(INOUT)         :: tmp         ! Aux array needed in ADM type
  TREAL, DIMENSION(jmax_total,*), INTENT(INOUT)         :: wrk1d
 
! -----------------------------------------------------------------------
  TINTEGER nxy, nxz

  TREAL, DIMENSION(:), POINTER :: p_org, p_dst

! #######################################################################
  nxy = nx*ny 
  nxz = nx*nz

! -------------------------------------------------------------------
! Make y-direction the last one
! -------------------------------------------------------------------
  IF ( nz .EQ. 1 ) THEN
     p_org => u
     p_dst => wrk3d
  ELSE
#ifdef USE_ESSL
     CALL DGETMO       (u, nxy, nxy, nz, wrk3d, nz)
#else
     CALL DNS_TRANSPOSE(u, nxy, nz, nxy, wrk3d, nz)
#endif
     p_org => wrk3d
     p_dst => u
  ENDIF

! ###################################################################
  IF      ( imode_filter .EQ. DNS_FILTER_COMPACT ) THEN; CALL FILT4C_KERNEL(jmax_total,nxz, p_org,p_dst, j1bc, bcs_jmin, bcs_jmax, wrk1d, cy)
     IF ( j1bc .EQ. 0 ) THEN
        CALL TRIDPFS(jmax_total,     wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5))
        CALL TRIDPSS(jmax_total,nxz, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5), p_dst, wrk2d)
     ELSE
        CALL TRIDFS(jmax_total,     wrk1d(1,1),wrk1d(1,2),wrk1d(1,3))
        CALL TRIDSS(jmax_total,nxz, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3), p_dst)
     ENDIF

  ELSE IF ( imode_filter .EQ. DNS_FILTER_6E      ) THEN; CALL FILT6E_KERNEL (nxz,jmax_total, j1bc, bcs_jmin, bcs_jmax, p_org, p_dst)
  ELSE IF ( imode_filter .EQ. DNS_FILTER_4E      ) THEN; CALL FILT4E_KERNEL (jmax_total,nxz, j1bc,                     p_org, p_dst, cy)
  ELSE IF ( imode_filter .EQ. DNS_FILTER_ADM     ) THEN; CALL FILTADM_KERNEL(jmax_total,nxz, j1bc,                     p_org, p_dst, tmp, cy)
  ENDIF

! -------------------------------------------------------------------
! Put arrays back in the order in which they came in
! -------------------------------------------------------------------
  IF ( nz .GT. 1 ) THEN
#ifdef USE_ESSL
     CALL DGETMO       (p_dst, nz, nz, nxy, p_org, nxy)
#else
     CALL DNS_TRANSPOSE(p_dst, nz, nxy, nz, p_org, nxy)
#endif
  ENDIF

  u(1:nx*ny*nz) = wrk3d(1:nx*ny*nz)

  NULLIFY(p_org,p_dst)

  RETURN
END SUBROUTINE OPR_FILTER_Y
