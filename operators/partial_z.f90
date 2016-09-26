#include "types.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

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
!# k1bc      In  grid structure: 0  periodic
!#                               1  non-periodic
!# bcs_kmin  In  BC derivative at kmin: 0  biased, non-zero
!#                                      1  forced to zero
!# bcs_kmax  In  BC derivative at kmax: 0  biased, non-zero
!#                                      1  forced to zero
!#
!########################################################################
SUBROUTINE PARTIAL_Z(imode_fdm, nx,ny,nz, k1bc, dz, u,up, bcs_kmin,bcs_kmax, wrk1d,wrk2d,wrk3d)

  USE DNS_GLOBAL, ONLY : g !kmax_total, inb_grid_1
#ifdef USE_MPI
  USE DNS_MPI
#endif

  IMPLICIT NONE

  TINTEGER imode_fdm, nx,ny,nz, k1bc, bcs_kmin, bcs_kmax
  TREAL, DIMENSION(*)        :: dz
  TREAL, DIMENSION(nx*ny*nz),    TARGET :: u, up, wrk3d
  TREAL, DIMENSION(*)                   :: wrk1d ! not used, to be removed
  TREAL, DIMENSION(nx*ny)               :: wrk2d

! -------------------------------------------------------------------
  TINTEGER nxy

  TREAL, DIMENSION(:), POINTER :: p_a, p_b

#ifdef USE_MPI
  TINTEGER id
#endif

! ###################################################################
#ifdef USE_MPI
  id  = DNS_MPI_K_PARTIAL
#endif

  IF ( g(3)%size .EQ. 1 ) THEN ! Set to zero in 2D case
  up = C_0_R
     
  ELSE
! ###################################################################
! -------------------------------------------------------------------
! Transposition
! -------------------------------------------------------------------
#ifdef USE_MPI         
  IF ( ims_npro_k .GT. 1 ) THEN
     CALL DNS_MPI_TRPF_K(u, up, ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))
     p_a => up
     p_b => wrk3d
     nxy = ims_size_k(id)
  ELSE
#endif
     p_a => u
     p_b => up
     nxy = nx*ny 
#ifdef USE_MPI         
  ENDIF
#endif

! ###################################################################
  CALL OPR_PARTIAL(imode_fdm, nxy, g(3), p_a,p_b, bcs_kmin,bcs_kmax, wrk2d)
  
! ! ###################################################################
! ! -------------------------------------------------------------------
! ! Periodic case
! ! -------------------------------------------------------------------
!   IF ( k1bc .EQ. 0 ) THEN
!      IF      ( imode_fdm .eq. FDM_COM4_JACOBIAN                                 ) THEN; CALL FDM_C1N4P_RHS(kmax_total,nxy, p_a, p_b)
!      ELSE IF ( imode_fdm .eq. FDM_COM6_JACOBIAN .OR. imode_fdm .EQ. FDM_COM6_DIRECT ) THEN; CALL FDM_C1N6P_RHS(kmax_total,nxy, p_a, p_b)
!      ELSE IF ( imode_fdm .eq. FDM_COM8_JACOBIAN                                 ) THEN; CALL FDM_C1N8P_RHS(kmax_total,nxy, p_a, p_b)
!      ENDIF
     
!      ip  = inb_grid_1 - 1
!      CALL TRIDPSS(kmax_total,nxy, dz(1,ip+1),dz(1,ip+2),dz(1,ip+3),dz(1,ip+4),dz(1,ip+5), p_b,wrk2d)
     
! ! -------------------------------------------------------------------
! ! Nonperiodic case
! ! -------------------------------------------------------------------
!   ELSE
!      IF      ( imode_fdm .eq. FDM_COM4_JACOBIAN ) THEN; CALL FDM_C1N4_RHS(kmax_total,nxy, bcs_kmin,bcs_kmax, p_a, p_b)
!      ELSE IF ( imode_fdm .eq. FDM_COM6_JACOBIAN ) THEN; CALL FDM_C1N6_RHS(kmax_total,nxy, bcs_kmin,bcs_kmax, p_a, p_b)
!      ELSE IF ( imode_fdm .eq. FDM_COM8_JACOBIAN ) THEN; CALL FDM_C1N8_RHS(kmax_total,nxy, bcs_kmin,bcs_kmax, p_a, p_b)
!      ELSE IF ( imode_fdm .eq. FDM_COM6_DIRECT   ) THEN; CALL FDM_C1N6_RHS(kmax_total,nxy, bcs_kmin,bcs_kmax, p_a, p_b) ! not yet implemented
!      ENDIF

!      ip = inb_grid_1 + (bcs_kmin + bcs_kmax*2)*3 - 1
!      CALL TRIDSS(kmax_total,nxy, dz(1,ip+1),dz(1,ip+2),dz(1,ip+3), p_b)
     
!   ENDIF
     
! ###################################################################
! -------------------------------------------------------------------
! Transposition
! -------------------------------------------------------------------
#ifdef USE_MPI         
  IF ( ims_npro_k .GT. 1 ) THEN
     CALL DNS_MPI_TRPB_K(p_b, up, ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))
  ENDIF
#endif

  NULLIFY(p_a,p_b)

  ENDIF

  RETURN
END SUBROUTINE PARTIAL_Z
