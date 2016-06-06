#include "types.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2007/11/12 - J.P. Mellado
!#              Created
!# 2016/05/27 - J.P. Mellado
!#              Parallelizing
!#
!########################################################################
!# DESCRIPTION
!#
!########################################################################
SUBROUTINE OPR_FILTER_Z(imode_filter, nx,ny,nz, k1bc, bcs_kmin,bcs_kmax, u, cz, tmp, wrk1d,wrk2d,wrk3d)

  USE DNS_GLOBAL, ONLY : kmax_total
#ifdef USE_MPI
  USE DNS_MPI
#endif

  IMPLICIT NONE

  TINTEGER,                       INTENT(IN)            :: imode_filter, nx,ny,nz, k1bc, bcs_kmin,bcs_kmax
  TREAL, DIMENSION(nx*ny*nz),     INTENT(INOUT), TARGET :: u, wrk3d    ! in-place operation
  TREAL, DIMENSION(kmax_total,6), INTENT(IN)            :: cz          ! Filter kernel information
  TREAL, DIMENSION(nx*ny),        INTENT(INOUT)         :: wrk2d       ! Aux arrays
  TREAL, DIMENSION(nx*ny*nz),     INTENT(INOUT)         :: tmp         ! Aux array needed in ADM type
  TREAL, DIMENSION(kmax_total,*), INTENT(INOUT)         :: wrk1d
 
! -------------------------------------------------------------------
  TINTEGER nxy

  TREAL, DIMENSION(:), POINTER :: p_a, p_b

#ifdef USE_MPI
  TINTEGER id
#endif

! ###################################################################
#ifdef USE_MPI         
  id = DNS_MPI_K_PARTIAL
#endif

! -------------------------------------------------------------------
! Transposition
! -------------------------------------------------------------------
#ifdef USE_MPI         
  IF ( ims_npro_k .GT. 1 ) THEN
     CALL DNS_MPI_TRPF_K(u, wrk3d, ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))
     p_a => wrk3d
     p_b => u
     nxy = ims_size_k(id)
  ELSE
#endif
     p_a => u
     p_b => wrk3d
     nxy = nx*ny
#ifdef USE_MPI         
  ENDIF
#endif

! ###################################################################
  IF      ( imode_filter .EQ. DNS_FILTER_COMPACT ) THEN; CALL FILT4C_KERNEL(kmax_total,nxy, p_a,p_b, k1bc, bcs_kmin,bcs_kmax, wrk1d, cz)
     IF ( k1bc .EQ. 0 ) THEN
        CALL TRIDPFS(kmax_total,     wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5))
        CALL TRIDPSS(kmax_total,nxy, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5), p_b, wrk2d)
     ELSE
        CALL TRIDFS(kmax_total,     wrk1d(1,1),wrk1d(1,2),wrk1d(1,3))
        CALL TRIDSS(kmax_total,nxy, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3), p_b)
     ENDIF

  ELSE IF ( imode_filter .EQ. DNS_FILTER_6E      ) THEN; CALL FILT6E_KERNEL (nxy,kmax_total, k1bc, bcs_kmin, bcs_kmax, p_a,p_b)
  ELSE IF ( imode_filter .EQ. DNS_FILTER_4E      ) THEN; CALL FILT4E_KERNEL (kmax_total,nxy, k1bc,                     p_a,p_b, cz)
  ELSE IF ( imode_filter .EQ. DNS_FILTER_ADM     ) THEN; CALL FILTADM_KERNEL(kmax_total,nxy, k1bc,                     p_a,p_b, tmp, cz)
  ENDIF

! ###################################################################
! -------------------------------------------------------------------
! Transposition
! -------------------------------------------------------------------
#ifdef USE_MPI         
  IF ( ims_npro_k .GT. 1 ) THEN
     CALL DNS_MPI_TRPB_K(p_b, p_a, ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))
  ENDIF
#endif

  u(1:nx*ny*nz) = wrk3d(1:nx*ny*nz)

  NULLIFY(p_a,p_b)

  RETURN
END SUBROUTINE OPR_FILTER_Z
