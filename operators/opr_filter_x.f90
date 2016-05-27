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
SUBROUTINE OPR_FILTER_X(imode_filter, nx,ny,nz, i1bc, bcs_imin,bcs_imax, u, cx, tmp, wrk1d,wrk2d,wrk3d)

  USE DNS_GLOBAL, ONLY : imax_total
#ifdef USE_MPI
  USE DNS_MPI
#endif

  IMPLICIT NONE

  TINTEGER,                       INTENT(IN)            :: imode_filter, nx,ny,nz, i1bc, bcs_imin,bcs_imax
  TREAL, DIMENSION(nx*ny*nz),     INTENT(INOUT), TARGET :: u, wrk3d    ! in-place operation
  TREAL, DIMENSION(imax_total,6), INTENT(IN)            :: cx          ! Filter kernel information
  TREAL, DIMENSION(ny*nz),        INTENT(INOUT)         :: wrk2d       ! Aux arrays
  TREAL, DIMENSION(nx*ny*nz),     INTENT(INOUT)         :: tmp         ! Aux array needed in ADM type
  TREAL, DIMENSION(imax_total,*), INTENT(INOUT)         :: wrk1d
 
! -------------------------------------------------------------------
  TINTEGER nyz

  TREAL, DIMENSION(:), POINTER :: p_a, p_b

#ifdef USE_MPI
  TINTEGER id
#endif

! ###################################################################
#ifdef USE_MPI         
  id = DNS_MPI_I_PARTIAL
#endif

! -------------------------------------------------------------------
! Transposition
! -------------------------------------------------------------------
#ifdef USE_MPI         
  IF ( ims_npro_i .GT. 1 ) THEN
     CALL DNS_MPI_TRPF_I(u, wrk3d, ims_ds_i(1,id), ims_dr_i(1,id), ims_ts_i(1,id), ims_tr_i(1,id))
     p_a => wrk3d
     p_b => u
     nyz = ims_size_i(id)
  ELSE
#endif
     p_a => u
     p_b => wrk3d
     nyz = ny*nz    
#ifdef USE_MPI         
  ENDIF
#endif

! -------------------------------------------------------------------
! Make  x  direction the last one
! -------------------------------------------------------------------
#ifdef USE_ESSL
  CALL DGETMO       (p_a, imax_total, imax_total, nyz,        p_b, nyz)
#else
  CALL DNS_TRANSPOSE(p_a, imax_total, nyz,        imax_total, p_b, nyz)
#endif

! ###################################################################
  IF      ( imode_filter .EQ. DNS_FILTER_COMPACT ) THEN; CALL FILT4C_KERNEL(imax_total,nyz, p_b,p_a, i1bc, bcs_imin,bcs_imax, wrk1d, cx)
     IF ( i1bc .EQ. 0 ) THEN
        CALL TRIDPFS(imax_total,     wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5))
        CALL TRIDPSS(imax_total,nyz, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5), p_a, wrk2d)
     ELSE
        CALL TRIDFS(imax_total,     wrk1d(1,1),wrk1d(1,2),wrk1d(1,3))
        CALL TRIDSS(imax_total,nyz, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3), p_a)
     ENDIF

  ELSE IF ( imode_filter .EQ. DNS_FILTER_6E      ) THEN; CALL FILT6E_KERNEL (nyz,imax_total, i1bc, bcs_imin, bcs_imax, p_b,p_a)
  ELSE IF ( imode_filter .EQ. DNS_FILTER_4E      ) THEN; CALL FILT4E_KERNEL (imax_total,nyz, i1bc,                     p_b,p_a, cx)
  ELSE IF ( imode_filter .EQ. DNS_FILTER_ADM     ) THEN; CALL FILTADM_KERNEL(imax_total,nyz, i1bc,                     p_b,p_a, tmp, cx)
  ENDIF

! ###################################################################
! -------------------------------------------------------------------
! Put arrays back in the order in which they came in
! -------------------------------------------------------------------
#ifdef USE_ESSL
  CALL DGETMO       (p_a, nyz, nyz,        imax_total, p_b, imax_total)
#else
  CALL DNS_TRANSPOSE(p_a, nyz, imax_total, nyz,        p_b, imax_total)
#endif

! -------------------------------------------------------------------
! Transposition
! -------------------------------------------------------------------
#ifdef USE_MPI         
  IF ( ims_npro_i .GT. 1 ) THEN
     CALL DNS_MPI_TRPB_I(p_b, p_a, ims_ds_i(1,id), ims_dr_i(1,id), ims_ts_i(1,id), ims_tr_i(1,id))
  ENDIF
#endif

  u(1:nx*ny*nz) = wrk3d(1:nx*ny*nz)

  NULLIFY(p_a,p_b)

  RETURN
END SUBROUTINE OPR_FILTER_X

!########################################################################
!########################################################################

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
SUBROUTINE OPR_FILTER_X_OLD(itype, nx, ny, nz, i1bc, bcs_imin, bcs_imax, u, cx, tmp, wrk1d,wrk2d,wrk3d)

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
END SUBROUTINE OPR_FILTER_X_OLD
