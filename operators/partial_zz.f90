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
!# k1bc       In  grid structure: 0  periodic
!#                                1  non-periodic
!# bcs2_kmin  In  BC derivative at kmin: 0  biased, non-zero
!#                                       1  forced to zero
!# bcs2_kmax  In  BC derivative at kmax: 0  biased, non-zero
!#                                       1  forced to zero
!#
!# up2        Out second derivative
!# up1        Out first derivative, if ifirst=1
!#
!########################################################################
SUBROUTINE PARTIAL_ZZ(ifirst,iunif,imode_fdm, nx,ny,nz, k1bc, dz, u, up2, &
     bcs1_kmin,bcs1_kmax, bcs2_kmin,bcs2_kmax, up1, wrk1d,wrk2d,wrk3d)

  USE DNS_GLOBAL, ONLY : g
#ifdef USE_MPI
  USE DNS_MPI
#endif

  IMPLICIT NONE

  TINTEGER ifirst, iunif
  TINTEGER imode_fdm, nx,ny,nz, k1bc, bcs1_kmin, bcs1_kmax, bcs2_kmin, bcs2_kmax
  TREAL, DIMENSION(nx*ny*nz),    TARGET :: u, up1, up2, wrk3d
  TREAL, DIMENSION(*)                   :: dz, wrk1d ! not used, to be removed
  TREAL, DIMENSION(nx*ny)               :: wrk2d

! -------------------------------------------------------------------
  TINTEGER nxy, bcs_min(2), bcs_max(2)!, ifirst_loc

  TREAL, DIMENSION(:), POINTER :: p_a, p_b, p_c

#ifdef USE_MPI
  TINTEGER, PARAMETER :: id = DNS_MPI_K_PARTIAL
#endif

! ###################################################################
  bcs_min(1) = bcs1_kmin; bcs_max(1) = bcs1_kmax
  bcs_min(2) = bcs2_kmin; bcs_max(2) = bcs2_kmax

  IF ( g(3)%size .EQ. 1 ) THEN ! Set to zero in 2D case
     up2 = C_0_R
     IF ( ifirst .EQ. 1 ) up1 = C_0_R
     
  ELSE
! ###################################################################
! -------------------------------------------------------------------
! MPI Transposition
! -------------------------------------------------------------------
#ifdef USE_MPI         
  IF ( ims_npro_k .GT. 1 ) THEN
     CALL DNS_MPI_TRPF_K(u, up2, ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))
     p_a => up2
     p_b => wrk3d
     p_c => up1
     nxy = ims_size_k(id)
 ELSE
#endif
    p_a => u
    p_b => up1
    p_c => up2
    nxy = nx*ny 
#ifdef USE_MPI         
  ENDIF
#endif

! ###################################################################
  CALL OPR_PARTIAL2(nxy, g(3), p_a,p_c, bcs_min,bcs_max, wrk2d,p_b)
  
! Check whether we need to calculate the 1. order derivative
  IF ( ifirst .EQ. 1 ) THEN
     IF ( g(3)%uniform .OR. imode_fdm .EQ. FDM_COM6_DIRECT ) THEN
        CALL OPR_PARTIAL1(nxy, g(3), p_a,p_b, bcs_min(1),bcs_max(1), wrk2d)
     ENDIF
  ENDIF
  
! ###################################################################
! Put arrays back in the order in which they came in
#ifdef USE_MPI         
  IF ( ims_npro_k .GT. 1 ) THEN
     CALL DNS_MPI_TRPB_K(p_c, up2, ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))
     IF ( ifirst .EQ. 1 ) THEN ! only if you really want first derivative back
        CALL DNS_MPI_TRPB_K(p_b, up1, ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))
     ENDIF
  ENDIF
#endif

  NULLIFY(p_a,p_b,p_c)

  ENDIF

  RETURN
END SUBROUTINE PARTIAL_ZZ
