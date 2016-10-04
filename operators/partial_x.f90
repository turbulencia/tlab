#include "types.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

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
!# bcs_imin  In  BC derivative at imin: 0  biased, non-zero
!#                                      1  forced to zero
!# bcs_imax  In  BC derivative at imax: 0  biased, non-zero
!#                                      1  forced to zero
!#
!########################################################################
SUBROUTINE PARTIAL_X(imode_fdm, nx,ny,nz, i1bc, dx, u,up, bcs_imin,bcs_imax, wrk1d,wrk2d,wrk3d)

  USE DNS_GLOBAL, ONLY : g
#ifdef USE_MPI
  USE DNS_MPI
#endif

  IMPLICIT NONE

  TINTEGER imode_fdm, nx,ny,nz, i1bc, bcs_imin, bcs_imax
  TREAL, DIMENSION(*)        :: dx
  TREAL, DIMENSION(nx*ny*nz),    TARGET :: u, up, wrk3d
  TREAL, DIMENSION(*)                   :: wrk1d ! not used, to be removed
  TREAL, DIMENSION(ny*nz)               :: wrk2d
 
! -------------------------------------------------------------------
  TINTEGER nyz

  TREAL, DIMENSION(:), POINTER :: p_a, p_b, p_c

#ifdef USE_MPI
  TINTEGER, PARAMETER :: id = DNS_MPI_I_PARTIAL
#endif

! ###################################################################

! -------------------------------------------------------------------
! MPI transposition
! -------------------------------------------------------------------
#ifdef USE_MPI         
  IF ( ims_npro_i .GT. 1 ) THEN
     CALL DNS_MPI_TRPF_I(u, up, ims_ds_i(1,id), ims_dr_i(1,id), ims_ts_i(1,id), ims_tr_i(1,id))
     p_a => up
     p_b => wrk3d
     p_c => up
     nyz = ims_size_i(id)
  ELSE
#endif
     p_a => u
     p_b => up
     p_c => wrk3d
     nyz = ny*nz    
#ifdef USE_MPI         
  ENDIF
#endif

! -------------------------------------------------------------------
! Local transposition: make x-direction the last one
! -------------------------------------------------------------------
#ifdef USE_ESSL
  CALL DGETMO       (p_a, g(1)%size, g(1)%size, nyz,       p_b, nyz)
#else
  CALL DNS_TRANSPOSE(p_a, g(1)%size, nyz,       g(1)%size, p_b, nyz)
#endif

! ###################################################################
  CALL OPR_PARTIAL1(nyz, g(1), p_b,p_c, bcs_imin,bcs_imax, wrk2d)
  
! ###################################################################
! Put arrays back in the order in which they came in
#ifdef USE_ESSL
  CALL DGETMO       (p_c, nyz, nyz,       g(1)%size, p_b, g(1)%size)
#else
  CALL DNS_TRANSPOSE(p_c, nyz, g(1)%size, nyz,       p_b, g(1)%size)
#endif

#ifdef USE_MPI         
  IF ( ims_npro_i .GT. 1 ) THEN
     CALL DNS_MPI_TRPB_I(p_b, up, ims_ds_i(1,id), ims_dr_i(1,id), ims_ts_i(1,id), ims_tr_i(1,id))
  ENDIF
#endif

  NULLIFY(p_a,p_b,p_c)

  RETURN
END SUBROUTINE PARTIAL_X
