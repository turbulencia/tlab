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
!# i1bc       In  grid structure: 0  periodic
!#                                1  non-periodic
!# bcs2_imin  In  BC derivative at imin: 0  biased, non-zero
!#                                       1  forced to zero
!# bcs2_imax  In  BC derivative at imax: 0  biased, non-zero
!#                                       1  forced to zero
!#
!# up2        Out second derivative
!# up1        Out first derivative, if ifirst=1
!#
!########################################################################
SUBROUTINE PARTIAL_XX(ifirst,iunif,imode_fdm, nx,ny,nz, i1bc, dx, u, up2, &
     bcs1_imin,bcs1_imax, bcs2_imin,bcs2_imax, up1, wrk1d,wrk2d,wrk3d)

  USE DNS_GLOBAL, ONLY : g
#ifdef USE_MPI
  USE DNS_MPI
#endif

  IMPLICIT NONE

  TINTEGER ifirst
  TINTEGER nx,ny,nz, bcs1_imin,bcs1_imax, bcs2_imin,bcs2_imax
  TREAL, DIMENSION(nx*ny*nz),    TARGET :: u, up1, up2, wrk3d
  TREAL, DIMENSION(ny*nz)               :: wrk2d

  TREAL, DIMENSION(*)                   :: dx, wrk1d ! not used, to be removed
  TINTEGER imode_fdm, iunif, i1bc

! -------------------------------------------------------------------
  TINTEGER nyz, bcs(2,2)

  TREAL, DIMENSION(:), POINTER :: p_a, p_b, p_c, p_d

#ifdef USE_MPI
  TINTEGER, PARAMETER :: id = DNS_MPI_I_PARTIAL
#endif

! ###################################################################
  bcs(1,1) = bcs1_imin; bcs(2,1) = bcs1_imax ! 1. order derivative
  bcs(1,2) = bcs2_imin; bcs(2,2) = bcs2_imax ! 2. order derivative

! -------------------------------------------------------------------
! MPI transposition
! -------------------------------------------------------------------
#ifdef USE_MPI         
  IF ( ims_npro_i .GT. 1 ) THEN
     CALL DNS_MPI_TRPF_I(u, up2, ims_ds_i(1,id), ims_dr_i(1,id), ims_ts_i(1,id), ims_tr_i(1,id))
     p_a => up2
     p_b => up1
     p_c => up2
     p_d => wrk3d
     nyz = ims_size_i(id)
  ELSE
#endif
     p_a => u
     p_b => up2
     p_c => wrk3d
     p_d => up1
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
  CALL OPR_PARTIAL2(nyz, bcs, g(1), p_b,p_d, wrk2d,p_c)
  
! Check whether we need to calculate the 1. order derivative
  IF ( ifirst .EQ. 1 ) THEN
     IF ( g(1)%uniform .OR. g(1)%mode_fdm .EQ. FDM_COM6_DIRECT ) THEN
        CALL OPR_PARTIAL1(nyz, bcs, g(1), p_b,p_c, wrk2d)
     ENDIF
  ENDIF  

! ###################################################################
! Put arrays back in the order in which they came in
#ifdef USE_ESSL
  CALL DGETMO(p_d, nyz, nyz,        g(1)%size, p_b, g(1)%size)
#else
  CALL DNS_TRANSPOSE(p_d, nyz, g(1)%size, nyz,        p_b, g(1)%size)
#endif

  IF ( ifirst .EQ. 1 ) THEN
#ifdef USE_ESSL
  CALL DGETMO(p_c, nyz, nyz,        g(1)%size, p_d, g(1)%size)
#else
  CALL DNS_TRANSPOSE(p_c, nyz, g(1)%size, nyz,        p_d, g(1)%size)
#endif
  ENDIF

#ifdef USE_MPI         
  IF ( ims_npro_i .GT. 1 ) THEN
     CALL DNS_MPI_TRPB_I(p_b, up2, ims_ds_i(1,id), ims_dr_i(1,id), ims_ts_i(1,id), ims_tr_i(1,id))
     IF ( ifirst .EQ. 1 ) THEN ! only if you really want first derivative back
        CALL DNS_MPI_TRPB_I(p_d, up1, ims_ds_i(1,id), ims_dr_i(1,id), ims_ts_i(1,id), ims_tr_i(1,id))
     ENDIF
  ENDIF
#endif

  NULLIFY(p_a,p_b,p_c,p_d)

  RETURN
END SUBROUTINE PARTIAL_XX
