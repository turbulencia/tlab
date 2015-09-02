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
!# i1bc      In  grid structure: 0  periodic
!#                               1  non-periodic
!# bcs_imin  In  BC derivative at imin: 0  biased, non-zero
!#                                      1  forced to zero
!# bcs_imax  In  BC derivative at imax: 0  biased, non-zero
!#                                      1  forced to zero
!#
!########################################################################
SUBROUTINE PARTIAL_X(imode_fdm, nx,ny,nz, i1bc, dx, u,up, bcs_imin,bcs_imax, wrk1d,wrk2d,wrk3d)

  USE DNS_GLOBAL, ONLY : imax_total, inb_grid_1
#ifdef USE_MPI
  USE DNS_MPI
#endif

  IMPLICIT NONE

  TINTEGER imode_fdm, nx,ny,nz, i1bc, bcs_imin, bcs_imax
  TREAL, DIMENSION(imax_total,*)        :: dx
  TREAL, DIMENSION(nx*ny*nz),    TARGET :: u, up, wrk3d
  TREAL, DIMENSION(*)                   :: wrk1d ! not used, to be removed
  TREAL, DIMENSION(ny*nz)               :: wrk2d
 
! -------------------------------------------------------------------
  TINTEGER ip, nyz

  TREAL, DIMENSION(:), POINTER :: p_a, p_b, p_c

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
! Make  x  direction the last one
! -------------------------------------------------------------------
#ifdef USE_ESSL
  CALL DGETMO(p_a, imax_total, imax_total, nyz,        p_b, nyz)
#else
  CALL DNS_TRANSPOSE(p_a, imax_total, nyz,        imax_total, p_b, nyz)
#endif

! ###################################################################
! -------------------------------------------------------------------
! Periodic case
! -------------------------------------------------------------------
  IF ( i1bc .EQ. 0 ) THEN
     IF      ( imode_fdm .eq. FDM_COM4_JACOBIAN                                 ) THEN; CALL FDM_C1N4P_RHS(imax_total,nyz, p_b, p_c)
     ELSE IF ( imode_fdm .eq. FDM_COM6_JACOBIAN .OR. imode_fdm .EQ. FDM_COM6_DIRECT ) THEN; CALL FDM_C1N6P_RHS(imax_total,nyz, p_b, p_c)
     ELSE IF ( imode_fdm .eq. FDM_COM8_JACOBIAN                                 ) THEN; CALL FDM_C1N8P_RHS(imax_total,nyz, p_b, p_c)
     ENDIF

     ip = inb_grid_1 - 1
     CALL TRIDPSS(imax_total,nyz, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3),dx(1,ip+4),dx(1,ip+5), p_c,wrk2d)

! -------------------------------------------------------------------
! Nonperiodic case
! -------------------------------------------------------------------
  ELSE
     IF      ( imode_fdm .eq. FDM_COM4_JACOBIAN ) THEN; CALL FDM_C1N4_RHS(imax_total,nyz, bcs_imin,bcs_imax, p_b, p_c)
     ELSE IF ( imode_fdm .eq. FDM_COM6_JACOBIAN ) THEN; CALL FDM_C1N6_RHS(imax_total,nyz, bcs_imin,bcs_imax, p_b, p_c)
     ELSE IF ( imode_fdm .eq. FDM_COM8_JACOBIAN ) THEN; CALL FDM_C1N8_RHS(imax_total,nyz, bcs_imin,bcs_imax, p_b, p_c)
     ELSE IF ( imode_fdm .eq. FDM_COM6_DIRECT   ) THEN; CALL FDM_C1N6_RHS(imax_total,nyz, bcs_imin,bcs_imax, p_b, p_c) ! not yet implemented
     ENDIF
     
     ip = inb_grid_1 + (bcs_imin + bcs_imax*2)*3 - 1
     CALL TRIDSS(imax_total,nyz, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3), p_c)

  ENDIF

! ###################################################################
! -------------------------------------------------------------------
! Put arrays back in the order in which they came in
! -------------------------------------------------------------------
#ifdef USE_ESSL
  CALL DGETMO(p_c, nyz, nyz,        imax_total, p_b, imax_total)
#else
  CALL DNS_TRANSPOSE(p_c, nyz, imax_total, nyz,        p_b, imax_total)
#endif

! -------------------------------------------------------------------
! Transposition
! -------------------------------------------------------------------
#ifdef USE_MPI         
  IF ( ims_npro_i .GT. 1 ) THEN
     CALL DNS_MPI_TRPB_I(p_b, up, ims_ds_i(1,id), ims_dr_i(1,id), ims_ts_i(1,id), ims_tr_i(1,id))
  ENDIF
#endif

  NULLIFY(p_a,p_b,p_c)

  RETURN
END SUBROUTINE PARTIAL_X
