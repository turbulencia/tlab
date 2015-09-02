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

  USE DNS_GLOBAL, ONLY : imax_total, inb_grid_1, inb_grid_2
#ifdef USE_MPI
  USE DNS_MPI
#endif

  IMPLICIT NONE

  TINTEGER ifirst, iunif
  TINTEGER imode_fdm, nx,ny,nz, i1bc, bcs1_imin,bcs1_imax, bcs2_imin,bcs2_imax
  TREAL, DIMENSION(imax_total,*)        :: dx
  TREAL, DIMENSION(nx*ny*nz),    TARGET :: u, up1, up2, wrk3d
  TREAL, DIMENSION(*)                   :: wrk1d ! not used, to be removed
  TREAL, DIMENSION(ny*nz)               :: wrk2d

! -------------------------------------------------------------------
  TINTEGER ip, nyz, ifirst_loc

  TREAL, DIMENSION(:), POINTER :: p_a, p_b, p_c, p_d

#ifdef USE_MPI
  TINTEGER id
#endif

! ###################################################################
#ifdef USE_MPI
  id  = DNS_MPI_I_PARTIAL
#endif

! Check whether we need to calculate the 1. order derivative
  ifirst_loc = ifirst
  IF ( iunif .NE. 0 ) THEN
     IF ( imode_fdm .eq. FDM_COM4_JACOBIAN .OR. &
          imode_fdm .eq. FDM_COM6_JACOBIAN .OR. &
          imode_fdm .eq. FDM_COM8_JACOBIAN      ) THEN; ifirst_loc = MAX(ifirst_loc,1)
     ENDIF
  ENDIF

! -------------------------------------------------------------------
! Transposition
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
     IF ( ifirst_loc .EQ. 1 ) THEN ! First derivative
     IF      ( imode_fdm .EQ. FDM_COM4_JACOBIAN                                 ) THEN; CALL FDM_C1N4P_RHS(imax_total,nyz, p_b, p_c)
     ELSE IF ( imode_fdm .EQ. FDM_COM6_JACOBIAN .OR. imode_fdm .EQ. FDM_COM6_DIRECT ) THEN; CALL FDM_C1N6P_RHS(imax_total,nyz, p_b, p_c)
     ELSE IF ( imode_fdm .EQ. FDM_COM8_JACOBIAN                                 ) THEN; CALL FDM_C1N8P_RHS(imax_total,nyz, p_b, p_c)
     ENDIF
     ip = inb_grid_1 - 1
     CALL TRIDPSS(imax_total,nyz, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3),dx(1,ip+4),dx(1,ip+5), p_c,wrk2d)
     ENDIF

     IF      ( imode_fdm .EQ. FDM_COM4_JACOBIAN                                 ) THEN; CALL FDM_C2N4P_RHS(imax_total,nyz, p_b, p_d)
     ELSE IF ( imode_fdm .EQ. FDM_COM6_JACOBIAN .OR. imode_fdm .EQ. FDM_COM6_DIRECT ) THEN; CALL FDM_C2N6P_RHS(imax_total,nyz, p_b, p_d)
     ELSE IF ( imode_fdm .EQ. FDM_COM8_JACOBIAN                                 ) THEN; CALL FDM_C2N6P_RHS(imax_total,nyz, p_b, p_d) !8th not yet developed
     ENDIF
     ip = inb_grid_2 - 1
     CALL TRIDPSS(imax_total,nyz, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3),dx(1,ip+4),dx(1,ip+5), p_d,wrk2d)

! -------------------------------------------------------------------
! Nonperiodic case
! -------------------------------------------------------------------
  ELSE
     IF ( ifirst_loc .EQ. 1 ) THEN ! First derivative
     IF      ( imode_fdm .eq. FDM_COM4_JACOBIAN ) THEN; CALL FDM_C1N4_RHS(imax_total,nyz, bcs1_imin,bcs1_imax, p_b, p_c)
     ELSE IF ( imode_fdm .eq. FDM_COM6_JACOBIAN ) THEN; CALL FDM_C1N6_RHS(imax_total,nyz, bcs1_imin,bcs1_imax, p_b, p_c)
     ELSE IF ( imode_fdm .eq. FDM_COM8_JACOBIAN ) THEN; CALL FDM_C1N8_RHS(imax_total,nyz, bcs1_imin,bcs1_imax, p_b, p_c)
     ELSE IF ( imode_fdm .eq. FDM_COM6_DIRECT   ) THEN; CALL FDM_C1N6_RHS(imax_total,nyz, bcs1_imin,bcs1_imax, p_b, p_c) ! not yet implemented
     ENDIF
     ip = inb_grid_1 + (bcs1_imin + bcs1_imax*2)*3 - 1
     CALL TRIDSS(imax_total,nyz, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3), p_c)
     ENDIF

     IF      ( imode_fdm .eq. FDM_COM4_JACOBIAN ) THEN; CALL FDM_C2N4_RHS(iunif, imax_total,nyz, bcs2_imin,bcs2_imax, dx, p_b, p_c, p_d)
     ELSE IF ( imode_fdm .eq. FDM_COM6_JACOBIAN ) THEN; CALL FDM_C2N6_RHS(iunif, imax_total,nyz, bcs2_imin,bcs2_imax, dx, p_b, p_c, p_d)
     ELSE IF ( imode_fdm .eq. FDM_COM8_JACOBIAN ) THEN; CALL FDM_C2N6_RHS(iunif, imax_total,nyz, bcs2_imin,bcs2_imax, dx, p_b, p_c, p_d) !8th not yet developed
     ELSE IF ( imode_fdm .eq. FDM_COM6_DIRECT   ) THEN; CALL FDM_C2N6N_RHS(imax_total,nyz, dx(1,inb_grid_2+3), p_b, p_d)
     ENDIF
     ip = inb_grid_2 + (bcs2_imin + bcs2_imax*2)*3 - 1
     CALL TRIDSS(imax_total,nyz, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3), p_d)

  ENDIF

! ###################################################################
! Put arrays back in the order in which they came in
#ifdef USE_ESSL
  CALL DGETMO(p_d, nyz, nyz,        imax_total, p_b, imax_total)
#else
  CALL DNS_TRANSPOSE(p_d, nyz, imax_total, nyz,        p_b, imax_total)
#endif

  IF ( ifirst .EQ. 1 ) THEN
#ifdef USE_ESSL
  CALL DGETMO(p_c, nyz, nyz,        imax_total, p_d, imax_total)
#else
  CALL DNS_TRANSPOSE(p_c, nyz, imax_total, nyz,        p_d, imax_total)
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
