#include "types.h"
#include "dns_const.h"
#include "dns_error.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

!########################################################################
!# Tool/Library PADE
!#
!########################################################################
!# HISTORY
!#
!# 2011/11/01 - J.P. Mellado
!#              Created
!# 2013/01/20 - J.P. Mellado
!#              Introducing direct formulation of non-uniform grid
!#
!########################################################################
!# DESCRIPTION
!#
!# Apply the non-linear operator N(u) = visc* d^2/dx^2 s - u d/dx s
!# Derived from OPR_PARTIAL_XX to avoid one MPI transposition
!#
!########################################################################
!# ARGUMENTS 
!#
!# i1bc       In   grid structure: 0  periodic
!#                                 1  non-periodic
!# bcs2_imin  In   BC derivative at imin: 0  biased, non-zero
!#                                        1  forced to zero
!# bcs2_imax  In   BC derivative at imax: 0  biased, non-zero
!#                                        1  forced to zero
!#
!# ivel       In   Flag indicating the array containing the velocity:
!#                    0 for velocity being the scalar itself
!#                    1 for velocity passed through u1, or u2 if transposed required
!# is         In   Scalar index; if 0, then velocity
!# result     Out  Result N(u)
!# tmp1       Out  Transpose velocity
!#
!########################################################################
SUBROUTINE OPR_BURGERS_X(ivel, is, iunif,imode_fdm, nx,ny,nz, i1bc, dx, s,u1,u2, result, &
     bcs1_imin,bcs1_imax, bcs2_imin,bcs2_imax, tmp1, wrk2d,wrk3d)

  USE DNS_GLOBAL, ONLY : imax_total, inb_grid_1, inb_grid_2, inb_grid_3
  USE DNS_CONSTANTS, ONLY : efile
#ifdef USE_MPI
  USE DNS_MPI
#endif

  IMPLICIT NONE

  TINTEGER ivel, is, iunif
  TINTEGER imode_fdm, nx,ny,nz, i1bc, bcs1_imin,bcs1_imax, bcs2_imin,bcs2_imax
  TREAL, DIMENSION(imax_total,*), INTENT(IN)    :: dx
  TREAL, DIMENSION(nx*ny*nz),     INTENT(IN)    :: s,u1,u2
  TREAL, DIMENSION(nx*ny*nz),     INTENT(OUT)   :: result
  TREAL, DIMENSION(nx*ny*nz),     INTENT(INOUT) :: tmp1, wrk3d
  TREAL, DIMENSION(ny*nz),        INTENT(INOUT) :: wrk2d

  TARGET s,u1,u2, tmp1, result, wrk3d

! -------------------------------------------------------------------
  TINTEGER ip, nyz, ij

  TREAL, DIMENSION(:), POINTER :: p_a,p_b,p_c,p_d, p_vel

#ifdef USE_MPI
  TINTEGER id
#endif

! ###################################################################
  IF ( bcs2_imin + bcs2_imax .GT. 0 ) THEN
     CALL IO_WRITE_ASCII(efile,'OPR_BURGERS_X. Only developed for biased BCs.')
     CALL DNS_STOP(DNS_ERROR_UNDEVELOP)
  ENDIF

#ifdef USE_MPI
  id  = DNS_MPI_I_PARTIAL
#endif

! -------------------------------------------------------------------
! MPI transposition
! -------------------------------------------------------------------
#ifdef USE_MPI         
  IF ( ims_npro_i .GT. 1 ) THEN
     CALL DNS_MPI_TRPF_I(s, result, ims_ds_i(1,id), ims_dr_i(1,id), ims_ts_i(1,id), ims_tr_i(1,id))
     p_a => result
     p_b => tmp1
     p_c => wrk3d
     p_d => result
     nyz = ims_size_i(id)
  ELSE
#endif
     p_a => s
     p_b => tmp1
     p_c => result
     p_d => wrk3d
     nyz = ny*nz    
#ifdef USE_MPI         
  ENDIF
#endif

! pointer to velocity
  IF ( ivel .EQ. 0 ) THEN; p_vel => p_b
  ELSE;                    p_vel => u2; ENDIF ! always transposed needed

! -------------------------------------------------------------------
! Local transposition: make  x  direction the last one
! -------------------------------------------------------------------
#ifdef USE_ESSL
  CALL DGETMO       (p_a, imax_total, imax_total, nyz,        p_b, nyz)
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

! Calculate second derivative
     IF      ( imode_fdm .eq. FDM_COM4_JACOBIAN                                 ) THEN; CALL FDM_C2N4P_RHS(imax_total,nyz, p_b, p_d)
     ELSE IF ( imode_fdm .eq. FDM_COM6_JACOBIAN .OR. imode_fdm .EQ. FDM_COM6_DIRECT ) THEN; CALL FDM_C2N6P_RHS(imax_total,nyz, p_b, p_d)
     ELSE IF ( imode_fdm .eq. FDM_COM8_JACOBIAN                                 ) THEN; CALL FDM_C2N6P_RHS(imax_total,nyz, p_b, p_d) !8th not yet developed
     ENDIF
     ip = inb_grid_3 + is*5 - 1 ! LU decomposition containing the diffusivity
     CALL TRIDPSS(imax_total,nyz, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3),dx(1,ip+4),dx(1,ip+5), p_d,wrk2d)

! -------------------------------------------------------------------
! Nonperiodic case
! -------------------------------------------------------------------
  ELSE
     IF      ( imode_fdm .eq. FDM_COM4_JACOBIAN ) THEN; CALL FDM_C1N4_RHS(imax_total,nyz, bcs1_imin,bcs1_imax, p_b, p_c)
     ELSE IF ( imode_fdm .eq. FDM_COM6_JACOBIAN ) THEN; CALL FDM_C1N6_RHS(imax_total,nyz, bcs1_imin,bcs1_imax, p_b, p_c)
     ELSE IF ( imode_fdm .eq. FDM_COM8_JACOBIAN ) THEN; CALL FDM_C1N8_RHS(imax_total,nyz, bcs1_imin,bcs1_imax, p_b, p_c)
     ELSE IF ( imode_fdm .eq. FDM_COM6_DIRECT   ) THEN; CALL FDM_C1N6_RHS(imax_total,nyz, bcs1_imin,bcs1_imax, p_b, p_c) ! not yet implemented
     ENDIF
     ip = inb_grid_1 + (bcs1_imin + bcs1_imax*2)*3 - 1
     CALL TRIDSS(imax_total,nyz, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3), p_c)

! Calculate second derivative
     IF      ( imode_fdm .eq. FDM_COM4_JACOBIAN ) THEN; CALL FDM_C2N4_RHS(iunif, imax_total,nyz, bcs2_imin,bcs2_imax, dx, p_b, p_c, p_d)
     ELSE IF ( imode_fdm .eq. FDM_COM6_JACOBIAN ) THEN; CALL FDM_C2N6_RHS(iunif, imax_total,nyz, bcs2_imin,bcs2_imax, dx, p_b, p_c, p_d)
     ELSE IF ( imode_fdm .eq. FDM_COM8_JACOBIAN ) THEN; CALL FDM_C2N6_RHS(iunif, imax_total,nyz, bcs2_imin,bcs2_imax, dx, p_b, p_c, p_d)
     ELSE IF ( imode_fdm .eq. FDM_COM6_DIRECT   ) THEN; CALL FDM_C2N6N_RHS(imax_total,nyz, dx(1,inb_grid_2+3), p_b, p_d)
     ENDIF
     ip = inb_grid_3 + is*5 - 1 ! LU decomposition containing the diffusivity
     CALL TRIDSS(imax_total,nyz, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3), p_d)

  ENDIF

! ###################################################################
! Operation
! ###################################################################
!$omp parallel default( shared ) private( ij )
!$omp do
  DO ij = 1,nx*ny*nz
     p_d(ij) = p_d(ij) - p_vel(ij)*p_c(ij) ! diffusivity included in 2.-order derivative
  ENDDO
!$omp end do
!$omp end parallel

! ###################################################################
! Put arrays back in the order in which they came in
#ifdef USE_ESSL
  CALL DGETMO       (p_d, nyz, nyz,        imax_total, p_c, imax_total)
#else
  CALL DNS_TRANSPOSE(p_d, nyz, imax_total, nyz,        p_c, imax_total)
#endif

#ifdef USE_MPI         
  IF ( ims_npro_i .GT. 1 ) THEN
     CALL DNS_MPI_TRPB_I(p_c, result, ims_ds_i(1,id), ims_dr_i(1,id), ims_ts_i(1,id), ims_tr_i(1,id))
  ENDIF
#endif

  NULLIFY(p_a,p_b,p_c,p_d, p_vel)

  RETURN
END SUBROUTINE OPR_BURGERS_X
