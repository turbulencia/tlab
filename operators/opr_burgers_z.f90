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
!# Apply the non-linear operator N(u) = visc* d^2/dz^2 s - u d/dz s
!# Derived from OPR_PARTIAL_ZZ to avoid one MPI transposition
!#
!########################################################################
!# ARGUMENTS 
!#
!# k1bc       In  grid structure: 0  periodic
!#                                1  non-periodic
!# bcs2_kmin  In  BC derivative at kmin: 0  biased, non-zero
!#                                       1  forced to zero
!# bcs2_kmax  In  BC derivative at kmax: 0  biased, non-zero
!#                                       1  forced to zero
!#
!# ivel       In   Flag indicating the array containing the velocity:
!#                    0 for velocity being the scalar itself
!#                    1 for velocity passed through u1, or u2 if transposed required
!# is         In   Scalar index; if 0, then velocity
!# result     Out  Result N(u)
!# tmp1       Out  Transpose velocity
!#
!########################################################################
SUBROUTINE OPR_BURGERS_Z(ivel, is, iunif,imode_fdm, nx,ny,nz, k1bc, dz, s,u1,u2, result, &
     bcs1_kmin,bcs1_kmax, bcs2_kmin,bcs2_kmax, tmp1, wrk2d,wrk3d)

  USE DNS_GLOBAL, ONLY : kmax_total, inb_grid_1, inb_grid_2, inb_grid_3
  USE DNS_CONSTANTS, ONLY : efile
#ifdef USE_MPI
  USE DNS_MPI
#endif

  IMPLICIT NONE

  TINTEGER ivel, is, iunif
  TINTEGER imode_fdm, nx,ny,nz, k1bc, bcs1_kmin, bcs1_kmax, bcs2_kmin, bcs2_kmax
  TREAL, DIMENSION(kmax_total,*), INTENT(IN)    :: dz
  TREAL, DIMENSION(nx*ny*nz),     INTENT(IN)    :: s,u1,u2
  TREAL, DIMENSION(nx*ny*nz),     INTENT(OUT)   :: result
  TREAL, DIMENSION(nx*ny*nz),     INTENT(INOUT) :: tmp1, wrk3d
  TREAL, DIMENSION(nx*ny),        INTENT(INOUT) :: wrk2d

  TARGET s,u1,u2, tmp1, result, wrk3d

! -------------------------------------------------------------------
  TINTEGER ip, nxy, ij

  TREAL, DIMENSION(:), POINTER :: p_a,p_b,p_c, p_vel

#ifdef USE_MPI
  TINTEGER id
#endif

! ###################################################################
  IF ( bcs2_kmin + bcs2_kmax .GT. 0 ) THEN
     CALL IO_WRITE_ASCII(efile,'OPR_BURGERS_Z. Only developed for biased BCs.')
     CALL DNS_STOP(DNS_ERROR_UNDEVELOP)
  ENDIF

#ifdef USE_MPI
  id  = DNS_MPI_K_PARTIAL
#endif

  IF ( kmax_total .EQ. 1 ) THEN ! Set to zero in 2D case
  result = C_0_R
     
  ELSE
! ###################################################################
! Transposition
! ###################################################################
#ifdef USE_MPI         
  IF ( ims_npro_k .GT. 1 ) THEN
     CALL DNS_MPI_TRPF_K(s, tmp1, ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))
     p_a => tmp1
     p_b => result
     p_c => wrk3d
     nxy = ims_size_k(id)
  ELSE
#endif
     p_a => s
     p_b => tmp1
     p_c => result
     nxy = nx*ny 
#ifdef USE_MPI         
  ENDIF
#endif

! pointer to velocity
  IF ( ivel .EQ. 0 ) THEN; p_vel => p_a
  ELSE
#ifdef USE_MPI         
     IF ( ims_npro_k .GT. 1 ) THEN ! I do     need the transposed
        p_vel => u2
     ELSE
#endif
        p_vel => u1                ! I do not need the transposed
#ifdef USE_MPI         
     ENDIF
#endif
  ENDIF

! ###################################################################
! -------------------------------------------------------------------
! Periodic case
! -------------------------------------------------------------------
  IF ( k1bc .EQ. 0 ) THEN
     IF      ( imode_fdm .eq. FDM_COM4_JACOBIAN                                 ) THEN; CALL FDM_C1N4P_RHS(kmax_total,nxy, p_a, p_b)
     ELSE IF ( imode_fdm .eq. FDM_COM6_JACOBIAN .OR. imode_fdm .EQ. FDM_COM6_DIRECT ) THEN; CALL FDM_C1N6P_RHS(kmax_total,nxy, p_a, p_b)
     ELSE IF ( imode_fdm .eq. FDM_COM8_JACOBIAN                                 ) THEN; CALL FDM_C1N8P_RHS(kmax_total,nxy, p_a, p_b)
     ENDIF
     ip = inb_grid_1 - 1
     CALL TRIDPSS(kmax_total,nxy, dz(1,ip+1),dz(1,ip+2),dz(1,ip+3),dz(1,ip+4),dz(1,ip+5), p_b,wrk2d)

! Calculate second derivative
     IF      ( imode_fdm .eq. FDM_COM4_JACOBIAN                                 ) THEN; CALL FDM_C2N4P_RHS(kmax_total,nxy, p_a, p_c)
     ELSE IF ( imode_fdm .eq. FDM_COM6_JACOBIAN .OR. imode_fdm .EQ. FDM_COM6_DIRECT ) THEN; CALL FDM_C2N6P_RHS(kmax_total,nxy, p_a, p_c)
     ELSE IF ( imode_fdm .eq. FDM_COM8_JACOBIAN                                 ) THEN; CALL FDM_C2N6P_RHS(kmax_total,nxy, p_a, p_c) ! 8th not yet developed
     ENDIF
     ip = inb_grid_3 + is*5 - 1 ! LU decomposition containing the diffusivity
     CALL TRIDPSS(kmax_total,nxy, dz(1,ip+1),dz(1,ip+2),dz(1,ip+3),dz(1,ip+4),dz(1,ip+5), p_c,wrk2d)

! -------------------------------------------------------------------
! Nonperiodic case
! -------------------------------------------------------------------
  ELSE
     IF      ( imode_fdm .eq. FDM_COM4_JACOBIAN ) THEN; CALL FDM_C1N4_RHS(kmax_total,nxy, bcs1_kmin,bcs1_kmax, p_a, p_b)
     ELSE IF ( imode_fdm .eq. FDM_COM6_JACOBIAN ) THEN; CALL FDM_C1N6_RHS(kmax_total,nxy, bcs1_kmin,bcs1_kmax, p_a, p_b)
     ELSE IF ( imode_fdm .eq. FDM_COM8_JACOBIAN ) THEN; CALL FDM_C1N8_RHS(kmax_total,nxy, bcs1_kmin,bcs1_kmax, p_a, p_b)
     ELSE IF ( imode_fdm .eq. FDM_COM6_DIRECT   ) THEN; CALL FDM_C1N6_RHS(kmax_total,nxy, bcs1_kmin,bcs1_kmax, p_a, p_b)
     ENDIF
     ip = inb_grid_1 + (bcs1_kmin + bcs1_kmax*2)*3 - 1
     CALL TRIDSS(kmax_total,nxy, dz(1,ip+1),dz(1,ip+2),dz(1,ip+3), p_b)

     IF      ( imode_fdm .eq. FDM_COM4_JACOBIAN ) THEN; CALL FDM_C2N4_RHS(iunif, kmax_total,nxy, bcs2_kmin,bcs2_kmax, dz, p_a,p_b,p_c)
     ELSE IF ( imode_fdm .eq. FDM_COM6_JACOBIAN ) THEN; CALL FDM_C2N6_RHS(iunif, kmax_total,nxy, bcs2_kmin,bcs2_kmax, dz, p_a,p_b,p_c)
     ELSE IF ( imode_fdm .eq. FDM_COM8_JACOBIAN ) THEN; CALL FDM_C2N6_RHS(iunif, kmax_total,nxy, bcs2_kmin,bcs2_kmax, dz, p_a,p_b,p_c)
     ELSE IF ( imode_fdm .eq. FDM_COM6_DIRECT   ) THEN; CALL FDM_C2N6N_RHS(kmax_total,nxy, dz(1,inb_grid_2+3), p_a, p_c)
     ENDIF
     ip = inb_grid_3 + is*5 - 1 ! LU decomposition containing the diffusivity
     CALL TRIDSS(kmax_total,nxy, dz(1,ip+1),dz(1,ip+2),dz(1,ip+3), p_c)
     
  ENDIF

! ###################################################################
! Operation
! ###################################################################
!$omp parallel default( shared ) private( ij )
!$omp do
  DO ij = 1,nx*ny*nz
     p_c(ij) = p_c(ij) - p_vel(ij)*p_b(ij) ! diffusivity included in 2.-order derivative
  ENDDO
!$omp end do
!$omp end parallel

! ###################################################################
! -------------------------------------------------------------------
! Transposition
! -------------------------------------------------------------------
#ifdef USE_MPI         
  IF ( ims_npro_k .GT. 1 ) THEN
     CALL DNS_MPI_TRPB_K(p_c, result, ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))
  ENDIF
#endif

  NULLIFY(p_a,p_b,p_c, p_vel)

  ENDIF

  RETURN
END SUBROUTINE OPR_BURGERS_Z
