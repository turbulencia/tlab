#include "types.h"
#include "dns_const.h"
#include "dns_error.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

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
!# along generic direction x to nyz lines of data
!#
!########################################################################
SUBROUTINE OPR_BURGERS(is, nlines, g, s,u, result, bcs_min,bcs_max, wrk2d,wrk3d)

  USE DNS_TYPES,     ONLY : grid_structure
  USE DNS_CONSTANTS, ONLY : efile
  
  IMPLICIT NONE

  TINTEGER,                        INTENT(IN)    :: is
  TINTEGER,                        INTENT(IN)    :: nlines     ! # of lines to be solved
  TYPE(grid_structure),            INTENT(IN)    :: g
  TREAL, DIMENSION(nlines*g%size), INTENT(IN)    :: s,u        ! argument field and velocity field
  TREAL, DIMENSION(nlines*g%size), INTENT(OUT)   :: result     ! N(u) applied to s
  TINTEGER,                        INTENT(IN)    :: bcs_min(2) ! BC derivative: 0 biased, non-zero
  TINTEGER,                        INTENT(IN)    :: bcs_max(2) !                1 forced to zero
  TREAL, DIMENSION(nlines*g%size), INTENT(INOUT) :: wrk3d      ! dsdx
  TREAL, DIMENSION(*),             INTENT(INOUT) :: wrk2d

! -------------------------------------------------------------------
  TINTEGER ip,ij

! ###################################################################
  IF ( bcs_min(2) + bcs_max(2) .GT. 0 ) THEN
     CALL IO_WRITE_ASCII(efile,'OPR_BURGERS. Only developed for biased BCs.')
     CALL DNS_STOP(DNS_ERROR_UNDEVELOP)
  ENDIF

! First derivative
  CALL OPR_PARTIAL1(nlines, g, s,wrk3d, bcs_min(1),bcs_max(1), wrk2d)

! Second derivative uses LE decomposition including diffusivity coefficient
! -------------------------------------------------------------------
! Periodic case
! -------------------------------------------------------------------
  IF ( g%periodic ) THEN 
     SELECT CASE( g%mode_fdm )
        
     CASE( FDM_COM4_JACOBIAN )
        CALL FDM_C2N4P_RHS(g%size,nlines, s, result)
        
     CASE( FDM_COM6_JACOBIAN, FDM_COM6_DIRECT ) ! Direct = Jacobian because uniform grid
        CALL FDM_C2N6P_RHS(g%size,nlines, s, result)
        
     CASE( FDM_COM8_JACOBIAN )                  ! Not yet implemented; default to 6. order
        CALL FDM_C2N6P_RHS(g%size,nlines, s, result)
        
     END SELECT

     ip = is*5 ! LU decomposition containing the diffusivity
     CALL TRIDPSS(g%size,nlines, g%lu2d(1,ip+1),g%lu2d(1,ip+2),g%lu2d(1,ip+3),g%lu2d(1,ip+4),g%lu2d(1,ip+5), result,wrk2d)

! -------------------------------------------------------------------
! Nonperiodic case
! -------------------------------------------------------------------
  ELSE
     SELECT CASE( g%mode_fdm )
        
     CASE( FDM_COM4_JACOBIAN )
        CALL FDM_C2N4_RHS(g%uniform, g%size,nlines, bcs_min(2),bcs_max(2), g%jac, s, wrk3d, result)

     CASE( FDM_COM6_JACOBIAN )
        CALL FDM_C2N6_RHS(g%uniform, g%size,nlines, bcs_min(2),bcs_max(2), g%jac, s, wrk3d, result)

     CASE( FDM_COM8_JACOBIAN ) ! Not yet implemented; defaulting to 6. order
        CALL FDM_C2N6_RHS(g%uniform, g%size,nlines, bcs_min(2),bcs_max(2), g%jac, s, wrk3d, result)

     CASE( FDM_COM6_DIRECT   )
        CALL FDM_C2N6N_RHS(g%size,nlines, g%lu2(1,4), u, result)

     END SELECT

     ip = is*3 ! LU decomposition containing the diffusivity
     CALL TRIDSS(g%size,nlines, g%lu2d(1,ip+1),g%lu2d(1,ip+2),g%lu2d(1,ip+3), result)

  ENDIF

! ###################################################################
! Operation
! ###################################################################
!$omp parallel default( shared ) private( ij )
!$omp do
  DO ij = 1,nlines*g%size
     result(ij) = result(ij) - u(ij)*wrk3d(ij) ! diffusivity included in 2.-order derivative
  ENDDO
!$omp end do
!$omp end parallel

  RETURN
END SUBROUTINE OPR_BURGERS

!########################################################################
!# Routines for different specific directions:
!#
!# ivel       In   Flag indicating the array containing the velocity:
!#                    0 for velocity being the scalar itself
!#                    1 for velocity passed through u1, or u2 if transposed required
!# is         In   Scalar index; if 0, then velocity
!# tmp1       Out  Transpose velocity
!########################################################################
SUBROUTINE OPR_BURGERS_X(ivel, is, nx,ny,nz, g, s,u1,u2, result, &
     bcs1_imin,bcs1_imax, bcs2_imin,bcs2_imax, tmp1, wrk2d,wrk3d)

  USE DNS_TYPES, ONLY : grid_structure
#ifdef USE_MPI
  USE DNS_MPI
#endif

  IMPLICIT NONE

  TINTEGER ivel, is, nx,ny,nz
  TINTEGER bcs1_imin,bcs1_imax, bcs2_imin,bcs2_imax
  TYPE(grid_structure),           INTENT(IN)    :: g
  TREAL, DIMENSION(nx*ny*nz),     INTENT(IN)    :: s,u1,u2
  TREAL, DIMENSION(nx*ny*nz),     INTENT(OUT)   :: result
  TREAL, DIMENSION(nx*ny*nz),     INTENT(INOUT) :: tmp1, wrk3d
  TREAL, DIMENSION(ny*nz),        INTENT(INOUT) :: wrk2d

  TARGET s,u1,u2, tmp1, result, wrk3d

! -------------------------------------------------------------------
  TINTEGER nyz, bcs_min(2), bcs_max(2)
  TREAL, DIMENSION(:), POINTER :: p_a,p_b,p_c,p_d, p_vel
#ifdef USE_MPI
  TINTEGER, PARAMETER :: id = DNS_MPI_I_PARTIAL
#endif

! ###################################################################
  bcs_min(1) = bcs1_imin; bcs_max(1) = bcs1_imax
  bcs_min(2) = bcs2_imin; bcs_max(2) = bcs2_imax

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
! Local transposition: make x-direction the last one
! -------------------------------------------------------------------
#ifdef USE_ESSL
  CALL DGETMO       (p_a, g%size, g%size, nyz,    p_b, nyz)
#else
  CALL DNS_TRANSPOSE(p_a, g%size, nyz,    g%size, p_b, nyz)
#endif

! ###################################################################
  CALL OPR_BURGERS(is, nyz, g, p_b, p_vel, p_d, bcs_min,bcs_max, wrk2d,p_c)
  
! ###################################################################
! Put arrays back in the order in which they came in
#ifdef USE_ESSL
  CALL DGETMO       (p_d, nyz, nyz,        g%size, p_c, g%size)
#else
  CALL DNS_TRANSPOSE(p_d, nyz, g%size, nyz,        p_c, g%size)
#endif

#ifdef USE_MPI         
  IF ( ims_npro_i .GT. 1 ) THEN
     CALL DNS_MPI_TRPB_I(p_c, result, ims_ds_i(1,id), ims_dr_i(1,id), ims_ts_i(1,id), ims_tr_i(1,id))
  ENDIF
#endif

  NULLIFY(p_a,p_b,p_c,p_d, p_vel)

  RETURN
END SUBROUTINE OPR_BURGERS_X

!########################################################################
!########################################################################
SUBROUTINE OPR_BURGERS_Y(ivel, is, nx,ny,nz, g, s,u1,u2, result, &
     bcs1_jmin,bcs1_jmax, bcs2_jmin,bcs2_jmax, tmp1, wrk2d,wrk3d)

  USE DNS_TYPES, ONLY : grid_structure
  IMPLICIT NONE

  TINTEGER ivel, is, nx,ny,nz
  TINTEGER bcs1_jmin,bcs1_jmax, bcs2_jmin,bcs2_jmax
  TYPE(grid_structure),           INTENT(IN)    :: g
  TREAL, DIMENSION(nx*ny*nz),     INTENT(IN)    :: s,u1,u2
  TREAL, DIMENSION(nx*ny*nz),     INTENT(OUT)   :: result
  TREAL, DIMENSION(nx*ny*nz),     INTENT(INOUT) :: tmp1, wrk3d
  TREAL, DIMENSION(nx*nz),        INTENT(INOUT) :: wrk2d

  TARGET s,u1,u2, tmp1, result, wrk3d

! -------------------------------------------------------------------
  TINTEGER nxy, nxz, bcs_min(2), bcs_max(2)
  TREAL, DIMENSION(:),   POINTER :: p_org, p_dst1, p_dst2, p_vel

! ###################################################################
  bcs_min(1) = bcs1_jmin; bcs_max(1) = bcs1_jmax
  bcs_min(2) = bcs2_jmin; bcs_max(2) = bcs2_jmax

  IF ( g%size .EQ. 1 ) THEN ! Set to zero in 2D case
     result = C_0_R
     
  ELSE
! ###################################################################
  nxy = nx*ny 
  nxz = nx*nz

! -------------------------------------------------------------------
! Local transposition: Make y-direction the last one
! -------------------------------------------------------------------
  IF ( nz .EQ. 1 ) THEN
     p_org  => s
     p_dst1 => tmp1
     p_dst2 => result
  ELSE
#ifdef USE_ESSL
     CALL DGETMO       (s, nxy, nxy, nz, tmp1, nz)
#else
     CALL DNS_TRANSPOSE(s, nxy, nz, nxy, tmp1, nz)
#endif
     p_org  => tmp1
     p_dst1 => result
     p_dst2 => wrk3d
  ENDIF

! pointer to velocity
  IF ( ivel .EQ. 0 ) THEN; p_vel => p_org
  ELSE;
     IF ( nz .EQ. 1 ) THEN; p_vel => u1         ! I do not need the transposed
     ELSE;                  p_vel => u2; ENDIF  ! I do     need the transposed
  ENDIF

! ###################################################################
  CALL OPR_BURGERS(is, nxz, g, p_org, p_vel, p_dst2, bcs_min,bcs_max, wrk2d,p_dst1)
  
! ###################################################################
! Put arrays back in the order in which they came in
  IF ( nz .GT. 1 ) THEN
#ifdef USE_ESSL
     CALL DGETMO       (p_dst2, nz, nz, nxy, result, nxy)
#else
     CALL DNS_TRANSPOSE(p_dst2, nz, nxy, nz, result, nxy)
#endif
  ENDIF

  NULLIFY(p_org,p_dst1,p_dst2, p_vel)

  ENDIF

  RETURN
END SUBROUTINE OPR_BURGERS_Y

!########################################################################
!########################################################################
SUBROUTINE OPR_BURGERS_Z(ivel, is, nx,ny,nz, g, s,u1,u2, result, &
     bcs1_kmin,bcs1_kmax, bcs2_kmin,bcs2_kmax, tmp1, wrk2d,wrk3d)

  USE DNS_TYPES, ONLY : grid_structure
#ifdef USE_MPI
  USE DNS_MPI
#endif

  IMPLICIT NONE

  TINTEGER ivel, is, nx,ny,nz
  TINTEGER bcs1_kmin, bcs1_kmax, bcs2_kmin, bcs2_kmax
  TYPE(grid_structure),           INTENT(IN)    :: g
  TREAL, DIMENSION(nx*ny*nz),     INTENT(IN)    :: s,u1,u2
  TREAL, DIMENSION(nx*ny*nz),     INTENT(OUT)   :: result
  TREAL, DIMENSION(nx*ny*nz),     INTENT(INOUT) :: tmp1, wrk3d
  TREAL, DIMENSION(nx*ny),        INTENT(INOUT) :: wrk2d

  TARGET s,u1,u2, tmp1, result, wrk3d

! -------------------------------------------------------------------
  TINTEGER nxy, bcs_min(2), bcs_max(2)  
  TREAL, DIMENSION(:), POINTER :: p_a,p_b,p_c, p_vel
#ifdef USE_MPI
  TINTEGER, PARAMETER :: id  = DNS_MPI_K_PARTIAL
#endif

! ###################################################################
  bcs_min(1) = bcs1_kmin; bcs_max(1) = bcs1_kmax
  bcs_min(2) = bcs2_kmin; bcs_max(2) = bcs2_kmax

  IF ( g%size .EQ. 1 ) THEN ! Set to zero in 2D case
     result = C_0_R
     
  ELSE
! ###################################################################
! -------------------------------------------------------------------
! MPI transposition
! -------------------------------------------------------------------
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
  CALL OPR_BURGERS(is, nxy, g, p_a, p_vel, p_c, bcs_min,bcs_max, wrk2d,p_b)

! ###################################################################
! Put arrays back in the order in which they came in
#ifdef USE_MPI         
  IF ( ims_npro_k .GT. 1 ) THEN
     CALL DNS_MPI_TRPB_K(p_c, result, ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))
  ENDIF
#endif

  NULLIFY(p_a,p_b,p_c, p_vel)

  ENDIF

  RETURN
END SUBROUTINE OPR_BURGERS_Z
