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
!# ARGUMENTS 
!#
!# bcs2_imin  In   BC derivative at imin: 0  biased, non-zero
!#                                        1  forced to zero
!# bcs2_imax  In   BC derivative at imax: 0  biased, non-zero
!#                                        1  forced to zero
!#
!########################################################################
SUBROUTINE OPR_BURGERS(imode_fdm, is, nyz, g, s,u, &
     bcs1_imin,bcs1_imax,bcs2_imin,bcs2_imax, result, wrk2d,wrk3d)

  USE DNS_TYPES,     ONLY : grid_structure
  USE DNS_CONSTANTS, ONLY : efile

  IMPLICIT NONE

  TINTEGER,                     INTENT(IN)    :: is, nyz
  TYPE(grid_structure),         INTENT(IN)    :: g
  TREAL, DIMENSION(nyz*g%size), INTENT(IN)    :: s,u    ! argument field and velocity field
  TREAL, DIMENSION(nyz*g%size), INTENT(OUT)   :: result ! result N(u) applied to s
  TREAL, DIMENSION(nyz*g%size), INTENT(INOUT) :: wrk3d  ! dsdx
  TREAL, DIMENSION(*),          INTENT(INOUT) :: wrk2d

  TINTEGER imode_fdm
  TINTEGER bcs1_imin,bcs1_imax, bcs2_imin,bcs2_imax

! -------------------------------------------------------------------
  TINTEGER iunif, ip,ij

! ###################################################################
  IF ( bcs2_imin + bcs2_imax .GT. 0 ) THEN
     CALL IO_WRITE_ASCII(efile,'OPR_BURGERS. Only developed for biased BCs.')
     CALL DNS_STOP(DNS_ERROR_UNDEVELOP)
  ENDIF

  IF ( g%uniform ) THEN; iunif = 0;
  ELSE;                  iunif = 1; ENDIF
       
! -------------------------------------------------------------------
! Periodic case
! -------------------------------------------------------------------
  IF ( g%periodic ) THEN 
     IF      ( imode_fdm .eq. FDM_COM4_JACOBIAN ) THEN
        CALL FDM_C1N4P_RHS(g%size,nyz, s, wrk3d)
     ELSE IF ( imode_fdm .eq. FDM_COM6_JACOBIAN ) THEN
        CALL FDM_C1N6P_RHS(g%size,nyz, s, wrk3d)
     ELSE IF ( imode_fdm .eq. FDM_COM8_JACOBIAN ) THEN
        CALL FDM_C1N8P_RHS(g%size,nyz, s, wrk3d)
     ELSE IF ( imode_fdm .EQ. FDM_COM6_DIRECT   ) THEN
        CALL FDM_C1N6P_RHS(g%size,nyz, s, wrk3d)
     ENDIF
     ! ip = inb_grid_1 - 1
     ! CALL TRIDPSS(g%size,nyz, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3),&
     !      dx(1,ip+4),dx(1,ip+5), wrk3d,wrk2d)
     CALL TRIDPSS(g%size,nyz, g%lu1(1,1),g%lu1(1,2),g%lu1(1,3),g%lu1(1,4),g%lu1(1,5), wrk3d,wrk2d)

! Calculate second derivative
     IF      ( imode_fdm .eq. FDM_COM4_JACOBIAN ) THEN
        CALL FDM_C2N4P_RHS(g%size,nyz, s, result)
     ELSE IF ( imode_fdm .eq. FDM_COM6_JACOBIAN ) THEN 
        CALL FDM_C2N6P_RHS(g%size,nyz, s, result)
     ELSE IF ( imode_fdm .eq. FDM_COM8_JACOBIAN ) THEN 
        CALL FDM_C2N6P_RHS(g%size,nyz, s, result) !8th not yet developed
     ELSE IF ( imode_fdm .EQ. FDM_COM6_DIRECT   ) THEN
        CALL FDM_C2N6P_RHS(g%size,nyz, s, result)
     ENDIF
     ! ip = inb_grid_3 + is*5 - 1 ! LU decomposition containing the diffusivity
     ! CALL TRIDPSS(g%size,nyz, dx(1,ip+1),dx(1,ip+2),&
     !      dx(1,ip+3),dx(1,ip+4),dx(1,ip+5), result,wrk2d)
     ip = is*5 ! LU decomposition containing the diffusivity
     CALL TRIDPSS(g%size,nyz, g%lu2d(1,ip+1),g%lu2d(1,ip+2),g%lu2d(1,ip+3),g%lu2d(1,ip+4),g%lu2d(1,ip+5), result,wrk2d)

! -------------------------------------------------------------------
! Nonperiodic case
! -------------------------------------------------------------------
  ELSE
     IF      ( imode_fdm .eq. FDM_COM4_JACOBIAN ) THEN;
        CALL FDM_C1N4_RHS(g%size,nyz, bcs1_imin,bcs1_imax, s, wrk3d)
     ELSE IF ( imode_fdm .eq. FDM_COM6_JACOBIAN ) THEN; 
        CALL FDM_C1N6_RHS(g%size,nyz, bcs1_imin,bcs1_imax, s, wrk3d)
     ELSE IF ( imode_fdm .eq. FDM_COM8_JACOBIAN ) THEN; 
        CALL FDM_C1N8_RHS(g%size,nyz, bcs1_imin,bcs1_imax, s, wrk3d)
     ELSE IF ( imode_fdm .eq. FDM_COM6_DIRECT   ) THEN; 
        CALL FDM_C1N6_RHS(g%size,nyz, bcs1_imin,bcs1_imax, s, wrk3d) ! not yet implemented
     ENDIF
     ! ip = inb_grid_1 + (bcs1_imin + bcs1_imax*2)*3 - 1
     ! CALL TRIDSS(g%size,nyz, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3), wrk3d)
     ip = (bcs1_imin + bcs1_imax*2)*3 
     CALL TRIDSS(g%size,nyz, g%lu1(1,ip+1),g%lu1(1,ip+2),g%lu1(1,ip+3), wrk3d)

! Calculate second derivative
     IF      ( imode_fdm .eq. FDM_COM4_JACOBIAN ) THEN; 
        CALL FDM_C2N4_RHS(iunif, g%size,nyz, bcs2_imin,bcs2_imax, g%aux, s, wrk3d, result)
     ELSE IF ( imode_fdm .eq. FDM_COM6_JACOBIAN ) THEN; 
        CALL FDM_C2N6_RHS(iunif, g%size,nyz, bcs2_imin,bcs2_imax, g%aux, s, wrk3d, result)
     ELSE IF ( imode_fdm .eq. FDM_COM8_JACOBIAN ) THEN; 
        CALL FDM_C2N6_RHS(iunif, g%size,nyz, bcs2_imin,bcs2_imax, g%aux, s, wrk3d, result)
     ELSE IF ( imode_fdm .eq. FDM_COM6_DIRECT   ) THEN; 
!        CALL FDM_C2N6N_RHS(g%size,nyz, dx(1,inb_grid_2+3), &
        CALL FDM_C2N6N_RHS(g%size,nyz, g%lu2(1,4), s, result)
     ENDIF
     ! ip = inb_grid_3 + is*5 - 1 ! LU decomposition containing the diffusivity
     ! CALL TRIDSS(g%size,nyz, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3), result)
     ip = is*5 ! LU decomposition containing the diffusivity
     CALL TRIDSS(g%size,nyz, g%lu2d(1,ip+1),g%lu2d(1,ip+2),g%lu2d(1,ip+3), result)

  ENDIF

! ###################################################################
! Operation
! ###################################################################
!$omp parallel default( shared ) private( ij )
!$omp do
  DO ij = 1,nyz*g%size
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
SUBROUTINE OPR_BURGERS_X(ivel, is, imode_fdm, nx,ny,nz, g, s,u1,u2, result, &
     bcs1_imin,bcs1_imax, bcs2_imin,bcs2_imax, tmp1, wrk2d,wrk3d)

  USE DNS_TYPES, ONLY : grid_structure
#ifdef USE_MPI
  USE DNS_MPI
#endif

  IMPLICIT NONE

  TINTEGER ivel, is, imode_fdm, nx,ny,nz
  TINTEGER bcs1_imin,bcs1_imax, bcs2_imin,bcs2_imax
  TYPE(grid_structure),           INTENT(IN)    :: g
  TREAL, DIMENSION(nx*ny*nz),     INTENT(IN)    :: s,u1,u2
  TREAL, DIMENSION(nx*ny*nz),     INTENT(OUT)   :: result
  TREAL, DIMENSION(nx*ny*nz),     INTENT(INOUT) :: tmp1, wrk3d
  TREAL, DIMENSION(ny*nz),        INTENT(INOUT) :: wrk2d

  TARGET s,u1,u2, tmp1, result, wrk3d

! -------------------------------------------------------------------
  TINTEGER nyz
  TREAL, DIMENSION(:), POINTER :: p_a,p_b,p_c,p_d, p_vel
#ifdef USE_MPI
  TINTEGER, PARAMETER :: id = DNS_MPI_I_PARTIAL
#endif

! ###################################################################
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
  CALL DGETMO       (p_a, g%size, g%size, nyz,        p_b, nyz)
#else
  CALL DNS_TRANSPOSE(p_a, g%size, nyz,        g%size, p_b, nyz)
#endif

! ###################################################################
  CALL OPR_BURGERS(imode_fdm, is, nyz, g, p_b, p_vel, &
       bcs1_imin,bcs1_imax, bcs2_imin,bcs2_imax, p_d, wrk2d,p_c)
  
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
SUBROUTINE OPR_BURGERS_Y(ivel, is, imode_fdm, nx,ny,nz, g, s,u1,u2, result, &
     bcs1_jmin,bcs1_jmax, bcs2_jmin,bcs2_jmax, tmp1, wrk2d,wrk3d)

  USE DNS_TYPES, ONLY : grid_structure
  IMPLICIT NONE

  TINTEGER ivel, is, imode_fdm, nx,ny,nz
  TINTEGER bcs1_jmin,bcs1_jmax, bcs2_jmin,bcs2_jmax
  TYPE(grid_structure),           INTENT(IN)    :: g
  TREAL, DIMENSION(nx*ny*nz),     INTENT(IN)    :: s,u1,u2
  TREAL, DIMENSION(nx*ny*nz),     INTENT(OUT)   :: result
  TREAL, DIMENSION(nx*ny*nz),     INTENT(INOUT) :: tmp1, wrk3d
  TREAL, DIMENSION(nx*nz),        INTENT(INOUT) :: wrk2d

  TARGET s,u1,u2, tmp1, result, wrk3d

! -------------------------------------------------------------------
  TINTEGER nxy, nxz
  TREAL, DIMENSION(:),   POINTER :: p_org, p_dst1, p_dst2, p_vel

! ###################################################################
  IF ( g%size .EQ. 1 ) THEN ! Set to zero in 2D case
     result = C_0_R
     
  ELSE
! ###################################################################
  nxy = nx*ny 
  nxz = nx*nz

! -------------------------------------------------------------------
! Local transposition: Make y direction the last one
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
  CALL OPR_BURGERS(imode_fdm, is, nxz, g, p_org, p_vel, &
       bcs1_jmin,bcs1_jmax, bcs2_jmin,bcs2_jmax, p_dst2, wrk2d,p_dst1)
  
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
SUBROUTINE OPR_BURGERS_Z(ivel, is, imode_fdm, nx,ny,nz, g, s,u1,u2, result, &
     bcs1_kmin,bcs1_kmax, bcs2_kmin,bcs2_kmax, tmp1, wrk2d,wrk3d)

  USE DNS_TYPES, ONLY : grid_structure
#ifdef USE_MPI
  USE DNS_MPI
#endif

  IMPLICIT NONE

  TINTEGER ivel, is, imode_fdm, nx,ny,nz
  TINTEGER bcs1_kmin, bcs1_kmax, bcs2_kmin, bcs2_kmax
  TYPE(grid_structure),           INTENT(IN)    :: g
  TREAL, DIMENSION(nx*ny*nz),     INTENT(IN)    :: s,u1,u2
  TREAL, DIMENSION(nx*ny*nz),     INTENT(OUT)   :: result
  TREAL, DIMENSION(nx*ny*nz),     INTENT(INOUT) :: tmp1, wrk3d
  TREAL, DIMENSION(nx*ny),        INTENT(INOUT) :: wrk2d

  TARGET s,u1,u2, tmp1, result, wrk3d

! -------------------------------------------------------------------
  TINTEGER nxy
  TREAL, DIMENSION(:),   POINTER :: p_a,p_b,p_c, p_vel
#ifdef USE_MPI
  TINTEGER, PARAMETER : id  = DNS_MPI_K_PARTIAL
#endif

! ###################################################################
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
  CALL OPR_BURGERS(imode_fdm, is, nxy, g, p_a, p_vel, &
       bcs1_kmin,bcs1_kmax, bcs2_kmin,bcs2_kmax, p_c, wrk2d,p_b)

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
