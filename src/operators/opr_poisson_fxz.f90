#include "types.h"
#include "dns_error.h"
#include "dns_const_mpi.h"

!########################################################################
!# DESCRIPTION
!#
!# Solve Lap p = a using Fourier in xOz planes, to rewritte the problem
!# as 
!#     \hat{p}''-\lambda \hat{p} = \hat{a}
!#
!# where \lambda = kx^2+kz^2
!#
!# The reference value of p at the lower boundary is set to zero
!#
!########################################################################
!# ARGUMENTS 
!#
!# ibc     In    BCs at j1/jmax: 0, for Dirichlet & Dirichlet
!#                               1, for Neumann   & Dirichlet 
!#                               2, for Dirichlet & Neumann   
!#                               3, for Neumann   & Neumann
!#
!# The global variable isize_txc_field defines the size of array txc equal
!# to (imax+2)*(jmax+2)*kmax, or larger if PARALLEL mode
!#
!########################################################################
SUBROUTINE OPR_POISSON_FXZ(flag, nx,ny,nz, g, ibc, &
     a,dpdy, tmp1,tmp2, bcs_hb,bcs_ht, aux, wrk1d,wrk3d)

  USE TLAB_TYPES,    ONLY : grid_dt
  USE TLAB_VARS,     ONLY : isize_txc_dimz
  USE TLAB_VARS,     ONLY : ivfilter, istagger, vfilter_param
#ifdef USE_MPI
  USE TLAB_MPI_VARS, ONLY : ims_offset_i, ims_offset_k
#endif

  IMPLICIT NONE

#include "integers.h"

  LOGICAL,                                  INTENT(IN)    :: flag
  TINTEGER,                                 INTENT(INOUT)    :: nx,ny,nz, ibc
  TYPE(grid_dt),                            INTENT(IN)    :: g(3)
  TREAL,    DIMENSION(nx,ny,nz),            INTENT(INOUT) :: a    ! Forcing term, ans solution field p
  TREAL,    DIMENSION(nx,ny,nz),            INTENT(INOUT) :: dpdy ! Derivative, flag .TRUE.
  TREAL,    DIMENSION(nx,nz),               INTENT(IN)    :: bcs_hb, bcs_ht   ! Boundary-condition fields
  TCOMPLEX, DIMENSION(isize_txc_dimz/2,nz), INTENT(INOUT) :: tmp1,tmp2, wrk3d
  TCOMPLEX, DIMENSION(ny,2),                INTENT(INOUT) :: aux
  TCOMPLEX, DIMENSION(ny,7),                INTENT(INOUT) :: wrk1d

! -----------------------------------------------------------------------
  TINTEGER i, j, k, iglobal, kglobal, ip, isize_line
  TREAL lambda, norm, w1, w2
  TCOMPLEX bcs(3)

! #######################################################################
  norm = C_1_R /M_REAL( g(1)%size *g(3)%size )

  isize_line = nx/2+1

! #######################################################################
! Fourier transform of forcing term; output of this section in array tmp1
! #######################################################################
  IF ( g(3)%size .GT. 1 ) THEN
     CALL OPR_FOURIER_F_X_EXEC(nx,ny,nz, a,bcs_hb,bcs_ht,tmp2, tmp1,wrk3d) 
     CALL OPR_FOURIER_F_Z_EXEC(tmp2,tmp1) ! tmp2 might be overwritten
  ELSE  
     CALL OPR_FOURIER_F_X_EXEC(nx,ny,nz, a,bcs_hb,bcs_ht,tmp1, tmp2,wrk3d)
  ENDIF 

! ###################################################################
! Solve FDE \hat{p}''-\lambda \hat{p} = \hat{a}
! ###################################################################
  DO k = 1,nz
#ifdef USE_MPI
  kglobal = k+ ims_offset_k
#else
  kglobal = k
#endif

  DO i = 1,isize_line
#ifdef USE_MPI
  iglobal = i + ims_offset_i/2
#else
  iglobal = i
#endif

! Define \lambda based on modified wavenumbers (real)
     IF ( g(3)%size .GT. 1 ) THEN; lambda = g(1)%mwn(iglobal,1) + g(3)%mwn(kglobal,1)
     ELSE;                         lambda = g(1)%mwn(iglobal,1); ENDIF

! forcing term
     DO j = 1,ny
        ip = (j-1)*isize_line + i; aux(j,1) = tmp1(ip,k)
     ENDDO

! BCs
     j = ny+1; ip = (j-1)*isize_line + i; bcs(1) = tmp1(ip,k) ! Dirichlet or Neumann
     j = ny+2; ip = (j-1)*isize_line + i; bcs(2) = tmp1(ip,k) ! Dirichlet or Neumann

! -----------------------------------------------------------------------
! Solve for each (kx,kz) a system of 1 complex equation as 2 independent real equations
! -----------------------------------------------------------------------
     SELECT CASE(ibc)

     CASE(3) ! Neumann   & Neumann   BCs
        IF ( istagger .EQ. 0 ) THEN
           IF ( kglobal .EQ. 1             .AND. (iglobal .EQ. 1 .OR. iglobal .EQ. g(1)%size/2+1) .OR.&
                kglobal .EQ. g(3)%size/2+1 .AND. (iglobal .EQ. 1 .OR. iglobal .EQ. g(1)%size/2+1)     )THEN
              CALL FDE_BVP_SINGULAR_NN(g(2)%mode_fdm, ny,i2, &
                    g(2)%jac, aux(1,2),aux(1,1), bcs, wrk1d(1,1), wrk1d(1,3))
           ELSE
              CALL FDE_BVP_REGULAR_NN(g(2)%mode_fdm, ny,i2, lambda, &
                    g(2)%jac, aux(1,2),aux(1,1), bcs, wrk1d(1,1), wrk1d(1,2))
           ENDIF
        ELSE ! In case of staggering only one singular mode + different modified wavenumbers
           IF ( kglobal .EQ. 1 .AND. iglobal .EQ. 1 )THEN
              CALL FDE_BVP_SINGULAR_NN(g(2)%mode_fdm, ny,i2, &
                    g(2)%jac, aux(1,2),aux(1,1), bcs, wrk1d(1,1), wrk1d(1,3))
           ELSE
              CALL FDE_BVP_REGULAR_NN(g(2)%mode_fdm, ny,i2, lambda, &
                    g(2)%jac, aux(1,2),aux(1,1), bcs, wrk1d(1,1), wrk1d(1,2))
           ENDIF
        ENDIF

     CASE(0) ! Dirichlet & Dirichlet BCs
        IF ( istagger .EQ. 0 ) THEN
           IF ( kglobal .EQ. 1             .AND. (iglobal .EQ. 1 .OR. iglobal .EQ. g(1)%size/2+1) .OR.&
                kglobal .EQ. g(3)%size/2+1 .AND. (iglobal .EQ. 1 .OR. iglobal .EQ. g(1)%size/2+1)     )THEN
              CALL FDE_BVP_SINGULAR_DD(g(2)%mode_fdm, ny,i2, &
                    g(2)%nodes,g(2)%jac, aux(1,2),aux(1,1), bcs, wrk1d(1,1), wrk1d(1,3))
           ELSE
              CALL FDE_BVP_REGULAR_DD(g(2)%mode_fdm, ny,i2, lambda, &
                    g(2)%jac, aux(1,2),aux(1,1), bcs, wrk1d(1,1), wrk1d(1,2))
           ENDIF
        ELSE ! In case of staggering only one singular mode + different modified wavenumbers
           IF ( kglobal .EQ. 1 .AND. iglobal .EQ. 1 )THEN
              CALL FDE_BVP_SINGULAR_DD(g(2)%mode_fdm, ny,i2, &
                    g(2)%nodes,g(2)%jac, aux(1,2),aux(1,1), bcs, wrk1d(1,1), wrk1d(1,3))
           ELSE
              CALL FDE_BVP_REGULAR_DD(g(2)%mode_fdm, ny,i2, lambda, &
                    g(2)%jac, aux(1,2),aux(1,1), bcs, wrk1d(1,1), wrk1d(1,2))
           ENDIF
        ENDIF

     END SELECT

   ! Vertical filtering of p and dpdy with a spectral erf filter
     IF (ivfilter .EQ. 1 ) THEN
        CALL FILTER_VERTICAL_1D(aux(:,2),   ny, vfilter_param, aux(1,1))
        CALL FILTER_VERTICAL_1D(wrk1d(:,1), ny, vfilter_param, aux(1,1))
     ENDIF
 
   ! Normalize
     DO j = 1,ny
        ip = (j-1)*isize_line + i
        tmp1(ip,k) = aux(j,2)  *norm ! solution
        tmp2(ip,k) = wrk1d(j,1)*norm ! Oy derivative
     ENDDO
      
   ENDDO
   ENDDO
 

! ###################################################################
! Fourier field p (based on array tmp1)
! ###################################################################
  IF ( g(3)%size .GT. 1 ) THEN
     CALL OPR_FOURIER_B_Z_EXEC(tmp1,wrk3d) 
     CALL OPR_FOURIER_B_X_EXEC(nx,ny,nz, wrk3d,a, tmp1) 
  ELSE
     CALL OPR_FOURIER_B_X_EXEC(nx,ny,nz, tmp1,a, wrk3d) 
  ENDIF
     
! ###################################################################
! Fourier derivatives (based on array tmp2)
! ###################################################################
  IF ( flag ) THEN
     IF ( g(3)%size .GT. 1 ) THEN
        CALL OPR_FOURIER_B_Z_EXEC(tmp2,wrk3d) 
        CALL OPR_FOURIER_B_X_EXEC(nx,ny,nz, wrk3d,dpdy, tmp2) 
     ELSE
        CALL OPR_FOURIER_B_X_EXEC(nx,ny,nz, tmp2,dpdy, wrk3d) 
     ENDIF
  ENDIF

  RETURN
END SUBROUTINE OPR_POISSON_FXZ

!########################################################################
!# DESCRIPTTION
!# 
!# Filter function
!#
!########################################################################
SUBROUTINE FILTER_VERTICAL_1D(a, n, lcut, wrk1)
  IMPLICIT NONE
  
#include "integers.h"
#ifdef USE_FFTW
#include "fftw3.f"
#endif 
  
  TCOMPLEX, DIMENSION(n), INTENT(INOUT) :: a       ! 1d-array to be filtered
  TINTEGER,               INTENT(INOUT) :: n       ! ny
  TREAL,                  INTENT(INOUT) :: lcut    ! filter parameter
  TCOMPLEX, DIMENSION(n), INTENT(INOUT) :: wrk1    ! aux 1d-array
  
! -------------------------------------------------------------------

  INTEGER(8)                            :: plan_f, plan_b  
  TREAL                                 :: norm 
  
! ###################################################################
  norm = C_1_R / M_REAL(n)
  
! Build plans
  CALL dfftw_plan_dft_1d(plan_f, n, a,    wrk1, FFTW_FORWARD,  FFTW_ESTIMATE)
  CALL dfftw_plan_dft_1d(plan_b, n, wrk1, a,    FFTW_BACKWARD, FFTW_ESTIMATE)

! Execute + Filtering
  CALL dfftw_execute_dft(plan_f, a,    wrk1)
  CALL FILTER_ERF_1D(wrk1, n, lcut)
  CALL dfftw_execute_dft(plan_b, wrk1, a   )

! Normalize
  a = a * norm

! Destroy plans  
  CALL dfftw_destroy_plan(plan_f)
  CALL dfftw_destroy_plan(plan_b)

  RETURN
END SUBROUTINE FILTER_VERTICAL_1D
!########################################################################
SUBROUTINE FILTER_ERF_1D(a, n, lcut)

  IMPLICIT NONE
 
#include "integers.h"
 
  TCOMPLEX, DIMENSION(n), INTENT(INOUT) :: a       ! 1d-array to be filtered
  TINTEGER,               INTENT(IN   ) :: n 
  TREAL,                  INTENT(IN   ) :: lcut    ! filter parameter
 
 ! -------------------------------------------------------------------
   
  TREAL                                 :: k_ref, k_rel
  TINTEGER                              :: i, j, n_ny
 
 ! ###################################################################
  n_ny  = n/2+1
  k_ref = C_2_R * n / lcut
 
  DO i=1,n_ny
    k_rel = (i-1) / k_ref
    a(i)  = a(i) * (ERF(C_8_R * (C_1_R - k_rel)) + C_1_R) / C_2_R
   !  write(*,*)'i,    krel = ', i, j, k_rel
  ENDDO

!   write(*,*)'======================================================='

  i = n_ny+1

  DO j=n_ny-1,n-(n-2),-1 
    k_rel = (j-1) / k_ref
    a(i)  = a(i) * (ERF(C_8_R * (C_1_R - k_rel)) + C_1_R) / C_2_R
   !  write(*,*)'i, j, krel = ', i, j, k_rel
    i = i + 1
  ENDDO

  RETURN
END SUBROUTINE FILTER_ERF_1D