#include "types.h"
#include "dns_const.h"
#include "dns_error.h"
#include "dns_const_mpi.h"

!########################################################################
!# HISTORY
!#
!# 2010/11/11 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Solve Lap p + \alpha p = a using Fourier in xOz planes, to rewritte
!# the problem as
!#
!#      \hat{p}''-(\lambda-\alpha) \hat{p} = \hat{a}
!#
!# where \lambda = kx^2+kz^2
!#
!########################################################################
!# ARGUMENTS
!#
!# a       In    Forcing term
!#         Out   Solution field p
!# ibc     In    BCs at j1/jmax: 0, for Dirichlet & Dirichlet
!#                               1, for Neumann   & Dirichlet
!#                               2, for Dirichlet & Neumann
!#                               3, for Neumann   & Neumann
!# bcs_??  In    BCs at j1/jmax
!#
!# The global variable isize_txc_field defines the size of array txc equal
!# to (imax+2)*(jmax+2)*kmax, or larger if PARALLEL mode
!#
!########################################################################
SUBROUTINE OPR_HELMHOLTZ_FXZ(nx,ny,nz, g, ibc, alpha,&
     a, tmp1,tmp2, bcs_hb,bcs_ht, aux, wrk1d,wrk3d)

  USE TLAB_CONSTANTS, ONLY : efile
  USE TLAB_TYPES,     ONLY : grid_dt
  USE TLAB_VARS,    ONLY : isize_txc_dimz
  USE TLAB_PROCS
#ifdef USE_MPI
  USE TLAB_MPI_VARS, ONLY : ims_offset_i, ims_offset_k
#endif

  IMPLICIT NONE

#include "integers.h"

  TINTEGER nx,ny,nz, ibc
  TYPE(grid_dt),      INTENT(IN)    :: g(3)
  TREAL alpha
  TREAL,    DIMENSION(nx,ny,nz)            :: a
  TREAL,    DIMENSION(nx,nz)               :: bcs_hb, bcs_ht
  TCOMPLEX, DIMENSION(isize_txc_dimz/2,nz) :: tmp1,tmp2, wrk3d
  TCOMPLEX, DIMENSION(ny,2)                :: aux
  TCOMPLEX, DIMENSION(ny,7)                :: wrk1d

! -----------------------------------------------------------------------
  TINTEGER i,j,k, iglobal,kglobal, ip, isize_line
  TREAL lambda, norm
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
! Solve FDE (\hat{p}')'-\lambda \hat{p} = \hat{a}
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

     lambda = lambda - alpha

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
     IF      ( ibc .EQ. 3 ) THEN ! Neumman BCs
        CALL FDE_BVP_REGULAR_NN(g(2)%mode_fdm, ny,i2, lambda, g(2)%jac, aux(1,2),aux(1,1), bcs, wrk1d(1,1), wrk1d(1,2))
     ELSE IF ( ibc .EQ. 0 ) THEN ! Dirichlet BCs
        CALL FDE_BVP_REGULAR_DD(g(2)%mode_fdm, ny,i2, lambda, g(2)%jac, aux(1,2),aux(1,1), bcs, wrk1d(1,1), wrk1d(1,2))
     ELSE
        CALL TLAB_WRITE_ASCII(efile,'OPR_HELMHOLT_FXZ. Undeveloped BCs.')
        CALL TLAB_STOP(DNS_ERROR_UNDEVELOP)
     ENDIF

! normalize
     DO j = 1,ny
        ip = (j-1)*isize_line + i
        tmp1(ip,k) = aux(j,2)  *norm ! solution
!        tmp2(ip,k) = wrk1d(j,1)*norm ! Oy derivative; not used in this routine
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

  RETURN
END SUBROUTINE OPR_HELMHOLTZ_FXZ

!########################################################################
!########################################################################
! Same, but using the direct mode of FDM
SUBROUTINE OPR_HELMHOLTZ_FXZ_2(nx,ny,nz, g, ibc, alpha,&
     a, tmp1,tmp2, bcs_hb,bcs_ht, aux, wrk1d,wrk3d)

  USE TLAB_CONSTANTS, ONLY : efile
  USE TLAB_TYPES,     ONLY : grid_dt
  USE TLAB_VARS,    ONLY : isize_field, isize_txc_dimz
  USE TLAB_PROCS
#ifdef USE_MPI
  USE TLAB_MPI_VARS, ONLY : ims_offset_i, ims_offset_k
#endif

  IMPLICIT NONE

#include "integers.h"

  TINTEGER nx,ny,nz, ibc
  TREAL alpha
  TYPE(grid_dt),      INTENT(IN)    :: g(3)
  TREAL,    DIMENSION(isize_field)         :: a
  TREAL,    DIMENSION(nx,nz)               :: bcs_hb, bcs_ht
  TCOMPLEX, DIMENSION(isize_txc_dimz/2,nz) :: tmp1,tmp2, wrk3d
  TCOMPLEX, DIMENSION(ny,2)                :: aux
  TREAL,    DIMENSION(ny,7)                :: wrk1d

! -----------------------------------------------------------------------
  TINTEGER i,j,k, iglobal,kglobal, ip, isize_line
  TREAL lambda, norm
  TCOMPLEX bcs(3)

! #######################################################################
  norm = C_1_R/M_REAL(g(1)%size*g(3)%size)

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
     IF ( g(3)%size .GT. 1 ) THEN; lambda = g(1)%mwn(iglobal,2) + g(3)%mwn(kglobal,2)
     ELSE;                         lambda = g(1)%mwn(iglobal,2); ENDIF

     lambda = lambda - alpha

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
     IF ( ibc .EQ. 0 ) THEN ! Dirichlet BCs
        IF      ( g(2)%mode_fdm .EQ. FDM_COM6_JACOBIAN .OR. g(2)%mode_fdm .EQ. FDM_COM6_JACPENTA ) THEN
           CALL INT_C2N6_LHS_E(ny,    g(2)%jac, lambda, &
                wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5), wrk1d(1,6),wrk1d(1,7))
           CALL INT_C2N6_RHS  (ny,i2, g(2)%jac, aux(1,1),aux(1,2))
        ELSE IF ( g(2)%mode_fdm .EQ. FDM_COM6_DIRECT   ) THEN
           wrk1d = C_0_R
           CALL INT_C2N6N_LHS_E(ny,    g(2)%lu2(1,8), g(2)%lu2(1,4), lambda, &
                wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5), wrk1d(1,6),wrk1d(1,7))
           CALL INT_C2N6N_RHS  (ny,i2, g(2)%lu2(1,8), aux(1,1),aux(1,2))
        ENDIF

        CALL PENTADFS(ny-2,    wrk1d(2,1),wrk1d(2,2),wrk1d(2,3),wrk1d(2,4),wrk1d(2,5))

        CALL PENTADSS(ny-2,i1, wrk1d(2,1),wrk1d(2,2),wrk1d(2,3),wrk1d(2,4),wrk1d(2,5), wrk1d(2,6))
        CALL PENTADSS(ny-2,i1, wrk1d(2,1),wrk1d(2,2),wrk1d(2,3),wrk1d(2,4),wrk1d(2,5), wrk1d(2,7))

        CALL PENTADSS(ny-2,i2, wrk1d(2,1),wrk1d(2,2),wrk1d(2,3),wrk1d(2,4),wrk1d(2,5), aux(2,2))

        aux(:,2) = aux(:,2) +  bcs(1)*wrk1d(:,6) + bcs(2)*wrk1d(:,7)

     ELSE
        CALL TLAB_WRITE_ASCII(efile,'OPR_HELMHOLT_FXZ_2. Undeveloped BCs.')
        CALL TLAB_STOP(DNS_ERROR_UNDEVELOP)
     ENDIF

! normalize
     DO j = 1,ny
        ip = (j-1)*isize_line + i
        tmp1(ip,k) = aux(j,2)  *norm ! solution
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

  RETURN
END SUBROUTINE OPR_HELMHOLTZ_FXZ_2

!########################################################################
!########################################################################
! Same, but for n fields
SUBROUTINE OPR_HELMHOLTZ_FXZ_2_N(nx,ny,nz, nfield, ibc, alpha, &
     a, tmp1,tmp2, bcs_hb,bcs_ht, aux, wrk1d,wrk3d)

  USE TLAB_CONSTANTS, ONLY : efile
  USE TLAB_TYPES,     ONLY : pointers_dt
  USE TLAB_VARS,    ONLY : isize_txc_dimz
  USE TLAB_VARS,    ONLY : g
  USE TLAB_PROCS
#ifdef USE_MPI
  USE TLAB_MPI_VARS, ONLY : ims_offset_i, ims_offset_k
#endif

  IMPLICIT NONE
#include "integers.h"


  TINTEGER ibc, nx,ny,nz,nfield
  TREAL alpha
  TYPE(pointers_dt), DIMENSION(nfield)            :: a
  TREAL,    DIMENSION(nx,nz,nfield)               :: bcs_hb, bcs_ht
  TCOMPLEX, DIMENSION(isize_txc_dimz/2,nz,nfield) :: tmp1
  TCOMPLEX, DIMENSION(isize_txc_dimz/2,nz)        :: tmp2, wrk3d
  TCOMPLEX, DIMENSION(nfield,ny,2)                :: aux
  TREAL,    DIMENSION(ny,7)                       :: wrk1d
! -----------------------------------------------------------------------
  TINTEGER i,j,k, iglobal, kglobal, ip, isize_line,ifield
  TREAL lambda, norm
  TCOMPLEX bcs(3,nfield)

! #######################################################################
  norm = C_1_R/M_REAL(g(1)%size*g(3)%size)

  isize_line = nx/2+1

! #######################################################################
! Fourier transform of forcing term; output of this section in array tmp1
! #######################################################################
  DO ifield=1,nfield
     IF ( g(3)%size .GT. 1 ) THEN
        CALL OPR_FOURIER_F_X_EXEC(nx,ny,nz, a(ifield)%field,&
             bcs_hb(1,1,ifield),bcs_ht(1,1,ifield),tmp2, tmp1(1,1,ifield),wrk3d)
        CALL OPR_FOURIER_F_Z_EXEC(tmp2,tmp1(1,1,ifield)) ! tmp2 might be overwritten
     ELSE
        CALL OPR_FOURIER_F_X_EXEC(nx,ny,nz, a(ifield)%field,bcs_hb(1,1,ifield),bcs_ht(1,1,ifield),tmp1(1,1,ifield), tmp2,wrk3d)
     ENDIF
  ENDDO

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
     IF ( g(3)%size .GT. 1 ) THEN; lambda = g(1)%mwn(iglobal,2) + g(3)%mwn(kglobal,2)
     ELSE;                         lambda = g(1)%mwn(iglobal,2); ENDIF

     lambda = lambda - alpha

! forcing term
     DO ifield=1,nfield
        DO j = 1,ny
           ip = (j-1)*isize_line + i; aux(ifield,j,1) = tmp1(ip,k,ifield)
        ENDDO
! BCs
        j = ny+1; ip = (j-1)*isize_line + i;
        bcs(ifield,1) = tmp1(ip,k,ifield)   ! Dirichlet or Neumann
        j = ny+2; ip = (j-1)*isize_line + i;
        bcs(ifield,2) = tmp1(ip,k,ifield)   ! Dirichlet or Neumann
     ENDDO

! -----------------------------------------------------------------------
! Solve for each (kx,kz) a system of 1 complex equation as 2 independent real equations
! -----------------------------------------------------------------------
     IF ( ibc .EQ. 0 ) THEN ! Dirichlet BCs
        IF      ( g(2)%mode_fdm .EQ. FDM_COM6_JACOBIAN .OR. g(2)%mode_fdm .EQ. FDM_COM6_JACPENTA ) THEN
           CALL INT_C2N6_LHS_E(ny,    g(2)%jac, lambda, &
                wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5), wrk1d(1,6),wrk1d(1,7))
           CALL INT_C2N6_RHS  (ny,i2, g(2)%jac, aux(1,1,1),aux(1,1,2))
        ELSE IF ( g(2)%mode_fdm .EQ. FDM_COM6_DIRECT   ) THEN
           wrk1d = C_0_R
           CALL INT_C2N6N_LHS_E(ny,    g(2)%lu2(1,8), g(2)%lu2(1,4), lambda, &
                wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5), wrk1d(1,6),wrk1d(1,7))
           CALL INT_C2N6N_RHS  (ny,i2, g(2)%lu2(1,8), aux(1,1,1),aux(1,1,2))
        ENDIF

        CALL PENTADFS(ny-2,    wrk1d(2,1),wrk1d(2,2),wrk1d(2,3),wrk1d(2,4),wrk1d(2,5))

        CALL PENTADSS(ny-2,i1,        wrk1d(2,1),wrk1d(2,2),wrk1d(2,3),wrk1d(2,4),wrk1d(2,5),wrk1d(2,6))
        CALL PENTADSS(ny-2,i1,        wrk1d(2,1),wrk1d(2,2),wrk1d(2,3),wrk1d(2,4),wrk1d(2,5),wrk1d(2,7))
        CALL PENTADSS(ny-2,i2*nfield, wrk1d(2,1),wrk1d(2,2),wrk1d(2,3),wrk1d(2,4),wrk1d(2,5),aux(1,2,2))

        DO ifield=1,nfield
           aux(ifield,1:ny,2) = aux(ifield,1:ny,2) &
                + bcs(ifield,1)*wrk1d(1:ny,6) + bcs(ifield,2)*wrk1d(1:ny,7)
        ENDDO

     ELSE
        CALL TLAB_WRITE_ASCII(efile,'OPR_HELMHOLT_FXZ_2_N. Undeveloped BCs.')
        CALL TLAB_STOP(DNS_ERROR_UNDEVELOP)
     ENDIF

! normalize
     DO ifield=1,nfield
        DO j = 1,ny
           ip = (j-1)*isize_line + i
           tmp1(ip,k,ifield) = aux(ifield,j,2)  *norm ! solution
        ENDDO
     ENDDO
  ENDDO
  ENDDO

! ###################################################################
! Fourier field p (based on array tmp1)
! ###################################################################
  DO ifield=1,nfield
     IF ( g(3)%size .GT. 1 ) THEN
        CALL OPR_FOURIER_B_Z_EXEC(tmp1(1,1,ifield),wrk3d)
        CALL OPR_FOURIER_B_X_EXEC(nx,ny,nz, wrk3d,a(ifield)%field, tmp1(1,1,ifield))
     ELSE
        CALL OPR_FOURIER_B_X_EXEC(nx,ny,nz, tmp1(1,1,ifield),a(ifield)%field, wrk3d)
     ENDIF
  ENDDO

  RETURN
END SUBROUTINE OPR_HELMHOLTZ_FXZ_2_N
