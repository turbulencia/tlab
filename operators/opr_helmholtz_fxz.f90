#include "types.h"
#include "dns_const.h"
#include "dns_error.h"
#include "dns_const_mpi.h"

!########################################################################
!# Tool/Library
!#
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
SUBROUTINE OPR_HELMHOLTZ_FXZ(nx,ny,nz, ibc, alpha,&
     dx,dy,dz, a, tmp1,tmp2, bcs_hb,bcs_ht, aux, wrk1d,wrk3d)

  USE DNS_GLOBAL, ONLY : isize_txc_dimz, imax_total,jmax_total,kmax_total, inb_grid_1
  USE DNS_GLOBAL, ONLY : imode_fdm
  USE DNS_CONSTANTS, ONLY : efile
#ifdef USE_MPI
  USE DNS_MPI, ONLY : ims_offset_i, ims_offset_k
#endif

  IMPLICIT NONE

#include "integers.h"

  TINTEGER nx,ny,nz, ibc
  TREAL alpha
  TREAL,    DIMENSION(imax_total,*)        :: dx
  TREAL,    DIMENSION(jmax_total,*)        :: dy
  TREAL,    DIMENSION(kmax_total,*)        :: dz
  TREAL,    DIMENSION(nx,ny,nz)            :: a
  TREAL,    DIMENSION(nx,nz)               :: bcs_hb, bcs_ht  
  TCOMPLEX, DIMENSION(isize_txc_dimz/2,nz) :: tmp1,tmp2, wrk3d
  TCOMPLEX, DIMENSION(ny,2)                :: aux
  TCOMPLEX, DIMENSION(ny,6)                :: wrk1d

! -----------------------------------------------------------------------
  TINTEGER i,j,k, iglobal,kglobal, ip, isize_line
  TREAL lambda, norm
  TCOMPLEX bcs(3)

! #######################################################################
  norm = C_1_R/M_REAL(imax_total*kmax_total)

  isize_line = nx/2+1

! #######################################################################
! Fourier transform of forcing term; output of this section in array tmp1
! #######################################################################
  IF ( kmax_total .GT. 1 ) THEN
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
     ip = inb_grid_1 + 5 ! pointer to position in arrays dx, dz
     IF ( kmax_total .GT. 1 ) THEN; lambda = dx(iglobal,ip) + dz(kglobal,ip)
     ELSE;                          lambda = dx(iglobal,ip); ENDIF

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
        CALL FDE_BVP_REGULAR_NN(imode_fdm, ny,i2, lambda, dy, aux(1,2),aux(1,1), bcs, wrk1d(1,1), wrk1d(1,2))
     ELSE IF ( ibc .EQ. 0 ) THEN ! Dirichlet BCs
        CALL FDE_BVP_REGULAR_DD(imode_fdm, ny,i2, lambda, dy, aux(1,2),aux(1,1), bcs, wrk1d(1,1), wrk1d(1,2))
     ELSE
        CALL IO_WRITE_ASCII(efile,'OPR_HELMHOLT_FXZ. Undeveloped BCs.')
        CALL DNS_STOP(DNS_ERROR_UNDEVELOP)
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
  IF ( kmax_total .GT. 1 ) THEN
     CALL OPR_FOURIER_B_Z_EXEC(tmp1,wrk3d) 
     CALL OPR_FOURIER_B_X_EXEC(nx,ny,nz, wrk3d,a, tmp1) 
  ELSE
     CALL OPR_FOURIER_B_X_EXEC(nx,ny,nz, tmp1,a, wrk3d) 
  ENDIF

  RETURN
END SUBROUTINE OPR_HELMHOLTZ_FXZ

!########################################################################
!########################################################################
SUBROUTINE OPR_HELMHOLTZ_FXZ_2(nx,ny,nz, ibc, alpha,&
     dx,dy,dz, a, tmp1,tmp2, bcs_hb,bcs_ht, aux, wrk1d,wrk3d)

  USE DNS_GLOBAL, ONLY : isize_field, isize_txc_dimz, imax_total,jmax_total,kmax_total, inb_grid_2
  USE DNS_GLOBAL, ONLY : imode_fdm
  USE DNS_CONSTANTS, ONLY : efile
#ifdef USE_MPI
  USE DNS_MPI, ONLY : ims_offset_i, ims_offset_k
#endif

  IMPLICIT NONE

#include "integers.h"

  TINTEGER nx,ny,nz, ibc
  TREAL alpha
  TREAL,    DIMENSION(imax_total,*)        :: dx
  TREAL,    DIMENSION(jmax_total,*)        :: dy
  TREAL,    DIMENSION(kmax_total,*)        :: dz
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
  norm = C_1_R/M_REAL(imax_total*kmax_total)

  isize_line = nx/2+1

! #######################################################################
! Fourier transform of forcing term; output of this section in array tmp1
! #######################################################################
  IF ( kmax_total .GT. 1 ) THEN
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
     ip = inb_grid_2 + 5 ! pointer to position in arrays dx, dz
     IF ( kmax_total .GT. 1 ) THEN; lambda = dx(iglobal,ip) + dz(kglobal,ip)
     ELSE;                          lambda = dx(iglobal,ip); ENDIF

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
        IF      ( imode_fdm .EQ. FDM_COM6_JACOBIAN ) THEN
           CALL INT_C2N6_LHS_E(ny,    dy, lambda, &
                wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5), wrk1d(1,6),wrk1d(1,7))
           CALL INT_C2N6_RHS  (ny,i2, dy, aux(1,1),aux(1,2))
        ELSE IF ( imode_fdm .EQ. FDM_COM6_DIRECT   ) THEN
           wrk1d = C_0_R
           CALL INT_C2N6N_LHS_E(ny,    dy(1,inb_grid_2+7),dy(1,inb_grid_2+3), lambda, &
                wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5), wrk1d(1,6),wrk1d(1,7))
           CALL INT_C2N6N_RHS  (ny,i2, dy(1,inb_grid_2+7), aux(1,1),aux(1,2))
        ENDIF

        CALL PENTADFS(ny-2,    wrk1d(2,1),wrk1d(2,2),wrk1d(2,3),wrk1d(2,4),wrk1d(2,5))

        CALL PENTADSS(ny-2,i1, wrk1d(2,1),wrk1d(2,2),wrk1d(2,3),wrk1d(2,4),wrk1d(2,5), wrk1d(2,6))
        CALL PENTADSS(ny-2,i1, wrk1d(2,1),wrk1d(2,2),wrk1d(2,3),wrk1d(2,4),wrk1d(2,5), wrk1d(2,7))

        CALL PENTADSS(ny-2,i2, wrk1d(2,1),wrk1d(2,2),wrk1d(2,3),wrk1d(2,4),wrk1d(2,5), aux(2,2))

        aux(:,2) = aux(:,2) +  bcs(1)*wrk1d(:,6) + bcs(2)*wrk1d(:,7)

     ELSE
        CALL IO_WRITE_ASCII(efile,'OPR_HELMHOLT_FXZ_2. Undeveloped BCs.')
        CALL DNS_STOP(DNS_ERROR_UNDEVELOP)
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
  IF ( kmax_total .GT. 1 ) THEN
     CALL OPR_FOURIER_B_Z_EXEC(tmp1,wrk3d) 
     CALL OPR_FOURIER_B_X_EXEC(nx,ny,nz, wrk3d,a, tmp1) 
  ELSE
     CALL OPR_FOURIER_B_X_EXEC(nx,ny,nz, tmp1,a, wrk3d) 
  ENDIF

  RETURN
END SUBROUTINE OPR_HELMHOLTZ_FXZ_2
