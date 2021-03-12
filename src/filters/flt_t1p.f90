#include "types.h"

!########################################################################
!# DESCRIPTION
!#
!# Top-hat filter
!# Trapezoidal rule
!# Periodic boundary conditions
!# 
!########################################################################

!########################################################################
! Uniform grid
!########################################################################
SUBROUTINE FLT_T1P(kmax, ijmax, nx, z1, zf1)

  IMPLICIT NONE

  TINTEGER kmax, ijmax
  TINTEGER nx
  TREAL z1(ijmax,*)
  TREAL zf1(ijmax,*)

! -------------------------------------------------------------------
  TREAL dmul05 
  TREAL dmul1
  TINTEGER k,ij
  TINTEGER ii,im

  dmul1  = C_1_R/M_REAL(nx)
  dmul05 = C_05_R/M_REAL(nx)

! ##################################################
! # First nx/2 points where periodicity is imposed #
! ##################################################
  DO k = 1,nx/2
! first point and last point
     ii =  k-nx/2
     ii = ii+kmax-1
     ii = MOD(ii,kmax)+1
     im = k+nx/2
     DO ij = 1,ijmax
!        zf1(ij,k) = z1(ij,ii)+z1(ij,im)
        zf1(ij,k) = dmul05*(z1(ij,ii)+z1(ij,im))
     ENDDO
! inner points         
     DO ii = k-nx/2+1,k+nx/2-1
        im = ii+kmax-1
        im = MOD(im,kmax)+1
        DO ij = 1,ijmax
!           zf1(ij,k) = dmul05*(zf1(ij,k)+C_2_R*z1(ij,im))
           zf1(ij,k) = zf1(ij,k)+dmul1*z1(ij,im)
        ENDDO
     ENDDO
  ENDDO

! #################################################
! # Interior points where boundary does not enter #
! #################################################
  DO k = 1+nx/2,kmax-nx/2
! first point and last point
     ii = k-nx/2
     im = k+nx/2
     DO ij = 1,ijmax
!        zf1(ij,k) = z1(ij,ii)+z1(ij,im)
        zf1(ij,k) = dmul05*(z1(ij,ii)+z1(ij,im))
     ENDDO
! inner points
     DO ii = k-nx/2+1,k+nx/2-1
        DO ij = 1,ijmax
!           zf1(ij,k) = dmul05*(zf1(ij,k)+C_2_R*z1(ij,ii))
           zf1(ij,k) = zf1(ij,k)+dmul1*z1(ij,ii)
        ENDDO
     ENDDO
  ENDDO

! #################################################
! # Last nx/2 points where periodicity is imposed #
! #################################################
  DO k = kmax-nx/2+1,kmax
! first point and last point
     ii = k-nx/2
     im = k+nx/2
     im = im+kmax-1
     im = MOD(im,kmax)+1
     DO ij = 1,ijmax
!        zf1(ij,k) = z1(ij,ii)+z1(ij,im)
        zf1(ij,k) = dmul05*(z1(ij,ii)+z1(ij,im))
     ENDDO
! inner points
     DO ii = k-nx/2+1,k+nx/2-1
        im = ii+kmax-1
        im = MOD(im,kmax)+1
        DO ij = 1,ijmax
!           zf1(ij,k) = dmul05*(zf1(ij,k)+C_2_R*z1(ij,im))
           zf1(ij,k) = zf1(ij,k)+dmul1*z1(ij,im)
        ENDDO
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE FLT_T1P

!########################################################################
! Particularized for a stencil of 3 points, nx=2
!########################################################################
SUBROUTINE FLT_T1PD2(kmax, ijmax, z1, zf1)
  
  IMPLICIT NONE

  TINTEGER kmax, ijmax
  TREAL z1(ijmax,*)
  TREAL zf1(ijmax,*)

! -------------------------------------------------------------------
  TINTEGER k,ij

! ############################################
! # First point where periodicity is imposed #
! ############################################
  DO ij = 1,ijmax
     zf1(ij,1) = C_025_R*(z1(ij,kmax) + z1(ij,2)) + C_05_R*z1(ij,1)
  ENDDO

! #################################################
! # Interior points where boundary does not enter #
! #################################################
  DO k = 2,kmax-1
     DO ij = 1,ijmax
        zf1(ij,k) = C_025_R*(z1(ij,k-1) + z1(ij,k+1)) + C_05_R*z1(ij,k)
     ENDDO
  ENDDO

! ###########################################
! # Last point where periodicity is imposed #
! ###########################################
  DO ij = 1,ijmax
     zf1(ij,kmax) = C_025_R*(z1(ij,kmax-1) + z1(ij,1)) + C_05_R*z1(ij,kmax)
  ENDDO

  RETURN
END SUBROUTINE FLT_T1PD2

!########################################################################
! Particularized for a stencil of 5 points, nx=4
!########################################################################
SUBROUTINE FLT_T1PD4(kmax, ijmax, z1, zf1)
  
  IMPLICIT NONE

  TINTEGER kmax, ijmax
  TREAL z1(ijmax,*)
  TREAL zf1(ijmax,*)

! -------------------------------------------------------------------
  TINTEGER k,ij

! ###############################################
! # First 2 points where periodicity is imposed #
! ###############################################
  DO ij = 1,ijmax
     zf1(ij,1) = C_0125_R*(z1(ij,kmax-1) + z1(ij,3)) + C_025_R*(z1(ij,kmax) + z1(ij,1) + z1(ij,2))
  ENDDO

  DO ij = 1,ijmax
     zf1(ij,2) = C_0125_R*(z1(ij,kmax) + z1(ij,4)) + C_025_R*(z1(ij,1) + z1(ij,2) + z1(ij,3))
  ENDDO

! #################################################
! # Interior points where boundary does not enter #
! #################################################
  DO k = 3,kmax-2
     DO ij = 1,ijmax
        zf1(ij,k) = C_0125_R*(z1(ij,k-2) + z1(ij,k+2)) + C_025_R*(z1(ij,k-1)+z1(ij,k)+z1(ij,k+1))
     ENDDO
  ENDDO

! ##############################################
! # Last 2 points where periodicity is imposed #
! ##############################################
  DO ij = 1,ijmax
     zf1(ij,kmax-1) = C_0125_R*(z1(ij,kmax-3) + z1(ij,1)) + &
          C_025_R*(z1(ij,kmax-2) + z1(ij,kmax-1) + z1(ij,kmax))
  ENDDO

  DO ij = 1,ijmax
     zf1(ij,kmax) = C_0125_R*(z1(ij,kmax-2) + z1(ij,2)) + &
          C_025_R*(z1(ij,kmax-1) + z1(ij,kmax) + z1(ij,1))
  ENDDO

  RETURN
END SUBROUTINE FLT_T1PD4

!########################################################################
! Non-uniform grid
!########################################################################
SUBROUTINE FLT_T1P_ND(kmax, ijmax, nx, cf, z1, zf1)

  IMPLICIT NONE

  TINTEGER kmax, ijmax
  TINTEGER nx
  TREAL cf(nx+1,kmax)
  TREAL z1(ijmax,*)
  TREAL zf1(ijmax,*)

! -------------------------------------------------------------------
  TINTEGER k,ij,ic
  TINTEGER ii,im

! ##################################################
! # First nx/2 points where periodicity is imposed #
! ##################################################
  DO k = 1,nx/2
     ii = k-nx/2
     im = ii+kmax-1
     im = MOD(im,kmax)+1
     ic = ii-k+nx/2+1
     DO ij = 1,ijmax
        zf1(ij,k) = z1(ij,im)*cf(ic,k)
     ENDDO
     !     DO ii = k-nx/2,k+nx/2
     DO ii = k-nx/2+1,k+nx/2
        im = ii+kmax-1
        im = MOD(im,kmax)+1
        ic = ii-k+nx/2+1
        DO ij = 1,ijmax
           zf1(ij,k) = zf1(ij,k)+z1(ij,im)*cf(ic,k)
        ENDDO
     ENDDO
  ENDDO

! #################################################
! # Interior points where boundary does not enter #
! #################################################
  DO k = 1+nx/2,kmax-nx/2
     ii = k-nx/2
     ic = ii-k+nx/2+1
     DO ij = 1,ijmax
        zf1(ij,k) = z1(ij,ii)*cf(ic,k)
     ENDDO
     !     DO ii = k-nx/2,k+nx/2
     DO ii = k-nx/2+1,k+nx/2
        ic = ii-k+nx/2+1
        DO ij = 1,ijmax
           zf1(ij,k) = zf1(ij,k)+z1(ij,ii)*cf(ic,k)
        ENDDO
     ENDDO
  ENDDO

! #################################################
! # Last nx/2 points where periodicity is imposed #
! #################################################
  DO k = kmax-nx/2+1,kmax
     ii = k-nx/2
     im = ii+kmax-1
     im = MOD(im,kmax)+1
     ic = ii-k+nx/2+1
     DO ij = 1,ijmax
        zf1(ij,k) = z1(ij,im)*cf(ic,k)
     ENDDO
     !     DO ii = k-nx/2,k+nx/2
     DO ii = k-nx/2+1,k+nx/2
        im = ii+kmax-1
        im = MOD(im,kmax)+1
        ic = ii-k+nx/2+1
        DO ij = 1,ijmax
           zf1(ij,k) = zf1(ij,k)+z1(ij,im)*cf(ic,k)
        ENDDO
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE FLT_T1P_ND
