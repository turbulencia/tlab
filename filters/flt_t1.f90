#include "types.h"

!########################################################################
!# DESCRIPTION
!#
!# Top-hat filter
!# Trapezoidal rule
!# Free boundary conditions
!# 
!########################################################################

!########################################################################
! Uniform grid
!########################################################################
SUBROUTINE FLT_T1(kmax, ijmax, nx, cf, z1, zf1)
        
  IMPLICIT NONE
  
  TINTEGER kmax, ijmax
  TINTEGER nx
  TREAL cf(2,nx/2)
  TREAL z1(ijmax,*)
  TREAL zf1(ijmax,*)

! -------------------------------------------------------------------
  TREAL r0, dmul05, dmul1
  TINTEGER i1, i2
  TINTEGER k,ij
  TINTEGER ii,im

  i1 = 1
  i2 = 2
  r0 = C_0_R

  dmul1  = C_1_R/M_REAL(nx)
  dmul05 = C_05_R/M_REAL(nx)

! ##################################################
! # First nx/2 points where periodicity is imposed #
! ##################################################
  DO k = 1,nx/2
! first 2 points reflect bc
     DO ij = 1,ijmax
        zf1(ij,k) = dmul1*(z1(ij,1)*cf(1,k)+z1(ij,2)*cf(2,k))
     ENDDO
! inner points
     DO ii = 3,k+nx/2-1
        DO ij = 1,ijmax
           zf1(ij,k) = zf1(ij,k)+dmul1*z1(ij,ii)
        ENDDO
     ENDDO
! last point
     IF ( nx .GT. 2 ) THEN
        ii =  k+nx/2
        DO ij = 1,ijmax
           zf1(ij,k) = zf1(ij,k)+dmul05*z1(ij,ii)
        ENDDO
     ENDIF
  ENDDO

! #################################################
! # Interior points where boundary does not enter #
! #################################################
  DO k = 1+nx/2,kmax-nx/2
! first point
     ii =  k-nx/2
     DO ij = 1,ijmax
        zf1(ij,k) = dmul05*z1(ij,ii)
     ENDDO
! inner points
     DO ii = k-nx/2+1,k+nx/2-1
        DO ij = 1,ijmax
           zf1(ij,k) = zf1(ij,k)+dmul1*z1(ij,ii)
        ENDDO
     ENDDO
! last point
     ii =  k+nx/2
     DO ij = 1,ijmax
        zf1(ij,k) = zf1(ij,k)+dmul05*z1(ij,ii)
     ENDDO
  ENDDO

! #################################################
! # Last nx/2 points where periodicity is imposed #
! #################################################
  DO k = kmax-nx/2+1,kmax
! last 2 points account for the bc
     im = kmax-k+i1
     DO ij = 1,ijmax
        zf1(ij,k) = dmul1*(z1(ij,kmax)*cf(1,im)+&
             z1(ij,kmax-1)*cf(2,im))
     ENDDO
! inner points
     DO ii = k-nx/2+1,kmax-2
        DO ij = 1,ijmax
           zf1(ij,k) = zf1(ij,k)+dmul1*z1(ij,ii)
        ENDDO
     ENDDO
! first point
     IF ( nx .GT. 2 ) THEN
        ii =  k-nx/2
        DO ij = 1,ijmax
           zf1(ij,k) = zf1(ij,k)+dmul05*z1(ij,ii)
        ENDDO
     ENDIF

  ENDDO

  RETURN
END SUBROUTINE FLT_T1

!########################################################################
! Non-uniform grid
!########################################################################
SUBROUTINE FLT_T1NU(kmax, ijmax, nx, cf, z1, zf1)

  IMPLICIT NONE
  
  TINTEGER kmax, ijmax
  TINTEGER nx
  TREAL cf(nx+1,kmax)
  TREAL z1(ijmax,*)
  TREAL zf1(ijmax,*)

  TREAL r0, dum
  TINTEGER i1, i2
  TINTEGER k,ij,ic
  TINTEGER ii, iiorg

  i1 = 1
  i2 = 2
  r0 = C_0_R

! ############################################
! # First nx/2 points where bc are imposed #
! ############################################
  DO k = 1,nx/2
     iiorg = nx/2-k+1
     ii = 1
     ic = iiorg + ii
     dum = cf(ic,k)
     DO ij = 1,ijmax
        zf1(ij,k) = z1(ij,ii)*dum
     ENDDO
     DO ii = 2,k+nx/2
        ic = iiorg + ii
        dum = cf(ic,k)
        DO ij = 1,ijmax
           zf1(ij,k) = zf1(ij,k)+z1(ij,ii)*dum
        ENDDO
     ENDDO
  ENDDO

! #################################################
! # Interior points where boundary does not enter #
! #################################################
  DO k = 1+nx/2,kmax-nx/2
     iiorg = nx/2-k+1
     ii = k-nx/2
     ic = ii + iiorg
     dum = cf(ic,k)
     DO ij = 1,ijmax
        zf1(ij,k) = z1(ij,ii)*dum
     ENDDO
     DO ii = k-nx/2+1,k+nx/2
        ic = ii + iiorg
        dum = cf(ic,k)
        DO ij = 1,ijmax
           zf1(ij,k) = zf1(ij,k)+z1(ij,ii)*dum
        ENDDO
     ENDDO
  ENDDO

! ###########################################
! # Last nx/2 points where bc are imposed #
! ###########################################
  DO k = kmax-nx/2+1,kmax
     iiorg = nx/2-k+1
     ii = k-nx/2
     ic = ii + iiorg
     dum = cf(ic,k)
     DO ij = 1,ijmax
        zf1(ij,k) = z1(ij,ii)*dum
     ENDDO
     DO ii = k-nx/2+1,kmax
        ic = ii + iiorg
        dum = cf(ic,k)
        DO ij = 1,ijmax
           zf1(ij,k) = zf1(ij,k)+z1(ij,ii)*dum
        ENDDO
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE FLT_T1NU

!########################################################################
! Non-uniform grid
! Particularized for a stencil of 3 points, nx=2
!########################################################################
SUBROUTINE FLT_T1NUN2(kmax, ijmax, cf, z1, zf1)
  
  IMPLICIT NONE

  TINTEGER kmax, ijmax
  TREAL cf(3,kmax)
  TREAL z1(ijmax,*)
  TREAL zf1(ijmax,*)

  TINTEGER ij, k
  TREAL a1, a2, a3

! ######################################
! # First point where bc are imposed #
! ######################################
  a2 = cf(2,1)
  a3 = cf(3,1)
  DO ij = 1,ijmax
     zf1(ij,1) = z1(ij,1)*a2 + z1(ij,2)*a3
  ENDDO

! #################################################
! # Interior points where boundary does not enter #
! #################################################
  DO k = 2,kmax-1
     a1 = cf(1,k)
     a2 = cf(2,k)
     a3 = cf(3,k)
     DO ij = 1,ijmax
        zf1(ij,k) = z1(ij,k-1)*a1 + z1(ij,k)*a2 + z1(ij,k+1)*a3
     ENDDO
  ENDDO

! #####################################
! # Last point where bc are imposed #
! #####################################
  a1 = cf(1,kmax)
  a2 = cf(2,kmax)
  DO ij = 1,ijmax
     zf1(ij,kmax) = z1(ij,kmax-1)*a1 + z1(ij,kmax)*a2
  ENDDO

  RETURN
END SUBROUTINE FLT_T1NUN2

!########################################################################
! Non-uniform grid
! Particularized for a stencil of 5 points, nx=4
!########################################################################
SUBROUTINE FLT_T1NUN4(kmax, ijmax, cf, z1, zf1)
  
  IMPLICIT NONE

  TINTEGER kmax, ijmax
  TREAL cf(5,kmax)
  TREAL z1(ijmax,*)
  TREAL zf1(ijmax,*)

  TINTEGER ij, k
  TREAL a1, a2, a3, a4, a5

! #########################################
! # First 2 points where bc are imposed #
! #########################################
  a3 = cf(3,1)
  a4 = cf(4,1)
  a5 = cf(5,1)
  DO ij = 1,ijmax
     zf1(ij,1) = z1(ij,1)*a3 + z1(ij,2)*a4 + z1(ij,3)*a5
  ENDDO

  a2 = cf(2,2)
  a3 = cf(3,2)
  a4 = cf(4,2)
  a5 = cf(5,2)
  DO ij = 1,ijmax
     zf1(ij,2) = z1(ij,1)*a2 + z1(ij,2)*a3 + z1(ij,3)*a4 + z1(ij,4)*a5
  ENDDO

! #################################################
! # Interior points where boundary does not enter #
! #################################################
  DO k = 3,kmax-2
     a1 = cf(1,k)
     a2 = cf(2,k)
     a3 = cf(3,k)
     a4 = cf(4,k)
     a5 = cf(5,k)
     DO ij = 1,ijmax
        zf1(ij,k) = z1(ij,k-2)*a1 + z1(ij,k-1)*a2 + z1(ij,k)*a3 + z1(ij,k+1)*a4 + z1(ij,k+2)*a5
     ENDDO
  ENDDO

! ########################################
! # Last 2 points where bc are imposed #
! ########################################
  a1 = cf(1,kmax-1)
  a2 = cf(2,kmax-1)
  a3 = cf(3,kmax-1)
  a4 = cf(4,kmax-1)
  DO ij = 1,ijmax
     zf1(ij,kmax-1) = z1(ij,kmax-3)*a1 + z1(ij,kmax-2)*a2 + z1(ij,kmax-1)*a3 + z1(ij,kmax)*a4
  ENDDO

  a1 = cf(1,kmax)
  a2 = cf(2,kmax)
  a3 = cf(3,kmax)
  DO ij = 1,ijmax
     zf1(ij,kmax) = z1(ij,kmax-2)*a1 + z1(ij,kmax-1)*a2 + z1(ij,kmax)*a3
  ENDDO

  RETURN
END SUBROUTINE FLT_T1NUN4

!########################################################################
! Non-uniform grid
! Particularized for a stencil of 7 points, nx=6
!########################################################################
SUBROUTINE FLT_T1NUN6(kmax, ijmax, cf, z1, zf1)
  
  IMPLICIT NONE

  TINTEGER kmax, ijmax
  TREAL cf(7,kmax)
  TREAL z1(ijmax,*)
  TREAL zf1(ijmax,*)

  TINTEGER ij, k
  TREAL a1, a2, a3, a4, a5, a6, a7

! #########################################
! # First 3 points where bc are imposed #
! #########################################
  a4 = cf(4,1)
  a5 = cf(5,1)
  a6 = cf(6,1)
  a7 = cf(7,1)
  DO ij = 1,ijmax
     zf1(ij,1) = z1(ij,1)*a4 + z1(ij,2)*a5 + z1(ij,3)*a6 + z1(ij,4)*a7
  ENDDO

  a3 = cf(3,2)
  a4 = cf(4,2)
  a5 = cf(5,2)
  a6 = cf(6,2)
  a7 = cf(7,2)
  DO ij = 1,ijmax
     zf1(ij,2) = z1(ij,1)*a3 + z1(ij,2)*a4 + z1(ij,3)*a5 + z1(ij,4)*a6 + z1(ij,5)*a7
  ENDDO

  a2 = cf(2,3)
  a3 = cf(3,3)
  a4 = cf(4,3)
  a5 = cf(5,3)
  a6 = cf(6,3)
  a7 = cf(7,3)
  DO ij = 1,ijmax
     zf1(ij,3) = z1(ij,1)*a2 + z1(ij,2)*a3 + z1(ij,3)*a4 + z1(ij,4)*a5 + z1(ij,5)*a6 + z1(ij,6)*a7
  ENDDO

! #################################################
! # Interior points where boundary does not enter #
! #################################################
  DO k = 4,kmax-3
     a1 = cf(1,k)
     a2 = cf(2,k)
     a3 = cf(3,k)
     a4 = cf(4,k)
     a5 = cf(5,k)
     a6 = cf(6,k)
     a7 = cf(7,k)
     DO ij = 1,ijmax
        zf1(ij,k) = z1(ij,k-3)*a1 + z1(ij,k-2)*a2 + z1(ij,k-1)*a3 + &
             z1(ij,k)*a4 + z1(ij,k+1)*a5 + z1(ij,k+2)*a6 + z1(ij,k+3)*a7
     ENDDO
  ENDDO

! ########################################
! # Last 3 points where bc are imposed #
! ########################################
  a1 = cf(1,kmax-2)
  a2 = cf(2,kmax-2)
  a3 = cf(3,kmax-2)
  a4 = cf(4,kmax-2)
  a5 = cf(5,kmax-2)
  a6 = cf(6,kmax-2)
  DO ij = 1,ijmax
     zf1(ij,kmax-2) = z1(ij,kmax-5)*a1 + z1(ij,kmax-4)*a2 +&
          z1(ij,kmax-3)*a3 + z1(ij,kmax-2)*a4 + z1(ij,kmax-1)*a5 + z1(ij,kmax)*a6
  ENDDO

  a1 = cf(1,kmax-1)
  a2 = cf(2,kmax-1)
  a3 = cf(3,kmax-1)
  a4 = cf(4,kmax-1)
  a5 = cf(5,kmax-1)
  DO ij = 1,ijmax
     zf1(ij,kmax-1) = z1(ij,kmax-4)*a1 + z1(ij,kmax-3)*a2 +&
          z1(ij,kmax-2)*a3 + z1(ij,kmax-1)*a4 + z1(ij,kmax)*a5
  ENDDO

  a1 = cf(1,kmax)
  a2 = cf(2,kmax)
  a3 = cf(3,kmax)
  a4 = cf(4,kmax)
  DO ij = 1,ijmax
     zf1(ij,kmax) = z1(ij,kmax-3)*a1 + z1(ij,kmax-2)*a2 + z1(ij,kmax-1)*a3 + z1(ij,kmax)*a4
  ENDDO

  RETURN
END SUBROUTINE FLT_T1NUN6
