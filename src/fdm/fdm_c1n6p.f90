!########################################################################
!# Tool/Library PADE
!#
!########################################################################
!# HISTORY
!#
!# 1999/01/01 - C. Pantano
!#              Created
!# 2007/10/01 - J.P. Mellado
!#              Fixed
!# 2009/12/09 - J.P. Mellado
!#              Splitting into two routines
!#
!########################################################################
!# DESCRIPTION
!#
!# Implementation of the first derivative finite difference with
!# 6th order tridiagonal compact scheme by JCP Lele 1992, nonperiodic.
!# Interior points according to Eq. 2.1.7 (\alpha=1/3).
!# System multiplied by 18/14 to eliminate one multiplication in the RHS
!#
!########################################################################
!# ARGUMENTS 
!#
!# u    In    function to be diferentiated
!# d    Out   right-hand side vector of the linear system
!#
!########################################################################
#include "types.h"

#define C_01D28_L 0.357142857142857d-1

! #######################################################################
! Left-hand side; tridiagonal matrix of the linear system
! #######################################################################
SUBROUTINE FDM_C1N6P_LHS(imax, dx, a,b,c)
  
  IMPLICIT NONE

  TINTEGER,                INTENT(IN) :: imax
  TREAL,   DIMENSION(imax),INTENT(IN) :: dx
  TREAL,   DIMENSION(imax),INTENT(OUT):: a,b,c

! -------------------------------------------------------------------
  TINTEGER i

! ###################################################################
  DO i = 1,imax
     a(i) = C_3_R /C_7_R
     b(i) = C_18_R/C_14_R
     c(i) = C_3_R /C_7_R
  ENDDO

! -------------------------------------------------------------------
! Jacobian Multiplication
! -------------------------------------------------------------------
  c(imax) = c(imax)*dx(1)
  b(1)    = b(1)   *dx(1)
  a(2)    = a(2)   *dx(1)

  DO i = 2,imax-1
     c(i-1) = c(i-1)*dx(i)
     b(i)   = b(i)  *dx(i)
     a(i+1) = a(i+1)*dx(i)
  ENDDO

  c(imax-1) = c(imax-1)*dx(imax)
  b(imax)   = b(imax)  *dx(imax)
  a(1)      = a(1)     *dx(imax)

  RETURN
END SUBROUTINE FDM_C1N6P_LHS

! #######################################################################
! Right-hand side; forcing term
! #######################################################################
SUBROUTINE FDM_C1N6P_RHS(imax,jkmax, u,d)
#ifdef USE_OPENMP
  USE OMP_LIB
#endif 
  
  IMPLICIT NONE

  TINTEGER,                      INTENT(IN) :: imax, jkmax
  TREAL,   DIMENSION(jkmax,imax),INTENT(IN) :: u
  TREAL,   DIMENSION(jkmax,imax),INTENT(OUT):: d

! -------------------------------------------------------------------
  TINTEGER i, jk, im2, im1, ip1, ip2

#ifdef USE_BLAS
  INTEGER ilen
  TREAL alpha
#else 
  TINTEGER                                  :: srt,end,siz, imm1
#endif

! #######################################################################
#ifdef USE_BLAS
  ilen = jkmax

  DO i = 1,imax
     im2 = i-2; im2=im2+imax-1; im2=MOD(im2,imax)+1
     im1 = i-1; im1=im1+imax-1; im1=MOD(im1,imax)+1
     ip1 = i+1; ip1=ip1+imax-1; ip1=MOD(ip1,imax)+1
     ip2 = i+2; ip2=ip2+imax-1; ip2=MOD(ip2,imax)+1

     !DVES is not part of BLAS but of ESSL - not supported on intel systems 
     !CALL DVES(ilen, u(1,ip1), 1, u(1,im1), 1, d(1,i), 1)    
     d(:,i) = -C_1_R*u(:,im1) 
     alpha=C_1_R
     CALL DAXPY(ilen, alpha, u(1,ip1), 1, d(1,i), 1)
     alpha = C_01D28_L
     CALL DAXPY(ilen, alpha, u(1,ip2), 1, d(1,i), 1)
     alpha =-C_01D28_L
     CALL DAXPY(ilen, alpha, u(1,im2), 1, d(1,i), 1)

  ENDDO

#else
!$omp parallel &
!$omp private(im1,ip1,im2,ip2,i,jk,srt,end,siz,imm1)  

  CALL DNS_OMP_PARTITION(imax,srt,end,siz) 
  imm1 = imax - 1 
  DO i = srt,end
     im2 = i-2; im2=im2+imm1; im2=MOD(im2,imax)+1
     im1 = i-1; im1=im1+imm1; im1=MOD(im1,imax)+1
     ip1 = i+1; ip1=ip1+imm1; ip1=MOD(ip1,imax)+1
     ip2 = i+2; ip2=ip2+imm1; ip2=MOD(ip2,imax)+1

     DO jk = 1,jkmax
        d(jk,i) = u(jk,ip1) - u(jk,im1) + C_01D28_L*(u(jk,ip2) - u(jk,im2))
     ENDDO

  ENDDO
!$omp end parallel 
#endif

  RETURN
END SUBROUTINE FDM_C1N6P_RHS


