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
!# 2010/03/09 - J.P. Mellado
!#              Splitting into two routines
!#
!########################################################################
!# DESCRIPTION
!#
!# Implementation of the first derivative finite difference with
!# 8th order tridiagonal compact scheme by JCP Lele 1992, nonperiodic.
!# Interior points according to Eq. 2.1.13 with \alpha=3/8.
!#
!########################################################################
!# ARGUMENTS 
!#
!# u    In    function to be diferentiated
!# d    Out   right-hand side vector of the linear system
!#
!########################################################################
#include "types.h"

! #######################################################################
! Left-hand side; tridiagonal matrix of the linear system
! #######################################################################
SUBROUTINE FDM_C1N8P_LHS(imax, dx, a,b,c)
  
  IMPLICIT NONE

  TINTEGER,                INTENT(IN) :: imax
  TREAL,   DIMENSION(imax),INTENT(IN) :: dx
  TREAL,   DIMENSION(imax),INTENT(OUT):: a,b,c

! -------------------------------------------------------------------
  TINTEGER i

! ###################################################################
  DO i = 1,imax
     a(i) = C_3_R/C_8_R
     b(i) = C_1_R
     c(i) = C_3_R/C_8_R
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
END SUBROUTINE FDM_C1N8P_LHS

! #######################################################################
! Right-hand side; forcing term
! #######################################################################
SUBROUTINE FDM_C1N8P_RHS(imax, jkmax, u,d)
#ifdef USE_OPENMP
  USE OMP_LIB
#endif

  IMPLICIT NONE

  TINTEGER,                      INTENT(IN) :: imax, jkmax
  TREAL,   DIMENSION(jkmax,imax),INTENT(IN) :: u
  TREAL,   DIMENSION(jkmax,imax),INTENT(OUT):: d

! -------------------------------------------------------------------
  TINTEGER i,jk, im3,im2,im1,ip1,ip2,ip3,imm1
  TINTEGER srt,siz,end
  TREAL c32dx, c20dx, c48dx


!$omp parallel default( none ) &
!$omp private( jk,i,im3,im2,im1,ip1,ip2,ip3,srt,end,siz,c32dx,c48dx,c20dx,imm1) &
!$omp shared(d,u,jkmax,imax)

  CALL DNS_OMP_PARTITION(imax,srt,end,siz) 

! ###################################################################
  c32dx = C_5_R*C_5_R/(C_4_R*C_8_R)
  c20dx = C_1_R      /(C_2_R*C_10_R)
  c48dx =-C_1_R      /(C_3_R*C_10_R*C_16_R)

  imm1 = imax-1

  DO i = srt,end
     im3 = i-3; im3=im3+imm1; im3=MOD(im3,imax)+1
     im2 = i-2; im2=im2+imm1; im2=MOD(im2,imax)+1
     im1 = i-1; im1=im1+imm1; im1=MOD(im1,imax)+1
     ip1 = i+1; ip1=ip1+imm1; ip1=MOD(ip1,imax)+1
     ip2 = i+2; ip2=ip2+imm1; ip2=MOD(ip2,imax)+1
     ip3 = i+3; ip3=ip3+imm1; ip3=MOD(ip3,imax)+1

     DO jk = 1,jkmax
        d(jk,i) = c32dx*(u(jk,ip1) - u(jk,im1)) + c20dx*(u(jk,ip2) - u(jk,im2)) &
                + c48dx*(u(jk,ip3) - u(jk,im3))
     ENDDO
  ENDDO

!$omp end parallel 

  RETURN
END SUBROUTINE FDM_C1N8P_RHS
