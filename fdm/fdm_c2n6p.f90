!########################################################################
!# Tool/Library PADE
!#
!########################################################################
!# HISTORY
!#
!# 1999/01/01 - C. Pantano
!#              Created
!# 2007/09/05 - J.P. Mellado
!#              Optimized
!# 2009/12/09 - J.P. Mellado
!#              Splitting into two routines
!#
!########################################################################
!# DESCRIPTION
!#
!# Implementation of the second derivative finite difference with
!# 6th order tridiagonal compact scheme by JCP, Lele 1992, periodic.
!# According to Eq. 2.2.7.
!# The linear system is normalized to reduce number of operations in RHS
!#
!########################################################################
!# ARGUMENTS 
!#
!# u    In    function to be diferentiated
!# d    Out   right-hand side vector of the linear system
!#
!########################################################################
#include "types.h"

#define C_01D06_L .166666666666667d+0
#define C_11D12_L .916666666666667d+0
#define C_17D08_L .212500000000000d+1
#define C_01D16_L .625000000000000d-1

! #######################################################################
! Left-hand side; tridiagonal matrix of the linear system
! #######################################################################
SUBROUTINE FDM_C2N6P_LHS(imax, dx, a,b,c)

  IMPLICIT NONE

  TINTEGER,                INTENT(IN) :: imax
  TREAL,   DIMENSION(imax),INTENT(IN) :: dx
  TREAL,   DIMENSION(imax),INTENT(OUT):: a,b,c

! -------------------------------------------------------------------
  TINTEGER i
  TREAL dx1, dxn, dxi

! #######################################################################
  dx1 = dx(1)   *dx(1)
  dxn = dx(imax)*dx(imax)

  a(1) = dxn*C_01D06_L 
  a(2) = dx1*C_01D06_L 
  b(1) = dx1*C_11D12_L 

  DO i = 2,imax-1
     dxi = dx(i)*dx(i)
     a(i+1) = dxi*C_01D06_L 
     b(i)   = dxi*C_11D12_L 
     c(i-1) = dxi*C_01D06_L 
  ENDDO

  b(imax)   = dxn*C_11D12_L 
  c(imax-1) = dxn*C_01D06_L 
  c(imax)   = dx1*C_01D06_L 

  RETURN
END SUBROUTINE FDM_C2N6P_LHS

! #######################################################################
! Right-hand side; forcing term
! #######################################################################
SUBROUTINE FDM_C2N6P_RHS(imax,jkmax, u,d)
#ifdef USE_OPENMP
  USE OMP_LIB
#endif

  IMPLICIT NONE

  TINTEGER,                    INTENT(IN) :: imax, jkmax
  TREAL, DIMENSION(jkmax,imax),INTENT(IN) :: u
  TREAL, DIMENSION(jkmax,imax),INTENT(OUT):: d

! -------------------------------------------------------------------
  TINTEGER i, jk, im2, im1, ip1, ip2

#ifdef USE_BLAS
  INTEGER ilen
  TREAL alpha
#else
  TINTEGER srt,end,siz,imm1
#endif

! #######################################################################
#ifdef USE_BLAS
  ilen = jkmax

  DO i = 1,imax
     im2 = i-2; im2=im2+imax-1; im2=MOD(im2,imax)+1
     im1 = i-1; im1=im1+imax-1; im1=MOD(im1,imax)+1
     ip1 = i+1; ip1=ip1+imax-1; ip1=MOD(ip1,imax)+1
     ip2 = i+2; ip2=ip2+imax-1; ip2=MOD(ip2,imax)+1

     !CALL DVEA(ilen, u(1,ip1), 1, u(1,im1), 1, d(1,i), 1)   
     !DVEA is not part of BLAS , but of ESSL -- not supported on intel systems 
     CALL DCOPY(ilen,u(1,ip1),1,d(1,i),1)  
     alpha=C_1_R
     CALL DAXPY(ilen,alpha,  u(1,im1), 1, d(1,i), x1) 
     alpha =-C_17D08_L
     CALL DAXPY(ilen, alpha, u(1,i  ), 1, d(1,i), 1)
     alpha = C_01D16_L
     CALL DAXPY(ilen, alpha, u(1,ip2), 1, d(1,i), 1)
     CALL DAXPY(ilen, alpha, u(1,im2), 1, d(1,i), 1)

  ENDDO

#else
!$omp parallel default ( none ) &
!$omp private( jk,i,im2,im1,ip1,ip2,srt,end,siz,imm1 ) &
!$omp shared(d,u,imax,jkmax) 
  
  CALL DNS_OMP_PARTITION(jkmax,srt,end,siz)
  imm1 = imax-1

  DO i = 1,imax
     im2 = i-2; im2=im2+imm1; im2=MOD(im2,imax)+1
     im1 = i-1; im1=im1+imm1; im1=MOD(im1,imax)+1
     ip1 = i+1; ip1=ip1+imm1; ip1=MOD(ip1,imax)+1
     ip2 = i+2; ip2=ip2+imm1; ip2=MOD(ip2,imax)+1

     DO jk = srt,end
        d(jk,i) = u(jk,ip1) + u(jk,im1) - C_17D08_L*u(jk,i)+ C_01D16_L*(u(jk,ip2)+u(jk,im2))
     ENDDO

  ENDDO
!$omp end parallel

#endif

  RETURN
END SUBROUTINE FDM_C2N6P_RHS

