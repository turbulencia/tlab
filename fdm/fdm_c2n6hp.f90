#include "types.h"

!########################################################################
!# HISTORY
!#
!# 2020/08/22 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Implementation of the second derivative finite difference with
!# 6th order tridiagonal compact scheme.
!# Interior points according to JCP, Lamballais et al. 2011, JCP 230:3270-3275
!# Eqs. 1,3 with kc = pi**2.
!# It adds one term in the RHS to Lele's Eq. 2.2.7 scheme to better match the
!# exact transfer function (slightly hyper- instead of strongly hypodiffusive).
!#
!# The linear system is normalized to reduce number of operations in RHS
!#
!########################################################################
!# ARGUMENTS
!#
!# u    In    function to be diferentiated
!# d    Out   right-hand side vector of the linear system
!#
!########################################################################

! #######################################################################
! Left-hand side; tridiagonal matrix of the linear system
! #######################################################################
#define C_LHS0_L .147077436439844d+1
#define C_LHS1_L .536071798487849d+0

SUBROUTINE FDM_C2N6HP_LHS(imax, dx, a,b,c)

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

  a(1) = dxn*C_LHS1_L
  a(2) = dx1*C_LHS1_L
  b(1) = dx1*C_LHS0_L

  DO i = 2,imax-1
    dxi = dx(i)*dx(i)
    a(i+1) = dxi*C_LHS1_L
    b(i)   = dxi*C_LHS0_L
    c(i-1) = dxi*C_LHS1_L
  ENDDO

  b(imax)   = dxn*C_LHS0_L
  c(imax-1) = dxn*C_LHS1_L
  c(imax)   = dx1*C_LHS1_L

  RETURN
END SUBROUTINE FDM_C2N6HP_LHS

! #######################################################################
! Right-hand side; forcing term
! #######################################################################
#define C_RHS0_L .281250399533387d+1
#define C_RHS1_L .422670003525653d+0
#define C_RHS2_L .164180058587192d-1

SUBROUTINE FDM_C2N6HP_RHS(imax,jkmax, u,d)
#ifdef USE_OPENMP
  USE OMP_LIB
#endif

  IMPLICIT NONE

  TINTEGER,                    INTENT(IN) :: imax, jkmax
  TREAL, DIMENSION(jkmax,imax),INTENT(IN) :: u
  TREAL, DIMENSION(jkmax,imax),INTENT(OUT):: d

  ! -------------------------------------------------------------------
  TINTEGER i, jk, im3, im2, im1, ip1, ip2, ip3,imm1

#ifdef USE_BLAS
  INTEGER ilen
  TREAL alpha
#else
  TINTEGER srt,end,siz
#endif

  ! #######################################################################
#ifdef USE_BLAS
  ilen = jkmax
  imm1 = imax-1

  DO i = 1,imax
    im3 = i-3; im3=im3+imm1; im3=MOD(im3,imax)+1
    im2 = i-2; im2=im2+imm1; im2=MOD(im2,imax)+1
    im1 = i-1; im1=im1+imm1; im1=MOD(im1,imax)+1
    ip1 = i+1; ip1=ip1+imm1; ip1=MOD(ip1,imax)+1
    ip2 = i+2; ip2=ip2+imm1; ip2=MOD(ip2,imax)+1
    ip3 = i+3; ip3=ip3+imm1; ip3=MOD(ip3,imax)+1

    !CALL DVEA(ilen, u(1,ip1), 1, u(1,im1), 1, d(1,i), 1)
    !DVEA is not part of BLAS , but of ESSL -- not supported on intel systems
    CALL DCOPY(ilen,u(1,ip1),1,d(1,i),1)
    alpha=C_1_R
    CALL DAXPY(ilen,alpha,  u(1,im1), 1, d(1,i), 1)
    alpha =-C_RHS0_L
    CALL DAXPY(ilen, alpha, u(1,i  ), 1, d(1,i), 1)
    alpha = C_RHS1_L
    CALL DAXPY(ilen, alpha, u(1,ip2), 1, d(1,i), 1)
    CALL DAXPY(ilen, alpha, u(1,im2), 1, d(1,i), 1)
    alpha =-C_RHS2_L
    CALL DAXPY(ilen, alpha, u(1,ip3), 1, d(1,i), 1)
    CALL DAXPY(ilen, alpha, u(1,im3), 1, d(1,i), 1)

  ENDDO

#else
  !$omp parallel default ( none ) &
  !$omp private( jk,i,im3,im2,im1,ip1,ip2,ip3,srt,end,siz,imm1 ) &
  !$omp shared(d,u,imax,jkmax)

  CALL DNS_OMP_PARTITION(jkmax,srt,end,siz)
  imm1 = imax-1

  DO i = 1,imax
    im3 = i-3; im3=im3+imm1; im3=MOD(im3,imax)+1
    im2 = i-2; im2=im2+imm1; im2=MOD(im2,imax)+1
    im1 = i-1; im1=im1+imm1; im1=MOD(im1,imax)+1
    ip1 = i+1; ip1=ip1+imm1; ip1=MOD(ip1,imax)+1
    ip2 = i+2; ip2=ip2+imm1; ip2=MOD(ip2,imax)+1
    ip3 = i+3; ip3=ip3+imm1; ip3=MOD(ip3,imax)+1

    DO jk = srt,end
      d(jk,i) = u(jk,ip1) + u(jk,im1) - C_RHS0_L*u(jk,i)+ C_RHS1_L*(u(jk,ip2)+u(jk,im2)) - C_RHS2_L*(u(jk,ip3)+u(jk,im3))
    ENDDO

  ENDDO
  !$omp end parallel

#endif

  RETURN
END SUBROUTINE FDM_C2N6HP_RHS
