!########################################################################
!# Tool/Library PADE
!#
!########################################################################
!# HISTORY
!#
!# 2021/12/23 - J. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Implementation of the interpolatory finite difference  first derivative
!# with 6th-order tridiagonal compact scheme by JCP Lele 1992, periodic.
!# Interior points according to Eq. B.1.1 (\alpha=9/62, \beta=0, c=0).
!# System multiplied by 62/63 to eliminate one multiplication in the RHS.
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
SUBROUTINE FDM_C1INT6P_LHS(imax, dx, a,b,c)
  
  IMPLICIT NONE

  TINTEGER,                  INTENT(IN ):: imax
  TREAL,    DIMENSION(imax), INTENT(IN ):: dx
  TREAL,    DIMENSION(imax), INTENT(OUT):: a,b,c

! -------------------------------------------------------------------
  TINTEGER i

! ###################################################################
  DO i = 1,imax
     a(i) = C_9_R  /  C_63_R ! 9/62
     b(i) = C_62_R /  C_63_R ! 1
     c(i) = C_9_R  /  C_63_R ! 9/62
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
END SUBROUTINE FDM_C1INT6P_LHS

! #######################################################################
! Right-hand side; forcing term 
! ==> interpolation from velocity to pressure grid
! #######################################################################
SUBROUTINE FDM_C1INTVP6P_RHS(imax,jkmax, u,d)
  
  IMPLICIT NONE

  TINTEGER,                        INTENT(IN ):: imax, jkmax
  TREAL,    DIMENSION(jkmax,imax), INTENT(IN ):: u
  TREAL,    DIMENSION(jkmax,imax), INTENT(OUT):: d

! -------------------------------------------------------------------
  TINTEGER                                    :: i, jk
  TINTEGER                                    :: im1, ip1, ip2, imm1
  TREAL                                       :: c17189

! #######################################################################

  c17189 = C_17_R / C_189_R

  imm1 = imax - 1 
  DO i = 1,imax
     im1 = i-1; im1=im1+imm1; im1=MOD(im1,imax)+1
     ip1 = i+1; ip1=ip1+imm1; ip1=MOD(ip1,imax)+1
     ip2 = i+2; ip2=ip2+imm1; ip2=MOD(ip2,imax)+1

     DO jk = 1,jkmax
        d(jk,i) = (u(jk,ip1) - u(jk,i)) + c17189*(u(jk,ip2) - u(jk,im1))
     ENDDO

  ENDDO

  RETURN
END SUBROUTINE FDM_C1INTVP6P_RHS

! #######################################################################
! Right-hand side; forcing term
! ==> interpolation from pressure to velocity grid
! #######################################################################
SUBROUTINE FDM_C1INTPV6P_RHS(imax,jkmax, u,d)
  
  IMPLICIT NONE

  TINTEGER,                        INTENT(IN ):: imax, jkmax
  TREAL,    DIMENSION(jkmax,imax), INTENT(IN ):: u
  TREAL,    DIMENSION(jkmax,imax), INTENT(OUT):: d

! -------------------------------------------------------------------
  TINTEGER                                    :: i, jk
  TINTEGER                                    :: im1, ip1, im2, imm1
  TREAL                                       :: c17189

! #######################################################################

  c17189 = C_17_R / C_189_R

  imm1 = imax - 1 
  DO i = 1,imax
     im1 = i-1; im1=im1+imm1; im1=MOD(im1,imax)+1
     im2 = i-2; im2=im2+imm1; im2=MOD(im2,imax)+1
     ip1 = i+1; ip1=ip1+imm1; ip1=MOD(ip1,imax)+1

     DO jk = 1,jkmax
        d(jk,i) = (u(jk,i) - u(jk,im1)) + c17189*(u(jk,ip1) - u(jk,im2))
     ENDDO

  ENDDO

  RETURN
END SUBROUTINE FDM_C1INTPV6P_RHS