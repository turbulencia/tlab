!########################################################################
!# Tool/Library PADE
!#
!########################################################################
!# HISTORY
!#
!# 2009/12/13 - J. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Implementation of the first derivative finite difference with
!# 6th order pentadiagonal compact scheme by JCP Lele 1992, periodic.
!# Interior points according to Eq. 2.1.10. Similar truncation error like 
!# Eq. 2.1.7 with (\alpha=1/3). Here alpha value (\alpha=0.6047306974511406)
!# is chosen such, that no inflection point in w'(w) appears.
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
! Left-hand side; pentadiagonal matrix of the linear system
! #######################################################################
SUBROUTINE FDM_C1N6MP_LHS(imax, dx, a,b,c,d,e)
  
  IMPLICIT NONE

  TINTEGER,                INTENT(IN) :: imax
  TREAL,   DIMENSION(imax),INTENT(IN) :: dx
  TREAL,   DIMENSION(imax),INTENT(OUT):: a,b,c,d,e

! -------------------------------------------------------------------
  TINTEGER                            :: i
  TREAL                               :: alpha, beta

! ###################################################################
! LHS coefficients of the modified pentadiagonal scheme
  alpha = 0.604730585697398
  beta  = (C_2_R/C_5_R)*(alpha - C_1_R/C_3_R)

  DO i = 1,imax
    a(i) = beta 
    b(i) = alpha
    c(i) = 1
    d(i) = alpha
    e(i) = beta 
  ENDDO

! -------------------------------------------------------------------
! Jacobian Multiplication
! -------------------------------------------------------------------
  e(imax-1) = e(imax-1) * dx(1)
  d(imax  ) = d(imax  ) * dx(1)
  c(1     ) = c(1     ) * dx(1)
  b(2     ) = b(2     ) * dx(1)
  a(3     ) = a(3     ) * dx(1)
  !
  e(imax  ) = e(imax  ) * dx(2)
  d(1     ) = d(1     ) * dx(2)
  c(2     ) = c(2     ) * dx(2)
  b(3     ) = b(3     ) * dx(2)
  a(4     ) = a(4     ) * dx(2)
  !
  Do i = 3, imax-2
    e(i-2) = e(i-2) * dx(i)
    d(i-1) = d(i-1) * dx(i)
    c(i  ) = c(i  ) * dx(i)
    b(i+1) = b(i+1) * dx(i)
    a(i+2) = a(i+2) * dx(i)
  ENDDO
  !
  e(imax-3) = e(imax-3) * dx(imax-1)
  d(imax-2) = d(imax-2) * dx(imax-1)
  c(imax-1) = c(imax-1) * dx(imax-1)
  b(imax  ) = b(imax  ) * dx(imax-1)
  a(1     ) = a(1     ) * dx(imax-1)
  !
  e(imax-2) = e(imax-2) * dx(imax)
  d(imax-1) = d(imax-1) * dx(imax)
  c(imax  ) = c(imax  ) * dx(imax)
  b(1     ) = b(1     ) * dx(imax)
  a(2     ) = a(2     ) * dx(imax)

! -------------------------------------------------------------------
! Jacobian Multiplication - short version for only uniform grids
! -------------------------------------------------------------------
  ! DO i = 1,imax
  !    a(i) = beta  * dx(1) 
  !    b(i) = alpha * dx(1) 
  !    c(i) =         dx(1) 
  !    d(i) = alpha * dx(1) 
  !    e(i) = beta  * dx(1) 
  ! ENDDO

  RETURN
END SUBROUTINE FDM_C1N6MP_LHS

! #######################################################################
! Right-hand side; forcing term
! #######################################################################
SUBROUTINE FDM_C1N6P_RHS(imax,jkmax, u,d)
  
  IMPLICIT NONE

  TINTEGER,                      INTENT(IN) :: imax, jkmax
  TREAL,   DIMENSION(jkmax,imax),INTENT(IN) :: u
  TREAL,   DIMENSION(jkmax,imax),INTENT(OUT):: d

! -------------------------------------------------------------------
  TINTEGER                                  :: i, jk, im2, im1, ip1, ip2
  TINTEGER                                  :: srt,end,siz, imm1

! #######################################################################


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


  RETURN
END SUBROUTINE FDM_C1N6P_RHS


