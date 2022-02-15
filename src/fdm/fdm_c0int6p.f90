!########################################################################
!# Tool/Library PADE
!#
!########################################################################
!# HISTORY
!#
!# 2021/12/22 - J. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Implementation of the finite difference interpolation with
!# 6th order tridiagonal compact scheme by JCP Lele 1992, periodic.
!# Interior points according to Eq.C.1.4 (\alpha=3/10, \beta=0, c=0).
!# System multiplied by 4/3 to eliminate one multiplication in the RHS.
!#
!########################################################################
!# ARGUMENTS 
!#
!# u    In    function to be interpolated
!# d    Out   right-hand side vector of the linear system
!#
!########################################################################
#include "types.h"

! #######################################################################
! Left-hand side; tridiagonal matrix of the linear system
! #######################################################################
SUBROUTINE FDM_C0INT6P_LHS(imax, a,b,c)
  
  IMPLICIT NONE

  TINTEGER,                  INTENT(IN ):: imax
  TREAL,    DIMENSION(imax), INTENT(OUT):: a,b,c

! -------------------------------------------------------------------
  TINTEGER i

! ###################################################################
  DO i = 1,imax
     a(i) = C_2_R / C_5_R ! 3/10
     b(i) = C_4_R / C_3_R ! 1
     c(i) = C_2_R / C_5_R ! 3/10
  ENDDO

  RETURN
END SUBROUTINE FDM_C0INT6P_LHS

! #######################################################################
! Right-hand side; forcing term
! ==> interpolation from velocity to pressure grid
! #######################################################################
SUBROUTINE FDM_C0INTVP6P_RHS(imax,jkmax, u,d)
  
  IMPLICIT NONE

  TINTEGER,                        INTENT(IN ):: imax, jkmax
  TREAL,    DIMENSION(jkmax,imax), INTENT(IN ):: u
  TREAL,    DIMENSION(jkmax,imax), INTENT(OUT):: d

! -------------------------------------------------------------------
  TINTEGER                                    :: i, jk
  TINTEGER                                    :: im1, ip1, ip2, imm1
  TREAL                                       :: c0115

! #######################################################################

  c0115 = C_1_R / C_15_R

  imm1 = imax - 1 
  DO i = 1,imax
     im1 = i-1; im1=im1+imm1; im1=MOD(im1,imax)+1
     ip1 = i+1; ip1=ip1+imm1; ip1=MOD(ip1,imax)+1
     ip2 = i+2; ip2=ip2+imm1; ip2=MOD(ip2,imax)+1

     DO jk = 1,jkmax
        d(jk,i) = u(jk,ip1) + u(jk,i) + c0115*(u(jk,ip2) + u(jk,im1))
     ENDDO

  ENDDO

  RETURN
END SUBROUTINE FDM_C0INTVP6P_RHS

! #######################################################################
! Right-hand side; forcing term 
! ==> interpolation from pressure to velocity grid
! #######################################################################
SUBROUTINE FDM_C0INTPV6P_RHS(imax,jkmax, u,d)
  
  IMPLICIT NONE

  TINTEGER,                        INTENT(IN ):: imax, jkmax
  TREAL,    DIMENSION(jkmax,imax), INTENT(IN ):: u
  TREAL,    DIMENSION(jkmax,imax), INTENT(OUT):: d

! -------------------------------------------------------------------
  TINTEGER                                    :: i, jk
  TINTEGER                                    :: im1, ip1, im2, imm1
  TREAL                                       :: c0115

! #######################################################################

  c0115 = C_1_R / C_15_R

  imm1 = imax - 1 
  DO i = 1,imax
     im1 = i-1; im1=im1+imm1; im1=MOD(im1,imax)+1
     im2 = i-2; im2=im2+imm1; im2=MOD(im2,imax)+1
     ip1 = i+1; ip1=ip1+imm1; ip1=MOD(ip1,imax)+1

     DO jk = 1,jkmax
        d(jk,i) = u(jk,i) + u(jk,im1) + c0115*(u(jk,ip1) + u(jk,im2))
     ENDDO

  ENDDO

  RETURN
END SUBROUTINE FDM_C0INTPV6P_RHS