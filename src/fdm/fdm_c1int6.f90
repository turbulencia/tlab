!########################################################################
!# Tool/Library PADE
!#
!########################################################################
!# HISTORY
!#
!# 2021/01/14 - J. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Implementation of the interpolatory finite difference  first derivative
!# with 6th-order tridiagonal compact scheme by JCP Lele 1992, non-periodic.
!# Interior points according to Eq. B.1.1 (\alpha=9/62, \beta=0, c=0).
!# System for this scheme is multiplied by 4/3 to eliminate one 
!# multiplication in the RHS.
!#
!# Different boundary closures can be found in Albin 2010
!# (https://doi.org/10.1002/fld.2520) table 6 (typo in Sd2_3i).
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
! ==> interpolation to the right
! #######################################################################
SUBROUTINE FDM_C1INT6R_LHS(imax, imaxp, dx, a,b,c)
  
  IMPLICIT NONE

  TINTEGER,                   INTENT(IN ):: imax  ! velocity grid 
  TINTEGER,                   INTENT(IN ):: imaxp ! pressure grid (imaxp==imax-1)
  TREAL,    DIMENSION(imax ), INTENT(IN ):: dx
  TREAL,    DIMENSION(imaxp), INTENT(OUT):: a,b,c

! -------------------------------------------------------------------
  TINTEGER                               :: i

! ###################################################################
! ! third-order biased (explicit)
!   a(1)     = C_0_R
!   b(1)     = C_1_R
!   c(1)     = C_0_R
!   !
!   a(imaxp) = C_0_R
!   b(imaxp) = C_1_R
!   c(imaxp) = C_0_R

! third-order biased (implicit)
  a(1)     =  C_0_R
  b(1)     =  C_1_R
  c(1)     = -C_1_R
  !
  a(imaxp) = -C_1_R
  b(imaxp) =  C_1_R
  c(imaxp) =  C_0_R

! sixth-order (implicit)
  DO i = 2,imaxp-1
     a(i) = C_9_R  /  C_63_R ! 9/62
     b(i) = C_62_R /  C_63_R ! 1
     c(i) = C_9_R  /  C_63_R ! 9/62
  ENDDO

! ! forth-order (implicit) - for adjacent boundary nodes 
!   a(2)       = C_1_R / C_22_R
!   b(2)       = C_1_R 
!   c(2)       = C_1_R / C_22_R
!   !
!   a(imaxp-1) = C_1_R / C_22_R
!   b(imaxp-1) = C_1_R
!   c(imaxp-1) = C_1_R / C_22_R

! -------------------------------------------------------------------
! Jacobian Multiplication
! -------------------------------------------------------------------
! just for uniform grids !!!!!!
  DO i = 1,imaxp
    a(i) = a(i)*dx(1)
    b(i) = b(i)*dx(1)
    c(i) = c(i)*dx(1)
  ENDDO

  ! c(imaxp) = c(imaxp)*dx(1)
  ! b(1)     = b(1)    *dx(1)
  ! a(2)     = a(2)    *dx(1)

  ! DO i = 2,imaxp-1
  !   c(i-1) = c(i-1)*dx(i)
  !   b(i)   = b(i)  *dx(i)
  !   a(i+1) = a(i+1)*dx(i)
  ! ENDDO

  ! c(imaxp-1) = c(imaxp-1)*dx(imaxp)
  ! b(imaxp)   = b(imaxp)  *dx(imaxp)
  ! a(1)       = a(1)      *dx(imaxp)
  
  RETURN
END SUBROUTINE FDM_C1INT6R_LHS

! #######################################################################
! Left-hand side; tridiagonal matrix of the linear system  
! ==> interpolation to the left
! #######################################################################
SUBROUTINE FDM_C1INT6L_LHS(imax, imaxp, dx, a,b,c)
  
  IMPLICIT NONE

  TINTEGER,                  INTENT(IN ):: imax  ! velocity grid 
  TINTEGER,                  INTENT(IN ):: imaxp ! pressure grid (imaxp==imax-1)
  TREAL,    DIMENSION(imax), INTENT(IN ):: dx
  TREAL,    DIMENSION(imax), INTENT(OUT):: a,b,c

! -------------------------------------------------------------------
  TINTEGER                              :: i

! ###################################################################
! ! third-order biased (explicit)
!   a(1)    = C_0_R
!   b(1)    = C_1_R
!   c(1)    = C_0_R
!   !
!   a(imax) = C_0_R
!   b(imax) = C_1_R
!   c(imax) = C_0_R

! third-order biased (implicit)
  a(1)      = C_0_R
  b(1)      = C_1_R
  c(1)      = C_23_R
  !
  a(imax)   = C_23_R
  b(imax)   = C_1_R
  c(imax)   = C_0_R

  ! forth-order (implicit) - for adjacent boundary nodes 
  a(2)      = C_1_R / C_22_R
  b(2)      = C_1_R 
  c(2)      = C_1_R / C_22_R
  !
  a(imax-1) = C_1_R / C_22_R
  b(imax-1) = C_1_R
  c(imax-1) = C_1_R / C_22_R
  
! sixth-order (implicit)
  DO i = 3,imax-2
     a(i) = C_9_R  /  C_63_R ! 9/62
     b(i) = C_62_R /  C_63_R ! 1
     c(i) = C_9_R  /  C_63_R ! 9/62
  ENDDO

! -------------------------------------------------------------------
! Jacobian Multiplication
! -------------------------------------------------------------------
! just for uniform grids !!!!!!
  DO i = 1,imax
    a(i) = a(i)*dx(1)
    b(i) = b(i)*dx(1)
    c(i) = c(i)*dx(1)
  ENDDO

  ! c(imax) = c(imax)*dx(1)
  ! b(1)    = b(1)   *dx(1)
  ! a(2)    = a(2)   *dx(1)

  ! DO i = 2,imax-1
  !   c(i-1) = c(i-1)*dx(i)
  !   b(i)   = b(i)  *dx(i)
  !   a(i+1) = a(i+1)*dx(i)
  ! ENDDO

  ! c(imax-1) = c(imax-1)*dx(imax)
  ! b(imax)   = b(imax)  *dx(imax)
  ! a(1)      = a(1)     *dx(imax)

  RETURN
END SUBROUTINE FDM_C1INT6L_LHS

! #######################################################################
! Right-hand side; forcing term ==> interpolation to the right
! #######################################################################
SUBROUTINE FDM_C1INT6R_RHS(imax,imaxp,jkmax, u,d)
  
  IMPLICIT NONE

  TINTEGER,                         INTENT(IN ):: imax, imaxp, jkmax
  TREAL,    DIMENSION(jkmax,imax),  INTENT(IN ):: u
  TREAL,    DIMENSION(jkmax,imaxp), INTENT(OUT):: d

! -------------------------------------------------------------------
  TINTEGER                                     :: i, jk
  TREAL                                        :: c17189, c0120
  ! TREAL                                        :: c0108, c0708, c0124, c2324, c1211

! #######################################################################
  c17189 = C_17_R / C_189_R
  c0120  = C_1_R  / C_20_R

  ! c0108  = C_1_R  / C_8_R
  ! c0708  = C_7_R  / C_8_R
  ! c0124  = C_1_R  / C_24_R
  ! c2324  = C_23_R / C_24_R
  ! c1211  = C_12_R / C_11_R

  DO jk =1,jkmax
    ! ! third-order biased (explicit)
    ! d(jk,1)       = - c2324*u(jk,1)    + c0708*u(jk,2)      + c0108*u(jk,3)      - c0124*u(jk,4)
    ! d(jk,imaxp)   =   c2324*u(jk,imax) - c0708*u(jk,imax-1) - c0108*u(jk,imax-2) + c0124*u(jk,imax-3)
    ! forth-order biased (implicit)
    d(jk,1)       = - u(jk,1)    + C_2_R*u(jk,2)      - u(jk,3)
    d(jk,imaxp)   =   u(jk,imax) - C_2_R*u(jk,imax-1) + u(jk,imax-2)
    ! ! forth-order (implicit)
    ! d(jk,2)       =   c1211 * (u(jk,3)       - u(jk,2))
    ! d(jk,imaxp-1) = - c1211 * (u(jk,imax-2) - u(jk,imax-1))
  ENDDO

  ! sixth-order (implicit)  
  DO i = 2,imaxp-1
  ! DO i = 3,imaxp-2 ! for forth order (implict)
    DO jk = 1,jkmax
      d(jk,i) = (u(jk,i+1) - u(jk,i)) + c17189*(u(jk,i+2) - u(jk,i-1))
    ENDDO
  ENDDO

  RETURN
END SUBROUTINE FDM_C1INT6R_RHS

! #######################################################################
! Right-hand side; forcing term ==> interpolation to the left
! #######################################################################
SUBROUTINE FDM_C1INT6L_RHS(imax,imaxp,jkmax, u,d)
  
  IMPLICIT NONE

  TINTEGER,                         INTENT(IN ):: imax, imaxp, jkmax
  TREAL,    DIMENSION(jkmax,imaxp), INTENT(IN ):: u
  TREAL,    DIMENSION(jkmax,imax),  INTENT(OUT):: d

! -------------------------------------------------------------------
  TINTEGER                                     :: i, jk
  TREAL                                        :: c17189, c0120, c1211
  ! TREAL                                        :: c0108, c0708, c0124, c2324

! #######################################################################
  c17189 = C_17_R / C_189_R
  c0120  = C_1_R  / C_20_R
  c1211  = C_12_R / C_11_R

  ! c2324  = C_23_R / C_24_R
  ! c0708  = C_7_R  / C_8_R
  ! c0108  = C_1_R  / C_8_R
  ! c0124  = C_1_R  / C_24_R

  DO jk =1,jkmax
    ! ! third-order biased (explicit)
    ! d(jk,1)      = - c2324*u(jk,1)     + c0708*u(jk,2)       + c0108*u(jk,3)       - c0124*u(jk,4)
    ! d(jk,imax)   =   c2324*u(jk,imaxp) - c0708*u(jk,imaxp-1) - c0108*u(jk,imaxp-2) + c0124*u(jk,imaxp-3)
    ! forth-order biased (implicit)
    d(jk,1)      = - C_25_R*u(jk,1)     + C_26_R*u(jk,2)       - u(jk,3)
    d(jk,imax)   =   C_25_R*u(jk,imaxp) - C_26_R*u(jk,imaxp-1) + u(jk,imaxp-2)
    ! forth-order (implicit)
    d(jk,2)      =   c1211 * (u(jk,2)       - u(jk,1))
    d(jk,imax-1) = - c1211 * (u(jk,imaxp-1) - u(jk,imaxp))
  ENDDO

  ! sixth-order (implicit)  
  DO i = 3,imax-2
  ! DO i = 3,imaxp-2 ! for forth order (implict)
    DO jk = 1,jkmax
      d(jk,i) = (u(jk,i) - u(jk,i-1)) + c17189*(u(jk,i+1) - u(jk,i-2))
    ENDDO
  ENDDO

  RETURN
END SUBROUTINE FDM_C1INT6L_RHS