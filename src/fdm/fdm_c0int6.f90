!########################################################################
!# Tool/Library PADE
!#
!########################################################################
!# HISTORY
!#
!# 2022/01/14 - J. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Implementation of the finite difference interpolation with
!# 6th order tridiagonal compact scheme by JCP Lele 1992, non-periodic.
!# Interior points according to Eq.C.1.4 (\alpha=3/10, \beta=0, c=0).
!#
!# Different boundary closures can be found in Albin 2010
!# (https://doi.org/10.1002/fld.2520) table 6 (typos in Si2_5i and Si2_6i).
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
! ==> interpolation from velocity to pressure grid
! #######################################################################
SUBROUTINE FDM_C0INTVP6_LHS(imaxp, a,b,c)
  
  IMPLICIT NONE

  TINTEGER,                   INTENT(IN ):: imaxp ! pressure grid (imaxp==imax-1)
  TREAL,    DIMENSION(imaxp), INTENT(OUT):: a,b,c

! -------------------------------------------------------------------
  TINTEGER                              :: i

! ###################################################################
! forth-order biased (implicit)
  a(1)  = C_0_R
  b(1)  = C_1_R
  c(1)  = C_1_R
  !
  a(imaxp) = C_1_R
  b(imaxp) = C_1_R
  c(imaxp) = C_0_R

! sixth-order (implicit)
  DO i = 2,imaxp-1
    a(i) = C_3_R / C_10_R
    b(i) = C_1_R 
    c(i) = C_3_R / C_10_R
  ENDDO

! ! forth-order (implicit) - for adjacent boundary nodes 
!   a(2)      = C_1_R / C_6_R
!   b(2)      = C_1_R 
!   c(2)      = C_1_R / C_6_R
!   !
!   a(imaxp-1) = C_1_R / C_6_R
!   b(imaxp-1) = C_1_R
!   c(imaxp-1) = C_1_R / C_6_R
  
  RETURN
END SUBROUTINE FDM_C0INTVP6_LHS

! #######################################################################
! Left-hand side; tridiagonal matrix of the linear system
! ==> interpolation from pressure to velocity grid
! #######################################################################
SUBROUTINE FDM_C0INTPV6_LHS(imax, a,b,c)
  
  IMPLICIT NONE

  TINTEGER,                  INTENT(IN ):: imax ! velocity grid
  TREAL,    DIMENSION(imax), INTENT(OUT):: a,b,c

! -------------------------------------------------------------------
  TINTEGER                              :: i

! ###################################################################
! third-order biased (implicit)
  a(1)  = C_0_R
  b(1)  = C_1_R
  c(1)  = C_3_R
  !
  a(imax) = C_3_R
  b(imax) = C_1_R
  c(imax) = C_0_R
  
! forth-order (implicit)
  a(2)      = C_1_R / C_6_R
  b(2)      = C_1_R 
  c(2)      = C_1_R / C_6_R
  !
  a(imax-1) = C_1_R / C_6_R
  b(imax-1) = C_1_R
  c(imax-1) = C_1_R / C_6_R

! sixth-order (implicit)
  DO i = 3,imax-2
    a(i) = C_3_R / C_10_R
    b(i) = C_1_R 
    c(i) = C_3_R / C_10_R
  ENDDO

  RETURN
END SUBROUTINE FDM_C0INTPV6_LHS

! #######################################################################
! Right-hand side; forcing term
! ==> interpolation from velocity to pressure grid
! #######################################################################
SUBROUTINE FDM_C0INTVP6_RHS(imax,imaxp,jkmax, u,d)
  
  IMPLICIT NONE

  TINTEGER,                         INTENT(IN ):: imax, imaxp, jkmax
  TREAL,    DIMENSION(jkmax,imax),  INTENT(IN ):: u
  TREAL,    DIMENSION(jkmax,imaxp), INTENT(OUT):: d

! -------------------------------------------------------------------
  TINTEGER                                     :: i, jk
  TREAL                                        :: c0302, c0104, c0304, c0120 !, c0203

! #######################################################################
  c0302 = C_3_R / C_2_R
  c0104 = C_1_R / C_4_R
  c0304 = C_3_R / C_4_R
  c0120 = C_1_R / C_20_R
  ! c0203 = C_2_R / C_3_R

  DO jk =1,jkmax
    ! forth-order biased (implicit)
    d(jk,1)       = c0302*u(jk,2)      + c0104*(u(jk,1)    + u(jk,3)     )
    d(jk,imaxp)   = c0302*u(jk,imax-1) + c0104*(u(jk,imax) + u(jk,imax-2))
    ! forth-order (implicit)
    ! d(jk,2)       = c0203*(u(jk,3)      + u(jk,2)     )
    ! d(jk,imaxp-1) = c0203*(u(jk,imax-2) + u(jk,imax-1))
  ENDDO
  
! sixth-order (implicit)  
  DO i = 2,imaxp-1
  ! DO i = 3,imaxp-2 ! for forth order (implict)
    DO jk = 1,jkmax
      d(jk,i) = c0304*(u(jk,i+1) + u(jk,i)) + c0120*(u(jk,i+2) + u(jk,i-1))
    ENDDO
  ENDDO

  RETURN
END SUBROUTINE FDM_C0INTVP6_RHS

! #######################################################################
! Right-hand side; forcing term
! ==> interpolation from pressure to velocity grid
! #######################################################################
SUBROUTINE FDM_C0INTPV6_RHS(imax,imaxp,jkmax, u,d)
  
  IMPLICIT NONE

  TINTEGER,                         INTENT(IN ):: imax, imaxp, jkmax
  TREAL,    DIMENSION(jkmax,imaxp), INTENT(IN ):: u
  TREAL,    DIMENSION(jkmax,imax),  INTENT(OUT):: d

! -------------------------------------------------------------------
  TINTEGER                                     :: i, jk
  TREAL                                        :: c0304, c0120, c0203

! #######################################################################
  c0304 = C_3_R / C_4_R
  c0120 = C_1_R / C_20_R
  c0203 = C_2_R / C_3_R

  DO jk =1,jkmax
    ! forth-order biased (implicit)
    d(jk,1)      = C_3_R*u(jk,1)      + u(jk,2)   
    d(jk,imax)   = C_3_R*u(jk,imaxp)  + u(jk,imaxp-1)
    ! forth-order (implicit)
    d(jk,2)      = c0203*(u(jk,1)     + u(jk,2)      )
    d(jk,imax-1) = c0203*(u(jk,imaxp) + u(jk,imaxp-1))
  ENDDO

! sixth-order (implicit)
  DO i = 3,imax-2
    DO jk = 1,jkmax
      d(jk,i) = c0304*(u(jk,i) + u(jk,i-1)) + c0120*(u(jk,i+1) + u(jk,i-2))
    ENDDO
  ENDDO

  RETURN
END SUBROUTINE FDM_C0INTPV6_RHS