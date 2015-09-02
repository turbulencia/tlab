!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2008/08/07 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Solve LU of nonadiagonal system
!#
!# First equation is divided by d1 to get 1 in the first diagonal elem.
!#
!########################################################################
!# ARGUMENTS 
!#
!# nmax    In     Size of the pentadiagonal system
!# len     In     Number of simultaneous systems to be solved
!# frc     In     Array with forcing term
!#         Out    Array with solution
!#
!########################################################################
#include "types.h"

! #######################################################################
! LU factorization stage
! #######################################################################
SUBROUTINE NONADFS(nmax, inorm, a, b, c, d, e, f, g, h, i)

  IMPLICIT NONE

  TINTEGER nmax, inorm
  TREAL, DIMENSION(nmax), INTENT(INOUT) :: a,b,c,d,e,f,g,h,i

! -----------------------------------------------------------------------
  TINTEGER n

! #######################################################################
  a(1) = C_1_R ! padding
  b(1) = C_1_R ! padding
  c(1) = C_1_R ! padding
  IF ( inorm .EQ. 1 ) THEN
     d(1) = C_1_R/e(1) ! padding, and used in nonadss to normalize 1st eqn.
     e(1) = C_1_R
     f(1) = f(1)*d(1)
     g(1) = g(1)*d(1)
     h(1) = h(1)*d(1)
     i(1) = i(1)*d(1)
  ELSE
     d(1) = C_1_R ! padding
  ENDIF

!  n = 2 
  a(2) = C_1_R ! padding
  b(2) = C_1_R ! padding
  c(2) = C_1_R ! padding
  d(2) = d(2)             /e(1)
  e(2) = e(2) - d(2)*f(1)
  f(2) = f(2) - d(2)*g(1)
  g(2) = g(2) - d(2)*h(1)
  h(2) = h(2) - d(2)*i(1)

  n = 3 
  a(n) = C_1_R ! padding
  b(n) = C_1_R ! padding
  c(n) = c(n)                              /e(n-2)
  d(n) =(d(n)               - c(n)*f(n-2) )/e(n-1)
  e(n) = e(n) - d(n)*f(n-1) - c(n)*g(n-2) 
  f(n) = f(n) - d(n)*g(n-1) - c(n)*h(n-2) 
  g(n) = g(n) - d(n)*h(n-1) - c(n)*i(n-2)
  h(n) = h(n) - d(n)*i(n-1)

  n = 4
  a(n) = C_1_R ! padding
  b(n) = b(n)                                            /e(n-3)
  c(n) =(c(n)                             - b(n)*f(n-3) )/e(n-2)
  d(n) =(d(n)               - c(n)*f(n-2) - b(n)*g(n-3) )/e(n-1)
  e(n) = e(n) - d(n)*f(n-1) - c(n)*g(n-2) - b(n)*h(n-3) 
  f(n) = f(n) - d(n)*g(n-1) - c(n)*h(n-2) - b(n)*i(n-3)
  g(n) = g(n) - d(n)*h(n-1) - c(n)*i(n-2)
  h(n) = h(n) - d(n)*i(n-1)

! -----------------------------------------------------------------------
  DO n = 5,nmax-3
     a(n) = a(n)                                                         /e(n-4)
     b(n) =(b(n)                                           - a(n)*f(n-4))/e(n-3)
     c(n) =(c(n)                             - b(n)*f(n-3) - a(n)*g(n-4))/e(n-2)
     d(n) =(d(n)               - c(n)*f(n-2) - b(n)*g(n-3) - a(n)*h(n-4))/e(n-1)
     e(n) = e(n) - d(n)*f(n-1) - c(n)*g(n-2) - b(n)*h(n-3) - a(n)*i(n-4)
     f(n) = f(n) - d(n)*g(n-1) - c(n)*h(n-2) - b(n)*i(n-3)
     g(n) = g(n) - d(n)*h(n-1) - c(n)*i(n-2)
     h(n) = h(n) - d(n)*i(n-1)
  ENDDO
  i(n-1) = C_1_R ! padding

! -----------------------------------------------------------------------
  n = nmax-2
  a(n) = a(n)                                                         /e(n-4)
  b(n) =(b(n)                                           - a(n)*f(n-4))/e(n-3)
  c(n) =(c(n)                             - b(n)*f(n-3) - a(n)*g(n-4))/e(n-2)
  d(n) =(d(n)               - c(n)*f(n-2) - b(n)*g(n-3) - a(n)*h(n-4))/e(n-1)
  e(n) = e(n) - d(n)*f(n-1) - c(n)*g(n-2) - b(n)*h(n-3) - a(n)*i(n-4)
  f(n) = f(n) - d(n)*g(n-1) - c(n)*h(n-2) - b(n)*i(n-3)
  g(n) = g(n) - d(n)*h(n-1) - c(n)*i(n-2)
  h(n) = C_1_R ! padding
  i(n) = C_1_R ! padding

  n = nmax-1
  a(n) = a(n)                                                         /e(n-4)
  b(n) =(b(n)                                           - a(n)*f(n-4))/e(n-3)
  c(n) =(c(n)                             - b(n)*f(n-3) - a(n)*g(n-4))/e(n-2)
  d(n) =(d(n)               - c(n)*f(n-2) - b(n)*g(n-3) - a(n)*h(n-4))/e(n-1)
  e(n) = e(n) - d(n)*f(n-1) - c(n)*g(n-2) - b(n)*h(n-3) - a(n)*i(n-4)
  f(n) = f(n) - d(n)*g(n-1) - c(n)*h(n-2) - b(n)*i(n-3)
  g(n) = C_1_R ! padding
  h(n) = C_1_R ! padding
  i(n) = C_1_R ! padding

  n = nmax
  a(n) = a(n)                                                         /e(n-4)
  b(n) =(b(n)                                           - a(n)*f(n-4))/e(n-3)
  c(n) =(c(n)                             - b(n)*f(n-3) - a(n)*g(n-4))/e(n-2)
  d(n) =(d(n)               - c(n)*f(n-2) - b(n)*g(n-3) - a(n)*h(n-4))/e(n-1)
  e(n) = e(n) - d(n)*f(n-1) - c(n)*g(n-2) - b(n)*h(n-3) - a(n)*i(n-4)
  f(n) = C_1_R ! padding
  g(n) = C_1_R ! padding
  h(n) = C_1_R ! padding
  i(n) = C_1_R ! padding

  RETURN
END SUBROUTINE NONADFS

! #######################################################################
! Backward substitution step in the Thomas algorith
! #######################################################################
SUBROUTINE NONADSS(nmax, len, a, b, c, d, e, f, g, h, i, frc)

  IMPLICIT NONE

  TINTEGER nmax, len
  TREAL, DIMENSION(nmax),     INTENT(IN)    :: a,b,c,d,e,f,g,h,i
  TREAL, DIMENSION(len,nmax), INTENT(INOUT) :: frc

! -----------------------------------------------------------------------
  TINTEGER n, ij

! #######################################################################
! -----------------------------------------------------------------------
! Solve Ly=frc, forward
! -----------------------------------------------------------------------
  DO ij = 1,len
     frc(ij,1) =            frc(ij,1)*d(1) ! Normalize first eqn. See NONADFS
     frc(ij,2) = frc(ij,2) -frc(ij,1)*d(2)
     frc(ij,3) = frc(ij,3) -frc(ij,2)*d(3) -frc(ij,1)*c(3)
     frc(ij,4) = frc(ij,4) -frc(ij,3)*d(4) -frc(ij,2)*c(4) -frc(ij,1)*b(4)
  ENDDO

  DO n = 5,nmax
     DO ij = 1,len
        frc(ij,n) = frc(ij,n) -frc(ij,n-1)*d(n) -frc(ij,n-2)*c(n) -frc(ij,n-3)*b(n) -frc(ij,n-4)*a(n)
     ENDDO
  ENDDO

! -----------------------------------------------------------------------
! Solve Ux=y, backward
! -----------------------------------------------------------------------
  DO ij = 1,len
     frc(ij,nmax  ) = frc(ij,nmax  ) &
           /e(nmax  )
     frc(ij,nmax-1) =(frc(ij,nmax-1) -frc(ij,nmax  )*f(nmax-1) &
          )/e(nmax-1)
     frc(ij,nmax-2) =(frc(ij,nmax-2) -frc(ij,nmax-1)*f(nmax-2) -frc(ij,nmax  )*g(nmax-2) &
          )/e(nmax-2)
     frc(ij,nmax-3) =(frc(ij,nmax-3) -frc(ij,nmax-2)*f(nmax-3) -frc(ij,nmax-1)*g(nmax-3) -frc(ij,nmax)*h(nmax-3) &
          )/e(nmax-3)
  ENDDO

  DO n = nmax-4,1,-1
     DO ij = 1,len
        frc(ij,n) =(frc(ij,n) -frc(ij,n+1)*f(n) -frc(ij,n+2)*g(n) -frc(ij,n+3)*h(n) -frc(ij,n+4)*i(n))/e(n)
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE NONADSS
