      subroutine fpsysy(a,n,g)
! subroutine fpsysy solves a linear n x n symmetric system
!(a) * (b) = (g)
! on input, vector g contains the right hand side ; on output it will
! contain the solution (b).
!  ..

      IMPLICIT NONE

#include "types.h"

!  ..scalar arguments..
      TINTEGER n
!  ..array arguments..
      TREAL a(6,6),g(6)
!  ..local scalars..
      TREAL fac
      TINTEGER i,i1,j,k
!  ..
      g(1) = g(1)/a(1,1)
      if(n.eq.1) return
!  decomposition of the symmetric matrix (a) = (l) * (d) *(l)
!  with (l) a unit lower triangular matrix and (d) a diagonal
!  matrix
      do 10 k=2,n
         a(k,1) = a(k,1)/a(1,1)
  10  continue
      do 40 i=2,n
         i1 = i-1
         do 30 k=i,n
            fac = a(k,i)
            do 20 j=1,i1
               fac = fac-a(j,j)*a(k,j)*a(i,j)
  20        continue
            a(k,i) = fac
            if(k.gt.i) a(k,i) = fac/a(i,i)
  30     continue
  40  continue
!  solve the system (l)*(d)*(l)*(b) = (g).
!  first step : solve (l)*(d)*(c) = (g).
      do 60 i=2,n
         i1 = i-1
         fac = g(i)
         do 50 j=1,i1
            fac = fac-g(j)*a(j,j)*a(i,j)
  50     continue
         g(i) = fac/a(i,i)
  60  continue
!  second step : solve (l)*(b) = (c)
      i = n
      do 80 j=2,n
         i1 = i
         i = i-1
         fac = g(i)
         do 70 k=i1,n
            fac = fac-g(k)*a(k,i)
  70     continue
         g(i) = fac
  80  continue
      return
      end
