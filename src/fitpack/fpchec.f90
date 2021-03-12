      subroutine fpchec(x,m,t,n,k,ier)
!  subroutine fpchec verifies the number and the position of the knots
!  t(j),j=1,2,...,n of a spline of degree k, in relation to the number
!  and the position of the data points x(i),i=1,2,...,m. if all of the
!  following conditions are fulfilled, the error parameter ier is set
!  to zero. if one of the conditions is violated ier is set to ten.
!  1) k+1 <= n-k-1 <= m
!  2) t(1) <= t(2) <= ... <= t(k+1)
! t(n-k) <= t(n-k+1) <= ... <= t(n)
!  3) t(k+1) < t(k+2) < ... < t(n-k)
!  4) t(k+1) <= x(i) <= t(n-k)
!  5) the conditions specified by schoenberg and whitney must hold
! for at least one subset of data points, i.e. there must be a
! subset of data points y(j) such that
!     t(j) < y(j) < t(j+k+1), j=1,2,...,n-k-1
!  ..

      IMPLICIT NONE

#include "types.h"

!  ..scalar arguments..
      TINTEGER m,n,k,ier
!  ..array arguments..
      TREAL x(m),t(n)
!  ..local scalars..
      TINTEGER i,j,k1,k2,l,nk1,nk2,nk3
      TREAL tj,tl
!  ..
      k1 = k+1
      k2 = k1+1
      nk1 = n-k1
      nk2 = nk1+1
      ier = 10
!  check condition no 1
      if(nk1.lt.k1 .or. nk1.gt.m) go to 80
!  check condition no 2
      j = n
      do 20 i=1,k
        if(t(i).gt.t(i+1)) go to 80
        if(t(j).lt.t(j-1)) go to 80
        j = j-1
  20  continue
!  check condition no 3
      do 30 i=k2,nk2
        if(t(i).le.t(i-1)) go to 80
  30  continue
!  check condition no 4
      if(x(1).lt.t(k1) .or. x(m).gt.t(nk2)) go to 80
!  check condition no 5
      if(x(1).ge.t(k2) .or. x(m).le.t(nk1)) go to 80
      i = 1
      l = k2
      nk3 = nk1-1
      if(nk3.lt.2) go to 70
      do 60 j=2,nk3
        tj = t(j)
        l = l+1
        tl = t(l)
  40    i = i+1
        if(i.ge.m) go to 80
        if(x(i).le.tj) go to 40
        if(x(i).ge.tl) go to 80
  60  continue
  70  ier = 0
  80  return
      end
