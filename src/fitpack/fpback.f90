      subroutine fpback(a,z,n,k,c,nest)
!  subroutine fpback calculates the solution of the system of
!  equations a*c = z with a a n x n upper triangular matrix
!  of bandwidth k.
!  ..

      IMPLICIT NONE

#include "types.h"

!  ..scalar arguments..
      TINTEGER n,k,nest
!  ..array arguments..
      TREAL a(nest,k),z(n),c(n)
!  ..local scalars..
      TREAL store
      TINTEGER i,i1,j,k1,l,m
!  ..
      k1 = k-1
      c(n) = z(n)/a(n,1)
      i = n-1
      if(i.eq.0) go to 30
      do 20 j=2,n
        store = z(i)
        i1 = k1
        if(j.le.k1) i1 = j-1
        m = i
        do 10 l=1,i1
          m = m+1
          store = store-c(m)*a(i,l+1)
  10    continue
        c(i) = store/a(i,1)
        i = i-1
  20  continue
  30  return
      end
