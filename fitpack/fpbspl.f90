      subroutine fpbspl(t,n,k,x,l,h)
!  subroutine fpbspl evaluates the (k+1) non-zero b-splines of
!  degree k at t(l) <= x < t(l+1) using the stable recurrence
!  relation of de boor and cox.
!  ..

      IMPLICIT NONE

#include "types.h"

!  ..scalar arguments..
      TREAL x
      TINTEGER n,k,l
!  ..array arguments..
      TREAL t(n),h(6)
!  ..local scalars..
      TREAL f,one
      TINTEGER i,j,li,lj
!  ..local arrays..
      TREAL hh(5)
!  ..
      one = C_1_R
      h(1) = one
      do 20 j=1,k
        do 10 i=1,j
          hh(i) = h(i)
  10    continue
        h(1) = C_0_R
        do 20 i=1,j
          li = l+i
          lj = li-j
          f = hh(i)/(t(li)-t(lj))
          h(i) = h(i)+f*(t(li)-x)
          h(i+1) = f*(x-t(lj))
  20  continue
      return
      end
