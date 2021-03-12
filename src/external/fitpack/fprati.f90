      real function fprati(p1,f1,p2,f2,p3,f3)
!  given three points (p1,f1),(p2,f2) and (p3,f3), function fprati
!  gives the value of p such that the rational interpolating function
!  of the form r(p) = (u*p+v)/(p+w) equals zero at p.
!  ..

      IMPLICIT NONE

#include "types.h"

!  ..scalar arguments..
      TREAL p1,f1,p2,f2,p3,f3
!  ..local scalars..
      TREAL h1,h2,h3,p
!  ..
      if(p3.gt.0.) go to 10
!  value of p in case p3 = infinity.
      p = (p1*(f1-f3)*f2-p2*(f2-f3)*f1)/((f1-f2)*f3)
      go to 20
!  value of p in case p3 ^= infinity.
  10  h1 = f1*(f2-f3)
      h2 = f2*(f3-f1)
      h3 = f3*(f1-f2)
      p = -(p1*p2*h3+p2*p3*h1+p3*p1*h2)/(p1*h1+p2*h2+p3*h3)
!  adjust the value of p1,f1,p3 and f3 such that f1 > 0 and f3 < 0.
  20  if(f2.lt.0.) go to 30
      p1 = p2
      f1 = f2
      go to 40
  30  p3 = p2
      f3 = f2
  40  fprati = p
      return
      end
