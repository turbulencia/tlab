      subroutine splev(t,n,c,k,x,y,m,ier)

      IMPLICIT NONE

#include "types.h"

!  subroutine splev evaluates in a number of points x(i),i=1,2,...,m
!  a spline s(x) of degree k, given in its b-spline representation.
!
!  calling sequence:
! call splev(t,n,c,k,x,y,m,ier)
!
!  input parameters:
!t    : array,length n, which contains the position of the knots.
!n    : integer, giving the total number of knots of s(x).
!c    : array,length n, which contains the b-spline coefficients.
!k    : integer, giving the degree of s(x).
!x    : array,length m, which contains the points where s(x) must
!   be evaluated.
!m    : integer, giving the number of points where s(x) must be
!   evaluated.
!
!  output parameter:
!y    : array,length m, giving the value of s(x) at the different
!   points.
!ier  : error flag
!  ier = 0 : normal return
!  ier =10 : invalid input data (see restrictions)
!
!  restrictions:
!m >= 1
!t(k+1) <= x(i) <= x(i+1) <= t(n-k) , i=1,2,...,m-1.
!
!  other subroutines required: fpbspl.
!
!  references :
!de boor c  : on calculating with b-splines, j. approximation theory
!         6 (1972) 50-62.
!cox m.g.   : the numerical evaluation of b-splines, j. inst. maths
!         applics 10 (1972) 134-149.
!dierckx p. : curve and surface fitting with splines, monographs on
!         numerical analysis, oxford university press, 1993.
!
!  author :
!p.dierckx
!dept. computer science, k.u.leuven
!celestijnenlaan 200a, b-3001 heverlee, belgium.
!e-mail : Paul.Dierckx@cs.kuleuven.ac.be
!
!  latest update : march 1987
!
!  ..scalar arguments..
      TINTEGER n,k,m,ier
!  ..array arguments..
      TREAL t(n),c(n),x(m),y(m)
!  ..local scalars..
      TINTEGER i,j,k1,l,ll,l1,nk1
      TREAL arg,sp,tb,te
!  ..local array..
      TREAL h(6)
!  ..
!  before starting computations a data check is made. if the input data
!  are invalid control is immediately repassed to the calling program.
      ier = 10
      if(m-1) 100,30,10
  10  do 20 i=2,m
        if(x(i).lt.x(i-1)) then
           ier = 101
           go to 100
        endif
  20  continue
  30  ier = 0
!  fetch tb and te, the boundaries of the approximation interval.
      k1 = k+1
      nk1 = n-k1
      tb = t(k1)
      te = t(nk1+1)
      l = k1
      l1 = l+1
!  main loop for the different points.
      do 80 i=1,m
!  fetch a new x-value arg.
        arg = x(i)
        if(arg.lt.tb) arg = tb
        if(arg.gt.te) arg = te
!  search for knot interval t(l) <= arg < t(l+1)
  40    if(arg.lt.t(l1) .or. l.eq.nk1) go to 50
        l = l1
        l1 = l+1
        go to 40
!  evaluate the non-zero b-splines at arg.
  50    call fpbspl(t,n,k,arg,l,h)
!  find the value of s(x) at x=arg.
        sp = 0.
        ll = l-k1
        do 60 j=1,k1
          ll = ll+1
          sp = sp+c(ll)*h(j)
  60    continue
        y(i) = sp
  80  continue
 100  return
      end
