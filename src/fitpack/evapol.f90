      function evapol(tu,nu,tv,nv,c,rad,x,y)

      IMPLICIT NONE

#include "types.h"

!  function program evacir evaluates the function f(x,y) = s(u,v),
!  defined through the transformation
!  x = u*rad(v)*cos(v)    y = u*rad(v)*sin(v)
!  and where s(u,v) is a bicubic spline ( 0<=u<=1 , -pi<=v<=pi ), given
!  in its standard b-spline representation.
!
!  calling sequence:
! f = evapol(tu,nu,tv,nv,c,rad,x,y)
!
!  input parameters:
!   tu    : real array, length nu, which contains the position of the
!   knots in the u-direction.
!   nu    : integer, giving the total number of knots in the u-direction
!   tv    : real array, length nv, which contains the position of the
!   knots in the v-direction.
!   nv    : integer, giving the total number of knots in the v-direction
!   c     : real array, length (nu-4)*(nv-4), which contains the
!   b-spline coefficients.
!   rad   : real function subprogram, defining the boundary of the
!   approximation domain. must be declared external in the
!   calling (sub)-program
!   x,y   : real values.
!   before entry x and y must be set to the co-ordinates of
!   the point where f(x,y) must be evaluated.
!
!  output parameter:
!   f     : real
!   on exit f contains the value of f(x,y)
!
!  other subroutines required:
!bispev,fpbisp,fpbspl
!
!  references :
!de boor c : on calculating with b-splines, j. approximation theory
!        6 (1972) 50-62.
!cox m.g.  : the numerical evaluation of b-splines, j. inst. maths
!        applics 10 (1972) 134-149.
!dierckx p. : curve and surface fitting with splines, monographs on
!         numerical analysis, oxford university press, 1993.
!
!  author :
!p.dierckx
!dept. computer science, k.u.leuven
!celestijnenlaan 200a, b-3001 heverlee, belgium.
!e-mail : Paul.Dierckx@cs.kuleuven.ac.be
!
!  latest update : march 1989
!
!  ..scalar arguments..
      TINTEGER nu,nv
      TREAL x,y
!  ..array arguments..
      TREAL tu(nu),tv(nv),c((nu-4)*(nv-4))
!  ..user specified function
      TREAL rad
!  ..local scalars..
      TINTEGER ier
      TREAL u,v,r,f,one,dist
!  ..local arrays
      TREAL wrk(8)
      TINTEGER iwrk(2)
!  ..function references
      TREAL atan2,sqrt
      TREAL evapol
      TINTEGER i1, i2, i3, i8
!  ..
!  calculate the (u,v)-coordinates of the given point.
      one = 1
      u = 0.
      v = 0.
      dist = x**2+y**2
      if(dist.le.0.) go to 10
      v = atan2(y,x)
      r = rad(v)
      if(r.le.0.) go to 10
      u = sqrt(dist)/r
      if(u.gt.one) u = one
!  evaluate s(u,v)
      i1 = 1
      i2 = 2
      i3 = 3
      i8 = 8
  10  call bispev(tu,nu,tv,nv,c,i3,i3,u,i1,v,i1,f,wrk,i8,iwrk,i2,ier)
      evapol = f
      return
      end

