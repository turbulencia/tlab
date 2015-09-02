      subroutine fpcuro(a,b,c,d,x,n)
!  subroutine fpcuro finds the real zeros of a cubic polynomial
!  p(x) = a*x**3+b*x**2+c*x+d.
!
!  calling sequence:
! call fpcuro(a,b,c,d,x,n)
!
!  input parameters:
!a,b,c,d: real values, containing the coefficients of p(x).
!
!  output parameters:
!x      : real array,length 3, which contains the real zeros of p(x)
!n      : integer, giving the number of real zeros of p(x).
!  ..
!  ..scalar arguments..
      real a,b,c,d
      integer n
!  ..array argument..
      real x(3)
!  ..local scalars..
      integer i
      real a1,b1,c1,df,disc,d1,e3,f,four,half,ovfl,pi3,p3,q,r,&
       step,tent,three,two,u,u1,u2,y
!  ..function references..
      real abs,amax1,atan,atan2,cos,sign,sqrt
!  set constants
      two = 0.2e+01
      three = 0.3e+01
      four = 0.4e+01
      ovfl =0.1e+05
      half = 0.5e+0
      tent = 0.1e+0
      e3 = tent/0.3e0
      pi3 = atan(0.1e+01)/0.75e0
      a1 = abs(a)
      b1 = abs(b)
      c1 = abs(c)
      d1 = abs(d)
!  test whether p(x) is a third degree polynomial.
      if(amax1(b1,c1,d1).lt.a1*ovfl) go to 300
!  test whether p(x) is a second degree polynomial.
      if(amax1(c1,d1).lt.b1*ovfl) go to 200
!  test whether p(x) is a first degree polynomial.
      if(d1.lt.c1*ovfl) go to 100
!  p(x) is a constant function.
      n = 0
      go to 800
!  p(x) is a first degree polynomial.
 100  n = 1
      x(1) = -d/c
      go to 500
!  p(x) is a second degree polynomial.
 200  disc = c*c-four*b*d
      n = 0
      if(disc.lt.0.) go to 800
      n = 2
      u = sqrt(disc)
      b1 = b+b
      x(1) = (-c+u)/b1
      x(2) = (-c-u)/b1
      go to 500
!  p(x) is a third degree polynomial.
 300  b1 = b/a*e3
      c1 = c/a
      d1 = d/a
      q = c1*e3-b1*b1
      r = b1*b1*b1+(d1-b1*c1)*half
      disc = q*q*q+r*r
      if(disc.gt.0.) go to 400
      u = sqrt(abs(q))
      if(r.lt.0.) u = -u
      p3 = atan2(sqrt(-disc),abs(r))*e3
      u2 = u+u
      n = 3
      x(1) = -u2*cos(p3)-b1
      x(2) = u2*cos(pi3-p3)-b1
      x(3) = u2*cos(pi3+p3)-b1
      go to 500
 400  u = sqrt(disc)
      u1 = -r+u
      u2 = -r-u
      n = 1
      x(1) = sign(abs(u1)**e3,u1)+sign(abs(u2)**e3,u2)-b1
!  apply a newton iteration to improve the accuracy of the roots.
 500  do 700 i=1,n
        y = x(i)
        f = ((a*y+b)*y+c)*y+d
        df = (three*a*y+two*b)*y+c
        step = 0.
        if(abs(f).lt.abs(df)*tent) step = f/df
        x(i) = y-step
 700  continue
 800  return
      end
