      subroutine splder(t,n,c,k,nu,x,y,m,wrk,ier)
!  subroutine splder evaluates in a number of points x(i),i=1,2,...,m
!  the derivative of order nu of a spline s(x) of degree k,given in
!  its b-spline representation.
!
!  calling sequence:
! call splder(t,n,c,k,nu,x,y,m,wrk,ier)
!
!  input parameters:
!t    : array,length n, which contains the position of the knots.
!n    : integer, giving the total number of knots of s(x).
!c    : array,length n, which contains the b-spline coefficients.
!k    : integer, giving the degree of s(x).
!nu   : integer, specifying the order of the derivative. 0<=nu<=k
!x    : array,length m, which contains the points where the deriv-
!   ative of s(x) must be evaluated.
!m    : integer, giving the number of points where the derivative
!   of s(x) must be evaluated
!wrk  : real array of dimension n. used as working space.
!
!  output parameters:
!y    : array,length m, giving the value of the derivative of s(x)
!   at the different points.
!ier  : error flag
!  ier = 0 : normal return
!  ier =10 : invalid input data (see restrictions)
!
!  restrictions:
!0 <= nu <= k
!m >= 1
!t(k+1) <= x(i) <= x(i+1) <= t(n-k) , i=1,2,...,m-1.
!
!  other subroutines required: fpbspl
!
!  references :
!de boor c : on calculating with b-splines, j. approximation theory
!        6 (1972) 50-62.
!cox m.g.  : the numerical evaluation of b-splines, j. inst. maths
!        applics 10 (1972) 134-149.
!   dierckx p. : curve and surface fitting with splines, monographs on
!        numerical analysis, oxford university press, 1993.
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
      integer n,k,nu,m,ier
!  ..array arguments..
      real t(n),c(n),x(m),y(m),wrk(n)
!  ..local scalars..
      integer i,j,kk,k1,k2,l,ll,l1,l2,nk1,nk2,nn
      real ak,arg,fac,sp,tb,te
!  ..local arrays ..
      real h(6)
!  before starting computations a data check is made. if the input data
!  are invalid control is immediately repassed to the calling program.
      ier = 10
      if(nu.lt.0 .or. nu.gt.k) go to 200
      if(m-1) 200,30,10
  10  do 20 i=2,m
        if(x(i).lt.x(i-1)) go to 200
  20  continue
  30  ier = 0
!  fetch tb and te, the boundaries of the approximation interval.
      k1 = k+1
      nk1 = n-k1
      tb = t(k1)
      te = t(nk1+1)
!  the derivative of order nu of a spline of degree k is a spline of
!  degree k-nu,the b-spline coefficients wrk(i) of which can be found
!  using the recurrence scheme of de boor.
      l = 1
      kk = k
      nn = n
      do 40 i=1,nk1
         wrk(i) = c(i)
  40  continue
      if(nu.eq.0) go to 100
      nk2 = nk1
      do 60 j=1,nu
         ak = kk
         nk2 = nk2-1
         l1 = l
         do 50 i=1,nk2
            l1 = l1+1
            l2 = l1+kk
            fac = t(l2)-t(l1)
            if(fac.le.0.) go to 50
            wrk(i) = ak*(wrk(i+1)-wrk(i))/fac
  50     continue
         l = l+1
         kk = kk-1
  60  continue
      if(kk.ne.0) go to 100
!  if nu=k the derivative is a piecewise constant function
      j = 1
      do 90 i=1,m
         arg = x(i)
  70     if(arg.lt.t(l+1) .or. l.eq.nk1) go to 80
         l = l+1
         j = j+1
         go to 70
  80     y(i) = wrk(j)
  90  continue
      go to 200
 100  l = k1
      l1 = l+1
      k2 = k1-nu
!  main loop for the different points.
      do 180 i=1,m
!  fetch a new x-value arg.
        arg = x(i)
        if(arg.lt.tb) arg = tb
        if(arg.gt.te) arg = te
!  search for knot interval t(l) <= arg < t(l+1)
 140    if(arg.lt.t(l1) .or. l.eq.nk1) go to 150
        l = l1
        l1 = l+1
        go to 140
!  evaluate the non-zero b-splines of degree k-nu at arg.
 150    call fpbspl(t,n,kk,arg,l,h)
!  find the value of the derivative at x=arg.
        sp = 0.
        ll = l-k1
        do 160 j=1,k2
          ll = ll+1
          sp = sp+wrk(ll)*h(j)
 160    continue
        y(i) = sp
 180  continue
 200  return
      end
