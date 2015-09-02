      subroutine curev(idim,t,n,c,nc,k,u,m,x,mx,ier)
!  subroutine curev evaluates in a number of points u(i),i=1,2,...,m
!  a spline curve s(u) of degree k and dimension idim, given in its
!  b-spline representation.
!
!  calling sequence:
! call curev(idim,t,n,c,nc,k,u,m,x,mx,ier)
!
!  input parameters:
!idim : integer, giving the dimension of the spline curve.
!t    : array,length n, which contains the position of the knots.
!n    : integer, giving the total number of knots of s(u).
!c    : array,length nc, which contains the b-spline coefficients.
!nc   : integer, giving the total number of coefficients of s(u).
!k    : integer, giving the degree of s(u).
!u    : array,length m, which contains the points where s(u) must
!   be evaluated.
!m    : integer, giving the number of points where s(u) must be
!   evaluated.
!mx   : integer, giving the dimension of the array x. mx >= m*idim
!
!  output parameters:
!x    : array,length mx,giving the value of s(u) at the different
!   points. x(idim*(i-1)+j) will contain the j-th coordinate
!   of the i-th point on the curve.
!ier  : error flag
!  ier = 0 : normal return
!  ier =10 : invalid input data (see restrictions)
!
!  restrictions:
!m >= 1
!mx >= m*idim
!t(k+1) <= u(i) <= u(i+1) <= t(n-k) , i=1,2,...,m-1.
!
!  other subroutines required: fpbspl.
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
!  latest update : march 1987
!
!  ..scalar arguments..
      integer idim,n,nc,k,m,mx,ier
!  ..array arguments..
      real t(n),c(nc),u(m),x(mx)
!  ..local scalars..
      integer i,j,jj,j1,k1,l,ll,l1,mm,nk1
      real arg,sp,tb,te
!  ..local array..
      real h(6)
!  ..
!  before starting computations a data check is made. if the input data
!  are invalid control is immediately repassed to the calling program.
      ier = 10
      if(m-1) 100,30,10
  10  do 20 i=2,m
        if(u(i).lt.u(i-1)) go to 100
  20  continue
  30  if(mx.lt.(m*idim)) go to 100
      ier = 0
!  fetch tb and te, the boundaries of the approximation interval.
      k1 = k+1
      nk1 = n-k1
      tb = t(k1)
      te = t(nk1+1)
      l = k1
      l1 = l+1
!  main loop for the different points.
      mm = 0
      do 80 i=1,m
!  fetch a new u-value arg.
        arg = u(i)
        if(arg.lt.tb) arg = tb
        if(arg.gt.te) arg = te
!  search for knot interval t(l) <= arg < t(l+1)
  40    if(arg.lt.t(l1) .or. l.eq.nk1) go to 50
        l = l1
        l1 = l+1
        go to 40
!  evaluate the non-zero b-splines at arg.
  50    call fpbspl(t,n,k,arg,l,h)
!  find the value of s(u) at u=arg.
        ll = l-k1
        do 70 j1=1,idim
          jj = ll
          sp = 0.
          do 60 j=1,k1
            jj = jj+1
            sp = sp+c(jj)*h(j)
  60      continue
          mm = mm+1
          x(mm) = sp
          ll = ll+n
  70    continue
  80  continue
 100  return
      end
