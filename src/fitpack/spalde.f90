      subroutine spalde(t,n,c,k1,x,d,ier)
!  subroutine spalde evaluates at a point x all the derivatives
!      (j-1)
!  d(j) = s     (x) , j=1,2,...,k1
!  of a spline s(x) of order k1 (degree k=k1-1), given in its b-spline
!  representation.
!
!  calling sequence:
! call spalde(t,n,c,k1,x,d,ier)
!
!  input parameters:
!t    : array,length n, which contains the position of the knots.
!n    : integer, giving the total number of knots of s(x).
!c    : array,length n, which contains the b-spline coefficients.
!k1   : integer, giving the order of s(x) (order=degree+1)
!x    : real, which contains the point where the derivatives must
!   be evaluated.
!
!  output parameters:
!d    : array,length k1, containing the derivative values of s(x).
!ier  : error flag
!  ier = 0 : normal return
!  ier =10 : invalid input data (see restrictions)
!
!  restrictions:
!t(k1) <= x <= t(n-k1+1)
!
!  further comments:
!if x coincides with a knot, right derivatives are computed
!( left derivatives if x = t(n-k1+1) ).
!
!  other subroutines required: fpader.
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
      integer n,k1,ier
      real x
!  ..array arguments..
      real t(n),c(n),d(k1)
!  ..local scalars..
      integer l,nk1
!  ..
!  before starting computations a data check is made. if the input data
!  are invalid control is immediately repassed to the calling program.
      ier = 10
      nk1 = n-k1
      if(x.lt.t(k1) .or. x.gt.t(nk1+1)) go to 300
!  search for knot interval t(l) <= x < t(l+1)
      l = k1
 100  if(x.lt.t(l+1) .or. l.eq.nk1) go to 200
      l = l+1
      go to 100
 200  if(t(l).ge.t(l+1)) go to 300
      ier = 0
!  calculate the derivatives.
      call fpader(t,n,c,k1,x,l,d)
 300  return
      end
