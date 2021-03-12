      subroutine insert(iopt,t,n,c,k,x,tt,nn,cc,nest,ier)
!  subroutine insert inserts a new knot x into a spline function s(x)
!  of degree k and calculates the b-spline representation of s(x) with
!  respect to the new set of knots. in addition, if iopt.ne.0, s(x)
!  will be considered as a periodic spline with period per=t(n-k)-t(k+1)
!  satisfying the boundary constraints
!   t(i+n-2*k-1) = t(i)+per  ,i=1,2,...,2*k+1
!   c(i+n-2*k-1) = c(i)      ,i=1,2,...,k
!  in that case, the knots and b-spline coefficients returned will also
!  satisfy these boundary constraints, i.e.
!   tt(i+nn-2*k-1) = tt(i)+per  ,i=1,2,...,2*k+1
!   cc(i+nn-2*k-1) = cc(i)      ,i=1,2,...,k
!
!  calling sequence:
! call insert(iopt,t,n,c,k,x,tt,nn,cc,nest,ier)
!
!  input parameters:
!iopt : integer flag, specifying whether (iopt.ne.0) or not (iopt=0)
!   the given spline must be considered as being periodic.
!t    : array,length nest, which contains the position of the knots.
!n    : integer, giving the total number of knots of s(x).
!c    : array,length nest, which contains the b-spline coefficients.
!k    : integer, giving the degree of s(x).
!x    : real, which gives the location of the knot to be inserted.
!nest : integer specifying the dimension of the arrays t,c,tt and cc
!   nest > n.
!
!  output parameters:
!tt   : array,length nest, which contains the position of the knots
!   after insertion.
!nn   : integer, giving the total number of knots after insertion
!cc   : array,length nest, which contains the b-spline coefficients
!   of s(x) with respect to the new set of knots.
!ier  : error flag
!  ier = 0 : normal return
!  ier =10 : invalid input data (see restrictions)
!
!  restrictions:
!nest > n
!t(k+1) <= x <= t(n-k)
!in case of a periodic spline (iopt.ne.0) there must be
!   either at least k interior knots t(j) satisfying t(k+1)<t(j)<=x
!   or at least k interior knots t(j) satisfying x<=t(j)<t(n-k)
!
!  other subroutines required: fpinst.
!
!  further comments:
!   subroutine insert may be called as follows
!call insert(iopt,t,n,c,k,x,t,n,c,nest,ier)
!   in which case the new representation will simply replace the old one
!
!  references :
!boehm w : inserting new knots into b-spline curves. computer aided
!      design 12 (1980) 199-201.
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
      integer iopt,n,k,nn,nest,ier
      real x
!  ..array arguments..
      real t(nest),c(nest),tt(nest),cc(nest)
!  ..local scalars..
      integer kk,k1,l,nk,nk1
!  ..
!  before starting computations a data check is made. if the input data
!  are invalid control is immediately repassed to the calling program.
      ier = 10
      if(nest.le.n) go to 40
      k1 = k+1
      nk = n-k
      if(x.lt.t(k1) .or. x.gt.t(nk)) go to 40
!  search for knot interval t(l) <= x < t(l+1).
      nk1 = nk-1
      l = k1
  10  if(x.lt.t(l+1) .or. l.eq.nk1) go to 20
      l = l+1
      go to 10
  20  if(t(l).ge.t(l+1)) go to 40
      if(iopt.eq.0) go to 30
      kk = 2*k
      if(l.le.kk .and. l.ge.(n-kk)) go to 40
  30  ier = 0
!  insert the new knot.
      call fpinst(iopt,t,n,c,k,x,l,tt,nn,cc,nest)
  40  return
      end
