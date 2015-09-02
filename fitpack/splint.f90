      real function splint(t,n,c,k,a,b,wrk)
!  function splint calculates the integral of a spline function s(x)
!  of degree k, which is given in its normalized b-spline representation
!
!  calling sequence:
! aint = splint(t,n,c,k,a,b,wrk)
!
!  input parameters:
!t    : array,length n,which contains the position of the knots
!   of s(x).
!n    : integer, giving the total number of knots of s(x).
!c    : array,length n, containing the b-spline coefficients.
!k    : integer, giving the degree of s(x).
!a,b  : real values, containing the end points of the integration
!   interval. s(x) is considered to be identically zero outside
!   the interval (t(k+1),t(n-k)).
!
!  output parameter:
!aint : real, containing the integral of s(x) between a and b.
!wrk  : real array, length n.  used as working space
!   on output, wrk will contain the integrals of the normalized
!   b-splines defined on the set of knots.
!
!  other subroutines required: fpintb.
!
!  references :
!gaffney p.w. : the calculation of indefinite integrals of b-splines
!           j. inst. maths applics 17 (1976) 37-41.
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
      real a,b
      integer n,k
!  ..array arguments..
      real t(n),c(n),wrk(n)
!  ..local scalars..
      integer i,nk1
!  ..
      nk1 = n-k-1
!  calculate the integrals wrk(i) of the normalized b-splines
!  ni,k+1(x), i=1,2,...nk1.
      call fpintb(t,n,wrk,nk1,a,b)
!  calculate the integral of s(x).
      splint = 0.
      do 10 i=1,nk1
        splint = splint+c(i)*wrk(i)
  10  continue
      return
      end
