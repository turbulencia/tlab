      subroutine percur(iopt,m,x,y,w,k,s,nest,n,t,c,fp,&
       wrk,lwrk,iwrk,ier)
      
      IMPLICIT NONE

#include "types.h"

!  given the set of data points (x(i),y(i)) and the set of positive
!  numbers w(i),i=1,2,...,m-1, subroutine percur determines a smooth
!  periodic spline approximation of degree k with period per=x(m)-x(1).
!  if iopt=-1 percur calculates the weighted least-squares periodic
!  spline according to a given set of knots.
!  if iopt>=0 the number of knots of the spline s(x) and the position
!  t(j),j=1,2,...,n is chosen automatically by the routine. the smooth-
!  ness of s(x) is then achieved by minimalizing the discontinuity
!  jumps of the k-th derivative of s(x) at the knots t(j),j=k+2,k+3,...,
!  n-k-1. the amount of smoothness is determined by the condition that
!  f(p)=sum((w(i)*(y(i)-s(x(i))))**2) be <= s, with s a given non-
!  negative constant, called the smoothing factor.
!  the fit s(x) is given in the b-spline representation (b-spline coef-
!  ficients c(j),j=1,2,...,n-k-1) and can be evaluated by means of
!  subroutine splev.
!
!  calling sequence:
! call percur(iopt,m,x,y,w,k,s,nest,n,t,c,fp,wrk,
!* lwrk,iwrk,ier)
!
!  parameters:
!   iopt  : integer flag. on entry iopt must specify whether a weighted
!   least-squares spline (iopt=-1) or a smoothing spline (iopt=
!   0 or 1) must be determined. if iopt=0 the routine will start
!   with an initial set of knots t(i)=x(1)+(x(m)-x(1))*(i-k-1),
!   i=1,2,...,2*k+2. if iopt=1 the routine will continue with
!   the knots found at the last call of the routine.
!   attention: a call with iopt=1 must always be immediately
!   preceded by another call with iopt=1 or iopt=0.
!   unchanged on exit.
!   m     : integer. on entry m must specify the number of data points.
!   m > 1. unchanged on exit.
!   x     : real array of dimension at least (m). before entry, x(i)
!   must be set to the i-th value of the independent variable x,
!   for i=1,2,...,m. these values must be supplied in strictly
!   ascending order. x(m) only indicates the length of the
!   period of the spline, i.e per=x(m)-x(1).
!   unchanged on exit.
!   y     : real array of dimension at least (m). before entry, y(i)
!   must be set to the i-th value of the dependent variable y,
!   for i=1,2,...,m-1. the element y(m) is not used.
!   unchanged on exit.
!   w     : real array of dimension at least (m). before entry, w(i)
!   must be set to the i-th value in the set of weights. the
!   w(i) must be strictly positive. w(m) is not used.
!   see also further comments. unchanged on exit.
!   k     : integer. on entry k must specify the degree of the spline.
!   1<=k<=5. it is recommended to use cubic splines (k=3).
!   the user is strongly dissuaded from choosing k even,together
!   with a small s-value. unchanged on exit.
!   s     : real.on entry (in case iopt>=0) s must specify the smoothing
!   factor. s >=0. unchanged on exit.
!   for advice on the choice of s see further comments.
!   nest  : integer. on entry nest must contain an over-estimate of the
!   total number of knots of the spline returned, to indicate
!   the storage space available to the routine. nest >=2*k+2.
!   in most practical situation nest=m/2 will be sufficient.
!   always large enough is nest=m+2*k,the number of knots needed
!   for interpolation (s=0). unchanged on exit.
!   n     : integer.
!   unless ier = 10 (in case iopt >=0), n will contain the
!   total number of knots of the spline approximation returned.
!   if the computation mode iopt=1 is used this value of n
!   should be left unchanged between subsequent calls.
!   in case iopt=-1, the value of n must be specified on entry.
!   t     : real array of dimension at least (nest).
!   on succesful exit, this array will contain the knots of the
!   spline,i.e. the position of the interior knots t(k+2),t(k+3)
!   ...,t(n-k-1) as well as the position of the additional knots
!   t(1),t(2),...,t(k+1)=x(1) and t(n-k)=x(m),..,t(n) needed for
!   the b-spline representation.
!   if the computation mode iopt=1 is used, the values of t(1),
!   t(2),...,t(n) should be left unchanged between subsequent
!   calls. if the computation mode iopt=-1 is used, the values
!   t(k+2),...,t(n-k-1) must be supplied by the user, before
!   entry. see also the restrictions (ier=10).
!   c     : real array of dimension at least (nest).
!   on succesful exit, this array will contain the coefficients
!   c(1),c(2),..,c(n-k-1) in the b-spline representation of s(x)
!   fp    : real. unless ier = 10, fp contains the weighted sum of
!   squared residuals of the spline approximation returned.
!   wrk   : real array of dimension at least (m*(k+1)+nest*(8+5*k)).
!   used as working space. if the computation mode iopt=1 is
!   used, the values wrk(1),...,wrk(n) should be left unchanged
!   between subsequent calls.
!   lwrk  : integer. on entry,lwrk must specify the actual dimension of
!   the array wrk as declared in the calling (sub)program. lwrk
!   must not be too small (see wrk). unchanged on exit.
!   iwrk  : integer array of dimension at least (nest).
!   used as working space. if the computation mode iopt=1 is
!   used,the values iwrk(1),...,iwrk(n) should be left unchanged
!   between subsequent calls.
!   ier   : integer. unless the routine detects an error, ier contains a
!   non-positive value on exit, i.e.
!ier=0  : normal return. the spline returned has a residual sum of
!     squares fp such that abs(fp-s)/s <= tol with tol a relat-
!     ive tolerance set to 0.001 by the program.
!ier=-1 : normal return. the spline returned is an interpolating
!     periodic spline (fp=0).
!ier=-2 : normal return. the spline returned is the weighted least-
!     squares constant. in this extreme case fp gives the upper
!     bound fp0 for the smoothing factor s.
!ier=1  : error. the required storage space exceeds the available
!     storage space, as specified by the parameter nest.
!     probably causes : nest too small. if nest is already
!     large (say nest > m/2), it may also indicate that s is
!     too small
!     the approximation returned is the least-squares periodic
!     spline according to the knots t(1),t(2),...,t(n). (n=nest)
!     the parameter fp gives the corresponding weighted sum of
!     squared residuals (fp>s).
!ier=2  : error. a theoretically impossible result was found during
!     the iteration proces for finding a smoothing spline with
!     fp = s. probably causes : s too small.
!     there is an approximation returned but the corresponding
!     weighted sum of squared residuals does not satisfy the
!     condition abs(fp-s)/s < tol.
!ier=3  : error. the maximal number of iterations maxit (set to 20
!     by the program) allowed for finding a smoothing spline
!     with fp=s has been reached. probably causes : s too small
!     there is an approximation returned but the corresponding
!     weighted sum of squared residuals does not satisfy the
!     condition abs(fp-s)/s < tol.
!ier=10 : error. on entry, the input data are controlled on validity
!     the following restrictions must be satisfied.
!     -1<=iopt<=1, 1<=k<=5, m>1, nest>2*k+2, w(i)>0,i=1,...,m-1
!     x(1)<x(2)<...<x(m), lwrk>=(k+1)*m+nest*(8+5*k)
!     if iopt=-1: 2*k+2<=n<=min(nest,m+2*k)
!                 x(1)<t(k+2)<t(k+3)<...<t(n-k-1)<x(m)
!               the schoenberg-whitney conditions, i.e. there
!               must be a subset of data points xx(j) with
!               xx(j) = x(i) or x(i)+(x(m)-x(1)) such that
!                 t(j) < xx(j) < t(j+k+1), j=k+1,...,n-k-1
!     if iopt>=0: s>=0
!                 if s=0 : nest >= m+2*k
!     if one of these conditions is found to be violated,control
!     is immediately repassed to the calling program. in that
!     case there is no approximation returned.
!
!  further comments:
!   by means of the parameter s, the user can control the tradeoff
!   between closeness of fit and smoothness of fit of the approximation.
!   if s is too large, the spline will be too smooth and signal will be
!   lost ; if s is too small the spline will pick up too much noise. in
!   the extreme cases the program will return an interpolating periodic
!   spline if s=0 and the weighted least-squares constant if s is very
!   large. between these extremes, a properly chosen s will result in
!   a good compromise between closeness of fit and smoothness of fit.
!   to decide whether an approximation, corresponding to a certain s is
!   satisfactory the user is highly recommended to inspect the fits
!   graphically.
!   recommended values for s depend on the weights w(i). if these are
!   taken as 1/d(i) with d(i) an estimate of the standard deviation of
!   y(i), a good s-value should be found in the range (m-sqrt(2*m),m+
!   sqrt(2*m)). if nothing is known about the statistical error in y(i)
!   each w(i) can be set equal to one and s determined by trial and
!   error, taking account of the comments above. the best is then to
!   start with a very large value of s ( to determine the least-squares
!   constant and the corresponding upper bound fp0 for s) and then to
!   progressively decrease the value of s ( say by a factor 10 in the
!   beginning, i.e. s=fp0/10, fp0/100,...and more carefully as the
!   approximation shows more detail) to obtain closer fits.
!   to economize the search for a good s-value the program provides with
!   different modes of computation. at the first call of the routine, or
!   whenever he wants to restart with the initial set of knots the user
!   must set iopt=0.
!   if iopt=1 the program will continue with the set of knots found at
!   the last call of the routine. this will save a lot of computation
!   time if percur is called repeatedly for different values of s.
!   the number of knots of the spline returned and their location will
!   depend on the value of s and on the complexity of the shape of the
!   function underlying the data. but, if the computation mode iopt=1
!   is used, the knots returned may also depend on the s-values at
!   previous calls (if these were smaller). therefore, if after a number
!   of trials with different s-values and iopt=1, the user can finally
!   accept a fit as satisfactory, it may be worthwhile for him to call
!   percur once more with the selected value for s but now with iopt=0.
!   indeed, percur may then return an approximation of the same quality
!   of fit but with fewer knots and therefore better if data reduction
!   is also an important objective for the user.
!
!  other subroutines required:
!fpbacp,fpbspl,fpchep,fpperi,fpdisc,fpgivs,fpknot,fprati,fprota
!
!  references:
!   dierckx p. : algorithms for smoothing data with periodic and
!        parametric splines, computer graphics and image
!        processing 20 (1982) 171-184.
!   dierckx p. : algorithms for smoothing data with periodic and param-
!        etric splines, report tw55, dept. computer science,
!        k.u.leuven, 1981.
!   dierckx p. : curve and surface fitting with splines, monographs on
!        numerical analysis, oxford university press, 1993.
!
!  author:
!p.dierckx
!dept. computer science, k.u. leuven
!celestijnenlaan 200a, b-3001 heverlee, belgium.
!e-mail : Paul.Dierckx@cs.kuleuven.ac.be
!
!  creation date : may 1979
!  latest update : march 1987
!
!  ..
!  ..scalar arguments..
      TREAL s,fp
      TINTEGER iopt,m,k,nest,n,lwrk,ier
!  ..array arguments..
      TREAL x(m),y(m),w(m),t(nest),c(nest),wrk(lwrk)
      TINTEGER iwrk(nest)
!  ..local scalars..
      TREAL per,tol
      TINTEGER i,ia1,ia2,ib,ifp,ig1,ig2,iq,iz,i1,i2,j1,j2,k1,k2,&
           lwest, maxit,m1,nmin
!  ..subroutine references..
!perper,pcheck
!  ..
!  we set up the parameters tol and maxit
      maxit = 20
      tol = C_1EM3_R
!  before starting computations a data check is made. if the input data
!  are invalid, control is immediately repassed to the calling program.
      ier = 10
      if(k.le.0 .or. k.gt.5) go to 50
      k1 = k+1
      k2 = k1+1
      if(iopt.lt.(-1) .or. iopt.gt.1) go to 50
      nmin = 2*k1
      if(m.lt.2 .or. nest.lt.nmin) go to 50
      lwest = m*k1+nest*(8+5*k)
      if(lwrk.lt.lwest) go to 50
      m1 = m-1
      do 10 i=1,m1
         if(x(i).ge.x(i+1) .or. w(i).le.0.) go to 50
  10  continue
      if(iopt.ge.0) go to 30
      if(n.le.nmin .or. n.gt.nest) go to 50
      per = x(m)-x(1)
      j1 = k1
      t(j1) = x(1)
      i1 = n-k
      t(i1) = x(m)
      j2 = j1
      i2 = i1
      do 20 i=1,k
         i1 = i1+1
         i2 = i2-1
         j1 = j1+1
         j2 = j2-1
         t(j2) = t(i2)-per
         t(i1) = t(j1)+per
  20  continue
      call fpchep(x,m,t,n,k,ier)
      if(ier) 50,40,50
  30  if(s.lt.0.) go to 50
      if(s.eq.0. .and. nest.lt.(m+2*k)) go to 50
      ier = 0
! we partition the working space and determine the spline approximation.
  40  ifp = 1
      iz = ifp+nest
      ia1 = iz+nest
      ia2 = ia1+nest*k1
      ib = ia2+nest*k
      ig1 = ib+nest*k2
      ig2 = ig1+nest*k2
      iq = ig2+nest*k1
      call fpperi(iopt,x,y,w,m,k,s,nest,tol,maxit,k1,k2,n,t,c,fp,&
       wrk(ifp),wrk(iz),wrk(ia1),wrk(ia2),wrk(ib),wrk(ig1),wrk(ig2),&
       wrk(iq),iwrk,ier)
  50  return
      end
