      subroutine concur(iopt,idim,m,u,mx,x,xx,w,ib,db,nb,ie,de,ne,k,s,&
       nest,n,t,nc,c,np,cp,fp,wrk,lwrk,iwrk,ier)
!  given the ordered set of m points x(i) in the idim-dimensional space
!  and given also a corresponding set of strictly increasing values u(i)
!  and the set of positive numbers w(i),i=1,2,...,m, subroutine concur
!  determines a smooth approximating spline curve s(u), i.e.
!  x1 = s1(u)
!  x2 = s2(u)      ub = u(1) <= u <= u(m) = ue
!  .........
!  xidim = sidim(u)
!  with sj(u),j=1,2,...,idim spline functions of odd degree k with
!  common knots t(j),j=1,2,...,n.
!  in addition these splines will satisfy the following boundary
!  constraints        (l)
!  if ib > 0 :  sj   (u(1)) = db(idim*l+j) ,l=0,1,...,ib-1
!  and                (l)
!  if ie > 0 :  sj   (u(m)) = de(idim*l+j) ,l=0,1,...,ie-1.
!  if iopt=-1 concur calculates the weighted least-squares spline curve
!  according to a given set of knots.
!  if iopt>=0 the number of knots of the splines sj(u) and the position
!  t(j),j=1,2,...,n is chosen automatically by the routine. the smooth-
!  ness of s(u) is then achieved by minimalizing the discontinuity
!  jumps of the k-th derivative of s(u) at the knots t(j),j=k+2,k+3,...,
!  n-k-1. the amount of smoothness is determined by the condition that
!  f(p)=sum((w(i)*dist(x(i),s(u(i))))**2) be <= s, with s a given non-
!  negative constant, called the smoothing factor.
!  the fit s(u) is given in the b-spline representation and can be
!  evaluated by means of subroutine curev.
!
!  calling sequence:
! call concur(iopt,idim,m,u,mx,x,xx,w,ib,db,nb,ie,de,ne,k,s,nest,n,
!* t,nc,c,np,cp,fp,wrk,lwrk,iwrk,ier)
!
!  parameters:
!   iopt  : integer flag. on entry iopt must specify whether a weighted
!   least-squares spline curve (iopt=-1) or a smoothing spline
!   curve (iopt=0 or 1) must be determined.if iopt=0 the routine
!   will start with an initial set of knots t(i)=ub,t(i+k+1)=ue,
!   i=1,2,...,k+1. if iopt=1 the routine will continue with the
!   knots found at the last call of the routine.
!   attention: a call with iopt=1 must always be immediately
!   preceded by another call with iopt=1 or iopt=0.
!   unchanged on exit.
!   idim  : integer. on entry idim must specify the dimension of the
!   curve. 0 < idim < 11.
!   unchanged on exit.
!   m     : integer. on entry m must specify the number of data points.
!   m > k-max(ib-1,0)-max(ie-1,0). unchanged on exit.
!   u     : real array of dimension at least (m). before entry,
!   u(i) must be set to the i-th value of the parameter variable
!   u for i=1,2,...,m. these values must be supplied in
!   strictly ascending order and will be unchanged on exit.
!   mx    : integer. on entry mx must specify the actual dimension of
!   the arrays x and xx as declared in the calling (sub)program
!   mx must not be too small (see x). unchanged on exit.
!   x     : real array of dimension at least idim*m.
!   before entry, x(idim*(i-1)+j) must contain the j-th coord-
!   inate of the i-th data point for i=1,2,...,m and j=1,2,...,
!   idim. unchanged on exit.
!   xx    : real array of dimension at least idim*m.
!   used as working space. on exit xx contains the coordinates
!   of the data points to which a spline curve with zero deriv-
!   ative constraints has been determined.
!   if the computation mode iopt =1 is used xx should be left
!   unchanged between calls.
!   w     : real array of dimension at least (m). before entry, w(i)
!   must be set to the i-th value in the set of weights. the
!   w(i) must be strictly positive. unchanged on exit.
!   see also further comments.
!   ib    : integer. on entry ib must specify the number of derivative
!   constraints for the curve at the begin point. 0<=ib<=(k+1)/2
!   unchanged on exit.
!   db    : real array of dimension nb. before entry db(idim*l+j) must
!   contain the l-th order derivative of sj(u) at u=u(1) for
!   j=1,2,...,idim and l=0,1,...,ib-1 (if ib>0).
!   unchanged on exit.
!   nb    : integer, specifying the dimension of db. nb>=max(1,idim*ib)
!   unchanged on exit.
!   ie    : integer. on entry ie must specify the number of derivative
!   constraints for the curve at the end point. 0<=ie<=(k+1)/2
!   unchanged on exit.
!   de    : real array of dimension ne. before entry de(idim*l+j) must
!   contain the l-th order derivative of sj(u) at u=u(m) for
!   j=1,2,...,idim and l=0,1,...,ie-1 (if ie>0).
!   unchanged on exit.
!   ne    : integer, specifying the dimension of de. ne>=max(1,idim*ie)
!   unchanged on exit.
!   k     : integer. on entry k must specify the degree of the splines.
!   k=1,3 or 5.
!   unchanged on exit.
!   s     : real.on entry (in case iopt>=0) s must specify the smoothing
!   factor. s >=0. unchanged on exit.
!   for advice on the choice of s see further comments.
!   nest  : integer. on entry nest must contain an over-estimate of the
!   total number of knots of the splines returned, to indicate
!   the storage space available to the routine. nest >=2*k+2.
!   in most practical situation nest=m/2 will be sufficient.
!   always large enough is nest=m+k+1+max(0,ib-1)+max(0,ie-1),
!   the number of knots needed for interpolation (s=0).
!   unchanged on exit.
!   n     : integer.
!   unless ier = 10 (in case iopt >=0), n will contain the
!   total number of knots of the smoothing spline curve returned
!   if the computation mode iopt=1 is used this value of n
!   should be left unchanged between subsequent calls.
!   in case iopt=-1, the value of n must be specified on entry.
!   t     : real array of dimension at least (nest).
!   on succesful exit, this array will contain the knots of the
!   spline curve,i.e. the position of the interior knots t(k+2),
!   t(k+3),..,t(n-k-1) as well as the position of the additional
!   t(1)=t(2)=...=t(k+1)=ub and t(n-k)=...=t(n)=ue needed for
!   the b-spline representation.
!   if the computation mode iopt=1 is used, the values of t(1),
!   t(2),...,t(n) should be left unchanged between subsequent
!   calls. if the computation mode iopt=-1 is used, the values
!   t(k+2),...,t(n-k-1) must be supplied by the user, before
!   entry. see also the restrictions (ier=10).
!   nc    : integer. on entry nc must specify the actual dimension of
!   the array c as declared in the calling (sub)program. nc
!   must not be too small (see c). unchanged on exit.
!   c     : real array of dimension at least (nest*idim).
!   on succesful exit, this array will contain the coefficients
!   in the b-spline representation of the spline curve s(u),i.e.
!   the b-spline coefficients of the spline sj(u) will be given
!   in c(n*(j-1)+i),i=1,2,...,n-k-1 for j=1,2,...,idim.
!   cp    : real array of dimension at least 2*(k+1)*idim.
!   on exit cp will contain the b-spline coefficients of a
!   polynomial curve which satisfies the boundary constraints.
!   if the computation mode iopt =1 is used cp should be left
!   unchanged between calls.
!   np    : integer. on entry np must specify the actual dimension of
!   the array cp as declared in the calling (sub)program. np
!   must not be too small (see cp). unchanged on exit.
!   fp    : real. unless ier = 10, fp contains the weighted sum of
!   squared residuals of the spline curve returned.
!   wrk   : real array of dimension at least m*(k+1)+nest*(6+idim+3*k).
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
!ier=0  : normal return. the curve returned has a residual sum of
!     squares fp such that abs(fp-s)/s <= tol with tol a relat-
!     ive tolerance set to 0.001 by the program.
!ier=-1 : normal return. the curve returned is an interpolating
!     spline curve, satisfying the constraints (fp=0).
!ier=-2 : normal return. the curve returned is the weighted least-
!     squares polynomial curve of degree k, satisfying the
!     constraints. in this extreme case fp gives the upper
!     bound fp0 for the smoothing factor s.
!ier=1  : error. the required storage space exceeds the available
!     storage space, as specified by the parameter nest.
!     probably causes : nest too small. if nest is already
!     large (say nest > m/2), it may also indicate that s is
!     too small
!     the approximation returned is the least-squares spline
!     curve according to the knots t(1),t(2),...,t(n). (n=nest)
!     the parameter fp gives the corresponding weighted sum of
!     squared residuals (fp>s).
!ier=2  : error. a theoretically impossible result was found during
!     the iteration proces for finding a smoothing spline curve
!     with fp = s. probably causes : s too small.
!     there is an approximation returned but the corresponding
!     weighted sum of squared residuals does not satisfy the
!     condition abs(fp-s)/s < tol.
!ier=3  : error. the maximal number of iterations maxit (set to 20
!     by the program) allowed for finding a smoothing curve
!     with fp=s has been reached. probably causes : s too small
!     there is an approximation returned but the corresponding
!     weighted sum of squared residuals does not satisfy the
!     condition abs(fp-s)/s < tol.
!ier=10 : error. on entry, the input data are controlled on validity
!     the following restrictions must be satisfied.
!     -1<=iopt<=1, k = 1,3 or 5, m>k-max(0,ib-1)-max(0,ie-1),
!     nest>=2k+2, 0<idim<=10, lwrk>=(k+1)*m+nest*(6+idim+3*k),
!     nc >=nest*idim ,u(1)<u(2)<...<u(m),w(i)>0 i=1,2,...,m,
!     mx>=idim*m,0<=ib<=(k+1)/2,0<=ie<=(k+1)/2,nb>=1,ne>=1,
!     nb>=ib*idim,ne>=ib*idim,np>=2*(k+1)*idim,
!     if iopt=-1:2*k+2<=n<=min(nest,mmax) with mmax = m+k+1+
!                max(0,ib-1)+max(0,ie-1)
!                u(1)<t(k+2)<t(k+3)<...<t(n-k-1)<u(m)
!               the schoenberg-whitney conditions, i.e. there
!               must be a subset of data points uu(j) such that
!                 t(j) < uu(j) < t(j+k+1), j=1+max(0,ib-1),...
!                                           ,n+k-1-max(0,ie-1)
!     if iopt>=0: s>=0
!                 if s=0 : nest >=mmax (see above)
!     if one of these conditions is found to be violated,control
!     is immediately repassed to the calling program. in that
!     case there is no approximation returned.
!
!  further comments:
!   by means of the parameter s, the user can control the tradeoff
!   between closeness of fit and smoothness of fit of the approximation.
!   if s is too large, the curve will be too smooth and signal will be
!   lost ; if s is too small the curve will pick up too much noise. in
!   the extreme cases the program will return an interpolating curve if
!   s=0 and the least-squares polynomial curve of degree k if s is
!   very large. between these extremes, a properly chosen s will result
!   in a good compromise between closeness of fit and smoothness of fit.
!   to decide whether an approximation, corresponding to a certain s is
!   satisfactory the user is highly recommended to inspect the fits
!   graphically.
!   recommended values for s depend on the weights w(i). if these are
!   taken as 1/d(i) with d(i) an estimate of the standard deviation of
!   x(i), a good s-value should be found in the range (m-sqrt(2*m),m+
!   sqrt(2*m)). if nothing is known about the statistical error in x(i)
!   each w(i) can be set equal to one and s determined by trial and
!   error, taking account of the comments above. the best is then to
!   start with a very large value of s ( to determine the least-squares
!   polynomial curve and the upper bound fp0 for s) and then to
!   progressively decrease the value of s ( say by a factor 10 in the
!   beginning, i.e. s=fp0/10, fp0/100,...and more carefully as the
!   approximating curve shows more detail) to obtain closer fits.
!   to economize the search for a good s-value the program provides with
!   different modes of computation. at the first call of the routine, or
!   whenever he wants to restart with the initial set of knots the user
!   must set iopt=0.
!   if iopt=1 the program will continue with the set of knots found at
!   the last call of the routine. this will save a lot of computation
!   time if concur is called repeatedly for different values of s.
!   the number of knots of the spline returned and their location will
!   depend on the value of s and on the complexity of the shape of the
!   curve underlying the data. but, if the computation mode iopt=1 is
!   used, the knots returned may also depend on the s-values at previous
!   calls (if these were smaller). therefore, if after a number of
!   trials with different s-values and iopt=1, the user can finally
!   accept a fit as satisfactory, it may be worthwhile for him to call
!   concur once more with the selected value for s but now with iopt=0.
!   indeed, concur may then return an approximation of the same quality
!   of fit but with fewer knots and therefore better if data reduction
!   is also an important objective for the user.
!
!   the form of the approximating curve can strongly be affected by
!   the choice of the parameter values u(i). if there is no physical
!   reason for choosing a particular parameter u, often good results
!   will be obtained with the choice
!v(1)=0, v(i)=v(i-1)+q(i), i=2,...,m, u(i)=v(i)/v(m), i=1,..,m
!   where
!q(i)= sqrt(sum j=1,idim (xj(i)-xj(i-1))**2 )
!   other possibilities for q(i) are
!q(i)= sum j=1,idim (xj(i)-xj(i-1))**2
!q(i)= sum j=1,idim abs(xj(i)-xj(i-1))
!q(i)= max j=1,idim abs(xj(i)-xj(i-1))
!q(i)= 1
!
!  other subroutines required:
!fpback,fpbspl,fpched,fpcons,fpdisc,fpgivs,fpknot,fprati,fprota
!curev,fppocu,fpadpo,fpinst
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
      real s,fp
      integer iopt,idim,m,mx,ib,nb,ie,ne,k,nest,n,nc,np,lwrk,ier
!  ..array arguments..
      real u(m),x(mx),xx(mx),db(nb),de(ne),w(m),t(nest),c(nc),wrk(lwrk)
      real cp(np)
      integer iwrk(nest)
!  ..local scalars..
      real tol !,dist
      integer i,ib1,ie1,ja,jb,jfp,jg,jq,jz,j,k1,k2,lwest,maxit,nmin,&
       ncc,kk,mmin,nmax,mxx
! ..function references
      integer max0
!  ..
!  we set up the parameters tol and maxit
      maxit = 20
      tol = 0.1e-02
!  before starting computations a data check is made. if the input data
!  are invalid, control is immediately repassed to the calling program.
      ier = 10
      if(iopt.lt.(-1) .or. iopt.gt.1) go to 90
      if(idim.le.0 .or. idim.gt.10) go to 90
      if(k.le.0 .or. k.gt.5) go to 90
      k1 = k+1
      kk = k1/2
      if(kk*2.ne.k1) go to 90
      k2 = k1+1
      if(ib.lt.0 .or. ib.gt.kk) go to 90
      if(ie.lt.0 .or. ie.gt.kk) go to 90
      nmin = 2*k1
      ib1 = max0(0,ib-1)
      ie1 = max0(0,ie-1)
      mmin = k1-ib1-ie1
      if(m.lt.mmin .or. nest.lt.nmin) go to 90
      if(nb.lt.(idim*ib) .or. ne.lt.(idim*ie)) go to 90
      if(np.lt.(2*k1*idim)) go to 90
      mxx = m*idim
      ncc = nest*idim
      if(mx.lt.mxx .or. nc.lt.ncc) go to 90
      lwest = m*k1+nest*(6+idim+3*k)
      if(lwrk.lt.lwest) go to 90
      if(w(1).le.0.) go to 90
      do 10 i=2,m
         if(u(i-1).ge.u(i) .or. w(i).le.0.) go to 90
  10  continue
      if(iopt.ge.0) go to 30
      if(n.lt.nmin .or. n.gt.nest) go to 90
      j = n
      do 20 i=1,k1
         t(i) = u(1)
         t(j) = u(m)
         j = j-1
  20  continue
      call fpched(u,m,t,n,k,ib,ie,ier)
      if(ier) 90,40,90
  30  if(s.lt.0.) go to 90
      nmax = m+k1+ib1+ie1
      if(s.eq.0. .and. nest.lt.nmax) go to 90
      ier = 0
      if(iopt.gt.0) go to 70
!  we determine a polynomial curve satisfying the boundary constraints.
  40  call fppocu(idim,k,u(1),u(m),ib,db,nb,ie,de,ne,cp,np)
!  we generate new data points which will be approximated by a spline
!  with zero derivative constraints.
      j = nmin
      do 50 i=1,k1
        wrk(i) = u(1)
        wrk(j) = u(m)
        j = j-1
  50  continue
!  evaluate the polynomial curve
      call curev(idim,wrk,nmin,cp,np,k,u,m,xx,mxx,ier)
!  substract from the old data, the values of the polynomial curve
      do 60 i=1,mxx
        xx(i) = x(i)-xx(i)
  60  continue
! we partition the working space and determine the spline curve.
  70  jfp = 1
      jz = jfp+nest
      ja = jz+ncc
      jb = ja+nest*k1
      jg = jb+nest*k2
      jq = jg+nest*k2
      call fpcons(iopt,idim,m,u,mxx,xx,w,ib,ie,k,s,nest,tol,maxit,k1,&
       k2,n,t,ncc,c,fp,wrk(jfp),wrk(jz),wrk(ja),wrk(jb),wrk(jg),wrk(jq),&
       iwrk,ier)
!  add the polynomial curve to the calculated spline.
      call fpadpo(idim,t,n,c,ncc,k,cp,np,wrk(jz),wrk(ja),wrk(jb))
  90  return
      end
