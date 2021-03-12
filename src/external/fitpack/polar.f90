      subroutine polar(iopt,m,x,y,z,w,rad,s,nuest,nvest,eps,nu,tu,&
        nv,tv,u,v,c,fp,wrk1,lwrk1,wrk2,lwrk2,iwrk,kwrk,ier)
!  subroutine polar fits a smooth function f(x,y) to a set of data
!  points (x(i),y(i),z(i)) scattered arbitrarily over an approximation
!  domain  x**2+y**2 <= rad(atan(y/x))**2. through the transformation
!x = u*rad(v)*cos(v) , y = u*rad(v)*sin(v)
!  the approximation problem is reduced to the determination of a bi-
!  cubic spline s(u,v) fitting a corresponding set of data points
!  (u(i),v(i),z(i)) on the rectangle 0<=u<=1,-pi<=v<=pi.
!  in order to have continuous partial derivatives
!      i+j
!     d   f(0,0)
!g(i,j) = ----------
!        i   j
!      dx  dy
!
!  s(u,v)=f(x,y) must satisfy the following conditions
!
!(1) s(0,v) = g(0,0)   -pi <=v<= pi.
!
!d s(0,v)
!(2) -------- = rad(v)*(cos(v)*g(1,0)+sin(v)*g(0,1))
!d u
!                                            -pi <=v<= pi
! 2
!d s(0,v)         2       2             2
!(3) -------- = rad(v)*(cos(v)*g(2,0)+sin(v)*g(0,2)+sin(2*v)*g(1,1))
!   2
!d u                                         -pi <=v<= pi
!
!  moreover, s(u,v) must be periodic in the variable v, i.e.
!
! j            j
!d s(u,-pi)   d s(u,pi)
!(4) ---------- = ---------   0 <=u<= 1, j=0,1,2
!   j           j
!d v         d v
!
!  if iopt(1) < 0 circle calculates a weighted least-squares spline
!  according to a given set of knots in u- and v- direction.
!  if iopt(1) >=0, the number of knots in each direction and their pos-
!  ition tu(j),j=1,2,...,nu ; tv(j),j=1,2,...,nv are chosen automatical-
!  ly by the routine. the smoothness of s(u,v) is then achieved by mini-
!  malizing the discontinuity jumps of the derivatives of the spline
!  at the knots. the amount of smoothness of s(u,v) is determined  by
!  the condition that fp = sum((w(i)*(z(i)-s(u(i),v(i))))**2) be <= s,
!  with s a given non-negative constant.
!  the bicubic spline is given in its standard b-spline representation
!  and the corresponding function f(x,y) can be evaluated by means of
!  function program evapol.
!
! calling sequence:
! call polar(iopt,m,x,y,z,w,rad,s,nuest,nvest,eps,nu,tu,
!*  nv,tv,u,v,wrk1,lwrk1,wrk2,lwrk2,iwrk,kwrk,ier)
!
! parameters:
!  iopt  : integer array of dimension 3, specifying different options.
!  unchanged on exit.
!  iopt(1):on entry iopt(1) must specify whether a weighted
!  least-squares polar spline (iopt(1)=-1) or a smoothing
!  polar spline (iopt(1)=0 or 1) must be determined.
!  if iopt(1)=0 the routine will start with an initial set of
!  knots tu(i)=0,tu(i+4)=1,i=1,...,4;tv(i)=(2*i-9)*pi,i=1,...,8.
!  if iopt(1)=1 the routine will continue with the set of knots
!  found at the last call of the routine.
!  attention: a call with iopt(1)=1 must always be immediately
!  preceded by another call with iopt(1) = 1 or iopt(1) = 0.
!  iopt(2):on entry iopt(2) must specify the requested order of conti-
!  nuity for f(x,y) at the origin.
!  if iopt(2)=0 only condition (1) must be fulfilled,
!  if iopt(2)=1 conditions (1)+(2) must be fulfilled and
!  if iopt(2)=2 conditions (1)+(2)+(3) must be fulfilled.
!  iopt(3):on entry iopt(3) must specify whether (iopt(3)=1) or not
!  (iopt(3)=0) the approximation f(x,y) must vanish at the
!  boundary of the approximation domain.
!  m     : integer. on entry m must specify the number of data points.
!  m >= 4-iopt(2)-iopt(3) unchanged on exit.
!  x     : real array of dimension at least (m).
!  y     : real array of dimension at least (m).
!  z     : real array of dimension at least (m).
!  before entry, x(i),y(i),z(i) must be set to the co-ordinates
!  of the i-th data point, for i=1,...,m. the order of the data
!  points is immaterial. unchanged on exit.
!  w     : real array of dimension at least (m). before entry, w(i) must
!  be set to the i-th value in the set of weights. the w(i) must
!  be strictly positive. unchanged on exit.
!  rad   : real function subprogram defining the boundary of the approx-
!  imation domain, i.e   x = rad(v)*cos(v) , y = rad(v)*sin(v),
!  -pi <= v <= pi.
!  must be declared external in the calling (sub)program.
!  s     : real. on entry (in case iopt(1) >=0) s must specify the
!  smoothing factor. s >=0. unchanged on exit.
!  for advice on the choice of s see further comments
!  nuest : integer. unchanged on exit.
!  nvest : integer. unchanged on exit.
!  on entry, nuest and nvest must specify an upper bound for the
!  number of knots required in the u- and v-directions resp.
!  these numbers will also determine the storage space needed by
!  the routine. nuest >= 8, nvest >= 8.
!  in most practical situation nuest = nvest = 8+sqrt(m/2) will
!  be sufficient. see also further comments.
!  eps   : real.
!  on entry, eps must specify a threshold for determining the
!  effective rank of an over-determined linear system of equat-
!  ions. 0 < eps < 1.  if the number of decimal digits in the
!  computer representation of a real number is q, then 10**(-q)
!  is a suitable value for eps in most practical applications.
!  unchanged on exit.
!  nu    : integer.
!  unless ier=10 (in case iopt(1) >=0),nu will contain the total
!  number of knots with respect to the u-variable, of the spline
!  approximation returned. if the computation mode iopt(1)=1
!  is used, the value of nu should be left unchanged between
!  subsequent calls.
!  in case iopt(1)=-1,the value of nu must be specified on entry
!  tu    : real array of dimension at least nuest.
!  on succesful exit, this array will contain the knots of the
!  spline with respect to the u-variable, i.e. the position
!  of the interior knots tu(5),...,tu(nu-4) as well as the
!  position of the additional knots tu(1)=...=tu(4)=0 and
!  tu(nu-3)=...=tu(nu)=1 needed for the b-spline representation
!  if the computation mode iopt(1)=1 is used,the values of
!  tu(1),...,tu(nu) should be left unchanged between subsequent
!  calls. if the computation mode iopt(1)=-1 is used,the values
!  tu(5),...tu(nu-4) must be supplied by the user, before entry.
!  see also the restrictions (ier=10).
!  nv    : integer.
!  unless ier=10 (in case iopt(1)>=0), nv will contain the total
!  number of knots with respect to the v-variable, of the spline
!  approximation returned. if the computation mode iopt(1)=1
!  is used, the value of nv should be left unchanged between
!  subsequent calls. in case iopt(1)=-1, the value of nv should
!  be specified on entry.
!  tv    : real array of dimension at least nvest.
!  on succesful exit, this array will contain the knots of the
!  spline with respect to the v-variable, i.e. the position of
!  the interior knots tv(5),...,tv(nv-4) as well as the position
!  of the additional knots tv(1),...,tv(4) and tv(nv-3),...,
!  tv(nv) needed for the b-spline representation.
!  if the computation mode iopt(1)=1 is used, the values of
!  tv(1),...,tv(nv) should be left unchanged between subsequent
!  calls. if the computation mode iopt(1)=-1 is used,the values
!  tv(5),...tv(nv-4) must be supplied by the user, before entry.
!  see also the restrictions (ier=10).
!  u     : real array of dimension at least (m).
!  v     : real array of dimension at least (m).
!  on succesful exit, u(i),v(i) contains the co-ordinates of
!  the i-th data point with respect to the transformed rectan-
!  gular approximation domain, for i=1,2,...,m.
!  if the computation mode iopt(1)=1 is used the values of
!  u(i),v(i) should be left unchanged between subsequent calls.
!  c     : real array of dimension at least (nuest-4)*(nvest-4).
!  on succesful exit, c contains the coefficients of the spline
!  approximation s(u,v).
!  fp    : real. unless ier=10, fp contains the weighted sum of
!  squared residuals of the spline approximation returned.
!  wrk1  : real array of dimension (lwrk1). used as workspace.
!  if the computation mode iopt(1)=1 is used the value of
!  wrk1(1) should be left unchanged between subsequent calls.
!  on exit wrk1(2),wrk1(3),...,wrk1(1+ncof) will contain the
!  values d(i)/max(d(i)),i=1,...,ncof=1+iopt(2)*(iopt(2)+3)/2+
!  (nv-7)*(nu-5-iopt(2)-iopt(3)) with d(i) the i-th diagonal el-
!  ement of the triangular matrix for calculating the b-spline
!  coefficients.it includes those elements whose square is < eps
!  which are treated as 0 in the case of rank deficiency(ier=-2)
!  lwrk1 : integer. on entry lwrk1 must specify the actual dimension of
!  the array wrk1 as declared in the calling (sub)program.
!  lwrk1 must not be too small. let
!    k = nuest-7, l = nvest-7, p = 1+iopt(2)*(iopt(2)+3)/2,
!    q = k+2-iopt(2)-iopt(3) then
!  lwrk1 >= 129+10*k+21*l+k*l+(p+l*q)*(1+8*l+p)+8*m
!  wrk2  : real array of dimension (lwrk2). used as workspace, but
!  only in the case a rank deficient system is encountered.
!  lwrk2 : integer. on entry lwrk2 must specify the actual dimension of
!  the array wrk2 as declared in the calling (sub)program.
!  lwrk2 > 0 . a save upper bound  for lwrk2 = (p+l*q+1)*(4*l+p)
!  +p+l*q where p,l,q are as above. if there are enough data
!  points, scattered uniformly over the approximation domain
!  and if the smoothing factor s is not too small, there is a
!  good chance that this extra workspace is not needed. a lot
!  of memory might therefore be saved by setting lwrk2=1.
!  (see also ier > 10)
!  iwrk  : integer array of dimension (kwrk). used as workspace.
!  kwrk  : integer. on entry kwrk must specify the actual dimension of
!  the array iwrk as declared in the calling (sub)program.
!  kwrk >= m+(nuest-7)*(nvest-7).
!  ier   : integer. unless the routine detects an error, ier contains a
!  non-positive value on exit, i.e.
!   ier=0  : normal return. the spline returned has a residual sum of
!    squares fp such that abs(fp-s)/s <= tol with tol a relat-
!    ive tolerance set to 0.001 by the program.
!   ier=-1 : normal return. the spline returned is an interpolating
!    spline (fp=0).
!   ier=-2 : normal return. the spline returned is the weighted least-
!    squares constrained polynomial . in this extreme case
!    fp gives the upper bound for the smoothing factor s.
!   ier<-2 : warning. the coefficients of the spline returned have been
!    computed as the minimal norm least-squares solution of a
!    (numerically) rank deficient system. (-ier) gives the rank.
!    especially if the rank deficiency which can be computed as
!    1+iopt(2)*(iopt(2)+3)/2+(nv-7)*(nu-5-iopt(2)-iopt(3))+ier
!    is large the results may be inaccurate.
!    they could also seriously depend on the value of eps.
!   ier=1  : error. the required storage space exceeds the available
!    storage space, as specified by the parameters nuest and
!    nvest.
!    probably causes : nuest or nvest too small. if these param-
!    eters are already large, it may also indicate that s is
!    too small
!    the approximation returned is the weighted least-squares
!    polar spline according to the current set of knots.
!    the parameter fp gives the corresponding weighted sum of
!    squared residuals (fp>s).
!   ier=2  : error. a theoretically impossible result was found during
!    the iteration proces for finding a smoothing spline with
!    fp = s. probably causes : s too small or badly chosen eps.
!    there is an approximation returned but the corresponding
!    weighted sum of squared residuals does not satisfy the
!    condition abs(fp-s)/s < tol.
!   ier=3  : error. the maximal number of iterations maxit (set to 20
!    by the program) allowed for finding a smoothing spline
!    with fp=s has been reached. probably causes : s too small
!    there is an approximation returned but the corresponding
!    weighted sum of squared residuals does not satisfy the
!    condition abs(fp-s)/s < tol.
!   ier=4  : error. no more knots can be added because the dimension
!    of the spline 1+iopt(2)*(iopt(2)+3)/2+(nv-7)*(nu-5-iopt(2)
!    -iopt(3)) already exceeds the number of data points m.
!    probably causes : either s or m too small.
!    the approximation returned is the weighted least-squares
!    polar spline according to the current set of knots.
!    the parameter fp gives the corresponding weighted sum of
!    squared residuals (fp>s).
!   ier=5  : error. no more knots can be added because the additional
!    knot would (quasi) coincide with an old one.
!    probably causes : s too small or too large a weight to an
!    inaccurate data point.
!    the approximation returned is the weighted least-squares
!    polar spline according to the current set of knots.
!    the parameter fp gives the corresponding weighted sum of
!    squared residuals (fp>s).
!   ier=10 : error. on entry, the input data are controlled on validity
!    the following restrictions must be satisfied.
!    -1<=iopt(1)<=1 , 0<=iopt(2)<=2 , 0<=iopt(3)<=1 ,
!    m>=4-iopt(2)-iopt(3) , nuest>=8 ,nvest >=8, 0<eps<1,
!    0<=teta(i)<=pi, 0<=phi(i)<=2*pi, w(i)>0, i=1,...,m
!    lwrk1 >= 129+10*k+21*l+k*l+(p+l*q)*(1+8*l+p)+8*m
!    kwrk >= m+(nuest-7)*(nvest-7)
!    if iopt(1)=-1:9<=nu<=nuest,9+iopt(2)*(iopt(2)+1)<=nv<=nvest
!                  0<tu(5)<tu(6)<...<tu(nu-4)<1
!                  -pi<tv(5)<tv(6)<...<tv(nv-4)<pi
!    if iopt(1)>=0: s>=0
!    if one of these conditions is found to be violated,control
!    is immediately repassed to the calling program. in that
!    case there is no approximation returned.
!   ier>10 : error. lwrk2 is too small, i.e. there is not enough work-
!    space for computing the minimal least-squares solution of
!    a rank deficient system of linear equations. ier gives the
!    requested value for lwrk2. there is no approximation re-
!    turned but, having saved the information contained in nu,
!    nv,tu,tv,wrk1,u,v and having adjusted the value of lwrk2
!    and the dimension of the array wrk2 accordingly, the user
!    can continue at the point the program was left, by calling
!    polar with iopt(1)=1.
!
! further comments:
!  by means of the parameter s, the user can control the tradeoff
!   between closeness of fit and smoothness of fit of the approximation.
!   if s is too large, the spline will be too smooth and signal will be
!   lost ; if s is too small the spline will pick up too much noise. in
!   the extreme cases the program will return an interpolating spline if
!   s=0 and the constrained weighted least-squares polynomial if s is
!   very large. between these extremes, a properly chosen s will result
!   in a good compromise between closeness of fit and smoothness of fit.
!   to decide whether an approximation, corresponding to a certain s is
!   satisfactory the user is highly recommended to inspect the fits
!   graphically.
!   recommended values for s depend on the weights w(i). if these are
!   taken as 1/d(i) with d(i) an estimate of the standard deviation of
!   z(i), a good s-value should be found in the range (m-sqrt(2*m),m+
!   sqrt(2*m)). if nothing is known about the statistical error in z(i)
!   each w(i) can be set equal to one and s determined by trial and
!   error, taking account of the comments above. the best is then to
!   start with a very large value of s ( to determine the least-squares
!   polynomial and the corresponding upper bound fp0 for s) and then to
!   progressively decrease the value of s ( say by a factor 10 in the
!   beginning, i.e. s=fp0/10, fp0/100,...and more carefully as the
!   approximation shows more detail) to obtain closer fits.
!   to choose s very small is strongly discouraged. this considerably
!   increases computation time and memory requirements. it may also
!   cause rank-deficiency (ier<-2) and endager numerical stability.
!   to economize the search for a good s-value the program provides with
!   different modes of computation. at the first call of the routine, or
!   whenever he wants to restart with the initial set of knots the user
!   must set iopt(1)=0.
!   if iopt(1)=1 the program will continue with the set of knots found
!   at the last call of the routine. this will save a lot of computation
!   time if polar is called repeatedly for different values of s.
!   the number of knots of the spline returned and their location will
!   depend on the value of s and on the complexity of the shape of the
!   function underlying the data. if the computation mode iopt(1)=1
!   is used, the knots returned may also depend on the s-values at
!   previous calls (if these were smaller). therefore, if after a number
!   of trials with different s-values and iopt(1)=1,the user can finally
!   accept a fit as satisfactory, it may be worthwhile for him to call
!   polar once more with the selected value for s but now with iopt(1)=0
!   indeed, polar may then return an approximation of the same quality
!   of fit but with fewer knots and therefore better if data reduction
!   is also an important objective for the user.
!   the number of knots may also depend on the upper bounds nuest and
!   nvest. indeed, if at a certain stage in polar the number of knots
!   in one direction (say nu) has reached the value of its upper bound
!   (nuest), then from that moment on all subsequent knots are added
!   in the other (v) direction. this may indicate that the value of
!   nuest is too small. on the other hand, it gives the user the option
!   of limiting the number of knots the routine locates in any direction
!
!  other subroutines required:
!fpback,fpbspl,fppola,fpdisc,fpgivs,fprank,fprati,fprota,fporde,
!fprppo
!
!  references:
!   dierckx p.: an algorithm for fitting data over a circle using tensor
!       product splines,j.comp.appl.maths 15 (1986) 161-173.
!   dierckx p.: an algorithm for fitting data on a circle using tensor
!       product splines, report tw68, dept. computer science,
!       k.u.leuven, 1984.
!   dierckx p.: curve and surface fitting with splines, monographs on
!       numerical analysis, oxford university press, 1993.
!
!  author:
!p.dierckx
!dept. computer science, k.u. leuven
!celestijnenlaan 200a, b-3001 heverlee, belgium.
!e-mail : Paul.Dierckx@cs.kuleuven.ac.be
!
!  creation date : june 1984
!  latest update : march 1989
!
!  ..
!  ..scalar arguments..
      real s,eps,fp
      integer m,nuest,nvest,nu,nv,lwrk1,lwrk2,kwrk,ier
!  ..array arguments..
      real x(m),y(m),z(m),w(m),tu(nuest),tv(nvest),u(m),v(m),&
       c((nuest-4)*(nvest-4)),wrk1(lwrk1),wrk2(lwrk2)
      integer iopt(3),iwrk(kwrk)
!  ..user specified function
      real rad
!  ..local scalars..
      real tol,pi,dist,r,one
      integer i,ib1,ib3,ki,kn,kwest,la,lbu,lcc,lcs,lro,j,&
       lbv,lco,lf,lff,lfp,lh,lq,lsu,lsv,lwest,maxit,ncest,ncc,nuu,&
       nvv,nreg,nrint,nu4,nv4,iopt1,iopt2,iopt3,ipar,nvmin
!  ..function references..
      real atan2,sqrt
      external rad
!  ..subroutine references..
!fppola
!  ..
!  set up constants
      one = 1
!  we set up the parameters tol and maxit.
      maxit = 20
      tol = 0.1e-02
!  before starting computations a data check is made. if the input data
!  are invalid,control is immediately repassed to the calling program.
      ier = 10
      if(eps.le.0. .or. eps.ge.1.) go to 60
      iopt1 = iopt(1)
      if(iopt1.lt.(-1) .or. iopt1.gt.1) go to 60
      iopt2 = iopt(2)
      if(iopt2.lt.0 .or. iopt2.gt.2) go to 60
      iopt3 = iopt(3)
      if(iopt3.lt.0 .or. iopt3.gt.1) go to 60
      if(m.lt.(4-iopt2-iopt3)) go to 60
      if(nuest.lt.8 .or. nvest.lt.8) go to 60
      nu4 = nuest-4
      nv4 = nvest-4
      ncest = nu4*nv4
      nuu = nuest-7
      nvv = nvest-7
      ipar = 1+iopt2*(iopt2+3)/2
      ncc = ipar+nvv*(nuest-5-iopt2-iopt3)
      nrint = nuu+nvv
      nreg = nuu*nvv
      ib1 = 4*nvv
      ib3 = ib1+ipar
      lwest = ncc*(1+ib1+ib3)+2*nrint+ncest+m*8+ib3+5*nuest+12*nvest
      kwest = m+nreg
      if(lwrk1.lt.lwest .or. kwrk.lt.kwest) go to 60
      if(iopt1.gt.0) go to 40
      do 10 i=1,m
        if(w(i).le.0.) go to 60
        dist = x(i)**2+y(i)**2
        u(i) = 0.
        v(i) = 0.
        if(dist.le.0.) go to 10
        v(i) = atan2(y(i),x(i))
        r = rad(v(i))
        if(r.le.0.) go to 60
        u(i) = sqrt(dist)/r
        if(u(i).gt.one) go to 60
  10  continue
      if(iopt1.eq.0) go to 40
      nuu = nu-8
      if(nuu.lt.1 .or. nu.gt.nuest) go to 60
      tu(4) = 0.
      do 20 i=1,nuu
         j = i+4
         if(tu(j).le.tu(j-1) .or. tu(j).ge.one) go to 60
  20  continue
      nvv = nv-8
      nvmin = 9+iopt2*(iopt2+1)
      if(nv.lt.nvmin .or. nv.gt.nvest) go to 60
      pi = atan2(0.,-one)
      tv(4) = -pi
      do 30 i=1,nvv
         j = i+4
         if(tv(j).le.tv(j-1) .or. tv(j).ge.pi) go to 60
  30  continue
      go to 50
  40  if(s.lt.0.) go to 60
  50  ier = 0
!  we partition the working space and determine the spline approximation
      kn = 1
      ki = kn+m
      lq = 2
      la = lq+ncc*ib3
      lf = la+ncc*ib1
      lff = lf+ncc
      lfp = lff+ncest
      lco = lfp+nrint
      lh = lco+nrint
      lbu = lh+ib3
      lbv = lbu+5*nuest
      lro = lbv+5*nvest
      lcc = lro+nvest
      lcs = lcc+nvest
      lsu = lcs+nvest*5
      lsv = lsu+m*4
      call fppola(iopt1,iopt2,iopt3,m,u,v,z,w,rad,s,nuest,nvest,eps,tol,&
       maxit,ib1,ib3,ncest,ncc,nrint,nreg,nu,tu,nv,tv,c,fp,wrk1(1),&
       wrk1(lfp),wrk1(lco),wrk1(lf),wrk1(lff),wrk1(lro),wrk1(lcc),&
       wrk1(lcs),wrk1(la),wrk1(lq),wrk1(lbu),wrk1(lbv),wrk1(lsu),&
       wrk1(lsv),wrk1(lh),iwrk(ki),iwrk(kn),wrk2,lwrk2,ier)
  60  return
      end

