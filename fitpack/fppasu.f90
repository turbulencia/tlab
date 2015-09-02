      subroutine fppasu(iopt,ipar,idim,u,mu,v,mv,z,mz,s,nuest,nvest,&
       tol,maxit,nc,nu,tu,nv,tv,c,fp,fp0,fpold,reducu,reducv,fpintu,&
       fpintv,lastdi,nplusu,nplusv,nru,nrv,nrdatu,nrdatv,wrk,lwrk,ier)
!  ..
!  ..scalar arguments..
      real s,tol,fp,fp0,fpold,reducu,reducv
      integer iopt,idim,mu,mv,mz,nuest,nvest,maxit,nc,nu,nv,lastdi,&
       nplusu,nplusv,lwrk,ier
!  ..array arguments..
      real u(mu),v(mv),z(mz*idim),tu(nuest),tv(nvest),c(nc*idim),&
       fpintu(nuest),fpintv(nvest),wrk(lwrk)
      integer ipar(2),nrdatu(nuest),nrdatv(nvest),nru(mu),nrv(mv)
!  ..local scalars
      real acc,fpms,f1,f2,f3,p,p1,p2,p3,rn,one,con1,con9,con4,&
       peru,perv,ub,ue,vb,ve
      integer i,ich1,ich3,ifbu,ifbv,ifsu,ifsv,iter,j,lau1,lav1,laa,&
       l,lau,lav,lbu,lbv,lq,lri,lsu,lsv,l1,l2,l3,l4,mm,mpm,mvnu,ncof,&
       nk1u,nk1v,nmaxu,nmaxv,nminu,nminv,nplu,nplv,npl1,nrintu,&
       nrintv,nue,nuk,nve,nuu,nvv
!  ..function references..
      real abs,fprati
      integer max0,min0
!  ..subroutine references..
!fpgrpa,fpknot
!  ..
!   set constants
      one = 1
      con1 = 0.1e0
      con9 = 0.9e0
      con4 = 0.4e-01
!  set boundaries of the approximation domain
      ub = u(1)
      ue = u(mu)
      vb = v(1)
      ve = v(mv)
!  we partition the working space.
      lsu = 1
      lsv = lsu+mu*4
      lri = lsv+mv*4
      mm = max0(nuest,mv)
      lq = lri+mm*idim
      mvnu = nuest*mv*idim
      lau = lq+mvnu
      nuk = nuest*5
      lbu = lau+nuk
      lav = lbu+nuk
      nuk = nvest*5
      lbv = lav+nuk
      laa = lbv+nuk
      lau1 = lau
      if(ipar(1).eq.0) go to 10
      peru = ue-ub
      lau1 = laa
      laa = laa+4*nuest
  10  lav1 = lav
      if(ipar(2).eq.0) go to 20
      perv = ve-vb
      lav1 = laa
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! part 1: determination of the number of knots and their position.     c
! ****************************************************************     c
!  given a set of knots we compute the least-squares spline sinf(u,v), c
!  and the corresponding sum of squared residuals fp=f(p=inf).         c
!  if iopt=-1  sinf(u,v) is the requested approximation.               c
!  if iopt=0 or iopt=1 we check whether we can accept the knots:       c
!if fp <=s we will continue with the current set of knots.         c
!if fp > s we will increase the number of knots and compute the    c
!   corresponding least-squares spline until finally fp<=s.        c
!the initial choice of knots depends on the value of s and iopt.   c
!if s=0 we have spline interpolation; in that case the number of   c
!knots equals nmaxu = mu+4+2*ipar(1) and  nmaxv = mv+4+2*ipar(2)   c
!if s>0 and                                                        c
! *iopt=0 we first compute the least-squares polynomial            c
!  nu=nminu=8 and nv=nminv=8                                   c
! *iopt=1 we start with the knots found at the last call of the    c
!  routine, except for the case that s > fp0; then we can compute  c
!  the least-squares polynomial directly.                          c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  determine the number of knots for polynomial approximation.
  20  nminu = 8
      nminv = 8
      if(iopt.lt.0) go to 100
!  acc denotes the absolute tolerance for the root of f(p)=s.
      acc = tol*s
!  find nmaxu and nmaxv which denote the number of knots in u- and v-
!  direction in case of spline interpolation.
      nmaxu = mu+4+2*ipar(1)
      nmaxv = mv+4+2*ipar(2)
!  find nue and nve which denote the maximum number of knots
!  allowed in each direction
      nue = min0(nmaxu,nuest)
      nve = min0(nmaxv,nvest)
      if(s.gt.0.) go to 60
!  if s = 0, s(u,v) is an interpolating spline.
      nu = nmaxu
      nv = nmaxv
!  test whether the required storage space exceeds the available one.
      if(nv.gt.nvest .or. nu.gt.nuest) go to 420
!  find the position of the interior knots in case of interpolation.
!  the knots in the u-direction.
      nuu = nu-8
      if(nuu.eq.0) go to 40
      i = 5
      j = 3-ipar(1)
      do 30 l=1,nuu
        tu(i) = u(j)
        i = i+1
        j = j+1
  30  continue
!  the knots in the v-direction.
  40  nvv = nv-8
      if(nvv.eq.0) go to 60
      i = 5
      j = 3-ipar(2)
      do 50 l=1,nvv
        tv(i) = v(j)
        i = i+1
        j = j+1
  50  continue
      go to 100
!  if s > 0 our initial choice of knots depends on the value of iopt.
  60  if(iopt.eq.0) go to 90
      if(fp0.le.s) go to 90
!  if iopt=1 and fp0 > s we start computing the least- squares spline
!  according to the set of knots found at the last call of the routine.
!  we determine the number of grid coordinates u(i) inside each knot
!  interval (tu(l),tu(l+1)).
      l = 5
      j = 1
      nrdatu(1) = 0
      mpm = mu-1
      do 70 i=2,mpm
        nrdatu(j) = nrdatu(j)+1
        if(u(i).lt.tu(l)) go to 70
        nrdatu(j) = nrdatu(j)-1
        l = l+1
        j = j+1
        nrdatu(j) = 0
  70  continue
!  we determine the number of grid coordinates v(i) inside each knot
!  interval (tv(l),tv(l+1)).
      l = 5
      j = 1
      nrdatv(1) = 0
      mpm = mv-1
      do 80 i=2,mpm
        nrdatv(j) = nrdatv(j)+1
        if(v(i).lt.tv(l)) go to 80
        nrdatv(j) = nrdatv(j)-1
        l = l+1
        j = j+1
        nrdatv(j) = 0
  80  continue
      go to 100
!  if iopt=0 or iopt=1 and s>=fp0, we start computing the least-squares
!  polynomial (which is a spline without interior knots).
  90  nu = nminu
      nv = nminv
      nrdatu(1) = mu-2
      nrdatv(1) = mv-2
      lastdi = 0
      nplusu = 0
      nplusv = 0
      fp0 = 0.
      fpold = 0.
      reducu = 0.
      reducv = 0.
 100  mpm = mu+mv
      ifsu = 0
      ifsv = 0
      ifbu = 0
      ifbv = 0
      p = -one
!  main loop for the different sets of knots.mpm=mu+mv is a save upper
!  bound for the number of trials.
      do 250 iter=1,mpm
        if(nu.eq.nminu .and. nv.eq.nminv) ier = -2
!  find nrintu (nrintv) which is the number of knot intervals in the
!  u-direction (v-direction).
        nrintu = nu-nminu+1
        nrintv = nv-nminv+1
!  find ncof, the number of b-spline coefficients for the current set
!  of knots.
        nk1u = nu-4
        nk1v = nv-4
        ncof = nk1u*nk1v
!  find the position of the additional knots which are needed for the
!  b-spline representation of s(u,v).
        if(ipar(1).ne.0) go to 110
        i = nu
        do 105 j=1,4
          tu(j) = ub
          tu(i) = ue
          i = i-1
 105    continue
        go to 120
 110    l1 = 4
        l2 = l1
        l3 = nu-3
        l4 = l3
        tu(l2) = ub
        tu(l3) = ue
        do 115 j=1,3
          l1 = l1+1
          l2 = l2-1
          l3 = l3+1
          l4 = l4-1
          tu(l2) = tu(l4)-peru
          tu(l3) = tu(l1)+peru
 115    continue
 120    if(ipar(2).ne.0) go to 130
        i = nv
        do 125 j=1,4
          tv(j) = vb
          tv(i) = ve
          i = i-1
 125    continue
        go to 140
 130    l1 = 4
        l2 = l1
        l3 = nv-3
        l4 = l3
        tv(l2) = vb
        tv(l3) = ve
        do 135 j=1,3
          l1 = l1+1
          l2 = l2-1
          l3 = l3+1
          l4 = l4-1
          tv(l2) = tv(l4)-perv
          tv(l3) = tv(l1)+perv
 135    continue
!  find the least-squares spline sinf(u,v) and calculate for each knot
!  interval tu(j+3)<=u<=tu(j+4) (tv(j+3)<=v<=tv(j+4)) the sum
!  of squared residuals fpintu(j),j=1,2,...,nu-7 (fpintv(j),j=1,2,...
!  ,nv-7) for the data points having their absciss (ordinate)-value
!  belonging to that interval.
!  fp gives the total sum of squared residuals.
 140    call fpgrpa(ifsu,ifsv,ifbu,ifbv,idim,ipar,u,mu,v,mv,z,mz,tu,&
        nu,tv,nv,p,c,nc,fp,fpintu,fpintv,mm,mvnu,wrk(lsu),wrk(lsv),&
        wrk(lri),wrk(lq),wrk(lau),wrk(lau1),wrk(lav),wrk(lav1),&
        wrk(lbu),wrk(lbv),nru,nrv)
        if(ier.eq.(-2)) fp0 = fp
!  test whether the least-squares spline is an acceptable solution.
        if(iopt.lt.0) go to 440
        fpms = fp-s
        if(abs(fpms) .lt. acc) go to 440
!  if f(p=inf) < s, we accept the choice of knots.
        if(fpms.lt.0.) go to 300
!  if nu=nmaxu and nv=nmaxv, sinf(u,v) is an interpolating spline.
        if(nu.eq.nmaxu .and. nv.eq.nmaxv) go to 430
!  increase the number of knots.
!  if nu=nue and nv=nve we cannot further increase the number of knots
!  because of the storage capacity limitation.
        if(nu.eq.nue .and. nv.eq.nve) go to 420
        ier = 0
!  adjust the parameter reducu or reducv according to the direction
!  in which the last added knots were located.
        if(lastdi) 150,170,160
 150    reducu = fpold-fp
        go to 170
 160    reducv = fpold-fp
!  store the sum of squared residuals for the current set of knots.
 170    fpold = fp
!  find nplu, the number of knots we should add in the u-direction.
        nplu = 1
        if(nu.eq.nminu) go to 180
        npl1 = nplusu*2
        rn = nplusu
        if(reducu.gt.acc) npl1 = rn*fpms/reducu
        nplu = min0(nplusu*2,max0(npl1,nplusu/2,1))
!  find nplv, the number of knots we should add in the v-direction.
 180    nplv = 1
        if(nv.eq.nminv) go to 190
        npl1 = nplusv*2
        rn = nplusv
        if(reducv.gt.acc) npl1 = rn*fpms/reducv
        nplv = min0(nplusv*2,max0(npl1,nplusv/2,1))
 190    if(nplu-nplv) 210,200,230
 200    if(lastdi.lt.0) go to 230
 210    if(nu.eq.nue) go to 230
!  addition in the u-direction.
        lastdi = -1
        nplusu = nplu
        ifsu = 0
        do 220 l=1,nplusu
!  add a new knot in the u-direction
          call fpknot(u,mu,tu,nu,fpintu,nrdatu,nrintu,nuest,1)
!  test whether we cannot further increase the number of knots in the
!  u-direction.
          if(nu.eq.nue) go to 250
 220    continue
        go to 250
 230    if(nv.eq.nve) go to 210
!  addition in the v-direction.
        lastdi = 1
        nplusv = nplv
        ifsv = 0
        do 240 l=1,nplusv
!  add a new knot in the v-direction.
          call fpknot(v,mv,tv,nv,fpintv,nrdatv,nrintv,nvest,1)
!  test whether we cannot further increase the number of knots in the
!  v-direction.
          if(nv.eq.nve) go to 250
 240    continue
!  restart the computations with the new set of knots.
 250  continue
!  test whether the least-squares polynomial is a solution of our
!  approximation problem.
 300  if(ier.eq.(-2)) go to 440
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! part 2: determination of the smoothing spline sp(u,v)                c
! *****************************************************                c
!  we have determined the number of knots and their position. we now   c
!  compute the b-spline coefficients of the smoothing spline sp(u,v).  c
!  this smoothing spline varies with the parameter p in such a way thatc
!  f(p)=suml=1,idim(sumi=1,mu(sumj=1,mv((z(i,j,l)-sp(u(i),v(j),l))**2) c
!  is a continuous, strictly decreasing function of p. moreover the    c
!  least-squares polynomial corresponds to p=0 and the least-squares   c
!  spline to p=infinity. iteratively we then have to determine the     c
!  positive value of p such that f(p)=s. the process which is proposed c
!  here makes use of rational interpolation. f(p) is approximated by a c
!  rational function r(p)=(u*p+v)/(p+w); three values of p (p1,p2,p3)  c
!  with corresponding values of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s)c
!  are used to calculate the new value of p such that r(p)=s.          c
!  convergence is guaranteed by taking f1 > 0 and f3 < 0.              c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  initial value for p.
      p1 = 0.
      f1 = fp0-s
      p3 = -one
      f3 = fpms
      p = one
      ich1 = 0
      ich3 = 0
!  iteration process to find the root of f(p)=s.
      do 350 iter = 1,maxit
!  find the smoothing spline sp(u,v) and the corresponding sum of
!  squared residuals fp.
        call fpgrpa(ifsu,ifsv,ifbu,ifbv,idim,ipar,u,mu,v,mv,z,mz,tu,&
        nu,tv,nv,p,c,nc,fp,fpintu,fpintv,mm,mvnu,wrk(lsu),wrk(lsv),&
        wrk(lri),wrk(lq),wrk(lau),wrk(lau1),wrk(lav),wrk(lav1),&
        wrk(lbu),wrk(lbv),nru,nrv)
!  test whether the approximation sp(u,v) is an acceptable solution.
        fpms = fp-s
        if(abs(fpms).lt.acc) go to 440
!  test whether the maximum allowable number of iterations has been
!  reached.
        if(iter.eq.maxit) go to 400
!  carry out one more step of the iteration process.
        p2 = p
        f2 = fpms
        if(ich3.ne.0) go to 320
        if((f2-f3).gt.acc) go to 310
!  our initial choice of p is too large.
        p3 = p2
        f3 = f2
        p = p*con4
        if(p.le.p1) p = p1*con9 + p2*con1
        go to 350
 310    if(f2.lt.0.) ich3 = 1
 320    if(ich1.ne.0) go to 340
        if((f1-f2).gt.acc) go to 330
!  our initial choice of p is too small
        p1 = p2
        f1 = f2
        p = p/con4
        if(p3.lt.0.) go to 350
        if(p.ge.p3) p = p2*con1 + p3*con9
        go to 350
!  test whether the iteration process proceeds as theoretically
!  expected.
 330    if(f2.gt.0.) ich1 = 1
 340    if(f2.ge.f1 .or. f2.le.f3) go to 410
!  find the new value of p.
        p = fprati(p1,f1,p2,f2,p3,f3)
 350  continue
!  error codes and messages.
 400  ier = 3
      go to 440
 410  ier = 2
      go to 440
 420  ier = 1
      go to 440
 430  ier = -1
      fp = 0.
 440  return
      end
