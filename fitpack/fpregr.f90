      subroutine fpregr(iopt,x,mx,y,my,z,mz,xb,xe,yb,ye,kx,ky,s,&
       nxest,nyest,tol,maxit,nc,nx,tx,ny,ty,c,fp,fp0,fpold,reducx,&
       reducy,fpintx,fpinty,lastdi,nplusx,nplusy,nrx,nry,nrdatx,nrdaty,&
       wrk,lwrk,ier)
!  ..

      IMPLICIT NONE

#include "types.h"

!  ..scalar arguments..
      TREAL xb,xe,yb,ye,s,tol,fp,fp0,fpold,reducx,reducy
      TINTEGER iopt,mx,my,mz,kx,ky,nxest,nyest,maxit,nc,nx,ny,lastdi,&
       nplusx,nplusy,lwrk,ier
!  ..array arguments..
      TREAL x(mx),y(my),z(mz),tx(nxest),ty(nyest),c(nc),fpintx(nxest),&
       fpinty(nyest),wrk(lwrk)
      TINTEGER nrdatx(nxest),nrdaty(nyest),nrx(mx),nry(my)
!  ..local scalars
      TREAL acc,fpms,f1,f2,f3,p,p1,p2,p3,rn,one,half,con1,con9,con4
      TINTEGER i,ich1,ich3,ifbx,ifby,ifsx,ifsy,iter,j,kx1,kx2,ky1,ky2,&
       k3,l,lax,lay,lbx,lby,lq,lri,lsx,lsy,mk1,mm,mpm,mynx,ncof,&
       nk1x,nk1y,nmaxx,nmaxy,nminx,nminy,nplx,nply,npl1,nrintx,&
       nrinty,nxe,nxk,nye
!  ..function references..
      TREAL fprati
!  TINTEGER max0,min0
!  ..subroutine references..
!fpgrre,fpknot
!  ..
!   set constants
      one = C_1_R
      half = C_05_R
      con1 = C_01_R
#ifdef SINGLE_PREC
      con9 = 0.9e0
      con4 = 0.4e-01
#else
      con9 = 0.9d0
      con4 = 0.4d-01
#endif
     
!  we partition the working space.
      kx1 = kx+1
      ky1 = ky+1
      kx2 = kx1+1
      ky2 = ky1+1
      lsx = 1
      lsy = lsx+mx*kx1
      lri = lsy+my*ky1
      mm = max0(nxest,my)
      lq = lri+mm
      mynx = nxest*my
      lax = lq+mynx
      nxk = nxest*kx2
      lbx = lax+nxk
      lay = lbx+nxk
      lby = lay+nyest*ky2
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! part 1: determination of the number of knots and their position.     c
! ****************************************************************     c
!  given a set of knots we compute the least-squares spline sinf(x,y), c
!  and the corresponding sum of squared residuals fp=f(p=inf).         c
!  if iopt=-1  sinf(x,y) is the requested approximation.               c
!  if iopt=0 or iopt=1 we check whether we can accept the knots:       c
!if fp <=s we will continue with the current set of knots.         c
!if fp > s we will increase the number of knots and compute the    c
!   corresponding least-squares spline until finally fp<=s.        c
!the initial choice of knots depends on the value of s and iopt.   c
!if s=0 we have spline interpolation; in that case the number of   c
!knots equals nmaxx = mx+kx+1  and  nmaxy = my+ky+1.               c
!if s>0 and                                                        c
! *iopt=0 we first compute the least-squares polynomial of degree  c
!  kx in x and ky in y; nx=nminx=2*kx+2 and ny=nymin=2*ky+2.       c
! *iopt=1 we start with the knots found at the last call of the    c
!  routine, except for the case that s > fp0; then we can compute  c
!  the least-squares polynomial directly.                          c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  determine the number of knots for polynomial approximation.
      nminx = 2*kx1
      nminy = 2*ky1
      if(iopt.lt.0) go to 120
!  acc denotes the absolute tolerance for the root of f(p)=s.
      acc = tol*s
!  find nmaxx and nmaxy which denote the number of knots in x- and y-
!  direction in case of spline interpolation.
      nmaxx = mx+kx1
      nmaxy = my+ky1
!  find nxe and nye which denote the maximum number of knots
!  allowed in each direction
      nxe = min0(nmaxx,nxest)
      nye = min0(nmaxy,nyest)
      if(s.gt.0.) go to 100
!  if s = 0, s(x,y) is an interpolating spline.
      nx = nmaxx
      ny = nmaxy
!  test whether the required storage space exceeds the available one.
      if(ny.gt.nyest .or. nx.gt.nxest) go to 420
!  find the position of the interior knots in case of interpolation.
!  the knots in the x-direction.
      mk1 = mx-kx1
      if(mk1.eq.0) go to 60
      k3 = kx/2
      i = kx1+1
      j = k3+2
      if(k3*2.eq.kx) go to 40
      do 30 l=1,mk1
        tx(i) = x(j)
        i = i+1
        j = j+1
  30  continue
      go to 60
  40  do 50 l=1,mk1
        tx(i) = (x(j)+x(j-1))*half
        i = i+1
        j = j+1
  50  continue
!  the knots in the y-direction.
  60  mk1 = my-ky1
      if(mk1.eq.0) go to 120
      k3 = ky/2
      i = ky1+1
      j = k3+2
      if(k3*2.eq.ky) go to 80
      do 70 l=1,mk1
        ty(i) = y(j)
        i = i+1
        j = j+1
  70  continue
      go to 120
  80  do 90 l=1,mk1
        ty(i) = (y(j)+y(j-1))*half
        i = i+1
        j = j+1
  90  continue
      go to 120
!  if s > 0 our initial choice of knots depends on the value of iopt.
 100  if(iopt.eq.0) go to 115
      if(fp0.le.s) go to 115
!  if iopt=1 and fp0 > s we start computing the least- squares spline
!  according to the set of knots found at the last call of the routine.
!  we determine the number of grid coordinates x(i) inside each knot
!  interval (tx(l),tx(l+1)).
      l = kx2
      j = 1
      nrdatx(1) = 0
      mpm = mx-1
      do 105 i=2,mpm
        nrdatx(j) = nrdatx(j)+1
        if(x(i).lt.tx(l)) go to 105
        nrdatx(j) = nrdatx(j)-1
        l = l+1
        j = j+1
        nrdatx(j) = 0
 105  continue
!  we determine the number of grid coordinates y(i) inside each knot
!  interval (ty(l),ty(l+1)).
      l = ky2
      j = 1
      nrdaty(1) = 0
      mpm = my-1
      do 110 i=2,mpm
        nrdaty(j) = nrdaty(j)+1
        if(y(i).lt.ty(l)) go to 110
        nrdaty(j) = nrdaty(j)-1
        l = l+1
        j = j+1
        nrdaty(j) = 0
 110  continue
      go to 120
!  if iopt=0 or iopt=1 and s>=fp0, we start computing the least-squares
!  polynomial of degree kx in x and ky in y (which is a spline without
!  interior knots).
 115  nx = nminx
      ny = nminy
      nrdatx(1) = mx-2
      nrdaty(1) = my-2
      lastdi = 0
      nplusx = 0
      nplusy = 0
      fp0 = 0.
      fpold = 0.
      reducx = 0.
      reducy = 0.
 120  mpm = mx+my
      ifsx = 0
      ifsy = 0
      ifbx = 0
      ifby = 0
      p = -one
!  main loop for the different sets of knots.mpm=mx+my is a save upper
!  bound for the number of trials.
      do 250 iter=1,mpm
        if(nx.eq.nminx .and. ny.eq.nminy) ier = -2
!  find nrintx (nrinty) which is the number of knot intervals in the
!  x-direction (y-direction).
        nrintx = nx-nminx+1
        nrinty = ny-nminy+1
!  find ncof, the number of b-spline coefficients for the current set
!  of knots.
        nk1x = nx-kx1
        nk1y = ny-ky1
        ncof = nk1x*nk1y
!  find the position of the additional knots which are needed for the
!  b-spline representation of s(x,y).
        i = nx
        do 130 j=1,kx1
          tx(j) = xb
          tx(i) = xe
          i = i-1
 130    continue
        i = ny
        do 140 j=1,ky1
          ty(j) = yb
          ty(i) = ye
          i = i-1
 140    continue
!  find the least-squares spline sinf(x,y) and calculate for each knot
!  interval tx(j+kx)<=x<=tx(j+kx+1) (ty(j+ky)<=y<=ty(j+ky+1)) the sum
!  of squared residuals fpintx(j),j=1,2,...,nx-2*kx-1 (fpinty(j),j=1,2,
!  ...,ny-2*ky-1) for the data points having their absciss (ordinate)-
!  value belonging to that interval.
!  fp gives the total sum of squared residuals.
        call fpgrre(ifsx,ifsy,ifbx,ifby,x,mx,y,my,z,mz,kx,ky,tx,nx,ty,&
        ny,p,c,nc,fp,fpintx,fpinty,mm,mynx,kx1,kx2,ky1,ky2,wrk(lsx),&
        wrk(lsy),wrk(lri),wrk(lq),wrk(lax),wrk(lay),wrk(lbx),wrk(lby),&
        nrx,nry)
        if(ier.eq.(-2)) fp0 = fp
!  test whether the least-squares spline is an acceptable solution.
        if(iopt.lt.0) go to 440
        fpms = fp-s
        if(abs(fpms) .lt. acc) go to 440
!  if f(p=inf) < s, we accept the choice of knots.
        if(fpms.lt.0.) go to 300
!  if nx=nmaxx and ny=nmaxy, sinf(x,y) is an interpolating spline.
        if(nx.eq.nmaxx .and. ny.eq.nmaxy) go to 430
!  increase the number of knots.
!  if nx=nxe and ny=nye we cannot further increase the number of knots
!  because of the storage capacity limitation.
        if(nx.eq.nxe .and. ny.eq.nye) go to 420
        ier = 0
!  adjust the parameter reducx or reducy according to the direction
!  in which the last added knots were located.
        if(lastdi) 150,170,160
 150    reducx = fpold-fp
        go to 170
 160    reducy = fpold-fp
!  store the sum of squared residuals for the current set of knots.
 170    fpold = fp
!  find nplx, the number of knots we should add in the x-direction.
        nplx = 1
        if(nx.eq.nminx) go to 180
        npl1 = nplusx*2
        rn = nplusx
        if(reducx.gt.acc) npl1 = rn*fpms/reducx
        nplx = min0(nplusx*2,max0(npl1,nplusx/2,1))
!  find nply, the number of knots we should add in the y-direction.
 180    nply = 1
        if(ny.eq.nminy) go to 190
        npl1 = nplusy*2
        rn = nplusy
        if(reducy.gt.acc) npl1 = rn*fpms/reducy
        nply = min0(nplusy*2,max0(npl1,nplusy/2,1))
 190    if(nplx-nply) 210,200,230
 200    if(lastdi.lt.0) go to 230
 210    if(nx.eq.nxe) go to 230
!  addition in the x-direction.
        lastdi = -1
        nplusx = nplx
        ifsx = 0
        do 220 l=1,nplusx
!  add a new knot in the x-direction
          call fpknot(x,mx,tx,nx,fpintx,nrdatx,nrintx,nxest,1)
!  test whether we cannot further increase the number of knots in the
!  x-direction.
          if(nx.eq.nxe) go to 250
 220    continue
        go to 250
 230    if(ny.eq.nye) go to 210
!  addition in the y-direction.
        lastdi = 1
        nplusy = nply
        ifsy = 0
        do 240 l=1,nplusy
!  add a new knot in the y-direction.
          call fpknot(y,my,ty,ny,fpinty,nrdaty,nrinty,nyest,1)
!  test whether we cannot further increase the number of knots in the
!  y-direction.
          if(ny.eq.nye) go to 250
 240    continue
!  restart the computations with the new set of knots.
 250  continue
!  test whether the least-squares polynomial is a solution of our
!  approximation problem.
 300  if(ier.eq.(-2)) go to 440
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! part 2: determination of the smoothing spline sp(x,y)                c
! *****************************************************                c
!  we have determined the number of knots and their position. we now   c
!  compute the b-spline coefficients of the smoothing spline sp(x,y).  c
!  this smoothing spline varies with the parameter p in such a way thatc
!f(p) = sumi=1,mx(sumj=1,my((z(i,j)-sp(x(i),y(j)))**2)             c
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
!  find the smoothing spline sp(x,y) and the corresponding sum of
!  squared residuals fp.
        call fpgrre(ifsx,ifsy,ifbx,ifby,x,mx,y,my,z,mz,kx,ky,tx,nx,ty,&
        ny,p,c,nc,fp,fpintx,fpinty,mm,mynx,kx1,kx2,ky1,ky2,wrk(lsx),&
        wrk(lsy),wrk(lri),wrk(lq),wrk(lax),wrk(lay),wrk(lbx),wrk(lby),&
        nrx,nry)
!  test whether the approximation sp(x,y) is an acceptable solution.
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

