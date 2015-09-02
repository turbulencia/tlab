      subroutine parder(tx,nx,ty,ny,c,kx,ky,nux,nuy,x,mx,y,my,z,&
       wrk,lwrk,iwrk,kwrk,ier)
!  subroutine parder evaluates on a grid (x(i),y(j)),i=1,...,mx; j=1,...
!  ,my the partial derivative ( order nux,nuy) of a bivariate spline
!  s(x,y) of degrees kx and ky, given in the b-spline representation.
!
!  calling sequence:
! call parder(tx,nx,ty,ny,c,kx,ky,nux,nuy,x,mx,y,my,z,wrk,lwrk,
!* iwrk,kwrk,ier)
!
!  input parameters:
!   tx    : real array, length nx, which contains the position of the
!   knots in the x-direction.
!   nx    : integer, giving the total number of knots in the x-direction
!   ty    : real array, length ny, which contains the position of the
!   knots in the y-direction.
!   ny    : integer, giving the total number of knots in the y-direction
!   c     : real array, length (nx-kx-1)*(ny-ky-1), which contains the
!   b-spline coefficients.
!   kx,ky : integer values, giving the degrees of the spline.
!   nux   : integer values, specifying the order of the partial
!   nuy     derivative. 0<=nux<kx, 0<=nuy<ky.
!   x     : real array of dimension (mx).
!   before entry x(i) must be set to the x co-ordinate of the
!   i-th grid point along the x-axis.
!   tx(kx+1)<=x(i-1)<=x(i)<=tx(nx-kx), i=2,...,mx.
!   mx    : on entry mx must specify the number of grid points along
!   the x-axis. mx >=1.
!   y     : real array of dimension (my).
!   before entry y(j) must be set to the y co-ordinate of the
!   j-th grid point along the y-axis.
!   ty(ky+1)<=y(j-1)<=y(j)<=ty(ny-ky), j=2,...,my.
!   my    : on entry my must specify the number of grid points along
!   the y-axis. my >=1.
!   wrk   : real array of dimension lwrk. used as workspace.
!   lwrk  : integer, specifying the dimension of wrk.
!   lwrk >= mx*(kx+1-nux)+my*(ky+1-nuy)+(nx-kx-1)*(ny-ky-1)
!   iwrk  : integer array of dimension kwrk. used as workspace.
!   kwrk  : integer, specifying the dimension of iwrk. kwrk >= mx+my.
!
!  output parameters:
!   z     : real array of dimension (mx*my).
!   on succesful exit z(my*(i-1)+j) contains the value of the
!   specified partial derivative of s(x,y) at the point
!   (x(i),y(j)),i=1,...,mx;j=1,...,my.
!   ier   : integer error flag
!ier=0 : normal return
!ier=10: invalid input data (see restrictions)
!
!  restrictions:
!   mx >=1, my >=1, 0 <= nux < kx, 0 <= nuy < ky, kwrk>=mx+my
!   lwrk>=mx*(kx+1-nux)+my*(ky+1-nuy)+(nx-kx-1)*(ny-ky-1),
!   tx(kx+1) <= x(i-1) <= x(i) <= tx(nx-kx), i=2,...,mx
!   ty(ky+1) <= y(j-1) <= y(j) <= ty(ny-ky), j=2,...,my
!
!  other subroutines required:
!fpbisp,fpbspl
!
!  references :
!de boor c : on calculating with b-splines, j. approximation theory
!        6 (1972) 50-62.
!   dierckx p. : curve and surface fitting with splines, monographs on
!        numerical analysis, oxford university press, 1993.
!
!  author :
!p.dierckx
!dept. computer science, k.u.leuven
!celestijnenlaan 200a, b-3001 heverlee, belgium.
!e-mail : Paul.Dierckx@cs.kuleuven.ac.be
!
!  latest update : march 1989
!
!  ..scalar arguments..
      integer nx,ny,kx,ky,nux,nuy,mx,my,lwrk,kwrk,ier
!  ..array arguments..
      integer iwrk(kwrk)
      real tx(nx),ty(ny),c((nx-kx-1)*(ny-ky-1)),x(mx),y(my),z(mx*my),&
       wrk(lwrk)
!  ..local scalars..
      integer i,iwx,iwy,j,kkx,kky,kx1,ky1,lx,ly,lwest,l1,l2,m,m0,m1,&
       nc,nkx1,nky1,nxx,nyy
      real ak,fac
!  ..
!  before starting computations a data check is made. if the input data
!  are invalid control is immediately repassed to the calling program.
      ier = 10
      kx1 = kx+1
      ky1 = ky+1
      nkx1 = nx-kx1
      nky1 = ny-ky1
      nc = nkx1*nky1
      if(nux.lt.0 .or. nux.ge.kx) go to 400
      if(nuy.lt.0 .or. nuy.ge.ky) go to 400
      lwest = nc +(kx1-nux)*mx+(ky1-nuy)*my
      if(lwrk.lt.lwest) go to 400
      if(kwrk.lt.(mx+my)) go to 400
      if(mx-1) 400,30,10
  10  do 20 i=2,mx
        if(x(i).lt.x(i-1)) go to 400
  20  continue
  30  if(my-1) 400,60,40
  40  do 50 i=2,my
        if(y(i).lt.y(i-1)) go to 400
  50  continue
  60  ier = 0
      nxx = nkx1
      nyy = nky1
      kkx = kx
      kky = ky
!  the partial derivative of order (nux,nuy) of a bivariate spline of
!  degrees kx,ky is a bivariate spline of degrees kx-nux,ky-nuy.
!  we calculate the b-spline coefficients of this spline
      do 70 i=1,nc
        wrk(i) = c(i)
  70  continue
      if(nux.eq.0) go to 200
      lx = 1
      do 100 j=1,nux
        ak = kkx
        nxx = nxx-1
        l1 = lx
        m0 = 1
        do 90 i=1,nxx
          l1 = l1+1
          l2 = l1+kkx
          fac = tx(l2)-tx(l1)
          if(fac.le.0.) go to 90
          do 80 m=1,nyy
            m1 = m0+nyy
            wrk(m0) = (wrk(m1)-wrk(m0))*ak/fac
            m0  = m0+1
  80      continue
  90    continue
        lx = lx+1
        kkx = kkx-1
 100  continue
 200  if(nuy.eq.0) go to 300
      ly = 1
      do 230 j=1,nuy
        ak = kky
        nyy = nyy-1
        l1 = ly
        do 220 i=1,nyy
          l1 = l1+1
          l2 = l1+kky
          fac = ty(l2)-ty(l1)
          if(fac.le.0.) go to 220
          m0 = i
          do 210 m=1,nxx
            m1 = m0+1
            wrk(m0) = (wrk(m1)-wrk(m0))*ak/fac
            m0  = m0+nky1
 210      continue
 220    continue
        ly = ly+1
        kky = kky-1
 230  continue
      m0 = nyy
      m1 = nky1
      do 250 m=2,nxx
        do 240 i=1,nyy
          m0 = m0+1
          m1 = m1+1
          wrk(m0) = wrk(m1)
 240    continue
        m1 = m1+nuy
 250  continue
!  we partition the working space and evaluate the partial derivative
 300  iwx = 1+nxx*nyy
      iwy = iwx+mx*(kx1-nux)
      call fpbisp(tx(nux+1),nx-2*nux,ty(nuy+1),ny-2*nuy,wrk,kkx,kky,&
       x,mx,y,my,z,wrk(iwx),wrk(iwy),iwrk(1),iwrk(mx+1))
 400  return
      end

