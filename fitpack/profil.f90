      subroutine profil(iopt,tx,nx,ty,ny,c,kx,ky,u,nu,cu,ier)
!  if iopt=0 subroutine profil calculates the b-spline coefficients of
!  the univariate spline f(y) = s(u,y) with s(x,y) a bivariate spline of
!  degrees kx and ky, given in the b-spline representation.
!  if iopt = 1 it calculates the b-spline coefficients of the univariate
!  spline g(x) = s(x,u)
!
!  calling sequence:
! call profil(iopt,tx,nx,ty,ny,c,kx,ky,u,nu,cu,ier)
!
!  input parameters:
!   iopt  : integer flag, specifying whether the profile f(y) (iopt=0)
!   or the profile g(x) (iopt=1) must be determined.
!   tx    : real array, length nx, which contains the position of the
!   knots in the x-direction.
!   nx    : integer, giving the total number of knots in the x-direction
!   ty    : real array, length ny, which contains the position of the
!   knots in the y-direction.
!   ny    : integer, giving the total number of knots in the y-direction
!   c     : real array, length (nx-kx-1)*(ny-ky-1), which contains the
!   b-spline coefficients.
!   kx,ky : integer values, giving the degrees of the spline.
!   u     : real value, specifying the requested profile.
!   tx(kx+1)<=u<=tx(nx-kx), if iopt=0.
!   ty(ky+1)<=u<=ty(ny-ky), if iopt=1.
!   nu    : on entry nu must specify the dimension of the array cu.
!   nu >= ny if iopt=0, nu >= nx if iopt=1.
!
!  output parameters:
!   cu    : real array of dimension (nu).
!   on succesful exit this array contains the b-spline
!   ier   : integer error flag
!ier=0 : normal return
!ier=10: invalid input data (see restrictions)
!
!  restrictions:
!   if iopt=0 : tx(kx+1) <= u <= tx(nx-kx), nu >=ny.
!   if iopt=1 : ty(ky+1) <= u <= ty(ny-ky), nu >=nx.
!
!  other subroutines required:
!fpbspl
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
      integer iopt,nx,ny,kx,ky,nu,ier
      real u
!  ..array arguments..
      real tx(nx),ty(ny),c((nx-kx-1)*(ny-ky-1)),cu(nu)
!  ..local scalars..
      integer i,j,kx1,ky1,l,l1,m,m0,nkx1,nky1
      real sum
!  ..local array
      real h(6)
!  ..
!  before starting computations a data check is made. if the input data
!  are invalid control is immediately repassed to the calling program.
      kx1 = kx+1
      ky1 = ky+1
      nkx1 = nx-kx1
      nky1 = ny-ky1
      ier = 10
      if(iopt.ne.0) go to 200
      if(nu.lt.ny) go to 300
      if(u.lt.tx(kx1) .or. u.gt.tx(nkx1+1)) go to 300
!  the b-splinecoefficients of f(y) = s(u,y).
      ier = 0
      l = kx1
      l1 = l+1
 110  if(u.lt.tx(l1) .or. l.eq.nkx1) go to 120
      l = l1
      l1 = l+1
      go to 110
 120  call fpbspl(tx,nx,kx,u,l,h)
      m0 = (l-kx1)*nky1+1
      do 140 i=1,nky1
        m = m0
        sum = 0.
        do 130 j=1,kx1
          sum = sum+h(j)*c(m)
          m = m+nky1
 130    continue
        cu(i) = sum
        m0 = m0+1
 140  continue
      go to 300
 200  if(nu.lt.nx) go to 300
      if(u.lt.ty(ky1) .or. u.gt.ty(nky1+1)) go to 300
!  the b-splinecoefficients of g(x) = s(x,u).
      ier = 0
      l = ky1
      l1 = l+1
 210  if(u.lt.ty(l1) .or. l.eq.nky1) go to 220
      l = l1
      l1 = l+1
      go to 210
 220  call fpbspl(ty,ny,ky,u,l,h)
      m0 = l-ky
      do 240 i=1,nkx1
        m = m0
        sum = 0.
        do 230 j=1,ky1
          sum = sum+h(j)*c(m)
          m = m+1
 230    continue
        cu(i) = sum
        m0 = m0+nky1
 240  continue
 300  return
      end

