      real function dblint(tx,nx,ty,ny,c,kx,ky,xb,xe,yb,ye,wrk)
!  function dblint calculates the double integral
! / xe  / ye
!|     |      s(x,y) dx dy
!xb /  yb /
!  with s(x,y) a bivariate spline of degrees kx and ky, given in the
!  b-spline representation.
!
!  calling sequence:
! aint = dblint(tx,nx,ty,ny,c,kx,ky,xb,xe,yb,ye,wrk)
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
!   xb,xe : real values, containing the boundaries of the integration
!   yb,ye   domain. s(x,y) is considered to be identically zero out-
!   side the rectangle (tx(kx+1),tx(nx-kx))*(ty(ky+1),ty(ny-ky))
!
!  output parameters:
!   aint  : real , containing the double integral of s(x,y).
!   wrk   : real array of dimension at least (nx+ny-kx-ky-2).
!   used as working space.
!   on exit, wrk(i) will contain the integral
!        / xe
!       | ni,kx+1(x) dx , i=1,2,...,nx-kx-1
!   xb /
!   with ni,kx+1(x) the normalized b-spline defined on
!   the knots tx(i),...,tx(i+kx+1)
!   wrk(j+nx-kx-1) will contain the integral
!        / ye
!       | nj,ky+1(y) dy , j=1,2,...,ny-ky-1
!   yb /
!   with nj,ky+1(y) the normalized b-spline defined on
!   the knots ty(j),...,ty(j+ky+1)
!
!  other subroutines required: fpintb
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
!  latest update : march 1989
!
!  ..scalar arguments..
      integer nx,ny,kx,ky
      real xb,xe,yb,ye
!  ..array arguments..
      real tx(nx),ty(ny),c((nx-kx-1)*(ny-ky-1)),wrk(nx+ny-kx-ky-2)
!  ..local scalars..
      integer i,j,l,m,nkx1,nky1
      real res
!  ..
      nkx1 = nx-kx-1
      nky1 = ny-ky-1
!  we calculate the integrals of the normalized b-splines ni,kx+1(x)
      call fpintb(tx,nx,wrk,nkx1,xb,xe)
!  we calculate the integrals of the normalized b-splines nj,ky+1(y)
      call fpintb(ty,ny,wrk(nkx1+1),nky1,yb,ye)
!  calculate the integral of s(x,y)
      dblint = 0.
      do 200 i=1,nkx1
        res = wrk(i)
        if(res.eq.0.) go to 200
        m = (i-1)*nky1
        l = nkx1
        do 100 j=1,nky1
          m = m+1
          l = l+1
          dblint = dblint+res*wrk(l)*c(m)
 100    continue
 200  continue
      return
      end
