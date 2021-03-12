      subroutine bispev2(tx,nx,ty,ny,c,kx,ky,xy,mxy,z,ier)
!  subroutine bispev evaluates on a grid (x(i),y(j)),i=1,...,mx; j=1,...
!  ,my a bivariate spline s(x,y) of degrees kx and ky, given in the
!  b-spline representation.
!
!  calling sequence:
! call bispev(tx,nx,ty,ny,c,kx,ky,x,mx,y,my,z,wrk,lwrk,
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
!
!  output parameters:
!   z     : real array of dimension (mx*my).
!   on succesful exit z(my*(i-1)+j) contains the value of s(x,y)
!   at the point (x(i),y(j)),i=1,...,mx;j=1,...,my.
!   ier   : integer error flag
!ier=0 : normal return
!ier=10: invalid input data (see restrictions)
!
!  restrictions:
!   mx >=1, my >=1, lwrk>=mx*(kx+1)+my*(ky+1), kwrk>=mx+my
!   tx(kx+1) <= x(i-1) <= x(i) <= tx(nx-kx), i=2,...,mx
!   ty(ky+1) <= y(j-1) <= y(j) <= ty(ny-ky), j=2,...,my
!
!  other subroutines required:
!fpbisp,fpbspl
!
!  references :
!de boor c : on calculating with b-splines, j. approximation theory
!        6 (1972) 50-62.
!cox m.g.  : the numerical evaluation of b-splines, j. inst. maths
!        applics 10 (1972) 134-149.
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

      IMPLICIT NONE

#include "types.h"

!  ..scalar arguments..
      TINTEGER nx,ny,kx,ky,mxy,ier
! ..array arguments..
      TREAL tx(nx),ty(ny),c((nx-kx-1)*(ny-ky-1)),xy(2,mxy),z(mxy)
! ..local scalars..
!  ..
! before starting computations a data check is made. if the input data
! are invalid control is immediately repassed to the calling program.
      ier = 0
      call fpbisp2(tx,nx,ty,ny,c,kx,ky,xy,mxy,z)
      return
      end
