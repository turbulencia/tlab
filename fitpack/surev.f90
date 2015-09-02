      subroutine surev(idim,tu,nu,tv,nv,c,u,mu,v,mv,f,mf,wrk,lwrk,&
       iwrk,kwrk,ier)
!  subroutine surev evaluates on a grid (u(i),v(j)),i=1,...,mu; j=1,...
!  ,mv a bicubic spline surface of dimension idim, given in the
!  b-spline representation.
!
!  calling sequence:
! call surev(idim,tu,nu,tv,nv,c,u,mu,v,mv,f,mf,wrk,lwrk,
!* iwrk,kwrk,ier)
!
!  input parameters:
!   idim  : integer, specifying the dimension of the spline surface.
!   tu    : real array, length nu, which contains the position of the
!   knots in the u-direction.
!   nu    : integer, giving the total number of knots in the u-direction
!   tv    : real array, length nv, which contains the position of the
!   knots in the v-direction.
!   nv    : integer, giving the total number of knots in the v-direction
!   c     : real array, length (nu-4)*(nv-4)*idim, which contains the
!   b-spline coefficients.
!   u     : real array of dimension (mu).
!   before entry u(i) must be set to the u co-ordinate of the
!   i-th grid point along the u-axis.
!   tu(4)<=u(i-1)<=u(i)<=tu(nu-3), i=2,...,mu.
!   mu    : on entry mu must specify the number of grid points along
!   the u-axis. mu >=1.
!   v     : real array of dimension (mv).
!   before entry v(j) must be set to the v co-ordinate of the
!   j-th grid point along the v-axis.
!   tv(4)<=v(j-1)<=v(j)<=tv(nv-3), j=2,...,mv.
!   mv    : on entry mv must specify the number of grid points along
!   the v-axis. mv >=1.
!   mf    : on entry, mf must specify the dimension of the array f.
!   mf >= mu*mv*idim
!   wrk   : real array of dimension lwrk. used as workspace.
!   lwrk  : integer, specifying the dimension of wrk.
!   lwrk >= 4*(mu+mv)
!   iwrk  : integer array of dimension kwrk. used as workspace.
!   kwrk  : integer, specifying the dimension of iwrk. kwrk >= mu+mv.
!
!  output parameters:
!   f     : real array of dimension (mf).
!   on succesful exit f(mu*mv*(l-1)+mv*(i-1)+j) contains the
!   l-th co-ordinate of the bicubic spline surface at the
!   point (u(i),v(j)),l=1,...,idim,i=1,...,mu;j=1,...,mv.
!   ier   : integer error flag
!ier=0 : normal return
!ier=10: invalid input data (see restrictions)
!
!  restrictions:
!   mu >=1, mv >=1, lwrk>=4*(mu+mv), kwrk>=mu+mv , mf>=mu*mv*idim
!   tu(4) <= u(i-1) <= u(i) <= tu(nu-3), i=2,...,mu
!   tv(4) <= v(j-1) <= v(j) <= tv(nv-3), j=2,...,mv
!
!  other subroutines required:
!fpsuev,fpbspl
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
!
!  ..scalar arguments..
      integer idim,nu,nv,mu,mv,mf,lwrk,kwrk,ier
!  ..array arguments..
      integer iwrk(kwrk)
      real tu(nu),tv(nv),c((nu-4)*(nv-4)*idim),u(mu),v(mv),f(mf),&
       wrk(lwrk)
!  ..local scalars..
      integer i,muv
!  ..
!  before starting computations a data check is made. if the input data
!  are invalid control is immediately repassed to the calling program.
      ier = 10
      if(mf.lt.mu*mv*idim) go to 100
      muv = mu+mv
      if(lwrk.lt.4*muv) go to 100
      if(kwrk.lt.muv) go to 100
      if(mu-1) 100,30,10
  10  do 20 i=2,mu
        if(u(i).lt.u(i-1)) go to 100
  20  continue
  30  if(mv-1) 100,60,40
  40  do 50 i=2,mv
        if(v(i).lt.v(i-1)) go to 100
  50  continue
  60  ier = 0
      call fpsuev(idim,tu,nu,tv,nv,c,u,mu,v,mv,f,wrk(1),wrk(4*mu+1),&
       iwrk(1),iwrk(mu+1))
 100  return
      end
