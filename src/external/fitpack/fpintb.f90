      subroutine fpintb(t,n,bint,nk1,x,y)
!  subroutine fpintb calculates integrals of the normalized b-splines
!  nj,k+1(x) of degree k, defined on the set of knots t(j),j=1,2,...n.
!  it makes use of the formulae of gaffney for the calculation of
!  indefinite integrals of b-splines.
!
!  calling sequence:
! call fpintb(t,n,bint,nk1,x,y)
!
!  input parameters:
!t    : real array,length n, containing the position of the knots.
!n    : integer value, giving the number of knots.
!nk1  : integer value, giving the number of b-splines of degree k,
!   defined on the set of knots ,i.e. nk1 = n-k-1.
!x,y  : real values, containing the end points of the integration
!   interval.
!  output parameter:
!bint : array,length nk1, containing the integrals of the b-splines.
!  ..
!  ..scalars arguments..
      integer n,nk1
      real x,y
!  ..array arguments..
      real t(n),bint(nk1)
!  ..local scalars..
      integer i,ia,ib,it,j,j1,k,k1,l,li,lj,lk,l0,min
      real a,ak,arg,b,f,one
!  ..local arrays..
      real aint(6),h(6),h1(6)
!  initialization.
      one = 0.1e+01
      k1 = n-nk1
      ak = k1
      k = k1-1
      do 10 i=1,nk1
        bint(i) = 0.
  10  continue
!  the integration limits are arranged in increasing order.
      a = x
      b = y
      min = 0
      if(a-b) 30,160,20
  20  a = y
      b = x
      min = 1
  30  if(a.lt.t(k1)) a = t(k1)
      if(b.gt.t(nk1+1)) b = t(nk1+1)
!  using the expression of gaffney for the indefinite integral of a
!  b-spline we find that
!  bint(j) = (t(j+k+1)-t(j))*(res(j,b)-res(j,a))/(k+1)
!where for t(l) <= x < t(l+1)
!res(j,x) = 0, j=1,2,...,l-k-1
!     = 1, j=l+1,l+2,...,nk1
!     = aint(j+k-l+1), j=l-k,l-k+1,...,l
!       = sumi((x-t(j+i))*nj+i,k+1-i(x)/(t(j+k+1)-t(j+i)))
!         i=0,1,...,k
      l = k1
      l0 = l+1
!  set arg = a.
      arg = a
      do 90 it=1,2
!  search for the knot interval t(l) <= arg < t(l+1).
  40    if(arg.lt.t(l0) .or. l.eq.nk1) go to 50
        l = l0
        l0 = l+1
        go to 40
!  calculation of aint(j), j=1,2,...,k+1.
!  initialization.
  50    do 55 j=1,k1
          aint(j) = 0.
  55    continue
        aint(1) = (arg-t(l))/(t(l+1)-t(l))
        h1(1) = one
        do 70 j=1,k
!  evaluation of the non-zero b-splines of degree j at arg,i.e.
!h(i+1) = nl-j+i,j(arg), i=0,1,...,j.
          h(1) = 0.
          do 60 i=1,j
            li = l+i
            lj = li-j
            f = h1(i)/(t(li)-t(lj))
            h(i) = h(i)+f*(t(li)-arg)
            h(i+1) = f*(arg-t(lj))
  60      continue
!  updating of the integrals aint.
          j1 = j+1
          do 70 i=1,j1
            li = l+i
            lj = li-j1
            aint(i) = aint(i)+h(i)*(arg-t(lj))/(t(li)-t(lj))
            h1(i) = h(i)
  70    continue
        if(it.eq.2) go to 100
!  updating of the integrals bint
        lk = l-k
        ia = lk
        do 80 i=1,k1
          bint(lk) = -aint(i)
          lk = lk+1
  80    continue
!  set arg = b.
        arg = b
  90  continue
!  updating of the integrals bint.
 100  lk = l-k
      ib = lk-1
      do 110 i=1,k1
        bint(lk) = bint(lk)+aint(i)
        lk = lk+1
 110  continue
      if(ib.lt.ia) go to 130
      do 120 i=ia,ib
        bint(i) = bint(i)+one
 120  continue
!  the scaling factors are taken into account.
 130  f = one/ak
      do 140 i=1,nk1
        j = i+k1
        bint(i) = bint(i)*(t(j)-t(i))*f
 140  continue
!  the order of the integration limits is taken into account.
      if(min.eq.0) go to 160
      do 150 i=1,nk1
        bint(i) = -bint(i)
 150  continue
 160  return
      end
